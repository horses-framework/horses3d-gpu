#include "Includes.h"

module ActuatorDisk
!
!  ///////////////////////////////////////////////////////////////////////////////////////
!
!  Actuator disk model for ducted-fan / BLI simulations.
!
!  This module extends the actuator-line (AL) infrastructure of
!  Solver/src/libs/sources/ActuatorLine.f90 to a stationary actuator disk, following the
!  design note "An Actuator Disk Model for Ducted-Fan Simulations in HORSES3D-GPU". Two
!  forcing strategies are implemented, selected with the control-file key
!  "actuator disk model":
!
!    - "pressure"     (Model A) a simple disk with a prescribed pressure/enthalpy rise
!                      (geometry-only radial force law, no blade/airfoil data), Eq. (18).
!    - "bladeelement"  (Model B) a local blade-element disk: Cl, Cd are evaluated from the
!                      locally sampled flow exactly as in the AL blade-element solve, but the
!                      resulting force is smeared uniformly in the azimuthal direction
!                      instead of following a small number of discrete rotating blades,
!                      Eq. (19).
!
!  As in the AL model, the disk is realized as a volumetric momentum source deposited with
!  an isotropic Gaussian kernel (Sorensen & Shen). Unlike the AL model, the disk points are
!  stationary in space (they do not rotate), so the per-partition element search performed
!  once in ConstructDiskFarm remains valid for the whole run: no per-time-step point
!  relocation/MPI dance is required. In addition, a volumetric ENERGY source is deposited
!  (added to S_NS(IRHOE,...) under NAVIERSTOKES), consistent with the momentum source: the
!  power delivered to the fluid at each actuator point is simply the dot product of the
!  local force (already expressed as "force on the fluid") with the locally sampled
!  velocity, Eq. (16). This term is required for a fan (which adds energy to the flow),
!  unlike the wind-turbine extraction case the AL model was written for.
!
!  Sign convention: both models return, at every ring point, the AXIAL and TANGENTIAL
!  (theta) force delivered TO THE FLUID (not to a blade):
!    - Model A: a positive "delta_p0" (prescribed pressure/stagnation-enthalpy rise) always
!      pushes the fluid forward (+x, i.e. increases its axial momentum), matching the
!      classical actuator-disk convention for a fan/propeller.
!    - Model B: the same blade-element equations as FarmGetLocalForces are used verbatim
!      (lift/drag on a virtual blade element), and the force on the fluid is the reaction
!      of that force, exactly as in the AL model. Whether the disk then behaves as a fan
!      (adds energy) or a turbine (extracts energy) follows from the sign of the blade
!      twist/pitch supplied in the actuator-disk definition file, exactly as for any BEM
!      code; the equations themselves are generic.
!
!  Only the axis-aligned configuration already assumed by the AL model is supported: the
!  disk axis is the +x direction, and the rotor/azimuthal plane is y-z (the "normal_x/y/z"
!  entry is kept for interface completeness/future generalization but is not currently used
!  to rotate the geometry, exactly mirroring the AL model's own convention).
!
!  ///////////////////////////////////////////////////////////////////////////////////////
!
#if defined(NAVIERSTOKES) || defined(INCNS) || defined(MULTIPHASE)
    use SMConstants
    use MPI_Process_Info
#ifdef _HAS_MPI_
    use mpi
#endif
    implicit none

private
public diskFarm, ConstructDiskFarm, DestructDiskFarm, UpdateDiskFarm, ForcesDiskFarm, WriteDiskFarmForces
!
!   Model selector
!   --------------
    integer, parameter :: AD_MODEL_PRESSURE      = 0   ! Model A: prescribed pressure/enthalpy rise
    integer, parameter :: AD_MODEL_BLADE_ELEMENT = 1   ! Model B: local blade-element (Cl/Cd)
    integer, parameter :: MAX_AIRFOIL_FILES_AD   = 10
!
!  ****************************************
!  DEFINE AIRFOIL AND DISK (ACTUATOR RING)
!  ****************************************
!
    type airfoil_disk_t
    integer                         :: num_aoa
    integer                         :: num_Re
    real(KIND=RP), allocatable      :: aoa(:)  ! in rad
    real(KIND=RP), allocatable      :: Re(:)   ! reynolds number
    real(KIND=RP), allocatable      :: cl(:,:)
    real(KIND=RP), allocatable      :: cd(:,:)
    end type

    type disk_t
!
!       geometry
!       --------
        real(KIND=RP)                   :: hub_cood_x, hub_cood_y, hub_cood_z
        real(KIND=RP)                   :: normal_x, normal_y, normal_z  ! kept for interface completeness, see header note
        real(KIND=RP)                   :: radius, hub_radius
        integer                         :: n_radial, n_azimuthal
        real(KIND=RP), allocatable      :: r_R(:)              ! (0:n_radial), r_R(0) = hub_radius
        real(KIND=RP), allocatable      :: theta(:)            ! (n_azimuthal), fixed azimuthal positions
        real(KIND=RP), allocatable      :: point_xyz_loc(:,:,:)! (n_radial,n_azimuthal,3), fixed in space
!
!       Model A - prescribed pressure/enthalpy rise, size (n_radial)
!       -------------------------------------------------------------
        real(KIND=RP), allocatable      :: delta_p0(:)
        real(KIND=RP), allocatable      :: swirl_q(:)
        logical                         :: use_swirl = .false.
!
!       Model B - local blade-element (Cl/Cd), size (n_radial) unless noted
!       --------------------------------------------------------------------
        integer                         :: num_blades
        real(KIND=RP)                   :: rot_speed        ! rad/s
        real(KIND=RP)                   :: blade_pitch      ! rad
        real(KIND=RP), allocatable      :: chord(:)         ! m
        real(KIND=RP), allocatable      :: twist(:)         ! rad
        integer, allocatable            :: num_airfoils(:)
        CHARACTER(LEN=30), allocatable  :: airfoil_files(:,:)
        type(airfoil_disk_t), allocatable :: airfoil_t(:)
        real(KIND=RP)                   :: tip_c1, tip_c2
        logical                         :: compressibility_correction = .false.
!
!       numerics
!       --------
        real(KIND=RP)                   :: gauss_epsil
!
!       point search (performed once; points are stationary)
!       -----------------------------------------------------
        integer, allocatable            :: eID(:,:)             ! (n_radial,n_azimuthal)
        real(KIND=RP), allocatable      :: xi_loc(:,:,:)         ! (n_radial,n_azimuthal,3)
        logical, allocatable            :: elementFound(:,:)
!
!       per-time-step sampled state and forces (non-projection strategy)
!       -------------------------------------------------------------------
        real(KIND=RP), allocatable      :: Q_disk(:,:,:)         ! (n_radial,n_azimuthal,NCONS)
        real(KIND=RP), allocatable      :: Q_disk_all(:,:,:)
!
!       per-point forces/diagnostics, size (n_radial,n_azimuthal)
!       -------------------------------------------------------------
        real(KIND=RP), allocatable      :: local_thrust(:,:)       ! axial force ON THE FLUID [N]
        real(KIND=RP), allocatable      :: local_rotor_force(:,:)  ! tangential force ON THE FLUID [N]
        real(KIND=RP), allocatable      :: local_power(:,:)        ! power delivered TO THE FLUID [W]
        real(KIND=RP), allocatable      :: local_velocity(:,:)     ! local relative velocity [m/s]
        real(KIND=RP), allocatable      :: local_angle(:,:)        ! local flow angle [rad]
        real(KIND=RP), allocatable      :: local_aoa(:,:)          ! local angle of attack [rad] (0 for Model A)
        real(KIND=RP), allocatable      :: local_Re(:,:)           ! local Reynolds number
!
!       projection-strategy accumulators (weighted sums, reset every UpdateDiskFarm call)
!       -------------------------------------------------------------------------------------
        real(KIND=RP), allocatable      :: local_thrust_temp(:,:)
        real(KIND=RP), allocatable      :: local_rotor_force_temp(:,:)
        real(KIND=RP), allocatable      :: local_power_temp(:,:)
        real(KIND=RP), allocatable      :: local_velocity_temp(:,:)
        real(KIND=RP), allocatable      :: local_angle_temp(:,:)
        real(KIND=RP), allocatable      :: local_aoa_temp(:,:)
        real(KIND=RP), allocatable      :: local_Re_temp(:,:)
        real(KIND=RP), allocatable      :: local_gaussian_sum(:,:)
!
!       whole-disk integrated outputs
!       ------------------------------
        real(KIND=RP)                   :: total_thrust   ! N
        real(KIND=RP)                   :: total_torque   ! N.m
        real(KIND=RP)                   :: total_power    ! W
        real(KIND=RP)                   :: Cp, Ct
!
!       time-averaged radial profile (r, U, AoA, Re, tangential force, axial force), azimuth index 1
!       -------------------------------------------------------------------------------------------------
        real(KIND=RP), allocatable      :: average_conditions(:,:)
    end type

    type DiskFarm_t
    integer                        :: num_disks = 0
    integer                        :: model = AD_MODEL_PRESSURE
    type(disk_t), allocatable      :: disk_t(:)
    real(KIND=RP)                  :: gauss_epsil
    real(KIND=RP)                  :: time
    real(KIND=RP)                  :: tolerance_factor
    integer                        :: epsilon_type
    logical                        :: calculate_with_projection = .false.
    logical                        :: active = .false.
    logical                        :: save_average = .false.
    logical                        :: save_instant = .false.
    logical                        :: verbose = .false.
    character(len=LINE_LENGTH)     :: file_name
    integer                        :: number_iterations
    integer                        :: save_iterations
    end type

    type(DiskFarm_t)                :: diskFarm

    integer, dimension(:), allocatable  :: elementsActuatedAD, diskOfElement, numElementsPerDiskAD
    !$acc declare create(elementsActuatedAD, diskOfElement)

!  ========
contains
!  ========
!
!///////////////////////////////////////////////////////////////////////////////////////
!
   subroutine ConstructDiskFarm(self, controlVariables, t0, mesh)
       use FTValueDictionaryClass
       use mainKeywordsModule, only: solutionFileNameKey
       use FileReadingUtilities      , only: getFileName
       use PhysicsStorage
       use HexMeshClass
       use fluiddata
       use MPI_Process_Info
       implicit none
       type(DiskFarm_t) , intent(inout)             :: self
       TYPE(FTValueDictionary), intent(in)          :: controlVariables
       real(kind=RP), intent(in)                    :: t0
       type(HexMesh), intent(in)                    :: mesh
!        ---------------
!        Local variables
!        ---------------
!
         integer     ::  i, j, k, ii, fid, n_aoa, n_airfoil
         character(LEN=LINE_LENGTH) :: arg, char1, model_str
         character(LEN=LINE_LENGTH) :: solution_file
         character(LEN=5)           :: file_id
         integer        :: nelem, eID, eIndex, tip_flag_int
         real(kind=RP)  :: tolerance, r_square
         real(kind=RP)  :: delta, delta_temp, tip_c1_common, tip_c2_common, dtheta
         real(kind=RP), dimension(NDIM)  :: x, xi
         integer        :: delta_count, delta_paritions, ierr
         logical                    :: found, allfound

    if (.not. controlVariables % logicalValueForKey("use actuatordisk")) return

    self % time = t0

    self % epsilon_type = controlVariables % getValueOrDefault("actuator disk epsilon type", 0)
    self % calculate_with_projection = controlVariables % getValueOrDefault("actuator disk calculate with projection", .false.)
    self % save_average = controlVariables % getValueOrDefault("actuator disk save average", .false.)
    self % save_instant = controlVariables % getValueOrDefault("actuator disk save instant", .false.)
    self % save_iterations = controlVariables % getValueOrDefault("actuator disk save iteration", 1)
    self % verbose = controlVariables % getValueOrDefault("actuator disk verbose", .false.)
!
!   For a thin disk, "tolerance" also sets the axial thickness of the selected element
!   region (see design note, Sec. 6); a smaller default than the AL model (0.2) is used.
!   ------------------------------------------------------------------------------------
    self % tolerance_factor = controlVariables % getValueOrDefault("actuator disk tolerance", 0.1_RP)
!
!   Model selector: "pressure" (A, default) or "bladeelement" (B)
!   ---------------------------------------------------------------
    if (controlVariables % containsKey("actuator disk model")) then
        model_str = controlVariables % stringValueForKey("actuator disk model", requestedLength = LINE_LENGTH)
    else
        model_str = "pressure"
    end if
    model_str = adjustl(model_str)
    do i = 1, len_trim(model_str)
        if (model_str(i:i) >= 'A' .and. model_str(i:i) <= 'Z') model_str(i:i) = achar(iachar(model_str(i:i))+32)
    end do
    if ( index(trim(model_str),"blade") .gt. 0 .or. trim(model_str) .eq. "b" ) then
        self % model = AD_MODEL_BLADE_ELEMENT
    else
        self % model = AD_MODEL_PRESSURE
    end if

    arg='./ActuatorDef/Act_ActuatorDiskDef.dat'
    OPEN( newunit = fid,file=trim(arg),status="old",action="read")

    READ(fid,'(A132)') char1
    READ(fid,'(A132)') char1
    READ(fid,'(A132)') char1
    READ(fid,'(A132)') char1
    READ(fid,*) self%num_disks

    if (self % verbose .and. MPI_Process % isRoot) then
        print *,'-------------------------'
        print *,achar(27)//'[34m READING ACTUATOR DISK FARM DEFINITION'
        write(*,*) "Number of disks in farm:", self%num_disks
        select case (self % model)
        case (AD_MODEL_PRESSURE)
            write(*,*) "Model: A (prescribed pressure/enthalpy rise)"
        case (AD_MODEL_BLADE_ELEMENT)
            write(*,*) "Model: B (local blade-element)"
        end select
    endif

    allocate(self%disk_t(self%num_disks))

    READ(fid,'(A132)') char1
    do k = 1, self%num_disks
       READ(fid,*) self%disk_t(k)%hub_cood_x, self%disk_t(k)%hub_cood_y, self%disk_t(k)%hub_cood_z
    end do

    READ(fid,'(A132)') char1
    do k = 1, self%num_disks
       READ(fid,*) self%disk_t(k)%normal_x, self%disk_t(k)%normal_y, self%disk_t(k)%normal_z
    end do

    READ(fid,'(A132)') char1
    do k = 1, self%num_disks
       READ(fid,*) self%disk_t(k)%radius
    end do

    READ(fid,'(A132)') char1
    do k = 1, self%num_disks
       READ(fid,*) self%disk_t(k)%hub_radius
    end do

    READ(fid,'(A132)') char1
    do k = 1, self%num_disks
       READ(fid,*) self%disk_t(k)%n_radial, self%disk_t(k)%n_azimuthal
    end do

    READ(fid,'(A132)') char1
    READ(fid,'(A132)') char1
    READ(fid,'(A132)') char1

    select case (self % model)
    case (AD_MODEL_PRESSURE)
!
!      ----------------------------------------------
!      Model A: r, delta_p0(r), swirl_q(r) per station
!      ----------------------------------------------
!
       do k = 1, self%num_disks
         associate ( nr => self%disk_t(k)%n_radial )
            allocate( self%disk_t(k)%r_R(0:nr), self%disk_t(k)%delta_p0(nr), self%disk_t(k)%swirl_q(nr) )
            self%disk_t(k)%r_R(0) = self%disk_t(k)%hub_radius
            READ(fid,'(A132)') char1  ! per-disk comment line
            do i = 1, nr
               READ(fid,*) self%disk_t(k)%r_R(i), self%disk_t(k)%delta_p0(i), self%disk_t(k)%swirl_q(i)
            end do
            self%disk_t(k)%use_swirl = any(abs(self%disk_t(k)%swirl_q) .gt. 0.0_RP)
         end associate
       end do

    case default ! AD_MODEL_BLADE_ELEMENT
!
!      -----------------------------------------------------------------
!      Model B: num_blades, rot_speed, blade_pitch, r-chord-twist table
!      -----------------------------------------------------------------
!
       READ(fid,'(A132)') char1
       do k = 1, self%num_disks
          READ(fid,*) self%disk_t(k)%num_blades
       end do

       READ(fid,'(A132)') char1
       do k = 1, self%num_disks
          READ(fid,*) self%disk_t(k)%rot_speed, self%disk_t(k)%blade_pitch
       end do

       READ(fid,'(A132)') char1
       do k = 1, self%num_disks
         associate ( nr => self%disk_t(k)%n_radial )
            allocate( self%disk_t(k)%r_R(0:nr), self%disk_t(k)%chord(nr), self%disk_t(k)%twist(nr), &
                      self%disk_t(k)%num_airfoils(nr), self%disk_t(k)%airfoil_files(nr,MAX_AIRFOIL_FILES_AD), &
                      self%disk_t(k)%airfoil_t(nr) )
            self%disk_t(k)%r_R(0) = self%disk_t(k)%hub_radius
            do i = 1, nr
               self%disk_t(k)%airfoil_files(i,:) = ' '
            end do

            READ(fid,'(A132)') char1  ! per-disk comment line
            do i = 1, nr
               READ(fid,*) self%disk_t(k)%r_R(i), self%disk_t(k)%chord(i), self%disk_t(k)%twist(i), self%disk_t(k)%num_airfoils(i)
               self%disk_t(k)%airfoil_t(i)%num_Re = self%disk_t(k)%num_airfoils(i)
               do n_airfoil = 1, self%disk_t(k)%num_airfoils(i)
                  READ(fid,*) self%disk_t(k)%airfoil_files(i,n_airfoil)
               end do
            end do
         end associate
       end do

    end select

! read numerical parameters
     READ(fid,'(A132)') char1
     READ(fid,'(A132)') char1
     READ(fid,'(A132)') char1
     READ(fid,'(A132)') char1

     READ(fid,*) self%gauss_epsil

     if (self % model .eq. AD_MODEL_BLADE_ELEMENT) then
        READ(fid,'(A132)') char1
        READ(fid,*) tip_c1_common, tip_c2_common
        do k = 1, self%num_disks
           self%disk_t(k)%tip_c1 = tip_c1_common
           self%disk_t(k)%tip_c2 = tip_c2_common
        end do

        READ(fid,'(A132)') char1
        READ(fid,*) tip_flag_int
        do k = 1, self%num_disks
           self%disk_t(k)%compressibility_correction = (tip_flag_int .ne. 0)
        end do
     end if

    close(fid)

    if (self % verbose .and. MPI_Process % isRoot) then
        print *,achar(27)//'[34m END OF READING ACTUATOR DISK FARM DEFINITION'
    end if

    if (MPI_Process % isRoot) then
        call Subsection_Header("Actuator Disk")
        write(STD_OUT,'(30X,A,A28,I0)') "->", "Number of disks: ", self % num_disks
        select case (self % epsilon_type)
        case (0)
            write(STD_OUT,'(30X,A,A28,ES10.3)') "->", 'Fixed Epsilon value: ',self%gauss_epsil
        case (2)
            write(STD_OUT,'(30X,A)') 'Epsilon calculated based on element size and polynomial order'
            write(STD_OUT,'(30X,A,A28,F10.3)') "->", 'Constant for Epsilon: ',self%gauss_epsil
        case default
            write(STD_OUT,'(30X,A,A28,ES10.3)') "->", 'Fixed Epsilon value: ',self%gauss_epsil
        end select
        write(STD_OUT,'(30X,A,A28,L1)') "->", "Projection formulation: ", self % calculate_with_projection
        write(STD_OUT,'(30X,A,A28,L1)') "->", "Save disk average values: ", self % save_average
    end if
!
!   Read airfoil polar files (model B only), one table per disk per radial station
!   -----------------------------------------------------------------------------------
    if (self % model .eq. AD_MODEL_BLADE_ELEMENT) then
      do k = 1, self%num_disks
        do i = 1, self%disk_t(k)%n_radial
            arg=trim('./ActuatorDef/'//trim(self%disk_t(k)%airfoil_files(i,1)))
            OPEN( newunit = fid,file=trim(arg),status="old",action="read")
            READ(fid,'(A132)') char1
            READ(fid,*) self%disk_t(k)%airfoil_t(i)%num_aoa
            close(fid)

            associate (num_re => self%disk_t(k)%airfoil_t(i)%num_Re, num_aoa => self%disk_t(k)%airfoil_t(i)%num_aoa)

               allocate( self%disk_t(k)%airfoil_t(i)%aoa(num_aoa), self%disk_t(k)%airfoil_t(i)%cl(num_aoa,num_re), &
                         self%disk_t(k)%airfoil_t(i)%cd(num_aoa,num_re), self%disk_t(k)%airfoil_t(i)%Re(num_re) )

               do n_airfoil = 1, num_re

                  arg=trim('./ActuatorDef/'//trim(self%disk_t(k)%airfoil_files(i,n_airfoil)))
                  OPEN( newunit = fid,file=trim(arg),status="old",action="read")
                  READ(fid,'(A132)') char1

                  READ(fid,*) n_aoa
                  if ( n_aoa .ne. num_aoa ) then
                      print *, "Error: not same number of AoA in all files for same disk radial station, file: ", trim(arg)
                      call exit(99)
                  end if

                  READ(fid,'(A132)') char1
                  READ(fid,*) self%disk_t(k)%airfoil_t(i)%Re(n_airfoil)

                  READ(fid,'(A132)') char1

                  do ii = 1, num_aoa
                       READ(fid,*) self%disk_t(k)%airfoil_t(i)%aoa(ii), self%disk_t(k)%airfoil_t(i)%cl(ii,n_airfoil), &
                                   self%disk_t(k)%airfoil_t(i)%cd(ii,n_airfoil)
                       ! file is in deg, convert to rad
                       self%disk_t(k)%airfoil_t(i)%aoa(ii) = self%disk_t(k)%airfoil_t(i)%aoa(ii) * PI / 180.0_RP
                  end do

                  close(fid)
               end do ! number of airfoil files

            end associate
        end do ! number of radial stations
      end do ! number of disks
    end if
!
!   Build the (fixed, stationary) ring-point geometry
!   ----------------------------------------------------
    do k = 1, self%num_disks
      associate ( d => self%disk_t(k), nr => self%disk_t(k)%n_radial, na_ => self%disk_t(k)%n_azimuthal )
         allocate( d%theta(na_), d%point_xyz_loc(nr,na_,3) )
         dtheta = 2.0_RP*PI / real(na_,kind=RP)
         do j = 1, na_
            d%theta(j) = real(j-1,kind=RP) * dtheta
         end do
         do i = 1, nr
            do j = 1, na_
               d%point_xyz_loc(i,j,1) = d%hub_cood_x
               d%point_xyz_loc(i,j,2) = d%hub_cood_y + d%r_R(i)*cos(d%theta(j))
               d%point_xyz_loc(i,j,3) = d%hub_cood_z + d%r_R(i)*sin(d%theta(j))
            end do
         end do
!
!        Allocate per-point state arrays
!        --------------------------------
         allocate( d%eID(nr,na_), d%xi_loc(nr,na_,3), d%elementFound(nr,na_), &
                   d%Q_disk(nr,na_,NCONS), &
                   d%local_thrust(nr,na_), d%local_rotor_force(nr,na_), d%local_power(nr,na_), &
                   d%local_velocity(nr,na_), d%local_angle(nr,na_), d%local_aoa(nr,na_), d%local_Re(nr,na_) )
         d%local_thrust = 0.0_RP ; d%local_rotor_force = 0.0_RP ; d%local_power = 0.0_RP
         d%local_velocity = 0.0_RP ; d%local_angle = 0.0_RP ; d%local_aoa = 0.0_RP ; d%local_Re = 0.0_RP

         if (self % calculate_with_projection) then
            allocate( d%local_thrust_temp(nr,na_), d%local_rotor_force_temp(nr,na_), d%local_power_temp(nr,na_), &
                      d%local_velocity_temp(nr,na_), d%local_angle_temp(nr,na_), d%local_aoa_temp(nr,na_), &
                      d%local_Re_temp(nr,na_), d%local_gaussian_sum(nr,na_) )
         else
            if ( MPI_Process % doMPIAction ) then
               allocate( d%Q_disk_all(nr,na_,NCONS) )
            end if
         end if
      end associate
    end do
!
!   Get the solution file name
!   --------------------------
    solution_file = controlVariables % stringValueForKey( solutionFileNameKey, requestedLength = LINE_LENGTH )
    solution_file = trim(getFileName(solution_file))
    self % file_name = trim(solution_file)
!
!   Get the elements that are relevant (axis-aligned disk: axis = +x, plane = y-z)
!   ---------------------------------------------------------------------------------
    allocate(numElementsPerDiskAD(self%num_disks))
    numElementsPerDiskAD = 0
    element_loop:do eID = 1, mesh%no_of_elements
        do k=1, self%num_disks
            tolerance = self%tolerance_factor*self%disk_t(k)%radius
            r_square = minval(POW2(mesh%elements(eID)%geom%x(2,:,:,:)-self%disk_t(k)%hub_cood_y)) + minval(POW2(mesh%elements(eID)%geom%x(3,:,:,:)-self%disk_t(k)%hub_cood_z))
            if( r_square <= POW2(self%disk_t(k)%radius+tolerance) &
                .and. minval(mesh%elements(eID)%geom%x(1,:,:,:)) < self%disk_t(k)%hub_cood_x+tolerance &
                .and. maxval(mesh%elements(eID)%geom%x(1,:,:,:)) >self%disk_t(k)%hub_cood_x-tolerance) then
                numElementsPerDiskAD(k) = numElementsPerDiskAD(k) + 1
                cycle element_loop
            end if
        end do
    end do element_loop
    nelem = sum(numElementsPerDiskAD)
    allocate(elementsActuatedAD(nelem),diskOfElement(nelem))
    eIndex = 0
    element_loop2:do eID = 1, mesh%no_of_elements
        do k=1, self%num_disks
            tolerance = self%tolerance_factor*self%disk_t(k)%radius
            r_square = minval(POW2(mesh%elements(eID)%geom%x(2,:,:,:)-self%disk_t(k)%hub_cood_y)) + minval(POW2(mesh%elements(eID)%geom%x(3,:,:,:)-self%disk_t(k)%hub_cood_z))
            if( r_square <= POW2(self%disk_t(k)%radius+tolerance) &
                .and. minval(mesh%elements(eID)%geom%x(1,:,:,:)) < self%disk_t(k)%hub_cood_x+tolerance &
                .and. maxval(mesh%elements(eID)%geom%x(1,:,:,:)) >self%disk_t(k)%hub_cood_x-tolerance) then
                eIndex = eIndex + 1
                elementsActuatedAD(eIndex) = eID
                diskOfElement(eIndex) = k
                cycle element_loop2
            end if
        end do
    end do element_loop2
!
   select case (self % epsilon_type)
        case (0)
            do k = 1, self%num_disks
               self % disk_t(k) % gauss_epsil = self % gauss_epsil
            end do
        case (2)
            ! eps = k*delta; k is in gauss_epsil of farm; delta precalculated using element 1
            if (MPI_Process % doMPIAction) then
                if (nelem .gt. 0) then
                  delta_temp = (mesh % elements(elementsActuatedAD(1)) % geom % Volume / product(mesh % elements(elementsActuatedAD(1)) % Nxyz + 1)) ** (1.0_RP / 3.0_RP)
                  delta_count = 1
                else
                  delta_temp = 0.0_RP
                  delta_count = 0
                end if
#ifdef _HAS_MPI_
                call mpi_allreduce(delta_temp, delta, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD, ierr)
                call mpi_allreduce(delta_count, delta_paritions, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD, ierr)
#endif
                delta = delta / real(delta_paritions,kind=RP)
            else
                delta = (mesh % elements(elementsActuatedAD(1)) % geom % Volume / product(mesh % elements(elementsActuatedAD(1)) % Nxyz + 1)) ** (1.0_RP / 3.0_RP)
            end if
            do k = 1, self%num_disks
               self % disk_t(k) % gauss_epsil = self % gauss_epsil * delta
            end do
        case default
            do k = 1, self%num_disks
               self % disk_t(k) % gauss_epsil = self % gauss_epsil
            end do
    end select
!
!   One-time search for the element/reference coordinates of every ring point.
!   Points are stationary, so (unlike the AL model) this search is never repeated.
!   -----------------------------------------------------------------------------
    do k = 1, self%num_disks
      associate ( d => self%disk_t(k) )
      do j = 1, d%n_azimuthal
         do i = 1, d%n_radial
             x = d%point_xyz_loc(i,j,:)
             call FindActuatorPointElementAD(mesh, x, eID, xi, found)
             if ( (MPI_Process % doMPIAction) ) then
#ifdef _HAS_MPI_
                 call mpi_allreduce(found, allfound, 1, MPI_LOGICAL, MPI_LOR, MPI_COMM_WORLD, ierr)
#endif
             else
                 allfound = found
             end if

             if (allfound) then
                 d%eID(i,j) = eID
                 d%xi_loc(i,j,:) = xi
                 d%elementFound(i,j) = found
             else
                 print*, "Actuator disk point not found in mesh, x: ", x
                 print *, "i,j,k: ", i,j,k
                 call exit(99)
             end if
         end do
      end do
      end associate
    end do
!
    if (MPI_Process % isRoot) print*, "I allocate the AD device data"
!$acc update device(elementsActuatedAD)
!$acc update device(diskOfElement)
!$acc enter data copyin(self)
!$acc enter data copyin(self%disk_t)
    do k=1, self % num_disks
        !$acc enter data copyin(self%disk_t(k)%gauss_epsil)
        !$acc enter data copyin(self%disk_t(k)%theta)
        !$acc enter data copyin(self%disk_t(k)%point_xyz_loc)
        !$acc enter data copyin(self%disk_t(k)%local_thrust)
        !$acc enter data copyin(self%disk_t(k)%local_rotor_force)
        !$acc enter data copyin(self%disk_t(k)%local_power)
    end do
!
!   Create output files
!   -------------------
    if (MPI_Process % isRoot) then
      do k=1, self%num_disks
        write(file_id, '(I3.3)') k

        write(arg , '(A,A,A,A)') trim(self%file_name), "_Actuator_Disk_Forces_disk_", trim(file_id) , ".dat"
        open ( newunit = fID , file = trim(arg) , status = "unknown" , action = "write" )
        write(fid,'(4(2X,A24))') "time", "thrust", "torque", "power"
        close(fid)

        write(arg , '(A,A,A,A)') trim(self%file_name), "_Actuator_Disk_CP_CT_disk_", trim(file_id) , ".dat"
        open ( newunit = fID , file = trim(arg) , status = "unknown" , action = "write" )
        write(fid,'(3(2X,A24))') "time", "Cp(power_coef)", "Ct(thust_coef)"
        close(fid)

        if (self % save_average) then
            write(arg , '(A,A,A,A)') trim(self%file_name), "_Actuator_Disk_average_disk_", trim(file_id) , ".dat"
            open ( newunit = fID , file = trim(arg) , status = "unknown" , action = "write" )
            write(fid,'(6(2X,A24))') "R", "U", "AoA", "Re", "Tangential_Force", "Axial_Force"
            close(fid)
        end if
      end do
    end if
!
       self % number_iterations = 0
       self % active = .true.
!
   end subroutine ConstructDiskFarm
!
!///////////////////////////////////////////////////////////////////////////////////////
!
   subroutine DestructDiskFarm(self)
   implicit none
   type(DiskFarm_t), intent(inout)       :: self

   safedeallocate(self%disk_t)

   end subroutine DestructDiskFarm
!
!///////////////////////////////////////////////////////////////////////////////////////
!
   subroutine UpdateDiskFarm(self,time, mesh)
   use fluiddata
   use HexMeshClass
   use PhysicsStorage
   use MPI_Process_Info
   implicit none

   type(DiskFarm_t), intent(inout)   :: self
   real(kind=RP), intent(in)         :: time
   type(HexMesh), intent(in)         :: mesh

   !local variables
   integer                           :: ii, jj, kk, ierr

   if (.not. self % active) return

   self % time = time

   if (self % calculate_with_projection) then
!
!      Forces are (re-)computed lazily, node by node, inside ForcesDiskFarm; here we only
!      reset the accumulators used for the (Gaussian-weighted) diagnostic output.
!      -------------------------------------------------------------------------------------
       do kk = 1, self%num_disks
          self%disk_t(kk)%local_thrust_temp(:,:) = 0.0_RP
          self%disk_t(kk)%local_rotor_force_temp(:,:) = 0.0_RP
          self%disk_t(kk)%local_power_temp(:,:) = 0.0_RP
          self%disk_t(kk)%local_velocity_temp(:,:) = 0.0_RP
          self%disk_t(kk)%local_angle_temp(:,:) = 0.0_RP
          self%disk_t(kk)%local_aoa_temp(:,:) = 0.0_RP
          self%disk_t(kk)%local_Re_temp(:,:) = 0.0_RP
          self%disk_t(kk)%local_gaussian_sum(:,:) = 0.0_RP
       end do

   else
!
!      Non-projection: sample Q once per (stationary) ring point and get its local force.
!      Since the points never move, eID/xi found in ConstructDiskFarm remain valid.
!      -------------------------------------------------------------------------------------
       do kk = 1, self%num_disks
         associate ( d => self%disk_t(kk) )
         do jj = 1, d%n_azimuthal
            do ii = 1, d%n_radial
               if (d%elementFound(ii,jj)) then
                  d%Q_disk(ii,jj,:) = interpolateQAD(mesh, d%eID(ii,jj), d%xi_loc(ii,jj,:))
               else
                  d%Q_disk(ii,jj,:) = 0.0_RP
               end if
            end do
         end do
         end associate
       end do

       if ( (MPI_Process % doMPIAction) ) then
         do kk = 1, self%num_disks
           associate ( d => self%disk_t(kk) )
#ifdef _HAS_MPI_
           call mpi_allreduce(d%Q_disk(:,:,:), d%Q_disk_all(:,:,:), NCONS*d%n_radial*d%n_azimuthal, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD, ierr)
#endif
           d%Q_disk = d%Q_disk_all
           end associate
         end do
       end if

       do kk = 1, self%num_disks
         associate ( d => self%disk_t(kk) )
         do jj = 1, d%n_azimuthal
            do ii = 1, d%n_radial
               call DiskGetLocalForces(self, ii, jj, kk, d%Q_disk(ii,jj,:), 1.0_RP, d%local_angle(ii,jj), &
                                        d%local_aoa(ii,jj), d%local_velocity(ii,jj), d%local_Re(ii,jj), d%local_thrust(ii,jj), &
                                        d%local_rotor_force(ii,jj), d%local_power(ii,jj))
            end do
         end do
         end associate
       end do
!
       do kk=1, self % num_disks
           !$acc update device(self%disk_t(kk)%local_thrust)
           !$acc update device(self%disk_t(kk)%local_rotor_force)
           !$acc update device(self%disk_t(kk)%local_power)
       end do

   end if
!$acc wait
!
   end subroutine UpdateDiskFarm
!
!///////////////////////////////////////////////////////////////////////////////////////
!
   subroutine ForcesDiskFarm(self, mesh, time)
   use PhysicsStorage
   use HexMeshClass
   use fluiddata
   implicit none

   type(DiskFarm_t) , intent(inout)  :: self
   type(HexMesh), intent(in)         :: mesh
   real(kind=RP),intent(in)          :: time

! local vars
   real(kind=RP)                     :: Non_dimensional, Non_dimensional_energy, t, interp
   integer                           :: ii,jj, kk
   integer                           :: i,j, k
   integer                           :: eID, eIndex
   real(kind=RP), dimension(NDIM)    :: actuator_source
   real(kind=RP)                     :: energy_source
   real(kind=RP)                     :: actuator_source_x_iijj
   real(kind=RP)                     :: actuator_source_y_iijj
   real(kind=RP)                     :: actuator_source_z_iijj
   real(kind=RP)                     :: energy_source_iijj
   real(kind=RP), dimension(NDIM)    :: xlocal
   real(kind=RP)                     :: local_angle, local_aoa, local_velocity, local_Re, local_thrust, local_rotor_force, local_power
   real(kind=RP)                     :: theta_p, epsil
#if defined(MULTIPHASE)
   real(kind=RP)                     :: invSqrtRho
#endif

    if (.not. self % active) return

    Non_dimensional = POW2(refValues % V) * refValues % rho / Lref
    Non_dimensional_energy = Non_dimensional * refValues % V
    t = time * Lref / refValues % V

    if (self % calculate_with_projection) then
!
!       ---------------------------------------------------------------------------
!       Projection: loop over every mesh quadrature node in every actuated element,
!       re-evaluate the local force law using that node's own Q (CPU only, as in AL)
!       ---------------------------------------------------------------------------
!
        do eIndex = 1, size(elementsActuatedAD)
            eID = elementsActuatedAD(eIndex)
            kk = diskOfElement(eIndex)

            do k = 0, mesh % elements(eID) % Nxyz(3)   ; do j = 0, mesh % elements(eID) % Nxyz(2) ; do i = 0, mesh % elements(eID) % Nxyz(1)
                actuator_source(:) = 0.0_RP
                energy_source = 0.0_RP

                associate ( d => self % disk_t(kk) )
                do jj = 1, d % n_azimuthal
                    do ii = 1, d % n_radial
                        interp = GaussianInterpolationAD( mesh % elements(eID) % geom % x(:,i,j,k), d%point_xyz_loc(ii,jj,:), d%gauss_epsil)
                        call DiskGetLocalForces(self, ii, jj, kk, mesh%elements(eID)%storage%Q(:,i,j,k), interp, local_angle, local_aoa, local_velocity, local_Re, local_thrust, local_rotor_force, local_power)

                        theta_p = d % theta(jj)
!
!                       local_thrust/local_rotor_force are already "force on the fluid"
!                       ------------------------------------------------------------------
                        actuator_source(1) = actuator_source(1) + local_thrust
                        actuator_source(2) = actuator_source(2) - local_rotor_force*sin(theta_p)
                        actuator_source(3) = actuator_source(3) + local_rotor_force*cos(theta_p)
                        energy_source = energy_source + local_power

                        d%local_thrust_temp(ii,jj)      = d%local_thrust_temp(ii,jj)      + local_thrust
                        d%local_rotor_force_temp(ii,jj) = d%local_rotor_force_temp(ii,jj) + local_rotor_force
                        d%local_power_temp(ii,jj)       = d%local_power_temp(ii,jj)       + local_power
                        d%local_velocity_temp(ii,jj)    = d%local_velocity_temp(ii,jj)    + local_velocity*interp
                        d%local_angle_temp(ii,jj)       = d%local_angle_temp(ii,jj)       + local_angle*interp
                        d%local_aoa_temp(ii,jj)         = d%local_aoa_temp(ii,jj)         + local_aoa*interp
                        d%local_Re_temp(ii,jj)          = d%local_Re_temp(ii,jj)          + local_Re*interp
                        d%local_gaussian_sum(ii,jj)     = d%local_gaussian_sum(ii,jj)     + interp
                    end do
                end do
                end associate

                actuator_source = actuator_source / Non_dimensional
                energy_source = energy_source / Non_dimensional_energy

#if defined(NAVIERSTOKES)
                mesh % elements(eID) % storage % S_NS(IRHOU,i,j,k) = actuator_source(1)
                mesh % elements(eID) % storage % S_NS(IRHOV,i,j,k) = actuator_source(2)
                mesh % elements(eID) % storage % S_NS(IRHOW,i,j,k) = actuator_source(3)
                mesh % elements(eID) % storage % S_NS(IRHOE,i,j,k) = energy_source
#elif defined(INCNS)
                mesh % elements(eID) % storage % S_NS(INSRHOU,i,j,k) = actuator_source(1)
                mesh % elements(eID) % storage % S_NS(INSRHOV,i,j,k) = actuator_source(2)
                mesh % elements(eID) % storage % S_NS(INSRHOW,i,j,k) = actuator_source(3)
#elif defined(MULTIPHASE)
                invSqrtRho = 1.0_RP / sqrt(mesh % elements(eID) % storage % rho(i,j,k))
                mesh % elements(eID) % storage % S_NS(IMSQRHOU,i,j,k) = actuator_source(1)*invSqrtRho
                mesh % elements(eID) % storage % S_NS(IMSQRHOV,i,j,k) = actuator_source(2)*invSqrtRho
                mesh % elements(eID) % storage % S_NS(IMSQRHOW,i,j,k) = actuator_source(3)*invSqrtRho
#endif
            end do                  ; end do                ; end do
        end do

    else ! no projection
!
!       ---------------------------------------------------------------------------------
!       Non-projection: local_thrust/local_rotor_force/local_power were already computed
!       once per point in UpdateDiskFarm; here we only deposit them with the Gaussian
!       kernel. As in the AL model, the collapsed loop below assumes a uniform, isotropic
!       polynomial order across the actuated elements (mesh % Nx(1) is used for every
!       direction/element so the OpenACC collapse clause stays rectangular).
!       ---------------------------------------------------------------------------------
!
        !$acc parallel loop gang collapse(4) present(self,mesh) private(xlocal)
        do eIndex = 1, size(elementsActuatedAD) ;  do k = 0, mesh % Nx(1) ;  do j = 0, mesh % Nx(1) ; do i = 0, mesh % Nx(1)

            eID = elementsActuatedAD(eIndex)
            kk = diskOfElement(eIndex)

                actuator_source_x_iijj = 0.0_RP
                actuator_source_y_iijj = 0.0_RP
                actuator_source_z_iijj = 0.0_RP
                energy_source_iijj = 0.0_RP

                xlocal = mesh % elements(eID) % geom % x(:,i,j,k)
                epsil = self%disk_t(kk) % gauss_epsil
                !$acc loop vector collapse(2) reduction(+:actuator_source_x_iijj, actuator_source_y_iijj, actuator_source_z_iijj, energy_source_iijj)
                do jj = 1, self % disk_t(kk) % n_azimuthal ; do ii = 1, self % disk_t(kk) % n_radial

                    theta_p = self%disk_t(kk)%theta(jj)
                    local_rotor_force = self%disk_t(kk)%local_rotor_force(ii,jj)
                    local_thrust = self%disk_t(kk)%local_thrust(ii,jj)
                    local_power = self%disk_t(kk)%local_power(ii,jj)

                    interp = exp( -(POW2(xlocal(1) - self%disk_t(kk)%point_xyz_loc(ii,jj,1)) &
                                  + POW2(xlocal(2) - self%disk_t(kk)%point_xyz_loc(ii,jj,2)) &
                                  + POW2(xlocal(3) - self%disk_t(kk)%point_xyz_loc(ii,jj,3))) &
                                  / POW2(epsil) ) / ( POW3(epsil) * pi**(3.0_RP/2.0_RP) )

                    actuator_source_x_iijj = actuator_source_x_iijj + local_thrust * interp
                    actuator_source_y_iijj = actuator_source_y_iijj - local_rotor_force*sin(theta_p) * interp
                    actuator_source_z_iijj = actuator_source_z_iijj + local_rotor_force*cos(theta_p) * interp
                    energy_source_iijj     = energy_source_iijj     + local_power * interp
                end do                  ; end do

                actuator_source_x_iijj = actuator_source_x_iijj / Non_dimensional
                actuator_source_y_iijj = actuator_source_y_iijj / Non_dimensional
                actuator_source_z_iijj = actuator_source_z_iijj / Non_dimensional
                energy_source_iijj     = energy_source_iijj     / Non_dimensional_energy

#if defined(NAVIERSTOKES)
                mesh % elements(eID) % storage % S_NS(IRHOU,i,j,k) = actuator_source_x_iijj
                mesh % elements(eID) % storage % S_NS(IRHOV,i,j,k) = actuator_source_y_iijj
                mesh % elements(eID) % storage % S_NS(IRHOW,i,j,k) = actuator_source_z_iijj
                mesh % elements(eID) % storage % S_NS(IRHOE,i,j,k) = energy_source_iijj
#elif defined(INCNS)
                mesh % elements(eID) % storage % S_NS(INSRHOU,i,j,k) = actuator_source_x_iijj
                mesh % elements(eID) % storage % S_NS(INSRHOV,i,j,k) = actuator_source_y_iijj
                mesh % elements(eID) % storage % S_NS(INSRHOW,i,j,k) = actuator_source_z_iijj
#elif defined(MULTIPHASE)
                invSqrtRho = 1.0_RP / sqrt(mesh % elements(eID) % storage % rho(i,j,k))
                mesh % elements(eID) % storage % S_NS(IMSQRHOU,i,j,k) = actuator_source_x_iijj*invSqrtRho
                mesh % elements(eID) % storage % S_NS(IMSQRHOV,i,j,k) = actuator_source_y_iijj*invSqrtRho
                mesh % elements(eID) % storage % S_NS(IMSQRHOW,i,j,k) = actuator_source_z_iijj*invSqrtRho
#endif
            end do                  ; end do                ; end do
        end do
!$acc end parallel loop
    endif
!
   end subroutine  ForcesDiskFarm
!
!///////////////////////////////////////////////////////////////////////////////////////
!
   subroutine WriteDiskFarmForces(self,time,iter,last)
   use fluiddata
   use PhysicsStorage
   use MPI_Process_Info
   implicit none

   type(DiskFarm_t), intent(inout)  :: self
   real(kind=RP),intent(in)      :: time
   integer, intent(in)           :: iter
   logical, optional             :: last
   integer                       :: fid
   CHARACTER(LEN=LINE_LENGTH)    :: arg
   real(kind=RP)                 :: t
   integer                       :: ii, jj, kk, ierr
   logical                       :: isLast
   logical                       :: save_instant
   real(kind=RP), dimension(:,:), allocatable :: temp2d
   CHARACTER(LEN=5)           :: file_id

   if (.not. self % active) return

   if (present(last)) then
       isLast = last
   else
       isLast = .false.
   end if

   if (isLast) then
      if ( .not. self % save_average ) return
      if ( .not. MPI_Process % isRoot ) return
      do kk=1, self%num_disks
          write(file_id, '(I3.3)') kk
          write(arg , '(A,A,A,A)') trim(self%file_name), "_Actuator_Disk_average_disk_", trim(file_id) , ".dat"
          open( newunit = fID , file = trim(arg) , action = "write" , access = "append" , status = "old" )
          do ii = 1, self % disk_t(kk) % n_radial
            write(fid,"(6(2X,ES24.16))") self%disk_t(kk)%r_R(ii), self%disk_t(kk)%average_conditions(ii,:)
          end do
          close(fid)
      end do
      return
    end if

   save_instant = self%save_instant .and. ( mod(iter,self % save_iterations) .eq. 0 )
   t = time * Lref / refValues%V
!
!  In projection mode, normalize the Gaussian-weighted accumulators before computing
!  the integrated (whole-disk) diagnostics, exactly as the AL model does.
!  ---------------------------------------------------------------------------------------
   if (self%calculate_with_projection) then

     if ( (MPI_Process % doMPIAction) ) then
       do kk = 1, self%num_disks
         associate ( d => self%disk_t(kk) )
           allocate( temp2d(d%n_radial,d%n_azimuthal) )
#ifdef _HAS_MPI_
           temp2d = d%local_thrust_temp
           call mpi_allreduce(temp2d, d%local_thrust_temp, d%n_radial*d%n_azimuthal, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD, ierr)
           temp2d = d%local_rotor_force_temp
           call mpi_allreduce(temp2d, d%local_rotor_force_temp, d%n_radial*d%n_azimuthal, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD, ierr)
           temp2d = d%local_power_temp
           call mpi_allreduce(temp2d, d%local_power_temp, d%n_radial*d%n_azimuthal, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD, ierr)
           temp2d = d%local_velocity_temp
           call mpi_allreduce(temp2d, d%local_velocity_temp, d%n_radial*d%n_azimuthal, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD, ierr)
           temp2d = d%local_angle_temp
           call mpi_allreduce(temp2d, d%local_angle_temp, d%n_radial*d%n_azimuthal, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD, ierr)
           temp2d = d%local_aoa_temp
           call mpi_allreduce(temp2d, d%local_aoa_temp, d%n_radial*d%n_azimuthal, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD, ierr)
           temp2d = d%local_Re_temp
           call mpi_allreduce(temp2d, d%local_Re_temp, d%n_radial*d%n_azimuthal, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD, ierr)
           temp2d = d%local_gaussian_sum
           call mpi_allreduce(temp2d, d%local_gaussian_sum, d%n_radial*d%n_azimuthal, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD, ierr)
#endif
           deallocate(temp2d)
         end associate
       end do
     end if

     do kk = 1, self%num_disks
        associate ( d => self%disk_t(kk) )
        do jj = 1, d%n_azimuthal
           do ii = 1, d%n_radial
              if (d%local_gaussian_sum(ii,jj) .gt. 0.0_RP) then
                 d%local_thrust(ii,jj)      = d%local_thrust_temp(ii,jj)      / d%local_gaussian_sum(ii,jj)
                 d%local_rotor_force(ii,jj) = d%local_rotor_force_temp(ii,jj) / d%local_gaussian_sum(ii,jj)
                 d%local_power(ii,jj)       = d%local_power_temp(ii,jj)       / d%local_gaussian_sum(ii,jj)
                 d%local_velocity(ii,jj)    = d%local_velocity_temp(ii,jj)    / d%local_gaussian_sum(ii,jj)
                 d%local_angle(ii,jj)       = d%local_angle_temp(ii,jj)       / d%local_gaussian_sum(ii,jj)
                 d%local_aoa(ii,jj)         = d%local_aoa_temp(ii,jj)         / d%local_gaussian_sum(ii,jj)
                 d%local_Re(ii,jj)          = d%local_Re_temp(ii,jj)          / d%local_gaussian_sum(ii,jj)
              end if
           end do
        end do
        end associate
     end do

   end if

   if ( .not. MPI_Process % isRoot ) return
!
!  Save in memory the time step integrated (whole-disk) forces
!  -----------------------------------------------------------
   call DiskUpdateIntegratedForces(self)

   do kk=1, self%num_disks
     write(file_id, '(I3.3)') kk
     write(arg , '(A,A,A,A)') trim(self%file_name), "_Actuator_Disk_Forces_disk_", trim(file_id) , ".dat"

     open( newunit = fID , file = trim(arg) , action = "write" , access = "append" , status = "old" )
     write(fid,"(4(2X,ES24.16))") t, self%disk_t(kk)%total_thrust, self%disk_t(kk)%total_torque, self%disk_t(kk)%total_power
     close(fid)

     write(arg , '(A,A,A,A)') trim(self%file_name), "_Actuator_Disk_CP_CT_disk_", trim(file_id) , ".dat"
     open( newunit = fID , file = trim(arg) , action = "write" , access = "append" , status = "old" )
     write(fid,"(3(2X,ES24.16))") t, self%disk_t(kk)%Cp, self%disk_t(kk)%Ct
     close(fid)

     if (save_instant) then
       write(arg , '(2A,I10.10,A,I3.3,A)') trim(self%file_name), "_Actuator_Disk_instant_",iter, "_disk_", kk, ".dat"
       open ( newunit = fID , file = trim(arg) , status = "unknown" , action = "write" )
       write(fid,'(6(2X,A24))') "R", "U", "AoA", "Re", "Tangential_Force", "Axial_Force"
       do ii = 1, self % disk_t(kk) % n_radial
         write(fid,"(6(2X,ES24.16))") self%disk_t(kk)%r_R(ii), &
              self%disk_t(kk)%local_velocity(ii,1), &
              self%disk_t(kk)%local_aoa(ii,1) * 180.0_RP / PI, &
              self%disk_t(kk)%local_Re(ii,1), &
              sum(self%disk_t(kk)%local_rotor_force(ii,:)), &
              sum(self%disk_t(kk)%local_thrust(ii,:))
       end do
       close(fid)
     end if

     if (self % save_average) then
        associate ( d => self%disk_t(kk) )
        if (.not. allocated(d%average_conditions)) then
           allocate( d%average_conditions(d%n_radial,5) )
           d%average_conditions = 0.0_RP
        end if
        d%average_conditions = d%average_conditions * real(self % number_iterations,RP)
        do ii = 1, d%n_radial
           d%average_conditions(ii,1) = d%average_conditions(ii,1) + d%local_velocity(ii,1)
           d%average_conditions(ii,2) = d%average_conditions(ii,2) + d%local_aoa(ii,1) * 180.0_RP / PI
           d%average_conditions(ii,3) = d%average_conditions(ii,3) + d%local_Re(ii,1)
           d%average_conditions(ii,4) = d%average_conditions(ii,4) + sum(d%local_rotor_force(ii,:))
           d%average_conditions(ii,5) = d%average_conditions(ii,5) + sum(d%local_thrust(ii,:))
        end do
        end associate
     end if
   end do

   if (self % save_average) self % number_iterations = self % number_iterations + 1
   if (self % save_average) then
      do kk = 1, self%num_disks
         self % disk_t(kk) % average_conditions = self % disk_t(kk) % average_conditions / real(self % number_iterations,RP)
      end do
   end if

end subroutine WriteDiskFarmForces
!
!///////////////////////////////////////////////////////////////////////////////////////
!
    Subroutine DiskGetLocalForces(self, ii, jj, kk, Q, interp, local_angle, local_aoa, local_velocity, local_Re, local_thrust, local_rotor_force, local_power)
        use PhysicsStorage
        use fluiddata
#if defined(NAVIERSTOKES)
        use VariableConversion, only: Temperature, SutherlandsLaw, Pressure
#endif
        implicit none
        type(DiskFarm_t)                              :: self
        integer, intent(in)                           :: ii, jj, kk
        real(kind=RP), dimension(NCONS), intent(in)   :: Q
        real(kind=RP), intent(in)                     :: interp
        real(kind=RP), intent(out)                    :: local_angle
        real(kind=RP), intent(out)                    :: local_aoa           ! angle of attack [rad] (0 for Model A)
        real(kind=RP), intent(out)                    :: local_velocity
        real(kind=RP), intent(out)                    :: local_Re
        real(kind=RP), intent(out)                    :: local_thrust        ! axial force ON THE FLUID [N]
        real(kind=RP), intent(out)                    :: local_rotor_force   ! tangential force ON THE FLUID [N]
        real(kind=RP), intent(out)                    :: local_power         ! power delivered TO THE FLUID [W]

        !local variables
        real(kind=RP)                                 :: density, Cl, Cd, aoa, g1_func, tip_correct
        real(kind=RP)                                 :: wind_speed_axial, wind_speed_rot
        real(kind=RP)                                 :: lift_force, drag_force
        real(kind=RP)                                 :: T, muL, theta_p, r, dr, dtheta, scale_B
#if defined(NAVIERSTOKES)
        real(kind=RP)                                 :: soundSpeed, localMach, betaPG
#endif
#if defined(MULTIPHASE)
        real(kind=RP)                                 :: rho, invSqrtRho
#endif

        theta_p = self % disk_t(kk) % theta(jj)
!
!       -----------------------------
!       get flow related variables
!       -----------------------------
!
#if defined(NAVIERSTOKES)
        density = Q(IRHO) * refValues % rho
        wind_speed_axial = (Q(IRHOU)/Q(IRHO)) * refValues % V
        wind_speed_rot = ( -Q(IRHOV)*sin(theta_p) + Q(IRHOW)*cos(theta_p) ) / Q(IRHO) * refValues % V

#elif defined(INCNS)
        density = Q(INSRHO) * refValues % rho
        wind_speed_axial = (Q(INSRHOU)/Q(INSRHO)) * refValues % V
        wind_speed_rot = ( -Q(INSRHOV)*sin(theta_p) + Q(INSRHOW)*cos(theta_p) ) / Q(INSRHO) * refValues % V
#elif defined(MULTIPHASE)
        rho = dimensionless % rho(2) + (dimensionless % rho(1)-dimensionless % rho(2)) * Q(IMC)
        rho = min(max(rho, dimensionless % rho_min),dimensionless % rho_max)
        invSqrtRho = 1.0_RP/sqrt(rho)
        density = rho * refValues % rho

        wind_speed_axial = (Q(IMSQRHOU)*invSqrtRho) * refValues % V
        wind_speed_rot = ( -Q(IMSQRHOV)*sin(theta_p) + Q(IMSQRHOW)*cos(theta_p) )*invSqrtRho * refValues % V
#endif

        r  = self % disk_t(kk) % r_R(ii)
        dr = self % disk_t(kk) % r_R(ii) - self % disk_t(kk) % r_R(ii-1)
        dtheta = 2.0_RP*PI / real(self % disk_t(kk) % n_azimuthal, kind=RP)

        select case (self % model)
        case (AD_MODEL_PRESSURE)
!
!          -----------------------------------------------
!          Model A: prescribed pressure/enthalpy rise, Eq. (18)
!          -----------------------------------------------
!
           local_thrust = self % disk_t(kk) % delta_p0(ii) * r * dr * dtheta * interp

           if (self % disk_t(kk) % use_swirl) then
              local_rotor_force = self % disk_t(kk) % swirl_q(ii) * r * dr * dtheta * interp
           else
              local_rotor_force = 0.0_RP
           end if

           local_velocity = sqrt(POW2(wind_speed_axial) + POW2(wind_speed_rot))
           local_angle    = 0.0_RP
           local_aoa      = 0.0_RP
           local_Re       = 0.0_RP

        case default ! AD_MODEL_BLADE_ELEMENT
!
!          --------------------------------------------------------------------------
!          Model B: local blade-element (Cl-Cd), Eqs. (7)-(13), smeared via Eq. (19)
!          --------------------------------------------------------------------------
!
           tip_correct = 1.0_RP
           aoa = 0.0_RP

#if defined(NAVIERSTOKES)
           T   = Temperature(Q)
           muL = SutherlandsLaw(T) * refValues % mu
#else
           muL = refValues % mu
#endif

           local_velocity = sqrt( POW2(self % disk_t(kk) % rot_speed * r - wind_speed_rot) + POW2(wind_speed_axial) )
           local_angle    = atan( wind_speed_axial / (self % disk_t(kk) % rot_speed * r - wind_speed_rot) )

           aoa      = local_angle - (self % disk_t(kk) % twist(ii) + self % disk_t(kk) % blade_pitch)
           local_aoa = aoa
           local_Re = local_velocity * self % disk_t(kk) % chord(ii) * density / muL

           call Get_Cl_Cd_from_airfoil_data_AD(self % disk_t(kk) % airfoil_t(ii), aoa, local_Re, Cl, Cd)

#if defined(NAVIERSTOKES)
           if (self % disk_t(kk) % compressibility_correction) then
              soundSpeed = refValues % V * sqrt( thermodynamics % gamma * Pressure(Q) / Q(IRHO) )
              localMach  = local_velocity / soundSpeed
              betaPG     = sqrt(max(1.0_RP - POW2(min(localMach,0.95_RP)), 0.0975_RP))
              Cl = Cl / betaPG
           end if
#endif
!
!          tip-loss correction, Eqs. (12)-(13)
!          --------------------------------------
           g1_func = exp( -self % disk_t(kk) % tip_c1 * ( real(self % disk_t(kk) % num_blades,kind=RP) * self % disk_t(kk) % rot_speed * self % disk_t(kk) % radius / refValues % V &
                          - self % disk_t(kk) % tip_c2) ) + 0.1_RP

           tip_correct = 2.0_RP/PI * ( acos( exp( -g1_func * real(self % disk_t(kk) % num_blades,kind=RP) * (self % disk_t(kk) % radius - r) / &
                                       abs(2.0_RP * r * sin(local_angle)) ) ) )

           scale_B = real(self % disk_t(kk) % num_blades, kind=RP) / real(self % disk_t(kk) % n_azimuthal, kind=RP)

           lift_force = 0.5_RP * density * Cl * tip_correct * POW2(local_velocity) * self % disk_t(kk) % chord(ii) * dr * interp * scale_B
           drag_force = 0.5_RP * density * Cd * tip_correct * POW2(local_velocity) * self % disk_t(kk) % chord(ii) * dr * interp * scale_B
!
!          action-reaction: lift_force/drag_force above act ON the virtual blade element;
!          the force delivered TO THE FLUID is the opposite, exactly as in the AL model.
!          ----------------------------------------------------------------------------------
           local_thrust      = -( lift_force * cos(local_angle) + drag_force * sin(local_angle) )
           local_rotor_force = -( lift_force * sin(local_angle) - drag_force * cos(local_angle) )

        end select
!
!       Power delivered to the fluid at this point, Eq. (16)
!       -------------------------------------------------------
        local_power = local_thrust * wind_speed_axial + local_rotor_force * wind_speed_rot

    End Subroutine DiskGetLocalForces
!
!///////////////////////////////////////////////////////////////////////////////////////
!
    Subroutine DiskUpdateIntegratedForces(self)
        use fluiddata
        Implicit None

        type(DiskFarm_t), intent(inout)      :: self
        !local variables
        integer                           :: kk

    do kk = 1, self%num_disks
      associate ( d => self%disk_t(kk) )
      d%total_thrust = sum(d%local_thrust)
      d%total_torque = sum(d%local_rotor_force * spread(d%r_R(1:d%n_radial),2,d%n_azimuthal))
      d%total_power  = sum(d%local_power)

      d%Ct = 2.0_RP * d%total_thrust / (refValues%rho * POW2(refValues%V) * pi * POW2(d%radius))
      d%Cp = 2.0_RP * d%total_power  / (refValues%rho * POW3(refValues%V) * pi * POW2(d%radius))
      end associate
    end do
!
    End Subroutine DiskUpdateIntegratedForces
!
!///////////////////////////////////////////////////////////////////////////////////////
!
    Function GaussianInterpolationAD(x, x_point, gauss_epsil)
        implicit none
        !$acc routine seq
        real(kind=RP), intent(in)               :: x(NDIM)
        real(kind=RP), intent(in)               :: x_point(NDIM)
        real(kind=RP), intent(in)               :: gauss_epsil
        real(kind=RP)                           :: GaussianInterpolationAD

        GaussianInterpolationAD = exp( -(POW2(x(1) - x_point(1)) + POW2(x(2) - x_point(2)) + POW2(x(3) - x_point(3))) / POW2(gauss_epsil) ) / ( POW3(gauss_epsil) * pi**(3.0_RP/2.0_RP) )

    End Function GaussianInterpolationAD
!
!///////////////////////////////////////////////////////////////////////////////////////
!
! based on HexMesh_FindPointWithCoords, without curvature only in the precalculated list
    Subroutine FindActuatorPointElementAD(mesh, x, eID, xi, success)
       use HexMeshClass
       Implicit None

       type(HexMesh), intent(in)                     :: mesh
       real(kind=RP), dimension(NDIM), intent(in)    :: x       ! physical space
       integer, intent(out)                          :: eID
       real(kind=RP), dimension(NDIM), intent(out)   :: xi      ! computational space
       logical, intent(out)                          :: success
       !
       logical                                       :: found
       integer                                       :: eIndex

       success = .false.
       found = .false.
!
!      Search in linear (not curved) mesh (faster and safer); AD points are expected linear
!      -----------------------------------------------------------------------------------------
       do eIndex = 1, size(elementsActuatedAD)
          eID = elementsActuatedAD(eIndex)
          found = mesh % elements(eID) % FindPointInLinElement(x, mesh % nodes)
          if ( found ) exit
       end do

       if (found) then
           success = mesh % elements(eID) % FindPointWithCoords(x, mesh % dir2D_ctrl, xi)
       else
          eID = 0
       end if
!
    End Subroutine FindActuatorPointElementAD
!
!///////////////////////////////////////////////////////////////////////////////////////
!
    subroutine Get_Cl_Cd_from_airfoil_data_AD(airfoil, aoa, Re, Cl_out, Cd_out)
         implicit none

         type (airfoil_disk_t), intent(in)   :: airfoil
         real(KIND=RP), intent(in)      :: aoa, Re
         real(KIND=RP), intent(out)     :: Cl_out, Cd_out
         integer                        :: i,k
         real(kind=RP), dimension(2)    :: Cl_inter, Cd_inter

         Cl_out=0.0_RP
         Cd_out=0.0_RP

         if (airfoil%num_Re .eq. 1) then
             do i=1, airfoil % num_aoa-1
                 if (airfoil%aoa(i+1)>=aoa .and. airfoil%aoa(i)<=aoa ) then
                    Cl_out=InterpolateAirfoilDataAD(airfoil%aoa(i),airfoil%aoa(i+1),airfoil%cl(i,1),airfoil%cl(i+1,1),aoa)
                    Cd_out=InterpolateAirfoilDataAD(airfoil%aoa(i),airfoil%aoa(i+1),airfoil%cd(i,1),airfoil%cd(i+1,1),aoa)
                    exit
                 endif
             end do

         else

             do k=1, airfoil % num_Re-1
                 if (airfoil%Re(k+1)>=Re .and. airfoil%Re(k)<=Re ) then
                     do i=1, airfoil % num_aoa-1
                         if (airfoil%aoa(i+1)>=aoa .and. airfoil%aoa(i)<=aoa ) then

                            Cl_inter(1) = InterpolateAirfoilDataAD(airfoil%aoa(i),airfoil%aoa(i+1),airfoil%cl(i,k),airfoil%cl(i+1,k),aoa)
                            Cl_inter(2) = InterpolateAirfoilDataAD(airfoil%aoa(i),airfoil%aoa(i+1),airfoil%cl(i,k+1),airfoil%cl(i+1,k+1),aoa)
                            Cl_out=InterpolateAirfoilDataAD(airfoil%Re(k),airfoil%Re(k+1),Cl_inter(1),Cl_inter(2),Re)

                            Cd_inter(1) = InterpolateAirfoilDataAD(airfoil%aoa(i),airfoil%aoa(i+1),airfoil%cd(i,k),airfoil%cd(i+1,k),aoa)
                            Cd_inter(2) = InterpolateAirfoilDataAD(airfoil%aoa(i),airfoil%aoa(i+1),airfoil%cd(i,k+1),airfoil%cd(i+1,k+1),aoa)
                            Cd_out=InterpolateAirfoilDataAD(airfoil%Re(k),airfoil%Re(k+1),Cd_inter(1),Cd_inter(2),Re)

                            exit
                         endif
                     end do
                 endif
             end do

         end if

    end subroutine Get_Cl_Cd_from_airfoil_data_AD

! linear interpolation given two points; returns y for new_x following line coefs (a,b) with y=ax+b
function InterpolateAirfoilDataAD(x1,x2,y1,y2,new_x)
   implicit none

   real(KIND=RP), intent(in)    :: x1, x2, y1, y2, new_x
   real(KIND=RP)                :: a, b, InterpolateAirfoilDataAD

    if(abs(x1-x2)<1.0e-6_RP) then
      a=100.0_RP
   else
      a=(y1- y2)/(x1- x2)
   endif
    b= y1-a*x1;
    InterpolateAirfoilDataAD=a*new_x+b
end function

! high order interpolation of Q
Function interpolateQAD(mesh,eID,xi) result(Qe)
   use HexMeshClass
   use PhysicsStorage
   use NodalStorageClass

   Implicit None
   type(HexMesh), intent(in)    :: mesh
   integer, intent(in)          :: eID
   real(kind=RP), dimension(NDIM), intent(in) :: xi
   real(kind=RP), dimension(NCONS)   :: Qe

   integer                        :: k, j, i
   integer, dimension(NDIM)       :: Nxyz
   type(NodalStorage_t), pointer  :: spAxi, spAeta, spAzeta
   real(kind=RP), allocatable     :: lxi(:) , leta(:), lzeta(:) !interpolants

     Nxyz = mesh % elements(eID) % Nxyz

     spAxi   => NodalStorage(Nxyz(1))
     spAeta   => NodalStorage(Nxyz(2))
     spAzeta   => NodalStorage(Nxyz(3))

     allocate( lxi(0:Nxyz(1)), leta(0:Nxyz(2)), lzeta(0:Nxyz(3)) )

     lxi = spAxi % lj(xi(1))
     leta = spAeta % lj(xi(2))
     lzeta = spAzeta % lj(xi(3))

     !$acc update self(mesh%elements(eID)%Storage%Q) async(eid)
     !$acc wait(eid)

     Qe = 0.0_RP
     do k = 0, Nxyz(3)    ; do j = 0, Nxyz(2)  ; do i = 0, Nxyz(1)
         Qe = Qe + mesh % elements(eID) % Storage % Q(:,i,j,k) * lxi(i) * leta(j) * lzeta(k)
     end do               ; end do             ; end do

     deallocate(lxi,leta,lzeta)
     nullify(spAxi,spAeta,spAzeta)

END Function interpolateQAD

#endif
end module ActuatorDisk
