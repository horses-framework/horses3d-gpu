#include "Includes.h"
#ifdef FLOW
module ProbeClass
   use SMConstants
   use HexMeshClass
   use MonitorDefinitions
   use PhysicsStorage
   use VariableConversion
   use MPI_Process_Info
   use FluidData
   use FileReadingUtilities, only: getRealArrayFromString
   use NodalStorageClass   , only: NodalStorage
#ifdef _HAS_MPI_
   use mpi
#endif
   implicit none
   
   private
   public   Probe_t
!
!  **********************
!  Probe class definition
!  **********************
!
   type Probe_t
      logical                         :: active
      integer                         :: rank
      integer                         :: ID
      integer                         :: eID
      integer                         :: nVars
      real(kind=RP)                   :: x(NDIM)
      real(kind=RP)                   :: xi(NDIM)
      real(kind=RP), allocatable      :: values(:,:)
      real(kind=RP), allocatable      :: lxi(:) , leta(:), lzeta(:)
      real(kind=RP), allocatable      :: var(:,:,:)
      logical                         :: isFileProbe = .false.
      integer                         :: fileUnit = -1
      real(kind=RP)                   :: saveTimestep
      real(kind=RP)                   :: lastSavedTime
      character(len=STR_LEN_MONITORS) :: fileName
      character(len=STR_LEN_MONITORS) :: monitorName
      character(len=STR_LEN_MONITORS), allocatable :: variableNames(:)
      contains
         procedure   :: Initialization => Probe_Initialization
         procedure   :: Update         => Probe_Update
         procedure   :: ComputeLocal   => Probe_ComputeLocal
         procedure   :: WriteLabel     => Probe_WriteLabel
         procedure   :: WriteValues    => Probe_WriteValue
         procedure   :: WriteToFile    => Probe_WriteToFile
         procedure   :: LookInOtherPartitions => Probe_LookInOtherPartitions
         procedure   :: destruct       => Probe_Destruct
         procedure   :: copy           => Probe_Assign
         generic     :: assignment(=)  => copy
   end type Probe_t

   contains

      subroutine Probe_Initialization(self, mesh, ID, solution_file, FirstCall, x_in, variables_in, name_in, isFileProbe_in, outputFormat_in, eID_hint)
         use ParamfileRegions
         use MPI_Process_Info
         use Utilities, only: toLower
         implicit none
         class(Probe_t)          :: self
         class(HexMesh)          :: mesh
         integer                 :: ID
         character(len=*)        :: solution_file
         logical, intent(in)     :: FirstCall
         real(kind=RP),     intent(in), optional :: x_in(NDIM)
         character(len=*),  intent(in), optional :: variables_in(:)
         character(len=*),  intent(in), optional :: name_in
         logical,           intent(in), optional :: isFileProbe_in
         character(len=*),  intent(in), optional :: outputFormat_in
         integer,           intent(in), optional :: eID_hint
!
!        ---------------
!        Local variables
!        ---------------
!
         integer                          :: i, j, k, v, fid
         character(len=STR_LEN_MONITORS)  :: in_label
         character(len=STR_LEN_MONITORS)  :: fileName
         character(len=STR_LEN_MONITORS)  :: paramFile
         character(len=STR_LEN_MONITORS)  :: coordinates
         character(len=STR_LEN_MONITORS)  :: variable
         character(len=STR_LEN_MONITORS)  :: outputFormat

         self % isFileProbe = .false.
         if ( present(isFileProbe_in) ) self % isFileProbe = isFileProbe_in

         outputFormat = "ASCII"
         if ( present(outputFormat_in) ) outputFormat = trim(outputFormat_in)


         if (FirstCall) then
!
!           Get monitor ID
!           --------------
            self % ID = ID
!
!           Get the probe definition, either from a probes file (x_in/variables_in
!           provided) or from a "#define probe" block in the case file
!           ------------------------------------------------------------------------
            if ( present(x_in) ) then
               if ( .not. present(name_in) .or. .not. present(variables_in) ) then
                  write(*,'(A)') "ERROR: Probe_Initialization: x_in requires name_in and variables_in."
                  stop
               end if
               self % monitorName = name_in
               self % x           = x_in
               self % nVars       = size(variables_in)
               allocate( self % variableNames(self % nVars) )
               self % variableNames = variables_in

            else
!
!              Search for the parameters in the case file
!              ------------------------------------------
               write(in_label , '(A,I0)') "#define probe " , self % ID

               call get_command_argument(1, paramFile)
               call readValueInRegion(trim(paramFile), "name"    , self % monitorName, in_label, "# end" )
               call readValueInRegion(trim(paramFile), "variable", variable          , in_label, "# end" )
               call readValueInRegion(trim(paramFile), "position", coordinates       , in_label, "# end" )
!
!              Get the coordinates
!              -------------------
               self % x = getRealArrayFromString(coordinates)

               self % nVars = 1
               allocate( self % variableNames(1) )
               self % variableNames(1) = variable
            end if
!
!           Allocate memory
!           ---------------
            if (self % isFileProbe) then
               allocate ( self % values(self % nVars, 1) )
            else
               allocate ( self % values(self % nVars, BUFFER_SIZE) )
            end if
            self % saveTimestep  = 0.0_RP
            self % lastSavedTime = -huge(self % lastSavedTime)
!
!           Check the variables
!           --------------------
            do v = 1, self % nVars
            call tolower(self % variableNames(v))

            select case ( trim(self % variableNames(v)) )
#ifdef NAVIERSTOKES
            case ("pressure")
            case ("velocity")
            case ("u")
            case ("v")
            case ("w")
            case ("mach")
            case ("k")
            case ("rho")
            case default
               print*, 'Probe variable "',trim(self % variableNames(v)),'" not implemented.'
               print*, "Options available are:"
               print*, "   * pressure"
               print*, "   * velocity"
               print*, "   * u"
               print*, "   * v"
               print*, "   * w"
               print*, "   * Mach"
               print*, "   * K"
               print*, "   * rho"
            end select
#endif
#ifdef INCNS
            case ("pressure")
            case ("velocity")
            case ("u")
            case ("v")
            case ("w")
            case ("rho")
            case default
               print*, 'Probe variable "',trim(self % variableNames(v)),'" not implemented.'
               print*, "Options available are:"
               print*, "   * pressure"
               print*, "   * velocity"
               print*, "   * u"
               print*, "   * v"
               print*, "   * w"
               print*, "   * rho"
            end select
#endif
#ifdef MULTIPHASE
            case ("static-pressure")

            case default
               print*, 'Probe variable "',trim(self % variableNames(v)),'" not implemented.'
               print*, "Options available are:"
               print*, "   * static-pressure"

            end select
#endif
#ifdef ACOUSTIC
            case ("pressure")
            case ("density")
            case ("u")
            case ("v")
            case ("w")
            case default
               print*, 'Probe variable "',trim(self % variableNames(v)),'" not implemented.'
               print*, "Options available are:"
               print*, "   * pressure"
               print*, "   * velocity"
               print*, "   * u"
               print*, "   * v"
               print*, "   * w"
            end select
#endif
            end do

!
!           Find the requested point in the mesh (use hint for local search if available)
!           ---------------------------------------------------------------------------
            self % active = mesh % FindPointWithCoords(self % x, self % eID, self % xi, eID_hint=eID_hint)
!
!           Check whether the probe is located in other partition
!           -----------------------------------------------------
            call self % LookInOtherPartitions
!
!           Disable the probe if the point is not found
!           -------------------------------------------
            if ( .not. self % active ) then
               if ( MPI_Process % isRoot ) then
                  write(STD_OUT,'(A,I0,A)') "Probe ", ID, " was not successfully initialized."
                  print*, "Probe is set to inactive."
               end if

               return
            end if
!
!           Set the fileName
!           ----------------
            write( self % fileName , '(A,A,A,A)') trim(solution_file) , "." , &
                                               trim(self % monitorName) , ".probe"

            if ( MPI_Process % isRoot .and. .not. self % isFileProbe ) then
               write(STD_OUT,'(/,30X,A,I0)') "** Probe ", self % ID
               write(STD_OUT,'(30X,A,A28,A)') "   ->", "Name: ", trim(self % monitorName)
               write(STD_OUT,'(30X,A,A28,3ES12.4)') "   ->", "Position: ", self % x
               write(STD_OUT,'(30X,A,A28)',advance="no") "   ->", "Variable(s): "
               do v = 1, self % nVars
                  write(STD_OUT,'(A,A)',advance="no") trim(self % variableNames(v)), "  "
               end do
               write(STD_OUT,*)
            end if
         end if
!
!
!           Return if the process does not contain the partition
!           ----------------------------------------------------
         if ( self % rank .ne. MPI_Process % rank ) then
            self % eID = 1
            return
         end if
         
!
!        If this is not the first call, just reload the reference frame coordinates
!        --------------------------------------------------------------------------
         if (.not. firstCall) self % active = mesh % elements(self % eID) % FindPointWithCoords(self % x,mesh % dir2D_ctrl, self % xi)
!
!        Get the Lagrange interpolants
!        -----------------------------
         associate(e => mesh % elements(self % eID))
         associate( spAxi   => NodalStorage(e % Nxyz(1)), &
                    spAeta  => NodalStorage(e % Nxyz(2)), &
                    spAzeta => NodalStorage(e % Nxyz(3)) )
         safedeallocate(self % lxi  ) ; allocate( self % lxi(0 : e % Nxyz(1)) )
         safedeallocate(self % leta ) ; allocate( self % leta(0 : e % Nxyz(2)) )
         safedeallocate(self % lzeta) ; allocate( self % lzeta(0 : e % Nxyz(3)) )
         self % lxi = spAxi % lj(self % xi(1))
         self % leta = spAeta % lj(self % xi(2))
         self % lzeta = spAzeta % lj(self % xi(3))
!
!        Allocate storage for the probe
!        -----------------------------   
         safedeallocate(self % var  ) ; allocate( self % var(0 : e % Nxyz(1),0 : e % Nxyz(2),0 : e % Nxyz(3)) )
         self % var = 0.0_RP
         
         ! File-probes use SoA arrays in Monitor_t (Monitor_InitFileProbesGPU);
         ! skip per-probe GPU copyin to avoid O(N) individual transfer overhead.
         if ( .not. self % isFileProbe ) then
            !$acc enter data copyin(self)
            !$acc enter data copyin(self % eiD)
            !$acc enter data copyin(self % id)
            !$acc enter data copyin(self % var)
            !$acc enter data copyin(self % lxi)
            !$acc enter data copyin(self % leta)
            !$acc enter data copyin(self % lzeta)
         end if
!
!        ****************
!        Prepare the file
!        ****************
!
!        Create file (skip for file-probes using HDF5 output)
!        -----------------------------------------------------
         if (FirstCall .and. .not. (self % isFileProbe .and. trim(outputFormat) .eq. "HDF5")) then
            open ( newunit = fID , file = trim(self % fileName) , status = "unknown" , action = "write" )
!
!        Write the file headers
!        ----------------------
            write( fID , '(A20,A  )') "Monitor name:      ", trim(self % monitorName)
            write( fID , '(A25,ES24.10,2(4X,ES24.10))') "x, y, z coordinates: ", self % x(1), self % x(2), self % x(3)

            write( fID , * )
            write( fID , '(A10,2X,A24)' , advance = "no") "Iteration" , "Time"
            do v = 1 , self % nVars
               write( fID , '(2X,A24)' , advance = "no") trim(self % variableNames(v))
            end do
            write( fID , * )

            if ( self % isFileProbe ) then
               ! Keep file open for low-overhead per-timestep appends
               self % fileUnit = fID
            else
               close ( fID )
            end if
         end if
         end associate
         end associate
      end subroutine Probe_Initialization

      subroutine Probe_Update(self, mesh, bufferPosition)
         use Physics
         use MPI_Process_Info
         implicit none
         class(Probe_t) :: self
         type(HexMesh)  :: mesh
         integer        :: bufferPosition
!
!        ---------------
!        Local variables
!        ---------------
!
         integer        :: i, j, k, v, ierr
         real(kind=RP)  :: value

         if ( .not. self % active ) return

         if ( MPI_Process % rank .eq. self % rank ) then

!
!           Update the probe
!           ----------------
            do v = 1, self % nVars

            select case (trim(self % variableNames(v)))
#ifdef NAVIERSTOKES
            case("pressure")
               !$acc parallel loop collapse(3) present(mesh,self) async(self % ID)
               do k = 0, mesh % elements(self % eID) % Nxyz(3) ; do j = 0, mesh % elements(self % eID) % Nxyz(2)  ; do i = 0, mesh % elements(self % eID) % Nxyz(1) 
                  self % var(i,j,k) = Pressure(mesh % elements(self % eID) % storage % Q(:,i,j,k))
               end do            ; end do             ; end do
               !$acc end parallel loop
   
            case("velocity")
               !$acc parallel loop collapse(3) present(mesh, self) async(self % ID)
               do k = 0, mesh % elements(self % eID) % Nxyz(3) ; do j = 0, mesh % elements(self % eID) % Nxyz(2)  ; do i = 0, mesh % elements(self % eID) % Nxyz(1) 
                  self % var(i,j,k) = sqrt(POW2(mesh % elements(self % eID) % storage % Q(IRHOU,i,j,k)) + &
                                           POW2(mesh % elements(self % eID) % storage % Q(IRHOV,i,j,k)) + &
                                           POW2(mesh % elements(self % eID) % storage % Q(IRHOW,i,j,k)))/mesh % elements(self % eID) % storage % Q(IRHO,i,j,k)
               end do         ; end do         ; end do
               !$acc end parallel loop
   
            case("u")
               !$acc parallel loop collapse(3) present(mesh,self) async(self % ID)
               do k = 0, mesh % elements(self % eID) % Nxyz(3) ; do j = 0, mesh % elements(self % eID) % Nxyz(2)  ; do i = 0, mesh % elements(self % eID) % Nxyz(1) 
                  self % var(i,j,k) = mesh % elements(self % eID) % storage % Q(IRHOU,i,j,k) / mesh % elements(self % eID) % storage % Q(IRHO,i,j,k)
               end do            ; end do             ; end do
               !$acc end parallel loop

            case("v")
               !$acc parallel loop collapse(3) present(mesh,self) async(self % ID)
               do k = 0, mesh % elements(self % eID) % Nxyz(3) ; do j = 0, mesh % elements(self % eID) % Nxyz(2)  ; do i = 0, mesh % elements(self % eID) % Nxyz(1) 
                  self % var(i,j,k) = mesh % elements(self % eID) % storage % Q(IRHOV,i,j,k) / mesh % elements(self % eID) % storage % Q(IRHO,i,j,k)
               end do            ; end do             ; end do
               !$acc end parallel loop
   
            case("w")
               !$acc parallel loop collapse(3) present(mesh,self) async(self % ID)
               do k = 0, mesh % elements(self % eID) % Nxyz(3) ; do j = 0, mesh % elements(self % eID) % Nxyz(2)  ; do i = 0, mesh % elements(self % eID) % Nxyz(1) 
                  self % var(i,j,k) = mesh % elements(self % eID) % storage % Q(IRHOW,i,j,k) / mesh % elements(self % eID) % storage % Q(IRHO,i,j,k)
               end do            ; end do             ; end do
               !$acc end parallel loop
   
            case("mach")
               !$acc parallel loop collapse(3) present(mesh,self) async(self % ID)
               do k = 0, mesh % elements(self % eID) % Nxyz(3) ; do j = 0, mesh % elements(self % eID) % Nxyz(2)  ; do i = 0, mesh % elements(self % eID) % Nxyz(1) 
                  self % var(i,j,k) = (POW2(mesh % elements(self % eID) % storage % Q(IRHOU,i,j,k)) + &
                                       POW2(mesh % elements(self % eID) % storage % Q(IRHOV,i,j,k)) + &
                                       POW2(mesh % elements(self % eID) % storage % Q(IRHOW,i,j,k)))/POW2(mesh % elements(self % eID) % storage % Q(IRHO,i,j,k))   ! Vabs**2
                  self % var(i,j,k) = sqrt( self % var(i,j,k) / ( thermodynamics % gamma*(thermodynamics % gamma-1.0_RP)*&
                                           (mesh % elements(self % eID) % storage % Q(IRHOE,i,j,k)/mesh % elements(self % eID) % storage % Q(IRHO,i,j,k)-0.5_RP * self % var(i,j,k)) ) )
               end do         ; end do         ; end do
               !$acc end parallel loop
      
            case("k")
               !$acc parallel loop collapse(3) present(mesh,self) async(self % ID)
               do k = 0, mesh % elements(self % eID) % Nxyz(3) ; do j = 0, mesh % elements(self % eID) % Nxyz(2)  ; do i = 0, mesh % elements(self % eID) % Nxyz(1)
                  self % var(i,j,k) = 0.5_RP * (POW2(mesh % elements(self % eID) % storage % Q(IRHOU,i,j,k)) + &
                                                POW2(mesh % elements(self % eID) % storage % Q(IRHOV,i,j,k)) + &
                                                POW2(mesh % elements(self % eID) % storage % Q(IRHOW,i,j,k)))/mesh % elements(self % eID) % storage % Q(IRHO,i,j,k)
               end do         ; end do         ; end do
               !$acc end parallel loop

            case("rho")
               !$acc parallel loop collapse(3) present(mesh,self) async(self % ID)
               do k = 0, mesh % elements(self % eID) % Nxyz(3) ; do j = 0, mesh % elements(self % eID) % Nxyz(2)  ; do i = 0, mesh % elements(self % eID) % Nxyz(1)
                  self % var(i,j,k) = mesh % elements(self % eID) % storage % Q(IRHO,i,j,k)
               end do         ; end do         ; end do
               !$acc end parallel loop
#endif
#ifdef INCNS

            case("pressure")
               !$acc parallel loop collapse(3) present(mesh,self) async(self % ID)
               do k = 0, mesh % elements(self % eID) % Nxyz(3) ; do j = 0, mesh % elements(self % eID) % Nxyz(2)  ; do i = 0, mesh % elements(self % eID) % Nxyz(1) 
                  self % var(i,j,k) = mesh % elements(self % eID) % storage % Q(INSP,i,j,k)
               end do            ; end do             ; end do
               !$acc end parallel loop
   
            case("velocity")
               !$acc parallel loop collapse(3) present(mesh,self) async(self % ID)
               do k = 0, mesh % elements(self % eID) % Nxyz(3) ; do j = 0, mesh % elements(self % eID) % Nxyz(2)  ; do i = 0, mesh % elements(self % eID) % Nxyz(1) 
                  self % var(i,j,k) = sqrt(POW2(mesh % elements(self % eID) % storage % Q(INSRHOU,i,j,k)) + &
                                           POW2(mesh % elements(self % eID) % storage % Q(INSRHOV,i,j,k)) + &
                                           POW2( mesh % elements(self % eID) % storage % Q(INSRHOW,i,j,k)))/mesh % elements(self % eID) % storage % Q(INSRHO,i,j,k)
               end do         ; end do         ; end do
               !$acc end parallel loop
   
            case("u")
               !$acc parallel loop collapse(3) present(mesh,self) async(self % ID)
               do k = 0, mesh % elements(self % eID) % Nxyz(3) ; do j = 0, mesh % elements(self % eID) % Nxyz(2)  ; do i = 0, mesh % elements(self % eID) % Nxyz(1) 
                  self % var(i,j,k) = mesh % elements(self % eID) % storage % Q(INSRHOU,i,j,k) / mesh % elements(self % eID) % storage % Q(INSRHO,i,j,k)
               end do            ; end do             ; end do
               !$acc end parallel loop
   
            case("v")
               !$acc parallel loop collapse(3) present(mesh,self) async(self % ID)
               do k = 0, mesh % elements(self % eID) % Nxyz(3) ; do j = 0, mesh % elements(self % eID) % Nxyz(2)  ; do i = 0, mesh % elements(self % eID) % Nxyz(1) 
                  self % var(i,j,k) = mesh % elements(self % eID) % storage % Q(INSRHOV,i,j,k) / mesh % elements(self % eID) % storage % Q(INSRHO,i,j,k)
               end do            ; end do             ; end do
               !$acc end parallel loop
   
            case("w")
               !$acc parallel loop collapse(3) present(mesh,self) async(self % ID)
               do k = 0, mesh % elements(self % eID) % Nxyz(3) ; do j = 0, mesh % elements(self % eID) % Nxyz(2)  ; do i = 0, mesh % elements(self % eID) % Nxyz(1)
                  self % var(i,j,k) = mesh % elements(self % eID) % storage % Q(INSRHOW,i,j,k) / mesh % elements(self % eID) % storage % Q(INSRHO,i,j,k)
               end do            ; end do             ; end do
               !$acc end parallel loop

            case("rho")
               !$acc parallel loop collapse(3) present(mesh,self) async(self % ID)
               do k = 0, mesh % elements(self % eID) % Nxyz(3) ; do j = 0, mesh % elements(self % eID) % Nxyz(2)  ; do i = 0, mesh % elements(self % eID) % Nxyz(1)
                  self % var(i,j,k) = mesh % elements(self % eID) % storage % Q(INSRHO,i,j,k)
               end do            ; end do             ; end do
               !$acc end parallel loop
#endif
#ifdef MULTIPHASE
            case("static-pressure")
               !$acc parallel loop collapse(3) present(mesh,self) async(self % ID)
               do k = 0, mesh % elements(self % eID) % Nxyz(3) ; do j = 0, mesh % elements(self % eID) % Nxyz(2)  ; do i = 0, mesh % elements(self % eID) % Nxyz(1) 
                  self % var(i,j,k) = mesh % elements(self % eID) % storage % Q(IMP,i,j,k) + mesh % elements(self % eID) % storage % Q(IMC,i,j,k)*mesh % elements(self % eID) % storage % mu(1,i,j,k) &
                               - 12.0_RP*multiphase%sigma*multiphase%invEps*(POW2(mesh % elements(self % eID) % storage % Q(IMC,i,j,k)*(1.0_RP-mesh % elements(self % eID) % storage % Q(IMC,i,j,k)))) &
                               - 0.25_RP*3.0_RP*multiphase % sigma * multiphase % eps * (POW2(mesh % elements(self % eID) % storage % c_x(1,i,j,k))+POW2(mesh % elements(self % eID) % storage % c_y(1,i,j,k))+POW2(mesh % elements(self % eID) % storage % c_z(1,i,j,k)))
               end do         ; end do         ; end do
               !$acc end parallel loop
#endif 
#ifdef ACOUSTIC
            case("pressure")
               do k = 0, mesh % elements(self % eID) % Nxyz(3) ; do j = 0, mesh % elements(self % eID) % Nxyz(2)  ; do i = 0, mesh % elements(self % eID) % Nxyz(1) 
                  self % var(i,j,k) = Q(ICAAP,i,j,k)
               end do            ; end do             ; end do

            case("density")
               do k = 0, mesh % elements(self % eID) % Nxyz(3) ; do j = 0, mesh % elements(self % eID) % Nxyz(2)  ; do i = 0, mesh % elements(self % eID) % Nxyz(1) 
                  self % var(i,j,k) = Q(ICAARHO,i,j,k)
               end do            ; end do             ; end do
   
            case("u")
               do k = 0, mesh % elements(self % eID) % Nxyz(3) ; do j = 0, mesh % elements(self % eID) % Nxyz(2)  ; do i = 0, mesh % elements(self % eID) % Nxyz(1) 
                  self % var(i,j,k) = Q(ICAAU,i,j,k)
               end do            ; end do             ; end do
   
            case("v")
               do k = 0, mesh % elements(self % eID) % Nxyz(3) ; do j = 0, mesh % elements(self % eID) % Nxyz(2)  ; do i = 0, mesh % elements(self % eID) % Nxyz(1) 
                  self % var(i,j,k) = Q(ICAAV,i,j,k)
               end do            ; end do             ; end do
   
            case("w")
               do k = 0, mesh % elements(self % eID) % Nxyz(3) ; do j = 0, mesh % elements(self % eID) % Nxyz(2)  ; do i = 0, mesh % elements(self % eID) % Nxyz(1) 
                  self % var(i,j,k) = Q(ICAAW,i,j,k)
               end do            ; end do             ; end do
#endif
            end select
   
            value = 0.0_RP
            !$acc parallel loop collapse(3) present(mesh, self) reduction(+:value) async(self % ID)
            do k = 0, mesh % elements(self % eID) % Nxyz(3) ; do j = 0, mesh % elements(self % eID) % Nxyz(2)  ; do i = 0, mesh % elements(self % eID) % Nxyz(1)
               value = value + self % var(i,j,k) * self % lxi(i) * self % leta(j) * self % lzeta(k)
            end do               ; end do             ; end do
            !$acc end parallel loop

            !$acc wait

            self % values(v, bufferPosition) = value

            end do
!
#ifdef _HAS_MPI_
            if ( MPI_Process % doMPIAction ) then
!
!              Share the result with the rest of the processes
!              -----------------------------------------------
               call mpi_bcast(self % values(:,bufferPosition), self % nVars, MPI_DOUBLE, self % rank, MPI_COMM_WORLD, ierr)

            end if
#endif
         else
!
!           Receive the result from the rank that contains the probe
!           --------------------------------------------------------
#ifdef _HAS_MPI_
            if ( MPI_Process % doMPIAction ) then
               call mpi_bcast(self % values(:,bufferPosition), self % nVars, MPI_DOUBLE, self % rank, MPI_COMM_WORLD, ierr)
            end if
#endif
         end if
      end subroutine Probe_Update

      subroutine Probe_ComputeLocal(self, mesh, bufferPosition)
!
!        *************************************************************
!           CPU-only computation for file-probes, no MPI.
!           Uses fused extract+Lagrange loops (single pass over Q)
!           to avoid writing to self%var and re-reading it.
!           Non-owning ranks set values to 0 so a caller can
!           accumulate results with a single MPI_Allreduce(SUM).
!        *************************************************************
!
         use Physics
         use MPI_Process_Info
         implicit none
         class(Probe_t) :: self
         type(HexMesh)  :: mesh
         integer        :: bufferPosition
!
!        ---------------
!        Local variables
!        ---------------
!
         integer        :: i, j, k, v
         real(kind=RP)  :: value, w
#ifdef NAVIERSTOKES
         real(kind=RP)  :: u2
#endif

         if ( .not. self % active ) then
            self % values(:, bufferPosition) = 0.0_RP
            return
         end if

         if ( MPI_Process % rank .ne. self % rank ) then
            self % values(:, bufferPosition) = 0.0_RP
            return
         end if

         do v = 1, self % nVars
            value = 0.0_RP
            select case (trim(self % variableNames(v)))
#ifdef NAVIERSTOKES
            case("pressure")
               do k = 0, mesh % elements(self % eID) % Nxyz(3) ; do j = 0, mesh % elements(self % eID) % Nxyz(2) ; do i = 0, mesh % elements(self % eID) % Nxyz(1)
                  w = self % lxi(i) * self % leta(j) * self % lzeta(k)
                  value = value + w * Pressure(mesh % elements(self % eID) % storage % Q(:,i,j,k))
               end do ; end do ; end do
            case("velocity")
               do k = 0, mesh % elements(self % eID) % Nxyz(3) ; do j = 0, mesh % elements(self % eID) % Nxyz(2) ; do i = 0, mesh % elements(self % eID) % Nxyz(1)
                  w = self % lxi(i) * self % leta(j) * self % lzeta(k)
                  value = value + w * sqrt(POW2(mesh % elements(self % eID) % storage % Q(IRHOU,i,j,k)) + &
                                           POW2(mesh % elements(self % eID) % storage % Q(IRHOV,i,j,k)) + &
                                           POW2(mesh % elements(self % eID) % storage % Q(IRHOW,i,j,k))) &
                                    / mesh % elements(self % eID) % storage % Q(IRHO,i,j,k)
               end do ; end do ; end do
            case("u")
               do k = 0, mesh % elements(self % eID) % Nxyz(3) ; do j = 0, mesh % elements(self % eID) % Nxyz(2) ; do i = 0, mesh % elements(self % eID) % Nxyz(1)
                  w = self % lxi(i) * self % leta(j) * self % lzeta(k)
                  value = value + w * mesh % elements(self % eID) % storage % Q(IRHOU,i,j,k) / mesh % elements(self % eID) % storage % Q(IRHO,i,j,k)
               end do ; end do ; end do
            case("v")
               do k = 0, mesh % elements(self % eID) % Nxyz(3) ; do j = 0, mesh % elements(self % eID) % Nxyz(2) ; do i = 0, mesh % elements(self % eID) % Nxyz(1)
                  w = self % lxi(i) * self % leta(j) * self % lzeta(k)
                  value = value + w * mesh % elements(self % eID) % storage % Q(IRHOV,i,j,k) / mesh % elements(self % eID) % storage % Q(IRHO,i,j,k)
               end do ; end do ; end do
            case("w")
               do k = 0, mesh % elements(self % eID) % Nxyz(3) ; do j = 0, mesh % elements(self % eID) % Nxyz(2) ; do i = 0, mesh % elements(self % eID) % Nxyz(1)
                  w = self % lxi(i) * self % leta(j) * self % lzeta(k)
                  value = value + w * mesh % elements(self % eID) % storage % Q(IRHOW,i,j,k) / mesh % elements(self % eID) % storage % Q(IRHO,i,j,k)
               end do ; end do ; end do
            case("mach")
               do k = 0, mesh % elements(self % eID) % Nxyz(3) ; do j = 0, mesh % elements(self % eID) % Nxyz(2) ; do i = 0, mesh % elements(self % eID) % Nxyz(1)
                  w  = self % lxi(i) * self % leta(j) * self % lzeta(k)
                  u2 = (POW2(mesh % elements(self % eID) % storage % Q(IRHOU,i,j,k)) + &
                        POW2(mesh % elements(self % eID) % storage % Q(IRHOV,i,j,k)) + &
                        POW2(mesh % elements(self % eID) % storage % Q(IRHOW,i,j,k))) / &
                        POW2(mesh % elements(self % eID) % storage % Q(IRHO,i,j,k))
                  value = value + w * sqrt( u2 / ( thermodynamics % gamma*(thermodynamics % gamma-1.0_RP) * &
                             (mesh % elements(self % eID) % storage % Q(IRHOE,i,j,k) / &
                              mesh % elements(self % eID) % storage % Q(IRHO,i,j,k) - 0.5_RP*u2) ) )
               end do ; end do ; end do
            case("k")
               do k = 0, mesh % elements(self % eID) % Nxyz(3) ; do j = 0, mesh % elements(self % eID) % Nxyz(2) ; do i = 0, mesh % elements(self % eID) % Nxyz(1)
                  w = self % lxi(i) * self % leta(j) * self % lzeta(k)
                  value = value + w * 0.5_RP * (POW2(mesh % elements(self % eID) % storage % Q(IRHOU,i,j,k)) + &
                                                 POW2(mesh % elements(self % eID) % storage % Q(IRHOV,i,j,k)) + &
                                                 POW2(mesh % elements(self % eID) % storage % Q(IRHOW,i,j,k))) / &
                                                 mesh % elements(self % eID) % storage % Q(IRHO,i,j,k)
               end do ; end do ; end do
            case("rho")
               do k = 0, mesh % elements(self % eID) % Nxyz(3) ; do j = 0, mesh % elements(self % eID) % Nxyz(2) ; do i = 0, mesh % elements(self % eID) % Nxyz(1)
                  w = self % lxi(i) * self % leta(j) * self % lzeta(k)
                  value = value + w * mesh % elements(self % eID) % storage % Q(IRHO,i,j,k)
               end do ; end do ; end do
#endif
#ifdef INCNS
            case("pressure")
               do k = 0, mesh % elements(self % eID) % Nxyz(3) ; do j = 0, mesh % elements(self % eID) % Nxyz(2) ; do i = 0, mesh % elements(self % eID) % Nxyz(1)
                  w = self % lxi(i) * self % leta(j) * self % lzeta(k)
                  value = value + w * mesh % elements(self % eID) % storage % Q(INSP,i,j,k)
               end do ; end do ; end do
            case("velocity")
               do k = 0, mesh % elements(self % eID) % Nxyz(3) ; do j = 0, mesh % elements(self % eID) % Nxyz(2) ; do i = 0, mesh % elements(self % eID) % Nxyz(1)
                  w = self % lxi(i) * self % leta(j) * self % lzeta(k)
                  value = value + w * sqrt(POW2(mesh % elements(self % eID) % storage % Q(INSRHOU,i,j,k)) + &
                                           POW2(mesh % elements(self % eID) % storage % Q(INSRHOV,i,j,k)) + &
                                           POW2(mesh % elements(self % eID) % storage % Q(INSRHOW,i,j,k))) &
                                    / mesh % elements(self % eID) % storage % Q(INSRHO,i,j,k)
               end do ; end do ; end do
            case("u")
               do k = 0, mesh % elements(self % eID) % Nxyz(3) ; do j = 0, mesh % elements(self % eID) % Nxyz(2) ; do i = 0, mesh % elements(self % eID) % Nxyz(1)
                  w = self % lxi(i) * self % leta(j) * self % lzeta(k)
                  value = value + w * mesh % elements(self % eID) % storage % Q(INSRHOU,i,j,k) / mesh % elements(self % eID) % storage % Q(INSRHO,i,j,k)
               end do ; end do ; end do
            case("v")
               do k = 0, mesh % elements(self % eID) % Nxyz(3) ; do j = 0, mesh % elements(self % eID) % Nxyz(2) ; do i = 0, mesh % elements(self % eID) % Nxyz(1)
                  w = self % lxi(i) * self % leta(j) * self % lzeta(k)
                  value = value + w * mesh % elements(self % eID) % storage % Q(INSRHOV,i,j,k) / mesh % elements(self % eID) % storage % Q(INSRHO,i,j,k)
               end do ; end do ; end do
            case("w")
               do k = 0, mesh % elements(self % eID) % Nxyz(3) ; do j = 0, mesh % elements(self % eID) % Nxyz(2) ; do i = 0, mesh % elements(self % eID) % Nxyz(1)
                  w = self % lxi(i) * self % leta(j) * self % lzeta(k)
                  value = value + w * mesh % elements(self % eID) % storage % Q(INSRHOW,i,j,k) / mesh % elements(self % eID) % storage % Q(INSRHO,i,j,k)
               end do ; end do ; end do
            case("rho")
               do k = 0, mesh % elements(self % eID) % Nxyz(3) ; do j = 0, mesh % elements(self % eID) % Nxyz(2) ; do i = 0, mesh % elements(self % eID) % Nxyz(1)
                  w = self % lxi(i) * self % leta(j) * self % lzeta(k)
                  value = value + w * mesh % elements(self % eID) % storage % Q(INSRHO,i,j,k)
               end do ; end do ; end do
#endif
#ifdef MULTIPHASE
            case("static-pressure")
               do k = 0, mesh % elements(self % eID) % Nxyz(3) ; do j = 0, mesh % elements(self % eID) % Nxyz(2) ; do i = 0, mesh % elements(self % eID) % Nxyz(1)
                  w = self % lxi(i) * self % leta(j) * self % lzeta(k)
                  value = value + w * ( mesh % elements(self % eID) % storage % Q(IMP,i,j,k) &
                               + mesh % elements(self % eID) % storage % Q(IMC,i,j,k)*mesh % elements(self % eID) % storage % mu(1,i,j,k) &
                               - 12.0_RP*multiphase%sigma*multiphase%invEps*(POW2(mesh % elements(self % eID) % storage % Q(IMC,i,j,k)*(1.0_RP-mesh % elements(self % eID) % storage % Q(IMC,i,j,k)))) &
                               - 0.25_RP*3.0_RP*multiphase%sigma*multiphase%eps*(POW2(mesh % elements(self % eID) % storage % c_x(1,i,j,k))+POW2(mesh % elements(self % eID) % storage % c_y(1,i,j,k))+POW2(mesh % elements(self % eID) % storage % c_z(1,i,j,k))) )
               end do ; end do ; end do
#endif
#ifdef ACOUSTIC
            case("pressure")
               do k = 0, mesh % elements(self % eID) % Nxyz(3) ; do j = 0, mesh % elements(self % eID) % Nxyz(2) ; do i = 0, mesh % elements(self % eID) % Nxyz(1)
                  w = self % lxi(i) * self % leta(j) * self % lzeta(k)
                  value = value + w * Q(ICAAP,i,j,k)
               end do ; end do ; end do
            case("density")
               do k = 0, mesh % elements(self % eID) % Nxyz(3) ; do j = 0, mesh % elements(self % eID) % Nxyz(2) ; do i = 0, mesh % elements(self % eID) % Nxyz(1)
                  w = self % lxi(i) * self % leta(j) * self % lzeta(k)
                  value = value + w * Q(ICAARHO,i,j,k)
               end do ; end do ; end do
            case("u")
               do k = 0, mesh % elements(self % eID) % Nxyz(3) ; do j = 0, mesh % elements(self % eID) % Nxyz(2) ; do i = 0, mesh % elements(self % eID) % Nxyz(1)
                  w = self % lxi(i) * self % leta(j) * self % lzeta(k)
                  value = value + w * Q(ICAAU,i,j,k)
               end do ; end do ; end do
            case("v")
               do k = 0, mesh % elements(self % eID) % Nxyz(3) ; do j = 0, mesh % elements(self % eID) % Nxyz(2) ; do i = 0, mesh % elements(self % eID) % Nxyz(1)
                  w = self % lxi(i) * self % leta(j) * self % lzeta(k)
                  value = value + w * Q(ICAAV,i,j,k)
               end do ; end do ; end do
            case("w")
               do k = 0, mesh % elements(self % eID) % Nxyz(3) ; do j = 0, mesh % elements(self % eID) % Nxyz(2) ; do i = 0, mesh % elements(self % eID) % Nxyz(1)
                  w = self % lxi(i) * self % leta(j) * self % lzeta(k)
                  value = value + w * Q(ICAAW,i,j,k)
               end do ; end do ; end do
#endif
            end select
            self % values(v, bufferPosition) = value
         end do

      end subroutine Probe_ComputeLocal

      subroutine Probe_WriteLabel ( self )
!
!        *************************************************************
!              This subroutine writes the label for the probe,
!           when invoked from the time integrator Display
!           procedure.
!        *************************************************************
!
         implicit none
         class(Probe_t)             :: self

         write(STD_OUT , '(3X,A10)' , advance = "no") trim(self % monitorName(1 : MONITOR_LENGTH))

      end subroutine Probe_WriteLabel

      subroutine Probe_WriteValue ( self , bufferLine ) 
!
!        *************************************************************
!              This subroutine writes the monitor value for the time
!           integrator Display procedure.
!        *************************************************************
!
         implicit none
         class(Probe_t) :: self
         integer                 :: bufferLine

         write(STD_OUT , '(1X,A,1X,ES10.3)' , advance = "no") "|" , self % values ( 1 , bufferLine )

      end subroutine Probe_WriteValue

      subroutine Probe_WriteToFile ( self , iter , t , no_of_lines)
!
!        *************************************************************
!              This subroutine writes the buffer to the file.
!        *************************************************************
!
         implicit none  
         class(Probe_t) :: self
         integer                 :: iter(:)
         real(kind=RP)           :: t(:)
         integer                 :: no_of_lines
!
!        ---------------
!        Local variables
!        ---------------
!
         integer                    :: i, v
         integer                    :: fID

         if ( .not. self % active ) then
            if ( no_of_lines .ne. 0 ) self % values(:,1) = self % values(:,no_of_lines)
            return
         end if

         if ( MPI_Process % isRoot ) then
            do i = 1 , no_of_lines
               if ( self % saveTimestep > 0.0_RP ) then
                  if ( t(i) < self % lastSavedTime + self % saveTimestep ) cycle
               end if
               if ( self % fileUnit >= 0 ) then
                  fID = self % fileUnit
               else
                  open( newunit = fID , file = trim ( self % fileName ) , action = "write" , access = "append" , status = "old" )
               end if
               write( fID , '(I10,2X,ES24.16)' , advance = "no" ) iter(i) , t(i)
               do v = 1 , self % nVars
                  write( fID , '(2X,ES24.16)' , advance = "no" ) self % values(v,i)
               end do
               write( fID , * )
               if ( self % fileUnit < 0 ) close ( fID )
               self % lastSavedTime = t(i)
            end do
         end if


         if ( no_of_lines .ne. 0 ) self % values(:,1) = self % values(:,no_of_lines)
      
      end subroutine Probe_WriteToFile

      subroutine Probe_LookInOtherPartitions(self)
         use MPI_Process_Info
         implicit none
         class(Probe_t)    :: self
         integer           :: allActives(MPI_Process % nProcs)
         integer           :: active, ierr

         if ( MPI_Process % doMPIAction ) then
#ifdef _HAS_MPI_
!
!           Cast the logicals onto integers
!           -------------------------------
            if ( self % active ) then
               active = 1
            else
               active = 0
            end if
!
!           Gather all data from all processes
!           ----------------------------------
            call mpi_allgather(active, 1, MPI_INT, allActives, 1, MPI_INT, MPI_COMM_WORLD, ierr)
!
!           Check if any of them found the probe
!           ------------------------------------
            if ( any(allActives .eq. 1) ) then
!
!              Assign the domain of the partition that contains the probe
!              ----------------------------------------------------------
               self % active = .true.
               self % rank = maxloc(allActives, dim = 1) - 1

            else
!
!              Disable the probe
!              -----------------
               self % active = .false.
               self % rank   = -1

            end if
#endif
         else
!
!           Without MPI select the rank 0 as default
!           ----------------------------------------
            self % rank = 0

         end if

      end subroutine Probe_LookInOtherPartitions
      
      elemental subroutine Probe_Destruct (self)
         implicit none
         class(Probe_t), intent(inout) :: self
         
         safedeallocate (self % values)
         safedeallocate (self % lxi)
         safedeallocate (self % leta)
         safedeallocate (self % lzeta)
         safedeallocate (self % variableNames)
      end subroutine Probe_Destruct
      
      elemental subroutine Probe_Assign (to, from)
         implicit none
         class(Probe_t), intent(inout) :: to
         type(Probe_t) , intent(in) :: from
         
         to % active = from % active
         to % rank = from % rank
         to % ID = from %  ID
         to % eID = from % eID
         to % nVars = from % nVars
         to % x = from % x
         to % xi = from % xi

         if ( allocated(from % values) ) then
            safedeallocate ( to % values )
            allocate ( to % values ( size(from % values, 1), size(from % values, 2) ) )
            to % values = from % values
         end if

         if ( allocated(from % lxi) ) then
            safedeallocate ( to % lxi )
            allocate ( to % lxi ( size(from % lxi) ) )
            to % lxi = from % lxi
         end if

         if ( allocated(from % leta) ) then
            safedeallocate ( to % leta )
            allocate ( to % leta ( size(from % leta) ) )
            to % leta = from % leta
         end if

         if ( allocated(from % lzeta) ) then
            safedeallocate ( to % lzeta )
            allocate ( to % lzeta ( size(from % lzeta) ) )
            to % lzeta = from % lzeta
         end if
         
         to % isFileProbe   = from % isFileProbe
         to % fileUnit      = from % fileUnit
         to % saveTimestep  = from % saveTimestep
         to % lastSavedTime = from % lastSavedTime
         to % fileName = from % fileName
         to % monitorName = from % monitorName

         safedeallocate ( to % variableNames )
         allocate ( to % variableNames ( size(from % variableNames) ) )
         to % variableNames = from % variableNames

      end subroutine Probe_Assign
      
end module ProbeClass
#endif
