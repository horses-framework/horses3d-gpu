#include "Includes.h"
Module ChannelForcing  
#if defined(NAVIERSTOKES)
    use SMConstants
    use MPI_Process_Info
    use PhysicsStorage
#ifdef _HAS_MPI_
    use mpi
#endif
    Implicit None

private
public initializeChannel, updateChannel, channelSource

!
!=====================================================================
!     Module variables 
!=====================================================================
!
      real(kind=RP) :: flowDir(NDIM)            ! Flow direction
      real(kind=RP) :: f_x(2), f_y(2), f_z(2)   ! Correction parameters
      real(kind=RP) :: dpdx, dpdy, dpdz         ! Driving pressure gradient
      
      logical :: channelIsActive

      ! --- Global rho-extremes accumulated over all RK substeps ---
      real(kind=RP) :: rho_min_global, rho_max_global

      ! parameters
      real(kind=RP), parameter :: alpha = 0.9_RP

!$acc declare create(f_x, f_y, f_z)

contains
!=====================================================================
!
!---------------------------------------------------------------------
!   INITIALIZATION
!---------------------------------------------------------------------
!
     Subroutine initializeChannel(controlVariables)
       use FTValueDictionaryClass
       use FluidData, only: refValues
       use Headers
       implicit none
       type(FTValueDictionary), intent(in) :: controlVariables
       integer :: ierr

       channelIsActive = controlVariables % logicalValueForKey("use channel")
       if (.not. channelIsActive) return

! --- flow direction ---
       flowDir(1) = cos(refValues % AoATheta*pi/180._RP) * cos(refValues % AoAPhi*pi/180._RP)
       flowDir(2) = sin(refValues % AoAtheta*pi/180._RP) * cos(refValues % AoAphi*pi/180._RP)
       flowDir(3) = sin(refValues % AoAphi*pi/180._RP)

! --- initial dp forcing ---
       dpdx  = 3.08e-3_RP
       dpdy  = 0.0_RP 
       dpdz  = 3.08e-5_RP

! --- initial turbulent forcing parameters ---
       f_x(1) = 0.0_RP
       f_x(2) = (1.0_RP - alpha) * dpdx

       f_y = 0.0_RP
       f_z = 0.0_RP


! --- initialize rho extrema ---
       rho_min_global = huge(1.0_RP)
       rho_max_global = -huge(1.0_RP)

!$acc update device(f_x, f_y, f_z)

       if (.not. MPI_Process % isRoot) return
       call Subsection_Header("Channel")
       write(STD_OUT,'(30X,A,A28)') "->","Channel activated"
     End Subroutine initializeChannel

!=====================================================================
!
!---------------------------------------------------------------------
!   UPDATE CHANNEL FORCING  (called once per timestep)
!---------------------------------------------------------------------
!
     Subroutine updateChannel(mesh, dt, monitors)
        use HexMeshClass
        use MonitorsClass
        implicit none

        type(HexMesh)               :: mesh
        real(kind=RP)               :: dt
        type(Monitor_t), intent(in) :: monitors

        real(kind=RP) :: meanMomentum(NDIM)
        real(kind=RP) :: S_min, S_max
        real(kind=RP) :: dp_vals(3)
        integer       :: ierr, dp_ind

        if (.not. channelIsActive) return

! --- collect min/max rho across all partitions ---
#ifdef _HAS_MPI_
        if (MPI_Process % doMPIAction) then
            call MPI_Allreduce(MPI_IN_PLACE, rho_min_global, 1, MPI_DOUBLE_PRECISION, MPI_MIN, MPI_COMM_WORLD, ierr)
            call MPI_Allreduce(MPI_IN_PLACE, rho_max_global, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr)
        end if
#endif

! --- compute mean momentum ---
        meanMomentum = monitors % volumeMonitors(1) % getLast()

! --- compute forcing parameters ---
        f_x(1) = (flowDir(1) - meanMomentum(1)) * alpha / dt
        f_x(2) = (1._RP - alpha) * dpdx

        f_y(1) = (flowDir(2) - meanMomentum(2)) * alpha / dt
        f_y(2) = (1._RP - alpha) * dpdy

        f_z(1) = (flowDir(3) - meanMomentum(3)) * alpha / dt
        f_z(2) = (1._RP - alpha) * dpdz

! --- update dpdx using extreme rho values ---
        S_min = f_x(1)*rho_min_global + f_x(2)
        S_max = f_x(1)*rho_max_global + f_x(2)
        dp_vals = [dpdx, S_max, S_min]
        dp_ind = maxloc(abs(dp_vals), 1)
        dpdx = dp_vals(dp_ind)

! --- update dpdy (here zero for now but kept consistent) ---
        S_min = f_y(1)*rho_min_global + f_y(2)
        S_max = f_y(1)*rho_max_global + f_y(2)
        dp_vals = [dpdy, S_max, S_min]
        dp_ind = maxloc(abs(dp_vals), 1)
        dpdy = dp_vals(dp_ind)

! --- update dpdz ---
        S_min = f_z(1)*rho_min_global + f_z(2)
        S_max = f_z(1)*rho_max_global + f_z(2)
        dp_vals = [dpdz, S_max, S_min]
        dp_ind = maxloc(abs(dp_vals), 1)
        dpdz = dp_vals(dp_ind)

! --- send to device ---
!$acc update device(f_x, f_y, f_z)

! --- reset extrema for next timestep ---
        rho_min_global = huge(1.0_RP)
        rho_max_global = -huge(1.0_RP)

     End Subroutine updateChannel

!=====================================================================
!
!---------------------------------------------------------------------
!   CHANNEL SOURCE (called each RK substep)
!---------------------------------------------------------------------
!
     Subroutine channelSource(mesh)
         use HexMeshClass
         implicit none

         type(HexMesh), intent(inout) :: mesh

         real(kind=RP) :: rho, Su, Sv, Sw, Se
         ! Local variables for reduction over all elements
         real(kind=RP) :: rho_min, rho_max
         ! Local variables for element min/max
         real(kind=RP) :: local_rho_min, local_rho_max
         integer :: i,j,k,eID

         if (.not. channelIsActive) return
          
         ! Initialize reduction variables with current global values
         rho_min = rho_min_global
         rho_max = rho_max_global

!$acc parallel loop gang present(mesh) num_gangs(9700) reduction(max:rho_max) reduction(min:rho_min)         
         do eID = 1, mesh % no_of_elements
            local_rho_min = huge(1.0_RP)
            local_rho_max = -huge(1.0_RP)

            !$acc loop vector collapse(3) reduction(min:local_rho_min) reduction(max:local_rho_max)
            do k = 0, mesh % elements(eID)%Nxyz(3)
               do j = 0, mesh % elements(eID)%Nxyz(2)
                  do i = 0, mesh % elements(eID)%Nxyz(1)

                     rho = mesh%elements(eID)%storage%Q(IRHO,i,j,k)

                     Su = f_x(1)*rho + f_x(2)
                     Sv = 0.0_RP
                     Sw = f_z(1)*rho + f_z(2)

                     Se = ( mesh%elements(eID)%storage%Q(IRHOU,i,j,k)*Su + &
                            mesh%elements(eID)%storage%Q(IRHOV,i,j,k)*Sv + &
                            mesh%elements(eID)%storage%Q(IRHOW,i,j,k)*Sw ) / rho

                     mesh%elements(eID)%storage%S_NS(IRHOU,i,j,k) = Su
                     mesh%elements(eID)%storage%S_NS(IRHOV,i,j,k) = Sv
                     mesh%elements(eID)%storage%S_NS(IRHOW,i,j,k) = Sw
                     mesh%elements(eID)%storage%S_NS(IRHOE,i,j,k) = Se

                     local_rho_min = min(local_rho_min, rho)
                     local_rho_max = max(local_rho_max, rho)

                  end do
               end do
            end do

             rho_min = min(rho_min, local_rho_min)
             rho_max = max(rho_max, local_rho_max)

         end do
!$acc end parallel loop

          ! Update global values
          rho_min_global = rho_min
          rho_max_global = rho_max

     End Subroutine channelSource

#endif
End Module ChannelForcing