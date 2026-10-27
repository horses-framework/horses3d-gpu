!
!////////////////////////////////////////////////////////////////////////
!
!      RK integrators for DG approximation to conservation
!      laws in 3D
!
!////////////////////////////////////////////////////////////////////////
!
!  GPU-port status (this branch):
!    TakeExplicitEulerStep ..... GPU-ported, VALIDATED
!    TakeRK3Step ............... GPU-ported, VALIDATED  
!    TakeRK5Step ............... GPU-ported, VALIDATED
!    TakeRKOptStep ............. GPU-ported, pending validation
!    TakeSSPRK33Step ........... GPU-ported (unlimited path), VALIDATED
!    TakeSSPRK43Step ........... GPU-ported (unlimited path), VALIDATED
!    stage_limiter ............. GPU-ported, VALIDATED
!    TakeEulerRK3Step .......... NOT GPU-ported 
!    TakeExplicitBDFStep ....... NOT GPU-ported

#include "Includes.h"
MODULE ExplicitMethods
   USE SMConstants
   use HexMeshClass
   use TimeIntegratorDefinitions
   use DGSEMClass, only: ComputeTimeDerivative_f
   use ParticlesClass
   use PhysicsStorage, only: NCONS, CTD_IGNORE_MODE
   IMPLICIT NONE

   private
   public   EULER_NAME, RK3_NAME, RK5_NAME, OPTRK_NAME, SSPRK33_NAME, SSPRK43_NAME, EULER_RK3_NAME
   public   EULER_KEY, RK3_KEY, RK5_KEY, OPTRK_KEY, SSPRK33_KEY, SSPRK43_KEY, EULER_RK3_KEY
   public   TakeExplicitEulerStep, TakeRK3Step, TakeRK5Step, TakeRKOptStep
   public   TakeSSPRK33Step, TakeSSPRK43Step, TakeEulerRK3Step
   public   Enable_CTD_AFTER_STEPS, Enable_limiter, CTD_AFTER_STEPS, LIMITED, LIMITER_MIN

   integer,  protected :: eBDF_order = 3
   logical,  protected :: CTD_AFTER_STEPS = .false.
   logical,  protected :: LIMITED = .false.
   real(RP), protected :: LIMITER_MIN = 1e-13_RP
!
!  Implemented integration methods
!  -------------------------------
!  Remember to add them to `Enable_limiter` if they support limiting (or any kind of stage
!  callbacks if that is ever implemented)
!  ---------------------------------------------------------------------------------------
   character(len=*), parameter :: EULER_NAME   = "euler"
   character(len=*), parameter :: RK3_NAME     = "rk3"
   character(len=*), parameter :: RK5_NAME     = "rk5"
   character(len=*), parameter :: OPTRK_NAME   = "optimal rk"
   character(len=*), parameter :: SSPRK33_NAME = "ssprk33"
   character(len=*), parameter :: SSPRK43_NAME = "ssprk43"
   character(len=*), parameter :: EULER_RK3_NAME = "euler rk3"

   integer, parameter :: EULER_KEY   = 1
   integer, parameter :: RK3_KEY     = 2
   integer, parameter :: RK5_KEY     = 3
   integer, parameter :: OPTRK_KEY   = 4
   integer, parameter :: SSPRK33_KEY = 5
   integer, parameter :: SSPRK43_KEY = 6
   integer, parameter :: EULER_RK3_KEY = 7
!========
 CONTAINS
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!  ------------------------------
!  Routine for taking an explicit Euler - RK3 step depending on the polynomial order of each element.
!  NOTE: NOT GPU-ported. Host OpenMP only. Uses LO_Elements/HO_Elements subsets.
!  ------------------------------
 SUBROUTINE TakeEulerRK3Step( mesh, particles, t, deltaT, ComputeTimeDerivative, dt_vec, dts, global_dt, iter)
!
!     ----------------------------------
!     Williamson's 3rd order Runge-Kutta
!     ----------------------------------
!
      IMPLICIT NONE
!
!     -----------------
!     Input parameters:
!     -----------------
!
      type(HexMesh)      :: mesh
#ifdef FLOW
      type(Particles_t)  :: particles
#else
      logical            :: particles
#endif
      REAL(KIND=RP)   :: t, deltaT, tk
      real(kind=RP), allocatable, dimension(:), intent(in), optional :: dt_vec
      procedure(ComputeTimeDerivative_f)    :: ComputeTimeDerivative
      logical, intent(in), optional :: dts
      real(kind=RP), intent(in), optional :: global_dt
      integer, intent(in), optional :: iter
!
!     ---------------
!     Local variables
!     ---------------
!
      REAL(KIND=RP), DIMENSION(3) :: a = (/0.0_RP       , -5.0_RP /9.0_RP , -153.0_RP/128.0_RP/)
      REAL(KIND=RP), DIMENSION(3) :: b = (/0.0_RP       ,  1.0_RP /3.0_RP ,    3.0_RP/4.0_RP  /)
      REAL(KIND=RP), DIMENSION(3) :: c = (/1.0_RP/3.0_RP,  15.0_RP/16.0_RP,    8.0_RP/15.0_RP /)


      INTEGER :: i, j, k, id, eID
      integer :: interval

      interval = 10 !Compute time derivative each interval for those elements whose pmax = 1

      if (present(dt_vec)) then

         do k = 1,3
            tk = t + b(k)*deltaT
            if ((k==1) .and. (mod(iter, interval) == 0)) then
               call ComputeTimeDerivative( mesh, particles, tk, CTD_IGNORE_MODE)
            else
               call ComputeTimeDerivative( mesh, particles, tk, CTD_IGNORE_MODE, .true.)
            endif
            if ( present(dts) ) then
               if (dts) call ComputePseudoTimeDerivative(mesh, tk, global_dt)
            end if

            if (k==1) then
!$omp parallel do schedule(runtime) private(eID)
               do id = 1, SIZE( mesh % LO_Elements )
                  ! Explicit Euler
                  eID = mesh % LO_Elements(id)
#ifdef FLOW
                  mesh % elements(eID) % storage % Q = mesh % elements(eID) % storage % Q + dt_vec(eID)*mesh % elements(eID) % storage % QDot
#endif

#if (defined(CAHNHILLIARD)) && (!defined(FLOW))
                  mesh % elements(eID) % storage % c = mesh % elements(eID) % storage % c + dt_vec(eID)*mesh % elements(eID) % storage % cDot
#endif
               end do ! id
!$omp end parallel do
            end if

!$omp parallel do schedule(runtime) private(eID)
            do id = 1, SIZE( mesh % HO_Elements )
               ! Runge-Kutta 3
               eID = mesh % HO_Elements(id)
#ifdef FLOW
               mesh % elements(eID) % storage % G_NS = a(k)* mesh % elements(eID) % storage % G_NS  +              mesh % elements(eID) % storage % QDot
               mesh % elements(eID) % storage % Q =       mesh % elements(eID) % storage % Q  + c(k)*dt_vec(eID)* mesh % elements(eID) % storage % G_NS
#endif

#if (defined(CAHNHILLIARD)) && (!defined(FLOW))
               mesh % elements(eID) % storage % G_CH = a(k)*mesh % elements(eID) % storage % G_CH + mesh % elements(eID) % storage % cDot
               mesh % elements(eID) % storage % c    = mesh % elements(eID) % storage % c         + c(k)*dt_vec(eID)* mesh % elements(eID) % storage % G_CH
#endif
            end do ! id
!$omp end parallel do

         end do ! k

      else !not present dt_vec

         do k = 1,3
            tk = t + b(k)*deltaT
            if ((k==1) .and. (mod(iter, interval) == 0)) then
               call ComputeTimeDerivative( mesh, particles, tk, CTD_IGNORE_MODE)
            else
               call ComputeTimeDerivative( mesh, particles, tk, CTD_IGNORE_MODE, .true.)
            endif
            if ( present(dts) ) then
               if (dts) call ComputePseudoTimeDerivative(mesh, tk, global_dt)
            end if

            if (k==1) then
!$omp parallel do schedule(runtime) private(eID)
               do id = 1, SIZE( mesh % LO_Elements )
                  ! Explicit Euler
                  eID = mesh % LO_Elements(id)
#ifdef FLOW
                  mesh % elements(eID) % storage % Q = mesh % elements(eID) % storage % Q + deltaT*mesh % elements(eID) % storage % QDot
#endif

#if (defined(CAHNHILLIARD)) && (!defined(FLOW))
                  mesh % elements(eID) % storage % c = mesh % elements(eID) % storage % c + deltaT*mesh % elements(eID) % storage % cDot
#endif
               end do ! id
!$omp end parallel do
            end if

!$omp parallel do schedule(runtime) private(eID)
            do id = 1, SIZE( mesh % HO_Elements )
               ! Runge-Kutta 3
               eID = mesh % HO_Elements(id)
#ifdef FLOW
               mesh % elements(eID) % storage % G_NS = a(k)* mesh % elements(eID) % storage % G_NS  +              mesh % elements(eID) % storage % QDot
               mesh % elements(eID) % storage % Q =       mesh % elements(eID) % storage % Q  + c(k)*deltaT* mesh % elements(eID) % storage % G_NS
#endif

#if (defined(CAHNHILLIARD)) && (!defined(FLOW))
               mesh % elements(eID) % storage % G_CH = a(k)*mesh % elements(eID) % storage % G_CH + mesh % elements(eID) % storage % cDot
               mesh % elements(eID) % storage % c    = mesh % elements(eID) % storage % c         + c(k)*deltaT* mesh % elements(eID) % storage % G_CH
#endif
            end do ! id
!$omp end parallel do

         end do ! k

      end if
!
!     To obtain the updated residuals
      if ( CTD_AFTER_STEPS ) CALL ComputeTimeDerivative( mesh, particles, t+deltaT, CTD_IGNORE_MODE)

      call checkForNan(mesh, t)

   END SUBROUTINE TakeEulerRK3Step

!  ------------------------------
!  Routine for taking a RK3 step.
!  ------------------------------
   SUBROUTINE TakeRK3Step( mesh, particles, t, deltaT, ComputeTimeDerivative, dt_vec, dts, global_dt, iter)
!
!     ----------------------------------
!     Williamson's 3rd order Runge-Kutta
!     ----------------------------------
!
      IMPLICIT NONE
!
!     -----------------
!     Input parameters:
!     -----------------
!
      type(HexMesh)      :: mesh
#ifdef FLOW
      type(Particles_t)  :: particles
#else
      logical            :: particles
#endif
      REAL(KIND=RP)   :: t, deltaT, tk
      real(kind=RP), allocatable, dimension(:), intent(in), optional :: dt_vec
      procedure(ComputeTimeDerivative_f)    :: ComputeTimeDerivative
      logical, intent(in), optional :: dts
      real(kind=RP), intent(in), optional :: global_dt
      integer, intent(in), optional :: iter
!
!     ---------------
!     Local variables
!     ---------------
!
      REAL(KIND=RP), DIMENSION(3) :: a = (/0.0_RP       , -5.0_RP /9.0_RP , -153.0_RP/128.0_RP/)
      REAL(KIND=RP), DIMENSION(3) :: b = (/0.0_RP       ,  1.0_RP /3.0_RP ,    3.0_RP/4.0_RP  /)
      REAL(KIND=RP), DIMENSION(3) :: c = (/1.0_RP/3.0_RP,  15.0_RP/16.0_RP,    8.0_RP/15.0_RP /)

      INTEGER :: i, j, k, m, id, stage

      if (present(dt_vec)) then

         do stage = 1,3
            tk = t + b(stage)*deltaT
            call ComputeTimeDerivative( mesh, particles, tk, CTD_IGNORE_MODE)
            if ( present(dts) ) then
               if (dts) call ComputePseudoTimeDerivative(mesh, tk, global_dt)
            end if

!$omp parallel do schedule(runtime)
            do id = 1, SIZE( mesh % elements )
#ifdef FLOW
                  mesh % elements(id) % storage % G_NS = a(stage)* mesh % elements(id) % storage % G_NS  +              mesh % elements(id) % storage % QDot
                  mesh % elements(id) % storage % Q =       mesh % elements(id) % storage % Q  + c(stage)*dt_vec(id)* mesh % elements(id) % storage % G_NS
#endif
#if (defined(CAHNHILLIARD)) && (!defined(FLOW))
                  mesh % elements(id) % storage % G_CH = a(stage)*mesh % elements(id) % storage % G_CH + mesh % elements(id) % storage % cDot
                  mesh % elements(id) % storage % c    = mesh % elements(id) % storage % c         + c(stage)*dt_vec(id)* mesh % elements(id) % storage % G_CH
#endif
            end do ! id
!$omp end parallel do

         end do ! stage

      else

         !$acc data copyin(deltaT, a, b, c)

         do stage = 1,3
            tk = t + b(stage)*deltaT
            call ComputeTimeDerivative( mesh, particles, tk, CTD_IGNORE_MODE)
            if ( present(dts) ) then
               if (dts) call ComputePseudoTimeDerivative(mesh, tk, global_dt)
            end if

            !$acc parallel loop gang present(mesh,a,b,c,deltaT)
            do id = 1, SIZE( mesh % elements )
               !$acc loop vector collapse(3)
               do k = 0, mesh % elements(id) % Nxyz(3)
                   do j = 0, mesh % elements(id) % Nxyz(2)
                      do i = 0, mesh % elements(id) % Nxyz(1)
#ifdef FLOW
                        !$acc loop seq
                        do m = 1, NCONS
                           mesh % elements(id) % storage % G_NS(m,i,j,k) = a(stage)* mesh % elements(id) % storage % G_NS(m,i,j,k)  +                mesh % elements(id) % storage % QDot(m,i,j,k)
                           mesh % elements(id) % storage % Q(m,i,j,k)    =           mesh % elements(id) % storage % Q(m,i,j,k)     + c(stage)*deltaT* mesh % elements(id) % storage % G_NS(m,i,j,k)
                        enddo
#endif
#if (defined(CAHNHILLIARD)) && (!defined(FLOW))
                        !$acc loop seq
                        do m = 1, NCOMP
                           mesh % elements(id) % storage % G_CH(m,i,j,k) = a(stage)* mesh % elements(id) % storage % G_CH(m,i,j,k) + mesh % elements(id) % storage % cDot(m,i,j,k)
                           mesh % elements(id) % storage % c(m,i,j,k)    =           mesh % elements(id) % storage % c(m,i,j,k)    + c(stage)*deltaT* mesh % elements(id) % storage % G_CH(m,i,j,k)
                        end do
#endif
                     end do
                  end do
               end do
            end do ! id
            !$acc end parallel loop

         end do ! stage

         !$acc end data

      end if
!
!     To obtain the updated residuals
      if ( CTD_AFTER_STEPS ) CALL ComputeTimeDerivative( mesh, particles, t+deltaT, CTD_IGNORE_MODE)

      call checkForNan(mesh, t)

   END SUBROUTINE TakeRK3Step

!  ------------------------------
!  Routine for taking a RK5 step.   
!  ------------------------------
   SUBROUTINE TakeRK5Step( mesh, particles, t, deltaT, ComputeTimeDerivative , dt_vec, dts, global_dt, iter)
!
!        *****************************************************************************************
!           2N-storage 4th-order, 5-stage RK (Carpenter & Kennedy).
!        *****************************************************************************************
!
      implicit none
      type(HexMesh)                   :: mesh
#ifdef FLOW
      type(Particles_t)  :: particles
#else
      logical            :: particles
#endif
      REAL(KIND=RP)                   :: t, deltaT, tk
      procedure(ComputeTimeDerivative_f)      :: ComputeTimeDerivative
      real(kind=RP), allocatable, dimension(:), intent(in), optional :: dt_vec
      logical, intent(in), optional :: dts
      real(kind=RP), intent(in), optional :: global_dt
      integer, intent(in), optional :: iter
!
!     ---------------
!     Local variables
!     ---------------
!
      integer                    :: id, i, j, k, m, stage
      integer, parameter         :: N_STAGES = 5
      real(kind=RP), parameter  :: a(N_STAGES) = [0.0_RP , -0.4178904745_RP, -1.192151694643_RP ,     -1.697784692471_RP , -1.514183444257_RP ]
      real(kind=RP), parameter  :: b(N_STAGES) = [0.0_RP , 0.1496590219993_RP , 0.3704009573644_RP , 0.6222557631345_RP , 0.9582821306748_RP ]
      real(kind=RP), parameter  :: c(N_STAGES) = [0.1496590219993_RP , 0.3792103129999_RP , 0.8229550293869_RP , 0.6994504559488_RP , 0.1530572479681_RP]


      if (present(dt_vec)) then

      DO stage = 1, N_STAGES

         tk = t + b(stage)*deltaT
         CALL ComputeTimeDerivative( mesh, particles, tk, CTD_IGNORE_MODE)
         if ( present(dts) ) then
            if (dts) call ComputePseudoTimeDerivative(mesh, tk, global_dt)
         end if

!$omp parallel do schedule(runtime)
         DO id = 1, SIZE( mesh % elements )
#ifdef FLOW
             mesh % elements(id) % storage % G_NS = a(stage)* mesh % elements(id) % storage % G_NS  +              mesh % elements(id) % storage % QDot
             mesh % elements(id) % storage % Q =       mesh % elements(id) % storage % Q  + c(stage)*dt_vec(id)* mesh % elements(id) % storage % G_NS
#endif
#if (defined(CAHNHILLIARD)) && (!defined(FLOW))
            mesh % elements(id) % storage % G_CH = a(stage)*mesh % elements(id) % storage % G_CH + mesh % elements(id) % storage % cDot
            mesh % elements(id) % storage % c    = mesh % elements(id) % storage % c         + c(stage)*dt_vec(id)* mesh % elements(id) % storage % G_CH
#endif
         END DO
!$omp end parallel do

      END DO

      else

         !$acc data copyin(deltaT, a, b, c)

         DO stage = 1, N_STAGES

            tk = t + b(stage)*deltaT
            CALL ComputeTimeDerivative( mesh, particles, tk, CTD_IGNORE_MODE)
            if ( present(dts) ) then
               if (dts) call ComputePseudoTimeDerivative(mesh, tk, global_dt)
            end if

            !$acc parallel loop gang present(mesh, a, b, c, deltaT)
            do id = 1, SIZE( mesh % elements )
               !$acc loop vector collapse(3)
               do k = 0, mesh % elements(id) % Nxyz(3)
                  do j = 0, mesh % elements(id) % Nxyz(2)
                     do i = 0, mesh % elements(id) % Nxyz(1)
#ifdef FLOW
                        !$acc loop seq
                        do m = 1, NCONS
                           mesh % elements(id) % storage % G_NS(m,i,j,k) = a(stage)* mesh % elements(id) % storage % G_NS(m,i,j,k)  +                mesh % elements(id) % storage % QDot(m,i,j,k)
                           mesh % elements(id) % storage % Q(m,i,j,k)    =           mesh % elements(id) % storage % Q(m,i,j,k)     + c(stage)*deltaT* mesh % elements(id) % storage % G_NS(m,i,j,k)
                        end do
#endif
#if (defined(CAHNHILLIARD)) && (!defined(FLOW))
                        !$acc loop seq
                        do m = 1, NCOMP
                           mesh % elements(id) % storage % G_CH(m,i,j,k) = a(stage)* mesh % elements(id) % storage % G_CH(m,i,j,k) + mesh % elements(id) % storage % cDot(m,i,j,k)
                           mesh % elements(id) % storage % c(m,i,j,k)    =           mesh % elements(id) % storage % c(m,i,j,k)    + c(stage)*deltaT* mesh % elements(id) % storage % G_CH(m,i,j,k)
                        end do
#endif
                     end do
                  end do
               end do
            end do
            !$acc end parallel loop

         END DO

         !$acc end data

      end if

      call checkForNan(mesh, t)
!
!     To obtain the updated residuals
      if ( CTD_AFTER_STEPS ) CALL ComputeTimeDerivative( mesh, particles, t+deltaT, CTD_IGNORE_MODE)

   end subroutine TakeRK5Step
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!  ------------------------------
!  Routine for taking a SSP-RK 3-stage 3rd-order step.
!  ------------------------------
   subroutine TakeSSPRK33Step( mesh, particles, t, deltaT, ComputeTimeDerivative, dt_vec, dts, global_dt, iter)
      implicit none
!
!     ----------------
!     Input parameters
!     ----------------
!
      type(HexMesh)                               :: mesh
#if defined(FLOW)
      type(Particles_t)                           :: particles
#else
      logical                                     :: particles
#endif
      real(RP)                                    :: t
      real(RP)                                    :: deltaT
      procedure(ComputeTimeDerivative_f)          :: ComputeTimeDerivative
      real(RP), allocatable, optional, intent(in) :: dt_vec(:)
      logical,               optional, intent(in) :: dts
      real(RP),              optional, intent(in) :: global_dt
      integer,               optional, intent(in) :: iter
!
!     ---------------
!     Local variables
!     ---------------
!
      real(RP), parameter :: a(3) = [1.0_RP, 3.0_RP/4.0_RP, 1.0_RP/3.0_RP]
      real(RP), parameter :: b(3) = [0.0_RP, 1.0_RP/4.0_RP, 2.0_RP/3.0_RP]
      real(RP), parameter :: c(3) = [1.0_RP, 1.0_RP/4.0_RP, 2.0_RP/3.0_RP]
      real(RP), parameter :: d(3) = [0.0_RP, 1.0_RP,        0.5_RP]
      real(RP) :: tk
      integer  :: i, j, k, m, id, stage


      if (present(dt_vec)) then

!$omp parallel do
         do id = 1, size(mesh % elements)
#if defined(FLOW)
            mesh % elements(id) % storage % G_NS = mesh % elements(id) % storage % Q
#elif defined(CAHNHILLIARD)
            mesh % elements(id) % storage % G_CH = mesh % elements(id) % storage % c
#endif
         end do
!$omp end parallel do

         do stage = 1, 3
            tk = t + d(stage)*deltaT
            call ComputeTimeDerivative( mesh, particles, tk, CTD_IGNORE_MODE)
            if ( present(dts) ) then
               if (dts) call ComputePseudoTimeDerivative(mesh, tk, global_dt)
            end if

!$omp parallel do
            do id = 1, size( mesh % elements )
#if defined(FLOW)
               mesh % elements(id) % storage % Q = a(stage) * mesh % elements(id) % storage % G_NS &
                                                 + b(stage) * mesh % elements(id) % storage % Q    &
                                                 + c(stage) * dt_vec(id) * mesh % elements(id) % storage % Qdot
#elif defined(CAHNHILLIARD)
               mesh % elements(id) % storage % c = a(stage) * mesh % elements(id) % storage % G_CH &
                                                 + b(stage) * mesh % elements(id) % storage % c    &
                                                 + c(stage) * dt_vec(id) * mesh % elements(id) % storage % cDot
#endif
            end do ! id
!$omp end parallel do

            if (LIMITED) then
               call stage_limiter(mesh)
            end if

         end do ! stage

      else

         !$acc data copyin(deltaT, a, b, c)
!
!        Save initial solution Q^n into G_NS
!        -----------------------------------
         !$acc parallel loop gang present(mesh)
         do id = 1, size(mesh % elements)
            !$acc loop vector collapse(3)
            do k = 0, mesh % elements(id) % Nxyz(3) ; do j = 0, mesh % elements(id) % Nxyz(2) ; do i = 0, mesh % elements(id) % Nxyz(1)
               !$acc loop seq
#if defined(FLOW)
               do m = 1, NCONS
                  mesh % elements(id) % storage % G_NS(m,i,j,k) = mesh % elements(id) % storage % Q(m,i,j,k)
               end do
#elif defined(CAHNHILLIARD)
               do m = 1, NCOMP
                  mesh % elements(id) % storage % G_CH(m,i,j,k) = mesh % elements(id) % storage % c(m,i,j,k)
               end do
#endif
            end do                ; end do                ; end do
         end do
         !$acc end parallel loop

         do stage = 1, 3
            tk = t + d(stage) * deltaT
            call ComputeTimeDerivative( mesh, particles, tk, CTD_IGNORE_MODE)
            if ( present(dts) ) then
               if (dts) call ComputePseudoTimeDerivative(mesh, tk, global_dt)
            end if

            !$acc parallel loop gang present(mesh, a, b, c, deltaT)
            do id = 1, size( mesh % elements )
               !$acc loop vector collapse(3)
               do k = 0, mesh % elements(id) % Nxyz(3) ; do j = 0, mesh % elements(id) % Nxyz(2) ; do i = 0, mesh % elements(id) % Nxyz(1)
                  !$acc loop seq
#if defined(FLOW)
                  do m = 1, NCONS
                     mesh % elements(id) % storage % Q(m,i,j,k) = a(stage) * mesh % elements(id) % storage % G_NS(m,i,j,k) &
                                                                + b(stage) * mesh % elements(id) % storage % Q(m,i,j,k)    &
                                                                + c(stage) * deltaT * mesh % elements(id) % storage % Qdot(m,i,j,k)
                  end do
#elif defined(CAHNHILLIARD)
                  do m = 1, NCOMP
                     mesh % elements(id) % storage % c(m,i,j,k) = a(stage) * mesh % elements(id) % storage % G_CH(m,i,j,k) &
                                                                + b(stage) * mesh % elements(id) % storage % c(m,i,j,k)    &
                                                                + c(stage) * deltaT * mesh % elements(id) % storage % cDot(m,i,j,k)
                  end do
#endif
               end do                ; end do                ; end do
            end do ! id
            !$acc end parallel loop

            if (LIMITED) then
               call stage_limiter(mesh)
            end if

         end do ! stage

         !$acc end data

      end if
!
!     To obtain the updated residuals
      if ( CTD_AFTER_STEPS ) CALL ComputeTimeDerivative( mesh, particles, t+deltaT, CTD_IGNORE_MODE)

      call checkForNan(mesh, t)

   END SUBROUTINE TakeSSPRK33Step
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!  ------------------------------
!  Routine for taking a SSP-RK 4-stage 3rd-order step.
!  ------------------------------
   subroutine TakeSSPRK43Step( mesh, particles, t, deltaT, ComputeTimeDerivative, dt_vec, dts, global_dt, iter)
      implicit none
!
!     ----------------
!     Input parameters
!     ----------------
!
      type(HexMesh)                               :: mesh
#if defined(FLOW)
      type(Particles_t)                           :: particles
#else
      logical                                     :: particles
#endif
      real(RP)                                    :: t
      real(RP)                                    :: deltaT
      procedure(ComputeTimeDerivative_f)          :: ComputeTimeDerivative
      real(RP), allocatable, optional, intent(in) :: dt_vec(:)
      logical,               optional, intent(in) :: dts
      real(RP),              optional, intent(in) :: global_dt
      integer,               optional, intent(in) :: iter

!
!     ---------------
!     Local variables
!     ---------------
!
      real(RP), parameter :: a(4) = [1.0_RP, 0.0_RP, 2.0_RP/3.0_RP, 0.0_RP]
      real(RP), parameter :: b(4) = [0.0_RP, 1.0_RP, 1.0_RP/3.0_RP, 1.0_RP]
      real(RP), parameter :: c(4) = [0.5_RP, 0.5_RP, 1.0_RP/6.0_RP, 0.5_RP]
      real(RP), parameter :: d(4) = [0.0_RP, 0.5_RP, 1.0_RP,        0.5_RP]
      real(RP) :: tk
      integer  :: i, j, k, m, id, stage


      if (present(dt_vec)) then

!$omp parallel do
         do id = 1, size(mesh % elements)
#if defined(FLOW)
            mesh % elements(id) % storage % G_NS = mesh % elements(id) % storage % Q
#elif defined(CAHNHILLIARD)
            mesh % elements(id) % storage % G_CH = mesh % elements(id) % storage % c
#endif
         end do
!$omp end parallel do

         do stage = 1, 4
            tk = t + d(stage)*deltaT
            call ComputeTimeDerivative( mesh, particles, tk, CTD_IGNORE_MODE)
            if ( present(dts) ) then
               if (dts) call ComputePseudoTimeDerivative(mesh, tk, global_dt)
            end if

!$omp parallel do
            do id = 1, size( mesh % elements )
#if defined(FLOW)
               mesh % elements(id) % storage % Q = a(stage) * mesh % elements(id) % storage % G_NS &
                                                 + b(stage) * mesh % elements(id) % storage % Q    &
                                                 + c(stage) * dt_vec(id) * mesh % elements(id) % storage % Qdot
#elif defined(CAHNHILLIARD)
               mesh % elements(id) % storage % c = a(stage) * mesh % elements(id) % storage % G_CH &
                                                 + b(stage) * mesh % elements(id) % storage % c    &
                                                 + c(stage) * dt_vec(id) * mesh % elements(id) % storage % cDot
#endif
            end do ! id
!$omp end parallel do

            if (LIMITED) then
               call stage_limiter(mesh)
            end if

         end do ! stage

      else

         !$acc data copyin(deltaT, a, b, c)
!
!        Save initial solution Q^n into G_NS
!        -----------------------------------
         !$acc parallel loop gang present(mesh)
         do id = 1, size(mesh % elements)
            !$acc loop vector collapse(3)
            do k = 0, mesh % elements(id) % Nxyz(3) ; do j = 0, mesh % elements(id) % Nxyz(2) ; do i = 0, mesh % elements(id) % Nxyz(1)
               !$acc loop seq
#if defined(FLOW)
               do m = 1, NCONS
                  mesh % elements(id) % storage % G_NS(m,i,j,k) = mesh % elements(id) % storage % Q(m,i,j,k)
               end do
#elif defined(CAHNHILLIARD)
               do m = 1, NCOMP
                  mesh % elements(id) % storage % G_CH(m,i,j,k) = mesh % elements(id) % storage % c(m,i,j,k)
               end do
#endif
            end do                ; end do                ; end do
         end do
         !$acc end parallel loop

         do stage = 1, 4
            tk = t + d(stage) * deltaT
            call ComputeTimeDerivative( mesh, particles, tk, CTD_IGNORE_MODE)
            if ( present(dts) ) then
               if (dts) call ComputePseudoTimeDerivative(mesh, tk, global_dt)
            end if

            !$acc parallel loop gang present(mesh, a, b, c, deltaT)
            do id = 1, size( mesh % elements )
               !$acc loop vector collapse(3)
               do k = 0, mesh % elements(id) % Nxyz(3) ; do j = 0, mesh % elements(id) % Nxyz(2) ; do i = 0, mesh % elements(id) % Nxyz(1)
                  !$acc loop seq
#if defined(FLOW)
                  do m = 1, NCONS
                     mesh % elements(id) % storage % Q(m,i,j,k) = a(stage) * mesh % elements(id) % storage % G_NS(m,i,j,k) &
                                                                + b(stage) * mesh % elements(id) % storage % Q(m,i,j,k)    &
                                                                + c(stage) * deltaT * mesh % elements(id) % storage % Qdot(m,i,j,k)
                  end do
#elif defined(CAHNHILLIARD)
                  do m = 1, NCOMP
                     mesh % elements(id) % storage % c(m,i,j,k) = a(stage) * mesh % elements(id) % storage % G_CH(m,i,j,k) &
                                                                + b(stage) * mesh % elements(id) % storage % c(m,i,j,k)    &
                                                                + c(stage) * deltaT * mesh % elements(id) % storage % cDot(m,i,j,k)
                  end do
#endif
               end do                ; end do                ; end do
            end do ! id
            !$acc end parallel loop

            if (LIMITED) then
               call stage_limiter(mesh)
            end if

         end do ! stage

         !$acc end data

      end if
!
!     To obtain the updated residuals
      if ( CTD_AFTER_STEPS ) CALL ComputeTimeDerivative( mesh, particles, t+deltaT, CTD_IGNORE_MODE)

      call checkForNan(mesh, t)

   END SUBROUTINE TakeSSPRK43Step
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!  ------------------------------
!  Explicit forward Euler step.
!  ------------------------------
   SUBROUTINE TakeExplicitEulerStep( mesh, particles, t, deltaT, ComputeTimeDerivative , dt_vec, dts, global_dt, iter)
      implicit none
      type(HexMesh)                   :: mesh
#ifdef FLOW
      type(Particles_t)  :: particles
#else
      logical            :: particles
#endif
      REAL(KIND=RP)                   :: t, deltaT, tk
      procedure(ComputeTimeDerivative_f)      :: ComputeTimeDerivative
      real(kind=RP), allocatable, dimension(:), intent(in), optional :: dt_vec
      logical, intent(in), optional :: dts
      real(kind=RP), intent(in), optional :: global_dt
      integer, intent(in), optional :: iter
!
!     ---------------
!     Local variables
!     ---------------
!
      integer :: id, i, j, k, m

      tk = t

      CALL ComputeTimeDerivative( mesh, particles, t, CTD_IGNORE_MODE)
      if ( present(dts) ) then
         if (dts) call ComputePseudoTimeDerivative(mesh, t, global_dt)
      end if

      if (present(dt_vec)) then
!$omp parallel do schedule(runtime)
         DO id = 1, SIZE( mesh % elements )
#ifdef FLOW
            mesh % elements(id) % storage % Q = mesh % elements(id) % storage % Q + dt_vec(id)*mesh % elements(id) % storage % QDot
#endif
#if (defined(CAHNHILLIARD)) && (!defined(FLOW))
            mesh % elements(id) % storage % c = mesh % elements(id) % storage % c + dt_vec(id)*mesh % elements(id) % storage % cDot
#endif
         END DO
!$omp end parallel do
      else

         !$acc data copyin(deltaT)
         !$acc parallel loop gang present(mesh, deltaT)
         do id = 1, SIZE( mesh % elements )
            !$acc loop vector collapse(3)
            do k = 0, mesh % elements(id) % Nxyz(3)
               do j = 0, mesh % elements(id) % Nxyz(2)
                  do i = 0, mesh % elements(id) % Nxyz(1)
#ifdef FLOW
                     !$acc loop seq
                     do m = 1, NCONS
                        mesh % elements(id) % storage % Q(m,i,j,k) =   mesh % elements(id) % storage % Q(m,i,j,k) &
                                                                     + deltaT * mesh % elements(id) % storage % QDot(m,i,j,k)
                     end do
#endif
#if (defined(CAHNHILLIARD)) && (!defined(FLOW))
                     !$acc loop seq
                     do m = 1, NCOMP
                        mesh % elements(id) % storage % c(m,i,j,k) =   mesh % elements(id) % storage % c(m,i,j,k) &
                                                                     + deltaT * mesh % elements(id) % storage % cDot(m,i,j,k)
                     end do
#endif
                  end do
               end do
            end do
         end do
         !$acc end parallel loop
         !$acc end data

      end if

      call checkForNan(mesh, t)
!
!     To obtain the updated residuals
      if ( CTD_AFTER_STEPS ) CALL ComputeTimeDerivative( mesh, particles, t+deltaT, CTD_IGNORE_MODE)

   end subroutine TakeExplicitEulerStep

!  ------------------------------
!  Explicit BDF step.
!  ------------------------------
   subroutine TakeExplicitBDFStep(mesh, particles, t, deltaT, ComputeTimeDerivative)
      implicit none
      type(HexMesh)                      :: mesh
#ifdef FLOW
      type(Particles_t)                  :: particles
#else
      logical                            :: particles
#endif
      REAL(KIND=RP)                      :: t, deltaT, tk
      procedure(ComputeTimeDerivative_f) :: ComputeTimeDerivative
!
!     ---------------
!     Local variables
!     ---------------
!
      integer                  :: id
      real(kind=RP), parameter :: invGamma1 = 1.0_RP, invGamma2 = 2.0_RP/3.0_RP, invGamma3 = 6.0_RP / 11.0_RP
      logical, save            :: isFirst = .true., isSecond = .false., isThird = .false.

      if (isThird) then
!
!        Perform the third order stages
!        ------------------------------
!$omp parallel do schedule(runtime)
         do id = 1, size(mesh % elements)
!           Set y^{*,n+1} in Q and downgrade y^n and y^{n-1}
!           ------------------------------------------------
            mesh % elements(id) % storage % QDot = mesh % elements(id) % storage % prevQ(2) % Q
            mesh % elements(id) % storage % prevQ(2) % Q = mesh % elements(id) % storage % prevQ(1) % Q
            mesh % elements(id) % storage % prevQ(1) % Q = mesh % elements(id) % storage % Q
            mesh % elements(id) % storage % Q =   3.0_RP * mesh % elements(id) % storage % prevQ(1) % Q &
                                                  - 3.0_RP * mesh % elements(id) % storage % prevQ(2) % Q &
                                                  + mesh % elements(id) % storage % QDot
         end do
!$omp end parallel do
!
!        Compute QDot
!        ------------
         CALL ComputeTimeDerivative( mesh, particles, tk, CTD_IGNORE_MODE)
!
!        Perform the time-step
!        ---------------------
!$omp parallel do schedule(runtime)
         do id = 1, size(mesh % elements)
            mesh % elements(id) % storage % Q =   2.0_RP * mesh % elements(id) % storage % prevQ(1) % Q &
                                                - 0.5_RP * mesh % elements(id) % storage % prevQ(2) % Q &
                                                + (1.0_RP/3.0_RP) * mesh % elements(id) % storage % Q &
                                                + deltaT * mesh % elements(id) % storage % QDot
            mesh % elements(id) % storage % Q = invGamma3 * mesh % elements(id) % storage % Q
         end do
!$omp end parallel do


      elseif (isSecond) then
!
!        Perform the second order stages
!        -------------------------------
         if (eBDF_ORDER > 2) then
!$omp parallel do schedule(runtime)
            do id = 1, size(mesh % elements)
!
!              Set for the previous solution
!              -----------------------------
               mesh % elements(id) % storage % prevQ(2) % Q = mesh % elements(id) % storage % prevQ(1) % Q
            end do
!$omp end parallel do
         end if

!$omp parallel do schedule(runtime)
         do id = 1, size(mesh % elements)
!
!           Set y^{*,n+1} in Q
!           ------------------
            mesh % elements(id) % storage % Q = 2.0_RP*mesh % elements(id) % storage % Q -mesh % elements(id) % storage % prevQ(1) % Q
!
!           Set y^{n} in prevQ
!           --------------------
            mesh % elements(id) % storage % prevQ(1) % Q = 0.5_RP*(mesh % elements(id) % storage % Q + mesh % elements(id) % storage % prevQ(1) % Q)

         end do
!$omp end parallel do
!
!        Compute QDot
!        ------------
         CALL ComputeTimeDerivative( mesh, particles, tk, CTD_IGNORE_MODE)
!
!        Perform the time-step
!        ---------------------
!$omp parallel do schedule(runtime)
         do id = 1, size(mesh % elements)
            mesh % elements(id) % storage % Q =   mesh % elements(id) % storage % prevQ(1) % Q + 0.5_RP*mesh % elements(id) % storage % Q  &
                                                + deltaT * mesh % elements(id) % storage % QDot
            mesh % elements(id) % storage % Q = invGamma2 * mesh % elements(id) % storage % Q

         end do
!$omp end parallel do

         if (eBDF_ORDER > 1 ) then
!
!           Move to third order
!           -------------------
            isFirst = .false.
            isSecond = .false.
            isThird = .true.
         end if

      elseif (isFirst) then
!
!        Perform the first order stages
!        ------------------------------
         CALL ComputeTimeDerivative( mesh, particles, tk, CTD_IGNORE_MODE)

!$omp parallel do schedule(runtime)
         DO id = 1, SIZE( mesh % elements )
            if (eBDF_ORDER > 1 ) then
!
!              Set for the previous solution
!              -----------------------------
               mesh % elements(id) % storage % prevQ(1) % Q = mesh % elements(id) % storage % Q

            end if

            mesh % elements(id) % storage % Q = mesh % elements(id) % storage % Q + deltaT*mesh % elements(id) % storage % QDot

         END DO
!$omp end parallel do

            if (eBDF_ORDER > 1 ) then
!
!              Move to second order
!              --------------------
               isFirst = .false.
               isSecond = .true.
            end if

      end if

   end subroutine TakeExplicitBDFStep

   subroutine Enable_CTD_AFTER_STEPS()
      implicit none
      CTD_AFTER_STEPS = .true.
   end subroutine Enable_CTD_AFTER_STEPS
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!  ------------------------------
!  Optimal RK (Bassi 2009).
!  ------------------------------
   SUBROUTINE TakeRKOptStep( mesh, particles, t, deltaT, ComputeTimeDerivative , N_STAGES, dt_vec, dts, global_dt, iter)
!
!        *****************************************************************************************
!       Optimal RK coefficients from Bassi2009
!        *****************************************************************************************
!
      implicit none
      type(HexMesh)                   :: mesh
#ifdef FLOW
      type(Particles_t)  :: particles
#else
      logical            :: particles
#endif
      REAL(KIND=RP)                   :: t, deltaT, tk
      procedure(ComputeTimeDerivative_f)      :: ComputeTimeDerivative
      integer, intent(in)         :: N_STAGES
      real(kind=RP), allocatable, dimension(:), intent(in), optional :: dt_vec
      logical, intent(in), optional :: dts
      real(kind=RP), intent(in), optional :: global_dt
      integer, intent(in), optional :: iter
!
!     ---------------
!     Local variables
!     ---------------
!
      integer                    :: id, i, j, k, m, stage
      real(kind=RP), dimension(6,7) :: Am, Bm
      real(kind=RP) :: a(N_STAGES), b(N_STAGES)

      Bm(1,:) = (/ 0.9998596842476678, 0.5, 0.0, 0.0, 0.0, 0.0, 0.0 /)
      Bm(2,:) = (/ 0.528003175664866, 0.5193233361621609, 0.3209132144853066, 0.0, 0.0, 0.0, 0.0 /)
      Bm(3,:) = (/ 0.4062766523561424, 0.3590274668186006, 0.2782786562184366, 0.3031000737788218, 0.0, 0.0, 0.0 /)
      Bm(4,:) = (/ 0.32451232607547, 0.2850381110111294, 0.2299189950459751, 0.3245118697485674, 0.1925289659886354, 0.0, 0.0 /)
      Bm(5,:) = (/ 0.2712469842000987, 0.2506886071464794, 0.1571621623659122, 0.2281484031198761, 0.2523511867383585, 0.1918765315676819, 0.0 /)
      Bm(6,:) = (/ 0.2328811281838825, 0.2007437175473038, 0.1577268986576998, 0.2052094755795549, 0.2138853585222901, 0.192146045128098, 0.1428487872093191 /)

      Am(1,:) = (/ 0.0 / 1.0, -0.4998596842476678 / 0.5, 0.0, 0.0, 0.0, 0.0, 0.0 /)
      Am(2,:) = (/ 0.0 / 1.0, -0.1645541151603977 / 0.315637725010225, -0.20368561115193584 / 0.3209132144853066, 0.0, 0.0, 0.0, 0.0 /)
      Am(3,:) = (/ 0.0 / 1.0, -0.11076800268308828 / 0.1937832814991212, -0.1652441853194794 / 0.207607995049003, &
            -0.0706706611694336 / 0.3031000737788218, 0.0, 0.0, 0.0 /)
      Am(4,:) = (/ 0.0 / 1.0, -0.11517541944711737 / 0.1896693024156952, -0.0953688085954342 / 0.2159361297115306, &
            -0.01398286533444451 / 0.1925286952557861, -0.1319831744927813 / 0.1925289659886354, 0.0, 0.0 /)
      Am(5,:) = (/ 0.0 / 1.0, -0.056511008554826186 / 0.1229547684000589, -0.1277338387464205 / 0.1094339401033201, &
            -0.047728222262592115 / 0.1824888900957738, -0.045659513024102316 / 0.1785098941878928,-0.0738412925504657 / 0.1918765315676819, 0.0 /)
      Am(6,:) = (/ 0.0 / 1.0, -0.044038551801784204 / 0.1267393084745009, -0.0740044090728029 / 0.1061061571001407, &
            -0.051620741557559094 / 0.1525740388611983, -0.052635436718356604 / 0.1650271578073516,-0.04885820071493849 / 0.1178619741653911, &
             -0.0742840709627069 / 0.1428487872093191 /)

      a = Am(N_STAGES-1,1:N_STAGES)
      b = Bm(N_STAGES-1,1:N_STAGES)

      tk = t + deltaT

      if (present(dt_vec)) then

      DO stage = 1, N_STAGES

         CALL ComputeTimeDerivative( mesh, particles, tk, CTD_IGNORE_MODE)
         if ( present(dts) ) then
            if (dts) call ComputePseudoTimeDerivative(mesh, tk, global_dt)
         end if

!$omp parallel do schedule(runtime)
         DO id = 1, SIZE( mesh % elements )
#ifdef FLOW
             mesh % elements(id) % storage % G_NS = a(stage)* mesh % elements(id) % storage % G_NS  +  dt_vec(id) * mesh % elements(id) % storage % QDot
             mesh % elements(id) % storage % Q =       mesh % elements(id) % storage % Q   +  b(stage) * mesh % elements(id) % storage % G_NS
#endif
#if (defined(CAHNHILLIARD)) && (!defined(FLOW))
            mesh % elements(id) % storage % G_CH = a(stage)*mesh % elements(id) % storage % G_CH  +  dt_vec(id) * mesh % elements(id) % storage % cDot
            mesh % elements(id) % storage % c    = mesh % elements(id) % storage % c          +  b(stage) * mesh % elements(id) % storage % G_CH
#endif
         END DO
!$omp end parallel do

      END DO

      else

         !$acc data copyin(deltaT, a, b)

         DO stage = 1, N_STAGES

            CALL ComputeTimeDerivative( mesh, particles, tk, CTD_IGNORE_MODE)
            if ( present(dts) ) then
               if (dts) call ComputePseudoTimeDerivative(mesh, tk, global_dt)
            end if

            !$acc parallel loop gang present(mesh, a, b, deltaT)
            do id = 1, SIZE( mesh % elements )
               !$acc loop vector collapse(3)
               do k = 0, mesh % elements(id) % Nxyz(3)
                  do j = 0, mesh % elements(id) % Nxyz(2)
                     do i = 0, mesh % elements(id) % Nxyz(1)
#ifdef FLOW
                        !$acc loop seq
                        do m = 1, NCONS
                           mesh % elements(id) % storage % G_NS(m,i,j,k) = a(stage)* mesh % elements(id) % storage % G_NS(m,i,j,k)  + deltaT * mesh % elements(id) % storage % QDot(m,i,j,k)
                           mesh % elements(id) % storage % Q(m,i,j,k)    =           mesh % elements(id) % storage % Q(m,i,j,k)     + b(stage)* mesh % elements(id) % storage % G_NS(m,i,j,k)
                        end do
#endif
#if (defined(CAHNHILLIARD)) && (!defined(FLOW))
                        !$acc loop seq
                        do m = 1, NCOMP
                           mesh % elements(id) % storage % G_CH(m,i,j,k) = a(stage)* mesh % elements(id) % storage % G_CH(m,i,j,k) + deltaT * mesh % elements(id) % storage % cDot(m,i,j,k)
                           mesh % elements(id) % storage % c(m,i,j,k)    =           mesh % elements(id) % storage % c(m,i,j,k)    + b(stage)* mesh % elements(id) % storage % G_CH(m,i,j,k)
                        end do
#endif
                     end do
                  end do
               end do
            end do
            !$acc end parallel loop

         END DO

         !$acc end data

      end if

      call checkForNan(mesh, t)
!
!     To obtain the updated residuals
      if ( CTD_AFTER_STEPS ) CALL ComputeTimeDerivative( mesh, particles, tk, CTD_IGNORE_MODE)

   end subroutine TakeRKOptStep
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   subroutine Enable_limiter(integrator, minimum)
!
!     ---------
!     Interface
!     ---------
      integer,            intent(in) :: integrator
      real(RP), optional, intent(in) :: minimum


      LIMITED = (integrator == SSPRK33_KEY) .or. (integrator == SSPRK43_KEY)

      if (present(minimum)) then
         LIMITER_MIN = minimum
      end if

   end subroutine Enable_limiter

#if defined(NAVIERSTOKES) || defined(SPALARTALMARAS)
   subroutine stage_limiter(mesh)
!
!     -------
!     Modules
!     -------
!
      use ElementClass,      only: Element
      use NodalStorageClass, only: NodalStorage
      use FluidData,         only: thermodynamics
!
!     ---------
!     Interface
!     ---------
!
      type(HexMesh), target, intent(inout) :: mesh
!
!     ---------------
!     Local variables
!     ---------------
      real(RP)               :: m
      real(RP)               :: rho, minrho, rho_limited
      real(RP)               :: q2, q3, q4, q5
      real(RP)               :: p, pavg, minp
      real(RP)               :: theta_rho, theta_p
      real(RP)               :: gm1, lim_min
      real(RP)               :: wJ
      real(RP)               :: qavg(5)
      real(RP)               :: q1sum, q2sum, q3sum, q4sum, q5sum
      real(RP)               :: vol
      integer                :: eID
      integer                :: i, j, k, mEq
      integer                :: Nx, Ny, Nz
      type(Element), pointer :: e
      real(RP), pointer, contiguous :: Q(:,:,:,:)

      logical                :: do_rho_limit, do_p_limit

      gm1     = thermodynamics % gammaMinus1
      lim_min = LIMITER_MIN

      !$acc parallel loop gang present(mesh)                            &
      !$acc& firstprivate(gm1, lim_min)                                 &
      !$acc& private( Nx, Ny, Nz, vol,                                  &
      !$acc&          q1sum, q2sum, q3sum, q4sum, q5sum, minrho, qavg,  &
      !$acc&          pavg, minp,                                       &
      !$acc&          theta_rho, theta_p, m,                            &
      !$acc&          do_rho_limit, do_p_limit,                         &
      !$acc&          rho, rho_limited, q2, q3, q4, q5, p, wJ, mEq, e, Q )
      do eID = 1, mesh % no_of_elements
         e  => mesh % elements(eID)
         Q => e % storage % Q
         Nx = e % Nxyz(1); Ny = e % Nxyz(2); Nz = e % Nxyz(3)

         vol = e % geom % volume

         ! Pass 1
         q1sum = 0.0_RP; q2sum = 0.0_RP; q3sum = 0.0_RP
         q4sum = 0.0_RP; q5sum = 0.0_RP
         minrho = huge(1.0_RP)

         !$acc loop vector collapse(3)                         &
         !$acc& reduction(+:q1sum, q2sum, q3sum, q4sum, q5sum) &
         !$acc& reduction(min:minrho)                          &
         !$acc& private(rho, wJ)
         do k = 0, Nz; do j = 0, Ny; do i = 0, Nx
            rho = Q(1,i,j,k)
            minrho = min(minrho, rho)
            wJ = NodalStorage(Nx) % w(i) * &
                 NodalStorage(Ny) % w(j) * &
                 NodalStorage(Nz) % w(k) * &
                 e % geom % jacobian(i,j,k)
            q1sum = q1sum + Q(1,i,j,k) * wJ
            q2sum = q2sum + Q(2,i,j,k) * wJ
            q3sum = q3sum + Q(3,i,j,k) * wJ
            q4sum = q4sum + Q(4,i,j,k) * wJ
            q5sum = q5sum + Q(5,i,j,k) * wJ
         end do; end do; end do

         qavg = [q1sum, q2sum, q3sum, q4sum, q5sum] / vol

         theta_rho   = 1.0_RP
         do_rho_limit = .false.

         if (qavg(1) /= minrho) then
            m = min(lim_min, qavg(1))
            theta_rho = abs((qavg(1) - m) / (qavg(1) - minrho))
            if (theta_rho <= 1.0_RP) do_rho_limit = .true.
         end if

         ! Pass 2. Reference is p(qavg), not the mean nodal pressure: p is concave in Q,
         ! so this theta keeps every node >= m = min(LIMITER_MIN, p(qavg)) (Zhang & Shu, 2010)
         pavg = gm1 * (qavg(5) - 0.5_RP * (qavg(2)*qavg(2) + qavg(3)*qavg(3) + qavg(4)*qavg(4)) / qavg(1))
         minp = huge(1.0_RP)

         if (do_rho_limit) then
            !$acc loop vector collapse(3)                &
            !$acc& reduction(min:minp)                   &
            !$acc& private(rho, rho_limited, q2, q3, q4, q5, p)
            do k = 0, Nz; do j = 0, Ny; do i = 0, Nx
               rho = Q(1,i,j,k)
               rho_limited = theta_rho * (rho - qavg(1)) + qavg(1)

               q2 = Q(2,i,j,k)
               q3 = Q(3,i,j,k)
               q4 = Q(4,i,j,k)
               q5 = Q(5,i,j,k)

               p = gm1 * (q5 - 0.5_RP * (q2*q2 + q3*q3 + q4*q4) / rho_limited)

               minp = min(p, minp)
            end do; end do; end do
         else
            !$acc loop vector collapse(3)                &
            !$acc& reduction(min:minp)                   &
            !$acc& private(rho, q2, q3, q4, q5, p)
            do k = 0, Nz ; do j = 0, Ny ; do i = 0, Nx
               rho = Q(1,i,j,k)
               q2  = Q(2,i,j,k)
               q3  = Q(3,i,j,k)
               q4  = Q(4,i,j,k)
               q5  = Q(5,i,j,k)
               p = gm1 * (q5 - 0.5_RP * (q2*q2 + q3*q3 + q4*q4) / rho)
               minp = min(p, minp)
            end do; end do; end do
         end if

         theta_p    = 1.0_RP
         do_p_limit = .false.

         if (pavg /= minp) then
            m = min(lim_min, pavg)
            theta_p = abs((pavg - m) / (pavg - minp))
            if (theta_p <= 1.0_RP) do_p_limit = .true.
         end if

         ! Pass 3
         if (do_rho_limit) then
            if (do_p_limit) then
               !$acc loop seq
               do k = 0, Nz ;do j = 0, Ny; do i = 0, Nx
                  rho = Q(1,i,j,k)
                  rho_limited = theta_rho * (rho - qavg(1)) + qavg(1)
                  Q(1,i,j,k) = theta_p * (rho_limited - qavg(1)) + qavg(1)
                  !$acc loop seq
                  do mEq = 2, 5
                     Q(mEq,i,j,k) = theta_p * (Q(mEq,i,j,k) - qavg(mEq)) + qavg(mEq)
                  end do
               end do; end do; end do
            else
               !$acc loop seq
               do k = 0, Nz ;do j = 0, Ny; do i = 0, Nx
                  rho = Q(1,i,j,k)
                  Q(1,i,j,k) = theta_rho * (rho - qavg(1)) + qavg(1)
               end do; end do; end do;
            end if
         else if (do_p_limit) then
            !$acc loop seq
            do k = 0, Nz ;do j = 0, Ny; do i = 0, Nx
               Q(1,i,j,k) = theta_p * (Q(1,i,j,k) - qavg(1)) + qavg(1)
               !$acc loop seq
               do mEq = 2, 5
                  Q(mEq,i,j,k) = theta_p * (Q(mEq,i,j,k) - qavg(mEq)) + qavg(mEq)
               end do
            end do; end do; end do
         end if
      end do
      !$acc end parallel loop

   end subroutine stage_limiter
#else
   subroutine stage_limiter(mesh)
      type(HexMesh), intent(in) :: mesh
   end subroutine stage_limiter
#endif
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   Subroutine checkForNan(mesh, t)
!
!        **************************************************************************************************************
!        Look if there is a nan in the solution. NOTE: NOT CURRENTLY WORKING  
!        **************************************************************************************************************
!
      use MPI_Process_Info
#ifdef _HAS_MPI_
      use mpi
#endif
      implicit none

      type(HexMesh), intent(in)     :: mesh
      real(kind=RP), intent(in)     :: t

      !local variables
      integer                    :: eID
      CHARACTER(len=LINE_LENGTH) :: FinalName      !  Final name for particular restart file
      logical                    :: NanNotFound, allNan
      integer                    :: ierr

      ! use not found instead of found as OMP reduction initialized the private value as true
      ! this is redundant for the OMP but left for non parallel compilations
      NanNotFound = .TRUE.

!$omp parallel do reduction(.AND.:NanNotFound) schedule(runtime)
      do eID=1, mesh % no_of_elements
         !if ( any(isnan(mesh % elements(eID) % storage % Q))) then
         !   NanNotFound = .FALSE.
         !endif
      end do
!$omp end parallel do

#ifdef _HAS_MPI_
      call mpi_allreduce(NanNotFound, allNan, 1, MPI_LOGICAL, MPI_LOR, MPI_COMM_WORLD, ierr)
      NanNotFound = allNan
#endif

      if (.not. NanNotFound) then
          if ( MPI_Process % isRoot ) print*, "Numerical divergence obtained in solver."
          WRITE(FinalName,'(A,ES11.5,A)')  'RESULTS/horses_divergence_',t,'.hsol'
          if ( MPI_Process % isRoot ) print *, "Writing failed solution: ", FinalName
          call mesh % SaveSolution(0, t, FinalName, .FALSE.)
      end if

#ifdef _HAS_MPI_
     call mpi_barrier(MPI_COMM_WORLD, ierr)
#endif

      if (.not. NanNotFound) call exit(99)


   End Subroutine checkForNan
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
END MODULE ExplicitMethods