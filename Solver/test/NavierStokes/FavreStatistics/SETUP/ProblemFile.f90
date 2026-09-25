!
!////////////////////////////////////////////////////////////////////////
!
!      The Problem File contains user defined procedures
!      that are used to "personalize" i.e. define a specific
!      problem to be solved. These procedures include initial conditions,
!      exact solutions (e.g. for tests), etc. and allow modifications 
!      without having to modify the main code.
!
!      The procedures, *even if empty* that must be defined are
!
!      UserDefinedSetUp
!      UserDefinedInitialCondition(mesh)
!      UserDefinedPeriodicOperation(mesh)
!      UserDefinedFinalize(mesh)
!      UserDefinedTermination
!
!//////////////////////////////////////////////////////////////////////// 
! 
#include "Includes.h"
#if defined(NAVIERSTOKES)
!
!////////////////////////////////////////////////////////////////////////
!
!     Synthetic samples for the statistics monitor regression test.
!
!     The initial condition is sample 1 and UserDefinedPeriodicOperation
!     (called before every time step) sets the next one. With dt = 1e-30
!     the time step leaves Q unchanged (to round-off), so the statistics
!     monitor samples exactly these states: IC (t=0) + one per time step.
!
!     Every node gets the samples shifted by a node-dependent offset
!     (u + x, v - y, w + 2z): the means shift accordingly and the
!     fluctuations are unchanged, so node/element indexing errors are
!     caught as well.
!
!     Case A (control files without "Sine"), 4 samples:
!        rho = (1,2,3,4), u = (4,3,2,1), v = (1,1,3,3), w = (2,1,1,3)
!     Case B ("Sine" control file), 16 samples over one period:
!        rho = 1 + sin(theta)/2, u = 2 + sin(theta), v = w = 0,
!        theta_n = 2*pi*(n-1)/16
!
!     This module only uses types and parameters from the solver
!     libraries: the problem file is a shared library and must not pull
!     in (and duplicate) solver objects with module variables.
!
!////////////////////////////////////////////////////////////////////////
!
module FavreStatisticsSamples
   use SMConstants
   use PhysicsStorage, only: IRHO, IRHOU, IRHOV, IRHOW, IRHOE
   use HexMeshClass,   only: HexMesh
   implicit none

   private
   public OFFSET, N_SIN, favreTestCase, getSamples, setSample

   real(kind=RP), parameter :: OFFSET(NDIM) = [1.0_RP, -1.0_RP, 2.0_RP]
   integer,       parameter :: N_SIN = 16

   contains
!
!     Which test the control file (first command line argument) runs:
!        isSine:         case B instead of case A
!        expectReynolds: .false. for 'averaging = Favre' ("FavreOnly")
!     ---------------------------------------------------------------
      subroutine favreTestCase(isSine, expectReynolds)
         logical, intent(out) :: isSine, expectReynolds
         character(len=LINE_LENGTH) :: controlFile

         call get_command_argument(1, controlFile)
         isSine         = (index(controlFile, "Sine")      .ne. 0)
         expectReynolds = (index(controlFile, "FavreOnly") .eq. 0)
      end subroutine favreTestCase

      subroutine getSamples(isSine, rho, u, v, w)
         logical,                    intent(in)  :: isSine
         real(kind=RP), allocatable, intent(out) :: rho(:), u(:), v(:), w(:)
         real(kind=RP) :: theta(N_SIN)
         integer       :: n

         if (isSine) then
            theta = [(2.0_RP * PI * n / N_SIN, n = 0, N_SIN-1)]
            rho = 1.0_RP + 0.5_RP * sin(theta)
            u   = 2.0_RP + sin(theta)
            v   = [(0.0_RP, n = 1, N_SIN)]
            w   = [(0.0_RP, n = 1, N_SIN)]
         else
            rho = [1.0_RP, 2.0_RP, 3.0_RP, 4.0_RP]
            u   = [4.0_RP, 3.0_RP, 2.0_RP, 1.0_RP]
            v   = [1.0_RP, 1.0_RP, 3.0_RP, 3.0_RP]
            w   = [2.0_RP, 1.0_RP, 1.0_RP, 3.0_RP]
         end if
      end subroutine getSamples
!
!     Set sample s (clipped to the last one) at every node
!     ----------------------------------------------------
      subroutine setSample(mesh, s)
         class(HexMesh)      :: mesh
         integer, intent(in) :: s
         real(kind=RP), allocatable :: rho(:), u(:), v(:), w(:)
         real(kind=RP) :: vel(NDIM)
         logical       :: isSine, expectReynolds
         integer       :: eID, i, j, k, sc

         call favreTestCase(isSine, expectReynolds)
         call getSamples(isSine, rho, u, v, w)
         sc = min(s, size(rho))

         do eID = 1, mesh % no_of_elements
            associate(e => mesh % elements(eID))
            do k = 0, e % Nxyz(3) ; do j = 0, e % Nxyz(2) ; do i = 0, e % Nxyz(1)
               vel = [u(sc), v(sc), w(sc)] + OFFSET * e % geom % x(:,i,j,k)
               e % storage % Q(IRHO,i,j,k)        = rho(sc)
               e % storage % Q(IRHOU:IRHOW,i,j,k) = rho(sc) * vel
               e % storage % Q(IRHOE,i,j,k)       = 2.5_RP + 0.5_RP * rho(sc) * sum(vel**2)
            end do                ; end do                ; end do
            end associate
            !$acc update device(mesh % elements(eID) % storage % Q) if_present
         end do
      end subroutine setSample

end module FavreStatisticsSamples
#endif
module ProblemFileFunctions
   implicit none

   abstract interface
      subroutine UserDefinedStartup_f
      end subroutine UserDefinedStartup_f
   
      SUBROUTINE UserDefinedFinalSetup_f(mesh &
#ifdef FLOW
                                     , thermodynamics_ &
                                     , dimensionless_  &
                                     , refValues_ & 
#endif
#ifdef CAHNHILLIARD
                                     , multiphase_ &
#endif
                                     )
         USE HexMeshClass
         use FluidData
         IMPLICIT NONE
         CLASS(HexMesh)                      :: mesh
#ifdef FLOW
         type(Thermodynamics_t), intent(in)  :: thermodynamics_
         type(Dimensionless_t),  intent(in)  :: dimensionless_
         type(RefValues_t),      intent(in)  :: refValues_
#endif
#ifdef CAHNHILLIARD
         type(Multiphase_t),     intent(in)  :: multiphase_
#endif
      END SUBROUTINE UserDefinedFinalSetup_f

      subroutine UserDefinedInitialCondition_f(mesh &
#ifdef FLOW
                                     , thermodynamics_ &
                                     , dimensionless_  &
                                     , refValues_ & 
#endif
#ifdef CAHNHILLIARD
                                     , multiphase_ &
#endif
                                     )
         use smconstants
         use physicsstorage
         use hexmeshclass
         use fluiddata
         implicit none
         class(hexmesh)                      :: mesh
#ifdef FLOW
         type(Thermodynamics_t), intent(in)  :: thermodynamics_
         type(Dimensionless_t),  intent(in)  :: dimensionless_
         type(RefValues_t),      intent(in)  :: refValues_
#endif
#ifdef CAHNHILLIARD
         type(Multiphase_t),     intent(in)  :: multiphase_
#endif
      end subroutine UserDefinedInitialCondition_f
#ifdef FLOW
      subroutine UserDefinedState_f(x, t, nHat, Q, thermodynamics_, dimensionless_, refValues_)
         use SMConstants
         use PhysicsStorage
         use FluidData
         implicit none
         real(kind=RP)  :: x(NDIM)
         real(kind=RP)  :: t
         real(kind=RP)  :: nHat(NDIM)
         real(kind=RP)  :: Q(NCONS)
         type(Thermodynamics_t), intent(in)  :: thermodynamics_
         type(Dimensionless_t),  intent(in)  :: dimensionless_
         type(RefValues_t),      intent(in)  :: refValues_
      end subroutine UserDefinedState_f

      subroutine UserDefinedGradVars_f(x, t, nHat, Q, U, thermodynamics_, dimensionless_, refValues_)
         use SMConstants
         use PhysicsStorage
         use FluidData
         implicit none
         real(kind=RP), intent(in)          :: x(NDIM)
         real(kind=RP), intent(in)          :: t
         real(kind=RP), intent(in)          :: nHat(NDIM)
         real(kind=RP), intent(in)          :: Q(NCONS)
         real(kind=RP), intent(inout)       :: U(NGRAD)
         type(Thermodynamics_t), intent(in) :: thermodynamics_
         type(Dimensionless_t),  intent(in) :: dimensionless_
         type(RefValues_t),      intent(in) :: refValues_
      end subroutine UserDefinedGradVars_f


      subroutine UserDefinedNeumann_f(x, t, nHat, Q, U_x, U_y, U_z, flux, thermodynamics_, dimensionless_, refValues_)
         use SMConstants
         use PhysicsStorage
         use FluidData
         implicit none
         real(kind=RP), intent(in)    :: x(NDIM)
         real(kind=RP), intent(in)    :: t
         real(kind=RP), intent(in)    :: nHat(NDIM)
         real(kind=RP), intent(in)    :: Q(NCONS)
         real(kind=RP), intent(in)    :: U_x(NGRAD)
         real(kind=RP), intent(in)    :: U_y(NGRAD)
         real(kind=RP), intent(in)    :: U_z(NGRAD)
         real(kind=RP), intent(inout) :: flux(NCONS)
         type(Thermodynamics_t), intent(in) :: thermodynamics_
         type(Dimensionless_t),  intent(in) :: dimensionless_
         type(RefValues_t),      intent(in) :: refValues_
      end subroutine UserDefinedNeumann_f

#endif
!
!//////////////////////////////////////////////////////////////////////// 
! 
      SUBROUTINE UserDefinedPeriodicOperation_f(mesh, time, dt, Monitors)
         use SMConstants
         USE HexMeshClass
         use MonitorsClass
         IMPLICIT NONE
         CLASS(HexMesh)               :: mesh
         REAL(KIND=RP)                :: time
         REAL(KIND=RP)                :: dt
         type(Monitor_t), intent(in) :: monitors
      END SUBROUTINE UserDefinedPeriodicOperation_f
!
!//////////////////////////////////////////////////////////////////////// 
! 
#ifdef FLOW
      subroutine UserDefinedSourceTermNS_f(x, Q, time, S, thermodynamics_, dimensionless_, refValues_ &
#ifdef CAHNHILLIARD
,multiphase_ &
#endif
)
         use SMConstants
         USE HexMeshClass
         use FluidData
         use PhysicsStorage
         IMPLICIT NONE
         real(kind=RP),             intent(in)  :: x(NDIM)
         real(kind=RP),             intent(in)  :: Q(NCONS)
         real(kind=RP),             intent(in)  :: time
         real(kind=RP),             intent(inout) :: S(NCONS)
         type(Thermodynamics_t), intent(in)  :: thermodynamics_
         type(Dimensionless_t),  intent(in)  :: dimensionless_
         type(RefValues_t),      intent(in)  :: refValues_
#ifdef CAHNHILLIARD
         type(Multiphase_t),     intent(in)  :: multiphase_
#endif
      end subroutine UserDefinedSourceTermNS_f
#endif
!
!//////////////////////////////////////////////////////////////////////// 
! 
      SUBROUTINE UserDefinedFinalize_f(mesh, time, iter, maxResidual &
#ifdef FLOW
                                                 , thermodynamics_ &
                                                 , dimensionless_  &
                                                 , refValues_ & 
#endif   
#ifdef CAHNHILLIARD
                                                 , multiphase_ &
#endif
                                                 , monitors, &
                                                   elapsedTime, &
                                                   CPUTime   )
         use SMConstants
         USE HexMeshClass
         use FluidData
         use MonitorsClass
         IMPLICIT NONE
         CLASS(HexMesh)                        :: mesh
         REAL(KIND=RP)                         :: time
         integer                               :: iter
         real(kind=RP)                         :: maxResidual
#ifdef FLOW
         type(Thermodynamics_t), intent(in)    :: thermodynamics_
         type(Dimensionless_t),  intent(in)    :: dimensionless_
         type(RefValues_t),      intent(in)    :: refValues_
#endif
#ifdef CAHNHILLIARD
         type(Multiphase_t),     intent(in)    :: multiphase_
#endif
         type(Monitor_t),        intent(in)    :: monitors
         real(kind=RP),             intent(in) :: elapsedTime
         real(kind=RP),             intent(in) :: CPUTime
      END SUBROUTINE UserDefinedFinalize_f

      SUBROUTINE UserDefinedTermination_f
         implicit none
      END SUBROUTINE UserDefinedTermination_f
   end interface
   
end module ProblemFileFunctions

         SUBROUTINE UserDefinedStartup
!
!        --------------------------------
!        Called before any other routines
!        --------------------------------
!
            IMPLICIT NONE  
         END SUBROUTINE UserDefinedStartup
!
!//////////////////////////////////////////////////////////////////////// 
! 
         SUBROUTINE UserDefinedFinalSetup(mesh &
#ifdef FLOW
                                        , thermodynamics_ &
                                        , dimensionless_  &
                                        , refValues_ & 
#endif
#ifdef CAHNHILLIARD
                                        , multiphase_ &
#endif
                                        )
!
!           ----------------------------------------------------------------------
!           Called after the mesh is read in to allow mesh related initializations
!           or memory allocations.
!           ----------------------------------------------------------------------
!
            USE HexMeshClass
            use PhysicsStorage
            use FluidData
            IMPLICIT NONE
            CLASS(HexMesh)                      :: mesh
#ifdef FLOW
            type(Thermodynamics_t), intent(in)  :: thermodynamics_
            type(Dimensionless_t),  intent(in)  :: dimensionless_
            type(RefValues_t),      intent(in)  :: refValues_
#endif
#ifdef CAHNHILLIARD
            type(Multiphase_t),     intent(in)  :: multiphase_
#endif
         END SUBROUTINE UserDefinedFinalSetup
!
!//////////////////////////////////////////////////////////////////////// 
! 
         subroutine UserDefinedInitialCondition(mesh &
#ifdef FLOW
                                        , thermodynamics_ &
                                        , dimensionless_  &
                                        , refValues_ & 
#endif
#ifdef CAHNHILLIARD
                                        , multiphase_ &
#endif
                                        )
!
!           ------------------------------------------------
!           called to set the initial condition for the flow
!              - by default it sets an uniform initial
!                 condition.
!           ------------------------------------------------
!
            use smconstants
            use physicsstorage
            use hexmeshclass
            use fluiddata
#if defined(NAVIERSTOKES)
            use FavreStatisticsSamples, only: setSample
#endif
            implicit none
            class(hexmesh)                      :: mesh
#ifdef FLOW
            type(Thermodynamics_t), intent(in)  :: thermodynamics_
            type(Dimensionless_t),  intent(in)  :: dimensionless_
            type(RefValues_t),      intent(in)  :: refValues_
#endif
#ifdef CAHNHILLIARD
            type(Multiphase_t),     intent(in)  :: multiphase_
#endif
!
!           ---------------
!           local variables
!           ---------------
!
            integer        :: eid, i, j, k
            real(kind=RP)  :: qq, u, v, w, p
#if defined(NAVIERSTOKES)
            real(kind=RP)  :: Q(NCONS), phi, theta
#endif

!
!           ---------------------------------------
!           Navier-Stokes default initial condition
!           ---------------------------------------
!
#if defined(NAVIERSTOKES)
!           Statistics test: the initial condition is the first sample
            call setSample(mesh, 1)
#endif
!
!           ------------------------------------------------------
!           Incompressible Navier-Stokes default initial condition
!           ------------------------------------------------------
!
#if defined(INCNS)
            do eID = 1, mesh % no_of_elements
               associate( Nx => mesh % elements(eID) % Nxyz(1), &
                          ny => mesh % elemeNts(eID) % nxyz(2), &
                          Nz => mesh % elements(eID) % Nxyz(3) )
               do k = 0, Nz;  do j = 0, Ny;  do i = 0, Nx 
                  mesh % elements(eID) % storage % q(:,i,j,k) = [1.0_RP, 1.0_RP,0.0_RP,0.0_RP,0.0_RP] 
               end do;        end do;        end do
               end associate
            end do
#endif

!
!           ---------------------------------------
!           Cahn-Hilliard default initial condition
!           ---------------------------------------
!
#ifdef CAHNHILLIARD
            call random_seed()
         
            do eid = 1, mesh % no_of_elements
               associate( Nx => mesh % elements(eid) % Nxyz(1), &
                          Ny => mesh % elements(eid) % Nxyz(2), &
                          Nz => mesh % elements(eid) % Nxyz(3) )
               associate(e => mesh % elements(eID) % storage)
               call random_number(e % c) 
               e % c = 2.0_RP * (e % c - 0.5_RP)
               end associate
               end associate
            end do
#endif

         end subroutine UserDefinedInitialCondition
#ifdef FLOW
         subroutine UserDefinedState1(x, t, nHat, Q, thermodynamics_, dimensionless_, refValues_)
            use SMConstants
            use PhysicsStorage
            use FluidData
            implicit none
            real(kind=RP), intent(in)     :: x(NDIM)
            real(kind=RP), intent(in)     :: t
            real(kind=RP), intent(in)     :: nHat(NDIM)
            real(kind=RP), intent(inout)  :: Q(NCONS)
            type(Thermodynamics_t),    intent(in)  :: thermodynamics_
            type(Dimensionless_t),     intent(in)  :: dimensionless_
            type(RefValues_t),         intent(in)  :: refValues_
         end subroutine UserDefinedState1

         subroutine UserDefinedGradVars1(x, t, nHat, Q, U, thermodynamics_, dimensionless_, refValues_)
            use SMConstants
            use PhysicsStorage
            use FluidData
            implicit none
            real(kind=RP), intent(in)          :: x(NDIM)
            real(kind=RP), intent(in)          :: t
            real(kind=RP), intent(in)          :: nHat(NDIM)
            real(kind=RP), intent(in)          :: Q(NCONS)
            real(kind=RP), intent(inout)       :: U(NGRAD)
            type(Thermodynamics_t), intent(in) :: thermodynamics_
            type(Dimensionless_t),  intent(in) :: dimensionless_
            type(RefValues_t),      intent(in) :: refValues_
         end subroutine UserDefinedGradVars1

         subroutine UserDefinedNeumann1(x, t, nHat, Q, U_x, U_y, U_z, flux, thermodynamics_, dimensionless_, refValues_)
            use SMConstants
            use PhysicsStorage
            use FluidData
            implicit none
            real(kind=RP), intent(in)    :: x(NDIM)
            real(kind=RP), intent(in)    :: t
            real(kind=RP), intent(in)    :: nHat(NDIM)
            real(kind=RP), intent(in)    :: Q(NCONS)
            real(kind=RP), intent(in)    :: U_x(NGRAD)
            real(kind=RP), intent(in)    :: U_y(NGRAD)
            real(kind=RP), intent(in)    :: U_z(NGRAD)
            real(kind=RP), intent(inout) :: flux(NCONS)
            type(Thermodynamics_t), intent(in) :: thermodynamics_
            type(Dimensionless_t),  intent(in) :: dimensionless_
            type(RefValues_t),      intent(in) :: refValues_
         end subroutine UserDefinedNeumann1
#endif
!
!//////////////////////////////////////////////////////////////////////// 
! 
         SUBROUTINE UserDefinedPeriodicOperation(mesh, time, dt, Monitors)
!
!           ----------------------------------------------------------
!           Called before every time-step to allow periodic operations
!           to be performed
!           ----------------------------------------------------------
!
            use SMConstants
            USE HexMeshClass
            use MonitorsClass
#if defined(NAVIERSTOKES)
            use FavreStatisticsSamples, only: setSample
#endif
            IMPLICIT NONE
            CLASS(HexMesh)               :: mesh
            REAL(KIND=RP)                :: time
            REAL(KIND=RP)                :: dt
            type(Monitor_t), intent(in) :: monitors
#if defined(NAVIERSTOKES)
!
!           Statistics test: set the next sample before the time step
!           (dt = 1e-30 leaves it unchanged, so the monitor samples it)
!           ------------------------------------------------------------
            call setSample(mesh, monitors % stats % no_of_samples + 1)
#endif

         END SUBROUTINE UserDefinedPeriodicOperation
!
!//////////////////////////////////////////////////////////////////////// 
! 
#ifdef FLOW
         subroutine UserDefinedSourceTermNS(x, Q, time, S, thermodynamics_, dimensionless_, refValues_ &
#ifdef CAHNHILLIARD
, multiphase_ &
#endif
)
!
!           --------------------------------------------
!           Called to apply source terms to the equation
!           --------------------------------------------
!
            use SMConstants
            USE HexMeshClass
            use PhysicsStorage
            use FluidData
            IMPLICIT NONE
            real(kind=RP),             intent(in)  :: x(NDIM)
            real(kind=RP),             intent(in)  :: Q(NCONS)
            real(kind=RP),             intent(in)  :: time
            real(kind=RP),             intent(inout) :: S(NCONS)
            type(Thermodynamics_t), intent(in)  :: thermodynamics_
            type(Dimensionless_t),  intent(in)  :: dimensionless_
            type(RefValues_t),      intent(in)  :: refValues_
#ifdef CAHNHILLIARD
            type(Multiphase_t),     intent(in)  :: multiphase_
#endif
!
!           ---------------
!           Local variables
!           ---------------
!
            integer  :: i, j, k, eID
!
!           Usage example
!           -------------
!           S(:) = x(1) + x(2) + x(3) + time
   
         end subroutine UserDefinedSourceTermNS
#endif
!
!//////////////////////////////////////////////////////////////////////// 
! 
         SUBROUTINE UserDefinedFinalize(mesh, time, iter, maxResidual &
#ifdef FLOW
                                                    , thermodynamics_ &
                                                    , dimensionless_  &
                                                    , refValues_ & 
#endif   
#ifdef CAHNHILLIARD
                                                    , multiphase_ &
#endif
                                                    , monitors, &
                                                      elapsedTime, &
                                                      CPUTime   )
!
!           ---------------------------------------------------------------------
!           Statistics monitor regression test (Reynolds and Favre averaging).
!           The samples fed during the run are described in the
!           FavreStatisticsSamples module (top of this file). Expected values:
!
!           Case A: rho_mean = 2.5
!              Reynolds: u_mean = 2.5, <u'u'> = 1.25, <u'v'> = -1.0
!              Favre:    u~ = 2.0, v~ = 2.4, w~ = 1.9
!                        u''u''~ = 1.0, u''v''~ = -0.8, ...
!           Case B: rho_mean = 1, u_mean = 2, u~ = 2.25, u''u''~ = 0.4375
!
!           Averages are sample averages (not dt-weighted). The number of
!           samples includes the initial condition (t=0): N time steps give
!           N+1 samples.
!           ---------------------------------------------------------------------
!
            use SMConstants
            use FTAssertions
            USE HexMeshClass
            use PhysicsStorage
            use FluidData
            use MonitorsClass
#if defined(NAVIERSTOKES)
            use StatisticsMonitor, only: S_U => U, S_V => V, S_W => W, S_UU => UU, S_VV => VV, &
                                         S_WW => WW, S_UV => UV, S_UW => UW, S_VW => VW, &
                                         NO_OF_VARIABLES_Sij, NO_OF_FAVRE_VARS
            use FavreStatisticsSamples
#endif
            IMPLICIT NONE
            CLASS(HexMesh)                        :: mesh
            REAL(KIND=RP)                         :: time
            integer                               :: iter
            real(kind=RP)                         :: maxResidual
#ifdef FLOW
            type(Thermodynamics_t), intent(in)    :: thermodynamics_
            type(Dimensionless_t),  intent(in)    :: dimensionless_
            type(RefValues_t),      intent(in)    :: refValues_
#endif
#ifdef CAHNHILLIARD
            type(Multiphase_t),     intent(in)    :: multiphase_
#endif
            type(Monitor_t),        intent(in)    :: monitors
            real(kind=RP),             intent(in) :: elapsedTime
            real(kind=RP),             intent(in) :: CPUTime
!
!           ---------------
!           Local variables
!           ---------------
!
            CHARACTER(LEN=29)                  :: testName = "Favre statistics"
            TYPE(FTAssertionsManager), POINTER :: sharedManager
#if defined(NAVIERSTOKES)
            real(kind=RP), parameter           :: TOL = 1.0e-12_RP
            character(len=*), parameter        :: COMP(6) = ["uu","vv","ww","uv","uw","vw"]
            integer,       parameter           :: PAIRS(2,6) = reshape([1,1, 2,2, 3,3, 1,2, 1,3, 2,3], [2,6])
            real(kind=RP), allocatable         :: rho(:), u(:), v(:), w(:)
            logical                            :: isSine, expectReynolds
            character(len=3)                   :: tag
            integer                            :: nSamples, nVars, nR, eID, i, j, k, c
            real(kind=RP)                      :: rhoMean, reyMean(NDIM), reyFluc(6), favMean(NDIM), favFluc(6)
            real(kind=RP)                      :: x(NDIM), rhoBar, reyU(NDIM), reyUU(6), favU(NDIM), favUU(6)
            real(kind=RP)                      :: errRho, errReyMean(NDIM), errReyFluc(6)
            real(kind=RP)                      :: errFavMean(NDIM), errFavFluc(6), errFavStress(6)

            call favreTestCase(isSine, expectReynolds)
            call getSamples(isSine, rho, u, v, w)
            nSamples = size(rho)
!
!           Exact values of the sample set
!           ------------------------------
            if (isSine) then
               tag     = "[B]"
               rhoMean = 1.0_RP
               reyMean = [2.0_RP, 0.0_RP, 0.0_RP]
               reyFluc = [0.5_RP, 0.0_RP, 0.0_RP, 0.0_RP, 0.0_RP, 0.0_RP]
               favMean = [2.25_RP, 0.0_RP, 0.0_RP]
               favFluc = [0.4375_RP, 0.0_RP, 0.0_RP, 0.0_RP, 0.0_RP, 0.0_RP]
            else
               tag     = "[A]"
               rhoMean = 2.5_RP
               reyMean = [2.5_RP, 2.0_RP, 1.75_RP]
               reyFluc = [1.25_RP, 1.0_RP, 0.6875_RP, -1.0_RP, -0.375_RP, 0.25_RP]
               favMean = [2.0_RP, 2.4_RP, 1.9_RP]
               favFluc = [1.0_RP, 0.84_RP, 0.89_RP, -0.8_RP, -0.6_RP, 0.34_RP]
            end if

            CALL initializeSharedAssertionsManager
            sharedManager => sharedAssertionsManager()
!
!           Number of samples: initial condition (t=0) + one per time step
!           --------------------------------------------------------------
            CALL FTAssertEqual(expectedValue = nSamples, actualValue = iter + 1, &
                               msg = tag // " Time steps + 1 = number of samples in the test case")
            CALL FTAssertEqual(expectedValue = iter + 1, actualValue = monitors % stats % no_of_samples, &
                               msg = tag // " Samples taken (N time steps + t=0)")
!
!           Storage layout: [Reynolds (9, optional)] [mean Q (NCONS)] [Favre (6)]
!           ---------------------------------------------------------------------
            nVars = size(mesh % elements(1) % storage % stats % data, 1)
            nR    = nVars - NCONS - NO_OF_FAVRE_VARS
            CALL FTAssertEqual(expectedValue = merge(NO_OF_VARIABLES_Sij, 0, expectReynolds) + NCONS + NO_OF_FAVRE_VARS, &
                               actualValue   = nVars, &
                               msg           = tag // " Number of statistics variables ('averaging' keyword)")

            errRho = 0.0_RP ; errReyMean = 0.0_RP ; errReyFluc = 0.0_RP
            errFavMean = 0.0_RP ; errFavFluc = 0.0_RP ; errFavStress = 0.0_RP

            do eID = 1, mesh % no_of_elements
               !$acc update self(mesh % elements(eID) % storage % stats % data) if_present
               associate(data => mesh % elements(eID) % storage % stats % data)
               do k = 0, mesh % elements(eID) % Nxyz(3) ; do j = 0, mesh % elements(eID) % Nxyz(2) ; do i = 0, mesh % elements(eID) % Nxyz(1)
                  x = OFFSET * mesh % elements(eID) % geom % x(:,i,j,k)
!
!                 Mean conservative variables: <rho>, <rho u_i>
                  rhoBar = data(nR+IRHO,i,j,k)
                  errRho = max(errRho, abs(rhoBar - rhoMean))
!
!                 Reynolds block: <u_i>, <u_i u_j>
                  if (nR .gt. 0) then
                     reyU  = data([S_U, S_V, S_W],i,j,k)
                     reyUU = data([S_UU, S_VV, S_WW, S_UV, S_UW, S_VW],i,j,k)
                     errReyMean = max(errReyMean, abs(reyU - (reyMean + x)))
                     do c = 1, 6
                        errReyFluc(c) = max(errReyFluc(c), abs(reyUU(c) - reyU(PAIRS(1,c))*reyU(PAIRS(2,c)) - reyFluc(c)))
                     end do
                  end if
!
!                 Favre: u_i~ = <rho u_i>/<rho>; the Favre block stores <rho u_i u_j>
                  favU  = data(nR+IRHOU:nR+IRHOW,i,j,k) / rhoBar
                  favUU = data(nR+NCONS+1:nR+NCONS+NO_OF_FAVRE_VARS,i,j,k)
                  errFavMean = max(errFavMean, abs(favU - (favMean + x)))
                  do c = 1, 6
!                    u_i''u_j''~ = <rho u_i u_j>/<rho> - u_i~ u_j~
                     errFavFluc(c)   = max(errFavFluc(c), &
                                           abs(favUU(c)/rhoBar - favU(PAIRS(1,c))*favU(PAIRS(2,c)) - favFluc(c)))
!                    Favre stress (as exported by horses2plt): <rho u_i''u_j''> = <rho u_i u_j> - <rho u_i><rho u_j>/<rho>
                     errFavStress(c) = max(errFavStress(c), &
                                           abs(favUU(c) - rhoBar*favU(PAIRS(1,c))*favU(PAIRS(2,c)) - rhoMean*favFluc(c)))
                  end do
               end do                ; end do                ; end do
               end associate
            end do
!
!           Only the maximum error over all nodes is asserted for each quantity
!           -------------------------------------------------------------------
            CALL FTAssertEqual(0.0_RP, errRho, TOL, tag // " <rho>")
            do c = 1, NDIM
               CALL FTAssertEqual(0.0_RP, errFavMean(c), TOL, tag // " Favre mean " // COMP(c)(1:1) // "~")
            end do
            do c = 1, 6
               CALL FTAssertEqual(0.0_RP, errFavFluc(c),   TOL, tag // " Favre " // COMP(c) // "''~")
               CALL FTAssertEqual(0.0_RP, errFavStress(c), TOL, tag // " Favre stress <rho " // COMP(c) // "''>")
            end do
            if (nR .gt. 0) then
               do c = 1, NDIM
                  CALL FTAssertEqual(0.0_RP, errReyMean(c), TOL, tag // " Reynolds mean " // COMP(c)(1:1))
               end do
               do c = 1, 6
                  CALL FTAssertEqual(0.0_RP, errReyFluc(c), TOL, tag // " Reynolds <" // COMP(c) // "'>")
               end do
            end if

            CALL sharedManager % summarizeAssertions(title = testName,iUnit = 6)
   
            IF ( sharedManager % numberOfAssertionFailures() == 0 )     THEN
               WRITE(6,*) testName, " ... Passed"
            ELSE
               WRITE(6,*) testName, " ... Failed"
               error stop 99
            END IF 
            WRITE(6,*)
            
            CALL finalizeSharedAssertionsManager
            CALL detachSharedAssertionsManager
#endif

         END SUBROUTINE UserDefinedFinalize
!
!//////////////////////////////////////////////////////////////////////// 
! 
      SUBROUTINE UserDefinedTermination
!
!        -----------------------------------------------
!        Called at the end of the main driver after 
!        everything else is done.
!        -----------------------------------------------
!
         IMPLICIT NONE  
      END SUBROUTINE UserDefinedTermination
      