!
!////////////////////////////////////////////////////////////////////////
!
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
module MMS
   use SMConstants
   use SolutionFile
   implicit none
   private
   public evalManufacturedSolution, compareLambFile

   contains
      pure subroutine evalManufacturedSolution(x, y, z, t, res)
         real(rp), intent(in)  :: x, y, z, t
         real(rp), intent(out) :: res(5)

         real(rp) :: rho_ms, u_ms, v_ms, w_ms, e_ms

         rho_ms = 3.0_rp + sin(x+y+z)*exp(t)
         u_ms = cos(x+y+z)*exp(-t)
         v_ms = sin(x+y+z)*exp(-t)
         w_ms = cos(x+y+z)*exp(-t)
         e_ms = rho_ms * (u_ms**2+v_ms**2+w_ms**2) / 2.0_rp + 1.0_rp

         res(1) = rho_ms
         res(2) = rho_ms * u_ms
         res(3) = rho_ms * v_ms
         res(4) = rho_ms * w_ms
         res(5) = rho_ms * e_ms
      end subroutine

      pure subroutine evalLambVector(x, y, z, t, res)
         real(rp), intent(in)  :: x, y, z, t
         real(rp), intent(out) :: res(NDIM)

         res(1) = (exp(-t)*sin(x+y+z)+exp(-t)*cos(x+y+z))*exp(-t)*sin(x+y+z)
         res(2) = (-2.0_rp*exp(-t)*sin(x+y+z)-2.0_rp*exp(-t)*cos(x+y+z))*exp(-t)*cos(x+y+z)
         res(3) = (exp(-t)*sin(x+y+z)+exp(-t)*cos(x+y+z))*exp(-t)*sin(x+y+z)
      end subroutine

      pure function evalLambVectorDivergence(x, y, z, t) result(res)
         real(rp), intent(in)  :: x, y, z, t
         real(rp)              :: res
!
!        The Lamb vector depends only on x+y+z, and Lx+Ly+Lz = -2 exp(-2t) cos(2(x+y+z))
!
         res = 4.0_rp*exp(-2.0_rp*t)*sin(2.0_rp*(x+y+z))
      end function
!
!     Reads a Lamb vector (nVars = NDIM) or divergence (nVars = 1) file record by record.
!     Its records must be the elements ids(:), in order, with the analytical values at time t.
!
      subroutine compareLambFile(fileName, nVars, elements, ids, t, tol, success)
         use ElementClass
         character(len=*), intent(in)    :: fileName
         integer,          intent(in)    :: nVars
         type(Element),    intent(in)    :: elements(:)
         integer,          intent(in)    :: ids(:)
         real(rp),         intent(in)    :: t, tol
         logical,          intent(inout) :: success

         integer               :: fid, n, i, j, k, rank, dims(4)
         real(rp)              :: x(NDIM), expected(NDIM), maxError
         real(rp), allocatable :: u_h(:,:,:,:)

         if ( getSolutionFileNoOfElements(fileName) /= size(ids) ) then
            print *, trim(fileName), ": number of elements ", getSolutionFileNoOfElements(fileName), &
                     ", expected ", size(ids)
            success = .false.
            return
         end if

         maxError = 0.0_rp
         open(newunit=fid, file=fileName, status="old", action="read", form="unformatted", access="stream")
         read(fid, pos=POS_INIT_DATA)
         do n = 1, size(ids)
            associate( e => elements(ids(n)) )
            read(fid) rank, dims
            if ( rank /= 4 .or. any(dims /= [nVars, e % Nxyz + 1]) ) then
               print *, trim(fileName), ": record ", n, " has dimensions ", dims
               success = .false.
               close(fid)
               return
            end if
            allocate(u_h(nVars, 0:e % Nxyz(1), 0:e % Nxyz(2), 0:e % Nxyz(3)))
            read(fid) u_h
            do k = 0, e % Nxyz(3) ; do j = 0, e % Nxyz(2) ; do i = 0, e % Nxyz(1)
               x = e % geom % x(:,i,j,k)
               if ( nVars == NDIM ) then
                  call evalLambVector(x(1), x(2), x(3), t, expected)
               else
                  expected(1) = evalLambVectorDivergence(x(1), x(2), x(3), t)
               end if
               maxError = max(maxError, maxval(abs(u_h(:,i,j,k) - expected(1:nVars))))
            end do                ; end do                ; end do
            deallocate(u_h)
            end associate
         end do
         close(fid)

         print *, trim(fileName), ": max. pointwise error ", maxError
         if ( maxError > tol ) success = .false.

      end subroutine compareLambFile
end module MMS

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
         use ElementConnectivityDefinitions
         implicit none
         real(kind=RP), intent(in)          :: x(NDIM)
         real(kind=RP), intent(in)          :: t
         real(kind=RP), intent(in)          :: nHat(NDIM)
         real(kind=RP), intent(inout)       :: Q(NCONS)
         type(Thermodynamics_t), intent(in) :: thermodynamics_
         type(Dimensionless_t),  intent(in) :: dimensionless_
         type(RefValues_t),      intent(in) :: refValues_
      end subroutine UserDefinedState_f

      subroutine UserDefinedGradVars_f(x, t, nHat, Q, U, GetGradients, thermodynamics_, dimensionless_, refValues_)
         use SMConstants
         use PhysicsStorage
         use FluidData
         use VariableConversion, only: GetGradientValues_f
         implicit none
         real(kind=RP), intent(in)          :: x(NDIM)
         real(kind=RP), intent(in)          :: t
         real(kind=RP), intent(in)          :: nHat(NDIM)
         real(kind=RP), intent(in)          :: Q(NCONS)
         real(kind=RP), intent(inout)       :: U(NGRAD)
         procedure(GetGradientValues_f)     :: GetGradients
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
      SUBROUTINE UserDefinedPeriodicOperation_f(mesh, time, dt, Monitors & 
#ifdef FLOW
         , thermodynamics_ &
         , dimensionless_  &
         , refValues_ & 
#endif   
#ifdef CAHNHILLIARD
         , multiphase_ &
#endif
      )
         use SMConstants
         USE HexMeshClass
         use MonitorsClass
         use fluiddata
         use physicsstorage
         IMPLICIT NONE
         CLASS(HexMesh)               :: mesh
         REAL(KIND=RP)                :: time
         REAL(KIND=RP)                :: dt
         type(Monitor_t), intent(in) :: monitors
#ifdef FLOW
         type(Thermodynamics_t), intent(in)    :: thermodynamics_
         type(Dimensionless_t),  intent(in)    :: dimensionless_
         type(RefValues_t),      intent(in)    :: refValues_
#endif
#ifdef CAHNHILLIARD
         type(Multiphase_t),     intent(in)    :: multiphase_
#endif
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
            use MMS
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
            do eID = 1, mesh % no_of_elements
               associate( e => mesh % elements(eID) )
               do k = 0, e % Nxyz(3);  do j = 0, e % Nxyz(2);  do i = 0, e % Nxyz(1)
                  call evalManufacturedSolution(e % geom % x(1,i,j,k), e % geom % x(2,i,j,k), e % geom % x(3,i,j,k), &
                                                0.0_RP, e % storage % Q(:,i,j,k))
               end do;        end do;        end do
               end associate
            end do
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
#ifdef MULTIPHASE
               do k = 0, Nz;  do j = 0, Ny;  do i = 0, Nx 
                  mesh % elements(eID) % storage % q(:,i,j,k) = [1.0_RP, 1.0_RP,0.0_RP,0.0_RP,0.0_RP] 
               end do;        end do;        end do
#endif
               end associate
               end associate
            end do
#endif
!
!           -----------------------------------
!           Acoustic default initial condition
!           -----------------------------------
!
#if defined(ACOUSTIC)
            do eID = 1, mesh % no_of_elements
               associate( Nx => mesh % elements(eID) % Nxyz(1), &
                          ny => mesh % elemeNts(eID) % nxyz(2), &
                          Nz => mesh % elements(eID) % Nxyz(3) )
               do k = 0, Nz;  do j = 0, Ny;  do i = 0, Nx 
                  ! mesh % elements(eID) % storage % q(:,i,j,k) = [0.0_RP, 0.0_RP,0.0_RP,0.0_RP,0.0_RP] 
                  mesh % elements(eID) % storage % q(:,i,j,k) = 0.0_RP 
               end do;        end do;        end do
               end associate
            end do
#endif
         end subroutine UserDefinedInitialCondition

         subroutine UserDefinedState1(x, t, nHat, Q & 
#ifdef FLOW
            , thermodynamics_ &
            , dimensionless_  &
            , refValues_ & 
#endif
#ifdef CAHNHILLIARD
            , multiphase_ &
#endif
            )
            use SMConstants
            use PhysicsStorage
            use FluidData
            implicit none
            real(kind=RP), intent(in)     :: x(NDIM)
            real(kind=RP), intent(in)     :: t
            real(kind=RP), intent(in)     :: nHat(NDIM)
            real(kind=RP), intent(inout)  :: Q(NCONS)
#ifdef FLOW
            type(Thermodynamics_t),    intent(in)  :: thermodynamics_
            type(Dimensionless_t),     intent(in)  :: dimensionless_
            type(RefValues_t),         intent(in)  :: refValues_
#endif
#ifdef CAHNHILLIARD
            type(Multiphase_t),     intent(in)  :: multiphase_
#endif
         end subroutine UserDefinedState1

#ifdef FLOW

         subroutine UserDefinedGradVars1(x, t, nHat, Q, U, GetGradients, thermodynamics_, dimensionless_, refValues_)
            use SMConstants
            use PhysicsStorage
            use FluidData
            use VariableConversion, only: GetGradientValues_f
            implicit none
            real(kind=RP), intent(in)          :: x(NDIM)
            real(kind=RP), intent(in)          :: t
            real(kind=RP), intent(in)          :: nHat(NDIM)
            real(kind=RP), intent(in)          :: Q(NCONS)
            real(kind=RP), intent(inout)       :: U(NGRAD)
            procedure(GetGradientValues_f)     :: GetGradients
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
         SUBROUTINE UserDefinedPeriodicOperation(mesh, time, dt, Monitors &
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
!           ----------------------------------------------------------
!           Called before every time-step to allow periodic operations
!           to be performed
!           ----------------------------------------------------------
!
            use SMConstants
            USE HexMeshClass
            use MonitorsClass
            use fluiddata
            use physicsstorage
            IMPLICIT NONE
            CLASS(HexMesh)               :: mesh
            REAL(KIND=RP)                :: time
            REAL(KIND=RP)                :: dt
            type(Monitor_t), intent(in) :: monitors
#ifdef FLOW
            type(Thermodynamics_t), intent(in)    :: thermodynamics_
            type(Dimensionless_t),  intent(in)    :: dimensionless_
            type(RefValues_t),      intent(in)    :: refValues_
#endif
#ifdef CAHNHILLIARD
            type(Multiphase_t),     intent(in)    :: multiphase_
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
!           --------------------------------------------------------
!           Called after the solution computed to allow, for example
!           error tests to be performed
!           --------------------------------------------------------
!
            use SMConstants
            USE HexMeshClass
            use PhysicsStorage
            use FluidData
            use MonitorsClass
            use FTAssertions
            use MMS
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
#if defined(NAVIERSTOKES)
            CHARACTER(LEN=29)                  :: testName = "Lamb vector cube"
            TYPE(FTAssertionsManager), POINTER :: sharedManager
            ! First time-based save, "save lamb vector interval" in LambVectorCube.control
            real(kind=RP), parameter           :: tSave = 1.5e-6_RP
            ! Same cube as "lamb vector cube" in LambVectorCube.control
            real(kind=RP), parameter           :: cube(6) = [0.0_RP, 2.0_RP, -1.0_RP, 1.5_RP, 0.5_RP, 3.0_RP]
            integer, allocatable               :: ids(:)
            integer                            :: eID
            logical                            :: success

            CALL initializeSharedAssertionsManager
            sharedManager => sharedAssertionsManager()

            call mesh % findElementsInsideCube(cube, ids)
            CALL FTAssertEqual(expectedValue = .true., &
                               actualValue   = size(ids) > 0 .and. size(ids) < size(mesh % elements), &
                               msg           = "The cube should contain some, but not all, elements.")

            success = .true.
            call compareLambFile("RESULTS/LambVectorCube_0000000002.Lamb.hsol", NDIM, mesh % elements, ids, &
                                 tSave, 0.01_RP, success)
            CALL FTAssertEqual(expectedValue = .true., &
                               actualValue   = success, &
                               msg           = "Error when comparing the Lamb vector saved in the cube.")

            success = .true.
            call compareLambFile("RESULTS/LambVectorCube_0000000002.divLamb.hsol", 1, mesh % elements, ids, &
                                 tSave, 0.05_RP, success)
            CALL FTAssertEqual(expectedValue = .true., &
                               actualValue   = success, &
                               msg           = "Error when comparing the Lamb vector divergence saved in the cube.")

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
