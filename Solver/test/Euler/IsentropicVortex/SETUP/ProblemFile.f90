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
#if defined(NAVIERSTOKES)
                                        , thermodynamics_ &
                                        , dimensionless_  &
                                        , refValues_ & 
#endif
#if defined(CAHNHILLIARD)
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
            class(HexMesh)                      :: mesh
#if defined(NAVIERSTOKES)
            type(Thermodynamics_t), intent(in)  :: thermodynamics_
            type(Dimensionless_t),  intent(in)  :: dimensionless_
            type(RefValues_t),      intent(in)  :: refValues_
#endif
#if defined(CAHNHILLIARD)
            type(Multiphase_t),     intent(in)  :: multiphase_
#endif
         END SUBROUTINE UserDefinedFinalSetup
!
!//////////////////////////////////////////////////////////////////////// 
! 
         subroutine UserDefinedInitialCondition(mesh &
#if defined(NAVIERSTOKES)
                                        , thermodynamics_ &
                                        , dimensionless_  &
                                        , refValues_ & 
#endif
#if defined(CAHNHILLIARD)
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
            implicit none
            class(hexmesh)                      :: mesh
#if defined(NAVIERSTOKES)
            type(Thermodynamics_t), intent(in)  :: thermodynamics_
            type(Dimensionless_t),  intent(in)  :: dimensionless_
            type(RefValues_t),      intent(in)  :: refValues_
#endif
#if defined(CAHNHILLIARD)
            type(Multiphase_t),     intent(in)  :: multiphase_
#endif
!
!           ---------------
!           Local variables
!           ---------------
!
            REAL(KIND=RP) :: x(3)        
            INTEGER       :: i, j, k, eID
            REAL(KIND=RP) :: rho , u , v , w , p
            REAL(KIND=RP) :: r, beta, T, Mach
            REAL(KIND=RP) :: x0, y0, xc, yc, dT, uInf
            integer       :: Nx, Ny, Nz

#if defined(NAVIERSTOKES)
            

            beta = 5.0_RP
            x0 = 0.0_RP
            y0 = 0.0_RP
            uInf = 1.0_RP

            associate( gamma => thermodynamics_ % gamma ) 
            DO eID = 1, SIZE(mesh % elements)
               Nx = mesh % elements(eID) % Nxyz(1)
               Ny = mesh % elements(eID) % Nxyz(2)
               Nz = mesh % elements(eID) % Nxyz(3)

               DO k = 0, Nz
                  DO j = 0, Ny
                     DO i = 0, Nx 

                         x = mesh % elements(eID) % geom % x(:,i,j,k)
                         
                         xc = x(1) - x0
                         yc = x(2) - y0
                         r = sqrt(xc*xc + yc*yc)
                         
                         dT = -(gamma - 1.0_RP) * beta * beta / (8.0_RP * gamma * PI * PI) * exp(1.0_RP - r*r)
                         T = 1.0_RP + dT
                         rho = T**(1.0_RP / (gamma - 1.0_RP))
                         u = uInf - beta * yc / (2.0_RP * PI) * exp(0.5_RP * (1.0_RP - r*r))
                         v = beta * xc / (2.0_RP * PI) * exp(0.5_RP * (1.0_RP - r*r))
                         w = 0.0_RP
                         p = rho**gamma

                         mesh % elements(eID) % storage % Q(1,i,j,k) = rho
                         mesh % elements(eID) % storage % Q(2,i,j,k) = rho*u
                         mesh % elements(eID) % storage % Q(3,i,j,k) = rho*v
                         mesh % elements(eID) % storage % Q(4,i,j,k) = rho*w
                         mesh % elements(eID) % storage % Q(5,i,j,k) = p / (gamma - 1.0_RP) + 0.5_RP * rho * (u*u + v*v + w*w)

                     END DO
                  END DO
               END DO 
               
            END DO 
            end associate
#endif

!
!           ---------------------------------------
!           Cahn-Hilliard default initial condition
!           ---------------------------------------
!
#if defined(CAHNHILLIARD)
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
#if defined(NAVIERSTOKES)
         subroutine UserDefinedState1(x, t, nHat, Q, thermodynamics_, dimensionless_, refValues_)
!
!           -------------------------------------------------
!           Used to define an user defined boundary condition
!           -------------------------------------------------
!
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

         subroutine UserDefinedNeumann1(x, t, nHat, U_x, U_y, U_z)
!
!           --------------------------------------------------------
!           Used to define a Neumann user defined boundary condition
!           --------------------------------------------------------
!
            use SMConstants
            use PhysicsStorage
            use FluidData
            implicit none
            real(kind=RP), intent(in)     :: x(NDIM)
            real(kind=RP), intent(in)     :: t
            real(kind=RP), intent(in)     :: nHat(NDIM)
            real(kind=RP), intent(inout)  :: U_x(NGRAD)
            real(kind=RP), intent(inout)  :: U_y(NGRAD)
            real(kind=RP), intent(inout)  :: U_z(NGRAD)
         end subroutine UserDefinedNeumann1
#endif
!
!//////////////////////////////////////////////////////////////////////// 
! 
         SUBROUTINE UserDefinedPeriodicOperation(mesh, time, dt, Monitors)
!
!           ----------------------------------------------------------
!           Called at the output interval to allow periodic operations
!           to be performed
!           ----------------------------------------------------------
!
            use SMConstants
            USE HexMeshClass
#if defined(NAVIERSTOKES)
            use MonitorsClass
#endif
            IMPLICIT NONE
            class(HexMesh)               :: mesh
            real(kind=RP)                :: time
            real(kind=RP)                :: dt
#if defined(NAVIERSTOKES)
            type(Monitor_t), intent(in) :: monitors
#else
            logical, intent(in) :: monitors
#endif
         END SUBROUTINE UserDefinedPeriodicOperation
!
!//////////////////////////////////////////////////////////////////////// 
! 
#if defined(NAVIERSTOKES)
         subroutine UserDefinedSourceTermNS(x, Q, time, S, thermodynamics_, dimensionless_, refValues_)
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
            type(Thermodynamics_t),    intent(in)  :: thermodynamics_
            type(Dimensionless_t),     intent(in)  :: dimensionless_
            type(RefValues_t),         intent(in)  :: refValues_
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
#if defined(NAVIERSTOKES)
                                                    , thermodynamics_ &
                                                    , dimensionless_  &
                                                    , refValues_ & 
#endif   
#if defined(CAHNHILLIARD)
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
            use FTAssertions
            USE HexMeshClass
            use PhysicsStorage
            use FluidData
            use MonitorsClass
#if defined(_HAS_MPI_)
            use mpi
#endif
            IMPLICIT NONE
            class(HexMesh)                        :: mesh
            REAL(KIND=RP)                         :: time
            integer                               :: iter
            real(kind=RP)                         :: maxResidual
#if defined(NAVIERSTOKES)
            type(Thermodynamics_t), intent(in)    :: thermodynamics_
            type(Dimensionless_t),  intent(in)    :: dimensionless_
            type(RefValues_t),      intent(in)    :: refValues_
#endif
#if defined(CAHNHILLIARD)
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
#if defined(NAVIERSTOKES)
            CHARACTER(LEN=47)                  :: testName           = "Translating isentropic vortex and entropy cons."
            REAL(KIND=RP)                      :: maxError
            REAL(KIND=RP), ALLOCATABLE         :: QExpected(:,:,:,:)
            INTEGER                            :: eID
            INTEGER                            :: i, j, k, N
            TYPE(FTAssertionsManager), POINTER :: sharedManager
            LOGICAL                            :: success
            integer                            :: rank
            integer                            :: nProcs
            real(kind=RP)                      :: errorLocal(5)
#if defined(_HAS_MPI_)
            integer                            :: ierr
#endif
            real(kind=RP), parameter           :: entropyRate = -2.18259E-12_RP
            real(kind=RP), parameter           :: EntropBalance = -2.18259E-12_RP            
!
!           Reference residuals after 100 iterations of the translating vortex
!           (res(5), the energy residual, is kept for reference but not asserted)
!
            real(kind=RP), parameter           :: res(5) = [0.4135553213840333_RP, &
                                                            0.5247316995389606_RP, &
                                                            0.7298627036149967_RP, &
                                                            1.5117E-12_RP, &
                                                            0.0_RP]

            real(KIND=RP), parameter :: wLGL4(0:4) = [0.100000000000000_RP, &
                                               0.544444444444444_RP, &
                                               0.711111111111111_RP, &
                                               0.544444444444444_RP, &
                                               0.100000000000000_RP]
            real(kind=RP) :: error(5), locErr(5), Q(5)
            real(kind=RP) :: x, y, z
            real(kind=RP) :: r, beta, T, rho, u, v, w, p, xc, yc, dT, x0, y0
            real(kind=RP) :: uInf
!
!           Periodic domain extents (see MESH/isentropicvortex2D.geo)
!
            real(kind=RP), parameter :: Lx = 30.0_RP, Ly = 30.0_RP

            beta = 5.0_RP
            uInf = 1.0_RP
!
!           The vortex is advected with the free stream along the x-axis
!
            x0 = 0.0_RP + uInf * time
            y0 = 0.0_RP
            error = 0.0_RP

            DO eID = 1, mesh % no_of_elements
            ASSOCIATE(e => mesh % elements(eID), gamma => thermodynamics_ % gamma)
            DO k = 0, e % Nxyz(3); DO j = 0, e % Nxyz(2); DO i = 0, e % Nxyz(1)
               x = e % geom % x(IX,i,j,k)
               y = e % geom % x(IY,i,j,k)
               z = e % geom % x(IZ,i,j,k)
               
!
!              Distance to the nearest periodic image of the vortex core
!
               xc = MODULO(x - x0 + 0.5_RP*Lx, Lx) - 0.5_RP*Lx
               yc = MODULO(y - y0 + 0.5_RP*Ly, Ly) - 0.5_RP*Ly
               r = SQRT(xc*xc + yc*yc)
               
               dT = -(gamma - 1.0_RP) * beta * beta / (8.0_RP * gamma * PI * PI) * EXP(1.0_RP - r*r)
               T = 1.0_RP + dT
               rho = T**(1.0_RP / (gamma - 1.0_RP))
               u = uInf - beta * yc / (2.0_RP * PI) * EXP(0.5_RP * (1.0_RP - r*r))
               v = beta * xc / (2.0_RP * PI) * EXP(0.5_RP * (1.0_RP - r*r))
               w = 0.0_RP
               p = rho**gamma
               
               Q(1) = rho
               Q(2) = rho*u
               Q(3) = rho*v
               Q(4) = rho*w
               Q(5) = p / (gamma - 1.0_RP) + 0.5_RP * rho * (u*u + v*v + w*w)
               
               locErr(1) = e % storage % Q(1,i,j,k) - Q(1)
               locErr(2) = e % storage % Q(2,i,j,k) - Q(2)
               locErr(3) = e % storage % Q(3,i,j,k) - Q(3)
               locErr(4) = e % storage % Q(4,i,j,k) - Q(4)
               locErr(5) = e % storage % Q(5,i,j,k) - Q(5)
               error = error + e % geom % jacobian(i,j,k) * wLGL4(i) * wLGL4(j) * wLGL4(k) * locErr**2
            END DO; END DO; END DO
            END ASSOCIATE
         END DO
         
         errorLocal = error
         rank    = 0
         nProcs  = 1
#if defined(_HAS_MPI_)
         call mpi_comm_rank(MPI_COMM_WORLD, rank, ierr)
         call mpi_comm_size(MPI_COMM_WORLD, nProcs, ierr)
         if ( nProcs .gt. 1 ) then
            call MPI_Allreduce(errorLocal, error, 5, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD, ierr)
         else
            error = errorLocal
         end if
#else
         error = errorLocal
#endif
         error = SQRT(error)
         if ( rank .eq. 0 ) then
            PRINT*, 'L2 Error (rho):   ', error(1)
            PRINT*, 'L2 Error (rho*u): ', error(2)
            PRINT*, 'L2 Error (rho*v): ', error(3)
            PRINT*, 'L2 Error (rho*w): ', error(4)
            PRINT*, 'L2 Error (rho*E): ', error(5)
         end if
!
!        Write exact solution to file
!
         DO eID = 1, mesh % no_of_elements
            ASSOCIATE(e => mesh % elements(eID), gamma => thermodynamics_ % gamma)
            DO k = 0, e % Nxyz(3); DO j = 0, e % Nxyz(2); DO i = 0, e % Nxyz(1)
               x = e % geom % x(IX,i,j,k)
               y = e % geom % x(IY,i,j,k)
               z = e % geom % x(IZ,i,j,k)
               
!
!              Distance to the nearest periodic image of the vortex core
!
               xc = MODULO(x - x0 + 0.5_RP*Lx, Lx) - 0.5_RP*Lx
               yc = MODULO(y - y0 + 0.5_RP*Ly, Ly) - 0.5_RP*Ly
               r = SQRT(xc*xc + yc*yc)
               
               dT = -(gamma - 1.0_RP) * beta * beta / (8.0_RP * gamma * PI * PI) * EXP(1.0_RP - r*r)
               T = 1.0_RP + dT
               rho = T**(1.0_RP / (gamma - 1.0_RP))
               u = uInf - beta * yc / (2.0_RP * PI) * EXP(0.5_RP * (1.0_RP - r*r))
               v = beta * xc / (2.0_RP * PI) * EXP(0.5_RP * (1.0_RP - r*r))
               w = 0.0_RP
               p = rho**gamma
               
               e % storage % Q(1,i,j,k) = rho
               e % storage % Q(2,i,j,k) = rho*u
               e % storage % Q(3,i,j,k) = rho*v
               e % storage % Q(4,i,j,k) = rho*w
               e % storage % Q(5,i,j,k) = p / (gamma - 1.0_RP) + 0.5_RP * rho * (u*u + v*v + w*w)
            END DO; END DO; END DO
            END ASSOCIATE
         END DO
         
         CALL mesh % SaveSolution(iter, time, "ExactSolution.hsol",.false.)

            CALL initializeSharedAssertionsManager
            sharedManager => sharedAssertionsManager()
            
            CALL FTAssertEqual(expectedValue = res(1) + 1.0_RP, &
                               actualValue   = monitors % residuals % values(1,1) + 1.0_RP, &
                               tol           = 1.0e-11_RP, &
                               msg           = "continuity residual + 1._RP")

            CALL FTAssertEqual(expectedValue = res(2) + 1.0_RP, &
                               actualValue   = monitors % residuals % values(2,1) + 1.0_RP, &
                               tol           = 1.0e-11_RP, &
                               msg           = "x-momentum residual + 1._RP")

            CALL FTAssertEqual(expectedValue = res(3) + 1.0_RP, &
                               actualValue   = monitors % residuals % values(3,1) + 1.0_RP, &
                               tol           = 1.0e-11_RP, &
                               msg           = "y-momentum residual + 1._RP")

            CALL FTAssertEqual(expectedValue = res(4) + 1.0_RP, &
                               actualValue   = monitors % residuals % values(4,1) + 1.0_RP, &
                               tol           = 1.0e-11_RP, &
                               msg           = "z-momentum residual + 1._RP")

            CALL FTAssertEqual(expectedValue = entropyRate + 1.0_RP, &
                               actualValue   = monitors % volumeMonitors(1) % values(1,1)+1.0_RP, &
                               tol           = 1.0e-11_RP, &
                               msg           = "Entropy Rate")

            CALL FTAssertEqual(expectedValue = EntropBalance + 1.0_RP, &
                   actualValue   = monitors % volumeMonitors(2) % values(1,1)+1.0_RP, &
                   tol           = 1.0e-11_RP, &
                   msg           = "Entropy balance")

            CALL sharedManager % summarizeAssertions(title = testName,iUnit = 6)
   
            IF ( sharedManager % numberOfAssertionFailures() == 0 )     THEN
               WRITE(6,*) testName, " ... Passed"
               WRITE(6,*) "This test case has no expected solution yet, only checks the residual after 100 iterations."
            ELSE
               WRITE(6,*) testName, " ... Failed"
               WRITE(6,*) "NOTE: Failure is expected when the max eigenvalue procedure is changed."
               WRITE(6,*) "      If that is done, re-compute the expected values and modify this procedure"
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
!        Called at the the end of the main driver after 
!        everything else is done.
!        -----------------------------------------------
!
         IMPLICIT NONE  
      END SUBROUTINE UserDefinedTermination
      