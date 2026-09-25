!
!//////////////////////////////////////////////////////////////////////////////////////
!
!     Unit test: consistency between the scalar wall shear stress / friction velocity
!     (u_tau, tau) and their vector counterparts (u_tau_vector, tau_x/y/z).
!
!     For a set of wall normals and velocity-gradient tensors, and for every gradient
!     variable set (STATE, ENERGY, ENTROPY), it checks that:
!
!        1) |u_tau_vector| = u_tau                (solver: getFrictionVelocityVector vs.
!                                                   getFrictionVelocityMagnitude)
!        2) tau_vector . n = 0                    (the vector is tangent to the wall)
!        3) rho * |u_tau_vector| * u_tau_vector = tangential traction computed here
!           independently from the velocity gradient
!        4) sqrt(tau_x^2+tau_y^2+tau_z^2) = tau   (horses2plt: TAUX_V/TAUY_V/TAUZ_V vs. Tauw_V)
!        5) sqrt(sqrt(u_t1^4 + u_t2^4)) = u_tau   (per-tangent getFrictionVelocity components,
!                                                   for an orthonormal tangent basis t1, t2)
!
!//////////////////////////////////////////////////////////////////////////////////////
!
program TestWallShearStress
   use SMConstants
   use FTValueDictionaryClass
   use PhysicsStorage_NS
   use FluidData_NS
   use VariableConversion_NS, only: Temperature, SutherlandsLaw
   use Physics_NS,            only: getFrictionVelocity, getFrictionVelocityVector, getFrictionVelocityMagnitude
   implicit none

   real(kind=RP), parameter :: TOL = 1.0e-10_RP
   integer,       parameter :: NO_OF_NORMALS   = 6
   integer,       parameter :: NO_OF_GRADIENTS = 5

   type(FTValueDictionary)  :: controlVariables
   real(kind=RP)            :: timeRef
   logical                  :: success
   integer                  :: iGradVars, iNormal, iGrad, numberOfFailures, numberOfChecks
   integer                  :: gradVarsList(3)
   character(len=8)         :: gradVarsNames(3)
   real(kind=RP)            :: normals(NDIM, NO_OF_NORMALS)
   real(kind=RP)            :: gradients(NDIM, NDIM, NO_OF_GRADIENTS)

   numberOfFailures = 0
   numberOfChecks   = 0
!
!  Set up the Navier-Stokes physics (air, Mach 0.3, Re 1000)
!  ---------------------------------------------------------
   call controlVariables % initWithSize(8)
   call controlVariables % addValueForKey("ns",     "flow equations")
   call controlVariables % addValueForKey("0.3",    "mach number")
   call controlVariables % addValueForKey("1000.0", "reynolds number")
   call ConstructPhysicsStorage_NS(controlVariables, 1.0_RP, timeRef, success)
   if ( .not. success ) then
      print*, "TestWallShearStress: could not construct the NS physics storage"
      error stop 1
   end if

   gradVarsList  = [GRADVARS_STATE, GRADVARS_ENERGY, GRADVARS_ENTROPY]
   gradVarsNames = ["STATE   ", "ENERGY  ", "ENTROPY "]
!
!  Wall normals: aligned with the axes and oblique
!  -----------------------------------------------
   normals(:,1) = [ 0.0_RP,  1.0_RP,  0.0_RP]
   normals(:,2) = [ 0.0_RP,  0.0_RP, -1.0_RP]
   normals(:,3) = [ 1.0_RP,  1.0_RP,  0.0_RP]
   normals(:,4) = [ 1.0_RP, -2.0_RP,  3.0_RP]
   normals(:,5) = [-0.3_RP,  0.7_RP,  0.2_RP]
   normals(:,6) = [ 0.1_RP,  0.2_RP, -0.9_RP]
   do iNormal = 1, NO_OF_NORMALS
      normals(:,iNormal) = normals(:,iNormal) / norm2(normals(:,iNormal))
   end do
!
!  Velocity gradients, grad(i,j) = du_i/dx_j
!  -----------------------------------------
!  1) Pure shear du/dy (the shear is aligned with one tangent for normal 1)
   gradients(:,:,1) = 0.0_RP
   gradients(1,2,1) = 1.0_RP
!  2) Shear in two directions: du/dy and dw/dy (crossflow, NOT aligned with any tangent)
   gradients(:,:,2) = 0.0_RP
   gradients(1,2,2) = 2.0_RP
   gradients(3,2,2) = 1.5_RP
!  3) Full, non-symmetric tensor with dilatation
   gradients(:,:,3) = reshape([ 0.3_RP, -1.2_RP,  0.8_RP, &
                                2.1_RP,  0.5_RP, -0.4_RP, &
                               -0.7_RP,  1.9_RP, -0.2_RP], [NDIM, NDIM])
!  4) Another full tensor
   gradients(:,:,4) = reshape([-1.0_RP,  0.4_RP,  2.2_RP, &
                                0.9_RP,  0.1_RP, -3.0_RP, &
                                1.3_RP, -0.6_RP,  0.7_RP], [NDIM, NDIM])
!  5) Rigid rotation: no viscous stress, u_tau must vanish
   gradients(:,:,5) = reshape([ 0.0_RP,  1.0_RP, -2.0_RP, &
                               -1.0_RP,  0.0_RP,  0.5_RP, &
                                2.0_RP, -0.5_RP,  0.0_RP], [NDIM, NDIM])

   do iGradVars = 1, size(gradVarsList)
      call SetGradientVariables(gradVarsList(iGradVars))
      do iNormal = 1, NO_OF_NORMALS ; do iGrad = 1, NO_OF_GRADIENTS
         call CheckWallShearStress(normals(:,iNormal), gradients(:,:,iGrad), &
                                   trim(gradVarsNames(iGradVars)), iNormal, iGrad)
      end do                        ; end do
   end do

   print "(A,I0,A,I0,A)", "TestWallShearStress: ", numberOfChecks - numberOfFailures, " of ", numberOfChecks, " checks passed"
   if ( numberOfFailures .gt. 0 ) then
      print*, "TestWallShearStress: FAILED"
      error stop 1
   end if
   print*, "TestWallShearStress: PASSED"

contains

   subroutine CheckWallShearStress(normal, grad, gradVarsName, iNormal, iGrad)
      implicit none
      real(kind=RP),    intent(in) :: normal(NDIM)
      real(kind=RP),    intent(in) :: grad(NDIM, NDIM)
      character(len=*), intent(in) :: gradVarsName
      integer,          intent(in) :: iNormal, iGrad
!
!     ---------------
!     Local variables
!     ---------------
!
      real(kind=RP) :: rho, p, vel(NDIM)
      real(kind=RP) :: Q(NCONS), Q_x(NGRAD), Q_y(NGRAD), Q_z(NGRAD)
      real(kind=RP) :: t1(NDIM), t2(NDIM)
      real(kind=RP) :: u_tau, u_tau_vec(NDIM), u_tau_t1, u_tau_t2
      real(kind=RP) :: mu, divV, S(NDIM, NDIM), traction(NDIM), tangentialTraction(NDIM)
      real(kind=RP) :: tau_scalar, tau_vec(NDIM), ref
      integer       :: i
      character(len=64) :: caseName

      write(caseName,'(A,A,I0,A,I0)') gradVarsName, " normal ", iNormal, " gradient ", iGrad
!
!     Flow state: non-trivial density, velocity and pressure
!     ------------------------------------------------------
      rho = 1.3_RP
      vel = [0.4_RP, -0.2_RP, 0.1_RP]
      p   = 1.1_RP / dimensionless % gammaM2
      Q(IRHO)        = rho
      Q(IRHOU:IRHOW) = rho * vel
      Q(IRHOE)       = p * thermodynamics % InvGammaMinus1 + 0.5_RP * rho * sum(vel*vel)
!
!     Gradients of the selected gradient variables that reproduce grad(:,:), with uniform rho and p
!     ---------------------------------------------------------------------------------------------
      Q_x = 0.0_RP ; Q_y = 0.0_RP ; Q_z = 0.0_RP
      select case (grad_vars)
      case (GRADVARS_STATE)      ! d(rho u)/dx = rho du/dx
         Q_x(IRHOU:IRHOW) = rho * grad(:,1)
         Q_y(IRHOU:IRHOW) = rho * grad(:,2)
         Q_z(IRHOU:IRHOW) = rho * grad(:,3)
      case (GRADVARS_ENERGY)     ! velocity gradient itself
         Q_x(IRHOU:IRHOW) = grad(:,1)
         Q_y(IRHOU:IRHOW) = grad(:,2)
         Q_z(IRHOU:IRHOW) = grad(:,3)
      case (GRADVARS_ENTROPY)    ! d(rho u/p)/dx = rho/p du/dx (with d(-rho/p)/dx = 0)
         Q_x(IRHOU:IRHOW) = rho / p * grad(:,1)
         Q_y(IRHOU:IRHOW) = rho / p * grad(:,2)
         Q_z(IRHOU:IRHOW) = rho / p * grad(:,3)
      end select
!
!     Orthonormal tangent basis
!     -------------------------
      if ( abs(normal(1)) .lt. 0.9_RP ) then
         t1 = [1.0_RP, 0.0_RP, 0.0_RP]
      else
         t1 = [0.0_RP, 1.0_RP, 0.0_RP]
      end if
      t1 = t1 - dot_product(t1, normal) * normal
      t1 = t1 / norm2(t1)
      t2 = [normal(2)*t1(3) - normal(3)*t1(2), normal(3)*t1(1) - normal(1)*t1(3), normal(1)*t1(2) - normal(2)*t1(1)]
!
!     Independent reference: tangential traction -(tau . n)_t from the Newtonian stress tensor
!     -----------------------------------------------------------------------------------------
      mu   = dimensionless % mu * SutherlandsLaw(Temperature(Q))
      divV = grad(1,1) + grad(2,2) + grad(3,3)
      S    = mu * (grad + transpose(grad))
      do i = 1, NDIM
         S(i,i) = S(i,i) - 2.0_RP / 3.0_RP * mu * divV
      end do
      traction           = -matmul(S, normal)
      tangentialTraction = traction - dot_product(traction, normal) * normal
      ref = max(norm2(tangentialTraction), 1.0_RP)
!
!     Solver quantities
!     -----------------
      call getFrictionVelocityMagnitude(Q, Q_x, Q_y, Q_z, normal, u_tau)
      call getFrictionVelocityVector   (Q, Q_x, Q_y, Q_z, normal, u_tau_vec)
      call getFrictionVelocity         (Q, Q_x, Q_y, Q_z, normal, t1, u_tau_t1)
      call getFrictionVelocity         (Q, Q_x, Q_y, Q_z, normal, t2, u_tau_t2)
!
!     horses2plt quantities (OutputVariables.f90: Tauw_V and TAUX_V/TAUY_V/TAUZ_V)
!     ---------------------------------------------------------------------------
      tau_scalar = rho * u_tau**2 * sign(1.0_RP, u_tau)
      tau_vec    = rho * norm2(u_tau_vec) * u_tau_vec
!
!     Checks
!     ------
!     All comparisons are made in stress units (rho * u_tau^2): the friction velocity scales as the
!     square root of the stress, so round-off in a vanishing wall shear stress (e.g. rigid rotation,
!     or a traction normal to the wall) would otherwise be amplified to O(sqrt(epsilon)).
      call Check(abs(rho * norm2(u_tau_vec)**2 - rho * u_tau**2) .le. TOL * ref, &
                 caseName, "|u_tau_vector| = u_tau", norm2(u_tau_vec), u_tau)
      call Check(abs(dot_product(tau_vec, normal)) .le. TOL * ref, &
                 caseName, "tau_vector . n = 0", dot_product(tau_vec, normal), 0.0_RP)
      call Check(norm2(tau_vec - tangentialTraction) .le. TOL * ref, &
                 caseName, "tau_vector = tangential traction", norm2(tau_vec), norm2(tangentialTraction))
      call Check(abs(norm2(tau_vec) - tau_scalar) .le. TOL * ref, &
                 caseName, "|tau_vector| = tau", norm2(tau_vec), tau_scalar)
      call Check(abs(rho * sqrt(u_tau_t1**4 + u_tau_t2**4) - rho * u_tau**2) .le. TOL * ref, &
                 caseName, "(u_t1^4 + u_t2^4)^(1/4) = u_tau", sqrt(sqrt(u_tau_t1**4 + u_tau_t2**4)), u_tau)

   end subroutine CheckWallShearStress

   subroutine Check(condition, caseName, checkName, obtained, expected)
      implicit none
      logical,          intent(in) :: condition
      character(len=*), intent(in) :: caseName, checkName
      real(kind=RP),    intent(in) :: obtained, expected

      numberOfChecks = numberOfChecks + 1
      if ( .not. condition ) then
         numberOfFailures = numberOfFailures + 1
         print "(A,A,A,A,A,ES24.16,A,ES24.16)", "   FAILED [", trim(caseName), "] ", checkName, &
                                               ": obtained ", obtained, ", expected ", expected
      end if

   end subroutine Check

end program TestWallShearStress
