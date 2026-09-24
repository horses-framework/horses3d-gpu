#include "Includes.h"
module VariableConversion
#if defined(NAVIERSTOKES) && (!(SPALARTALMARAS))
   use VariableConversion_NS
#elif defined(SPALARTALMARAS)
   USE VariableConversion_NSSA
#elif defined(INCNS)
   use VariableConversion_iNS
#elif defined(MULTIPHASE)
   use VariableConversion_MU
#endif
#if defined(CAHNHILLIARD)
   use VariableConversion_CH
#endif
#if defined(ACOUSTIC)
   use VariableConversion_CAA
#endif
   implicit none

   abstract interface
      subroutine GetGradientValues_f(nEqn, nGradEqn, Q, U, rho_)
         use SMConstants, only: RP
         implicit none
         integer, intent(in)                 :: nEqn, nGradEqn
         real(kind=RP), intent(in)           :: Q(nEqn)
         real(kind=RP), intent(out)          :: U(nGradEqn)
         real(kind=RP), intent(in), optional :: rho_
      end subroutine GetGradientValues_f
   end interface

   contains
!
!/////////////////////////////////////////////////////////////////////////////////////////////
!
!     ---------------------------------------------------------------------------------------
!     Computes the gradient variables, U, from the state vector, Q, for the requested set of
!     gradient variables (gradVars = GRADVARS_STATE, GRADVARS_ENTROPY, GRADVARS_ENERGY).
!     Sets that are not implemented for the current physics fall back to the state variables.
!
!     rho and mu are only used by the multiphase solver (density and chemical potential).
!     ---------------------------------------------------------------------------------------
!
      pure subroutine GradientVariables_Selector(gradVars, nEqn, nGradEqn, Q, U, rho, mu)
         !$acc routine seq
         use SMConstants, only: RP
         use PhysicsStorage
         implicit none
         integer,       intent(in)  :: gradVars
         integer,       intent(in)  :: nEqn, nGradEqn
         real(kind=RP), intent(in)  :: Q(nEqn)
         real(kind=RP), intent(out) :: U(nGradEqn)
         real(kind=RP), intent(in)  :: rho
         real(kind=RP), intent(in)  :: mu

         select case (gradVars)
#if defined(NAVIERSTOKES)
         case (GRADVARS_ENTROPY)
            call NSGradientVariables_ENTROPY(nEqn, nGradEqn, Q, U)

         case (GRADVARS_ENERGY)
            call NSGradientVariables_ENERGY(nEqn, nGradEqn, Q, U)
#elif defined(INCNS)
         case (GRADVARS_ENTROPY)
            call iNSGradientVariables(nEqn, nGradEqn, Q, U)
#elif defined(MULTIPHASE)
         case (GRADVARS_ENTROPY)
!
!           The multiphase solver needs the chemical potential as first entropy variable
!           ----------------------------------------------------------------------------
            call mGradientVariables(nEqn, nGradEqn, Q, U, rho)
            U(IGMU) = mu
#endif
         case default
            U(1:nGradEqn) = Q(1:nGradEqn)

         end select

      end subroutine GradientVariables_Selector

#if (defined(CAHNHILLIARD) && defined(NAVIERSTOKES))
      pure subroutine GetNSCHViscosity(phi, mu)
         use SMConstants, only: RP
         use FluidData
         implicit none
         real(kind=RP), intent(in)     :: phi
         real(kind=RP), intent(out)    :: mu
!
!        ---------------
!        Local variables         
!        ---------------
!
         real(kind=RP)  :: cIn01, p

         cIn01 = 0.5_RP * (phi + 1.0_RP)
         p = POW3(cIn01) * (6.0_RP * POW2(cIn01) - 15.0_RP * cIn01 + 10.0_RP)

         mu = dimensionless % mu * ( (1.0_RP - p) + (p)*multiphase % viscRatio)

      end subroutine GetNSCHViscosity
#endif

#if (defined(INCNS) && defined(CAHNHILLIARD))
      pure subroutine GetiNSCHViscosity(phi, mu)
         use SMConstants, only: RP
         use FluidData
         implicit none
         real(kind=RP), intent(in)     :: phi
         real(kind=RP), intent(out)    :: mu
!
!        ---------------
!        Local variables         
!        ---------------
!
         real(kind=RP)  :: cIn01, p

         cIn01 = 0.5_RP * (phi + 1.0_RP)
         p = POW3(cIn01) * (6.0_RP * POW2(cIn01) - 15.0_RP * cIn01 + 10.0_RP)

         mu = dimensionless % mu(1) * (1.0_RP - p) + (p)*dimensionless % mu(2)
      
      end subroutine GetiNSCHViscosity
#endif

end module VariableConversion
