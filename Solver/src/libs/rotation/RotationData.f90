#include "Includes.h"
module RotationData
!
!  ***********************************************************************
!    This moule implements two features:
!    1) Single Rotating Reference Frame (SRF)
!    2) Periodic revolution BCs
!
!  The rotation geometry (axis, center) is shared between both.
!
!  Equation-set-agnostic: this module is built once (see libs/rotation),
!  not once per equation set. The "SRF only supported in Navier-Stokes"
!  check therefore lives in the equation-set-specific caller (DGSEMClass),
!  not here.
!  ***********************************************************************
!
   use SMConstants
   implicit none

   private
   public :: RotationParams_t, rotationParams
   public :: SetSRF, RotationAxisIndex
   public :: SetPeriodicAngle, RotateVectorAroundAxis
   public :: RotateStateMomentum, RotateMomentumGradients
   public :: InitializeRotationFromControl

   type RotationParams_t

      ! Shared rotation geometry (used by both SRF and periodic)
      real(kind=RP) :: axis(NDIM) = [0.0_RP, 0.0_RP, 1.0_RP]    ! Unit rotation axis
      real(kind=RP) :: center(NDIM) = [0.0_RP, 0.0_RP, 0.0_RP]  ! Center

      ! For SRF
      logical       :: srfEnabled = .false.
      real(kind=RP) :: omega(NDIM) = 0.0_RP  ! Non-dimensional rotation vector for SRF

      ! For periodic revolution BCs
      logical       :: periodicEnabled = .false.
      real(kind=RP) :: periodicAngle = 0.0_RP  ! Sector (rad)
   end type RotationParams_t

   type(RotationParams_t), protected :: rotationParams
   !$acc declare create(rotationParams)

   contains

   subroutine SetSRF(omega_)
!
!     Set the non-dimensional SRF omega vector.
!     Done in spatial discretization, not in this class because it requires FluidData ref values
!
      implicit none
      real(kind=RP), intent(in) :: omega_(NDIM)

      rotationParams % omega = omega_
      !$acc update device(rotationParams)
   end subroutine SetSRF


   pure integer function RotationAxisIndex()
!     Return the index of the rotation axis if it is aligned with x, y, or z.
      real(kind=RP) :: axisAbs(NDIM)
      integer       :: idx

      axisAbs = abs(rotationParams % axis)
      idx = maxloc(axisAbs, dim = 1)
      if (axisAbs(idx) < (1.0_RP - 1.0e-8_RP)) then
         RotationAxisIndex = 0
      else
         RotationAxisIndex = idx
      end if
   end function RotationAxisIndex

   subroutine SetPeriodicAngle(angle_)
!
!     Store the sector angle (rad) computed during periodic face matching
!
      implicit none
      real(kind=RP), intent(in) :: angle_
      rotationParams % periodicAngle = angle_
      !$acc update device(rotationParams)
   end subroutine SetPeriodicAngle

   subroutine InitializeRotationFromControl(controlVariables)
!
!     Read rotation geometry (axis, center) and enable/disable flags
!     from the control file.  Called early during mesh construction so that
!     periodic face matching can use the rotation information.
!
      use FTValueDictionaryClass
      use FileReadingUtilities, only: getRealArrayFromString
      implicit none
      class(FTValueDictionary), intent(in) :: controlVariables
      real(kind=RP) :: axis_(NDIM), center_(NDIM)

      axis_   = rotationParams % axis
      center_ = rotationParams % center

      ! ---- rotation axis ----
      if (controlVariables % containsKey("rotation axis")) then
         axis_ = getRealArrayFromString(controlVariables % stringValueForKey("rotation axis", LINE_LENGTH))
      end if

      ! ---- rotation center ----
      if (controlVariables % containsKey("rotation center")) then
         center_ = getRealArrayFromString(controlVariables % stringValueForKey("rotation center", LINE_LENGTH))
      end if

      rotationParams % axis   = axis_ / norm2(axis_)
      rotationParams % center = center_
      !$acc update device(rotationParams)

      ! ---- SRF enabled ----
      if (controlVariables % containsKey("srf enabled")) then
         rotationParams % srfEnabled = controlVariables % logicalValueForKey("srf enabled")
         !$acc update device(rotationParams)
      end if

      ! ---- periodic enabled ----
      if (controlVariables % containsKey("periodic enabled")) then
         rotationParams % periodicEnabled = .false.
         if (controlVariables % logicalValueForKey("periodic enabled")) then
            if (RotationAxisIndex() /= 0) then
               rotationParams % periodicEnabled = .true.
            end if
         end if
         !$acc update device(rotationParams)
      end if

   end subroutine InitializeRotationFromControl

   subroutine BuildRotationMatrix(R, angle, axis)
     !$acc routine seq
     ! Construct the 3x3 rotation matrix for a given rotation angle and a UNIT axis
     implicit none
     real(kind=RP), intent(out) :: R(3,3)
     real(kind=RP), intent(in)  :: angle
     real(kind=RP), intent(in)  :: axis(3)
     real(kind=RP) :: c, s, t
     real(kind=RP) :: x, y, z

     x = axis(1)
     y = axis(2)
     z = axis(3)

     c = cos(angle)
     s = sin(angle)
     t = 1.0_rp - c

     R(1,1) = t*x*x + c
     R(1,2) = t*x*y - s*z
     R(1,3) = t*x*z + s*y

     R(2,1) = t*x*y + s*z
     R(2,2) = t*y*y + c
     R(2,3) = t*y*z - s*x

     R(3,1) = t*x*z - s*y
     R(3,2) = t*y*z + s*x
     R(3,3) = t*z*z + c

   end subroutine BuildRotationMatrix

   function RotateVectorAroundAxis(v, angle, axis) result(vRot)
     !$acc routine seq
     implicit none
     real(kind=RP), intent(in) :: v(3), axis(3)
     real(kind=RP), intent(in) :: angle
     real(kind=RP)             :: vRot(3)
     real(kind=RP)             :: R(3,3)

     call BuildRotationMatrix(R, angle, axis)
     vRot = matmul(R, v)

   end function

   subroutine RotateStateMomentum(Q_in, Q_out, rotAngle)
     !$acc routine seq
!
!        Rotate only the momentum components of Q by rotAngle.
!        Default angle should be -periodicAngle (right-to-left frame).
!        Use +periodicAngle when rotating the left state for MPI side-2.
!
      implicit none
      real(kind=RP), intent(in)  :: Q_in(:)
      real(kind=RP), intent(out) :: Q_out(:)
      real(kind=RP), intent(in)  :: rotAngle
      integer, parameter :: IRHOU = 2, IRHOV = 3, IRHOW = 4

      Q_out = Q_in
      Q_out(IRHOU:IRHOW) = RotateVectorAroundAxis(Q_in(IRHOU:IRHOW), &
                                         rotAngle, &
                                          rotationParams % axis)
   end subroutine RotateStateMomentum

   subroutine RotateMomentumGradients(U_x_in, U_y_in, U_z_in, &
                                      U_x_out, U_y_out, U_z_out, rotAngle)
     !$acc routine seq
   !
   !  Rotate the gradient of the momentum field.
   !  Let G = ∇U be the 3×3 gradient tensor of momentum components:
   !        G = [ U_x  U_y  U_z ]
   !  where:
   !        U_x = ∂U/∂x,  U_y = ∂U/∂y,  U_z = ∂U/∂z
   !  Under a rotation R, the tensor transforms as:
   !        G' = R * G * R^T
   !  This is implemented in two steps:
   !    1) Rotate momentum components (rows):
   !           G ← R * G
   !       → applied independently to each column:
   !           U_x, U_y, U_z
   !    2) Rotate spatial directions (columns):
   !           G ← G * R^T
   !       → applied by rotating, for each variable k:
   !           [G(k,1), G(k,2), G(k,3)]
   !  The combination yields the full tensor transformation G' = R G R^T.
      implicit none

      real(kind=RP), intent(in)  :: U_x_in(:), U_y_in(:), U_z_in(:)
      real(kind=RP), intent(out) :: U_x_out(:), U_y_out(:), U_z_out(:)
      real(kind=RP), intent(in)  :: rotAngle

      real(kind=RP) :: ax(NDIM)
      real(kind=RP) :: G_x(size(U_x_in)), G_y(size(U_y_in)), G_z(size(U_z_in))
      real(kind=RP) :: v(NDIM), vr(NDIM)
      real(kind=RP) :: R(3,3)

      integer :: k, nGrad
      integer, parameter :: IRHOU = 2, IRHOV = 3, IRHOW = 4

      ax    = rotationParams % axis
      nGrad = size(U_x_in)

      call BuildRotationMatrix(R, rotAngle, ax)

      !------------------------------------------------------------
      ! Step 1: G ← R * G  (rotate momentum components / rows)
      !------------------------------------------------------------
      G_x = U_x_in
      G_y = U_y_in
      G_z = U_z_in

      G_x(2:4) = matmul(R, U_x_in(2:4))
      G_y(2:4) = matmul(R, U_y_in(2:4))
      G_z(2:4) = matmul(R, U_z_in(2:4))

      !------------------------------------------------------------
      ! Step 2: G ← G * R^T  (rotate spatial directions / columns)
      !------------------------------------------------------------
      do k = 1, nGrad
         v  = [G_x(k), G_y(k), G_z(k)]
         vr = matmul(R, v)
         U_x_out(k) = vr(1)
         U_y_out(k) = vr(2)
         U_z_out(k) = vr(3)
      end do

   end subroutine RotateMomentumGradients


end module RotationData
