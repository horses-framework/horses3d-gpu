#include "Includes.h"
module ShockCapturing

   use SMConstants,                only: RP, STD_OUT, LINE_LENGTH, NDIM, IX, IY, IZ
   use PhysicsStorage,             only: NCONS, NGRAD
   use FluidData,                  only: dimensionless
   use Utilities,                  only: toLower
   use DGSEMClass,                 only: ComputeTimeDerivative_f, DGSem
   use HexMeshClass,               only: HexMesh
   use ElementClass,               only: Element
#if !defined(SPALARTALMARAS)
   use LESModels,                  only: Smagorinsky_t, Smagorinsky_ComputeViscosity
#endif
   use SpectralVanishingViscosity, only: SVV, InitializeSVV
   use DGIntegrals,                only: ScalarWeakIntegrals, ScalarStrongIntegrals

   use ShockCapturingKeywords
   use SCsensorClass, only: SCsensor_t, Set_SCsensor, Destruct_SCsensor

   implicit none

   public :: Initialize_ShockCapturing
   public :: ShockCapturingDriver

   type ArtificialViscosity_t

      integer,  private :: region           !< Sensor region where it acts (1 or 2)
      integer,  private :: updateMethod     !< Method to compute the viscosity coefficient
      real(RP), private :: mu1              !< First viscosity parameter (low)
      real(RP), private :: alpha1           !< Second viscosity parameter (low)
      real(RP), private :: mu2              !< First viscosity parameter (high)
      real(RP), private :: alpha2           !< Second viscosity parameter (high)
      real(RP), private :: mu2alpha         !< Ratio alpha/mu
      logical,  private :: alphaIsPropToMu  !< .true. if alpha/mu is defined
#if !defined (SPALARTALMARAS)
      type(Smagorinsky_t), private :: Smagorinsky  !< For automatic viscosity
#endif

      contains

         procedure :: Initialize => AV_initialize
         procedure :: Viscosity  => AV_viscosity
         procedure :: Describe   => AV_describe

   end type ArtificialViscosity_t

   type SCdriver_t

      logical  :: isActive = .false.  !< On/Off flag

      type(SCsensor_t),             private              :: sensor   !< Sensor
      class(ArtificialViscosity_t), private, allocatable :: method1  !< First SC method
      class(ArtificialViscosity_t), private, allocatable :: method2  !< Second SC method

      contains

         procedure :: Detect           => SC_detect
         procedure :: ComputeViscosity => SC_viscosity
         procedure :: Describe         => SC_describe

         final :: SC_destruct

   end type SCdriver_t

   type, extends(ArtificialViscosity_t) :: SC_NoSVV_t

      integer,                                 private :: fluxType
      procedure(Viscous_Int), nopass, pointer, private :: ViscousFlux => null()

      contains

         procedure :: Initialize => NoSVV_initialize
         procedure :: Viscosity  => NoSVV_viscosity
         procedure :: Describe   => NoSVV_describe

         final :: NoSVV_destruct

   end type SC_NoSVV_t

   type, extends(ArtificialViscosity_t) :: SC_SVV_t

      real(RP), private :: sqrt_mu1
      real(RP), private :: sqrt_alpha1
      real(RP), private :: sqrt_mu2
      real(RP), private :: sqrt_alpha2
      real(RP), private :: sqrt_mu2alpha

      contains

         procedure :: Initialize => SVV_initialize
         procedure :: Viscosity  => SVV_viscosity
         procedure :: Describe   => SVV_describe

   end type SC_SVV_t
!
!  Interfaces
!  ----------
   abstract interface
      pure subroutine Viscous_Int(nEqn, nGradEqn, Q, Q_x, Q_y, Q_z, mu, beta, kappa, F)
         import RP, NDIM
         integer,       intent(in)  :: nEqn
         integer,       intent(in)  :: nGradEqn
         real(kind=RP), intent(in)  :: Q   (1:nEqn    )
         real(kind=RP), intent(in)  :: Q_x (1:nGradEqn)
         real(kind=RP), intent(in)  :: Q_y (1:nGradEqn)
         real(kind=RP), intent(in)  :: Q_z (1:nGradEqn)
         real(kind=RP), intent(in)  :: mu
         real(kind=RP), intent(in)  :: beta
         real(kind=RP), intent(in)  :: kappa
         real(kind=RP), intent(out) :: F(1:nEqn, 1:NDIM)
      end subroutine Viscous_Int
   end interface

   type(SCdriver_t), allocatable :: ShockCapturingDriver
!
!  =====================================================================
!  SCDEV_DISPATCH -- device-side shock capturing: START HERE
!  =====================================================================
!
!  The host configuration above is polymorphic (class(ArtificialViscosity_t)
!  method1/method2) and holds a procedure pointer (ViscousFlux). Neither can be
!  dispatched from an OpenACC compute region, which is why the artificial
!  viscosity was left unwired when the solver was ported to the GPU.
!
!  The fix mirrors the GRADVARS_DISPATCH convention already used for the
!  gradient variables (see NSGradientVariables_selector in
!  VariableConversion_NS.f90): flatten the configuration into plain module
!  scalars that live on the device, and replace the procedure pointer with a
!  "select case" inside an "!$acc routine seq".
!
!  Arrays are indexed by sensor REGION (1 = first method, 2 = second method),
!  matching the "region" argument threaded through AV_initialize.
!
!  These are written ONCE, by SC_SyncDeviceConfig, at the end of
!  Initialize_ShockCapturing. Everything here is configuration, not state --
!  the only per-step quantity is the sensor, which is updated separately
!  (see SC_UpdateSensorOnDevice).
!
!  TO ADD A NEW VISCOUS FLUX TYPE: add a case to SC_ArtificialViscousFlux_0D
!  and to NoSVV_initialize. To add a new viscosity update method, add a case to
!  SC_ElementViscosity. Callers need no changes.
!
!        grep -rn "SCDEV_DISPATCH" Solver/src
!  =====================================================================
!
   logical  :: SCdev_isActive        = .false.
   logical  :: SCdev_onDevice        = .true.   !< .false. if any method needs the host path
   logical  :: SCdev_hasMethod(2)    = .false.
   integer  :: SCdev_update(2)       = SC_CONST_ID
   integer  :: SCdev_fluxType(2)     = SC_PHYS_ID
   real(RP) :: SCdev_mu1(2)          = 0.0_RP
   real(RP) :: SCdev_mu2(2)          = 0.0_RP
   real(RP) :: SCdev_smagC(2)        = 0.0_RP
   integer  :: SCdev_smagWallModel(2) = 0
!
!  copyin, not create: "create" would reserve the memory without initialising
!  it, and Initialize_ShockCapturing is not reached at all for an Euler run
!  (it sits inside "if (flowIsNavierStokes)"). copyin seeds the device with the
!  defaults above, so the guards below are always defined; SC_SyncDeviceConfig
!  then overwrites them when shock capturing is actually configured. Same
!  pattern as grad_vars in PhysicsStorage_NS.
!
   !$acc declare copyin(SCdev_isActive, SCdev_onDevice)
   !$acc declare copyin(SCdev_hasMethod, SCdev_update, SCdev_fluxType)
   !$acc declare copyin(SCdev_mu1, SCdev_mu2, SCdev_smagC, SCdev_smagWallModel)

   public :: SCdev_isActive, SCdev_onDevice
   public :: SC_ComputeElementAviscFlux, SC_ProlongAviscFluxToFaces
   public :: SC_UpdateSensorOnDevice
!
!  ========
   contains
!  ========
!
!
!///////////////////////////////////////////////////////////////////////////////
!
!     Initializer
!
!///////////////////////////////////////////////////////////////////////////////
!
   subroutine Initialize_ShockCapturing(self, controlVariables, sem, &
                                        TimeDerivative, TimeDerivativeIsolated)
!
!     -------
!     Modules
!     -------
      use FTValueDictionaryClass
!
!     ---------
!     Interface
!     ---------
      implicit none
      type(SCdriver_t), allocatable, intent(inout) :: self
      class(FTValueDictionary),      intent(in)    :: controlVariables
      class(DGSem),                  intent(inout) :: sem
      procedure(ComputeTimeDerivative_f)           :: TimeDerivative
      procedure(ComputeTimeDerivative_f)           :: TimeDerivativeIsolated
!
!     ---------------
!     Local variables
!     ---------------
      character(len=:), allocatable :: method
      integer                       :: minSteps

!
!     Check if shock-capturing is requested
!     -------------------------------------
      allocate(self)

      if (controlVariables % containsKey(SC_KEY)) then
         self % isActive = controlVariables % logicalValueForKey(SC_KEY)
      else
         self % isActive = .false.
      end if

      if (.not. self % isActive) then
!
!        Still has to reach the device: "!$acc declare create" reserves the
!        memory but does not initialise it, so without this the device copy of
!        SCdev_isActive is undefined and the guards in the time-derivative
!        kernels would branch on garbage.
!        ---------------------------------------------------------------------
         call SC_SyncDeviceConfig(self)
         return
      end if
!
!     Shock-capturing methods
!     -----------------------
      if (controlVariables % containsKey(SC_METHOD1_KEY)) then
         method = controlVariables % stringValueForKey(SC_METHOD1_KEY, LINE_LENGTH)
      else
         method = SC_NO_VAL
      end if
      call toLower(method)

      select case (trim(method))
      case (SC_NOSVV_VAL)
         allocate(SC_NoSVV_t :: self % method1)

      case (SC_SVV_VAL)
         allocate(SC_SVV_t :: self % method1)

      case (SC_NO_VAL)
         safedeallocate(self % method1)

      case default
         write(STD_OUT,*) 'ERROR. Unavailable first shock-capturing method. Options are:'
         write(STD_OUT,*) '   * ', SC_NO_VAL
         write(STD_OUT,*) '   * ', SC_NOSVV_VAL
         write(STD_OUT,*) '   * ', SC_SVV_VAL
         error stop

      end select

      if (controlVariables % containsKey(SC_METHOD2_KEY)) then
         method = controlVariables % stringValueForKey(SC_METHOD2_KEY, LINE_LENGTH)
      else
         method = SC_NO_VAL
      end if
      call toLower(method)

      select case (trim(method))
      case (SC_NOSVV_VAL)
         allocate(SC_NoSVV_t :: self % method2)

      case (SC_SVV_VAL)
         allocate(SC_SVV_t :: self % method2)

      case (SC_NO_VAL)
         safedeallocate(self % method2)

      case default
         write(STD_OUT,*) 'ERROR. Unavailable second shock-capturing method. Options are:'
         write(STD_OUT,*) '   * ', SC_NO_VAL
         write(STD_OUT,*) '   * ', SC_NOSVV_VAL
         write(STD_OUT,*) '   * ', SC_SVV_VAL
         error stop

      end select
!
!     Initialize viscous and hyperbolic terms
!     ---------------------------------------
      if (allocated(self % method1)) call self % method1 % Initialize(controlVariables, sem % mesh, 1)
      if (allocated(self % method2)) call self % method2 % Initialize(controlVariables, sem % mesh, 2)
!
!     Sensor 'inertia'
!     ----------------
      if (controlVariables % containsKey(SC_SENSOR_INERTIA_KEY)) then
         minSteps = controlVariables % realValueForKey(SC_SENSOR_INERTIA_KEY)
         if (minSteps < 1) then
            write(STD_OUT,*) 'ERROR. Sensor inertia must be at least 1.'
            error stop
         end if
      else
         minSteps = 1
      end if
!
!     Construct sensor
!     ----------------
      call Set_SCsensor(self % sensor, controlVariables, sem, minSteps, &
                        TimeDerivative, TimeDerivativeIsolated)
!
!     Flatten the configuration for the device
!     ----------------------------------------
      call SC_SyncDeviceConfig(self)

   end subroutine Initialize_ShockCapturing
!
!///////////////////////////////////////////////////////////////////////////////
!
   subroutine SC_SyncDeviceConfig(self)
!
!     ---------------------------------------------------------------------
!     SCDEV_DISPATCH -- copy the polymorphic host configuration into the flat
!     module scalars the device kernels read, then push them to the device.
!
!     Called once, from Initialize_ShockCapturing. See the SCDEV_DISPATCH
!     block at the top of this file for the convention.
!     ---------------------------------------------------------------------
!
      implicit none
      type(SCdriver_t), intent(in) :: self

      SCdev_isActive = self % isActive

      if (self % isActive) then
         if (allocated(self % method1)) call SC_FlattenMethod(self % method1, 1)
         if (allocated(self % method2)) call SC_FlattenMethod(self % method2, 2)
      end if

      !$acc update device(SCdev_isActive, SCdev_onDevice)
      !$acc update device(SCdev_hasMethod, SCdev_update, SCdev_fluxType)
      !$acc update device(SCdev_mu1, SCdev_mu2, SCdev_smagC, SCdev_smagWallModel)

   end subroutine SC_SyncDeviceConfig
!
!///////////////////////////////////////////////////////////////////////////////
!
   subroutine SC_FlattenMethod(method, region)
!
!     ---------------------------------------------------------------------
!     SCDEV_DISPATCH -- flatten one method into the region-indexed scalars.
!
!     SVV (filtered) needs a spectral filter matrix per element and was never
!     ported, so it cannot be flattened. That is not an error -- it just means
!     this run keeps using the original host routine (SC_viscosity); we record
!     that by clearing SCdev_onDevice, which TimeDerivative_ComputeArtificialViscosity
!     reads to pick a path.
!     ---------------------------------------------------------------------
!
      implicit none
      class(ArtificialViscosity_t), intent(in) :: method
      integer,                      intent(in) :: region

      select type (m => method)

      type is (SC_NoSVV_t)
         SCdev_hasMethod(region) = .true.
         SCdev_update(region)    = m % updateMethod
         SCdev_fluxType(region)  = m % fluxType
         SCdev_mu1(region)       = m % mu1
         SCdev_mu2(region)       = m % mu2
#if !defined (SPALARTALMARAS)
         SCdev_smagC(region)         = m % Smagorinsky % C
         SCdev_smagWallModel(region) = m % Smagorinsky % WallModel
#endif

      class default
         SCdev_hasMethod(region) = .true.
         SCdev_onDevice          = .false.

      end select

   end subroutine SC_FlattenMethod
!
!///////////////////////////////////////////////////////////////////////////////
!
   subroutine SC_UpdateSensorOnDevice(mesh)
!
!     ---------------------------------------------------------------------
!     SCDEV_DISPATCH -- push the freshly-computed sensor to the device.
!
!     The sensor is computed on the host (it clusters a scalar per element, so
!     a host round-trip is cheap compared with porting the clustering). The
!     element storage was copied to the device once at start-up, so without
!     this update every device kernel would read the sensor value from
!     initialisation forever.
!     ---------------------------------------------------------------------
!
      implicit none
      type(HexMesh), intent(in) :: mesh
      integer :: eID

      if (.not. SCdev_isActive) return

      do eID = 1, size(mesh % elements)
         !$acc update device(mesh % elements(eID) % storage % sensor) async(1)
      end do
      !$acc wait(1)

   end subroutine SC_UpdateSensorOnDevice
!
!///////////////////////////////////////////////////////////////////////////////
!
!     Base class
!
!///////////////////////////////////////////////////////////////////////////////
!
   subroutine SC_detect(self, sem, t)
!
!     ---------
!     Interface
!     ---------
      implicit none
      class(SCdriver_t), intent(inout) :: self
      type(DGSem),       intent(inout) :: sem
      real(RP),          intent(in)    :: t

!
!     The sensor runs on the host, but between solution files Q only exists on
!     the device -- the host copy is whatever was left over from the last save.
!     Without this pull the sensor clusters a stale solution and then gets
!     written out next to a fresh Q as if the two were contemporaneous.
!     ------------------------------------------------------------------------
      call sem % mesh % UpdateHostData()

      call self % sensor % Compute(sem, t)
!
!     Push the result back: the element storage was copied to the device once
!     at start-up, so the artificial-viscosity kernels would otherwise keep
!     reading the sensor value from initialisation.
!     ------------------------------------------------------------------------
      call SC_UpdateSensorOnDevice(sem % mesh)

   end subroutine SC_detect
!
!///////////////////////////////////////////////////////////////////////////////
!
   subroutine SC_viscosity(self, mesh, e, SCflux)
!
!     ---------
!     Interface
!     ---------
      implicit none
      class(SCdriver_t), intent(in)    :: self
      type(HexMesh),     intent(inout) :: mesh
      type(Element),     intent(inout) :: e
      real(RP),          intent(out)   :: SCflux(1:NCONS,       &
                                                 0:e % Nxyz(1), &
                                                 0:e % Nxyz(2), &
                                                 0:e % Nxyz(3), &
                                                 1:NDIM)
!
!     ---------------
!     Local variables
!     ---------------
      real(RP) :: switch
      logical  :: updated


      switch = e % storage % sensor
      updated = .false.

      if (switch >= 1.0_RP) then
         if (allocated(self % method2)) then
            call self % method2 % Viscosity(mesh, e, switch, SCflux)
            updated = .true.
         end if

      elseif (switch > 0.0_RP) then
         if (allocated(self % method1)) then
            call self % method1 % Viscosity(mesh, e, switch, SCflux)
            updated = .true.
         end if

      end if

      if (.not. updated) then
         SCflux = 0.0_RP
         e % storage % artificialDiss = 0.0_RP
      end if

   end subroutine SC_viscosity
!
!///////////////////////////////////////////////////////////////////////////////
!
   pure subroutine SC_ArtificialViscousFlux_0D(region, Q, Q_x, Q_y, Q_z, mu, beta, kappa, F)
!
!     ---------------------------------------------------------------------
!     SCDEV_DISPATCH -- device-side replacement for the "self % ViscousFlux"
!     procedure pointer of SC_NoSVV_t.
!
!     A procedure pointer cannot be called from an OpenACC compute region, so
!     the same choice is made here with a "select case" on the flattened
!     SCdev_fluxType. The Physical branch delegates to ViscousFlux_selector_0D
!     so the artificial viscosity automatically follows whichever gradient
!     variables the run uses (see GRADVARS_DISPATCH).
!     ---------------------------------------------------------------------
!
      !$acc routine seq
      use Physics, only: GuermondPopovFlux_ENTROPY, ViscousFlux_selector_0D
      implicit none
      integer,       intent(in)  :: region
      real(kind=RP), intent(in)  :: Q   (1:NCONS)
      real(kind=RP), intent(in)  :: Q_x (1:NGRAD)
      real(kind=RP), intent(in)  :: Q_y (1:NGRAD)
      real(kind=RP), intent(in)  :: Q_z (1:NGRAD)
      real(kind=RP), intent(in)  :: mu, beta, kappa
      real(kind=RP), intent(out) :: F   (1:NCONS, 1:NDIM)

      select case (SCdev_fluxType(region))
      case (SC_GP_ID)
         call GuermondPopovFlux_ENTROPY(NCONS, NGRAD, Q, Q_x, Q_y, Q_z, mu, beta, kappa, F)
      case default
         call ViscousFlux_selector_0D(NCONS, NGRAD, Q, Q_x, Q_y, Q_z, mu, beta, kappa, F)
      end select

   end subroutine SC_ArtificialViscousFlux_0D
!
!///////////////////////////////////////////////////////////////////////////////
!
   subroutine SC_ComputeElementAviscFlux(e)
!
!     ---------------------------------------------------------------------
!     SCDEV_DISPATCH -- artificial viscous flux of one element, in
!     contravariant form, written into e % storage % AviscContravariantFlux.
!
!     This is the device-ready equivalent of SC_viscosity + NoSVV_viscosity:
!     same mathematics, but the region is chosen with plain integers instead
!     of allocatable polymorphic components, and the flux with a select case
!     instead of a procedure pointer. Unlike its host counterpart it does NOT
!     prolong to the faces -- see SC_ProlongAviscFluxToFaces.
!
!     Both paths must stay in step: SC_viscosity is still what the isolated
!     time derivative (and hence the sensor) uses.
!     ---------------------------------------------------------------------
!
      !$acc routine vector
#if !defined (SPALARTALMARAS)
      use LESModels, only: Smagorinsky_ComputeViscosity
#endif
      implicit none
      type(Element), intent(inout) :: e
!
!     ---------------
!     Local variables
!     ---------------
      integer       :: i, j, k, eq, region
      real(kind=RP) :: switch, mu, kappa, delta
      real(kind=RP) :: covariantFlux(1:NCONS, 1:NDIM)

      switch = e % storage % sensor
!
!     Pick the region exactly as SC_viscosity does
!     --------------------------------------------
      if (switch >= 1.0_RP .and. SCdev_hasMethod(2)) then
         region = 2
      elseif (switch > 0.0_RP .and. SCdev_hasMethod(1)) then
         region = 1
      else
         region = 0
      end if

      if (region == 0) then

         !$acc loop vector collapse(3)
         do k = 0, e % Nxyz(3) ; do j = 0, e % Nxyz(2) ; do i = 0, e % Nxyz(1)
            !$acc loop seq
            do eq = 1, NCONS
               e % storage % AviscContravariantFlux(eq,i,j,k,IX) = 0.0_RP
               e % storage % AviscContravariantFlux(eq,i,j,k,IY) = 0.0_RP
               e % storage % AviscContravariantFlux(eq,i,j,k,IZ) = 0.0_RP
            end do
         end do                ; end do                ; end do
         e % storage % artificialDiss = 0.0_RP

      else

!        Written out rather than product(e % Nxyz + 1): array-valued intrinsics
!        on device routines are a portability hazard, and this is three terms.
         delta = ( e % geom % Volume                                            &
                 / real((e % Nxyz(1)+1)*(e % Nxyz(2)+1)*(e % Nxyz(3)+1), RP) )  &
                 ** (1.0_RP / 3.0_RP)
!
!     Viscosity and flux, node by node
!     --------------------------------
      !$acc loop vector collapse(3) private(covariantFlux, mu, kappa)
      do k = 0, e % Nxyz(3) ; do j = 0, e % Nxyz(2) ; do i = 0, e % Nxyz(1)

         select case (SCdev_update(region))
         case (SC_SENSOR_ID)
            mu = (SCdev_mu1(region) * (1.0_RP-switch) + SCdev_mu2(region) * switch) * e % hn

#if !defined (SPALARTALMARAS)
         case (SC_SMAG_ID)
            call Smagorinsky_ComputeViscosity(delta, e % geom % dWall(i,j,k), &
                                              e % storage % Q(:,i,j,k),       &
                                              e % storage % U_x(:,i,j,k),     &
                                              e % storage % U_y(:,i,j,k),     &
                                              e % storage % U_z(:,i,j,k),     &
                                              mu,                             &
                                              SCdev_smagC(region),            &
                                              SCdev_smagWallModel(region))
#endif

         case default   ! SC_CONST_ID
            if (switch >= 1.0_RP) then
               mu = SCdev_mu2(region) * e % hn
            else
               mu = SCdev_mu1(region) * e % hn
            end if

         end select

         kappa = dimensionless % mu_to_kappa * mu

         call SC_ArtificialViscousFlux_0D(region, e % storage % Q(:,i,j,k),   &
                                          e % storage % U_x(:,i,j,k),         &
                                          e % storage % U_y(:,i,j,k),         &
                                          e % storage % U_z(:,i,j,k),         &
                                          mu, 0.0_RP, kappa, covariantFlux)

         !$acc loop seq
         do eq = 1, NCONS
            e % storage % AviscContravariantFlux(eq,i,j,k,IX) =                       &
                 covariantFlux(eq,IX) * e % geom % jGradXi(IX,i,j,k)                  &
               + covariantFlux(eq,IY) * e % geom % jGradXi(IY,i,j,k)                  &
               + covariantFlux(eq,IZ) * e % geom % jGradXi(IZ,i,j,k)

            e % storage % AviscContravariantFlux(eq,i,j,k,IY) =                       &
                 covariantFlux(eq,IX) * e % geom % jGradEta(IX,i,j,k)                 &
               + covariantFlux(eq,IY) * e % geom % jGradEta(IY,i,j,k)                 &
               + covariantFlux(eq,IZ) * e % geom % jGradEta(IZ,i,j,k)

            e % storage % AviscContravariantFlux(eq,i,j,k,IZ) =                       &
                 covariantFlux(eq,IX) * e % geom % jGradZeta(IX,i,j,k)                &
               + covariantFlux(eq,IY) * e % geom % jGradZeta(IY,i,j,k)                &
               + covariantFlux(eq,IZ) * e % geom % jGradZeta(IZ,i,j,k)
         end do

      end do                ; end do                ; end do

      end if

   end subroutine SC_ComputeElementAviscFlux
!
!///////////////////////////////////////////////////////////////////////////////
!
   subroutine SC_ProlongAviscFluxToFaces(mesh)
!
!     ---------------------------------------------------------------------
!     SCDEV_DISPATCH -- send each element's artificial viscous flux to its six
!     faces, so the interface terms can average the two sides.
!
!     TODO(device): this is still a HOST loop. The interpolation it delegates
!     to (Face_AdaptAviscFluxToFace) reads the Tset projection matrices, which
!     are not in the device data region, so it cannot be given an
!     "!$acc routine" as it stands. Everything upstream of it
!     (SC_ComputeElementAviscFlux) is already device-ready, so porting this is
!     the one remaining step for a fully resident artificial viscosity --
!     mirror HexElement_ProlongSolToFaces / _GL, which solve exactly this
!     problem for the solution.
!
!     Until then the AviscContravariantFlux must be on the host when this
!     runs, and the face AviscFlux must be pushed back to the device after.
!     ---------------------------------------------------------------------
!
      implicit none
      type(HexMesh), intent(inout) :: mesh
      integer :: eID, fIDs(6)

      do eID = 1, size(mesh % elements)
         fIDs = mesh % elements(eID) % faceIDs
         call mesh % elements(eID) % ProlongAviscFluxToFaces(NCONS,                           &
                        mesh % elements(eID) % storage % AviscContravariantFlux,              &
                        mesh % faces(fIDs(1)), mesh % faces(fIDs(2)), mesh % faces(fIDs(3)),  &
                        mesh % faces(fIDs(4)), mesh % faces(fIDs(5)), mesh % faces(fIDs(6)))
      end do

   end subroutine SC_ProlongAviscFluxToFaces
!
!///////////////////////////////////////////////////////////////////////////////
!
   subroutine SC_describe(self)
!
!     -------
!     Modules
!     -------
      use MPI_Process_Info, only: MPI_Process
      use Headers,          only: Subsection_Header
!
!     ---------
!     Interface
!     ---------
      implicit none
      class(SCdriver_t), intent(in) :: self


      if (.not. MPI_Process % isRoot .or. .not. self % isActive) return

      write(STD_OUT, "(/)")
      call Subsection_Header("Shock-Capturing")

      call self % sensor % Describe()

      if (allocated(self % method1)) then
         write(STD_OUT,*) ""
         write(STD_OUT,"(30X,A,A30)") "=>", "First method"
         call self % method1 % Describe()
      end if

      if (allocated(self % method2)) then
         write(STD_OUT,*) ""
         write(STD_OUT,"(30X,A,A30)") "=>", "Second method"
         call self % method2 % Describe()
      end if

   end subroutine SC_describe
!
!///////////////////////////////////////////////////////////////////////////////
!
   subroutine SC_destruct(self)
!
!     ---------
!     Interface
!     ---------
      implicit none
      type(SCdriver_t), intent(inout) :: self


      self%isActive = .false.

      call Destruct_SCsensor(self % sensor)

      if (allocated(self % method1)) deallocate(self % method1)
      if (allocated(self % method2)) deallocate(self % method2)

   end subroutine SC_destruct
!
!///////////////////////////////////////////////////////////////////////////////
!
!     Shock-capturing methods base class
!
!///////////////////////////////////////////////////////////////////////////////
!
   subroutine AV_initialize(self, controlVariables, mesh, region)
!
!     -------
!     Modules
!     -------
      use FTValueDictionaryClass
!
!     ---------
!     Interface
!     ---------
      implicit none
      class(ArtificialViscosity_t), intent(inout) :: self
      type(FTValueDictionary),      intent(in)    :: controlVariables
      type(HexMesh),                intent(inout) :: mesh
      integer,                      intent(in)    :: region
!
!     ---------------
!     Local variables
!     ---------------
      character(len=:), allocatable :: update

!
!     Viscosity values (mu and alpha)
!     -------------------------------
      if (region == 1) then

         if (controlVariables % containsKey(SC_MU1_KEY)) then
            self % mu1 = controlVariables % doublePrecisionValueForKey(SC_MU1_KEY)
         else
            self % mu1 = 0.0_RP
         end if

         if (controlVariables % containsKey(SC_MU2_KEY)) then
            self % mu2 = controlVariables % doublePrecisionValueForKey(SC_MU2_KEY)
         else
            self % mu2 = self % mu1
         end if

      else

         if (controlVariables % containsKey(SC_MU2_KEY)) then
            self % mu2 = controlVariables % doublePrecisionValueForKey(SC_MU2_KEY)
         else
            self % mu2 = 0.0_RP
         end if
         self % mu1 = self % mu2

      end if

      if (controlVariables % containsKey(SC_ALPHA_MU_KEY)) then

         self % alphaIsPropToMu = .true.
         self % mu2alpha        = controlVariables % doublePrecisionValueForKey(SC_ALPHA_MU_KEY)
         self % alpha1          = self % mu2alpha * self % mu1
         self % alpha2          = self % mu2alpha * self % mu2

      else

         self % alphaIsPropToMu = .false.

         if (region == 1) then

            if (controlVariables % containsKey(SC_ALPHA1_KEY)) then
               self % alpha1 = controlVariables % doublePrecisionValueForKey(SC_ALPHA1_KEY)
            else
               self % alpha1 = 0.0_RP
            end if

            if (controlVariables % containsKey(SC_ALPHA2_KEY)) then
               self % alpha2 = controlVariables % doublePrecisionValueForKey(SC_ALPHA2_KEY)
            else
               self % alpha2 = self % alpha1
            end if

         else

            if (controlVariables % containsKey(SC_ALPHA2_KEY)) then
               self % alpha2 = controlVariables % doublePrecisionValueForKey(SC_ALPHA2_KEY)
            else
               self % alpha2 = 0.0_RP
            end if
            self % alpha1 = self % alpha2

         end if

      end if
!
!     Viscosity update method
!     -----------------------
      if (region == 1) then

         if (controlVariables % containsKey(SC_UPDATE_KEY)) then
            update = controlVariables % StringValueForKey(SC_UPDATE_KEY, LINE_LENGTH)
            call toLower(update)

            select case (trim(update))
            case (SC_CONST_VAL)
               self % updateMethod = SC_CONST_ID

            case (SC_SENSOR_VAL)
               self % updateMethod = SC_SENSOR_ID

#if !defined (SPALARTALMARAS)
            case (SC_SMAG_VAL)

               self % updateMethod = SC_SMAG_ID
               if (.not. self % alphaIsPropToMu) then
                  write(STD_OUT,*) 'ERROR. Alpha must be proportional to mu when using shock-capturing with LES.'
                  error stop
               end if

               ! TODO: Use the default constructor
               self % Smagorinsky % active = .true.
               self % Smagorinsky % requiresWallDistances = .false.
               self % Smagorinsky % WallModel = 0  ! No wall model
               self % Smagorinsky % C = self % mu1
#endif

            case default
               write(STD_OUT,*) 'ERROR. Unavailable shock-capturing update strategy. Options are:'
               write(STD_OUT,*) '   * ', SC_CONST_VAL
               write(STD_OUT,*) '   * ', SC_SENSOR_VAL
#if !defined (SPALARTALMARAS)
               write(STD_OUT,*) '   * ', SC_SMAG_VAL
#endif
               error stop

            end select

         else
            self % updateMethod = SC_CONST_ID

         end if

      else

         self % updateMethod = SC_CONST_ID

      end if

   end subroutine AV_initialize
!
!///////////////////////////////////////////////////////////////////////////////
!
   subroutine AV_viscosity(self, mesh, e, switch, SCflux)
!
!     ---------
!     Interface
!     ---------
      implicit none
      class(ArtificialViscosity_t), intent(in)    :: self
      type(HexMesh),                intent(inout) :: mesh
      type(Element),                intent(inout) :: e
      real(RP),                     intent(in)    :: switch
      real(RP),                     intent(out)   :: SCflux(1:NCONS,       &
                                                            0:e % Nxyz(1), &
                                                            0:e % Nxyz(2), &
                                                            0:e % Nxyz(3), &
                                                            1:NDIM)


      SCflux = 0.0_RP

   end subroutine AV_viscosity
!
!///////////////////////////////////////////////////////////////////////////////
!
   subroutine AV_describe(self)
!
!     ---------
!     Interface
!     ---------
      implicit none
      class(ArtificialViscosity_t), intent(in) :: self

   end subroutine AV_describe
!
!///////////////////////////////////////////////////////////////////////////////
!
!     Simple artificial viscosity
!
!///////////////////////////////////////////////////////////////////////////////
!
   subroutine NoSVV_initialize(self, controlVariables, mesh, region)
!
!     -------
!     Modules
!     -------
      use FTValueDictionaryClass
      use PhysicsStorage, only: grad_vars, GRADVARS_STATE, &
                                GRADVARS_ENTROPY, GRADVARS_ENERGY
      use Physics,        only: ViscousFlux_STATE, ViscousFlux_ENTROPY, &
                                ViscousFlux_ENERGY, GuermondPopovFlux_ENTROPY
!
!     ---------
!     Interface
!     ---------
      implicit none
      class(SC_NoSVV_t),       intent(inout) :: self
      type(FTValueDictionary), intent(in)    :: controlVariables
      type(HexMesh),           intent(inout) :: mesh
      integer,                 intent(in)    :: region
!
!     ---------------
!     Local variables
!     ---------------
      character(len=:), allocatable :: flux

!
!     Parent initializer
!     ------------------
      call self % ArtificialViscosity_t % Initialize(controlVariables, mesh, region)
      self % region = region
!
!     Set the flux type
!     -----------------
      if (region == 1 .and. controlVariables % containsKey(SC_VISC_FLUX1_KEY)) then
         flux = controlVariables % stringValueForKey(SC_VISC_FLUX1_KEY, LINE_LENGTH)
      elseif (region == 2 .and. controlVariables % containsKey(SC_VISC_FLUX2_KEY)) then
         flux = controlVariables % stringValueForKey(SC_VISC_FLUX2_KEY, LINE_LENGTH)
      else
         flux = SC_PHYS_VAL
      end if

      call toLower(flux)

      select case (trim(flux))
      case (SC_PHYS_VAL); self % fluxType = SC_PHYS_ID
      case (SC_GP_VAL);   self % fluxType = SC_GP_ID
      case default
         write(STD_OUT,'(A,I1,A)') 'ERROR. Artificial viscosity ', region, &
                                   ' not recognized. Options are:'
         write(STD_OUT,*) '   * ', SC_PHYS_VAL
         write(STD_OUT,*) '   * ', SC_GP_VAL
         error stop
      end select

      select case (self % fluxType)
      case (SC_PHYS_ID)
         select case (grad_vars)
         case (GRADVARS_STATE);   self % ViscousFlux => ViscousFlux_STATE
         case (GRADVARS_ENTROPY); self % ViscousFlux => ViscousFlux_ENTROPY
         case (GRADVARS_ENERGY);  self % ViscousFlux => ViscousFlux_ENERGY
         end select

      case (SC_GP_ID)
         select case (grad_vars)
         case (GRADVARS_ENTROPY); self % ViscousFlux => GuermondPopovFlux_ENTROPY
         case default
            write(STD_OUT,*) "ERROR. Guermond-Popov (2014) artificial ",  &
                              "viscosity is only configured for Entropy ", &
                              "gradient variables"
            error stop
         end select

      end select

   end subroutine NoSVV_initialize
!
!///////////////////////////////////////////////////////////////////////////////
!
   subroutine NoSVV_viscosity(self, mesh, e, switch, SCflux)
!
!     --------------------------------------------------------------------------
!     TODO: Introduce alpha viscosity, which probably means reimplementing here
!           all the viscous fluxes of `Physics_NS`...
!     --------------------------------------------------------------------------
!
!     ---------
!     Interface
!     ---------
      implicit none
      class(SC_NoSVV_t), intent(in)    :: self
      type(HexMesh),     intent(inout) :: mesh
      type(Element),     intent(inout) :: e
      real(RP),          intent(in)    :: switch
      real(RP),          intent(out)   :: SCflux(1:NCONS,       &
                                                 0:e % Nxyz(1), &
                                                 0:e % Nxyz(2), &
                                                 0:e % Nxyz(3), &
                                                 1:NDIM)
!
!     ---------------
!     Local variables
!     ---------------
      integer  :: i
      integer  :: j
      integer  :: k
      integer  :: fIDs(6)
      real(RP) :: delta
      real(RP) :: kappa
      real(RP) :: mu(0:e % Nxyz(1), 0:e % Nxyz(2), 0:e % Nxyz(3))
      real(RP) :: covariantFlux(1:NCONS, 1:NDIM)


      if (switch > 0.0_RP) then
!
!        Compute viscosity
!        -----------------
         select case (self % updateMethod)
         case (SC_CONST_ID)
            mu = merge(self % mu2, self % mu1, switch >= 1.0_RP) * e % hn

         case (SC_SENSOR_ID)
            mu = (self % mu1 * (1.0_RP-switch) + self % mu2 * switch) * e % hn

#if !defined (SPALARTALMARAS)
         case (SC_SMAG_ID)

            delta = (e % geom % Volume / product(e % Nxyz + 1)) ** (1.0_RP / 3.0_RP)
            do k = 0, e % Nxyz(3) ; do j = 0, e % Nxyz(2) ; do i = 0, e % Nxyz(1)
               call Smagorinsky_ComputeViscosity(delta, e % geom % dWall(i,j,k), &
                                                 e % storage % Q(:,i,j,k),       &
                                                 e % storage % U_x(:,i,j,k),     &
                                                 e % storage % U_y(:,i,j,k),     &
                                                 e % storage % U_z(:,i,j,k),     &
                                                 mu(i,j,k),                      &
                                                 self % smagorinsky % C,            &
                                                 self % smagorinsky % WallModel)
            end do                ; end do                ; end do
#endif

         end select
!
!        Compute the viscous flux
!        ------------------------
         do k = 0, e % Nxyz(3) ; do j = 0, e % Nxyz(2) ; do i = 0, e % Nxyz(1)

            kappa = dimensionless % mu_to_kappa * mu(i,j,k)
            call self % ViscousFlux(NCONS, NGRAD, e % storage % Q(:,i,j,k), &
                                    e % storage % U_x(:,i,j,k),             &
                                    e % storage % U_y(:,i,j,k),             &
                                    e % storage % U_z(:,i,j,k),             &
                                    mu(i,j,k), 0.0_RP, kappa,               &
                                    covariantflux)

            SCflux(:,i,j,k,IX) = covariantFlux(:,IX) * e % geom % jGradXi(IX,i,j,k) &
                               + covariantFlux(:,IY) * e % geom % jGradXi(IY,i,j,k) &
                               + covariantFlux(:,IZ) * e % geom % jGradXi(IZ,i,j,k)


            SCflux(:,i,j,k,IY) = covariantFlux(:,IX) * e % geom % jGradEta(IX,i,j,k) &
                               + covariantFlux(:,IY) * e % geom % jGradEta(IY,i,j,k) &
                               + covariantFlux(:,IZ) * e % geom % jGradEta(IZ,i,j,k)


            SCflux(:,i,j,k,IZ) = covariantFlux(:,IX) * e % geom % jGradZeta(IX,i,j,k) &
                               + covariantFlux(:,IY) * e % geom % jGradZeta(IY,i,j,k) &
                               + covariantFlux(:,IZ) * e % geom % jGradZeta(IZ,i,j,k)

         end do                ; end do                ; end do

      else

         e % storage % artificialDiss = 0.0_RP
         SCflux = 0.0_RP

      end if
!
!     Project to faces
!     ----------------
      fIDs = e % faceIDs
      call e % ProlongAviscFluxToFaces(NCONS, SCflux, mesh % faces(fIDs(1)), &
                                                      mesh % faces(fIDs(2)), &
                                                      mesh % faces(fIDs(3)), &
                                                      mesh % faces(fIDs(4)), &
                                                      mesh % faces(fIDs(5)), &
                                                      mesh % faces(fIDs(6))  )

   end subroutine NoSVV_viscosity
!
!///////////////////////////////////////////////////////////////////////////////
!
   subroutine NoSVV_describe(self)
!
!     -------
!     Modules
!     -------
      use MPI_Process_Info, only: MPI_Process
!
!     ---------
!     Interface
!     ---------
      implicit none
      class(SC_NoSVV_t), intent(in) :: self


      if (.not. MPI_Process % isRoot) return

      write(STD_OUT,"(30X,A,A30)", advance="no") "->", "Viscosity update method: "
      select case (self % updateMethod)
         case (SC_CONST_ID);  write(STD_OUT,"(A)") SC_CONST_VAL
         case (SC_SENSOR_ID); write(STD_OUT,"(A)") SC_SENSOR_VAL
         case (SC_SMAG_ID);   write(STD_OUT,"(A)") SC_SMAG_VAL
      end select

      if (self % updateMethod == SC_SMAG_ID) then
#if !defined (SPALARTALMARAS)
         write(STD_OUT,"(30X,A,A30,1pG10.3)") "->", "LES intensity (CS): ", self % Smagorinsky % C
#endif
      else
         if (self % region == 1) then
            write(STD_OUT,"(30X,A,A30,1pG10.3)") "->","Mu viscosity 1: ", self % mu1
            write(STD_OUT,"(30X,A,A30,1pG10.3)") "->","Mu viscosity 2: ", self % mu2
         else ! self % region == 2
            write(STD_OUT,"(30X,A,A30,1pG10.3)") "->","Mu viscosity: ", self % mu2
         end if
      end if

      if (self % alphaIsPropToMu) then
         write(STD_OUT,"(30X,A,A30,1pG10.3,A)") "->", "Alpha viscosity: ", self % mu2alpha, "x mu"
      else
         if (self % region == 1) then
            write(STD_OUT,"(30X,A,A30,1pG10.3)") "->","Alpha viscosity 1: ", self % alpha1
            write(STD_OUT,"(30X,A,A30,1pG10.3)") "->","Alpha viscosity 2: ", self % alpha2
         else ! self % region == 2
            write(STD_OUT,"(30X,A,A30,1pG10.3)") "->","Alpha viscosity: ", self % alpha2
         end if
      end if

      write(STD_OUT,"(30X,A,A30)", advance="no") "->", "Dissipation type: "
      select case (self % fluxType)
         case (SC_PHYS_ID); write(STD_OUT,"(A)") SC_PHYS_VAL
         case (SC_GP_ID);   write(STD_OUT,"(A)") SC_GP_VAL
      end select

   end subroutine NoSVV_describe
!
!///////////////////////////////////////////////////////////////////////////////
!
   pure subroutine NoSVV_destruct(self)
!
!     ---------
!     Interface
!     ---------
      implicit none
      type(SC_NoSVV_t), intent(inout) :: self


      if (associated(self % ViscousFlux)) nullify(self % ViscousFlux)

   end subroutine NoSVV_destruct
!
!///////////////////////////////////////////////////////////////////////////////
!
!     SVV filtered viscosity
!
!///////////////////////////////////////////////////////////////////////////////
!
   subroutine SVV_initialize(self, controlVariables, mesh, region)
!
!     -------
!     Modules
!     -------
      use FTValueDictionaryClass
!
!     ---------
!     Interface
!     ---------
      implicit none
      class(SC_SVV_t),         intent(inout) :: self
      type(FTValueDictionary), intent(in)    :: controlVariables
      type(HexMesh),           intent(inout) :: mesh
      integer,                 intent(in)    :: region


      ! TODO: Implement it also for region 2, but SVV does not seem very useful there...
      if (region == 2) then
         write(STD_OUT,*) "ERROR. SVV viscosity can be used only in the first region of the sensor."
         error stop
      end if
!
!     Parent initializer
!     ------------------
      call self % ArtificialViscosity_t % Initialize(controlVariables, mesh, region)
      self % region = region
!
!     Set the square root of the viscosities
!     --------------------------------------
      self % sqrt_mu1 = sqrt(self % mu1)
      self % sqrt_mu2 = sqrt(self % mu2)

      if (self % alphaIsPropToMu) then
         self % sqrt_mu2alpha = sqrt(self % mu2alpha)
      else
         self % sqrt_alpha1 = sqrt(self % alpha1)
         self % sqrt_alpha2 = sqrt(self % alpha2)
      end if
!
!     Start the SVV module
!     --------------------
      call InitializeSVV(SVV, controlVariables, mesh)

   end subroutine SVV_initialize
!
!///////////////////////////////////////////////////////////////////////////////
!
   subroutine SVV_viscosity(self, mesh, e, switch, SCflux)
!
!     ---------
!     Interface
!     ---------
      implicit none
      class(SC_SVV_t), intent(in)    :: self
      type(HexMesh),   intent(inout) :: mesh
      type(Element),   intent(inout) :: e
      real(RP),        intent(in)    :: switch
      real(RP),        intent(out)   :: SCflux(1:NCONS,       &
                                               0:e % Nxyz(1), &
                                               0:e % Nxyz(2), &
                                               0:e % Nxyz(3), &
                                               1:NDIM)
!
!     ---------------
!     Local variables
!     ---------------
      integer  :: i
      integer  :: j
      integer  :: k
      integer  :: fIDs(6)
      real(RP) :: delta
      real(RP) :: salpha
      real(RP) :: sqrt_mu(0:e % Nxyz(1), 0:e % Nxyz(2), 0:e % Nxyz(3))
      real(RP) :: sqrt_alpha(0:e % Nxyz(1), 0:e % Nxyz(2), 0:e % Nxyz(3))


      if (switch > 0.0_RP) then
!
!        Compute viscosities
!        -------------------
         select case (self % updateMethod)
         case (SC_CONST_ID)
            sqrt_mu = merge(self % sqrt_mu2,    self % sqrt_mu1,    switch >= 1.0_RP) * e % hn
            salpha  = merge(self % sqrt_alpha2, self % sqrt_alpha1, switch >= 1.0_RP) * e % hn

         case (SC_SENSOR_ID)
            sqrt_mu = (self % sqrt_mu1 * (1.0_RP-switch) + self % sqrt_mu2 * switch) * e % hn
            salpha  = (self % sqrt_alpha1 * (1.0_RP-switch) + self % sqrt_alpha2 * switch) * e % hn

#if !defined (SPALARTALMARAS)
         case (SC_SMAG_ID)

            delta = (e % geom % Volume / product(e % Nxyz + 1)) ** (1.0_RP / 3.0_RP)
            do k = 0, e % Nxyz(3) ; do j = 0, e % Nxyz(2) ; do i = 0, e % Nxyz(1)
               call Smagorinsky_ComputeViscosity(delta, e % geom % dWall(i,j,k), &
                                                 e % storage % Q(:,i,j,k),       &
                                                 e % storage % U_x(:,i,j,k),     &
                                                 e % storage % U_y(:,i,j,k),     &
                                                 e % storage % U_z(:,i,j,k),     &
                                                 sqrt_mu(i,j,k),                 &
                                                 self % smagorinsky % C,         &
                                                 self % smagorinsky % WallModel)
               sqrt_mu(i,j,k) = sqrt(sqrt_mu(i,j,k))
            end do                ; end do                ; end do
#endif

         end select

         if (self % alphaIsPropToMu) then
            sqrt_alpha = self % sqrt_mu2alpha * sqrt_mu
         else
            sqrt_alpha = salpha
         end if
!
!        Compute the viscous flux
!        ------------------------
         call SVV % ComputeInnerFluxes(e, sqrt_mu, sqrt_alpha, SCflux)

      else

         e % storage % artificialDiss = 0.0_RP
         SCflux = 0.0_RP

      end if
!
!     Project to faces
!     ----------------
      fIDs = e % faceIDs
      call e % ProlongAviscFluxToFaces(NCONS, SCflux, mesh % faces(fIDs(1)), &
                                                      mesh % faces(fIDs(2)), &
                                                      mesh % faces(fIDs(3)), &
                                                      mesh % faces(fIDs(4)), &
                                                      mesh % faces(fIDs(5)), &
                                                      mesh % faces(fIDs(6))  )

   end subroutine SVV_viscosity
!
!///////////////////////////////////////////////////////////////////////////////
!
   subroutine SVV_describe(self)
!
!     -------
!     Modules
!     -------
      use MPI_Process_Info, only: MPI_Process
!
!     ---------
!     Interface
!     ---------
      implicit none
      class(SC_SVV_t), intent(in) :: self


      if (.not. MPI_Process % isRoot) return

      write(STD_OUT,"(30X,A,A30)", advance="no") "->", "Viscosity update method: "
      select case (self % updateMethod)
         case (SC_CONST_ID);  write(STD_OUT,"(A)") SC_CONST_VAL
         case (SC_SENSOR_ID); write(STD_OUT,"(A)") SC_SENSOR_VAL
         case (SC_SMAG_ID);   write(STD_OUT,"(A)") SC_SMAG_VAL
      end select

      if (self % updateMethod == SC_SMAG_ID) then
#if !defined (SPALARTALMARAS)
         write(STD_OUT,"(30X,A,A30,1pG10.3)") "->", "LES intensity (CS): ", self % Smagorinsky % C
#endif
      else
         if (self % region == 1) then
            write(STD_OUT,"(30X,A,A30,1pG10.3)") "->","Mu viscosity 1: ", self % mu1
            write(STD_OUT,"(30X,A,A30,1pG10.3)") "->","Mu viscosity 2: ", self % mu2
         else ! self % region == 2
            write(STD_OUT,"(30X,A,A30,1pG10.3)") "->","Mu viscosity: ", self % mu2
         end if
      end if

      if (self % alphaIsPropToMu) then
         write(STD_OUT,"(30X,A,A30,1pG10.3,A)") "->", "Alpha viscosity: ", self % mu2alpha, "x mu"
      else
         if (self % region == 1) then
            write(STD_OUT,"(30X,A,A30,1pG10.3)") "->","Alpha viscosity 1: ", self % alpha1
            write(STD_OUT,"(30X,A,A30,1pG10.3)") "->","Alpha viscosity 2: ", self % alpha2
         else ! self % region == 2
            write(STD_OUT,"(30X,A,A30,1pG10.3)") "->","Alpha viscosity: ", self % alpha2
         end if
      end if

      call SVV % Describe()

   end subroutine SVV_describe

end module ShockCapturing
