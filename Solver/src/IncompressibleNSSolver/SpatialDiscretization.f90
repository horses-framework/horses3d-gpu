#include "Includes.h"
module SpatialDiscretization
      use SMConstants
      use HyperbolicDiscretizations
      use EllipticDiscretizations
      use DGIntegrals
      use MeshTypes
      use LESModels
      use HexMeshClass
      use ElementClass
      use PhysicsStorage
      use Physics
      use MPI_Face_Class
      use MPI_Process_Info
      use DGSEMClass
      use ParticlesClass
      use FluidData
      use VariableConversion, only: iNSGradientVariables, GetiNSOneFluidViscosity, GetiNSTwoFluidsViscosity
      use ProblemFileFunctions
      use BoundaryConditions, only: BCs
      use ProblemFileFunctions, only: UserDefinedSourceTermNS_f
#ifdef _HAS_MPI_
      use mpi
#endif

      private
      public   ComputeTimeDerivative, ComputeTimeDerivativeIsolated, viscousDiscretizationKey
      public   Initialize_SpaceAndTimeMethods, Finalize_SpaceAndTimeMethods,GetViscosity_selector

      abstract interface
         SUBROUTINE computeElementInterfaceFluxF(f)
            use FaceClass
            IMPLICIT NONE
            TYPE(Face)   , INTENT(inout) :: f   
         end subroutine computeElementInterfaceFluxF

         SUBROUTINE computeMPIFaceFluxF(f)
            use FaceClass
            IMPLICIT NONE
            TYPE(Face)   , INTENT(inout) :: f   
         end subroutine computeMPIFaceFluxF

         SUBROUTINE computeBoundaryFluxF(f, time)
            use SMConstants
            use FaceClass,  only: Face
            IMPLICIT NONE
            type(Face),    intent(inout) :: f
            REAL(KIND=RP)                :: time
         end subroutine computeBoundaryFluxF
      end interface
      
      character(len=LINE_LENGTH), parameter  :: viscousDiscretizationKey = "viscous discretization"
      procedure(GetViscosity_f), pointer, protected :: GetViscosity 
!
!     ========      
      CONTAINS 
!     ========      
!
!////////////////////////////////////////////////////////////////////////////////////////
!
      subroutine Initialize_SpaceAndTimeMethods(controlVariables, mesh)
         use FTValueDictionaryClass
         use Utilities, only: toLower
         use mainKeywordsModule
         use Headers
         use MPI_Process_Info
         implicit none
         class(FTValueDictionary),  intent(in)  :: controlVariables
         class(HexMesh)                         :: mesh
!
!        ---------------
!        Local variables
!        ---------------
!
         character(len=LINE_LENGTH)       :: inviscidDiscretizationName
         character(len=LINE_LENGTH)       :: viscousDiscretizationName
         
         if (.not. mesh % child) then ! If this is a child mesh, all these constructs were already initialized for the parent mesh
         
            if ( MPI_Process % isRoot ) then
               write(STD_OUT,'(/)')
               call Section_Header("Spatial discretization scheme")
               write(STD_OUT,'(/)')
            end if
   !
   !        Initialize inviscid discretization
   !        ----------------------------------
            inviscidDiscretizationName = controlVariables % stringValueForKey(inviscidDiscretizationKey,requestedLength = LINE_LENGTH)

            call toLower(inviscidDiscretizationName)
         
            select case ( trim(inviscidDiscretizationName) )
            case ( "standard" )
               if (.not. allocated(HyperbolicDiscretization)) allocate( StandardDG_t  :: HyperbolicDiscretization )

            case ( "split-form")
               if (.not. allocated(HyperbolicDiscretization)) allocate( SplitDG_t     :: HyperbolicDiscretization)

            case default
               write(STD_OUT,'(A,A,A)') 'Requested inviscid discretization "',trim(inviscidDiscretizationName),'" is not implemented.'
               write(STD_OUT,'(A)') "Implemented discretizations are:"
               write(STD_OUT,'(A)') "  * Standard"
               write(STD_OUT,'(A)') "  * Split-Form"
               errorMessage(STD_OUT)
               error stop 

            end select
               
            call HyperbolicDiscretization % Initialize(controlVariables)
   !
   !        Initialize viscous discretization
   !        ---------------------------------         
               if ( .not. controlVariables % ContainsKey(viscousDiscretizationKey) ) then
                  print*, "Input file is missing entry for keyword: viscous discretization"
                  errorMessage(STD_OUT)
                  error stop
               end if

               viscousDiscretizationName = controlVariables % stringValueForKey(viscousDiscretizationKey, requestedLength = LINE_LENGTH)
               call toLower(viscousDiscretizationName)
               
               select case ( trim(viscousDiscretizationName) )
               case("br1")
                  allocate(BassiRebay1_t     :: ViscousDiscretization)

               case("br2")
                  allocate(BassiRebay2_t     :: ViscousDiscretization)

               case("ip")
                  allocate(InteriorPenalty_t :: ViscousDiscretization)

               case default
                  write(STD_OUT,'(A,A,A)') 'Requested viscous discretization "',trim(viscousDiscretizationName),'" is not implemented.'
                  write(STD_OUT,'(A)') "Implemented discretizations are:"
                  write(STD_OUT,'(A)') "  * BR1"
                  write(STD_OUT,'(A)') "  * BR2"
                  write(STD_OUT,'(A)') "  * IP"
                  errorMessage(STD_OUT)
                  error stop 

               end select

               call ViscousDiscretization % Construct(controlVariables, ELLIPTIC_iNS)
               call ViscousDiscretization % Describe

               select case (thermodynamics % number_of_fluids)
               case(1)
                  GetViscosity => GetiNSOneFluidViscosity
               case(2)
                  GetViscosity => GetiNSTwoFluidsViscosity
               end select

!
!        Compute wall distances
!        ----------------------
         call mesh % ComputeWallDistances
!
!        Initialize models
!        -----------------
         call InitializeLESModel(LESModel, controlVariables)
         
         end if

      end subroutine Initialize_SpaceAndTimeMethods
!
!////////////////////////////////////////////////////////////////////////
!
      subroutine Finalize_SpaceAndTimeMethods
         implicit none
         IF ( ALLOCATED(HyperbolicDiscretization) ) DEALLOCATE( HyperbolicDiscretization )
         IF ( ALLOCATED(LESModel) )                 DEALLOCATE( LESModel )
      end subroutine Finalize_SpaceAndTimeMethods
!
!////////////////////////////////////////////////////////////////////////
!
      SUBROUTINE ComputeTimeDerivative( mesh, particles, time, mode, HO_Elements)
         IMPLICIT NONE 
!
!        ---------
!        Arguments
!        ---------
!
         TYPE(HexMesh), target           :: mesh
         type(Particles_t)               :: particles
         REAL(KIND=RP)                   :: time
         integer,             intent(in) :: mode
         logical, intent(in), optional   :: HO_Elements
         type(Element)                      :: e
!
!        ---------------
!        Local variables
!        ---------------
!
         INTEGER :: i, j, k, eID, nZones, zoneID
!
!        *******************************************************************
!        Construct the auxiliary state for the fluxes with density positivity
!        ******************************************************************* 
!$omp do schedule(runtime)
!$acc parallel loop gang vector_length(128) present(mesh) async(1)
         do eID = 1, size(mesh % elements)
            !$acc loop vector collapse(3)
            do k = 0, mesh % elements(eID) % Nxyz(3) ; do j = 0, mesh % elements(eID) % Nxyz(2) ; do i = 0, mesh % elements(eID) % Nxyz(1)
            mesh % elements(eID) % storage % rho(i,j,k) = mesh % elements(eID) % storage % Q(INSRHO,i,j,k)
            mesh % elements(eID) % storage % Q(INSRHO,i,j,k) = min(max(mesh % elements(eID) % storage % Q(INSRHO,i,j,k), thermodynamics % rho_min), &
                                                                   thermodynamics % rho_max)
            end do               ; end do                ; end do
         end do
!$acc end parallel loop 
!$omp end do nowait
!
!        -----------------------------------------
!        Prolongation of the solution to the faces
!        -----------------------------------------
!
!$omp parallel shared(mesh, time)
         call HexMesh_ProlongSolToFaces(mesh, NCONS)
!
!        ----------------
!        Update MPI Faces
!        ----------------
!
#ifdef _HAS_MPI_
!$omp single
         !$acc wait
         call HexMesh_UpdateMPIFacesSolution(mesh, NCONS)
!$omp end single
#endif

!        ------------------------------------------
!        Apply the Boundary conditions to the state
!        ------------------------------------------
!        This was done in the compute boundary flux before 
!        but it was called twice because we call it once in this file
!        and one in the Elliptic discretisation. So now we compute it
!        only once at the begining of time derivative and store it
! 

         nZones = size(mesh % zones)
         do zoneID=1, nZones
            CALL BCs(zoneID) % bc % FlowState(mesh, mesh % zones(zoneID))  
         enddo
!
!        -----------------
!        Compute gradients
!        -----------------
!
         if ( computeGradients ) then
            CALL DGSpatial_ComputeGradient(mesh , time)
         end if
         
         !$acc wait



#ifdef _HAS_MPI_
!$omp single
         !$acc wait
         call HexMesh_UpdateMPIFacesGradients(mesh, NGRAD)
!$omp end single
#endif
!        -----------------------
!        Compute time derivative
!        -----------------------
!
         call ComputeNSTimeDerivative(mesh = mesh , &
                                         particles = particles, &
                                         t    = time)
!
!        ***************************************
!        Return the density to its default value
!        ***************************************
!
!$omp do schedule(runtime)
!$acc parallel loop gang vector_length(128) present(mesh) async(1)
         do eID = 1, size(mesh % elements)
         !$acc loop vector collapse(3)
            do k = 0, mesh % elements(eID) % Nxyz(3) ; do j = 0, mesh % elements(eID) % Nxyz(2) ; do i = 0, mesh % elements(eID) % Nxyz(1)
                mesh % elements(eID) % storage % Q(INSRHO,i,j,k) = mesh % elements(eID) % storage % rho(i,j,k) 
            end do               ; end do                ; end do
         end do
!$acc end parallel loop 
!$omp end do

!$omp end parallel
!
      END SUBROUTINE ComputeTimeDerivative
!
!////////////////////////////////////////////////////////////////////////
!
!     This routine computes the time derivative element by element, without considering the Riemann Solvers
!     This is useful for estimating the isolated truncation error
!
      SUBROUTINE ComputeTimeDerivativeIsolated( mesh, particles, time, mode, HO_Elements)
         use EllipticDiscretizationClass
         IMPLICIT NONE 
!
!        ---------
!        Arguments
!        ---------
!
         TYPE(HexMesh), target           :: mesh
         type(Particles_t)               :: particles
         REAL(KIND=RP)                   :: time
         integer,             intent(in) :: mode
         logical, intent(in), optional   :: HO_Elements
!
!        ---------------
!        Local variables
!        ---------------
!
         INTEGER :: k
!
!        -----------------------------------------
!        Prolongation of the solution to the faces
!        -----------------------------------------
!
!$omp parallel shared(mesh, time)
         call mesh % ProlongSolutionToFaces(NCONS)
!
!        -----------------------------------------------------
!        Compute LOCAL gradients and prolong them to the faces
!        -----------------------------------------------------
!
         if ( computeGradients ) then
            CALL BaseClass_ComputeGradient( ViscousDiscretization, NCONS, NCONS, mesh , time, iNSGradientVariables)
!
!           The prolongation is usually done in the viscous methods, but not in the BaseClass
!           ---------------------------------------------------------------------------------
            !call mesh % ProlongGradientsToFaces(NCONS)
         end if

!
!        -----------------------
!        Compute time derivative
!        -----------------------
!
         call TimeDerivative_ComputeQDotIsolated(mesh = mesh , &
                                                 t    = time )
!$omp end parallel
!
      END SUBROUTINE ComputeTimeDerivativeIsolated
!
!////////////////////////////////////////////////////////////////////////////////////
!
!           Navier--Stokes procedures
!           -------------------------
!
!////////////////////////////////////////////////////////////////////////////////////
!
      subroutine ComputeNSTimeDerivative( mesh , particles, t )
         use SpongeClass, only: sponge
         implicit none
         type(HexMesh)              :: mesh
         type(Particles_t)          :: particles
         real(kind=RP)              :: t
         procedure(UserDefinedSourceTermNS_f) :: UserDefinedSourceTermNS
!
!        ---------------
!        Local variables
!        ---------------
!
         integer     :: eID , i, j, k, eq, ierr, fID, iFace, iEl
         real(kind=RP)  :: mu_smag, delta
         integer :: zoneID, zonefID
!
!        ****************
!        Volume integrals
!        ****************
!


         select type ( HyperbolicDiscretization )
         type is (StandardDG_t)
            call TimeDerivative_VolumetricContribution(mesh)         
         type is (SplitDG_t)
               call TimeDerivative_VolumetricContribution_Split(mesh)
         end select        


         if ( LESModel % active) then
            !$omp do schedule(runtime) private(i,j,k,delta,mu_smag)
                        !$acc parallel loop gang present(mesh, LESModel) async(1)
                        do eID = 1, size(mesh % elements)
                            delta = (mesh % elements(eID) % geom % Volume / product(mesh % elements(eID) % Nxyz + 1)) ** (1.0_RP / 3.0_RP)
                           !$acc loop vector collapse(3)
                            do k = 0, mesh % elements(eID) % Nxyz(3) ; do j = 0, mesh % elements(eID) % Nxyz(2) ; do i = 0, mesh % elements(eID) % Nxyz(1)
                               call LESModel_Selector(LESModel, delta, mesh % elements(eID) % geom % dWall(i,j,k), &
                                                                       mesh % elements(eID) % storage % Q(:,i,j,k),   &
                                                                       mesh % elements(eID) % storage % U_x(:,i,j,k), &
                                                                       mesh % elements(eID) % storage % U_y(:,i,j,k), &
                                                                       mesh % elements(eID) % storage % U_z(:,i,j,k), &
                                                                       mesh % elements(eID) % storage % mu_turb_NS(i,j,k) )
                                           
                               mesh % elements(eID) % storage % mu_NS(1,i,j,k) = mesh % elements(eID) % storage % mu_NS(1,i,j,k) + &
                                                                                 mesh % elements(eID) % storage % mu_turb_NS(i,j,k)
                            end do                ; end do                ; end do
                        end do
                        !$acc end parallel loop
            !$omp end do
         end if

!
!        Compute viscosity at interior and boundary faces
!        ------------------------------------------------
         call compute_viscosity_at_faces(size(mesh % faces_interior), 2, mesh % faces_interior, mesh)
         call compute_viscosity_at_faces(size(mesh % faces_boundary), 1, mesh % faces_boundary, mesh)
!
!        ******************************************
!        Compute Riemann solver of non-shared faces
!        ******************************************
!
!$omp do schedule(runtime) private(fID)
         !$acc parallel loop gang num_gangs(size(mesh % faces_interior)) present(mesh, mesh % faces, mesh % faces_interior) async(1)
         do iFace = 1, size(mesh % faces_interior)
            fID = mesh % faces_interior(iFace)
            call computeElementInterfaceFlux_INS(mesh % faces(fID))
         end do
         !$acc end parallel loop
!$omp end do nowait

               CALL computeBoundaryFlux_iNS(mesh, t) 
!
!        **************************************************************
!        Surface integrals and scaling of elements without shared faces
!        **************************************************************
! 
!$omp do schedule(runtime) private(i,j,k,eID)
         !$acc parallel loop gang present(mesh, mesh % elements, mesh % elements_sequential) copyin(t)  async(1)
         do iEl = 1, size(mesh % elements_sequential)
            eID = mesh % elements_sequential(iEl)
            call TimeDerivative_FacesContribution(mesh % elements(eID), t, mesh)
         end do
         !$acc end parallel loop 
!$omp end do
!
!        ****************************
!        Wait until messages are sent
!        ****************************
!
#ifdef _HAS_MPI_
         if ( MPI_Process % doMPIAction ) then
!$omp single
            !$acc wait
            call HexMesh_GatherMPIFacesGradients(mesh, NGRAD)
!$omp end single
!
!           Compute viscosity at MPI faces
!           ------------------------------
            call compute_viscosity_at_faces(size(mesh % faces_mpi), 2, mesh % faces_mpi, mesh)
!
!           **************************************
!           Compute Riemann solver of shared faces
!           **************************************
!
!$omp do schedule(runtime) private(fID)
!$acc parallel loop gang present(mesh) async(1)
         do iFace = 1, size(mesh % faces_mpi)
            fID = mesh % faces_mpi(iFace)
            call ComputeMPIFaceFlux_iNS(mesh % faces(fID))
         end do
!$acc end parallel loop
!$omp end do nowait
!
!           ***********************************************************
!           Surface integrals and scaling of elements with shared faces
!           ***********************************************************
! 
!$omp do schedule(runtime) private(i,j,k,eID)
!$acc parallel loop gang num_gangs(size(mesh % elements_mpi)) vector_length(128) present(mesh) async(1)
         do iEl = 1, size(mesh % elements_mpi)
            eID = mesh % elements_mpi(iEl)
            call TimeDerivative_FacesContribution(mesh % elements(eID), t, mesh)
         end do
!$acc end parallel loop 
!$omp end do
!
!           Add a MPI Barrier
!           -----------------
!$omp single
            call mpi_barrier(MPI_COMM_WORLD, ierr)
!$omp end single
         end if
#endif
!
!        ***********
!        Add gravity
!        ***********
!
!$omp do schedule(runtime) private(i,j,k)
            do eID = 1, size(mesh % elements)
               associate(e => mesh % elements(eID))
               do k = 0, e % Nxyz(3) ; do j = 0, e % Nxyz(2) ; do i = 0, e % Nxyz(1)
                  e % storage % QDot(INSRHOU:INSRHOW,i,j,k) = e % storage % QDot(INSRHOU:INSRHOW,i,j,k) + &
                                                        e % storage % Q(INSRHO,i,j,k) * &
                                    dimensionless % invFr2 * dimensionless % gravity_dir

               end do                ; end do                ; end do
               end associate
            end do
!$omp end do
!
!           ***************
!           Add source term
!           ***************
!$omp do schedule(runtime) private(i,j,k)
            do eID = 1, mesh % no_of_elements
               associate ( e => mesh % elements(eID) )
               do k = 0, e % Nxyz(3)   ; do j = 0, e % Nxyz(2) ; do i = 0, e % Nxyz(1)
                  call UserDefinedSourceTermNS(e % geom % x(:,i,j,k), e % storage % Q(:,i,j,k), t, e % storage % S_NS(:,i,j,k), thermodynamics, dimensionless, refValues)
               end do                  ; end do                ; end do
               end associate
            end do
!$omp end do

            ! for the sponge, loops are in the internal subroutine as values are precalculated
            !call sponge % addSource(mesh)

!$omp do schedule(runtime) private(i,j,k)
!$acc parallel loop gang vector_length(128) present(mesh) async(1)
         do eID = 1, mesh % no_of_elements
            !$acc loop vector collapse(4)
            do k = 0, mesh % elements(eID) % Nxyz(3)   ; do j = 0, mesh % elements(eID) % Nxyz(2) ; do i = 0, mesh % elements(eID) % Nxyz(1) ; do eq = 1, NCONS
               mesh % elements(eID) % storage % QDot(eq,i,j,k) = mesh % elements(eID) % storage % QDot(eq,i,j,k) + mesh % elements(eID) % storage % S_NS(eq,i,j,k)
            end do   ;  end do   ;  end do   ;  end do
         end do
!$acc end parallel loop
!$omp end do
!
!           ********************
!           Add Particles source
!           ********************
            if (.not. mesh % child) then
               if ( particles % active ) then             
!$omp do schedule(runtime)
                  do eID = 1, size(mesh % elements)
                  !   call particles % AddSource(mesh % elements(eID), t, thermodynamics, dimensionless, refValues)
                  end do
!$omp end do
               endif 
            end if



      end subroutine ComputeNSTimeDerivative
!
!////////////////////////////////////////////////////////////////////////
!
!     -------------------------------------------------------------------------------
!     This routine computes Qdot neglecting the interaction with neighboring elements
!     and boundaries. Therefore, the external states are not needed.
!     -------------------------------------------------------------------------------
      subroutine TimeDerivative_ComputeQDotIsolated( mesh , t )
         implicit none
         type(HexMesh)              :: mesh
         real(kind=RP)              :: t
!
!        ---------------
!        Local variables
!        ---------------
!
         integer     :: eID , i, j, k, fID
!
!        ****************
!        Volume integrals
!        ****************
!
!$omp do schedule(runtime) 
         do eID = 1 , size(mesh % elements)
            call TimeDerivative_StrongVolumetricContribution( mesh % elements(eID) , t)
         end do
!$omp end do
!
!        *******************
!        Scaling of elements
!        *******************
! 
!$omp do schedule(runtime) private(i,j,k)
         do eID = 1, size(mesh % elements) 
            associate(e => mesh % elements(eID))

            do k = 0, e % Nxyz(3) ; do j = 0, e % Nxyz(2) ; do i = 0, e % Nxyz(1) 
               e % storage % QDot(:,i,j,k) = e % storage % QDot(:,i,j,k) / e % geom % jacobian(i,j,k) 
            end do         ; end do          ; end do 
            end associate 
         end do
!$omp end do
         
      end subroutine TimeDerivative_ComputeQDotIsolated
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
      subroutine TimeDerivative_StrongVolumetricContribution( e , t )
         use HexMeshClass
         use ElementClass
         implicit none
         type(Element)      :: e
         real(kind=RP)      :: t

!
!        ---------------
!        Local variables
!        ---------------
!
         real(kind=RP) :: inviscidContravariantFlux ( 1:NCONS, 0:e%Nxyz(1) , 0:e%Nxyz(2) , 0:e%Nxyz(3), 1:NDIM ) 
         real(kind=RP) :: fSharp(1:NCONS, 0:e%Nxyz(1), 0:e%Nxyz(1), 0:e%Nxyz(2), 0:e%Nxyz(3))
         real(kind=RP) :: gSharp(1:NCONS, 0:e%Nxyz(2), 0:e%Nxyz(1), 0:e%Nxyz(2), 0:e%Nxyz(3))
         real(kind=RP) :: hSharp(1:NCONS, 0:e%Nxyz(3), 0:e%Nxyz(1), 0:e%Nxyz(2), 0:e%Nxyz(3))
         real(kind=RP) :: viscousContravariantFlux  ( 1:NCONS, 0:e%Nxyz(1) , 0:e%Nxyz(2) , 0:e%Nxyz(3), 1:NDIM ) 
         real(kind=RP) :: contravariantFlux         ( 1:NCONS, 0:e%Nxyz(1) , 0:e%Nxyz(2) , 0:e%Nxyz(3), 1:NDIM ) 
         integer       :: eID
!
!        *************************************
!        Compute interior contravariant fluxes
!        *************************************
!
!        Compute inviscid contravariant flux
!        -----------------------------------
         call HyperbolicDiscretization % ComputeInnerFluxes ( e , iEulerFlux, inviscidContravariantFlux ) 
!
!        Compute viscous contravariant flux
!        ----------------------------------
         call ViscousDiscretization  % ComputeInnerFluxes ( NCONS, NCONS, iViscousFlux, GetViscosity, e , viscousContravariantFlux) 
!
!        ************************
!        Perform volume integrals
!        ************************
!
         select type ( HyperbolicDiscretization )
         type is (StandardDG_t)
!
!           Compute the total Navier-Stokes flux
!           ------------------------------------
            contravariantFlux = inviscidContravariantFlux - viscousContravariantFlux 
!
!           Perform the Weak Volume Green integral
!           --------------------------------------
            e % storage % QDot = ScalarStrongIntegrals % StdVolumeGreen ( e , NCONS, contravariantFlux ) 

         type is (SplitDG_t)
            error stop ':: TimeDerivative_StrongVolumetricContribution not implemented for split form'
!~ !
!~ !           Compute sharp fluxes for skew-symmetric approximations
!~ !           ------------------------------------------------------
!~             call HyperbolicDiscretization % ComputeSplitFormFluxes(e, inviscidContravariantFlux, fSharp, gSharp, hSharp)
!~ !
!~ !           Perform the Weak volume green integral
!~ !           --------------------------------------
!~             viscousContravariantFlux = viscousContravariantFlux + SVVContravariantFlux

!~             e % storage % QDot = -ScalarWeakIntegrals % SplitVolumeDivergence( e, fSharp, gSharp, hSharp, viscousContravariantFlux)

         end select

      end subroutine TimeDerivative_StrongVolumetricContribution
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
      subroutine TimeDerivative_VolumetricContribution(mesh)
         use HexMeshClass
         use ElementClass
         use DGIntegrals
         implicit none
         type(HexMesh), intent (inout)           :: mesh

!
!        ---------------
!        Local variables
!        ---------------
!        
         integer       :: i, j, k, eq, eID
         real(kind=RP) :: inviscidFlux(1:NCONS, 1:NDIM)
         real(kind=RP) :: viscousFlux(1:NCONS, 1:NDIM)
         real(kind=RP) :: mu, beta, keppa
!
!        *************************************
!        Compute interior contravariant fluxes
!        *************************************
!
!        Compute inviscid contravariant flux
!        -----------------------------------
         !$omp do schedule(runtime)
         !$acc parallel loop gang vector_length(128) num_gangs(9700) present(mesh) async(1)
         do eID = 1 , size(mesh % elements)
         !$acc loop vector collapse(3) private(inviscidFlux, viscousFlux)
         do k = 0, mesh % elements(eID) % Nxyz(3) ; do j = 0, mesh % elements(eID) % Nxyz(2) ; do i = 0, mesh % elements(eID) % Nxyz(1)
                  
            call iEulerFlux(mesh % elements(eID) % storage % Q(:,i,j,k), inviscidFlux, mesh % elements(eID) % storage % rho(i,j,k))

            call GetViscosity_selector(mesh % elements(eID) % storage % Q(INSRHO,i,j,k), mu)

            call iViscousFlux( NCONS, NGRAD, mesh % elements(eID) % storage % Q(:,i,j,k) , mesh % elements(eID) % storage % U_x(:,i,j,k) , & 
                                    mesh % elements(eID) % storage % U_y(:,i,j,k) , mesh % elements(eID) % storage % U_z(:,i,j,k), mu, beta, keppa, viscousFlux)
            
         
            do eq =1, NCONS

               inviscidFlux(eq,:) = inviscidFlux(eq,:) - viscousFlux(eq,:)

               mesh % elements(eID) % storage % contravariantFlux(eq,i,j,k, IX) = inviscidFlux(eq,IX) * mesh % elements(eID) % geom % jGradXi(IX,i,j,k)  &
                                                  + inviscidFlux(eq,IY) * mesh % elements(eID) % geom % jGradXi(IY,i,j,k)  &
                                                  + inviscidFlux(eq,IZ) * mesh % elements(eID) % geom % jGradXi(IZ,i,j,k)

               mesh % elements(eID) % storage % contravariantFlux(eq,i,j,k, IY) = inviscidFlux(eq,IX) * mesh % elements(eID) % geom % jGradEta(IX,i,j,k)  &
                                                  + inviscidFlux(eq,IY) * mesh % elements(eID) % geom % jGradEta(IY,i,j,k)  &
                                                  + inviscidFlux(eq,IZ) * mesh % elements(eID) % geom % jGradEta(IZ,i,j,k)

               mesh % elements(eID) % storage % contravariantFlux(eq,i,j,k, IZ) = inviscidFlux(eq,IX) * mesh % elements(eID) % geom % jGradZeta(IX,i,j,k)  &
                                                  + inviscidFlux(eq,IY) * mesh % elements(eID) % geom % jGradZeta(IY,i,j,k)  &
                                                  + inviscidFlux(eq,IZ) * mesh % elements(eID) % geom % jGradZeta(IZ,i,j,k)
               !initialize to 0 to accumulate
               mesh % elements(eID) % storage % Qdot(eq,i,j,k)  = 0.0_RP
               end do
         end do               ; end do                ; end do

         call ScalarWeakIntegrals_StdVolumeGreen( mesh % elements(eID) % Nxyz, NCONS,&
                                                  mesh % elements(eID) % storage % contravariantFlux, &
                                                  mesh % elements(eID) % storage % QDot)

         end do
         !$acc end parallel loop 
         !$omp end do

      end subroutine TimeDerivative_VolumetricContribution
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
      subroutine TimeDerivative_VolumetricContribution_Split(mesh)
         use HexMeshClass
         use ElementClass
         use DGIntegrals
         use NodalStorageClass, only: NodalStorage
         use RiemannSolvers_iNS
         implicit none
         type(HexMesh), intent (inout)           :: mesh

!
!        ---------------
!        Local variables
!        ---------------
!        
         integer       :: i, j, k, l, eq, eID
         real(kind=RP) :: Flux(1:NCONS, 1:NDIM)
         real(kind=RP) :: mu, beta, keppa
!
!        *************************************
!        Compute interior contravariant fluxes
!        *************************************
!
!        Compute inviscid contravariant flux
!        -----------------------------------
         !$omp do schedule(runtime)
         !$acc parallel loop gang vector_length(128) num_gangs(9700) present(mesh) async(1)
         do eID = 1 , size(mesh % elements)
         !$acc loop vector collapse(3) private(Flux)
         do k = 0, mesh % elements(eID) % Nxyz(3) ; do j = 0, mesh % elements(eID) % Nxyz(2) ; do i = 0, mesh % elements(eID) % Nxyz(1)
                  
            call GetViscosity_selector(mesh % elements(eID) % storage % Q(INSRHO,i,j,k), mu)

            call iViscousFlux( NCONS, NGRAD, mesh % elements(eID) % storage % Q(:,i,j,k) , mesh % elements(eID) % storage % U_x(:,i,j,k) , & 
                                    mesh % elements(eID) % storage % U_y(:,i,j,k) , mesh % elements(eID) % storage % U_z(:,i,j,k), mu, beta, keppa, Flux)
            
         
            do eq =1, NCONS


               mesh % elements(eID) % storage % contravariantFlux(eq,i,j,k, IX) = -Flux(eq,IX) * mesh % elements(eID) % geom % jGradXi(IX,i,j,k)  &
                                                  - Flux(eq,IY) * mesh % elements(eID) % geom % jGradXi(IY,i,j,k)  &
                                                  - Flux(eq,IZ) * mesh % elements(eID) % geom % jGradXi(IZ,i,j,k)

               mesh % elements(eID) % storage % contravariantFlux(eq,i,j,k, IY) = -Flux(eq,IX) * mesh % elements(eID) % geom % jGradEta(IX,i,j,k)  &
                                                  - Flux(eq,IY) * mesh % elements(eID) % geom % jGradEta(IY,i,j,k)  &
                                                  - Flux(eq,IZ) * mesh % elements(eID) % geom % jGradEta(IZ,i,j,k)

               mesh % elements(eID) % storage % contravariantFlux(eq,i,j,k, IZ) = - Flux(eq,IX) * mesh % elements(eID) % geom % jGradZeta(IX,i,j,k)  &
                                                  - Flux(eq,IY) * mesh % elements(eID) % geom % jGradZeta(IY,i,j,k)  &
                                                  - Flux(eq,IZ) * mesh % elements(eID) % geom % jGradZeta(IZ,i,j,k)
               !initialize to 0 to accumulate
               mesh % elements(eID) % storage % Qdot(eq,i,j,k)  = 0.0_RP
               end do
         end do               ; end do                ; end do

         call ScalarWeakIntegrals_StdVolumeGreen( mesh % elements(eID) % Nxyz, NCONS,&
                                                  mesh % elements(eID) % storage % contravariantFlux, &
                                                  mesh % elements(eID) % storage % QDot)

         !$acc loop vector collapse(3) private(Flux)
         do k = 0, mesh % elements(eID) % Nxyz(3) ; do j = 0, mesh % elements(eID) % Nxyz(2) ; do i = 0, mesh % elements(eID) % Nxyz(1)
            !$acc loop seq
            do l = 0, mesh % elements(eID) % Nxyz(1)
                  call TwoPointFlux_Selector(mesh % elements(eID) % storage % Q(:,i,j,k), mesh % elements(eID) % storage % Q(:,l,j,k), mesh % elements(eID) % geom % jGradXi(:,i,j,k),  mesh % elements(eID) % geom % jGradXi(:,l,j,k), Flux(:,IX))
                  call TwoPointFlux_Selector(mesh % elements(eID) % storage % Q(:,i,j,k), mesh % elements(eID) % storage % Q(:,i,l,k), mesh % elements(eID) % geom % jGradEta(:,i,j,k), mesh % elements(eID) % geom % jGradEta(:,i,l,k), Flux(:,IY))
                  call TwoPointFlux_Selector(mesh % elements(eID) % storage % Q(:,i,j,k), mesh % elements(eID) % storage % Q(:,i,j,l), mesh % elements(eID) % geom % jGradZeta(:,i,j,k), mesh % elements(eID) % geom % jGradZeta(:,i,j,l), Flux(:,IZ))
                  !$acc loop seq
                   do eq = 1, NCONS
                      mesh % elements(eID) % storage % QDot(eq,i,j,k) = mesh % elements(eID) % storage % QDot(eq,i,j,k) &
                                                                      - NodalStorage(mesh % elements(eID) % Nxyz(1)) % sharpD(i,l) *  Flux(eq,IX) &
                                                                      - NodalStorage(mesh % elements(eID) % Nxyz(2)) % sharpD(j,l) *  Flux(eq,IY) &
                                                                      - NodalStorage(mesh % elements(eID) % Nxyz(3)) % sharpD(k,l) *  Flux(eq,IZ)
                   end do
            end do
         end do               ; end do                ; end do

         end do
         !$acc end parallel loop 
         !$omp end do

      end subroutine TimeDerivative_VolumetricContribution_Split

!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
      subroutine compute_viscosity_at_faces(no_of_faces, no_of_sides, face_ids, mesh)
         implicit none
         integer, intent(in)           :: no_of_faces
         integer, intent(in)           :: no_of_sides
         integer, intent(in)           :: face_ids(no_of_faces)
         class(HexMesh), intent(inout) :: mesh
!
!        ---------------
!        Local variables
!        ---------------
!
         integer       :: iFace, i, j, side, fID
         real(kind=RP) :: delta, mu_smag


         if ( LESModel % Active ) then
!$omp do schedule(runtime) private(i,j,delta,mu_smag)
            !$acc parallel loop gang present(mesh, LESModel) async(1)
            do iFace = 1, no_of_faces
               delta = sqrt(mesh % faces(face_ids(iFace)) % geom % surface / product(mesh % faces(face_ids(iFace)) % Nf + 1))
               !$acc loop vector collapse(3)
               do j = 0, mesh % faces(face_ids(iFace)) % Nf(2) ; do i = 0, mesh % faces(face_ids(iFace)) % Nf(1)
                  do side = 1, no_of_sides
                     call LESModel_Selector(LESModel, delta, mesh % faces(face_ids(iFace)) % geom % dWall(i,j), &
                                                             mesh % faces(face_ids(iFace)) % storage(side) % Q(:,i,j),   &
                                                             mesh % faces(face_ids(iFace)) % storage(side) % U_x(:,i,j), &
                                                             mesh % faces(face_ids(iFace)) % storage(side) % U_y(:,i,j), &
                                                             mesh % faces(face_ids(iFace)) % storage(side) % U_z(:,i,j), &
                                                                                    mu_smag)

                     mesh % faces(face_ids(iFace)) % storage(side) % mu_NS(1,i,j) = mesh % faces(face_ids(iFace)) % storage(side) % mu_NS(1,i,j) + mu_smag
                     !mesh % faces(face_ids(iFace)) % storage(side) % mu_NS(2,i,j) = mesh % faces(face_ids(iFace)) % storage(side) % mu_NS(2,i,j) + mu_smag * dimensionless % mu_to_kappa
                  end do
               end do              ; end do
            end do
            !$acc end parallel loop
!$omp end do
         end if




      end subroutine compute_viscosity_at_faces
!
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
      subroutine TimeDerivative_FacesContribution( e , t , mesh)
         !$acc routine vector
         use HexMeshClass
         implicit none
         type(Element)           :: e
         real(kind=RP)           :: t
         type(HexMesh)           :: mesh

         integer                 :: i,j,k,eID,eq

         call ScalarWeakIntegrals_StdFace( NCONS, e % Nxyz, &
                      mesh % faces(e % faceIDs(EFRONT))  % storage(e % faceSide(EFRONT))  % fStar, &
                      mesh % faces(e % faceIDs(EBACK))   % storage(e % faceSide(EBACK))   % fStar, &
                      mesh % faces(e % faceIDs(EBOTTOM)) % storage(e % faceSide(EBOTTOM)) % fStar, &
                      mesh % faces(e % faceIDs(ERIGHT))  % storage(e % faceSide(ERIGHT))  % fStar, &
                      mesh % faces(e % faceIDs(ETOP))    % storage(e % faceSide(ETOP))    % fStar, &
                      mesh % faces(e % faceIDs(ELEFT))   % storage(e % faceSide(ELEFT))   % fStar, &
                      -1, e % storage % QDot )

         !$acc loop vector collapse(3)
         do k = 0, e % Nxyz(3) ; do j = 0, e % Nxyz(2) ; do i = 0, e % Nxyz(1)
            !$acc loop seq
            do eq = 1,NCONS
            e % storage % QDot(eq,i,j,k) = e % storage % QDot(eq,i,j,k)  / e % geom % jacobian(i,j,k)
            enddo
         end do         ; end do          ; end do
         
      end subroutine TimeDerivative_FacesContribution
!
!///////////////////////////////////////////////////////////////////////////////////////////// 
! 
!        Riemann solver drivers 
!        ---------------------- 
! 
!///////////////////////////////////////////////////////////////////////////////////////////// 
! 
      SUBROUTINE computeElementInterfaceFlux_iNS(f)
         !$acc routine vector
         use FaceClass
         use RiemannSolvers_iNS
         use EllipticBR1
         IMPLICIT NONE
         TYPE(Face)   , INTENT(inout) :: f   
         integer       :: i, j, eq
         real(kind=RP) :: muL, muR

         !$acc loop vector collapse(2) 
               do j = 0, f % Nf(2) ; do i = 0, f % Nf(1)

                  call GetViscosity_selector(f % storage(1) % Q(INSRHO,i,j), muL)
                  call GetViscosity_selector(f % storage(2) % Q(INSRHO,i,j), muR)
                  !mu = 0.5_RP * (muL + muR)
!        
!                 --------------
!                 Viscous fluxes
!                 --------------
!        
                  !check  unstar unstar
                  call iViscousFlux( NCONS, NGRAD, f % storage(1) % Q(:,i,j) , f % storage(1) % U_x(:,i,j), & 
                                    f % storage(1) % U_y(:,i,j) , f % storage(1) % U_z(:,i,j), muL, 0.0_RP, 0.0_RP, f % storage(1) % unStar(:,:,i,j))

                  call iViscousFlux( NCONS, NGRAD, f % storage(2) % Q(:,i,j) , f % storage(2) % U_x(:,i,j), & 
                                    f % storage(2) % U_y(:,i,j) , f % storage(2) % U_z(:,i,j), muR, 0.0_RP, 0.0_RP, f % storage(2) % unStar(:,:,i,j))

            end do ; end do

            call BR1_RiemannSolver_acc(f, NCONS, NGRAD, f % storage(2) % FStar)

            call RiemannSolver_Selector(f % Nf(1), &
                                        f % Nf(2), &
                                        f % storage(1) % Q, &
                                        f % storage(2) % Q, &
                                        f % geom % normal, &
                                        f % geom % t1, &
                                        f % geom % t2, &
                                        f % storage(1) % FStar )

!        ------------------------
!        Multiply by the Jacobian
!        ------------------------
         !$acc loop vector collapse(2)
            do j = 0, f % Nf(2) ; do i = 0, f % Nf(1)
               !$acc loop seq
               do eq = 1, NCONS
                  f % storage(1) % FStar(eq,i,j) = (f % storage(1) % FStar(eq,i,j) - f % storage(2) % FStar(eq,i,j)) * f % geom % jacobian(i,j)
               enddo
            end do ;  end do
   !
   !        ---------------------------
   !        Return the flux to elements
   !        ---------------------------
   !       
           call Face_ProjectFluxToElements(f, NCONS, f % storage(1) % FStar, 1)
           call Face_ProjectFluxToElements(f, NCONS, f % storage(1) % FStar, 2)

   

      END SUBROUTINE computeElementInterfaceFlux_iNS

      SUBROUTINE computeMPIFaceFlux_iNS(f)
         !$acc routine vector
         use FaceClass
         use RiemannSolvers_iNS
         use EllipticBR1
         IMPLICIT NONE
         TYPE(Face)   , INTENT(inout) :: f   
         integer       :: i, j, eq
         integer       :: Sidearray, maxId
         real(kind=RP) :: mu
!
!        --------------
!        Invscid fluxes
!        --------------
!
         !$acc loop vector collapse(2) 
         DO j = 0, f % Nf(2)
            DO i = 0, f % Nf(1)
!      
!              --------------
!              Viscous fluxes
!              --------------
!      
               call GetViscosity_selector(f % storage(1) % Q(INSRHO,i,j), mu)

                  call iViscousFlux( NCONS, NGRAD, f % storage(1) % Q(:,i,j) , f % storage(1) % U_x(:,i,j), & 
                                    f % storage(1) % U_y(:,i,j) , f % storage(1) % U_z(:,i,j), mu, 0.0_RP, 0.0_RP, f % storage(1) % unStar(:,:,i,j))

                  call iViscousFlux( NCONS, NGRAD, f % storage(2) % Q(:,i,j) , f % storage(2) % U_x(:,i,j), & 
                                    f % storage(2) % U_y(:,i,j) , f % storage(2) % U_z(:,i,j), mu, 0.0_RP, 0.0_RP, f % storage(2) % unStar(:,:,i,j))

            end do ; end do

!
            call BR1_RiemannSolver_acc(f, NCONS, NGRAD, f % storage(2) % FStar)

            call RiemannSolver_Selector(f % Nf(1), &
                                        f % Nf(2), &
                                        f % storage(1) % Q, &
                                        f % storage(2) % Q, &
                                        f % geom % normal, &
                                        f % geom % t1, &
                                        f % geom % t2, &
                                        f % storage(1) % FStar )
!
!              Multiply by the Jacobian
!              ------------------------
         !$acc loop vector collapse(2)
            do j = 0, f % Nf(2) ; do i = 0, f % Nf(1)
               !$acc loop seq
               do eq = 1, NCONS
                  f % storage(1) % FStar(eq,i,j) = (f % storage(1) % FStar(eq,i,j) - f % storage(2) % FStar(eq,i,j)) * f % geom % jacobian(i,j)
               enddo
            end do ;  end do
!
!        ---------------------------
!        Return the flux to elements: The sign in eR % storage % FstarB has already been accouted.
!        ---------------------------
!
         maxId=MAXVAL(f % elementIDs)
         do i=1,SIZE(f % elementIDs)
             if(f % elementIDs(i)==maxId)THEN
               Sidearray=i
                 exit
             endif
         end do
         call Face_ProjectFluxToElements(f, NCONS, f % storage(1) % FStar, Sidearray)

      end subroutine ComputeMPIFaceFlux_iNS

      SUBROUTINE computeBoundaryFlux_iNS(mesh, time)
      USE ElementClass
      use FaceClass
      USE RiemannSolvers_iNS
      IMPLICIT NONE
!
!     ---------
!     Arguments
!     ---------
!
      type(HexMesh), intent(inout)    :: mesh
      REAL(KIND=RP)                :: time
!
!     ---------------
!     Local variables
!     ---------------
!
      INTEGER                         :: i, j, eq
      INTEGER                         :: nZones, zoneID, zonefID, fID
      real(kind=RP)                   :: mu
!
!     -------------------
!     Get external states
!     -------------------
!     

       nZones = size(mesh % zones)
       do zoneID=1, nZones


          !$acc parallel loop gang present(mesh) async(zoneID)
          do zonefID = 1, mesh % zones(zoneID) % no_of_faces
             fID =  mesh % zones(zoneID) % faces(zonefID)
     
             !$acc loop vector collapse(2)
             do j = 0, mesh % faces(fID) % Nf(2) ;  do i = 0, mesh % faces(fID) % Nf(1)

               call GetViscosity_selector(mesh % faces(fID) % storage(1) % Q(INSRHO,i,j), mu)
                
               call iViscousFlux( NCONS, NGRAD, mesh % faces(fID) % storage(1) % Q(:,i,j) , &
                                                mesh % faces(fID) % storage(1) % U_x(:,i,j), & 
                                                mesh % faces(fID) % storage(1) % U_y(:,i,j) , &
                                                mesh % faces(fID) % storage(1) % U_z(:,i,j), &
                                                mu, 0.0_RP, 0.0_RP, &
                                                mesh % faces(fID) % storage(1) % unStar(:,:,i,j))
             enddo ; enddo
            
             !$acc loop vector collapse(2)
             do j = 0, mesh % faces(fID) % Nf(2) ;  do i = 0, mesh % faces(fID) % Nf(1)
                !$acc loop seq
                do eq = 1, NCONS

                  mesh % faces(fID) % storage(2) % FStar(eq,i,j) = mesh % faces(fID) % storage(1) % unStar(eq,IX,i,j)* mesh % faces(fID) % geom % normal(IX,i,j) &
                                                                 + mesh % faces(fID) % storage(1) % unStar(eq,IY,i,j)* mesh % faces(fID) % geom % normal(IY,i,j) &
                                                                 + mesh % faces(fID) % storage(1) % unStar(eq,IZ,i,j)* mesh % faces(fID) % geom % normal(IZ,i,j)
                enddo
             enddo ; enddo
          enddo
          !$acc end parallel loop 

          !$acc wait

          CALL BCs(zoneID) % bc % FlowNeumann(mesh, mesh % zones(zoneID))                             
         
          !$acc wait

          !$acc parallel loop gang present(mesh) async(zoneID)
          do zonefID = 1, mesh % zones(zoneID) % no_of_faces
             fID =  mesh % zones(zoneID) % faces(zonefID)

            call RiemannSolver_Selector(Nx = mesh % faces(fID) % Nf(1), &
                                        Ny = mesh % faces(fID) % Nf(2), &
                                        QLeft  = mesh % faces(fID) % storage(1) % Q, &
                                        QRight = mesh % faces(fID) % storage(2) % Q, &
                                        nHat   = mesh % faces(fID) % geom % normal, &
                                        t1     = mesh % faces(fID) % geom % t1, &
                                        t2     = mesh % faces(fID) % geom % t2, &
                                        flux   = mesh % faces(fID) % storage(1) % FStar )
!           ------------------------
!           Multiply by the Jacobian
!           ------------------------
            !$acc loop vector collapse(2)
              do j = 0, mesh % faces(fID) % Nf(2) ; do i = 0, mesh % faces(fID) % Nf(1)
                 !$acc loop seq
                 do eq = 1, NCONS
                   mesh % faces(fID) % storage(1) % FStar(eq,i,j) = (mesh % faces(fID) % storage(1) % FStar(eq,i,j)  - &
                                                                     mesh % faces(fID) % storage(2) % FStar(eq,i,j)) * &
                                                                     mesh % faces(fID) % geom % jacobian(i,j)
                 enddo
             end do ;  end do
            
!
!           ---------------------------
!           Return the flux to elements
!           ---------------------------
!
             call Face_ProjectFluxToElements(mesh % faces(fID), NCONS, mesh % faces(fID) % storage(1) % FStar, 1)

          enddo
          !$acc end parallel loop 
       enddo

      END SUBROUTINE computeBoundaryFlux_iNS
!
!////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!              GRADIENT PROCEDURES
!              -------------------
!
!////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
      subroutine DGSpatial_ComputeGradient( mesh , time)
         use HexMeshClass
         use PhysicsStorage, only: NCONS
         implicit none
         type(HexMesh)                  :: mesh
         real(kind=RP),      intent(in) :: time

         call ViscousDiscretization % ComputeGradient( NCONS, NCONS, mesh , time , iNSGradientVariables)

      end subroutine DGSpatial_ComputeGradient

      subroutine DensityLimiter(N,Q)
         implicit none
         integer,       intent(in)    :: N(3)
         real(kind=RP), intent(inout) :: Q(1:NCONS,0:N(1),0:N(2),0:N(3))
!
!        ---------------
!        Local variables
!        ---------------
!
         integer  :: i, j, k 
         real(kind=RP) :: rhoIn01, p, rhomin, rhomax

         rhomin = thermodynamics % rho_min / refValues % rho
         rhomax = thermodynamics % rho_max / refValues % rho

         do k = 0, N(3) ; do j = 0, N(2) ; do i = 0, N(1)
            rhoIn01 = (Q(INSRHO,i,j,k)-rhomin)/(rhomax-rhomin)

            if ( rhoIn01 .ge. 1.0_RP ) then
               Q(INSRHO,i,j,k) = rhomax

            elseif ( rhoIn01 .le. 0.0_RP ) then
               Q(INSRHO,i,j,k) = rhomin

            else
               !p = POW3(rhoIn01)*(6.0_RP*POW2(rhoIn01)-15.0_RP*rhoIn01+10.0_RP)
               p = rhoIn01

               Q(INSRHO,i,j,k) = (rhomax-rhomin)*p + rhomin
         
            end if

         end do         ; end do         ; end do

      end subroutine DensityLimiter
!
!////////////////////////////////////////////////////////////////////////////////////////
!

      subroutine GetViscosity_selector(phi, mu)
         !$acc routine seq
         use VariableConversion
         implicit none
         real(kind=RP), intent(in)   :: phi
         real(kind=RP), intent(out)  :: mu
!
!        ---------------
!        Local variables
!        ---------------
!
         select case (thermodynamics % number_of_fluids)
         case(1)
            call GetiNSOneFluidViscosity(phi, mu)
         case(2)
            call GetiNSTwoFluidsViscosity(phi, mu)
         end select

      end subroutine GetViscosity_selector
!
!////////////////////////////////////////////////////////////////////////////////////////
!
end module SpatialDiscretization