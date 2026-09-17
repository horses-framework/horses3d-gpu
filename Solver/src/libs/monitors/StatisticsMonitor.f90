#include "Includes.h"
#if defined(NAVIERSTOKES) || defined(INCNS)
module StatisticsMonitor
   use SMConstants
   use HexMeshClass
   use StorageClass
   use StopwatchClass
   use Utilities, only: GreatestCommonDivisor, toLower
#ifdef _HAS_MPI_
   use mpi
#endif

   private
   public     StatisticsMonitor_t, U, V, W, UU, VV, WW, UV, UW, VW
   public     NO_OF_VARIABLES_Sij, NO_OF_VARIABLES, NO_OF_FAVRE_VARS
   public     stats_doReynolds, stats_doFavre
!
!  Commands for the parameter file
!  -------------------------------
#define START_COMMAND "@start"
#define PAUSE_COMMAND "@pause"
#define STOP_COMMAND "@stop"
#define RESET_COMMAND "@reset"
#define DUMP_COMMAND "@dump"
!
!  Statistics monitor states
!  -------------------------
#define OFF 0
#define WAITING_STATE 1
#define STOP_STATE 2
#define ACTIVE_STATE 3
!
!  Statistics monitor content
!  --------------------------
   integer            :: NO_OF_VARIABLES
   integer, parameter :: NO_OF_VARIABLES_Sij = 9
   integer, parameter :: NO_OF_FAVRE_VARS    = 6
   integer, parameter ::  U  = 1
   integer, parameter ::  V  = 2
   integer, parameter ::  W  = 3
   integer, parameter ::  UU = 4
   integer, parameter ::  VV = 5
   integer, parameter ::  WW = 6
   integer, parameter ::  UV = 7
   integer, parameter ::  UW = 8
   integer, parameter ::  VW = 9
   logical, save      :: stats_doReynolds = .true.
   logical, save      :: stats_doFavre    = .false.

   type StatisticsMonitor_t
      integer        :: state
      integer        :: sampling_interval
      integer        :: dump_interval
      integer        :: reset_interval
      integer        :: starting_iteration
      real(kind=RP)  :: starting_time
      integer        :: no_of_samples
      logical        :: saveGradients
      contains
         procedure   :: Construct    => StatisticsMonitor_Construct
         procedure   :: Update       => StatisticsMonitor_Update
         procedure   :: UpdateValues => StatisticsMonitor_UpdateValues
         procedure   :: GetState     => StatisticsMonitor_GetState
         procedure   :: WriteLabel   => StatisticsMonitor_WriteLabel
         procedure   :: WriteValue   => StatisticsMonitor_WriteValue
         procedure   :: WriteFile    => StatisticsMonitor_WriteFile
   end type StatisticsMonitor_t
!
!  ========
   contains
!  ========
!
!//////////////////////////////////////////////////////////////////////////////////
!
      subroutine StatisticsMonitor_Construct(self, mesh, saveGradients)
         use ParamfileRegions
         use PhysicsStorage, only: NCONS, NGRAD
         use HexMeshClass,   only: no_of_stats_variables, no_of_reynolds_variables, no_of_favre_variables
         implicit none
         class(StatisticsMonitor_t)    :: self
         class(HexMesh)                :: mesh
         logical, intent(in)           :: saveGradients
         integer, allocatable          :: Nsample, i0, Ndump, Nreset
         real(kind=RP), allocatable    :: t0
         integer                       :: eID
         character(len=LINE_LENGTH)    :: paramFile
         character(len=LINE_LENGTH)    :: averaging

         self % saveGradients = saveGradients
!
!        Search for the parameters in the case file
!        ------------------------------------------
         call get_command_argument(1, paramFile)
!
!        Read averaging mode (must be before the early-return below)
!        ------------------------------------------------------------
         stats_doReynolds = .true.
         stats_doFavre    = .false.
         averaging        = ""
         call readCharacterValueInRegion(trim(paramFile), "averaging", averaging, "#define statistics", "#end")
         call toLower(averaging)
         if (len_trim(averaging) > 0) then
            stats_doReynolds = (index(averaging, "reynolds") > 0)
            stats_doFavre    = (index(averaging, "favre")    > 0)
            if (.not. stats_doReynolds .and. .not. stats_doFavre) stats_doReynolds = .true.
         end if

         no_of_reynolds_variables = merge(NO_OF_VARIABLES_Sij, 0, stats_doReynolds)
         no_of_favre_variables    = merge(NO_OF_FAVRE_VARS,    0, stats_doFavre)

         NO_OF_VARIABLES = no_of_reynolds_variables + NCONS + no_of_favre_variables
         if (saveGradients) NO_OF_VARIABLES = NO_OF_VARIABLES + NGRAD * NDIM
         no_of_stats_variables = NO_OF_VARIABLES
!
!        Search if the #define statistics section is defined
!        ---------------------------------------------------
         if ( findIfStatsAreActive(paramFile) ) then
            self % state = WAITING_STATE

         else
            self % state = OFF
            return

         end if

         call readValueInRegion ( trim ( paramFile ) , "sampling interval"  , Nsample , "#define statistics" , "#end" )
         call readValueInRegion ( trim ( paramFile ) , "starting iteration" , i0      , "#define statistics" , "#end" )
         call readValueInRegion ( trim ( paramFile ) , "starting time"      , t0      , "#define statistics" , "#end" )
         call readValueInRegion ( trim ( paramFile ) , "dump interval"      , Ndump   , "#define statistics" , "#end" )
         call readValueInRegion ( trim ( paramFile ) , "reset interval"     , Nreset  , "#define statistics" , "#end" )
         
!
!        Check the input data
!        --------------------
         if ( allocated(i0) ) self % starting_iteration = i0      ! Both iteration and time can start the statistics calculation
         if ( allocated(t0) ) self % starting_time = t0           !

         if ( allocated(i0) .and. (.not. allocated(t0)) ) self % starting_time = huge(1.0_RP) ! Statistics controlled by iterations
         if ( allocated(t0) .and. (.not. allocated(i0)) ) self % starting_iteration = huge(1) ! Statistics controlled by time

         if ( (.not. allocated(i0)) .and. (.not. allocated(t0)) ) then  ! 
            self % starting_time = 0.0_RP                               ! Start since the beginning by default
            self % starting_iteration = 0                               !
         end if                                                         !

         if ( allocated(Nsample) ) then            !
            self % sampling_interval = Nsample     !  Set 1 by default
         else                                      !
            self % sampling_interval = 1           ! 
         end if     
         
         if ( allocated(Ndump) ) then              !
            self % dump_interval = Ndump           !  Set huge by default
         else                                      !
            self % dump_interval = huge(1)         ! 
         end if                                    !
         
         if ( allocated(Nreset) ) then                                        !
            if (allocated(Ndump) ) then                                       !
               self % dump_interval  = GreatestCommonDivisor(Nreset,Ndump)    !  Set huge by default
            else                                                              !
               self % dump_interval  = Nreset                                 !
            end if                                                            !
            self % reset_interval = Nreset                                    !
         else                                                                 !
            self % reset_interval = huge(1)                                   ! 
         end if                                                               !
!
!        Initial state: waiting
!        ----------------------
         self % state = WAITING_STATE
         self % no_of_samples = 0
!
!        Construct the data
!        ------------------
         !$acc enter data copyin(self) 
         do eID = 1, mesh % no_of_elements
            call mesh % elements(eID) % storage % stats % Construct(NO_OF_VARIABLES, mesh % elements(eID) % Nxyz)
         end do

      end subroutine StatisticsMonitor_Construct

      subroutine StatisticsMonitor_WriteFile(self, mesh, iter, t, solution_file)
         implicit none
         class(StatisticsMonitor_t) :: self
         class(HexMesh)             :: mesh
         integer, intent(in)        :: iter
         real(kind=RP), intent(in)  :: t
         character(len=*), intent(in)  :: solution_file
         character(len=LINE_LENGTH)    :: fileName

         write(fileName,'(A,A)') trim(solution_file),'.stats.hsol'
         if ( self % state .ne. OFF) call mesh % SaveStatistics(iter, t, trim(fileName), self % saveGradients)

      end subroutine StatisticsMonitor_WriteFile

      subroutine StatisticsMonitor_Update(self, mesh, iter, t, solution_file)
         implicit none
         class(StatisticsMonitor_t) :: self
         class(HexMesh)             :: mesh
         integer, intent(in)        :: iter
         real(kind=RP), intent(in)  :: t
         character(len=*), intent(in)  :: solution_file
!
!        ---------------
!        Local variables
!        ---------------
!
         logical     :: reset, dump
         character(len=LINE_LENGTH)    :: fileName
   
         if ( self % state .eq. OFF ) return
         if ( mod(iter, self % sampling_interval) .ne. 0) return
!
!        Read the parameter file to check the state
!        ------------------------------------------
         call self % GetState(reset, dump)
!
!        Dump the contents if requested
!        ------------------------------
         call Stopwatch % Pause("Solver") ! We dont want to measure the time of the statistics dump
         if ( dump .or. ( (mod(iter, self % dump_interval) == 0) .and. (iter > self % starting_iteration .or. t > self % starting_time) ) ) then
            write(fileName,'(A,A,I10.10,A)') trim(solution_file),'.stats.',iter,'.hsol'
            call mesh % SaveStatistics(iter, t, trim(fileName), self % saveGradients)
            if ( MPI_Process % isRoot ) write(STD_OUT,'(A,A,A)') '   *** Saving statistics file as "',trim(fileName),'".'
         end if
         call Stopwatch % Start("Solver") ! We dont want to measure the time of the statistics dump

!
!        Reset the statistics if requested
!        ---------------------------------
         if ( reset .or. ( (mod(iter, self % reset_interval) == 0) .and. (iter > self % starting_iteration) ) ) then
            call mesh % ResetStatistics
            self % no_of_samples = 0
         end if
!
!        Update the state if it is waiting
!        ---------------------------------
         if ( self % state .eq. WAITING_STATE ) then
            if ( (t .gt. self % starting_time) .or. (iter .ge. self % starting_iteration) ) then
               self % state = ACTIVE_STATE
            end if
         end if
            
!        Update the values if the state is "active"
!        ------------------------------------------
         if(self % state .eq. ACTIVE_STATE) call self % UpdateValues(mesh)

      end subroutine StatisticsMonitor_Update

      subroutine StatisticsMonitor_UpdateValues(self, mesh)
         use PhysicsStorage
         use HexMeshClass, only: no_of_reynolds_variables, no_of_favre_variables
         implicit none
         class(StatisticsMonitor_t)    :: self
         class(HexMesh)              :: mesh
!
!        ---------------
!        Local variables
!        ---------------
!
         integer  :: eID
         integer  :: i, j, k, eq
         real(RP) :: ratio, inv_nsamples_plus_1
         real(RP) :: rfactor1, rfactor2
         integer  :: nR, nF
         integer, dimension(5) :: limits

#ifdef NAVIERSTOKES
            nR = no_of_reynolds_variables
            nF = no_of_favre_variables

            inv_nsamples_plus_1 = 1.0_RP / (self % no_of_samples + 1)
            ratio = self % no_of_samples * inv_nsamples_plus_1

            if (stats_doReynolds) then
            !$acc parallel loop gang vector_length(128) present(mesh) firstprivate(inv_nsamples_plus_1, ratio) async(1)
            do eID = 1, size(mesh % elements)
               !$acc loop vector collapse(3) private(rfactor1, rfactor2)
               do k = 0, mesh % elements(eID) % Nxyz(3)   ; do j = 0, mesh % elements(eID) % Nxyz(2)    ; do i = 0, mesh % elements(eID) % Nxyz(1)
                  rfactor1 = inv_nsamples_plus_1 / mesh % elements(eID) % storage % Q(IRHO,i,j,k)
                  rfactor2 = inv_nsamples_plus_1 / POW2( mesh % elements(eID) % storage % Q(IRHO,i,j,k) )
                  mesh % elements(eID) % storage % stats % data(U,i,j,k)  = mesh % elements(eID) % storage % stats % data(U,i,j,k) &
                                                                          * ratio + mesh % elements(eID) % storage % Q(IRHOU,i,j,k) * rfactor1

                  mesh % elements(eID) % storage % stats % data(V,i,j,k)  = mesh % elements(eID) % storage % stats % data(V,i,j,k) &
                                                                          * ratio + mesh % elements(eID) % storage % Q(IRHOV,i,j,k) * rfactor1

                  mesh % elements(eID) % storage % stats % data(W,i,j,k)  = mesh % elements(eID) % storage % stats % data(W,i,j,k) &
                                                                          * ratio + mesh % elements(eID) % storage % Q(IRHOW,i,j,k) * rfactor1

                  mesh % elements(eID) % storage % stats % data(UU,i,j,k) = mesh % elements(eID) % storage % stats % data(UU,i,j,k) &
                                                                          * ratio + POW2( mesh % elements(eID) % storage % Q(IRHOU,i,j,k) ) * rfactor2

                  mesh % elements(eID) % storage % stats % data(VV,i,j,k) = mesh % elements(eID) % storage % stats % data(VV,i,j,k) &
                                                                          * ratio + POW2( mesh % elements(eID) % storage % Q(IRHOV,i,j,k) ) * rfactor2

                  mesh % elements(eID) % storage % stats % data(WW,i,j,k) = mesh % elements(eID) % storage % stats % data(WW,i,j,k) &
                                                                          * ratio + POW2( mesh % elements(eID) % storage % Q(IRHOW,i,j,k) ) * rfactor2

               end do                  ; end do                   ; end do
            end do
            !$acc end parallel loop

            !$acc parallel loop gang vector_length(128) present(mesh) firstprivate(inv_nsamples_plus_1, ratio) async(2)
            do eID = 1, size(mesh % elements)
               !$acc loop vector collapse(3) private(rfactor2)
               do k = 0, mesh % elements(eID) % Nxyz(3)   ; do j = 0, mesh % elements(eID) % Nxyz(2)    ; do i = 0, mesh % elements(eID) % Nxyz(1)
                  rfactor2 = inv_nsamples_plus_1 / POW2( mesh % elements(eID) % storage % Q(IRHO,i,j,k) )

                  mesh % elements(eID) % storage % stats % data(UV,i,j,k) = mesh % elements(eID) % storage % stats % data(UV,i,j,k) * ratio &
                                                                          + mesh % elements(eID) % storage % Q(IRHOU,i,j,k) * mesh % elements(eID) % storage % Q(IRHOV,i,j,k) * rfactor2

                  mesh % elements(eID) % storage % stats % data(UW,i,j,k) = mesh % elements(eID) % storage % stats % data(UW,i,j,k) * ratio &
                                                                          + mesh % elements(eID) % storage % Q(IRHOU,i,j,k) * mesh % elements(eID) % storage % Q(IRHOW,i,j,k) * rfactor2

                  mesh % elements(eID) % storage % stats % data(VW,i,j,k) = mesh % elements(eID) % storage % stats % data(VW,i,j,k) * ratio &
                                                                          + mesh % elements(eID) % storage % Q(IRHOV,i,j,k) * mesh % elements(eID) % storage % Q(IRHOW,i,j,k) * rfactor2

               end do                  ; end do                   ; end do
            end do
            !$acc end parallel loop
            end if ! stats_doReynolds

            !$acc parallel loop gang vector_length(128) present(mesh) firstprivate(inv_nsamples_plus_1, ratio, nR) async(3)
            do eID = 1, size(mesh % elements)
               !$acc loop vector collapse(4)
               do k = 0, mesh % elements(eID) % Nxyz(3)   ; do j = 0, mesh % elements(eID) % Nxyz(2)    ; do i = 0, mesh % elements(eID) % Nxyz(1)   ; do eq = 1, NCONS

                  mesh % elements(eID) % storage % stats % data(nR + eq,i,j,k) = mesh % elements(eID) % storage % stats % data(nR + eq,i,j,k) * ratio &
                                                                                           + mesh % elements(eID) % storage % Q(eq,i,j,k) * inv_nsamples_plus_1
               end do                  ; end do                   ; end do          ; end do
            end do
            !$acc end parallel loop

            if (stats_doFavre) then
            !$acc parallel loop gang vector_length(128) present(mesh) firstprivate(inv_nsamples_plus_1, ratio, nR) async(7)
            do eID = 1, size(mesh % elements)
               !$acc loop vector collapse(3) private(rfactor1)
               do k = 0, mesh % elements(eID) % Nxyz(3)   ; do j = 0, mesh % elements(eID) % Nxyz(2)    ; do i = 0, mesh % elements(eID) % Nxyz(1)
                  rfactor1 = inv_nsamples_plus_1 / mesh % elements(eID) % storage % Q(IRHO,i,j,k)
                  ! FUU = <(rho*u)^2 / rho> -> allows Favre stress <rho*u''u''> = FUU - <rho*u>^2/<rho>
                  mesh % elements(eID) % storage % stats % data(nR+NCONS+1,i,j,k) = mesh % elements(eID) % storage % stats % data(nR+NCONS+1,i,j,k) * ratio &
                     + POW2( mesh % elements(eID) % storage % Q(IRHOU,i,j,k) ) * rfactor1
                  ! FVV
                  mesh % elements(eID) % storage % stats % data(nR+NCONS+2,i,j,k) = mesh % elements(eID) % storage % stats % data(nR+NCONS+2,i,j,k) * ratio &
                     + POW2( mesh % elements(eID) % storage % Q(IRHOV,i,j,k) ) * rfactor1
                  ! FWW
                  mesh % elements(eID) % storage % stats % data(nR+NCONS+3,i,j,k) = mesh % elements(eID) % storage % stats % data(nR+NCONS+3,i,j,k) * ratio &
                     + POW2( mesh % elements(eID) % storage % Q(IRHOW,i,j,k) ) * rfactor1
                  ! FUV
                  mesh % elements(eID) % storage % stats % data(nR+NCONS+4,i,j,k) = mesh % elements(eID) % storage % stats % data(nR+NCONS+4,i,j,k) * ratio &
                     + mesh % elements(eID) % storage % Q(IRHOU,i,j,k) * mesh % elements(eID) % storage % Q(IRHOV,i,j,k) * rfactor1
                  ! FUW
                  mesh % elements(eID) % storage % stats % data(nR+NCONS+5,i,j,k) = mesh % elements(eID) % storage % stats % data(nR+NCONS+5,i,j,k) * ratio &
                     + mesh % elements(eID) % storage % Q(IRHOU,i,j,k) * mesh % elements(eID) % storage % Q(IRHOW,i,j,k) * rfactor1
                  ! FVW
                  mesh % elements(eID) % storage % stats % data(nR+NCONS+6,i,j,k) = mesh % elements(eID) % storage % stats % data(nR+NCONS+6,i,j,k) * ratio &
                     + mesh % elements(eID) % storage % Q(IRHOV,i,j,k) * mesh % elements(eID) % storage % Q(IRHOW,i,j,k) * rfactor1
               end do                  ; end do                   ; end do
            end do
            !$acc end parallel loop
            end if ! stats_doFavre

            if (self % saveGradients) then
            !$acc parallel loop gang vector_length(128) present(mesh) firstprivate(inv_nsamples_plus_1, ratio, nR, nF) async(4)
            do eID = 1, size(mesh % elements)
               !$acc loop vector collapse(4)
               do k = 0, mesh % elements(eID) % Nxyz(3)   ; do j = 0, mesh % elements(eID) % Nxyz(2)    ; do i = 0, mesh % elements(eID) % Nxyz(1)    ; do eq = 1, NGRAD

                     mesh % elements(eID) % storage % stats % data(nR+NCONS+nF+eq,i,j,k) = mesh % elements(eID) % storage % stats % data(nR+NCONS+nF+eq,i,j,k) * ratio &
                                                                                                + mesh % elements(eID) % storage % U_x(eq,i,j,k) * inv_nsamples_plus_1

               end do                  ; end do                   ; end do          ; end do
            end do
            !$acc end parallel loop

            !$acc parallel loop gang vector_length(128) present(mesh) firstprivate(inv_nsamples_plus_1, ratio, nR, nF) async(5)
            do eID = 1, size(mesh % elements)
               !$acc loop vector collapse(4)
               do k = 0, mesh % elements(eID) % Nxyz(3)   ; do j = 0, mesh % elements(eID) % Nxyz(2)    ; do i = 0, mesh % elements(eID) % Nxyz(1)    ; do eq = 1, NGRAD

                     mesh % elements(eID) % storage % stats % data(nR+NCONS+nF+NGRAD+eq,i,j,k) = mesh % elements(eID) % storage % stats % data(nR+NCONS+nF+NGRAD+eq,i,j,k) * ratio &
                                                                                                + mesh % elements(eID) % storage % U_y(eq,i,j,k) * inv_nsamples_plus_1
               end do                  ; end do                   ; end do          ; end do
            end do
            !$acc end parallel loop

            !$acc parallel loop gang vector_length(128) present(mesh) firstprivate(inv_nsamples_plus_1, ratio, nR, nF) async(6)
            do eID = 1, size(mesh % elements)
               !$acc loop vector collapse(4)
               do k = 0, mesh % elements(eID) % Nxyz(3)   ; do j = 0, mesh % elements(eID) % Nxyz(2)    ; do i = 0, mesh % elements(eID) % Nxyz(1)    ; do eq = 1, NGRAD

                     mesh % elements(eID) % storage % stats % data(nR+NCONS+nF+2*NGRAD+eq,i,j,k) = mesh % elements(eID) % storage % stats % data(nR+NCONS+nF+2*NGRAD+eq,i,j,k) * ratio &
                                                                                                + mesh % elements(eID) % storage % U_z(eq,i,j,k) * inv_nsamples_plus_1
               end do                  ; end do                   ; end do          ; end do
            end do
            !$acc end parallel loop

            end if

            !$acc wait
#endif

#ifdef INCNS
            nR = no_of_reynolds_variables
            nF = no_of_favre_variables

            inv_nsamples_plus_1 = 1.0_RP / (self % no_of_samples + 1)
            ratio = self % no_of_samples * inv_nsamples_plus_1

            do eID = 1, size(mesh % elements)
               associate(e    => mesh % elements(eID), &
                         data => mesh % elements(eID) % storage % stats % data)

               do k = 0, e % Nxyz(3)   ; do j = 0, e % Nxyz(2)    ; do i = 0, e % Nxyz(1)
                  rfactor1 = inv_nsamples_plus_1 / e % storage % Q(INSRHO,i,j,k)
                  rfactor2 = inv_nsamples_plus_1 / POW2( e % storage % Q(INSRHO,i,j,k) )
                  if (stats_doReynolds) then
                     data(U,i,j,k)  = data(U,i,j,k)  * ratio + e % storage % Q(INSRHOU,i,j,k) * rfactor1
                     data(V,i,j,k)  = data(V,i,j,k)  * ratio + e % storage % Q(INSRHOV,i,j,k) * rfactor1
                     data(W,i,j,k)  = data(W,i,j,k)  * ratio + e % storage % Q(INSRHOW,i,j,k) * rfactor1
                     data(UU,i,j,k) = data(UU,i,j,k) * ratio + POW2( e % storage % Q(INSRHOU,i,j,k) ) * rfactor2
                     data(VV,i,j,k) = data(VV,i,j,k) * ratio + POW2( e % storage % Q(INSRHOV,i,j,k) ) * rfactor2
                     data(WW,i,j,k) = data(WW,i,j,k) * ratio + POW2( e % storage % Q(INSRHOW,i,j,k) ) * rfactor2
                     data(UV,i,j,k) = data(UV,i,j,k) * ratio + e % storage % Q(INSRHOU,i,j,k) * e % storage % Q(INSRHOV,i,j,k) * rfactor2
                     data(UW,i,j,k) = data(UW,i,j,k) * ratio + e % storage % Q(INSRHOU,i,j,k) * e % storage % Q(INSRHOW,i,j,k) * rfactor2
                     data(VW,i,j,k) = data(VW,i,j,k) * ratio + e % storage % Q(INSRHOV,i,j,k) * e % storage % Q(INSRHOW,i,j,k) * rfactor2
                  end if
                  data(nR+1:nR+NCONS,i,j,k) = data(nR+1:nR+NCONS,i,j,k) * ratio + e % storage % Q(:,i,j,k) * inv_nsamples_plus_1
                  if (stats_doFavre) then
                     data(nR+NCONS+1,i,j,k) = data(nR+NCONS+1,i,j,k) * ratio + POW2( e % storage % Q(INSRHOU,i,j,k) ) * rfactor1
                     data(nR+NCONS+2,i,j,k) = data(nR+NCONS+2,i,j,k) * ratio + POW2( e % storage % Q(INSRHOV,i,j,k) ) * rfactor1
                     data(nR+NCONS+3,i,j,k) = data(nR+NCONS+3,i,j,k) * ratio + POW2( e % storage % Q(INSRHOW,i,j,k) ) * rfactor1
                     data(nR+NCONS+4,i,j,k) = data(nR+NCONS+4,i,j,k) * ratio + e % storage % Q(INSRHOU,i,j,k) * e % storage % Q(INSRHOV,i,j,k) * rfactor1
                     data(nR+NCONS+5,i,j,k) = data(nR+NCONS+5,i,j,k) * ratio + e % storage % Q(INSRHOU,i,j,k) * e % storage % Q(INSRHOW,i,j,k) * rfactor1
                     data(nR+NCONS+6,i,j,k) = data(nR+NCONS+6,i,j,k) * ratio + e % storage % Q(INSRHOV,i,j,k) * e % storage % Q(INSRHOW,i,j,k) * rfactor1
                  end if
                  if (self % saveGradients) then
                      data(nR+NCONS+nF+1:nR+NCONS+nF+NGRAD,i,j,k)           = data(nR+NCONS+nF+1:nR+NCONS+nF+NGRAD,i,j,k)           * ratio + e % storage % U_x(:,i,j,k) * inv_nsamples_plus_1
                      data(nR+NCONS+nF+NGRAD+1:nR+NCONS+nF+2*NGRAD,i,j,k)   = data(nR+NCONS+nF+NGRAD+1:nR+NCONS+nF+2*NGRAD,i,j,k)   * ratio + e % storage % U_y(:,i,j,k) * inv_nsamples_plus_1
                      data(nR+NCONS+nF+2*NGRAD+1:nR+NCONS+nF+3*NGRAD,i,j,k) = data(nR+NCONS+nF+2*NGRAD+1:nR+NCONS+nF+3*NGRAD,i,j,k) * ratio + e % storage % U_z(:,i,j,k) * inv_nsamples_plus_1
                  end if
               end do                  ; end do                   ; end do

               end associate
            end do
#endif

         self % no_of_samples = self % no_of_samples + 1
         
      end subroutine StatisticsMonitor_UpdateValues

      subroutine StatisticsMonitor_GetState(self, reset, dump)
         use ParamfileRegions
         use MPI_Process_Info
         implicit none
         class(StatisticsMonitor_t)    :: self
         logical, intent(out)          :: reset, dump
!
!        ---------------
!        Local variables
!        ---------------
!
         integer                       :: fID, io, ierr
         integer                       :: no_of_lines
         character(len=LINE_LENGTH)    :: paramFile
         character(len=LINE_LENGTH)    :: line
         logical                       :: isInside
         logical                       :: hasStart, hasPause, hasStop, hasDump, hasReset
         integer                       :: rstInt, dumpInt

         reset = .false.
         dump  = .false.
!
!        Only root will read, and it will broadcast the result to the rest
!        -----------------------------------------------------------------
         if ( MPI_Process % isRoot ) then
   
            call get_command_argument(1, paramFile)
!   
!           Open the file and get to the statistics zone
!           --------------------------------------------
            open( newunit = fID, file = trim(paramFile), status = "old", action = "read" )
      
            isInside = .false.
   
            no_of_lines = 0
   
            hasStart = .false.
            hasPause = .false.
            hasStop  = .false.
            hasReset = .false.
            hasDump  = .false.
   
            do
               read(fID,'(A)',iostat=io) line
               if ( io .ne. 0 ) exit
               no_of_lines = no_of_lines + 1 
!   
!              Check whether the line is in or out of the statistics region
!              ------------------------------------------------------------
               if ( index(getSquashedLine(trim(adjustl(line))), "#definestatistics" ) .ne. 0) then
                  isInside = .true.
                  cycle
               end if
   
               if ( isInside .and. (index(getSquashedLine(trim(adjustl(line))),"#end") .ne. 0) ) then
                  isInside = .false.
                  cycle
               end if
!   
!              Cycle if it is outside
!              ----------------------        
               if ( .not. isInside ) cycle
!   
!              It is inside, check whether the marks are present
!              -------------------------------------------------
               select case ( trim(adjustl(line)) )
               
               case(START_COMMAND)
                  if ( index(line,"*") .eq. 0 ) hasStart = .true.
   
               case(PAUSE_COMMAND)
                  if ( index(line,"*") .eq. 0 ) hasPause = .true.
   
               case(STOP_COMMAND)
                  if ( index(line,"*") .eq. 0 ) hasStop = .true.
   
               case(RESET_COMMAND)
                  if ( index(line,"*") .eq. 0 ) hasReset = .true.
   
               case(DUMP_COMMAND)
                  if ( index(line,"*") .eq. 0 ) hasDump = .true.
   
               end select
   
            end do
!   
!           Update the status of the statistics monitor
!           -------------------------------------------
            reset = hasReset
            dump  = hasDump
            
            if ( (self % state .eq. WAITING_STATE) .or. (self % state .eq. ACTIVE_STATE) ) then
               if ( hasStop ) then
                  dump = .true.
                  reset = .true.
                  self % state = STOP_STATE
               end if
   
               if ( hasPause ) then
                  self % state = STOP_STATE
               end if
   
               if ( hasStart ) then
                  self % state = ACTIVE_STATE
               end if
   
            elseif ( self % state .eq. STOP_STATE ) then
               if ( hasStart )  self % state = ACTIVE_STATE
           
            end if
   
            close(fID)
!   
!           Rewrite the control file to remove the controllers      
!           --------------------------------------------------
            if ( hasStart .or. hasPause .or. hasStop .or. hasReset .or. hasDump ) call rewriteControlFile(paramFile, no_of_lines)
         end if
!
!        Broadcast the result to the rest of the processes
!        -------------------------------------------------
         if ( MPI_Process % DoMPIAction ) then
#ifdef _HAS_MPI_
            if ( reset ) then
               rstInt = 1
            else
               rstInt = 0
            end if

            if ( dump ) then
               dumpInt = 1
            else
               dumpInt = 0
            end if

            call mpi_bcast(rstInt, 1, MPI_INT, 0, MPI_COMM_WORLD, ierr)
#endif
         end if
         

      end subroutine StatisticsMonitor_GetState

      logical function findIfStatsAreActive(paramFile)
         use ParamfileRegions
         implicit none
         character(len=*), intent(in)     :: paramFile
!
!        ---------------
!        Local variables
!        ---------------
!
         integer     :: fID, io
         character(len=LINE_LENGTH)    :: line

         open(newunit = fID, file = trim(paramFile), action = "read", status = "old" )
         
         findIfStatsAreActive = .false.

         do
            read(fID,'(A)', iostat = io) line
            if ( io .ne. 0 ) exit
      
            if ( index(getSquashedLine(line), '#definestatistics') .ne. 0 ) then
               findIfStatsAreActive = .true.
               close(fID)
               return
            end if
            
         end do

         close(fID)
         
      end function findIfStatsAreActive

      subroutine rewriteControlFile(paramFile, no_of_lines)
         implicit none
         character(len=*), intent(in)  :: paramFile
         integer,          intent(in)  :: no_of_lines
!
!        ---------------
!        Local variables
!        ---------------
!
         integer     :: fID,i
         character(len=LINE_LENGTH) :: lines(no_of_lines)

         open(newunit = fID, file = trim(paramFIle), action = "read", status = "old" )

         
         do i = 1, no_of_lines

            read(fID,'(A)') lines(i)
            
            select case (trim(adjustl(lines(i))))
            
            case(START_COMMAND)
               lines(i) = "   " // START_COMMAND // "*"
            
            case(PAUSE_COMMAND)
               lines(i) = "   " // PAUSE_COMMAND // "*"

            case(STOP_COMMAND)
               lines(i) = "   " // STOP_COMMAND // "*"

            case(DUMP_COMMAND)
               lines(i) = "   " // DUMP_COMMAND // "*"

            case(RESET_COMMAND)
               lines(i) = "   " // RESET_COMMAND // "*"

            end select
         end do

         close(fID)
   
         open(newunit = fID, file = trim(paramFile), action = "write", status = "old" )

         do i = 1, no_of_lines
            write(fID,'(A)') trim(lines(i))
         end do

         close(fID)

      end subroutine rewriteControlFile

      subroutine StatisticsMonitor_WriteLabel ( self )
         implicit none
         class(StatisticsMonitor_t)             :: self

         if ( self % state .eq. OFF ) return
         write(STD_OUT , '(3X,A10)' , advance = "no") "Stats."

      end subroutine StatisticsMonitor_WriteLabel

      subroutine StatisticsMonitor_WriteValue ( self ) 
!
!        *************************************************************
!              This subroutine writes the monitor value for the time
!           integrator Display procedure.
!        *************************************************************
!
         implicit none
         class(StatisticsMonitor_t) :: self
         integer                 :: bufferLine
         character(len=10)       :: state

         if ( self % state .eq. OFF ) return

         select case ( self % state)
         case(ACTIVE_STATE)
            state = "Active"
   
         case(WAITING_STATE)
            state = "Waiting"

         case(STOP_STATE)
            state = "Inactive"   

         end select
   
         write(STD_OUT , '(1X,A,1X,A10)' , advance = "no") "|" , trim(state)

      end subroutine StatisticsMonitor_WriteValue 

end module StatisticsMonitor
#endif