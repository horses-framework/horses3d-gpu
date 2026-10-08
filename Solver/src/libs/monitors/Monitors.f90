#include "Includes.h"
module MonitorsClass
   use SMConstants
   use NodalStorageClass
   use HexMeshClass
   use MonitorDefinitions
   use ResidualsMonitorClass
   use VolumeMonitorClass
   use LoadBalancingMonitorClass
   use FileReadingUtilities      , only: getFileName
#ifdef FLOW
   use ProbeClass
   use PhysicsStorage
   use VariableConversion
   use FluidData
#endif
#if defined(NAVIERSTOKES) || defined(INCNS)
   use StatisticsMonitor
   use SurfaceMonitorClass
#endif
#ifdef HAS_HDF5
   use HDF5
#endif
   implicit none
!

   private
   public      Monitor_t

   integer, parameter :: FPVAR_UNKNOWN    = 0
   integer, parameter :: FPVAR_PRESSURE   = 1
   integer, parameter :: FPVAR_VELOCITY   = 2
   integer, parameter :: FPVAR_U          = 3
   integer, parameter :: FPVAR_V          = 4
   integer, parameter :: FPVAR_W          = 5
   integer, parameter :: FPVAR_MACH       = 6
   integer, parameter :: FPVAR_K          = 7
   integer, parameter :: FPVAR_RHO        = 8
   integer, parameter :: FPVAR_STATICPRES = 9
   integer, parameter :: FPVAR_DENSITY    = 10
!
!  HDF5 persistent file handle for file-probe output.
!  Opening and closing the .probes.h5 file on every write is expensive
!  on network filesystems.  The file is created once in
!  Monitor_InitFileProbesHDF5 and kept open until Monitor_Destruct.
!
#ifdef HAS_HDF5
   integer(HID_T), save :: hdf5_fp_fid  = 0_HID_T
   logical,        save :: hdf5_fp_open = .false.
#endif
!
!  *****************************
!  Main monitor class definition
!  *****************************
!  
   type Monitor_t
      character(len=LINE_LENGTH)                 :: solution_file
      character(len=LINE_LENGTH)                 :: probes_solution_file = ""
      integer                                    :: no_of_probes
      integer                                    :: no_of_surfaceMonitors
      integer                                    :: no_of_volumeMonitors
      integer                                    :: no_of_loadBalancingMonitors
      integer                                    :: no_of_fileProbes = 0
      character(len=LINE_LENGTH)                 :: probesFileName = ""
      character(len=STR_LEN_MONITORS), allocatable :: probesVariables(:)
      real(kind=RP)                              :: probeFileSaveTimestep = 0.0_RP
      character(len=8)                           :: probeFileOutputFormat = "ASCII"
      real(kind=RP)                              :: fp_lastSavedTime = -huge(0.0_RP)
      integer                                    :: fp_nOwned = 0
      integer             , allocatable          :: fp_ownedIdx(:)
      real(kind=RP)       , allocatable          :: fp_buf(:)
!     fp_x/fp_fileUnit are the only per-file-probe state kept O(total) on
!     every rank - coordinates (needed by root for the HDF5 /coordinates
!     dataset and for ASCII headers) and the ASCII "keep file open" file
!     unit. Both are tiny (real/int arrays, not full Probe_t). Everything
!     else that used to live in Monitors % probes(:) for file-probes (one
!     Probe_t per probe, replicated on every rank regardless of ownership)
!     is now built directly by InitializeProbesFromFile into fp_nOwned/
!     fp_ownedIdx/fp_cpu_* (CPU) or fp_eID/fp_lxi/... (GPU) below, sized
!     O(owned probes on this rank) instead of O(total probes).
      real(kind=RP)       , allocatable          :: fp_x(:,:)
      integer             , allocatable          :: fp_fileUnit(:)
      logical             , allocatable          :: fp_active(:)
#ifndef _OPENACC
!     Compact SoA for CPU path: owned-probe data laid out contiguously to avoid
!     stride-112 scattered access into the main probes(:) array.
      integer             , allocatable          :: fp_cpu_eID(:)
      integer             , allocatable          :: fp_cpu_Nx(:)
      integer             , allocatable          :: fp_cpu_Ny(:)
      integer             , allocatable          :: fp_cpu_Nz(:)
      real(kind=RP)       , allocatable          :: fp_cpu_lxi(:,:)
      real(kind=RP)       , allocatable          :: fp_cpu_leta(:,:)
      real(kind=RP)       , allocatable          :: fp_cpu_lzeta(:,:)
      integer             , allocatable          :: fp_cpu_varCodes(:)
!     MPI_Gatherv infrastructure: replaces MPI_Allreduce(fp_buf) with a
!     many-to-one gather to root, eliminating the zero-fill of fp_buf and
!     reducing communication from all-to-all to many-to-one.
!     fp_owned_buf: compact result buffer, size fp_nOwned*nv (all ranks).
!     fp_gatherv_counts/displs: counts and byte-offsets for MPI_Gatherv (root).
!     fp_gatherv_perm: for each gathered entry k, the global probe index in fp_buf (root).
      real(kind=RP)       , allocatable          :: fp_owned_buf(:)
      integer             , allocatable          :: fp_gatherv_counts(:)
      integer             , allocatable          :: fp_gatherv_displs(:)
      integer             , allocatable          :: fp_gatherv_perm(:)
#endif
      integer                                    :: bufferLine
      integer                      , allocatable :: iter(:)
      integer                                    :: dt_restriction
      logical                                    :: write_dt_restriction
      real(kind=RP)                , allocatable :: t(:)
      real(kind=RP)                , allocatable :: SolverSimuTime(:)
      real(kind=RP)                , allocatable :: TotalSimuTime(:)
      type(Residuals_t)                          :: residuals
      class(VolumeMonitor_t)       , allocatable :: volumeMonitors(:)
      class(LoadBalancingMonitor_t), allocatable :: loadBalancingMonitors(:)
#ifdef FLOW
      class(Probe_t)               , allocatable :: probes(:)
#endif
#ifdef _OPENACC
      integer                                    :: fp_Nmax = 0
      integer             , allocatable          :: fp_eID(:)
      logical             , allocatable          :: fp_ownsProbe(:)
      real(kind=RP)       , allocatable          :: fp_lxi(:,:)
      real(kind=RP)       , allocatable          :: fp_leta(:,:)
      real(kind=RP)       , allocatable          :: fp_lzeta(:,:)
      integer             , allocatable          :: fp_varCodes(:)
      real(kind=RP)       , allocatable          :: fp_values_gpu(:,:)
#endif
#if defined(NAVIERSTOKES) || defined(INCNS)
      class(SurfaceMonitor_t)      , allocatable :: surfaceMonitors(:)
      type(StatisticsMonitor_t)                  :: stats
#endif
      contains
         procedure   :: Construct       => Monitors_Construct
         procedure   :: WriteLabel      => Monitor_WriteLabel
         procedure   :: WriteUnderlines => Monitor_WriteUnderlines
         procedure   :: WriteValues     => Monitor_WriteValues
         procedure   :: UpdateValues    => Monitor_UpdateValues
         procedure   :: WriteToFile     => Monitor_WriteToFile
         procedure   :: WritePostProcessingSummary => Monitor_WritePostProcessingSummary
         procedure   :: destruct        => Monitor_Destruct
         procedure   :: copy            => Monitor_Assign
         generic     :: assignment(=)   => copy
   end type Monitor_t
!
!  ========
   contains
!  ========
!
!///////////////////////////////////////////////////////////////////////////////////////
!
      subroutine Monitors_Construct( Monitors, mesh, controlVariables )
         use FTValueDictionaryClass
         use mainKeywordsModule
         use MPI_Process_Info
#ifdef _HAS_MPI_
         use mpi
#endif
         implicit none
         class(Monitor_t)                     :: Monitors
         class(HexMesh)                       :: mesh
         class(FTValueDictionary), intent(in) :: controlVariables
         
!
!        ---------------
!        Local variables
!        ---------------
!
         integer                         :: fID , io
         integer                         :: i
         character(len=STR_LEN_MONITORS) :: line
         character(len=STR_LEN_MONITORS) :: solution_file
         logical, save                   :: FirstCall = .TRUE.
         logical                         :: saveGradients
         character(len=LINE_LENGTH)      :: probesFileName
         character(len=LINE_LENGTH)      :: probesVariablesLine
         character(len=STR_LEN_MONITORS), allocatable :: probesVariables(:)
         integer                         :: no_of_fileProbes
         integer                         :: no_of_probesVariables
         real(kind=RP)                   :: probeFileSaveTimestep
         character(len=8)                :: probeFileOutputFormat
         character(len=LINE_LENGTH)      :: probes_solution_file
         character(len=LINE_LENGTH)      :: probes_dir
         integer                         :: last_slash, ierr
!
!        Setup the buffer
!        ----------------
         if (controlVariables % containsKey("monitors flush interval") ) then
            BUFFER_SIZE = controlVariables % integerValueForKey("monitors flush interval")
         end if
         
         if ( .not. allocated(Monitors % iter) ) then
            allocate ( Monitors % TotalSimuTime(BUFFER_SIZE), &
                       Monitors % SolverSimuTime(BUFFER_SIZE), &
                       Monitors % t(BUFFER_SIZE), &
                       Monitors % iter(BUFFER_SIZE) )
         end if
!
!        Get the solution file name
!        --------------------------
         solution_file = controlVariables % stringValueForKey( solutionFileNameKey, requestedLength = STR_LEN_MONITORS )
!
!        Remove the *.hsol termination
!        -----------------------------
         solution_file = trim(getFileName(solution_file))
         Monitors % solution_file = trim(solution_file)
!
!        Build the probes output subdirectory: <dir>/probes/<base>
!        ----------------------------------------------------------
         last_slash = index(trim(solution_file), '/', back=.true.)
         if (last_slash .gt. 0) then
            probes_dir           = trim(solution_file(1:last_slash)) // "probes"
            probes_solution_file = trim(solution_file(1:last_slash)) // "probes/" // trim(solution_file(last_slash+1:))
         else
            probes_dir           = "probes"
            probes_solution_file = "probes/" // trim(solution_file)
         end if
         if (MPI_Process % isRoot) then
            call execute_command_line('mkdir -p "' // trim(probes_dir) // '"', wait=.true.)
         end if
#ifdef _HAS_MPI_
         if ( MPI_Process % doMPIAction ) call MPI_Barrier(MPI_COMM_WORLD, ierr)
#endif
         Monitors % probes_solution_file = trim(probes_solution_file)
!
!        Search in case file for probes, surface monitors, and volume monitors
!        ---------------------------------------------------------------------
         no_of_fileProbes = 0
         if (mesh % child) then ! Return doing nothing if this is a child mesh
            Monitors % no_of_probes = 0
            Monitors % no_of_surfaceMonitors = 0
            Monitors % no_of_volumeMonitors = 0
            Monitors % no_of_loadBalancingMonitors = 0
         else
            call getNoOfMonitors( Monitors % no_of_probes, Monitors % no_of_surfaceMonitors, Monitors % no_of_volumeMonitors, Monitors % no_of_loadBalancingMonitors )
!
!           Check for an additional probes definition file (one probe per
!           line: "x y z variable1 [variable2 ...]"), allowing several
!           variables to be sampled and saved for the same probe location
!           ---------------------------------------------------------------
#ifdef FLOW
            call readProbesFileBlock( probesFileName, probesVariablesLine, no_of_probesVariables, probeFileSaveTimestep, probeFileOutputFormat )

            if ( len_trim(probesFileName) .gt. 0 ) then
!
!              Root-only count + broadcast: with O(1e3) MPI ranks all
!              opening/reading the same (possibly huge) probes file
!              independently, there is no guarantee every rank's count
!              agrees (seen in practice as a hung MPI_Allreduce deep
!              inside InitializeProbesFromFile, since a mismatched count
!              on even one rank breaks that collective). Reading it once
!              on root and broadcasting the result removes the ambiguity
!              entirely instead of hoping O(1e3) concurrent reads agree.
!              --------------------------------------------------------------
               if ( MPI_Process % isRoot ) then
                  call countProbesInFile( trim(probesFileName), no_of_fileProbes )
               end if
#ifdef _HAS_MPI_
               if ( MPI_Process % doMPIAction ) then
                  call MPI_Bcast(no_of_fileProbes, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
               end if
#endif

               if ( len_trim(probesVariablesLine) .gt. 0 ) then
                  call splitIntoTokens( trim(probesVariablesLine), probesVariables, no_of_probesVariables )
               else
                  write(STD_OUT,*) "Error: 'variables' must be specified inside the '#define probe file' block."
                  stop
               end if
            end if
#endif
         end if
!
!        Initialize the Monitors class in the GPU
!        ----------------------------------------
         !!$acc enter data copyin(Monitors)
         !$acc update device(Monitors)

!        Pro tip: This is necessary to avoid the compiler doing behind the scenes copies of the class.
!        Pro tip(cont): Its not really necessary for the class itself, but for the arrays inside the class.
!
!        Initialize
!        ----------
         call Monitors % residuals % Initialization( solution_file , FirstCall )

         allocate ( Monitors % volumeMonitors ( Monitors % no_of_volumeMonitors )  )
         !$acc update device(Monitors)
         do i = 1 , Monitors % no_of_volumeMonitors
            call Monitors % volumeMonitors(i) % Initialization ( mesh , i, solution_file , FirstCall  )
         end do

         allocate ( Monitors % loadBalancingMonitors ( Monitors % no_of_loadBalancingMonitors )  )
         !$acc update device(Monitors)
         do i = 1 , Monitors % no_of_loadBalancingMonitors
            call Monitors % loadBalancingMonitors(i) % Initialization ( mesh , i, solution_file , FirstCall )
         end do

#ifdef FLOW
         allocate ( Monitors % probes ( Monitors % no_of_probes )  )
         !$acc update device(Monitors)
         do i = 1 , Monitors % no_of_probes
            call Monitors % probes(i) % Initialization ( mesh , i, probes_solution_file , FirstCall )
         end do

!        File-probes (bulk "#define probe file" block) no longer get a
!        Probe_t entry each: InitializeProbesFromFile builds Monitors'
!        own SoA buffers (fp_x, fp_nOwned/fp_ownedIdx/fp_cpu_* or
!        fp_eID/fp_lxi/...) directly, sized O(owned probes on this rank)
!        instead of replicating a full Probe_t (two character(128)
!        strings plus several array descriptors each) for every one of
!        up to O(1e6) probes on every single MPI rank.
!        --------------------------------------------------------------
         if ( no_of_fileProbes .gt. 0 ) then
            Monitors % probesFileName         = trim(probesFileName)
            Monitors % probeFileSaveTimestep  = probeFileSaveTimestep
            Monitors % probeFileOutputFormat  = probeFileOutputFormat
            Monitors % fp_lastSavedTime       = -huge(0.0_RP)
            allocate( Monitors % probesVariables(size(probesVariables)) )
            Monitors % probesVariables = probesVariables

            call InitializeProbesFromFile( trim(probesFileName), Monitors, mesh, no_of_fileProbes )

#ifdef HAS_HDF5
            if ( trim(probeFileOutputFormat) .eq. "HDF5" ) then
               call Monitor_InitFileProbesHDF5( Monitors, no_of_fileProbes )
            end if
#endif
         end if

         Monitors % no_of_fileProbes = no_of_fileProbes
         Monitors % no_of_probes = Monitors % no_of_probes + no_of_fileProbes
#endif

#if defined(NAVIERSTOKES) || defined(INCNS)
         saveGradients    = controlVariables % logicalValueForKey(saveGradientsToSolutionKey)
         call Monitors % stats     % Construct(mesh, saveGradients)


         allocate ( Monitors % surfaceMonitors ( Monitors % no_of_surfaceMonitors )  )
         !$acc update device(Monitors)
         do i = 1 , Monitors % no_of_surfaceMonitors
            call Monitors % surfaceMonitors(i) % Initialization ( mesh , i, solution_file , FirstCall )
         end do
#endif

         Monitors % write_dt_restriction = controlVariables % logicalValueForKey( "write dt restriction" )
         
         Monitors % bufferLine = 0

         FirstCall = .FALSE.
!
!        Include the latest changes in the GPU
!        ----------------------------------------
         !$acc update device(Monitors)

      end subroutine Monitors_Construct

      subroutine Monitor_WriteLabel ( self )
!
!        ***************************************************
!           This subroutine prints the labels for the time
!         integrator Display procedure.
!        ***************************************************
!
         use MPI_Process_Info
         implicit none
         class(Monitor_t)              :: self
         integer                       :: i 
      
         if ( .not. MPI_Process % isRoot ) return
!
!        Write "Iteration" and "Time"
!        ----------------------------
         write ( STD_OUT , ' ( A10    ) ' , advance = "no" ) "Iteration"
         write ( STD_OUT , ' ( 3X,A10 ) ' , advance = "no" ) "Time"
!
!        Write residuals labels
!        ----------------------
         call self % residuals % WriteLabel
!
!        Write volume monitors labels
!        -----------------------------
         do i = 1 , self % no_of_volumeMonitors
            call self % volumeMonitors(i) % WriteLabel
         end do
!
!        Write load balancing monitor labels
!        ------------------------------------
         do i = 1 , self % no_of_loadBalancingMonitors
            call self % loadBalancingMonitors(i) % WriteLabel
         end do

#ifdef FLOW
!
!        Write probes labels (file-based probes excluded for readability)
!        ---------------------------------------------------------------
         do i = 1 , self % no_of_probes - self % no_of_fileProbes
            call self % probes(i) % WriteLabel
         end do
#endif

#if defined(NAVIERSTOKES) || defined(INCNS)
!
!        Write surface monitors labels
!        -----------------------------
         do i = 1 , self % no_of_surfaceMonitors
            call self % surfaceMonitors(i) % WriteLabel
         end do

         call self % stats % WriteLabel
#endif
!
!        Write label for dt restriction
!        ------------------------------
         if (self % write_dt_restriction) write ( STD_OUT , ' ( 3X,A10 ) ' , advance = "no" ) "dt restr."

         write(STD_OUT , *) 

      end subroutine Monitor_WriteLabel

      subroutine Monitor_WriteUnderlines( self ) 
!
!        ********************************************************
!              This subroutine displays the underlines for the
!           time integrator Display procedure.
!        ********************************************************
!
         use PhysicsStorage
         use MPI_Process_Info
         implicit none
         class(Monitor_t)                         :: self
!
!        ---------------
!        Local variables
!        ---------------
!
         integer                                  :: i, j
         character(len=MONITOR_LENGTH), parameter :: dashes = "----------"

         if ( .not. MPI_Process % isRoot ) return
!
!        Print dashes for "Iteration" and "Time"
!        ---------------------------------------
         write ( STD_OUT , ' ( A10    ) ' , advance = "no" ) trim ( dashes ) 
         write ( STD_OUT , ' ( 3X,A10 ) ' , advance = "no" ) trim ( dashes ) 
!
!        Print dashes for residuals
!        --------------------------
         do i = 1 , NCONS
            write(STD_OUT , '(3X,A10)' , advance = "no" ) trim(dashes)
         end do
!
!        Print dashes for volume monitors
!        --------------------------------
         do i = 1 , self % no_of_volumeMonitors  ; do j=1, size ( self % volumeMonitors(i) % values, 1 )
            write(STD_OUT , '(3X,A10)' , advance = "no" ) dashes(1 : min(10 , len_trim( self % volumeMonitors(i) % monitorName ) + 2 ) )
         end do                                  ; end do
!
!        Print dashes for load balancing monitor
!        --------------------------------------
         do i = 1 , self % no_of_loadBalancingMonitors ; do j=1, size ( self % loadBalancingMonitors(i) % values, 1 )
            write(STD_OUT , '(3X,A10)' , advance = "no" ) dashes(1 : min(10 , len_trim( self % loadBalancingMonitors(i) % monitorName ) + 2 ) )
         end do                                  ; end do

#ifdef FLOW
!
!        Print dashes for probes (file-based probes excluded for readability)
!        ----------------------------------------------------------------
         do i = 1 , self % no_of_probes - self % no_of_fileProbes
            if ( self % probes(i) % active ) then
               write(STD_OUT , '(3X,A10)' , advance = "no" ) dashes(1 : min(10 , len_trim( self % probes(i) % monitorName ) + 2 ) )
            end if
         end do
#endif

#if defined(NAVIERSTOKES) || defined(INCNS)
!
!        Print dashes for surface monitors
!        ---------------------------------
         do i = 1 , self % no_of_surfaceMonitors
            write(STD_OUT , '(3X,A10)' , advance = "no" ) dashes(1 : min(10 , len_trim( self % surfaceMonitors(i) % monitorName ) + 2 ) )
         end do

         if ( self % stats % state .ne. 0 ) write(STD_OUT,'(3X,A10)',advance="no") trim(dashes)
#endif

         
!
!        Print dashes for dt restriction
!        -------------------------------
         if (self % write_dt_restriction) write ( STD_OUT , ' ( 3X,A10 ) ' , advance = "no" ) trim ( dashes ) 
         
         write(STD_OUT , *) 

      end subroutine Monitor_WriteUnderlines

      subroutine Monitor_WriteValues ( self )
!
!        *******************************************************
!              This subroutine prints the values for the time
!           integrator Display procedure.
!        *******************************************************
!
         use MPI_Process_Info
         implicit none
         class(Monitor_t)           :: self
         integer                    :: i

         if ( .not. MPI_Process % isRoot ) return
!
!        Print iteration and time
!        ------------------------
         write ( STD_OUT , ' ( I10            ) ' , advance = "no" ) self % iter    ( self % bufferLine )
         write ( STD_OUT , ' ( 1X,A,1X,ES10.3 ) ' , advance = "no" ) "|" , self % t ( self % bufferLine ) 
!
!        Print residuals
!        ---------------
         call self % residuals % WriteValues( self % bufferLine )
!
!        Print volume monitors
!        ---------------------
         do i = 1 , self % no_of_volumeMonitors
            call self % volumeMonitors(i) % WriteValues ( self % bufferLine )
         end do
!
!        Print load balancing monitors
!        -----------------------------
         do i = 1 , self % no_of_loadBalancingMonitors
            call self % loadBalancingMonitors(i) % WriteValues ( self % bufferLine )
         end do

#ifdef FLOW
!
!        Print probes (file-based probes excluded for readability)
!        ----------------------------------------------------------
         do i = 1 , self % no_of_probes - self % no_of_fileProbes
            call self % probes(i) % WriteValues ( self % bufferLine )
         end do
#endif

#if defined(NAVIERSTOKES) || defined(INCNS)
!
!        Print surface monitors
!        ----------------------
         do i = 1 , self % no_of_surfaceMonitors
            call self % surfaceMonitors(i) % WriteValues ( self % bufferLine )
         end do

         call self % stats % WriteValue
#endif
!
!        Print dt restriction
!        --------------------
         if (self % write_dt_restriction) then
            select case (self % dt_restriction)
               case (DT_FIXED) ; write ( STD_OUT , ' ( 1X,A,1X,A10) ' , advance = "no" ) "|" , 'Fixed'
               case (DT_DIFF)  ; write ( STD_OUT , ' ( 1X,A,1X,A10) ' , advance = "no" ) "|" , 'Diffusive'
               case (DT_CONV)  ; write ( STD_OUT , ' ( 1X,A,1X,A10) ' , advance = "no" ) "|" , 'Convective'
            end select
         end if

         write(STD_OUT , *) 

      end subroutine Monitor_WriteValues

      subroutine Monitor_UpdateValues ( self, mesh, t , iter, maxResiduals, Autosave, dt )
!
!        ***************************************************************
!              This subroutine updates the values for the residuals,
!           for the probes, surface and volume monitors.
!        ***************************************************************
!        
         use PhysicsStorage
         use StopwatchClass
         use MPI_Process_Info
         implicit none
         class(Monitor_t)    :: self
         class(HexMesh)      :: mesh
         real(kind=RP)       :: t
         integer             :: iter
         real(kind=RP)       :: maxResiduals(NCONS), dt
         logical             :: Autosave
!
!        ---------------
!        Local variables
!        ---------------
!
         integer                       :: i

!
!        Move to next buffer line
!        ------------------------
         self % bufferLine = self % bufferLine + 1
!
!        Save time, iteration and CPU-time
!        -----------------------
         self % t       ( self % bufferLine )  = t
         self % iter    ( self % bufferLine )  = iter
         self % SolverSimuTime ( self % bufferLine )  = Stopwatch % ElapsedTime("Solver")
         self % TotalSimuTime ( self % bufferLine )  = Stopwatch % ElapsedTime("TotalTime")
!
!        Compute current residuals
!        -------------------------
         call self % residuals % Update( mesh, maxResiduals, self % bufferLine )
!
!        Update volume monitors
!        ----------------------
         do i = 1 , self % no_of_volumeMonitors
            call self % volumeMonitors(i) % Update( mesh , self % bufferLine )
         end do
!
!        Update load balancing monitors
!        ------------------------------
         do i = 1 , self % no_of_loadBalancingMonitors
            call self % loadBalancingMonitors(i) % Update( mesh , self % bufferLine )
         end do

#ifdef FLOW
!
!        Update standard probes (GPU path with per-probe MPI)
!        -----------------------------------------------------
         do i = 1 , self % no_of_probes - self % no_of_fileProbes
            call self % probes(i) % Update( mesh , self % bufferLine )
         end do
!
!        Update file probes: skip entirely if it is not yet time to save.
!        With O(1e6) probes and MPI, the Gatherv dominates; evaluating
!        every step when save_timestep >> dt wastes most of that cost.
!        ----------------------------------------------------------
         if ( self % no_of_fileProbes .gt. 0 ) then
            if ( self % probeFileSaveTimestep .le. 0.0_RP .or. &
                 t .ge. self % fp_lastSavedTime + self % probeFileSaveTimestep ) then
               call Monitor_UpdateFileProbes( self, mesh, 1 )
               call Monitor_FlushFileProbesNow( self, t, iter )
            end if
         end if
#endif

#if defined(NAVIERSTOKES) || defined(INCNS)
!
!        Update surface monitors
!        -----------------------
         do i = 1 , self % no_of_surfaceMonitors
            call self % surfaceMonitors(i) % Update( mesh , self % bufferLine, iter, autosave, dt )
         end do
!
!        Update statistics
!        -----------------
         call self % stats % Update(mesh, iter, t, trim(self % solution_file) )
#endif

!
!        Update dt restriction
!        ---------------------
         if (self % write_dt_restriction) self % dt_restriction = mesh % dt_restriction 
         
      end subroutine Monitor_UpdateValues

      subroutine Monitor_WriteToFile ( self , mesh, force) 
!
!        ******************************************************************
!              This routine has a double behaviour:
!           force = .true.  -> Writes to file and resets buffers
!           force = .false. -> Just writes to file if the buffer is full
!        ******************************************************************
!
         use MPI_Process_Info
         implicit none
         class(Monitor_t)        :: self
         class(HexMesh)          :: mesh
         logical, optional       :: force
!        ------------------------------------------------
         integer                 :: i 
         logical                 :: forceVal

         if ( present ( force ) ) then
            forceVal = force

         else
            forceVal = .false.

         end if

         if ( forceVal ) then 
!
!           In this case the monitors are exported to their files and the buffer is reset
!           -----------------------------------------------------------------------------
            call self % residuals % WriteToFile ( self % iter , self % t, self % TotalSimuTime, self % SolverSimuTime , self % bufferLine )
   
            do i = 1 , self % no_of_volumeMonitors
               call self % volumeMonitors(i) % WriteToFile ( self % iter , self % t , self % bufferLine )
            end do

            do i = 1 , self % no_of_loadBalancingMonitors
               call self % loadBalancingMonitors(i) % WriteToFile ( self % iter , self % t , self % bufferLine )
            end do

#ifdef FLOW
            do i = 1 , self % no_of_probes - self % no_of_fileProbes
               call self % probes(i) % WriteToFile ( self % iter , self % t , self % bufferLine )
            end do
#endif

#if defined(NAVIERSTOKES) || defined(INCNS)
            do i = 1 , self % no_of_surfaceMonitors
               call self % surfaceMonitors(i) % WriteToFile ( self % iter , self % t , self % bufferLine )
            end do
!
!              Write statistics
!              ----------------
            if ( self % bufferLine .eq. 0 ) then
               i = 1
            else
               i = self % bufferLine
            end if
            call self % stats % WriteFile(mesh, self % iter(i), self % t(i), self % solution_file)
#endif
!
!           Reset buffer
!           ------------
            self % bufferLine = 0

         else
!
!           The monitors are exported just if the buffer is full
!           ----------------------------------------------------
            if ( self % bufferLine .eq. BUFFER_SIZE ) then

               call self % residuals % WriteToFile ( self % iter , self % t, self % TotalSimuTime, self % SolverSimuTime, BUFFER_SIZE )

               do i = 1 , self % no_of_volumeMonitors
                  call self % volumeMonitors(i) % WriteToFile ( self % iter , self % t , self % bufferLine )
               end do

               do i = 1 , self % no_of_loadBalancingMonitors
                  call self % loadBalancingMonitors(i) % WriteToFile ( self % iter , self % t , self % bufferLine )
               end do

#ifdef FLOW
               do i = 1 , self % no_of_probes - self % no_of_fileProbes
                  call self % probes(i) % WriteToFile ( self % iter , self % t , self % bufferLine )
               end do
#endif

#if defined(NAVIERSTOKES) || defined(INCNS)
               do i = 1 , self % no_of_surfaceMonitors
                  call self % surfaceMonitors(i) % WriteToFile ( self % iter , self % t , self % bufferLine )
               end do
#endif
!
!              Reset buffer
!              ------------
               self % bufferLine = 0

            end if
         end if

      end subroutine Monitor_WriteToFile

      subroutine Monitor_WritePostProcessingSummary ( self )
!
!        ********************************************************************
!              Prints a unified post-processing summary section at startup,
!           covering probes, file probes, volume monitors, and surface
!           monitors.
!        ********************************************************************
!
         use Headers
         use MPI_Process_Info
         implicit none
         class(Monitor_t)        :: self
!
!        ---------------
!        Local variables
!        ---------------
!
         integer :: i, j
         integer :: no_of_stdProbes

         if ( .not. MPI_Process % isRoot ) return

         no_of_stdProbes = self % no_of_probes - self % no_of_fileProbes

         write(STD_OUT,'(/)')
         call Section_Header("Post-processing")
         write(STD_OUT,'(/)')

#ifdef FLOW
!
!        -- Standard probes (point probes defined in control file) -----------
!
         if ( no_of_stdProbes .gt. 0 ) then
            call SubSection_Header("Probes")
            write(STD_OUT,'(/)')
            write(STD_OUT,'(30X,A,A28,I0)') "->" , "Number of probes: " , no_of_stdProbes
            write(STD_OUT,'(/)')
         end if
!
!        -- File probes (bulk probe file) ------------------------------------
!
         if ( self % no_of_fileProbes .gt. 0 ) then
            call SubSection_Header("File probes")
            write(STD_OUT,'(/)')
            write(STD_OUT,'(30X,A,A28,A)')   "->" , "File: " , trim(self % probesFileName)
            write(STD_OUT,'(30X,A,A28,I0)')  "->" , "Number of probes: " , self % no_of_fileProbes
            write(STD_OUT,'(30X,A,A28)',advance="no") "->" , "Variables: "
            do j = 1 , size(self % probesVariables)
               write(STD_OUT,'(A)',advance="no") trim(self % probesVariables(j)) // " "
            end do
            write(STD_OUT,*)
            if ( self % probeFileSaveTimestep .gt. 0.0_RP ) then
               write(STD_OUT,'(30X,A,A28,ES14.6)') "->" , "Save timestep: " , self % probeFileSaveTimestep
            else
               write(STD_OUT,'(30X,A,A28,A)') "->" , "Save timestep: " , "every step"
            end if
            write(STD_OUT,'(30X,A,A28,A)') "->" , "Output format: " , trim(self % probeFileOutputFormat)
            write(STD_OUT,'(/)')
         end if
#endif

!
!        -- Volume monitors -------------------------------------------------
!
         if ( self % no_of_volumeMonitors .gt. 0 ) then
            call SubSection_Header("Volume monitors")
            write(STD_OUT,'(/)')
            do i = 1 , self % no_of_volumeMonitors
               write(STD_OUT,'(30X,A,I0,A,A,A,A)') "Monitor ", i, ": ", &
                  trim(self % volumeMonitors(i) % monitorName), " - ", trim(self % volumeMonitors(i) % variable)
            end do
            write(STD_OUT,'(/)')
         end if

#if defined(NAVIERSTOKES) || defined(INCNS)
!
!        -- Surface monitors ------------------------------------------------
!
         if ( self % no_of_surfaceMonitors .gt. 0 ) then
            call SubSection_Header("Surface monitors")
            write(STD_OUT,'(/)')
            do i = 1 , self % no_of_surfaceMonitors
               write(STD_OUT,'(30X,A,I0,A,A,A,A)') "Monitor ", i, ": ", &
                  trim(self % surfaceMonitors(i) % monitorName), " - ", trim(self % surfaceMonitors(i) % variable)
            end do
            write(STD_OUT,'(/)')
         end if
#endif

      end subroutine Monitor_WritePostProcessingSummary

      subroutine Monitor_Destruct (self)
         implicit none
         class(Monitor_t)        :: self
         integer                 :: i

         deallocate (self % iter)
         deallocate (self % t)
         deallocate (self % TotalSimuTime)
         deallocate (self % SolverSimuTime)
         
         call self % residuals % destruct
         
         call self % volumeMonitors % destruct
         safedeallocate(self % volumeMonitors)

         call self % loadBalancingMonitors % destruct
         safedeallocate(self % loadBalancingMonitors)
         
#ifdef FLOW
         if ( allocated(self % probes) ) then
            do i = 1, size(self % probes)
               if ( self % probes(i) % fileUnit >= 0 ) close( self % probes(i) % fileUnit )
            end do
         end if
         call self % probes % destruct
         safedeallocate (self % probes)

         if ( allocated(self % fp_fileUnit) ) then
            do i = 1, size(self % fp_fileUnit)
               if ( self % fp_fileUnit(i) .ge. 0 ) close( self % fp_fileUnit(i) )
            end do
         end if
         safedeallocate( self % fp_x )
         safedeallocate( self % fp_fileUnit )
         safedeallocate( self % fp_active )
#ifdef HAS_HDF5
         call Monitor_CloseHDF5FP()
#endif
#endif
#ifdef _OPENACC
         if ( allocated(self % fp_eID) ) then
            !$acc exit data delete(self % fp_eID, self % fp_ownsProbe)
            !$acc exit data delete(self % fp_lxi, self % fp_leta, self % fp_lzeta)
            !$acc exit data delete(self % fp_varCodes, self % fp_values_gpu)
            deallocate(self % fp_eID, self % fp_ownsProbe)
            deallocate(self % fp_lxi, self % fp_leta, self % fp_lzeta)
            deallocate(self % fp_varCodes, self % fp_values_gpu)
         end if
#else
         safedeallocate( self % fp_ownedIdx )
         safedeallocate( self % fp_buf )
         safedeallocate( self % fp_cpu_eID )
         safedeallocate( self % fp_cpu_Nx )
         safedeallocate( self % fp_cpu_Ny )
         safedeallocate( self % fp_cpu_Nz )
         safedeallocate( self % fp_cpu_lxi )
         safedeallocate( self % fp_cpu_leta )
         safedeallocate( self % fp_cpu_lzeta )
         safedeallocate( self % fp_cpu_varCodes )
         safedeallocate( self % fp_owned_buf )
         safedeallocate( self % fp_gatherv_counts )
         safedeallocate( self % fp_gatherv_displs )
         safedeallocate( self % fp_gatherv_perm )
#endif
         
#if defined(NAVIERSTOKES) || defined(INCNS)
         call self % surfaceMonitors % destruct
         safedeallocate (self % surfaceMonitors)
         
         !call self % stats % destruct
#endif         
      end subroutine
      
      impure elemental subroutine Monitor_Assign ( to, from )
         use MPI_Process_Info
         implicit none
         !-arguments--------------------------------------
         class(Monitor_t), intent(inout)  :: to
         type(Monitor_t) , intent(in)     :: from
         !-local-variables--------------------------------
         !------------------------------------------------

         to % solution_file               = from % solution_file
         to % probes_solution_file        = from % probes_solution_file
         to % no_of_probes                = from % no_of_probes
         to % no_of_surfaceMonitors       = from % no_of_surfaceMonitors
         to % no_of_volumeMonitors        = from % no_of_volumeMonitors
         to % no_of_loadBalancingMonitors = from % no_of_loadBalancingMonitors
         to % no_of_fileProbes            = from % no_of_fileProbes
         to % probesFileName              = from % probesFileName
         to % probeFileSaveTimestep       = from % probeFileSaveTimestep
         to % probeFileOutputFormat       = from % probeFileOutputFormat
         to % fp_lastSavedTime            = from % fp_lastSavedTime
         if ( allocated(from % probesVariables) ) then
            safedeallocate ( to % probesVariables )
            allocate ( to % probesVariables ( size(from % probesVariables) ) )
            to % probesVariables = from % probesVariables
         end if
         to % bufferLine                  = from % bufferLine
         
         safedeallocate ( to % iter )
         if ( allocated(from % iter) ) then
            allocate ( to % iter ( size(from % iter) ) )
            to % iter = from % iter
         end if
         
         to % dt_restriction        = from % dt_restriction
         to % write_dt_restriction  = from % write_dt_restriction
         
         safedeallocate (to % t)
         allocate (to % t (size (from % t) ) ) 
         to % t = from % t
         
         safedeallocate ( to % TotalSimuTime )
         allocate ( to % TotalSimuTime ( size(from % TotalSimuTime) ) )
         to % TotalSimuTime = from % TotalSimuTime
         
         safedeallocate ( to % SolverSimuTime )
         allocate ( to % SolverSimuTime ( size(from % SolverSimuTime) ) )
         to % SolverSimuTime = from % SolverSimuTime
         
         to % residuals = from % residuals
         
         safedeallocate ( to % volumeMonitors )
         allocate ( to % volumeMonitors ( size(from % volumeMonitors) ) )
         to % volumeMonitors = from % volumeMonitors

         safedeallocate ( to % loadBalancingMonitors )
         allocate ( to % loadBalancingMonitors ( size(from % loadBalancingMonitors) ) )
         to % loadBalancingMonitors = from % loadBalancingMonitors
      
#ifdef FLOW
         safedeallocate ( to % probes )
         allocate ( to % probes ( size(from % probes) ) )
         to % probes = from % probes
#endif

!
!        File-probe SoA buffers: fp_nOwned/fp_ownedIdx/fp_buf are used by
!        the CPU path (Monitor_ComputeFileProbesCPU, Monitor_UpdateFileProbes's
!        Allreduce); fp_cpu_*/fp_eID+friends are each build's own per-probe
!        arrays. None of this was copied before, so any 'to' object that
!        goes on to call the file-probe update path (e.g. a DGSem copy used
!        for truncation-error/load-balancing estimation) would hit
!        unallocated arrays there - a real crash, not just a stale-data bug.
!        --------------------------------------------------------------------
         if ( allocated(from % fp_x) ) then
            safedeallocate(to % fp_x) ; allocate(to % fp_x(size(from%fp_x,1),size(from%fp_x,2))) ; to % fp_x = from % fp_x
         end if
         if ( allocated(from % fp_fileUnit) ) then
            safedeallocate(to % fp_fileUnit) ; allocate(to % fp_fileUnit(size(from % fp_fileUnit))) ; to % fp_fileUnit = from % fp_fileUnit
         end if
         if ( allocated(from % fp_active) ) then
            safedeallocate(to % fp_active) ; allocate(to % fp_active(size(from % fp_active))) ; to % fp_active = from % fp_active
         end if

         to % fp_nOwned = from % fp_nOwned
         if ( allocated(from % fp_ownedIdx) ) then
            safedeallocate(to % fp_ownedIdx) ; allocate(to % fp_ownedIdx(size(from % fp_ownedIdx))) ; to % fp_ownedIdx = from % fp_ownedIdx
         end if
         if ( allocated(from % fp_buf) ) then
            safedeallocate(to % fp_buf) ; allocate(to % fp_buf(size(from % fp_buf))) ; to % fp_buf = from % fp_buf
         end if

#ifdef _OPENACC
         to % fp_Nmax = from % fp_Nmax
         if ( allocated(from % fp_eID) ) then
            safedeallocate(to % fp_eID)       ; allocate(to % fp_eID(size(from % fp_eID)))             ; to % fp_eID       = from % fp_eID
            safedeallocate(to % fp_ownsProbe) ; allocate(to % fp_ownsProbe(size(from % fp_ownsProbe))) ; to % fp_ownsProbe = from % fp_ownsProbe
            safedeallocate(to % fp_lxi)       ; allocate(to % fp_lxi(size(from%fp_lxi,1),size(from%fp_lxi,2)))        ; to % fp_lxi       = from % fp_lxi
            safedeallocate(to % fp_leta)      ; allocate(to % fp_leta(size(from%fp_leta,1),size(from%fp_leta,2)))     ; to % fp_leta      = from % fp_leta
            safedeallocate(to % fp_lzeta)     ; allocate(to % fp_lzeta(size(from%fp_lzeta,1),size(from%fp_lzeta,2))) ; to % fp_lzeta     = from % fp_lzeta
            safedeallocate(to % fp_varCodes)  ; allocate(to % fp_varCodes(size(from % fp_varCodes)))   ; to % fp_varCodes  = from % fp_varCodes
            safedeallocate(to % fp_values_gpu); allocate(to % fp_values_gpu(size(from%fp_values_gpu,1),size(from%fp_values_gpu,2))) ; to % fp_values_gpu = from % fp_values_gpu
         end if
#else
         if ( allocated(from % fp_owned_buf) ) then
            safedeallocate(to % fp_owned_buf) ; allocate(to % fp_owned_buf(size(from % fp_owned_buf))) ; to % fp_owned_buf = from % fp_owned_buf
         end if
         if ( allocated(from % fp_gatherv_counts) ) then
            safedeallocate(to % fp_gatherv_counts) ; allocate(to % fp_gatherv_counts(size(from % fp_gatherv_counts))) ; to % fp_gatherv_counts = from % fp_gatherv_counts
            safedeallocate(to % fp_gatherv_displs) ; allocate(to % fp_gatherv_displs(size(from % fp_gatherv_displs))) ; to % fp_gatherv_displs = from % fp_gatherv_displs
         end if
         if ( allocated(from % fp_gatherv_perm) ) then
            safedeallocate(to % fp_gatherv_perm) ; allocate(to % fp_gatherv_perm(size(from % fp_gatherv_perm))) ; to % fp_gatherv_perm = from % fp_gatherv_perm
         end if

         if ( allocated(from % fp_cpu_eID) ) then
            safedeallocate(to % fp_cpu_eID)      ; allocate(to % fp_cpu_eID(size(from % fp_cpu_eID)))           ; to % fp_cpu_eID      = from % fp_cpu_eID
            safedeallocate(to % fp_cpu_Nx)       ; allocate(to % fp_cpu_Nx(size(from % fp_cpu_Nx)))             ; to % fp_cpu_Nx       = from % fp_cpu_Nx
            safedeallocate(to % fp_cpu_Ny)       ; allocate(to % fp_cpu_Ny(size(from % fp_cpu_Ny)))             ; to % fp_cpu_Ny       = from % fp_cpu_Ny
            safedeallocate(to % fp_cpu_Nz)       ; allocate(to % fp_cpu_Nz(size(from % fp_cpu_Nz)))             ; to % fp_cpu_Nz       = from % fp_cpu_Nz
            safedeallocate(to % fp_cpu_lxi)      ; allocate(to % fp_cpu_lxi(size(from%fp_cpu_lxi,1),size(from%fp_cpu_lxi,2)))       ; to % fp_cpu_lxi      = from % fp_cpu_lxi
            safedeallocate(to % fp_cpu_leta)     ; allocate(to % fp_cpu_leta(size(from%fp_cpu_leta,1),size(from%fp_cpu_leta,2)))     ; to % fp_cpu_leta     = from % fp_cpu_leta
            safedeallocate(to % fp_cpu_lzeta)    ; allocate(to % fp_cpu_lzeta(size(from%fp_cpu_lzeta,1),size(from%fp_cpu_lzeta,2)))   ; to % fp_cpu_lzeta    = from % fp_cpu_lzeta
            safedeallocate(to % fp_cpu_varCodes) ; allocate(to % fp_cpu_varCodes(size(from % fp_cpu_varCodes))) ; to % fp_cpu_varCodes = from % fp_cpu_varCodes
         end if
#endif

#if defined(NAVIERSTOKES) || defined(INCNS)
         safedeallocate ( to % surfaceMonitors )
         allocate ( to % surfaceMonitors ( size(from % surfaceMonitors) ) )
         to % surfaceMonitors = from % surfaceMonitors
         
         to % stats = from % stats
#endif
         
      end subroutine Monitor_Assign
      
!
!//////////////////////////////////////////////////////////////////////////////
!
!        Auxiliars
!
!//////////////////////////////////////////////////////////////////////////////
!
   subroutine getNoOfMonitors(no_of_probes, no_of_surfaceMonitors, no_of_volumeMonitors, no_of_loadBalancingMonitors)
      use ParamfileRegions
      implicit none
      integer, intent(out)    :: no_of_probes
      integer, intent(out)    :: no_of_surfaceMonitors
      integer, intent(out)    :: no_of_volumeMonitors
      integer, intent(out)    :: no_of_loadBalancingMonitors
!
!     ---------------
!     Local variables
!     ---------------
!
      character(len=LINE_LENGTH) :: case_name, line
      integer                    :: fID
      integer                    :: io
!
!     Initialize
!     ----------
      no_of_probes = 0
      no_of_surfaceMonitors = 0
      no_of_volumeMonitors = 0
      no_of_loadBalancingMonitors = 0
!
!     Get case file name
!     ------------------
      call get_command_argument(1, case_name)

!
!     Open case file
!     --------------
      open ( newunit = fID , file = case_name , status = "old" , action = "read" )

!
!     Read the whole file to find monitors
!     ------------------------------------
readloop:do 
         read ( fID , '(A)' , iostat = io ) line

         if ( io .lt. 0 ) then
!
!           End of file
!           -----------
            line = ""
            exit readloop

         elseif ( io .gt. 0 ) then
!
!           Error
!           -----
            errorMessage(STD_OUT)
            error stop "Stopped."

         else
!
!           Succeeded
!           ---------
            line = getSquashedLine( line )

            if ( index ( line , '#defineprobefile' ) .gt. 0 ) then
!
!              The probe-file block is not an individual probe definition
!              -----------------------------------------------------------

            elseif ( index ( line , '#defineprobe' ) .gt. 0 ) then
               no_of_probes = no_of_probes + 1

            elseif ( index ( line , '#definesurfacemonitor' ) .gt. 0 ) then
               no_of_surfaceMonitors = no_of_surfaceMonitors + 1 

            elseif ( index ( line , '#definevolumemonitor' ) .gt. 0 ) then
               no_of_volumeMonitors = no_of_volumeMonitors + 1 

            elseif ( index ( line , '#defineloadbalancingmonitor' ) .gt. 0 ) then
               no_of_loadBalancingMonitors = no_of_loadBalancingMonitors + 1

            end if
            
         end if

      end do readloop
!
!     Close case file
!     ---------------
      close(fID)                             

end subroutine getNoOfMonitors

!
!///////////////////////////////////////////////////////////////////////////////////
!
!     Probes-from-file auxiliary routines
!
!     The probes file is a plain text file. Blank lines and lines starting
!     with "#" are ignored. Every remaining line holds the coordinates of
!     one probe:
!
!        x  y  z
!
!     The list of variables to sample (shared by every probe in the file)
!     is given through the "probes file variables" control-file keyword.
!
!     A separate output file "<solution_file>.probe_<N>.probe" is created for
!     each probe, with one column per variable, and one row written per
!     saved time step.
!
!///////////////////////////////////////////////////////////////////////////////////
!
   subroutine readProbesFileBlock(fileName, variablesLine, no_of_probesVariables, saveTimestep, outputFormat)
!
!     ******************************************************************
!        Reads the "#define probe file ... #end" block from the case
!     file, if present:
!
!        #define probe file
!           file                = Probe.dat
!           variables           = u
!           probe save timestep = 1.0E-3
!           output format       = HDF5
!        #end
!
!     Note: this is a hand-rolled parser (rather than readValueInRegion)
!     because readValueInRegion lower-cases every line it reads, which
!     would corrupt a case-sensitive file path.
!     ******************************************************************
!
      use ParamfileRegions, only: getSquashedLine
      implicit none
      character(len=LINE_LENGTH), intent(out) :: fileName
      character(len=LINE_LENGTH), intent(out) :: variablesLine
      integer,                    intent(out) :: no_of_probesVariables
      real(kind=RP),              intent(out) :: saveTimestep
      character(len=8),           intent(out) :: outputFormat
!
!     ---------------
!     Local variables
!     ---------------
!
      character(len=LINE_LENGTH) :: paramFile
      character(len=LINE_LENGTH) :: line, squashed, valStr
      integer                    :: fID, io, position

      logical                    :: inside

      fileName              = ""
      variablesLine         = ""
      no_of_probesVariables = 0
      saveTimestep          = 0.0_RP
      outputFormat          = "ASCII"
      inside                = .false.

      call get_command_argument(1, paramFile)

      open ( newunit = fID , file = trim(paramFile) , status = "old" , action = "read" )

      do
         read ( fID , '(A)' , iostat = io ) line
         if ( io .ne. 0 ) exit

         squashed = getSquashedLine(line)

         if ( squashed .eq. getSquashedLine("#define probe file") ) then
            inside = .true.
            cycle
         elseif ( squashed .eq. getSquashedLine("#end") ) then
            inside = .false.
            cycle
         end if

         if ( .not. inside ) cycle
!
!        Strip a trailing comment, keeping the value's original case
!        -------------------------------------------------------------
         position = index(line , '!')
         if ( position .gt. 0 ) line = line(1:position-1)

         position = max( index(line,'='), index(line,':') )
         if ( position .eq. 0 ) cycle

         if ( getSquashedLine(line(1:position-1)) .eq. getSquashedLine("file") ) then
            fileName = adjustl( removeQuotes( line(position+1:) ) )
         elseif ( getSquashedLine(line(1:position-1)) .eq. getSquashedLine("variables") ) then
            variablesLine = adjustl( removeQuotes( line(position+1:) ) )
         elseif ( getSquashedLine(line(1:position-1)) .eq. getSquashedLine("probe save timestep") ) then
            valStr = adjustl( removeQuotes( line(position+1:) ) )
            read( valStr , * ) saveTimestep
         elseif ( getSquashedLine(line(1:position-1)) .eq. getSquashedLine("output format") ) then
            valStr = adjustl( removeQuotes( line(position+1:) ) )
            valStr = adjustl( getSquashedLine(valStr) )
            if ( index(valStr, "HDF5") .gt. 0 ) then
               outputFormat = "HDF5"
            else
               outputFormat = "ASCII"
            end if
         end if
      end do

      close(fID)

   end subroutine readProbesFileBlock

   function removeQuotes(str) result(res)
      implicit none
      character(len=*), intent(in) :: str
      character(len=LINE_LENGTH)   :: res
      character(len=LINE_LENGTH)   :: auxstr
      integer                      :: i, j

      auxstr = trim(adjustl(str))
      res    = ""
      j      = 0

      do i = 1 , len_trim(auxstr)
         if ( auxstr(i:i) .eq. '"' .or. auxstr(i:i) .eq. "'" ) cycle
         j = j + 1
         res(j:j) = auxstr(i:i)
      end do

   end function removeQuotes

   subroutine countProbesInFile(fileName, n)
      implicit none
      character(len=*), intent(in)  :: fileName
      integer,          intent(out) :: n
!
!     ---------------
!     Local variables
!     ---------------
!
      integer                    :: fID, io
      character(len=LINE_LENGTH) :: line

      n = 0
      open ( newunit = fID , file = fileName , status = "old" , action = "read" )

      do
         read ( fID , '(A)' , iostat = io ) line
         if ( io .ne. 0 ) exit

         line = adjustl(line)
         if ( len_trim(line) .eq. 0 ) cycle
         if ( line(1:1) .eq. '#' )    cycle

         n = n + 1
      end do

      close(fID)

   end subroutine countProbesInFile

   subroutine splitIntoTokens(line, tokens, n)
      implicit none
      character(len=*),                              intent(in)  :: line
      character(len=STR_LEN_MONITORS), allocatable,   intent(out) :: tokens(:)
      integer,                                        intent(out) :: n
!
!     ---------------
!     Local variables
!     ---------------
!
      character(len=LINE_LENGTH)      :: auxline
      character(len=STR_LEN_MONITORS) :: buffer(64)
      integer                         :: pos

      auxline = adjustl(line)
      n = 0

      do while ( len_trim(auxline) .gt. 0 )
         if ( n .ge. size(buffer) ) then
            write(*,'(A)') "ERROR: splitIntoTokens: too many tokens (max 64). Truncating."
            exit
         end if
         pos = index(trim(auxline), " ")
         n = n + 1

         if ( pos .eq. 0 ) then
            buffer(n) = trim(auxline)
            auxline   = ""
         else
            buffer(n) = auxline(1:pos-1)
            auxline   = adjustl(auxline(pos+1:))
         end if
      end do

      allocate( tokens(n) )
      tokens(1:n) = buffer(1:n)

   end subroutine splitIntoTokens

#ifdef FLOW
   subroutine Monitor_FlushFileProbesNow(self, t_now, iter_now)
!
!     Writes the current (single-slot) file-probe values to disk,
!     applying the probeFileSaveTimestep filter.  Called every timestep
!     from Monitor_UpdateValues after Monitor_UpdateFileProbes.
!
      use MPI_Process_Info
      implicit none
      class(Monitor_t), intent(inout) :: self
      real(kind=RP),    intent(in)    :: t_now
      integer,          intent(in)    :: iter_now
!
!     ---------------
!     Local variables
!     ---------------
!
      integer        :: iter_arr(1)
      real(kind=RP)  :: t_arr(1)
      logical        :: do_write

      iter_arr(1) = iter_now
      t_arr(1)   = t_now

      ! Monitor-level timestep filter (shared by both output formats below):
      ! skip the write entirely outside the save-timestep window.
      do_write = .true.
      if ( self % probeFileSaveTimestep .gt. 0.0_RP ) then
         if ( t_now .lt. self % fp_lastSavedTime + self % probeFileSaveTimestep ) do_write = .false.
      end if

      if ( do_write ) then
         self % fp_lastSavedTime = t_now
#ifdef HAS_HDF5
         if ( trim(self % probeFileOutputFormat) .eq. "HDF5" ) then
            call Monitor_WriteFileProbesHDF5( self, iter_arr, t_arr, 1 )
         else
#endif
            call Monitor_WriteFileProbesASCII( self, iter_arr, t_arr, 1 )
#ifdef HAS_HDF5
         end if
#endif
      end if

   end subroutine Monitor_FlushFileProbesNow

   subroutine Monitor_WriteFileProbesASCII(self, iter, t, no_of_lines)
!
!     Appends one buffer of file-probe data to each probe's own ASCII
!     file. Root-only, like the old Probe_t % WriteToFile: fp_buf/
!     fp_values_gpu already hold the globally-reduced values on every
!     rank (via Monitor_UpdateFileProbes's Allreduce), so root alone can
!     write every probe regardless of which rank actually owns it. The
!     file is opened once (header written) and kept open across calls
!     via fp_fileUnit - no Probe_t needed to track this per probe.
!     Probes that were never found anywhere (fp_active=.false.) get no
!     file at all, matching the old behavior.
!     -------------------------------------------------------------------
      use MPI_Process_Info
      implicit none
      class(Monitor_t), intent(inout) :: self
      integer,          intent(in)    :: iter(:)
      real(kind=RP),    intent(in)    :: t(:)
      integer,          intent(in)    :: no_of_lines
!
!     ---------------
!     Local variables
!     ---------------
!
      integer                         :: i, v, l, fID, nfp, nv, offset
      character(len=LINE_LENGTH)      :: fname
      character(len=STR_LEN_MONITORS) :: pname

      if ( .not. MPI_Process % isRoot ) return

      nfp    = self % no_of_fileProbes
      nv     = size(self % probesVariables)
      offset = self % no_of_probes - nfp

      do i = 1, nfp
         if ( .not. self % fp_active(i) ) cycle

         if ( self % fp_fileUnit(i) .lt. 0 ) then
            write(pname,'(A,I0)') "probe_", offset + i
            write(fname,'(A,A,A,A)') trim(self % probes_solution_file), "." , trim(pname) , ".probe"
            open( newunit = fID , file = trim(fname) , status = "unknown" , action = "write" )

            write( fID , '(A20,A  )') "Monitor name:      ", trim(pname)
            write( fID , '(A25,ES24.10,2(4X,ES24.10))') "x, y, z coordinates: ", self % fp_x(1,i), self % fp_x(2,i), self % fp_x(3,i)
            write( fID , * )
            write( fID , '(A10,2X,A24)' , advance = "no") "Iteration" , "Time"
            do v = 1 , nv
               write( fID , '(2X,A24)' , advance = "no") trim(self % probesVariables(v))
            end do
            write( fID , * )

            self % fp_fileUnit(i) = fID
         end if

         fID = self % fp_fileUnit(i)
         do l = 1, no_of_lines
            write( fID , '(I10,2X,ES24.16)' , advance = "no" ) iter(l) , t(l)
            do v = 1 , nv
#ifdef _OPENACC
               write( fID , '(2X,ES24.16)' , advance = "no" ) self % fp_values_gpu(v,i)
#else
               write( fID , '(2X,ES24.16)' , advance = "no" ) self % fp_buf((i-1)*nv + v)
#endif
            end do
            write( fID , * )
         end do
      end do

   end subroutine Monitor_WriteFileProbesASCII

   subroutine Monitor_UpdateFileProbes(self, mesh, bufferPos)
      use MPI_Process_Info
#ifdef _HAS_MPI_
      use mpi
#endif
      implicit none
      class(Monitor_t), intent(inout) :: self
      class(HexMesh),   intent(in)    :: mesh
      integer,          intent(in)    :: bufferPos
!
!     ---------------
!     Local variables
!     ---------------
!
      integer        :: i, v, j, nfp, nv, ierr
#ifdef _OPENACC
      integer        :: Nm
#endif
!     --- timing ---
      integer(kind=8) :: t0, t1, t2, t3, rate
!     --- gatherv ---
      real(kind=RP), allocatable :: gather_tmp(:)
      integer :: k, jbuf_g

      nfp       = self % no_of_fileProbes
      nv        = size(self % probesVariables)

#ifdef _OPENACC
!
!     GPU path: parallel evaluation of all file-probes on device.
!     Lagrange weights and element IDs are pre-loaded in SoA arrays by InitializeProbesFromFile.
!     Non-owning ranks skip computation (fp_ownsProbe=.false.) and contribute 0 to MPI_Allreduce.
!     The kernel is in Monitor_FileProbeKernel so arrays arrive as dummy arguments — this avoids
!     NVFORTRAN accessing them through the host-side 'self' struct pointer on the GPU.
!
         Nm = self % fp_Nmax
         call Monitor_FileProbeKernel(self % fp_eID, self % fp_ownsProbe, &
                                      self % fp_lxi, self % fp_leta, self % fp_lzeta, &
                                      self % fp_varCodes, self % fp_values_gpu, &
                                      mesh, nfp, nv, Nm)
         !$acc update host(self % fp_values_gpu)
#ifdef _HAS_MPI_
         if ( MPI_Process % doMPIAction ) then
            call MPI_Allreduce(MPI_IN_PLACE, self % fp_values_gpu, nfp * nv, &
                               MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
         end if
#endif

#else
!
!     CPU path: compact SoA compute — reads eID/lxi/leta/lzeta from contiguous
!     fp_cpu_* arrays built once at init time (InitializeProbesFromFile),
!     writes directly into fp_buf.
!
      call system_clock(t0, rate)
!     (no zero-fill needed: ComputeFileProbesCPU writes compact fp_owned_buf)
      call system_clock(t1)
      call Monitor_ComputeFileProbesCPU(self, mesh, nv)
      call system_clock(t2)

#ifdef _HAS_MPI_
      if ( MPI_Process % doMPIAction ) then
!        MPI_Gatherv: each rank sends its compact fp_owned_buf to root.
!        Root unpacks using fp_gatherv_perm into fp_buf.
!        ALL ranks must call MPI_Gatherv (collective); only root uses the
!        recv-side arguments (recvcounts/recvdispls), so non-root passes
!        dummy arrays.  The previous guard "allocated(fp_gatherv_counts)"
!        was TRUE only on root, causing non-root to skip the collective and
!        leaving root blocked in MPI_Gatherv forever (deadlock).
         allocate( gather_tmp(merge(nfp * nv, 1, MPI_Process % isRoot)) )
         if ( MPI_Process % isRoot ) then
            call MPI_Gatherv(self % fp_owned_buf, self % fp_nOwned * nv, MPI_DOUBLE_PRECISION, &
                             gather_tmp, self % fp_gatherv_counts, self % fp_gatherv_displs, &
                             MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
         else
            block
               integer :: dummy_counts(1), dummy_displs(1)
               dummy_counts(1) = 0 ; dummy_displs(1) = 0
               call MPI_Gatherv(self % fp_owned_buf, self % fp_nOwned * nv, MPI_DOUBLE_PRECISION, &
                                gather_tmp, dummy_counts, dummy_displs, &
                                MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
            end block
         end if
         if ( MPI_Process % isRoot ) then
            do k = 1, nfp
               jbuf_g = (self % fp_gatherv_perm(k) - 1) * nv
               self % fp_buf(jbuf_g+1:jbuf_g+nv) = gather_tmp((k-1)*nv+1:k*nv)
            end do
         end if
         deallocate(gather_tmp)
      else
!        Single-rank fallback: copy compact buffer directly into fp_buf.
         do k = 1, self % fp_nOwned
            jbuf_g = (self % fp_ownedIdx(k) - 1) * nv
            self % fp_buf(jbuf_g+1:jbuf_g+nv) = self % fp_owned_buf((k-1)*nv+1:k*nv)
         end do
      end if
#else
      do k = 1, self % fp_nOwned
         jbuf_g = (self % fp_ownedIdx(k) - 1) * nv
         self % fp_buf(jbuf_g+1:jbuf_g+nv) = self % fp_owned_buf((k-1)*nv+1:k*nv)
      end do
#endif
      call system_clock(t3)

      if ( MPI_Process % isRoot ) then
         write(STD_OUT,'(/,30X,A)') "--- File-probe timing (rank 0) ---"
         write(STD_OUT,'(30X,A,ES12.4,A)') "  ComputeFileProbes  : ", real(t2-t1,8)/real(rate,8), " s"
         write(STD_OUT,'(30X,A,ES12.4,A)') "  MPI_Gatherv+unpack : ", real(t3-t2,8)/real(rate,8), " s"
         write(STD_OUT,'(30X,A,ES12.4,A)') "  Total              : ", real(t3-t0,8)/real(rate,8), " s"
      end if
#endif

   end subroutine Monitor_UpdateFileProbes

#ifndef _OPENACC
   subroutine Monitor_ComputeFileProbesCPU(self, mesh, nv)
!
!     Optimised SoA compute loop for file-probes on CPU.
!     Loop order: probe → node (kk,jj,ii) → variable.
!     Q(:,ii,jj,kk) is loaded once per node; all nv variables are accumulated
!     in a single pass, eliminating nv redundant cache-line fetches per node.
!     When probes are stored sorted by ascending eID (see InitializeProbesFromFile),
!     sequential probes touch the same or adjacent element Q arrays, further
!     improving cache reuse across probes.
!
      use Physics
      implicit none
      class(Monitor_t), intent(inout) :: self
      class(HexMesh),   intent(in)    :: mesh
      integer,          intent(in)    :: nv
!
!     Local variables
!
      integer        :: p, v, ii, jj, kk, eID, Nx, Ny, Nz, jbuf
      real(kind=RP)  :: w
      real(kind=RP)  :: acc(nv)
#ifdef NAVIERSTOKES
      real(kind=RP)  :: q_rho, q_rhou, q_rhov, q_rhow, q_rhoE, u2
#endif
#ifdef INCNS
      real(kind=RP)  :: q_rho, q_rhou, q_rhov, q_rhow, q_p
#endif
#ifdef MULTIPHASE
      real(kind=RP)  :: q_p, q_c, q_mu, q_cx, q_cy, q_cz
#endif

      do p = 1, self % fp_nOwned
         eID  = self % fp_cpu_eID(p)
         Nx   = self % fp_cpu_Nx(p)
         Ny   = self % fp_cpu_Ny(p)
         Nz   = self % fp_cpu_Nz(p)
         jbuf = (p - 1) * nv
         acc(1:nv) = 0.0_RP
         associate(Qe => mesh%elements(eID)%storage%Q)
         do kk = 0, Nz ; do jj = 0, Ny ; do ii = 0, Nx
            w = self%fp_cpu_lxi(ii,p)*self%fp_cpu_leta(jj,p)*self%fp_cpu_lzeta(kk,p)
!           Load conserved variables once for this node
#ifdef NAVIERSTOKES
            q_rho  = Qe(IRHO ,ii,jj,kk)
            q_rhou = Qe(IRHOU,ii,jj,kk)
            q_rhov = Qe(IRHOV,ii,jj,kk)
            q_rhow = Qe(IRHOW,ii,jj,kk)
            q_rhoE = Qe(IRHOE,ii,jj,kk)
#endif
#ifdef INCNS
            q_rho  = Qe(INSRHO ,ii,jj,kk)
            q_rhou = Qe(INSRHOU,ii,jj,kk)
            q_rhov = Qe(INSRHOV,ii,jj,kk)
            q_rhow = Qe(INSRHOW,ii,jj,kk)
            q_p    = Qe(INSP   ,ii,jj,kk)
#endif
#ifdef MULTIPHASE
            q_p  = Qe(IMP,ii,jj,kk)
            q_c  = Qe(IMC,ii,jj,kk)
            q_mu = mesh%elements(eID)%storage%mu(1,ii,jj,kk)
            q_cx = mesh%elements(eID)%storage%c_x(1,ii,jj,kk)
            q_cy = mesh%elements(eID)%storage%c_y(1,ii,jj,kk)
            q_cz = mesh%elements(eID)%storage%c_z(1,ii,jj,kk)
#endif
            do v = 1, nv
               select case (self % fp_cpu_varCodes(v))
#ifdef NAVIERSTOKES
               case(FPVAR_PRESSURE)
                  acc(v) = acc(v) + w * Pressure([q_rho,q_rhou,q_rhov,q_rhow,q_rhoE])
               case(FPVAR_VELOCITY)
                  acc(v) = acc(v) + w * sqrt(POW2(q_rhou)+POW2(q_rhov)+POW2(q_rhow)) / q_rho
               case(FPVAR_U)
                  acc(v) = acc(v) + w * q_rhou / q_rho
               case(FPVAR_V)
                  acc(v) = acc(v) + w * q_rhov / q_rho
               case(FPVAR_W)
                  acc(v) = acc(v) + w * q_rhow / q_rho
               case(FPVAR_MACH)
                  u2 = (POW2(q_rhou)+POW2(q_rhov)+POW2(q_rhow)) / POW2(q_rho)
                  acc(v) = acc(v) + w * sqrt( u2 / ( thermodynamics%gamma*(thermodynamics%gamma-1.0_RP) * &
                              (q_rhoE/q_rho - 0.5_RP*u2) ) )
               case(FPVAR_K)
                  acc(v) = acc(v) + w * 0.5_RP*(POW2(q_rhou)+POW2(q_rhov)+POW2(q_rhow)) / q_rho
               case(FPVAR_RHO)
                  acc(v) = acc(v) + w * q_rho
#endif
#ifdef INCNS
               case(FPVAR_PRESSURE)
                  acc(v) = acc(v) + w * q_p
               case(FPVAR_VELOCITY)
                  acc(v) = acc(v) + w * sqrt(POW2(q_rhou)+POW2(q_rhov)+POW2(q_rhow)) / q_rho
               case(FPVAR_U)
                  acc(v) = acc(v) + w * q_rhou / q_rho
               case(FPVAR_V)
                  acc(v) = acc(v) + w * q_rhov / q_rho
               case(FPVAR_W)
                  acc(v) = acc(v) + w * q_rhow / q_rho
               case(FPVAR_RHO)
                  acc(v) = acc(v) + w * q_rho
#endif
#ifdef MULTIPHASE
               case(FPVAR_STATICPRES)
                  acc(v) = acc(v) + w * ( q_p + q_c*q_mu &
                               - 12.0_RP*multiphase%sigma*multiphase%invEps*(POW2(q_c*(1.0_RP-q_c))) &
                               - 0.25_RP*3.0_RP*multiphase%sigma*multiphase%eps*(POW2(q_cx)+POW2(q_cy)+POW2(q_cz)) )
#endif
               end select
            end do
         end do ; end do ; end do
         end associate
         self % fp_owned_buf(jbuf+1:jbuf+nv) = acc(1:nv)
      end do

   end subroutine Monitor_ComputeFileProbesCPU
#endif

#ifdef _OPENACC
   subroutine Monitor_FileProbeKernel(l_eID, l_own, l_lxi, l_leta, l_lzeta, &
                                       l_varCodes, l_vals, mesh, nfp, nv, Nm)
      implicit none
      integer,          intent(in)    :: l_eID(:)
      logical,          intent(in)    :: l_own(:)
      real(kind=RP),    intent(in)    :: l_lxi(0:,1:)
      real(kind=RP),    intent(in)    :: l_leta(0:,1:)
      real(kind=RP),    intent(in)    :: l_lzeta(0:,1:)
      integer,          intent(in)    :: l_varCodes(:)
      real(kind=RP),    intent(inout) :: l_vals(:,:)
      class(HexMesh),   intent(in)    :: mesh
      integer,          intent(in)    :: nfp, nv, Nm
      integer        :: probe_idx, var_idx, ii, jj, kk, eID_loc
      real(kind=RP)  :: val, q_val

      !$acc parallel loop gang &
      !$acc& present(mesh, l_eID, l_own, l_lxi, l_leta, l_lzeta, l_varCodes, l_vals)
      do probe_idx = 1, nfp
         if ( .not. l_own(probe_idx) ) then
            do var_idx = 1, nv
               l_vals(var_idx, probe_idx) = 0.0_RP
            end do
         else
            eID_loc = l_eID(probe_idx)
            do var_idx = 1, nv
               val = 0.0_RP
               !$acc loop vector collapse(3) reduction(+:val)
               do kk = 0, Nm ; do jj = 0, Nm ; do ii = 0, Nm
                  select case (l_varCodes(var_idx))
#ifdef NAVIERSTOKES
                  case(FPVAR_PRESSURE)
                     q_val = Pressure(mesh % elements(eID_loc) % storage % Q(:,ii,jj,kk))
                  case(FPVAR_VELOCITY)
                     q_val = sqrt(POW2(mesh % elements(eID_loc) % storage % Q(IRHOU,ii,jj,kk)) + &
                                  POW2(mesh % elements(eID_loc) % storage % Q(IRHOV,ii,jj,kk)) + &
                                  POW2(mesh % elements(eID_loc) % storage % Q(IRHOW,ii,jj,kk))) / &
                             mesh % elements(eID_loc) % storage % Q(IRHO,ii,jj,kk)
                  case(FPVAR_U)
                     q_val = mesh % elements(eID_loc) % storage % Q(IRHOU,ii,jj,kk) / &
                             mesh % elements(eID_loc) % storage % Q(IRHO,ii,jj,kk)
                  case(FPVAR_V)
                     q_val = mesh % elements(eID_loc) % storage % Q(IRHOV,ii,jj,kk) / &
                             mesh % elements(eID_loc) % storage % Q(IRHO,ii,jj,kk)
                  case(FPVAR_W)
                     q_val = mesh % elements(eID_loc) % storage % Q(IRHOW,ii,jj,kk) / &
                             mesh % elements(eID_loc) % storage % Q(IRHO,ii,jj,kk)
                  case(FPVAR_MACH)
                     q_val = (POW2(mesh % elements(eID_loc) % storage % Q(IRHOU,ii,jj,kk)) + &
                              POW2(mesh % elements(eID_loc) % storage % Q(IRHOV,ii,jj,kk)) + &
                              POW2(mesh % elements(eID_loc) % storage % Q(IRHOW,ii,jj,kk))) / &
                             POW2(mesh % elements(eID_loc) % storage % Q(IRHO,ii,jj,kk))
                     q_val = sqrt( q_val / ( thermodynamics % gamma * (thermodynamics % gamma - 1.0_RP) * &
                             ( mesh % elements(eID_loc) % storage % Q(IRHOE,ii,jj,kk) / &
                               mesh % elements(eID_loc) % storage % Q(IRHO,ii,jj,kk) - 0.5_RP * q_val ) ) )
                  case(FPVAR_K)
                     q_val = 0.5_RP * (POW2(mesh % elements(eID_loc) % storage % Q(IRHOU,ii,jj,kk)) + &
                                       POW2(mesh % elements(eID_loc) % storage % Q(IRHOV,ii,jj,kk)) + &
                                       POW2(mesh % elements(eID_loc) % storage % Q(IRHOW,ii,jj,kk))) / &
                             mesh % elements(eID_loc) % storage % Q(IRHO,ii,jj,kk)
                  case(FPVAR_RHO)
                     q_val = mesh % elements(eID_loc) % storage % Q(IRHO,ii,jj,kk)
#endif
#ifdef INCNS
                  case(FPVAR_PRESSURE)
                     q_val = mesh % elements(eID_loc) % storage % Q(INSP,ii,jj,kk)
                  case(FPVAR_VELOCITY)
                     q_val = sqrt(POW2(mesh % elements(eID_loc) % storage % Q(INSRHOU,ii,jj,kk)) + &
                                  POW2(mesh % elements(eID_loc) % storage % Q(INSRHOV,ii,jj,kk)) + &
                                  POW2(mesh % elements(eID_loc) % storage % Q(INSRHOW,ii,jj,kk))) / &
                             mesh % elements(eID_loc) % storage % Q(INSRHO,ii,jj,kk)
                  case(FPVAR_U)
                     q_val = mesh % elements(eID_loc) % storage % Q(INSRHOU,ii,jj,kk) / &
                             mesh % elements(eID_loc) % storage % Q(INSRHO,ii,jj,kk)
                  case(FPVAR_V)
                     q_val = mesh % elements(eID_loc) % storage % Q(INSRHOV,ii,jj,kk) / &
                             mesh % elements(eID_loc) % storage % Q(INSRHO,ii,jj,kk)
                  case(FPVAR_W)
                     q_val = mesh % elements(eID_loc) % storage % Q(INSRHOW,ii,jj,kk) / &
                             mesh % elements(eID_loc) % storage % Q(INSRHO,ii,jj,kk)
                  case(FPVAR_RHO)
                     q_val = mesh % elements(eID_loc) % storage % Q(INSRHO,ii,jj,kk)
#endif
#ifdef MULTIPHASE
                  case(FPVAR_STATICPRES)
                     q_val = mesh % elements(eID_loc) % storage % Q(IMP,ii,jj,kk) &
                           + mesh % elements(eID_loc) % storage % Q(IMC,ii,jj,kk) &
                             * mesh % elements(eID_loc) % storage % mu(1,ii,jj,kk) &
                           - 12.0_RP*multiphase%sigma*multiphase%invEps &
                             * (mesh % elements(eID_loc) % storage % Q(IMC,ii,jj,kk) &
                                * (1.0_RP - mesh % elements(eID_loc) % storage % Q(IMC,ii,jj,kk)))**2 &
                           - 0.75_RP*multiphase%sigma*multiphase%eps &
                             * (POW2(mesh % elements(eID_loc) % storage % c_x(1,ii,jj,kk)) &
                              + POW2(mesh % elements(eID_loc) % storage % c_y(1,ii,jj,kk)) &
                              + POW2(mesh % elements(eID_loc) % storage % c_z(1,ii,jj,kk)))
#endif
#ifdef ACOUSTIC
                  case(FPVAR_PRESSURE)
                     q_val = mesh % elements(eID_loc) % storage % Q(ICAAP,ii,jj,kk)
                  case(FPVAR_DENSITY)
                     q_val = mesh % elements(eID_loc) % storage % Q(ICAARHO,ii,jj,kk)
                  case(FPVAR_U)
                     q_val = mesh % elements(eID_loc) % storage % Q(ICAAU,ii,jj,kk)
                  case(FPVAR_V)
                     q_val = mesh % elements(eID_loc) % storage % Q(ICAAV,ii,jj,kk)
                  case(FPVAR_W)
                     q_val = mesh % elements(eID_loc) % storage % Q(ICAAW,ii,jj,kk)
#endif
                  case default
                     q_val = 0.0_RP
                  end select
                  val = val + q_val * l_lxi(ii,probe_idx) &
                                    * l_leta(jj,probe_idx) &
                                    * l_lzeta(kk,probe_idx)
               end do ; end do ; end do
               l_vals(var_idx, probe_idx) = val
            end do
         end if
      end do
      !$acc end parallel loop

   end subroutine Monitor_FileProbeKernel
#endif

   subroutine InitializeProbesFromFile(fileName, Monitors, mesh, nfp)
      use MPI_Process_Info
#ifdef _HAS_MPI_
      use mpi
#endif
      implicit none
      character(len=*),   intent(in)    :: fileName
      class(Monitor_t),    intent(inout) :: Monitors
      class(HexMesh),     intent(inout) :: mesh
      integer,            intent(in)    :: nfp
!
!     ---------------
!     Local variables
!     ---------------
!
      integer                                      :: fID, io, nTok, i, nFound, ierr, prev_eID_local, nv
      character(len=LINE_LENGTH)                   :: line
      character(len=STR_LEN_MONITORS), allocatable  :: tokens(:)
      real(kind=RP),    allocatable :: allX(:,:)
      real(kind=RP),    allocatable :: xi_local(:,:)
      integer,          allocatable :: eID_local(:)
      logical,          allocatable :: foundLocal(:)
      integer,          allocatable :: ownerCandidate(:)
      integer,          allocatable :: globalOwner(:)
!     --- timing ---
      integer(kind=8) :: tinit_t0, tinit_t1, tinit_t2, tinit_t3, tinit_t4, tinit_t5, tinit_rate
!
!     nfp is the caller's own countProbesInFile result (Monitors_Construct),
!     not re-derived here: at O(1e3) MPI ranks every rank would otherwise
!     open and re-count this same file a second time over a shared/network
!     filesystem, and any inconsistency between the two counts (even from a
!     single rank, under heavy concurrent I/O) would size fp_x/fp_buf/etc.
!     below differently from what every other file-probe routine assumes
!     via Monitors % no_of_fileProbes - a silent out-of-bounds write.
!     --------------------------------------------------------------------
      if ( nfp .le. 0 ) return

      call system_clock(tinit_t0, tinit_rate)

      allocate( allX          (NDIM, nfp) )
      allocate( xi_local      (NDIM, nfp) )
      allocate( eID_local     (nfp)       )
      allocate( foundLocal    (nfp)       )
      allocate( ownerCandidate(nfp)       )
      allocate( globalOwner   (nfp)       )
!
!     Pass 1a: ONLY root reads the probes file and fills allX/nFound.
!     At O(1e3) MPI ranks, having every rank independently open and parse
!     the same (possibly huge, 1e6-line) file is not just wasteful - there
!     is no guarantee all of them count the same nFound under that much
!     concurrent I/O on a shared/network filesystem, and the MPI_Allreduce
!     a few lines below REQUIRES the same count on every rank or it hangs
!     forever (observed in practice on a 2240-rank MareNostrum 5 run, with
!     every rank's backtrace sitting inside that exact Allreduce). Reading
!     once on root and broadcasting removes the ambiguity by construction.
!     --------------------------------------------------------------------
      if ( MPI_Process % isRoot ) then
         nFound = 0
         open ( newunit = fID , file = fileName , status = "old" , action = "read" )

         do
            read ( fID , '(A)' , iostat = io ) line
            if ( io .ne. 0 ) exit

            line = adjustl(line)
            if ( len_trim(line) .eq. 0 ) cycle
            if ( line(1:1) .eq. '#' )    cycle

            call splitIntoTokens(line, tokens, nTok)

            if ( nTok .lt. 3 ) then
               deallocate(tokens)
               cycle
            end if

            nFound = nFound + 1
            read(tokens(1),*) allX(1,nFound)
            read(tokens(2),*) allX(2,nFound)
            read(tokens(3),*) allX(3,nFound)

            deallocate(tokens)
         end do

         close(fID)
      end if

#ifdef _HAS_MPI_
      if ( MPI_Process % doMPIAction ) then
         call MPI_Bcast(nFound, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
         call MPI_Bcast(allX, NDIM*nfp, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
      end if
#endif
!
      call system_clock(tinit_t1)
!
!     Pass 1b: every rank does its OWN local point search (no I/O, no MPI)
!     over the now-identical allX/nFound. The eID hint is cascaded from
!     one probe to the next only when the search ACTUALLY succeeded
!     locally, so it is always either -1 or a genuinely valid local
!     element on this rank.
!     --------------------------------------------------------------------
      if ( MPI_Process % isRoot ) &
         write(STD_OUT,'(30X,A,I0,A)') "-> Searching ", nFound, " probes in mesh..."

      prev_eID_local = -1
      do i = 1, nFound
         foundLocal(i) = mesh % FindPointWithCoords(allX(:,i), eID_local(i), &
                                                      xi_local(:,i), eID_hint=prev_eID_local)
         if ( foundLocal(i) ) prev_eID_local = eID_local(i)
         if ( MPI_Process % isRoot .and. mod(i, max(1,nFound/10)) .eq. 0 ) then
            write(STD_OUT,'(30X,A,I3,A)') "   ... ", (i*100)/nFound, "% done"
         end if
      end do
!
      call system_clock(tinit_t2)
!
!     Resolve, for EVERY file-probe at once, which rank owns it with a
!     single bulk MPI collective - not one mpi_allgather per probe (the
!     original Probe_LookInOtherPartitions path), which at O(1e6) probes
!     and O(1e3) ranks turned initialization into millions of tiny
!     sequential collectives instead of one.
!     --------------------------------------------------------------------
      do i = 1, nFound
         if ( foundLocal(i) ) then
            ownerCandidate(i) = MPI_Process % rank
         else
            ownerCandidate(i) = -1
         end if
      end do

#ifdef _HAS_MPI_
      if ( MPI_Process % doMPIAction ) then
         call MPI_Allreduce(ownerCandidate, globalOwner, nFound, MPI_INTEGER, MPI_MAX, MPI_COMM_WORLD, ierr)
      else
         globalOwner = ownerCandidate
      end if
#else
      globalOwner = ownerCandidate
#endif
!
      if ( MPI_Process % isRoot ) &
         write(STD_OUT,'(30X,A)') "-> Resolving probe ownership across ranks..."

      call system_clock(tinit_t3)
!
!     Pass 2: build Monitors' own SoA buffers directly from the LOCAL
!     eID_local/xi_local/globalOwner computed above - no Probe_t, no
!     per-probe allocation, no further communication. fp_x is the only
!     thing kept O(total probes) on every rank (just 3 reals/probe,
!     needed by root for the HDF5 /coordinates dataset and for ASCII
!     headers); everything else below is sized O(owned probes on this
!     rank).
!     --------------------------------------------------------------------
      nv = size(Monitors % probesVariables)

      allocate( Monitors % fp_x(NDIM, nFound) )
      Monitors % fp_x = allX(:, 1:nFound)

      allocate( Monitors % fp_active(nFound) )
      Monitors % fp_active = ( globalOwner(1:nFound) .ge. 0 )

      allocate( Monitors % fp_fileUnit(nFound) )
      Monitors % fp_fileUnit = -1

      allocate( Monitors % fp_buf(nFound * nv) )
      Monitors % fp_buf = 0.0_RP

#ifdef _OPENACC
      block
         integer :: ii, Nmax
         logical :: owns

         Nmax = 0
         do ii = 1, nFound
            if ( globalOwner(ii) .eq. MPI_Process % rank ) then
               Nmax = max(Nmax, maxval(mesh % elements(eID_local(ii)) % Nxyz))
            end if
         end do
#ifdef _HAS_MPI_
         if ( MPI_Process % doMPIAction ) then
            call MPI_Allreduce(MPI_IN_PLACE, Nmax, 1, MPI_INTEGER, MPI_MAX, MPI_COMM_WORLD, ierr)
         end if
#endif
         Monitors % fp_Nmax = Nmax

         allocate( Monitors % fp_varCodes(nv) )
         do ii = 1, nv
            select case (trim(Monitors % probesVariables(ii)))
            case("pressure")       ; Monitors % fp_varCodes(ii) = FPVAR_PRESSURE
            case("velocity")       ; Monitors % fp_varCodes(ii) = FPVAR_VELOCITY
            case("u")              ; Monitors % fp_varCodes(ii) = FPVAR_U
            case("v")              ; Monitors % fp_varCodes(ii) = FPVAR_V
            case("w")              ; Monitors % fp_varCodes(ii) = FPVAR_W
            case("mach")           ; Monitors % fp_varCodes(ii) = FPVAR_MACH
            case("k")              ; Monitors % fp_varCodes(ii) = FPVAR_K
            case("rho")            ; Monitors % fp_varCodes(ii) = FPVAR_RHO
            case("static-pressure"); Monitors % fp_varCodes(ii) = FPVAR_STATICPRES
            case("density")        ; Monitors % fp_varCodes(ii) = FPVAR_DENSITY
            case default           ; Monitors % fp_varCodes(ii) = FPVAR_UNKNOWN
            end select
         end do

         allocate( Monitors % fp_eID      (nFound) )
         allocate( Monitors % fp_ownsProbe(nFound) )
         allocate( Monitors % fp_lxi  (0:Nmax, nFound) )
         allocate( Monitors % fp_leta (0:Nmax, nFound) )
         allocate( Monitors % fp_lzeta(0:Nmax, nFound) )
         allocate( Monitors % fp_values_gpu(nv, nFound) )
         Monitors % fp_lxi        = 0.0_RP
         Monitors % fp_leta       = 0.0_RP
         Monitors % fp_lzeta      = 0.0_RP
         Monitors % fp_values_gpu = 0.0_RP

         do ii = 1, nFound
            owns = ( globalOwner(ii) .eq. MPI_Process % rank )
            Monitors % fp_eID(ii)       = eID_local(ii)
            Monitors % fp_ownsProbe(ii) = owns
            if ( owns ) then
               associate(eNxyz => mesh % elements(eID_local(ii)) % Nxyz)
               associate( spAxi   => NodalStorage(eNxyz(1)), &
                          spAeta  => NodalStorage(eNxyz(2)), &
                          spAzeta => NodalStorage(eNxyz(3)) )
               Monitors % fp_lxi  (0:eNxyz(1), ii) = spAxi   % lj(xi_local(1,ii))
               Monitors % fp_leta (0:eNxyz(2), ii) = spAeta  % lj(xi_local(2,ii))
               Monitors % fp_lzeta(0:eNxyz(3), ii) = spAzeta % lj(xi_local(3,ii))
               end associate
               end associate
            end if
         end do

         !$acc enter data copyin(Monitors % fp_eID)
         !$acc enter data copyin(Monitors % fp_ownsProbe)
         !$acc enter data copyin(Monitors % fp_lxi)
         !$acc enter data copyin(Monitors % fp_leta)
         !$acc enter data copyin(Monitors % fp_lzeta)
         !$acc enter data copyin(Monitors % fp_varCodes)
         !$acc enter data create(Monitors % fp_values_gpu)
      end block
#else
      block
         integer :: ii, Nmax_loc

         Monitors % fp_nOwned = 0
         Nmax_loc = 0
         do ii = 1, nFound
            if ( globalOwner(ii) .eq. MPI_Process % rank ) then
               Monitors % fp_nOwned = Monitors % fp_nOwned + 1
               Nmax_loc = max(Nmax_loc, maxval(mesh % elements(eID_local(ii)) % Nxyz))
            end if
         end do

         allocate( Monitors % fp_ownedIdx(Monitors % fp_nOwned) )
         allocate( Monitors % fp_cpu_eID  (Monitors % fp_nOwned) )
         allocate( Monitors % fp_cpu_Nx   (Monitors % fp_nOwned) )
         allocate( Monitors % fp_cpu_Ny   (Monitors % fp_nOwned) )
         allocate( Monitors % fp_cpu_Nz   (Monitors % fp_nOwned) )
         allocate( Monitors % fp_cpu_lxi  (0:Nmax_loc, Monitors % fp_nOwned) )
         allocate( Monitors % fp_cpu_leta (0:Nmax_loc, Monitors % fp_nOwned) )
         allocate( Monitors % fp_cpu_lzeta(0:Nmax_loc, Monitors % fp_nOwned) )
         allocate( Monitors % fp_cpu_varCodes(nv) )
         do ii = 1, nv
            select case (trim(Monitors % probesVariables(ii)))
            case("pressure")        ; Monitors % fp_cpu_varCodes(ii) = FPVAR_PRESSURE
            case("velocity")        ; Monitors % fp_cpu_varCodes(ii) = FPVAR_VELOCITY
            case("u")               ; Monitors % fp_cpu_varCodes(ii) = FPVAR_U
            case("v")               ; Monitors % fp_cpu_varCodes(ii) = FPVAR_V
            case("w")               ; Monitors % fp_cpu_varCodes(ii) = FPVAR_W
            case("mach")            ; Monitors % fp_cpu_varCodes(ii) = FPVAR_MACH
            case("k")               ; Monitors % fp_cpu_varCodes(ii) = FPVAR_K
            case("rho")             ; Monitors % fp_cpu_varCodes(ii) = FPVAR_RHO
            case("static-pressure") ; Monitors % fp_cpu_varCodes(ii) = FPVAR_STATICPRES
            case("density")         ; Monitors % fp_cpu_varCodes(ii) = FPVAR_DENSITY
            case default            ; Monitors % fp_cpu_varCodes(ii) = FPVAR_UNKNOWN
            end select
         end do

!        fp_ownedIdx now stores the LOCAL file-probe index (1..nFound)
!        directly - there is no global probes(:) index to translate back
!        from any more, so Monitor_ComputeFileProbesCPU's jbuf math below
!        drops the old "- fp_offset" adjustment.
         Monitors % fp_nOwned = 0
         do ii = 1, nFound
            if ( globalOwner(ii) .eq. MPI_Process % rank ) then
               Monitors % fp_nOwned = Monitors % fp_nOwned + 1
               Monitors % fp_ownedIdx(Monitors % fp_nOwned) = ii
               Monitors % fp_cpu_eID(Monitors % fp_nOwned)  = eID_local(ii)
               associate(eNxyz => mesh % elements(eID_local(ii)) % Nxyz)
               associate( spAxi   => NodalStorage(eNxyz(1)), &
                          spAeta  => NodalStorage(eNxyz(2)), &
                          spAzeta => NodalStorage(eNxyz(3)) )
               Monitors % fp_cpu_Nx(Monitors % fp_nOwned) = eNxyz(1)
               Monitors % fp_cpu_Ny(Monitors % fp_nOwned) = eNxyz(2)
               Monitors % fp_cpu_Nz(Monitors % fp_nOwned) = eNxyz(3)
               Monitors % fp_cpu_lxi  (0:eNxyz(1), Monitors % fp_nOwned) = spAxi   % lj(xi_local(1,ii))
               Monitors % fp_cpu_leta (0:eNxyz(2), Monitors % fp_nOwned) = spAeta  % lj(xi_local(2,ii))
               Monitors % fp_cpu_lzeta(0:eNxyz(3), Monitors % fp_nOwned) = spAzeta % lj(xi_local(3,ii))
               end associate
               end associate
            end if
         end do

      call system_clock(tinit_t4)
!        Sort owned probes by ascending eID so Monitor_ComputeFileProbesCPU
!        accesses mesh%elements(eID)%storage%Q sequentially, improving cache
!        reuse when many probes share the same or nearby elements.
!        Uses counting sort: O(fp_nOwned + no_of_elements).
         if (Monitors % fp_nOwned > 1) then
            block
               integer :: ne, cnt, pos_local, p_local, q_local, tmp_i, Nmax_sort
               integer, allocatable :: count_arr(:), prefix(:), sort_perm(:)
               integer, allocatable :: tmp_ownedIdx(:), tmp_eID(:), tmp_Nx(:), tmp_Ny(:), tmp_Nz(:)
               real(kind=RP), allocatable :: tmp_lxi(:,:), tmp_leta(:,:), tmp_lzeta(:,:)

               ne = mesh % no_of_elements
               allocate( count_arr(ne), prefix(ne+1) )
               count_arr = 0
               do p_local = 1, Monitors % fp_nOwned
                  count_arr(Monitors % fp_cpu_eID(p_local)) = &
                     count_arr(Monitors % fp_cpu_eID(p_local)) + 1
               end do
               prefix(1) = 1
               do cnt = 1, ne
                  prefix(cnt+1) = prefix(cnt) + count_arr(cnt)
               end do

               allocate( sort_perm(Monitors % fp_nOwned) )
               count_arr = 0
               do p_local = 1, Monitors % fp_nOwned
                  q_local   = Monitors % fp_cpu_eID(p_local)
                  pos_local = prefix(q_local) + count_arr(q_local)
                  sort_perm(pos_local) = p_local
                  count_arr(q_local) = count_arr(q_local) + 1
               end do
               deallocate(count_arr, prefix)

               Nmax_sort = ubound(Monitors % fp_cpu_lxi, 1)
               allocate( tmp_ownedIdx(Monitors % fp_nOwned) )
               allocate( tmp_eID(Monitors % fp_nOwned) )
               allocate( tmp_Nx (Monitors % fp_nOwned) )
               allocate( tmp_Ny (Monitors % fp_nOwned) )
               allocate( tmp_Nz (Monitors % fp_nOwned) )
               allocate( tmp_lxi  (0:Nmax_sort, Monitors % fp_nOwned) )
               allocate( tmp_leta (0:Nmax_sort, Monitors % fp_nOwned) )
               allocate( tmp_lzeta(0:Nmax_sort, Monitors % fp_nOwned) )

               do p_local = 1, Monitors % fp_nOwned
                  tmp_i = sort_perm(p_local)
                  tmp_ownedIdx(p_local)     = Monitors % fp_ownedIdx(tmp_i)
                  tmp_eID(p_local)          = Monitors % fp_cpu_eID(tmp_i)
                  tmp_Nx(p_local)           = Monitors % fp_cpu_Nx(tmp_i)
                  tmp_Ny(p_local)           = Monitors % fp_cpu_Ny(tmp_i)
                  tmp_Nz(p_local)           = Monitors % fp_cpu_Nz(tmp_i)
                  tmp_lxi  (:,p_local)      = Monitors % fp_cpu_lxi  (:,tmp_i)
                  tmp_leta (:,p_local)      = Monitors % fp_cpu_leta (:,tmp_i)
                  tmp_lzeta(:,p_local)      = Monitors % fp_cpu_lzeta(:,tmp_i)
               end do

               Monitors % fp_ownedIdx  = tmp_ownedIdx
               Monitors % fp_cpu_eID   = tmp_eID
               Monitors % fp_cpu_Nx    = tmp_Nx
               Monitors % fp_cpu_Ny    = tmp_Ny
               Monitors % fp_cpu_Nz    = tmp_Nz
               Monitors % fp_cpu_lxi   = tmp_lxi
               Monitors % fp_cpu_leta  = tmp_leta
               Monitors % fp_cpu_lzeta = tmp_lzeta

               deallocate(sort_perm, tmp_ownedIdx, tmp_eID, tmp_Nx, tmp_Ny, tmp_Nz, &
                          tmp_lxi, tmp_leta, tmp_lzeta)
            end block
         end if
      end block

!     Build MPI_Gatherv infrastructure: allows Monitor_UpdateFileProbes to use
!     MPI_Gatherv to root instead of MPI_Allreduce on the full fp_buf, avoiding
!     both the 43 MB zero-fill and the all-to-all reduction.
!     fp_owned_buf(p, v) = compact result buffer, size fp_nOwned*nv.
!     fp_gatherv_counts/displs: for MPI_Gatherv call (root receives from all ranks).
!     fp_gatherv_perm(k): global probe index of the k-th gathered entry (root only).
      allocate( Monitors % fp_owned_buf(max(1, Monitors % fp_nOwned * nv)) )
      Monitors % fp_owned_buf = 0.0_RP

#ifdef _HAS_MPI_
      if ( MPI_Process % doMPIAction ) then
         block
            integer :: nranks, my_owned_nv, ii, rank_offset, p_g
            integer, allocatable :: all_nowned(:), all_ownedIdx(:)
            integer, allocatable :: idx_counts(:), idx_displs(:)

            call MPI_Comm_size(MPI_COMM_WORLD, nranks, ierr)
            my_owned_nv = Monitors % fp_nOwned * nv

            if ( MPI_Process % isRoot ) then
               allocate( Monitors % fp_gatherv_counts(0:nranks-1) )
               allocate( Monitors % fp_gatherv_displs(0:nranks-1) )
            end if

            if ( MPI_Process % isRoot ) then
               call MPI_Gather(my_owned_nv, 1, MPI_INTEGER, &
                               Monitors % fp_gatherv_counts, 1, MPI_INTEGER, &
                               0, MPI_COMM_WORLD, ierr)
            else
               block
                  integer :: dummy_counts(1)
                  dummy_counts = 0
                  call MPI_Gather(my_owned_nv, 1, MPI_INTEGER, &
                                  dummy_counts, 1, MPI_INTEGER, &
                                  0, MPI_COMM_WORLD, ierr)
               end block
            end if

            if ( MPI_Process % isRoot ) then
               Monitors % fp_gatherv_displs(0) = 0
               do ii = 1, nranks-1
                  Monitors % fp_gatherv_displs(ii) = Monitors % fp_gatherv_displs(ii-1) + &
                                                     Monitors % fp_gatherv_counts(ii-1)
               end do
               allocate( Monitors % fp_gatherv_perm(nFound) )
               allocate( all_nowned(0:nranks-1) )
               all_nowned = Monitors % fp_gatherv_counts / max(1, nv)
               allocate( idx_counts(0:nranks-1), idx_displs(0:nranks-1) )
               idx_counts = all_nowned
               idx_displs(0) = 0
               do ii = 1, nranks-1
                  idx_displs(ii) = idx_displs(ii-1) + idx_counts(ii-1)
               end do
               allocate( all_ownedIdx(nFound) )
               call MPI_Gatherv(Monitors % fp_ownedIdx, Monitors % fp_nOwned, MPI_INTEGER, &
                                all_ownedIdx, idx_counts, idx_displs, MPI_INTEGER, &
                                0, MPI_COMM_WORLD, ierr)
               do p_g = 1, nFound
                  Monitors % fp_gatherv_perm(p_g) = all_ownedIdx(p_g)
               end do
               deallocate(all_nowned, all_ownedIdx, idx_counts, idx_displs)
            else
               block
                  integer :: dummy_recv(1), dummy_counts(1), dummy_displs(1)
                  dummy_recv = 0 ; dummy_counts = 0 ; dummy_displs = 0
                  call MPI_Gatherv(Monitors % fp_ownedIdx, Monitors % fp_nOwned, MPI_INTEGER, &
                                   dummy_recv, dummy_counts, dummy_displs, MPI_INTEGER, &
                                   0, MPI_COMM_WORLD, ierr)
               end block
            end if
         end block
      end if
#endif

      call system_clock(tinit_t5)

      if ( MPI_Process % isRoot ) then
         write(STD_OUT,'(/,30X,A)') "--- InitializeProbesFromFile timing (rank 0) ---"
         write(STD_OUT,'(30X,A,I0)')    "  Total probes         : ", nFound
         write(STD_OUT,'(30X,A,I0)')    "  Owned by this rank   : ", Monitors % fp_nOwned
         write(STD_OUT,'(30X,A,ES12.4,A)') "  File read + Bcast    : ", &
            real(tinit_t1-tinit_t0,8)/real(tinit_rate,8), " s"
         write(STD_OUT,'(30X,A,ES12.4,A)') "  FindPointWithCoords  : ", &
            real(tinit_t2-tinit_t1,8)/real(tinit_rate,8), " s"
         write(STD_OUT,'(30X,A,ES12.4,A)') "  Allreduce ownership  : ", &
            real(tinit_t3-tinit_t2,8)/real(tinit_rate,8), " s"
         write(STD_OUT,'(30X,A,ES12.4,A)') "  SoA build (pass 2)   : ", &
            real(tinit_t4-tinit_t3,8)/real(tinit_rate,8), " s"
         write(STD_OUT,'(30X,A,ES12.4,A)') "  Sort by eID          : ", &
            real(tinit_t5-tinit_t4,8)/real(tinit_rate,8), " s"
         write(STD_OUT,'(30X,A,ES12.4,A)') "  TOTAL                : ", &
            real(tinit_t5-tinit_t0,8)/real(tinit_rate,8), " s"
      end if
#endif

   end subroutine InitializeProbesFromFile

! ============================================================
!  HDF5 collective I/O for file-based probes
!  Compiled only when HAS_HDF5 is defined.
! ============================================================

#ifdef HAS_HDF5
   subroutine Monitor_InitFileProbesHDF5(self, no_of_fileProbes)
!
!     Creates the HDF5 output file for bulk file-probes with
!     the following structure:
!
!       /coordinates  (3, nProbes)   — fixed, written once
!       /time         (extendible)   — one value per saved step
!       /iteration    (extendible)   — one value per saved step
!       /<varname>    (nProbes, ext) — one row per saved step
!
      use HDF5
      use MPI_Process_Info
      implicit none
      class(Monitor_t), intent(inout) :: self
      integer,          intent(in)    :: no_of_fileProbes
!
!     ---------------
!     Local variables
!     ---------------
!
      integer(HID_T)   :: file_id, dset_id, dspace_id, dcpl_id
      integer(HSIZE_T) :: dims2(2), maxdims2(2), chunk2(2)
      integer(HSIZE_T) :: dims1(1), maxdims1(1), chunk1(1)
      integer          :: iError, v, nv
      character(len=LINE_LENGTH) :: fname

      if ( .not. MPI_Process % isRoot ) return

      nv = size(self % probesVariables)

      write(fname,'(A,A)') trim(self % probes_solution_file), ".probes.h5"

      call h5open_f(iError)
      call h5fcreate_f(trim(fname), H5F_ACC_TRUNC_F, file_id, iError)

      ! /coordinates  (3, nProbes) — fixed at creation, straight from fp_x
      ! (populated for every file-probe by InitializeProbesFromFile)
      dims2 = [ int(3, HSIZE_T), int(no_of_fileProbes, HSIZE_T) ]
      call h5screate_simple_f(2, dims2, dspace_id, iError)
      call h5dcreate_f(file_id, "coordinates", H5T_NATIVE_DOUBLE, dspace_id, dset_id, iError)
      call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, self % fp_x, dims2, iError)
      call h5dclose_f(dset_id, iError)
      call h5sclose_f(dspace_id, iError)

      ! /time  (extendible 1-D)
      dims1(1)    = 0
      maxdims1(1) = H5S_UNLIMITED_F
      chunk1(1)   = int(max(BUFFER_SIZE, 1), HSIZE_T)
      call h5screate_simple_f(1, dims1, dspace_id, iError, maxdims1)
      call h5pcreate_f(H5P_DATASET_CREATE_F, dcpl_id, iError)
      call h5pset_chunk_f(dcpl_id, 1, chunk1, iError)
      call h5dcreate_f(file_id, "time", H5T_NATIVE_DOUBLE, dspace_id, dset_id, iError, dcpl_id)
      call h5dclose_f(dset_id, iError)
      call h5sclose_f(dspace_id, iError)
      call h5pclose_f(dcpl_id, iError)

      ! /iteration  (extendible 1-D)
      call h5screate_simple_f(1, dims1, dspace_id, iError, maxdims1)
      call h5pcreate_f(H5P_DATASET_CREATE_F, dcpl_id, iError)
      call h5pset_chunk_f(dcpl_id, 1, chunk1, iError)
      call h5dcreate_f(file_id, "iteration", H5T_NATIVE_INTEGER, dspace_id, dset_id, iError, dcpl_id)
      call h5dclose_f(dset_id, iError)
      call h5sclose_f(dspace_id, iError)
      call h5pclose_f(dcpl_id, iError)

      ! /<varname>  (extendible × nProbes) — h5ls shows {nProbes, Inf}
      ! chunk2(1)=1: one time step per chunk avoids pre-allocation waste when
      ! the save-timestep filter makes n_write << BUFFER_SIZE per flush.
      dims2(1)    = 0
      dims2(2)    = int(no_of_fileProbes, HSIZE_T)
      maxdims2(1) = H5S_UNLIMITED_F
      maxdims2(2) = int(no_of_fileProbes, HSIZE_T)
      chunk2(1)   = int(1, HSIZE_T)
      chunk2(2)   = int(no_of_fileProbes, HSIZE_T)

      do v = 1, nv
         call h5screate_simple_f(2, dims2, dspace_id, iError, maxdims2)
         call h5pcreate_f(H5P_DATASET_CREATE_F, dcpl_id, iError)
         call h5pset_chunk_f(dcpl_id, 2, chunk2, iError)
         call h5dcreate_f(file_id, trim(self % probesVariables(v)), H5T_NATIVE_DOUBLE, &
                          dspace_id, dset_id, iError, dcpl_id)
         call h5dclose_f(dset_id, iError)
         call h5sclose_f(dspace_id, iError)
         call h5pclose_f(dcpl_id, iError)
      end do

      ! Keep the file open so Monitor_WriteFileProbesHDF5 can reuse the
      ! handle without paying a filesystem open/close on every probe save.
      hdf5_fp_fid  = file_id
      hdf5_fp_open = .true.

   end subroutine Monitor_InitFileProbesHDF5

   subroutine Monitor_WriteFileProbesHDF5(self, iter, t, no_of_lines)
!
!     Appends one buffer of file-probe data to the HDF5 output file.
!     Only lines that satisfy the saveTimestep filter are written.
!     Called from Monitor_WriteToFile in place of the per-probe ASCII loop.
!
      use HDF5
      use MPI_Process_Info
      implicit none
      class(Monitor_t), intent(inout) :: self
      integer,          intent(in)    :: iter(:)
      real(kind=RP),    intent(in)    :: t(:)
      integer,          intent(in)    :: no_of_lines
!
!     ---------------
!     Local variables
!     ---------------
!
      integer(HID_T)   :: file_id, dset_id, dspace_id, mspace_id
      integer(HSIZE_T) :: cur1(1), max1(1), off1(1), cnt1(1)
      integer(HSIZE_T) :: new1(1)
      integer(HSIZE_T) :: cur2(2), max2(2), off2(2), cnt2(2)
      integer(HSIZE_T) :: new2(2)
      integer          :: iError, i, j, k, v, nfp, nv, fp_offset, n_write
      character(len=LINE_LENGTH) :: fname
      logical,      allocatable :: wmask(:)
      integer,      allocatable :: ibuf(:)
      real(kind=RP), allocatable :: rbuf(:), vbuf(:)

      if ( .not. MPI_Process % isRoot ) return

      nfp       = self % no_of_fileProbes
      nv        = size(self % probesVariables)
      fp_offset = self % no_of_probes - nfp

      ! All no_of_lines slots are written: the caller (Monitor_FlushFileProbesNow)
      ! already applied the probeFileSaveTimestep filter and updated fp_lastSavedTime
      ! before invoking this routine.  No second filter here.
      allocate( wmask(no_of_lines) )
      wmask  = .true.
      n_write = no_of_lines

      ! Reuse the persistent file handle opened by Monitor_InitFileProbesHDF5.
      ! On restart (FirstCall=.false.) the init routine is skipped, so open
      ! lazily here if the handle is not yet valid.
      if ( hdf5_fp_open ) then
         file_id = hdf5_fp_fid
      else
         write(fname,'(A,A)') trim(self % probes_solution_file), ".probes.h5"
         call h5open_f(iError)
         call h5fopen_f(trim(fname), H5F_ACC_RDWR_F, file_id, iError)
         hdf5_fp_fid  = file_id
         hdf5_fp_open = .true.
      end if

      ! Query current extent of /time to get append offset
      call h5dopen_f(file_id, "time", dset_id, iError)
      call h5dget_space_f(dset_id, dspace_id, iError)
      call h5sget_simple_extent_dims_f(dspace_id, cur1, max1, iError)
      call h5sclose_f(dspace_id, iError)
      call h5dclose_f(dset_id, iError)
      ! cur1(1) = number of time steps already written

      cnt1(1) = int(n_write, HSIZE_T)
      off1(1) = cur1(1)
      new1(1) = cur1(1) + cnt1(1)

      ! Pack filtered iteration values
      allocate( ibuf(n_write) )
      j = 0
      do i = 1, no_of_lines
         if ( .not. wmask(i) ) cycle
         j = j + 1
         ibuf(j) = iter(i)
      end do

      call h5dopen_f(file_id, "iteration", dset_id, iError)
      call h5dextend_f(dset_id, new1, iError)
      call h5dget_space_f(dset_id, dspace_id, iError)
      call h5sselect_hyperslab_f(dspace_id, H5S_SELECT_SET_F, off1, cnt1, iError)
      call h5screate_simple_f(1, cnt1, mspace_id, iError)
      call h5dwrite_f(dset_id, H5T_NATIVE_INTEGER, ibuf, cnt1, iError, mspace_id, dspace_id)
      call h5sclose_f(mspace_id, iError)
      call h5sclose_f(dspace_id, iError)
      call h5dclose_f(dset_id, iError)
      deallocate(ibuf)

      ! Pack filtered time values
      allocate( rbuf(n_write) )
      j = 0
      do i = 1, no_of_lines
         if ( .not. wmask(i) ) cycle
         j = j + 1
         rbuf(j) = t(i)
      end do

      call h5dopen_f(file_id, "time", dset_id, iError)
      call h5dextend_f(dset_id, new1, iError)
      call h5dget_space_f(dset_id, dspace_id, iError)
      call h5sselect_hyperslab_f(dspace_id, H5S_SELECT_SET_F, off1, cnt1, iError)
      call h5screate_simple_f(1, cnt1, mspace_id, iError)
      call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, rbuf, cnt1, iError, mspace_id, dspace_id)
      call h5sclose_f(mspace_id, iError)
      call h5sclose_f(dspace_id, iError)
      call h5dclose_f(dset_id, iError)
      deallocate(rbuf)

      ! Write each variable: one extend + one write per flush (avoids chunk waste).
      ! vbuf2(nfp, n_write) packs all filtered timesteps into a contiguous 2D block.
      off2(1) = cur1(1)             ! append after existing time steps
      off2(2) = int(0,   HSIZE_T)
      cnt2(1) = int(n_write, HSIZE_T)
      cnt2(2) = int(nfp,    HSIZE_T)
      new2(1) = cur1(1) + int(n_write, HSIZE_T)
      new2(2) = int(nfp, HSIZE_T)

      allocate( vbuf(nfp * n_write) )

      do v = 1, nv
         call h5dopen_f(file_id, trim(self % probesVariables(v)), dset_id, iError)

         ! Pack into column-major 2D layout vbuf(n_write, nfp):
         ! time index j varies fastest (Fortran column-major with cnt2=[n_write,nfp])
         j = 0
         do i = 1, no_of_lines
            if ( .not. wmask(i) ) cycle
            j = j + 1
            do k = 1, nfp
#ifdef _OPENACC
               ! fp_values_gpu layout: (nv, nfp) -> variable first, probe second
               vbuf( j + (k-1)*n_write ) = self % fp_values_gpu(v, k)
#else
               ! fp_buf layout: (nfp*nv) -> probe first: (k-1)*nv + v
               vbuf( j + (k-1)*n_write ) = self % fp_buf( (k-1)*nv + v )
#endif
            end do
         end do

         call h5dextend_f(dset_id, new2, iError)
         call h5dget_space_f(dset_id, dspace_id, iError)
         call h5sselect_hyperslab_f(dspace_id, H5S_SELECT_SET_F, off2, cnt2, iError)
         call h5screate_simple_f(2, cnt2, mspace_id, iError)
         call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, vbuf, cnt2, iError, mspace_id, dspace_id)
         call h5sclose_f(mspace_id, iError)
         call h5sclose_f(dspace_id, iError)

         call h5dclose_f(dset_id, iError)
      end do

      deallocate(vbuf)
      deallocate(wmask)
      ! File stays open; closed by Monitor_Destruct via Monitor_CloseHDF5FP.

   end subroutine Monitor_WriteFileProbesHDF5

   subroutine Monitor_CloseHDF5FP()
!     Close the persistent file-probe HDF5 handle opened by
!     Monitor_InitFileProbesHDF5 (or lazily by Monitor_WriteFileProbesHDF5).
      implicit none
      integer :: iError
      if ( hdf5_fp_open ) then
         call h5fclose_f(hdf5_fp_fid, iError)
         call h5close_f(iError)
         hdf5_fp_open = .false.
      end if
   end subroutine Monitor_CloseHDF5FP

#endif   ! HAS_HDF5

#endif   ! FLOW

end module MonitorsClass
!
!///////////////////////////////////////////////////////////////////////////////////
!
