!
!///////////////////////////////////////////////////////////////////////////////////////
!
!  horses2probes - Extract solution variables at probe points from HORSES3D solution files
!
!  Usage (control file):
!     horses2probes control_file.convert
!
!  Usage (command line):
!     horses2probes mesh.hmesh solution.hsol probes.dat [--output-variables=p,u,v,w]
!
!  Control file format:
!     hmesh file= Mesh.hmesh
!     hsol file= Solution.hsol
!     probes file= probes.dat
!     output variables= p, u, v, w
!     flow equations= ns
!     boundary file= optional
!     partition file= optional
!
!  Probes file format (one probe per line, lines starting with # are comments):
!     probe_name  x  y  z
!
!  Probe location cache:
!     On the first run, probe locations (element ID + reference coordinates) are
!     saved to <probesfile>.pcache. Subsequent runs with the same mesh and probes
!     file load directly from cache, skipping the mesh search entirely.
!
!///////////////////////////////////////////////////////////////////////////////////////
!
program horses2probes
   use SMConstants
   use Storage
   use SolutionFile
   use SharedSpectralBasis
   use Headers
   use MPI_Process_Info
   use OutputVariables, only: outScale, hasVariablesFlag, askedVariables, Lreference, getNoOfCommas
   use Horses2probesModule
   implicit none

   integer                                 :: no_of_solutions
   character(len=LINE_LENGTH)              :: meshName
   character(len=LINE_LENGTH)              :: probesFileName
   character(len=LINE_LENGTH), allocatable :: solutionNames(:)
   integer                                 :: iSol
   type(Mesh_t)                            :: mesh
   type(ProbeData_t), allocatable          :: probes(:)
   integer                                 :: no_of_probes
   integer(kind=8)                         :: t_sol, t_end, sol_rate

   call MPI_Process % Init
   call Main_Header("HORSES Probes to TecPlot extraction utility", __DATE__, __TIME__)
!
!  Parse the control file or command line arguments
!  -------------------------------------------------
   call getHorses2probesArgs(meshName, probesFileName, no_of_solutions, solutionNames)
!
!  Construct spectral basis
!  ------------------------
   call ConstructSpectralBasis
!
!  Read mesh
!  ---------
   call mesh % ReadMesh(meshName)
!
!  Find probes in mesh (or load from cache) — done only once
!  ----------------------------------------------------------
   call FindAndCacheProbes(mesh, trim(meshName), trim(probesFileName), probes, no_of_probes)
!
!  Process each solution file
!  --------------------------
   if (no_of_probes .gt. 0) then
      do iSol = 1, no_of_solutions
         write(STD_OUT,'(/,/)')
         call Section_Header("Processing solution")
         call system_clock(t_sol, sol_rate)
         call mesh % ReadSolution(solutionNames(iSol))
         call system_clock(t_end)
         write(STD_OUT,'(30X,A,A)') "-> Solution: ", trim(solutionNames(iSol))
         write(STD_OUT,'(30X,A,F8.3,A)') "-> Read time:  ", &
            real(t_end-t_sol,8)/real(sol_rate,8), " s"
         call InterpolateAndWriteProbes(mesh, probes, no_of_probes, solutionNames(iSol), trim(probesFileName))
         call system_clock(t_end)
         write(STD_OUT,'(30X,A,F8.3,A)') "-> Total time: ", &
            real(t_end-t_sol,8)/real(sol_rate,8), " s"
      end do
   end if
!
!  Cleanup
!  -------
   if (allocated(probes)) deallocate(probes)
   call DestructSpectralBasis
   call MPI_Process % Close

   contains

   subroutine getHorses2probesArgs(meshName, probesFileName, no_of_solutions, solutionNames)
      use FTValueDictionaryClass, only: FTValueDictionary
      use FileReaders           , only: ReadControlFile
      use FileReadingUtilities  , only: getCharArrayFromString
      use Utilities             , only: toLower
      implicit none
      character(len=*),                        intent(out) :: meshName
      character(len=*),                        intent(out) :: probesFileName
      integer,                                 intent(out) :: no_of_solutions
      character(len=LINE_LENGTH), allocatable, intent(out) :: solutionNames(:)
!
!     ---------------
!     Local variables
!     ---------------
!
      integer                                              :: no_of_arguments, i, pos, pos2, fID, io
      character(len=LINE_LENGTH)                           :: arg, inputResultName
      type(FTValueDictionary)                              :: controlVariables
      logical                                              :: useControlFile

      no_of_arguments = command_argument_count()
      no_of_solutions = 0
!
!     Determine input mode: control file (1 arg) or command line (>1 args)
!     ---------------------------------------------------------------------
      useControlFile = (no_of_arguments .eq. 1)

      if (no_of_arguments .eq. 0) then
         write(STD_OUT,'(A)') "Usage: horses2probes control_file.convert"
         write(STD_OUT,'(A)') "       horses2probes mesh.hmesh solution.hsol probes.dat [--output-variables=p,u,v,w]"
         call exit(1)
      end if 
!
!     ---- Control file mode ----
!     ---------------------------
      if (useControlFile) then
         call controlVariables % initWithSize(16)
         call ReadControlFile(controlVariables)

         if (.not. controlVariables % containsKey("hmesh file")) then
            write(STD_OUT,'(A)') "ERROR: 'hmesh file' not specified in control file."
            call exit(1)
         end if
         meshName = controlVariables % stringValueForKey("hmesh file", LINE_LENGTH)

         if (.not. controlVariables % containsKey("probes file")) then
            write(STD_OUT,'(A)') "ERROR: 'probes file' not specified in control file."
            call exit(1)
         end if
         probesFileName = controlVariables % stringValueForKey("probes file", LINE_LENGTH)

         if (.not. controlVariables % containsKey("hsol file")) then
            write(STD_OUT,'(A)') "ERROR: 'hsol file' not specified in control file."
            call exit(1)
         end if

         inputResultName = controlVariables % stringValueForKey("hsol file", LINE_LENGTH)
         no_of_solutions = getNoOfCommas(trim(inputResultName)) + 1
         allocate(solutionNames(no_of_solutions))

         if (no_of_solutions .eq. 1) then
            solutionNames(1) = adjustl(trim(inputResultName))
         else
            pos = 0
            do i = 1, no_of_solutions - 1
               pos2 = index(trim(inputResultName(pos+1:)), ",") + pos
               solutionNames(i) = trim(adjustl(inputResultName(pos+1:pos2-1)))
               pos = pos2
            end do
            pos = index(trim(inputResultName), ",", BACK=.true.)
            solutionNames(no_of_solutions) = adjustl(trim(inputResultName(pos+1:)))
         end if

         hasVariablesFlag = controlVariables % containsKey("output variables")
         if (hasVariablesFlag) askedVariables = controlVariables % stringValueForKey("output variables", LINE_LENGTH)
         outScale   = .not. controlVariables % logicalValueForKey("dimensionless")
         Lreference = controlVariables % getValueOrDefault("reference length (m)", 1.0_RP)

         hasBoundaries = controlVariables % containsKey("boundary file")
         if (hasBoundaries) boundaryFileName = controlVariables % stringValueForKey("boundary file", LINE_LENGTH)
         hasMPIranks = controlVariables % containsKey("partition file")
         if (hasMPIranks) partitionFileName = controlVariables % stringValueForKey("partition file", LINE_LENGTH)

         if (controlVariables % containsKey("flow equations")) then
            flowEq = controlVariables % stringValueForKey("flow equations", LINE_LENGTH)
            call toLower(flowEq)
         else
            flowEq = "ns"
         end if
!
!     ---- Command line mode ----
!     ---------------------------
      else
         outScale         = .true.
         hasVariablesFlag = .false.
         hasBoundaries    = .false.
         hasMPIranks      = .false.
         flowEq           = "ns"
         Lreference       = 1.0_RP
         probesFileName   = ""

         do i = 1, no_of_arguments
            call get_command_argument(i, arg)
            open(newunit=fID, file=trim(arg), action="read", form="unformatted", &
                 access="stream", iostat=io)
            close(fID)
            if (io .ne. 0) cycle
            if (getSolutionFileType(arg) .eq. MESH_FILE .or. &
                getSolutionFileType(arg) .eq. ZONE_MESH_FILE) then
               meshName = trim(arg)
               exit
            end if
         end do

         no_of_solutions = 0
         do i = 1, no_of_arguments
            call get_command_argument(i, arg)
            open(newunit=fID, file=trim(arg), action="read", form="unformatted", &
                 access="stream", iostat=io)
            close(fID)
            if (io .ne. 0) cycle
            if (getSolutionFileType(arg) .ne. MESH_FILE .and. &
                getSolutionFileType(arg) .ne. ZONE_MESH_FILE) then
               no_of_solutions = no_of_solutions + 1
            end if
         end do

         allocate(solutionNames(max(1, no_of_solutions)))
         pos = 0
         do i = 1, no_of_arguments
            call get_command_argument(i, arg)
            if (arg(1:2) .eq. "--") then
               if (index(arg, "--output-variables=") .ne. 0) then
                  hasVariablesFlag = .true.
                  askedVariables = arg(len("--output-variables=")+1:)
               else if (index(arg, "--dimensionless") .ne. 0) then
                  outScale = .false.
               end if
               cycle
            end if

            open(newunit=fID, file=trim(arg), action="read", form="unformatted", &
                 access="stream", iostat=io)
            close(fID)

            if (io .eq. 0) then
               if (getSolutionFileType(arg) .ne. MESH_FILE .and. &
                   getSolutionFileType(arg) .ne. ZONE_MESH_FILE) then
                  pos = pos + 1
                  solutionNames(pos) = trim(arg)
               end if
            else
               open(newunit=fID, file=trim(arg), action="read", status="old", iostat=io)
               close(fID)
               if (io .eq. 0 .and. len_trim(probesFileName) .eq. 0) then
                  probesFileName = trim(arg)
               end if
            end if
         end do

         if (len_trim(probesFileName) .eq. 0) then
            write(STD_OUT,'(A)') "ERROR: No probes file found in arguments."
            write(STD_OUT,'(A)') "Usage: horses2probes mesh.hmesh solution.hsol probes.dat"
            call exit(1)
         end if
      end if

   end subroutine getHorses2probesArgs

end program horses2probes