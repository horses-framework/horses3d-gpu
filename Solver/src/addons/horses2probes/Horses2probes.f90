#include "Includes.h"
module Horses2probesModule
!
!  ///////////////////////////////////////////////////////////////////////////
!
!     This module implements the horses2probes functionality:
!     - Reads a probes file with (name, x, y, z) entries
!     - Finds each probe point in the mesh using Newton iteration
!       accelerated with element bounding boxes
!     - Caches probe locations (eID + xi) to a binary file for reuse
!     - Interpolates the solution at each probe position using Lagrange basis
!     - Writes output variables to a Tecplot .tec file
!
!  ///////////////////////////////////////////////////////////////////////////
!
   use SMConstants
   use Storage
   use SolutionFile
   use SharedSpectralBasis
   use NodalStorageClass
   use OutputVariables
   use PhysicsStorage
   use Headers
   implicit none

   private
   public ProbeData_t
   public ReadProbesFile
   public FindAndCacheProbes
   public InterpolateAndWriteProbes

   integer,       parameter :: MAX_NEWTON = 100
   real(kind=RP), parameter :: NEWTON_TOL = 1.0e-10_RP
   real(kind=RP), parameter :: INSIDE_TOL = 1.0e-6_RP
   real(kind=RP), parameter :: BBOX_TOL   = 1.0e-6_RP
   integer,       parameter :: MAX_WARN   = 20
#define PRECISION_FORMAT "(E18.10)"

   type ProbeData_t
      character(len=LINE_LENGTH) :: name
      real(kind=RP)              :: x(NDIM)
      real(kind=RP)              :: xi(NDIM)
      integer                    :: eID
      logical                    :: found
   end type ProbeData_t

   type ElementBBox_t
      real(kind=RP) :: xmin(NDIM)
      real(kind=RP) :: xmax(NDIM)
   end type ElementBBox_t

   contains

   pure function getBaseName(fullPath) result(base)
      character(len=*), intent(in) :: fullPath
      character(len=LINE_LENGTH)   :: base
      integer :: slashPos, dotPos

      slashPos = index(trim(fullPath), '/', BACK=.true.)
      base     = fullPath(slashPos+1:)
      dotPos   = index(trim(base), '.', BACK=.true.)
      if (dotPos > 1) base = base(1:dotPos-1)

   end function getBaseName

   pure function getDirName(fullPath) result(dir)
      character(len=*), intent(in) :: fullPath
      character(len=LINE_LENGTH)   :: dir
      integer :: slashPos

      slashPos = index(trim(fullPath), '/', BACK=.true.)
      if (slashPos > 0) then
         dir = fullPath(1:slashPos)
      else
         dir = ""
      end if

   end function getDirName

!
!/////////////////////////////////////////////////////////////////////////////////////
!
!     Read probe positions from file.
!     Format: one probe per line ->  name  x  y  z
!     Lines starting with # are comments. Name is optional.
!
!/////////////////////////////////////////////////////////////////////////////////////
!
   subroutine ReadProbesFile(fileName, probes, no_of_probes)
      implicit none
      character(len=*),               intent(in)  :: fileName
      type(ProbeData_t), allocatable, intent(out) :: probes(:)
      integer,                        intent(out) :: no_of_probes
!
!     ---------------
!     Local variables
!     ---------------
!
      integer            :: fid, io, lineCount, pID, nread
      character(len=512) :: line
      real(kind=RP)      :: x, y, z
      character(len=LINE_LENGTH) :: name

      write(STD_OUT,'(/)')
      call Section_Header("Reading probes file")
      write(STD_OUT,'(30X,A,A)') "-> File: ", trim(fileName)

      open(newunit=fid, file=trim(fileName), action="read", status="old", iostat=io)
      if (io .ne. 0) then
         write(STD_OUT,'(30X,A,A,A)') "ERROR: Could not open probes file '", trim(fileName), "'"
         no_of_probes = 0
         return
      end if

      lineCount = 0
      do
         read(fid, '(A)', iostat=io) line
         if (io .ne. 0) exit
         line = adjustl(line)
         if (len_trim(line) .eq. 0) cycle
         if (line(1:1) .eq. '#') cycle
         lineCount = lineCount + 1
      end do
      close(fid)

      no_of_probes = lineCount
      if (no_of_probes .eq. 0) return
      allocate(probes(no_of_probes))

      open(newunit=fid, file=trim(fileName), action="read", status="old")
      pID = 0
      do
         read(fid, '(A)', iostat=io) line
         if (io .ne. 0) exit
         line = adjustl(line)
         if (len_trim(line) .eq. 0) cycle
         if (line(1:1) .eq. '#') cycle
         pID = pID + 1
         read(line, *, iostat=nread) name, x, y, z
         if (nread .ne. 0) then
            read(line, *, iostat=nread) x, y, z
            if (nread .ne. 0) then
               write(STD_OUT,'(30X,A,I0,A)') "WARNING: Cannot parse probe line ", pID, ". Skipping."
               pID = pID - 1
               cycle
            end if
            write(probes(pID)%name, '(A,I0)') "probe_", pID
         else
            probes(pID)%name = trim(name)
         end if
         probes(pID)%x(1) = x
         probes(pID)%x(2) = y
         probes(pID)%x(3) = z
         probes(pID)%found = .false.
         probes(pID)%eID   = 0
         probes(pID)%xi    = 0.0_RP
      end do
      no_of_probes = pID
      close(fid)
      write(STD_OUT,'(30X,A,I0,A)') "-> Found ", no_of_probes, " probe(s)."

   end subroutine ReadProbesFile

!
!/////////////////////////////////////////////////////////////////////////////////////
!
!     Find all probes in the mesh (or load from cache if available).
!     Cache file is named <probesfile_base>.pcache and stores eID + xi
!     for each probe. It is invalidated if meshName or probesFileName change.
!
!/////////////////////////////////////////////////////////////////////////////////////
!
   subroutine FindAndCacheProbes(mesh, meshName, probesFileName, probes, no_of_probes)
      use FileReadingUtilities, only: getFileName
      implicit none
      type(Mesh_t),                    intent(in)  :: mesh
      character(len=*),                intent(in)  :: meshName
      character(len=*),                intent(in)  :: probesFileName
      type(ProbeData_t), allocatable,  intent(out) :: probes(:)
      integer,                         intent(out) :: no_of_probes
!
!     ---------------
!     Local variables
!     ---------------
!
      character(len=LINE_LENGTH)        :: cacheFileName
      logical                           :: cacheLoaded
      type(ElementBBox_t), allocatable  :: bbox(:)
      integer                           :: eID
      integer(kind=8)                   :: t0, t1, count_rate

      cacheFileName = trim(getFileName(probesFileName)) // ".pcache"

      write(STD_OUT,'(/)')
      call Section_Header("Probe location lookup")
      write(STD_OUT,'(30X,A,A)') "-> Cache file: ", trim(cacheFileName)

      call system_clock(t0, count_rate)
      call ReadProbeCache(cacheFileName, meshName, probesFileName, probes, no_of_probes, cacheLoaded)

      if (cacheLoaded) then
         call system_clock(t1)
         write(STD_OUT,'(30X,A,I0,A)') "-> Loaded ", no_of_probes, " probe locations from cache."
         write(STD_OUT,'(30X,A,F8.3,A)') "-> Cache load time: ", &
            real(t1-t0,RP)/real(count_rate,RP), " s"
         return
      end if

      write(STD_OUT,'(30X,A)') "-> Cache not found or stale. Searching mesh..."
!
!     Read probe positions from file
!     --------------------------------
      call ReadProbesFile(probesFileName, probes, no_of_probes)
      if (no_of_probes .eq. 0) return
!
!     Initialize spectral basis for mesh geometry orders
!     ---------------------------------------------------
      do eID = 1, mesh % no_of_elements
         associate(e => mesh % elements(eID))
         call addNewSpectralBasis(spA, e % Nmesh, mesh % nodeType)
         end associate
      end do
!
!     Build element bounding boxes
!     ----------------------------
      call BuildBoundingBoxes(mesh, bbox)
!
!     Find each probe in the mesh
!     ----------------------------
      call Section_Header("Finding probes in mesh")
      call system_clock(t0, count_rate)
      call FindAllProbesInMesh(mesh, bbox, probes, no_of_probes)
      call system_clock(t1)
      write(STD_OUT,'(30X,A,F8.3,A)') "-> Mesh search time: ", &
         real(t1-t0,RP)/real(count_rate,RP), " s"

      deallocate(bbox)
!
!     Save cache for future runs
!     --------------------------
      call SaveProbeCache(cacheFileName, meshName, probesFileName, probes, no_of_probes)

   end subroutine FindAndCacheProbes

!
!/////////////////////////////////////////////////////////////////////////////////////
!
!     Read probe location cache from binary file.
!
!/////////////////////////////////////////////////////////////////////////////////////
!
   subroutine ReadProbeCache(cacheFileName, meshName, probesFileName, probes, no_of_probes, success)
      implicit none
      character(len=*),               intent(in)  :: cacheFileName
      character(len=*),               intent(in)  :: meshName
      character(len=*),               intent(in)  :: probesFileName
      type(ProbeData_t), allocatable, intent(out) :: probes(:)
      integer,                        intent(out) :: no_of_probes
      logical,                        intent(out) :: success
!
!     ---------------
!     Local variables
!     ---------------
!
      integer :: fid, io, pID
      character(len=LINE_LENGTH) :: stored_mesh, stored_probes

      success      = .false.
      no_of_probes = 0

      open(newunit=fid, file=trim(cacheFileName), action="read", form="unformatted", &
           access="stream", status="old", iostat=io)
      if (io .ne. 0) return

      read(fid, iostat=io) stored_mesh  ;  if (io .ne. 0) then ; close(fid) ; return ; end if
      read(fid, iostat=io) stored_probes;  if (io .ne. 0) then ; close(fid) ; return ; end if

      if (trim(stored_mesh) .ne. trim(meshName) .or. &
          trim(stored_probes) .ne. trim(probesFileName)) then
         write(STD_OUT,'(30X,A)') "-> Cache is stale (mesh or probes file changed). Recomputing."
         close(fid)
         return
      end if

      read(fid, iostat=io) no_of_probes ; if (io .ne. 0) then ; close(fid) ; return ; end if

      allocate(probes(no_of_probes))
      do pID = 1, no_of_probes
         read(fid, iostat=io) probes(pID) % name, probes(pID) % x, &
                               probes(pID) % eID, probes(pID) % xi, probes(pID) % found
         if (io .ne. 0) then
            write(STD_OUT,'(30X,A)') "WARNING: Cache read error. Recomputing."
            deallocate(probes)
            no_of_probes = 0
            close(fid)
            return
         end if
      end do

      close(fid)
      success = .true.

   end subroutine ReadProbeCache

!
!/////////////////////////////////////////////////////////////////////////////////////
!
!     Save probe location cache to binary file.
!
!/////////////////////////////////////////////////////////////////////////////////////
!
   subroutine SaveProbeCache(cacheFileName, meshName, probesFileName, probes, no_of_probes)
      implicit none
      character(len=*),  intent(in) :: cacheFileName
      character(len=*),  intent(in) :: meshName
      character(len=*),  intent(in) :: probesFileName
      type(ProbeData_t), intent(in) :: probes(:)
      integer,           intent(in) :: no_of_probes
!
!     ---------------
!     Local variables
!     ---------------
!
      integer :: fid, io, pID
      character(len=LINE_LENGTH) :: stored_mesh, stored_probes

      stored_mesh   = trim(meshName)
      stored_probes = trim(probesFileName)

      open(newunit=fid, file=trim(cacheFileName), action="write", form="unformatted", &
           access="stream", status="replace", iostat=io)
      if (io .ne. 0) then
         write(STD_OUT,'(30X,A,A)') "WARNING: Could not write cache file: ", trim(cacheFileName)
         return
      end if

      write(fid) stored_mesh
      write(fid) stored_probes
      write(fid) no_of_probes
      do pID = 1, no_of_probes
         write(fid) probes(pID) % name, probes(pID) % x, &
                    probes(pID) % eID, probes(pID) % xi, probes(pID) % found
      end do
      close(fid)

      write(STD_OUT,'(30X,A,A)') "-> Cache saved: ", trim(cacheFileName)

   end subroutine SaveProbeCache

!
!/////////////////////////////////////////////////////////////////////////////////////
!
!     Build axis-aligned bounding boxes for all elements.
!
!/////////////////////////////////////////////////////////////////////////////////////
!
   subroutine BuildBoundingBoxes(mesh, bbox)
      implicit none
      type(Mesh_t),                     intent(in)  :: mesh
      type(ElementBBox_t), allocatable, intent(out) :: bbox(:)
!
!     ---------------
!     Local variables
!     ---------------
!
      integer :: e, i, j, k

      allocate(bbox(mesh % no_of_elements))

      do e = 1, mesh % no_of_elements
         associate(elem => mesh % elements(e))
         bbox(e) % xmin =  huge(1.0_RP)
         bbox(e) % xmax = -huge(1.0_RP)
         do k = 0, elem % Nmesh(3)
            do j = 0, elem % Nmesh(2)
               do i = 0, elem % Nmesh(1)
                  bbox(e) % xmin = min(bbox(e) % xmin, elem % x(:,i,j,k))
                  bbox(e) % xmax = max(bbox(e) % xmax, elem % x(:,i,j,k))
               end do
            end do
         end do
         end associate
      end do

   end subroutine BuildBoundingBoxes

!
!/////////////////////////////////////////////////////////////////////////////////////
!
!     Find all probes in the mesh using Newton iteration with bounding box filter.
!
!/////////////////////////////////////////////////////////////////////////////////////
!
   subroutine FindAllProbesInMesh(mesh, bbox, probes, no_of_probes)
      implicit none
      type(Mesh_t),        intent(in)    :: mesh
      type(ElementBBox_t), intent(in)    :: bbox(:)
      type(ProbeData_t),   intent(inout) :: probes(:)
      integer,             intent(in)    :: no_of_probes
!
!     ---------------
!     Local variables
!     ---------------
!
      integer         :: pID, found_count, warn_count, processed_count, new_count, milestone
      integer(kind=8) :: t0, t1, count_rate

      found_count     = 0
      warn_count      = 0
      processed_count = 0
      milestone       = max(1, no_of_probes / 20)   ! 5% steps

      write(STD_OUT,'(30X,A,I0,A)') "-> Searching ", no_of_probes, " probes..."

      call system_clock(t0, count_rate)

      !$OMP PARALLEL DO SCHEDULE(DYNAMIC,100) REDUCTION(+:found_count) DEFAULT(SHARED) PRIVATE(new_count)
      do pID = 1, no_of_probes
         call FindProbeInMesh(mesh, bbox, probes(pID) % x, &
                              probes(pID) % eID, probes(pID) % xi, probes(pID) % found)
         if (probes(pID) % found) found_count = found_count + 1
         !$OMP ATOMIC CAPTURE
         processed_count = processed_count + 1
         new_count       = processed_count
         !$OMP END ATOMIC
         if (mod(new_count, milestone) .eq. 0) then
            !$OMP CRITICAL (progress_bar)
            write(STD_OUT,'(30X,A,I0,A,I0,A,F5.1,A)') &
               "   [", new_count, "/", no_of_probes, &
               "]  ", 100.0_RP * new_count / no_of_probes, "%"
            !$OMP END CRITICAL (progress_bar)
         end if
      end do
      !$OMP END PARALLEL DO

      call system_clock(t1)
      write(STD_OUT,'(30X,A)') "   [done]"
      write(STD_OUT,'(30X,A,F8.3,A)') "-> Parallel search time: ", &
         real(t1-t0,RP)/real(count_rate,RP), " s"

      ! Print warnings sequentially after the parallel search
      do pID = 1, no_of_probes
         if (.not. probes(pID) % found) then
            warn_count = warn_count + 1
            if (warn_count .le. MAX_WARN) then
               write(STD_OUT,'(30X,A,A,A,3(ES12.4,1X),A)') "-> WARNING: Probe '", &
                  trim(probes(pID) % name), "' at (", probes(pID) % x, ") not found in mesh."
            end if
         end if
      end do

      if (warn_count .gt. MAX_WARN) then
         write(STD_OUT,'(30X,A,I0,A)') "   (", warn_count - MAX_WARN, " additional not-found probes suppressed)"
      end if
      write(STD_OUT,'(30X,A,I0,A,I0,A)') "-> ", found_count, " of ", no_of_probes, " probes found in mesh."

   end subroutine FindAllProbesInMesh

!
!/////////////////////////////////////////////////////////////////////////////////////
!
!     Search all elements for x_probe, using bounding box pre-filter.
!
!/////////////////////////////////////////////////////////////////////////////////////
!
   subroutine FindProbeInMesh(mesh, bbox, x_probe, eID, xi, found)
      implicit none
      type(Mesh_t),        intent(in)  :: mesh
      type(ElementBBox_t), intent(in)  :: bbox(:)
      real(kind=RP),       intent(in)  :: x_probe(NDIM)
      integer,             intent(out) :: eID
      real(kind=RP),       intent(out) :: xi(NDIM)
      logical,             intent(out) :: found
!
!     ---------------
!     Local variables
!     ---------------
!
      integer :: e

      found = .false.
      eID   = 0

      do e = 1, mesh % no_of_elements
         if (any(x_probe < bbox(e) % xmin - BBOX_TOL)) cycle
         if (any(x_probe > bbox(e) % xmax + BBOX_TOL)) cycle
         call FindPointInElement(mesh % elements(e), mesh % is2D, x_probe, xi, found)
         if (found) then
            eID = e
            return
         end if
      end do

   end subroutine FindProbeInMesh

!
!/////////////////////////////////////////////////////////////////////////////////////
!
!     Newton iteration to find reference coordinates (xi,eta,zeta) of x_probe
!     inside element e. Returns found=.true. if the point is inside the element.
!
!/////////////////////////////////////////////////////////////////////////////////////
!
   subroutine FindPointInElement(e, is2D, x_probe, xi, found)
      implicit none
      type(Element_t), intent(in)  :: e
      logical,         intent(in)  :: is2D
      real(kind=RP),   intent(in)  :: x_probe(NDIM)
      real(kind=RP),   intent(out) :: xi(NDIM)
      logical,         intent(out) :: found
!
!     ---------------
!     Local variables
!     ---------------
!
      real(kind=RP) :: Lxi   (0:e % Nmesh(1))
      real(kind=RP) :: Leta  (0:e % Nmesh(2))
      real(kind=RP) :: Lzeta (0:e % Nmesh(3))
      real(kind=RP) :: dLxi  (0:e % Nmesh(1))
      real(kind=RP) :: dLeta (0:e % Nmesh(2))
      real(kind=RP) :: dLzeta(0:e % Nmesh(3))
      real(kind=RP) :: x_calc(NDIM), F(NDIM), Jac(NDIM,NDIM), dxi_step(NDIM)
      integer        :: iter, i, j, k

      xi    = 0.0_RP
      found = .false.

      do iter = 1, MAX_NEWTON

         Lxi   = spA(e % Nmesh(1)) % lj(xi(1))
         Leta  = spA(e % Nmesh(2)) % lj(xi(2))
         dLxi  = spA(e % Nmesh(1)) % dlj(xi(1))
         dLeta = spA(e % Nmesh(2)) % dlj(xi(2))

         if (is2D) then
            Lzeta(0)  = 1.0_RP
            dLzeta(0) = 0.0_RP
         else
            Lzeta  = spA(e % Nmesh(3)) % lj(xi(3))
            dLzeta = spA(e % Nmesh(3)) % dlj(xi(3))
         end if

         x_calc = 0.0_RP
         Jac    = 0.0_RP
         do k = 0, e % Nmesh(3)
            do j = 0, e % Nmesh(2)
               do i = 0, e % Nmesh(1)
                  x_calc   = x_calc   + e % x(:,i,j,k) * Lxi(i) * Leta(j)  * Lzeta(k)
                  Jac(:,1) = Jac(:,1) + e % x(:,i,j,k) * dLxi(i)* Leta(j)  * Lzeta(k)
                  Jac(:,2) = Jac(:,2) + e % x(:,i,j,k) * Lxi(i) * dLeta(j) * Lzeta(k)
                  Jac(:,3) = Jac(:,3) + e % x(:,i,j,k) * Lxi(i) * Leta(j)  * dLzeta(k)
               end do
            end do
         end do

         F = x_calc - x_probe
         if (norm2(F) < NEWTON_TOL) exit

         dxi_step = Solve3x3(Jac, F)
         xi = xi - dxi_step

         xi(1) = max(min(xi(1),  2.0_RP), -2.0_RP)
         xi(2) = max(min(xi(2),  2.0_RP), -2.0_RP)
         if (.not. is2D) xi(3) = max(min(xi(3), 2.0_RP), -2.0_RP)

      end do

      if (is2D) then
         found = (abs(xi(1)) <= 1.0_RP + INSIDE_TOL) .and. &
                 (abs(xi(2)) <= 1.0_RP + INSIDE_TOL)
      else
         found = all(abs(xi) <= 1.0_RP + INSIDE_TOL)
      end if

   end subroutine FindPointInElement

!
!/////////////////////////////////////////////////////////////////////////////////////
!
!     Interpolate conservative variables Q at reference coordinates xi.
!
!/////////////////////////////////////////////////////////////////////////////////////
!
   function InterpolateQAtProbe(e, xi, nodeType) result(Q_probe)
      implicit none
      type(Element_t), intent(in) :: e
      real(kind=RP),   intent(in) :: xi(NDIM)
      integer,         intent(in) :: nodeType
      real(kind=RP)               :: Q_probe(NVARS)
!
!     ---------------
!     Local variables
!     ---------------
!
      real(kind=RP) :: Lxi  (0:e % Nsol(1))
      real(kind=RP) :: Leta (0:e % Nsol(2))
      real(kind=RP) :: Lzeta(0:e % Nsol(3))
      integer        :: i, j, k

      Lxi  = spA(e % Nsol(1)) % lj(xi(1))
      Leta = spA(e % Nsol(2)) % lj(xi(2))

      if (e % Nsol(3) .eq. 0) then
         Lzeta(0) = 1.0_RP
      else
         Lzeta = spA(e % Nsol(3)) % lj(xi(3))
      end if

      Q_probe = 0.0_RP
      do k = 0, e % Nsol(3)
         do j = 0, e % Nsol(2)
            do i = 0, e % Nsol(1)
               Q_probe = Q_probe + e % Q(:,i,j,k) * Lxi(i) * Leta(j) * Lzeta(k)
            end do
         end do
      end do

   end function InterpolateQAtProbe

!
!/////////////////////////////////////////////////////////////////////////////////////
!
!     For each solution: interpolate Q at all probe locations, compute output
!     variables and write a Tecplot .tec file. The mock element is allocated
!     once and reused across probes.
!
!/////////////////////////////////////////////////////////////////////////////////////
!
   subroutine InterpolateAndWriteProbes(mesh, probes, no_of_probes, solutionName, probesFileName)
      use FileReadingUtilities, only: getFileName
      implicit none
      type(Mesh_t),      intent(inout) :: mesh
      type(ProbeData_t), intent(in)    :: probes(:)
      integer,           intent(in)    :: no_of_probes
      character(len=*),  intent(in)    :: solutionName
      character(len=*),  intent(in)    :: probesFileName
!
!     ---------------
!     Local variables
!     ---------------
!
      type(Element_t)                :: mock_e
      real(kind=RP),     allocatable :: output(:,:,:,:)
      real(kind=RP),     allocatable :: outVars(:,:)
      real(kind=RP),     allocatable :: Q_all(:,:)
      integer                        :: pID, fid, eID, found_count
      integer(kind=8)                :: t_total, t_end, count_rate
      character(len=LINE_LENGTH)     :: outputFile, formatout
      character(len=1024)            :: title

      write(STD_OUT,'(/)')
      call Section_Header("Probes extraction")
      call system_clock(t_total, count_rate)
      found_count = count(probes(1:no_of_probes) % found)
      write(STD_OUT,'(30X,A,I0,A,I0,A)') "-> Extracting ", found_count, &
         " of ", no_of_probes, " probes..."
!
!     Initialize spectral basis for solution orders
!     -----------------------------------------------
      do eID = 1, mesh % no_of_elements
         associate(e => mesh % elements(eID))
         call addNewSpectralBasis(spA, e % Nsol, mesh % nodeType)
         end associate
      end do
!
!     Set up output variable metadata
!     ---------------------------------
      call getOutputVariables()
      allocate(output (no_of_outputVariables, 0:0, 0:0, 0:0))
      allocate(outVars(no_of_outputVariables, no_of_probes))
      outVars = 0.0_RP
!
!     Allocate mock element once — reused for every probe
!     -----------------------------------------------------
      mock_e % Nout = 0
      allocate(mock_e % Qout    (NVARS, 0:0, 0:0, 0:0))
      allocate(mock_e % xOut    (NDIM,  0:0, 0:0, 0:0))
      allocate(mock_e % QDot_out(NVARS, 0:0, 0:0, 0:0))
      allocate(mock_e % U_xout  (NVARS, 0:0, 0:0, 0:0))
      allocate(mock_e % U_yout  (NVARS, 0:0, 0:0, 0:0))
      allocate(mock_e % U_zout  (NVARS, 0:0, 0:0, 0:0))
      allocate(mock_e % statsout (1,    0:0, 0:0, 0:0))
      allocate(mock_e % mu_NSout (1,    0:0, 0:0, 0:0))
      allocate(mock_e % ut_NSout (1,    0:0, 0:0, 0:0))
      allocate(mock_e % wallYout (1,    0:0, 0:0, 0:0))
      allocate(mock_e % mu_sgsout(1,    0:0, 0:0, 0:0))
      allocate(mock_e % wallY    (1,    0:0, 0:0, 0:0))
      allocate(mock_e % ut_NS    (1,    0:0, 0:0, 0:0))

      mock_e % QDot_out  = 0.0_RP
      mock_e % U_xout    = 0.0_RP
      mock_e % U_yout    = 0.0_RP
      mock_e % U_zout    = 0.0_RP
      mock_e % wallY     = 0.0_RP
      mock_e % ut_NS     = 0.0_RP
      mock_e % statsout  = 0.0_RP
      mock_e % mu_NSout  = 0.0_RP
      mock_e % ut_NSout  = 0.0_RP
      mock_e % wallYout  = 0.0_RP
      mock_e % mu_sgsout = 0.0_RP
      mock_e % sensor    = 0.0_RP
!
!     Phase 1: interpolate Q at all probe locations (parallel)
!     ----------------------------------------------------------
      allocate(Q_all(NVARS, no_of_probes))
      Q_all = 0.0_RP

      !$OMP PARALLEL DO SCHEDULE(DYNAMIC,256) DEFAULT(SHARED)
      do pID = 1, no_of_probes
         if (.not. probes(pID) % found) cycle
         Q_all(:,pID) = InterpolateQAtProbe(mesh % elements(probes(pID) % eID), &
                                             probes(pID) % xi, mesh % nodeType)
      end do
      !$OMP END PARALLEL DO
!
!     Phase 2: compute output variables (sequential — shared mock element)
!     ---------------------------------------------------------------------
      do pID = 1, no_of_probes
         if (.not. probes(pID) % found) cycle

         mock_e % Qout(:,0,0,0) = Q_all(:,pID)
         mock_e % xOut(:,0,0,0) = probes(pID) % x

         call ComputeOutputVariables(no_of_outputVariables, outputVariableNames, &
                                      mock_e % Nout, mock_e, output, mesh % refs, &
                                      .false., .false., .false.)
         outVars(:,pID) = output(:,0,0,0)
      end do

      deallocate(Q_all)
!
!     Deallocate mock element
!     ------------------------
      deallocate(mock_e % Qout, mock_e % xOut, mock_e % QDot_out)
      deallocate(mock_e % U_xout, mock_e % U_yout, mock_e % U_zout)
      deallocate(mock_e % statsout, mock_e % mu_NSout, mock_e % ut_NSout)
      deallocate(mock_e % wallYout, mock_e % mu_sgsout)
      deallocate(mock_e % wallY, mock_e % ut_NS)
      deallocate(output)
!
!     Write output .tec file
!     -----------------------
      outputFile = trim(getDirName(solutionName)) // trim(getBaseName(solutionName)) // &
                   "." // trim(getBaseName(probesFileName)) // ".tec"
      write(STD_OUT,'(30X,A,A)') "-> Output file: ", trim(outputFile)

      open(newunit=fid, file=trim(outputFile), action="write", status="unknown")

      write(title,'(A,A,A)') '"Probes extracted from ', trim(solutionName), '"'
      write(fid,'(A,A)') "TITLE = ", trim(title)
      write(fid,'(A,A)') 'VARIABLES = "x","y","z"', trim(getOutputVariablesLabel())
      write(fid,'(A,I0,A)') 'ZONE T="Probes", I=', no_of_probes, ', J=1, K=1, F=POINT'

      write(formatout,'(A,I0,A,A)') "(", 3 + no_of_outputVariables, PRECISION_FORMAT, ")"

      do pID = 1, no_of_probes
         write(fid, trim(formatout)) probes(pID) % x, outVars(:,pID)
      end do

      close(fid)
      call system_clock(t_end)
      write(STD_OUT,'(30X,A,I0,A,I0,A,F8.3,A)') "-> Written ", found_count, &
         " of ", no_of_probes, " probes in ", &
         real(t_end-t_total,8)/real(count_rate,8), " s"

      deallocate(outVars)

   end subroutine InterpolateAndWriteProbes

!
!/////////////////////////////////////////////////////////////////////////////////////
!
!     Solve a 3x3 linear system A*x = b using Cramer's rule.
!
!/////////////////////////////////////////////////////////////////////////////////////
!
   pure function Solve3x3(A, b) result(x)
      implicit none
      real(kind=RP), intent(in) :: A(3,3), b(3)
      real(kind=RP)             :: x(3)
!
!     ---------------
!     Local variables
!     ---------------
!
      real(kind=RP) :: det, inv_det

      det = A(1,1)*(A(2,2)*A(3,3) - A(2,3)*A(3,2)) &
          - A(1,2)*(A(2,1)*A(3,3) - A(2,3)*A(3,1)) &
          + A(1,3)*(A(2,1)*A(3,2) - A(2,2)*A(3,1))

      if (abs(det) < 1.0e-30_RP) then
         x = 0.0_RP
         return
      end if

      inv_det = 1.0_RP / det

      x(1) = inv_det * ( b(1)*(A(2,2)*A(3,3) - A(2,3)*A(3,2)) &
                       - A(1,2)*(b(2)*A(3,3) - A(2,3)*b(3))   &
                       + A(1,3)*(b(2)*A(3,2) - A(2,2)*b(3)) )

      x(2) = inv_det * ( A(1,1)*(b(2)*A(3,3) - A(2,3)*b(3))   &
                       - b(1)*(A(2,1)*A(3,3) - A(2,3)*A(3,1)) &
                       + A(1,3)*(A(2,1)*b(3) - b(2)*A(3,1)) )

      x(3) = inv_det * ( A(1,1)*(A(2,2)*b(3) - b(2)*A(3,2))   &
                       - A(1,2)*(A(2,1)*b(3) - b(2)*A(3,1))   &
                       + b(1)*(A(2,1)*A(3,2) - A(2,2)*A(3,1)) )

   end function Solve3x3 

end module Horses2probesModule
