#include "Includes.h"
module Stats2PltModule
   use SMConstants
   use SolutionFile
   use Headers
   use InterpolationMatrices
   use FileReadingUtilities      , only: getFileName
   use Solution2PltModule        , only: WriteBoundaryToTecplot
   implicit none

   private
   public   Stats2Plt

#define PRECISION_FORMAT "(E13.5)"

   integer, parameter :: NSTATS_OUTVARS = 9
   integer, parameter :: NFAVRE_OUTVARS = 6
   ! Canonical output names for Reynolds and Favre stats variables (order matters)
   character(len=8), parameter :: STATS_OUT_NAMES(NSTATS_OUTVARS) = &
      ["Umean   ","Vmean   ","Wmean   ","Rxx     ","Ryy     ","Rzz     ","Rxy     ","Rxz     ","Ryz     "]
   character(len=4), parameter :: FAVRE_OUT_NAMES(NFAVRE_OUTVARS) = &
      ["Fxx ","Fyy ","Fzz ","Fxy ","Fxz ","Fyz "]
   ! Per-run output filter; set by buildStatsFilter() before each file write
   logical :: statsVarInclude(NSTATS_OUTVARS) = .true.
   logical :: favreVarInclude(NFAVRE_OUTVARS) = .true.

   contains
      subroutine Stats2Plt(meshName, solutionName, fixedOrder, basis, Nout, mode)
         use getTask
         implicit none
         character(len=*), intent(in)     :: meshName
         character(len=*), intent(in)     :: solutionName
         integer,          intent(in)     :: basis
         logical,          intent(in)     :: fixedOrder
         integer,          intent(in)     :: Nout(3)
         integer,          intent(in)     :: mode

         write(STD_OUT,'(/)')
         call SubSection_Header("Job description")

         select case (mode)
         case(MODE_FINITEELM)
            write(STD_OUT,'(30X,A3,A)') "->", " Output mode: Tecplot FE"
         case(MODE_MULTIZONE)
            write(STD_OUT,'(30X,A3,A)') "->", " Output mode: Tecplot Multi-Zone"
         end select

         select case ( basis )

         case(EXPORT_GAUSS)

            if ( fixedOrder ) then
               write(STD_OUT,'(30X,A3,A)') "->", " Export to Gauss points with fixed order"
               write(STD_OUT,'(30X,A,A30,I0,A,I0,A,I0,A)') "->" , "Output order: [",&
                                                Nout(1),",",Nout(2),",",Nout(3),"]."
               call Stats2Plt_GaussPoints_FixedOrder(meshName, solutionName, Nout, mode)

            else
               write(STD_OUT,'(30X,A3,A)') "->", " Export to Gauss points"
               call Stats2Plt_GaussPoints(meshName, solutionName, mode)

            end if

         case(EXPORT_HOMOGENEOUS)

            write(STD_OUT,'(30X,A3,A)') "->", " Export to homogeneous points"
            write(STD_OUT,'(30X,A,A30,I0,A,I0,A,I0,A)') "->" , "Output order: [",&
                                        Nout(1),",",Nout(2),",",Nout(3),"]."
            call Stats2Plt_Homogeneous(meshName, solutionName, Nout, mode)

         end select

      end subroutine Stats2Plt
!
!//////////////////////////////////////////////////////////////////////////////////////////
!
!     Gauss Points procedures
!     -----------------------
!
!//////////////////////////////////////////////////////////////////////////////////////////
!
      subroutine Stats2Plt_GaussPoints(meshName, solutionName, mode)
         use Storage
         use NodalStorageClass
         use SharedSpectralBasis
         use OutputVariables
         use getTask,          only: MODE_FINITEELM
         implicit none
         character(len=*), intent(in)     :: meshName
         character(len=*), intent(in)     :: solutionName
         integer,          intent(in)     :: mode
!
!        ---------------
!        Local variables
!        ---------------
!
         type(Mesh_t)                    :: mesh
         character(len=LINE_LENGTH)      :: meshPltName
         character(len=LINE_LENGTH)      :: solutionFile
         character(len=1024)             :: title
         integer                         :: no_of_elements, eID
         integer                         :: fid, bID
         integer                         :: Nmesh(4), Nsol(4)
!
!        Read the mesh and solution data
!        -------------------------------
         call mesh % ReadMesh(meshName)
         call mesh % ReadSolution(SolutionName)
         no_of_elements = mesh % no_of_elements
!
!        Transform zones to the output variables
!        ---------------------------------------
         do eID = 1, no_of_elements
            associate ( e => mesh % elements(eID) )
            e % Nout = e % Nsol
!
!           Construct spectral basis
!           ------------------------
            call addNewSpectralBasis(spA, e % Nmesh, mesh % nodeType)
            call addNewSpectralBasis(spA, e % Nsol, mesh % nodeType)
!
!           Project mesh and solution
!           -------------------------
            call ProjectStorageGaussPoints(e, spA, e % Nmesh, e % Nsol)

            end associate
         end do
!
!        Write the solution file name
!        ----------------------------
         solutionFile = trim(getFileName(solutionName)) // ".tec"
!
!        Create the file
!        ---------------
         open(newunit = fid, file = trim(solutionFile), action = "write", status = "unknown")
!
!        Add the title
!        -------------
         write(title,'(A,A,A,A,A)') '"Generated from ',trim(meshName),' and ',trim(solutionName),'"'
         write(fid,'(A,A)') "TITLE = ", trim(title)
!
!        Add the variables (filtered by output variables if set)
!        -------------------------------------------------------
         call buildStatsFilter()
         call printStatsOutputVariables()
         write(fid,'(A)') trim(buildVarsHeader())
!
!        Write each element zone
!        -----------------------
         if ( mode == MODE_FINITEELM ) then
            call WriteSingleFluidZoneToTecplotStats(fid, mesh)
         else
            do eID = 1, no_of_elements
               associate ( e => mesh % elements(eID) )
!
!              Write the tecplot file
!              ----------------------
               call WriteElementToTecplot(fid, e, mesh % refs)
               end associate
            end do
         end if
!
!        Write boundaries
!        ----------------
         if (hasBoundaries) then
            if ( mode == MODE_FINITEELM ) then
               do bID=1, size (mesh % boundaries)
                  call WriteSingleBoundaryZoneToTecplotStats(fid, mesh % boundaries(bID), mesh % elements)
               end do
            else
               do bID=1, size (mesh % boundaries)
                  call WriteBoundaryToTecplot(fid, mesh % boundaries(bID), mesh % elements)
               end do
            end if
         end if
!
!        Close the file
!        --------------
         close(fid)

      end subroutine Stats2Plt_GaussPoints

      subroutine ProjectStorageGaussPoints(e, spA, N1, N2)
         use Storage
         use NodalStorageClass
         use ProlongMeshAndSolution
         implicit none
         type(Element_t)     :: e
         type(NodalStorage_t), intent(in) :: spA(0:)
         integer           , intent(in) :: N1(3)
         integer           , intent(in) :: N2(3)

         e % Nout = e % Nsol
         if ( all(e % Nmesh .eq. e % Nout) ) then
            e % xOut(1:,0:,0:,0:) => e % x

         else
            allocate( e % xOut(1:3,0:e % Nout(1), 0:e % Nout(2), 0:e % Nout(3)) )
            call prolongMeshToGaussPoints(e, spA, N1, N2)

         end if

         if (NSTAT .gt. 0) e % statsout(1:,0:,0:,0:) => e % stats
         if (statsHasFavre) e % favreout(1:,0:,0:,0:) => e % favre

      end subroutine ProjectStorageGaussPoints
!
!//////////////////////////////////////////////////////////////////////////////////
!
!     Gauss points with fixed order procedures
!     ----------------------------------------
!
!//////////////////////////////////////////////////////////////////////////////////
!
      subroutine Stats2Plt_GaussPoints_FixedOrder(meshName, solutionName, Nout, mode)
         use Storage
         use NodalStorageClass
         use SharedSpectralBasis
         use OutputVariables
         use getTask,          only: MODE_FINITEELM
         implicit none
         character(len=*), intent(in)     :: meshName
         character(len=*), intent(in)     :: solutionName
         integer,          intent(in)     :: Nout(3)
         integer,          intent(in)     :: mode
!
!        ---------------
!        Local variables
!        ---------------
!
         type(Mesh_t)                               :: mesh
         character(len=LINE_LENGTH)                 :: meshPltName
         character(len=LINE_LENGTH)                 :: solutionFile
         character(len=1024)                        :: title
         integer                                    :: no_of_elements, eID
         integer                                    :: fid, bID
!
!        Read the mesh and solution data
!        -------------------------------
         call mesh % ReadMesh(meshName)
         call mesh % ReadSolution(SolutionName)
!
!        Allocate the output spectral basis
!        ----------------------------------
         call spA(Nout(1)) % Construct(GAUSS, Nout(1))
         call spA(Nout(2)) % Construct(GAUSS, Nout(2))
         call spA(Nout(3)) % Construct(GAUSS, Nout(3))
!
!        Write each element zone
!        -----------------------
         do eID = 1, mesh % no_of_elements
            associate ( e => mesh % elements(eID) )
            e % Nout = Nout
!
!           Construct spectral basis
!           ------------------------
            call addNewSpectralBasis(spA, e % Nmesh, mesh % nodeType)
            call addNewSpectralBasis(spA, e % Nsol , mesh % nodeType)
!
!           Construct interpolation matrices
!           --------------------------------
            associate( spAoutXi   => spA(Nout(1)), &
                       spAoutEta  => spA(Nout(2)), &
                       spAoutZeta => spA(Nout(3)) )
            call addNewInterpolationMatrix(Tset, e % Nsol(1), spA(e % Nsol(1)), e % Nout(1), spAoutXi   % x)
            call addNewInterpolationMatrix(Tset, e % Nsol(2), spA(e % Nsol(2)), e % Nout(2), spAoutEta  % x)
            call addNewInterpolationMatrix(Tset, e % Nsol(3), spA(e % Nsol(3)), e % Nout(3), spAoutZeta % x)
            end associate
!
!           Perform interpolation
!           ---------------------
            call ProjectStorageGaussPoints_FixedOrder(e, spA, e % Nmesh, e % Nsol, e % Nout, &
                                                                    Tset(e % Nout(1), e % Nsol(1)) % T, &
                                                                    Tset(e % Nout(2), e % Nsol(2)) % T, &
                                                                    Tset(e % Nout(3), e % Nsol(3)) % T    )

            end associate
         end do
!
!        Write the solution file name
!        ----------------------------
         solutionFile = trim(getFileName(solutionName)) // ".tec"
!
!        Create the file
!        ---------------
         open(newunit = fid, file = trim(solutionFile), action = "write", status = "unknown")
!
!        Add the title
!        -------------
         write(title,'(A,A,A,A,A)') '"Generated from ',trim(meshName),' and ',trim(solutionName),'"'
         write(fid,'(A,A)') "TITLE = ", trim(title)
!
!        Add the variables (filtered by output variables if set)
!        -------------------------------------------------------
         call buildStatsFilter()
         call printStatsOutputVariables()
         write(fid,'(A)') trim(buildVarsHeader())
!
!        Write elements
!        --------------
         if ( mode == MODE_FINITEELM ) then
            call WriteSingleFluidZoneToTecplotStats(fid, mesh)
         else
            do eID = 1, mesh % no_of_elements
               associate ( e => mesh % elements(eID) )

               call WriteElementToTecplot(fid, e, mesh % refs)
               end associate
            end do
         end if
!
!        Write boundaries
!        ----------------
         if (hasBoundaries) then
            if ( mode == MODE_FINITEELM ) then
               do bID=1, size (mesh % boundaries)
                  call WriteSingleBoundaryZoneToTecplotStats(fid, mesh % boundaries(bID), mesh % elements)
               end do
            else
               do bID=1, size (mesh % boundaries)
                  call WriteBoundaryToTecplot(fid, mesh % boundaries(bID), mesh % elements)
               end do
            end if
         end if

!
!        Close the file
!        --------------
         close(fid)

      end subroutine Stats2Plt_GaussPoints_FixedOrder

      subroutine ProjectStorageGaussPoints_FixedOrder(e, spA, NM, NS, Nout, Tx, Ty, Tz)
         use Storage
         use NodalStorageClass
         use ProlongMeshAndSolution
         implicit none
         type(Element_t)     :: e
         type(NodalStorage_t),  intent(in)  :: spA(0:)
         integer           ,  intent(in)  :: NM(3)
         integer           ,  intent(in)  :: NS(3)
         integer           ,  intent(in)  :: Nout(3)
         real(kind=RP),       intent(in)  :: Tx(0:e % Nout(1), 0:e % Nsol(1))
         real(kind=RP),       intent(in)  :: Ty(0:e % Nout(2), 0:e % Nsol(2))
         real(kind=RP),       intent(in)  :: Tz(0:e % Nout(3), 0:e % Nsol(3))
!
!        Project mesh
!        ------------
         if ( all(e % Nmesh .eq. e % Nout) ) then
            e % xOut(1:,0:,0:,0:) => e % x

         else
            allocate( e % xOut(1:3,0:e % Nout(1), 0:e % Nout(2), 0:e % Nout(3)) )
            call prolongMeshToGaussPoints(e, spA, NM, Nout)

         end if
!
!        Project the solution
!        --------------------
         if ( all( e % Nsol .eq. e % Nout ) ) then
            if (NSTAT .gt. 0) e % statsout(1:,0:,0:,0:) => e % stats
            if (statsHasFavre) e % favreout(1:,0:,0:,0:) => e % favre

         else
            if (NSTAT .gt. 0) then
               allocate( e % statsout(1:NSTAT,0:e % Nout(1), 0:e % Nout(2), 0:e % Nout(3)) )
               call prolongSolutionToGaussPoints(NSTAT, e % Nsol, e % stats, e % Nout, e % statsout, Tx, Ty, Tz)
            end if
            if (statsHasFavre) then
               allocate( e % favreout(1:NFAVRE_VARS,0:e % Nout(1), 0:e % Nout(2), 0:e % Nout(3)) )
               call prolongSolutionToGaussPoints(NFAVRE_VARS, e % Nsol, e % favre, e % Nout, e % favreout, Tx, Ty, Tz)
            end if

         end if

      end subroutine ProjectStorageGaussPoints_FixedOrder
!
!////////////////////////////////////////////////////////////////////////////
!
!     Homogeneous procedures
!     ----------------------
!
!////////////////////////////////////////////////////////////////////////////
!
      subroutine Stats2Plt_Homogeneous(meshName, solutionName, Nout, mode)
         use Storage
         use NodalStorageClass
         use SharedSpectralBasis
         use OutputVariables
         use getTask,          only: MODE_FINITEELM
         implicit none
         character(len=*), intent(in)     :: meshName
         character(len=*), intent(in)     :: solutionName
         integer,          intent(in)     :: Nout(3)
         integer,          intent(in)     :: mode
!
!        ---------------
!        Local variables
!        ---------------
!
         type(Mesh_t)                               :: mesh
         character(len=LINE_LENGTH)                 :: meshPltName
         character(len=LINE_LENGTH)                 :: solutionFile
         character(len=1024)                        :: title
         integer                                    :: no_of_elements, eID
         integer                                    :: fid, bID
         real(kind=RP)                              :: xi(0:Nout(1)), eta(0:Nout(2)), zeta(0:Nout(3))
         integer                                    :: i
!
!        Read the mesh and solution data
!        -------------------------------
         call mesh % ReadMesh(meshName)
         call mesh % ReadSolution(SolutionName)
!
!        Set homogeneous nodes
!        ---------------------
         xi   = RESHAPE( (/ (-1.0_RP + 2.0_RP*i/Nout(1),i=0,Nout(1)) /), (/ Nout(1)+1 /) )
         eta  = RESHAPE( (/ (-1.0_RP + 2.0_RP*i/Nout(2),i=0,Nout(2)) /), (/ Nout(2)+1 /) )
         zeta = RESHAPE( (/ (-1.0_RP + 2.0_RP*i/Nout(3),i=0,Nout(3)) /), (/ Nout(3)+1 /) )
!
!        Write each element zone
!        -----------------------
         do eID = 1, mesh % no_of_elements
            associate ( e => mesh % elements(eID) )
            e % Nout = Nout
!
!           Construct spectral basis for both mesh and solution
!           ---------------------------------------------------
            call addNewSpectralBasis(spA, e % Nmesh, mesh % nodeType)
            call addNewSpectralBasis(spA, e % Nsol , mesh % nodeType)
!
!           Construct interpolation matrices for the mesh
!           ---------------------------------------------
            call addNewInterpolationMatrix(Tset, e % Nmesh(1), spA(e % Nmesh(1)), e % Nout(1), xi)
            call addNewInterpolationMatrix(Tset, e % Nmesh(2), spA(e % Nmesh(2)), e % Nout(2), eta)   ! TODO: check why it was Nmesh(1)
            call addNewInterpolationMatrix(Tset, e % Nmesh(3), spA(e % Nmesh(3)), e % Nout(3), zeta)  ! TODO: check why it was Nmesh(1)

!
!           Construct interpolation matrices for the solution
!           -------------------------------------------------
            call addNewInterpolationMatrix(Tset, e % Nsol(1), spA(e % Nsol(1)), e % Nout(1), xi)
            call addNewInterpolationMatrix(Tset, e % Nsol(2), spA(e % Nsol(2)), e % Nout(2), eta)     ! TODO: check why it was Nsol(1)
            call addNewInterpolationMatrix(Tset, e % Nsol(3), spA(e % Nsol(3)), e % Nout(3), zeta)    ! TODO: check why it was Nsol(1)
!
!           Perform interpolation
!           ---------------------
            call ProjectStorageHomogeneousPoints(e, Tset(e % Nout(1), e % Nmesh(1)) % T, &
                                                    Tset(e % Nout(2), e % Nmesh(2)) % T, &
                                                    Tset(e % Nout(3), e % Nmesh(3)) % T, &
                                                     Tset(e % Nout(1), e % Nsol(1)) % T, &
                                                     Tset(e % Nout(2), e % Nsol(2)) % T, &
                                                     Tset(e % Nout(3), e % Nsol(3)) % T    )


            end associate
         end do
!
!        Write the solution file name
!        ----------------------------
         solutionFile = trim(getFileName(solutionName)) // ".tec"
!
!        Create the file
!        ---------------
         open(newunit = fid, file = trim(solutionFile), action = "write", status = "unknown")
!
!        Add the title
!        -------------
         write(title,'(A,A,A,A,A)') '"Generated from ',trim(meshName),' and ',trim(solutionName),'"'
         write(fid,'(A,A)') "TITLE = ", trim(title)
!
!        Add the variables (filtered by output variables if set)
!        -------------------------------------------------------
         call buildStatsFilter()
         call printStatsOutputVariables()
         write(fid,'(A)') trim(buildVarsHeader())
!
!        Write elements
!        --------------
         if ( mode == MODE_FINITEELM ) then
            call WriteSingleFluidZoneToTecplotStats(fid, mesh)
         else
            do eID = 1, mesh % no_of_elements
               associate ( e => mesh % elements(eID) )

               call WriteElementToTecplot(fid, e, mesh % refs)
               end associate
            end do
         end if
!
!        Write boundaries
!        ----------------
         if (hasBoundaries) then
            if ( mode == MODE_FINITEELM ) then
               do bID=1, size (mesh % boundaries)
                  call WriteSingleBoundaryZoneToTecplotStats(fid, mesh % boundaries(bID), mesh % elements)
               end do
            else
               do bID=1, size (mesh % boundaries)
                  call WriteBoundaryToTecplot(fid, mesh % boundaries(bID), mesh % elements)
               end do
            end if
         end if
!
!        Close the file
!        --------------
         close(fid)

      end subroutine Stats2Plt_Homogeneous

      subroutine ProjectStorageHomogeneousPoints(e, TxMesh, TyMesh, TzMesh, TxSol, TySol, TzSol)
         use Storage
         use NodalStorageClass
         implicit none
         type(Element_t)     :: e
         real(kind=RP),       intent(in)  :: TxMesh(0:e % Nout(1), 0:e % Nmesh(1))
         real(kind=RP),       intent(in)  :: TyMesh(0:e % Nout(2), 0:e % Nmesh(2))
         real(kind=RP),       intent(in)  :: TzMesh(0:e % Nout(3), 0:e % Nmesh(3))
         real(kind=RP),       intent(in)  :: TxSol(0:e % Nout(1), 0:e % Nsol(1))
         real(kind=RP),       intent(in)  :: TySol(0:e % Nout(2), 0:e % Nsol(2))
         real(kind=RP),       intent(in)  :: TzSol(0:e % Nout(3), 0:e % Nsol(3))
!
!        ---------------
!        Local variables
!        ---------------
!
         integer     :: i, j, k, iVar, l, m, n
!
!        Project mesh
!        ------------
         allocate( e % xOut(1:3,0:e % Nout(1), 0:e % Nout(2), 0:e % Nout(3)) )
         e % xOut = 0.0_RP

         do n = 0, e % Nmesh(3) ; do m = 0, e % Nmesh(2) ; do l = 0, e % Nmesh(1)
            do k = 0, e % Nout(3) ; do j = 0, e % Nout(2) ; do i = 0, e % Nout(1)
               e % xOut(:,i,j,k) = e % xOut(:,i,j,k) + e % x(:,l,m,n) * TxMesh(i,l) * TyMesh(j,m) * TzMesh(k,n)
            end do            ; end do            ; end do
         end do            ; end do            ; end do

!
!        Project the solution
!        --------------------
         if (NSTAT .gt. 0) then
            allocate( e % statsout(1:NSTAT,0:e % Nout(1), 0:e % Nout(2), 0:e % Nout(3)) )
            e % statsout = 0.0_RP
            do n = 0, e % Nsol(3) ; do m = 0, e % Nsol(2) ; do l = 0, e % Nsol(1)
               do k = 0, e % Nout(3) ; do j = 0, e % Nout(2) ; do i = 0, e % Nout(1)
                  e % statsout(:,i,j,k) = e % statsout(:,i,j,k) + e % stats(:,l,m,n) * TxSol(i,l) * TySol(j,m) * TzSol(k,n)
               end do            ; end do            ; end do
            end do            ; end do            ; end do
         end if
         if (statsHasFavre) then
            allocate( e % favreout(1:NFAVRE_VARS,0:e % Nout(1), 0:e % Nout(2), 0:e % Nout(3)) )
            e % favreout = 0.0_RP
            do n = 0, e % Nsol(3) ; do m = 0, e % Nsol(2) ; do l = 0, e % Nsol(1)
               do k = 0, e % Nout(3) ; do j = 0, e % Nout(2) ; do i = 0, e % Nout(1)
                  e % favreout(:,i,j,k) = e % favreout(:,i,j,k) + e % favre(:,l,m,n) * TxSol(i,l) * TySol(j,m) * TzSol(k,n)
               end do            ; end do            ; end do
            end do            ; end do            ; end do
         end if

      end subroutine ProjectStorageHomogeneousPoints
!
!/////////////////////////////////////////////////////////////////////////////
!
!     Write solution
!     --------------
!
!/////////////////////////////////////////////////////////////////////////////
!
      integer function countStatsOutputVars() result(nout_vars)
         use Storage, only: NSTAT, statsHasFavre
         implicit none
         integer :: var, vi

         vi = 0
         do var = 1, NSTATS_OUTVARS
            if (NSTAT .gt. 0 .and. statsVarInclude(var)) vi = vi + 1
         end do
         do var = 1, NFAVRE_OUTVARS
            if (statsHasFavre .and. favreVarInclude(var)) vi = vi + 1
         end do
         nout_vars = vi

      end function countStatsOutputVars

      subroutine ComputeStatsOutputVars(e, nout_vars)
         use Storage
         use StatisticsMonitor
         implicit none
         type(Element_t),    intent(inout) :: e
         integer,            intent(out)   :: nout_vars
!
!        ---------------
!        Local variables
!        ---------------
!
         integer                    :: i, j, k, var, vi
         integer                    :: statsIdx(NSTATS_OUTVARS), favreIdx(NFAVRE_OUTVARS)
         real(kind=RP)              :: stats9(NSTATS_OUTVARS)
!
!        Precompute output column indices (0 = not selected)
!        ---------------------------------------------------
         vi = 0
         do var = 1, NSTATS_OUTVARS
            if (NSTAT .gt. 0 .and. statsVarInclude(var)) then
               vi = vi + 1 ; statsIdx(var) = vi
            else
               statsIdx(var) = 0
            end if
         end do
         do var = 1, NFAVRE_OUTVARS
            if (statsHasFavre .and. favreVarInclude(var)) then
               vi = vi + 1 ; favreIdx(var) = vi
            else
               favreIdx(var) = 0
            end if
         end do
         nout_vars = vi

         allocate(e % outputVars(1:nout_vars, 0:e % Nout(1), 0:e % Nout(2), 0:e % Nout(3)))

         do k = 0, e % Nout(3) ; do j = 0, e % Nout(2) ; do i = 0, e % Nout(1)

            if (NSTAT .gt. 0) then
               stats9(1) = e % statsout(U, i,j,k)
               stats9(2) = e % statsout(V, i,j,k)
               stats9(3) = e % statsout(W, i,j,k)
               stats9(4) = e % statsout(UU,i,j,k) - POW2(e % statsout(U,i,j,k))
               stats9(5) = e % statsout(VV,i,j,k) - POW2(e % statsout(V,i,j,k))
               stats9(6) = e % statsout(WW,i,j,k) - POW2(e % statsout(W,i,j,k))
               stats9(7) = e % statsout(UV,i,j,k) - e % statsout(U,i,j,k)*e % statsout(V,i,j,k)
               stats9(8) = e % statsout(UW,i,j,k) - e % statsout(U,i,j,k)*e % statsout(W,i,j,k)
               stats9(9) = e % statsout(VW,i,j,k) - e % statsout(V,i,j,k)*e % statsout(W,i,j,k)
               do var = 1, NSTATS_OUTVARS
                  if (statsIdx(var) .gt. 0) &
                     e % outputVars(statsIdx(var), i,j,k) = stats9(var)
               end do
            end if

            if (statsHasFavre) then
               do var = 1, NFAVRE_OUTVARS
                  if (favreIdx(var) .gt. 0) &
                     e % outputVars(favreIdx(var), i,j,k) = e % favreout(var, i,j,k)
               end do
            end if

         end do ; end do ; end do

      end subroutine ComputeStatsOutputVars

      subroutine WriteElementToTecplot(fid, e, refs)
         use Storage
         implicit none
         integer,            intent(in)    :: fid
         type(Element_t),    intent(inout) :: e
         real(kind=RP),      intent(in)    :: refs(NO_OF_SAVED_REFS)
!
!        ---------------
!        Local variables
!        ---------------
!
         integer                    :: i, j, k, nout_vars
         character(len=LINE_LENGTH) :: formatout

         call ComputeStatsOutputVars(e, nout_vars)
!
!        Write zone header and data
!        --------------------------
         write(fid,'(A,I0,A,I0,A,I0,A)') "ZONE I=",e % Nout(1)+1,", J=",e % Nout(2)+1, &
                                            ", K=",e % Nout(3)+1,", F=POINT"

         formatout = getFormat(3 + nout_vars)

         do k = 0, e % Nout(3)   ; do j = 0, e % Nout(2)    ; do i = 0, e % Nout(1)
            write(fid,trim(formatout)) e % xOut(:,i,j,k), e % outputVars(:,i,j,k)
         end do               ; end do                ; end do

      end subroutine WriteElementToTecplot
!
!/////////////////////////////////////////////////////////////////////////////
!
!     Writes a single fluid/boundary zone using the FE Tecplot format
!     -> This format is more efficiently read by paraview and tecplot.
!
!/////////////////////////////////////////////////////////////////////////////
!
      subroutine WriteSingleFluidZoneToTecplotStats(fid, mesh)
         use Storage
         implicit none
         integer,      intent(in)    :: fid
         type(Mesh_t), intent(inout) :: mesh
!
!        ---------------
!        Local variables
!        ---------------
!
         integer :: numOfPoints, numOfFElems
         integer :: firstPoint(size(mesh % elements))
         integer :: eID, i, j, k, nout_vars
         integer :: corners(8), cornersFace(4)
         character(len=LINE_LENGTH) :: formatout

         nout_vars = countStatsOutputVars()
         formatout = getFormat(3 + nout_vars)
!
!        Count points and elements
!        -------------------------
         numOfPoints = product(mesh % elements(1) % Nout + 1)
         if (mesh % isSurface) then
            numOfFElems = product(mesh % elements(1) % Nout(1:2))
         else
            numOfFElems = product(mesh % elements(1) % Nout)
         end if
         firstPoint(1) = 1
         do eID = 2, size(mesh % elements)
            associate ( e => mesh % elements(eID) )
            firstPoint(eID) = numOfPoints + 1
            numOfPoints = numOfPoints + product(e % Nout + 1)
            if (mesh % isSurface) then
               numOfFElems = numOfFElems + product(e % Nout(1:2))
            else
               numOfFElems = numOfFElems + product(e % Nout)
            end if
            end associate
         end do

         if (mesh % isSurface) then
            write(fid,'(A,I0,A,I0,A)') 'ZONE T="FLUID" N=',numOfPoints,' E=',numOfFElems,' ET=QUADRILATERAL, F=FEPOINT'
         else
            write(fid,'(A,I0,A,I0,A)') 'ZONE T="FLUID" N=',numOfPoints,' E=',numOfFElems,' ET=BRICK, F=FEPOINT'
         end if
!
!        Write the points
!        ----------------
         do eID = 1, size(mesh % elements)
            associate ( e => mesh % elements(eID) )
            call ComputeStatsOutputVars(e, nout_vars)

            do k = 0, e % Nout(3) ; do j = 0, e % Nout(2) ; do i = 0, e % Nout(1)
               write(fid,trim(formatout)) e % xOut(:,i,j,k), e % outputVars(:,i,j,k)
            end do                ; end do                ; end do
            end associate
         end do
!
!        Write the elems connectivity
!        ----------------------------
         if (mesh % isSurface) then
            do eID = 1, size(mesh % elements)
               associate ( e => mesh % elements(eID) )

               do j = 0, e % Nout(2) - 1 ; do i = 0, e % Nout(1) - 1
                  cornersFace =  [ ij2localDOFStats(i,j,e%Nout(1:2)), ij2localDOFStats(i+1,j,e%Nout(1:2)), &
                                    ij2localDOFStats(i+1,j+1,e%Nout(1:2)), ij2localDOFStats(i,j+1,e%Nout(1:2)) ] + firstPoint(eID)
                  write(fid,*) cornersFace
               end do                  ; end do

               end associate
            end do
         else
            do eID = 1, size(mesh % elements)
               associate ( e => mesh % elements(eID) )

               do k = 0, e % Nout(3) - 1 ; do j = 0, e % Nout(2) - 1 ; do i = 0, e % Nout(1) - 1
                  corners =  [ ijk2localDOFStats(i,j,k  ,e%Nout), ijk2localDOFStats(i+1,j,k  ,e%Nout), &
                               ijk2localDOFStats(i+1,j+1,k  ,e%Nout), ijk2localDOFStats(i,j+1,k  ,e%Nout), &
                               ijk2localDOFStats(i,j,k+1,e%Nout), ijk2localDOFStats(i+1,j,k+1,e%Nout), &
                               ijk2localDOFStats(i+1,j+1,k+1,e%Nout), ijk2localDOFStats(i,j+1,k+1,e%Nout)  ] + firstPoint(eID)
                  write(fid,*) corners
               end do                    ; end do                    ; end do

               end associate
            end do
         end if

      end subroutine WriteSingleFluidZoneToTecplotStats

      function ijk2localDOFStats(i,j,k,Nout) result(idx)
         implicit none
         integer, intent(in)   :: i, j, k, Nout(3)
         integer               :: idx

         IF (i < 0 .OR. i > Nout(1))     error stop 'error in ijk2local, i has wrong value'
         IF (j < 0 .OR. j > Nout(2))     error stop 'error in ijk2local, j has wrong value'
         IF (k < 0 .OR. k > Nout(3))     error stop 'error in ijk2local, k has wrong value'

         idx = k*(Nout(1)+1)*(Nout(2)+1) + j*(Nout(1)+1) + i
      end function ijk2localDOFStats

      function ij2localDOFStats(i,j,Nout) result(idx)
         implicit none
         integer, intent(in)   :: i, j, Nout(2)
         integer               :: idx

         IF (i < 0 .OR. i > Nout(1))     error stop 'error in ijk2local, i has wrong value'
         IF (j < 0 .OR. j > Nout(2))     error stop 'error in ijk2local, j has wrong value'

         idx = j*(Nout(1)+1) + i
      end function ij2localDOFStats

      subroutine WriteSingleBoundaryZoneToTecplotStats(fd, boundary, elements)
         use Storage
         implicit none
         !-arguments-------------------------------------------
         integer         , intent(in) :: fd
         type(Boundary_t), intent(in) :: boundary
         type(Element_t) , intent(in) :: elements(:)
         !-local-variables-------------------------------------
         integer :: numOfPoints, numOfFElems
         integer :: fID, side
         integer :: corners(4)
         integer :: i,j,k
         integer :: N(3)
         integer :: firstPoint(boundary % no_of_faces)
         integer :: Nf      (2,boundary % no_of_faces)
         character(len=LINE_LENGTH) :: formatout
         integer :: nout_vars
         !-----------------------------------------------------

         nout_vars = countStatsOutputVars()
         formatout = getFormat(3 + nout_vars)
!
!        Count points and elements
!        -------------------------
         numOfPoints = 0
         numOfFElems = 0

         do fID = 1, boundary % no_of_faces
            associate (e => elements( boundary % elements(fID) ))
            side = boundary % elementSides(fID)

            select case (side)
               case(1,2) ; Nf(:,fID) = [e % Nout(1), e % Nout(3)]
               case(3,5) ; Nf(:,fID) = [e % Nout(1), e % Nout(2)]
               case(4,6) ; Nf(:,fID) = [e % Nout(2), e % Nout(3)]
            end select

            firstPoint(fID) = numOfPoints + 1
            numOfPoints     = numOfPoints + product(Nf(:,fID)+1)
            numOfFElems     = numOfFElems + product(Nf(:,fID)  )
            end associate
         end do

         ! don't write if boundary doesn't have elements associated, happens for periodic conditions
         if (numOfFElems .eq. 0) return

         write(fd,'(A,I0,A,I0,A,A,A)') "ZONE N=", numOfPoints,", E=", numOfFElems, &
                                                  ',ET=QUADRILATERAL, F=FEPOINT, T="boundary_', trim(boundary % Name), '"'
!
!        Write the points
!        ----------------
         do fID=1, boundary % no_of_faces

            associate (e => elements( boundary % elements(fID) ))
            side = boundary % elementSides(fID)
            N = e % Nout
            select case (side)

               case(1)
                  do k = 0, e % Nout(3)    ; do i = 0, e % Nout(1)
                     write(fd,trim(formatout)) e % xOut(:,i,0,k), e % outputVars(:,i,0,k)
                  end do                ; end do

               case(2)
                  do k = 0, e % Nout(3)    ; do i = 0, e % Nout(1)
                     write(fd,trim(formatout)) e % xOut(:,i,e % Nout(2),k), e % outputVars(:,i,e % Nout(2),k)
                  end do                ; end do

               case(3)
                  do j = 0, e % Nout(2)    ; do i = 0, e % Nout(1)
                     write(fd,trim(formatout)) e % xOut(:,i,j,0), e % outputVars(:,i,j,0)
                  end do                ; end do

               case(4)
                  do k = 0, e % Nout(3)    ; do j = 0, e % Nout(2)
                     write(fd,trim(formatout)) e % xOut(:,e % Nout(1),j,k), e % outputVars(:,e % Nout(1),j,k)
                  end do                ; end do

               case(5)
                  do j = 0, e % Nout(2)    ; do i = 0, e % Nout(1)
                     write(fd,trim(formatout)) e % xOut(:,i,j,e % Nout(3)), e % outputVars(:,i,j,e % Nout(3))
                  end do                ; end do

               case(6)
                  do k = 0, e % Nout(3)    ; do j = 0, e % Nout(2)
                     write(fd,trim(formatout)) e % xOut(:,0,j,k), e % outputVars(:,0,j,k)
                  end do                ; end do

            end select

            end associate
         end do
!
!        Write the elems connectivity
!        ----------------------------
         do fID = 1, boundary % no_of_faces

            do j = 0, Nf(2,fID) - 1 ; do i = 0, Nf(1,fID) - 1
               corners =  [ ij2localDOFStats(i,j,Nf(:,fID)), ij2localDOFStats(i+1,j,Nf(:,fID)), &
                            ij2localDOFStats(i+1,j+1,Nf(:,fID)), ij2localDOFStats(i,j+1,Nf(:,fID)) ] + firstPoint(fID)
               write(fd,*) corners
            end do                  ; end do

         end do

      end subroutine WriteSingleBoundaryZoneToTecplotStats

      character(len=LINE_LENGTH) function getFormat(ncols)
         implicit none
         integer, intent(in) :: ncols
         getFormat = ""
         write(getFormat,'(A,I0,A,A)') "(",ncols,PRECISION_FORMAT,")"
      end function getFormat
!
!/////////////////////////////////////////////////////////////////////////////
!
!     Variable filtering
!     ------------------
!
!/////////////////////////////////////////////////////////////////////////////
!
      subroutine buildStatsFilter()
         use OutputVariables, only: hasVariablesFlag, askedVariables, getNoOfCommas
         use Storage,         only: NSTAT, statsHasFavre
         implicit none
         integer            :: n, pos, pos2, i
         character(len=16)  :: tok

         ! Default: include all available variables (existing behaviour when no filter)
         statsVarInclude = (NSTAT .gt. 0)
         favreVarInclude = statsHasFavre

         if (.not. hasVariablesFlag) return

         ! User specified a list: start with nothing, include only what is requested
         statsVarInclude = .false.
         favreVarInclude = .false.

         n = getNoOfCommas(trim(askedVariables)) + 1
         pos = 0
         do i = 1, n
            tok = ""
            if (i .lt. n) then
               pos2 = index(trim(askedVariables(pos+1:)), ",") + pos
               read(askedVariables(pos+1:pos2), *, err=99, end=99) tok
               pos = pos2
            else
               read(askedVariables(pos+1:), *, err=99, end=99) tok
            end if
            tok = adjustl(trim(tok))
            call applyStatsToken(tok)
            99 continue
         end do

      end subroutine buildStatsFilter

      subroutine applyStatsToken(tok)
         implicit none
         character(len=*), intent(in) :: tok
         integer :: i

         ! Check against Reynolds stat names (Umean, Vmean, Wmean, Sxx, ..., Syz)
         do i = 1, NSTATS_OUTVARS
            if (trim(tok) .eq. trim(STATS_OUT_NAMES(i))) then
               statsVarInclude(i) = .true.
               return
            end if
         end do

         ! Check against Favre stat names (FUU, FVV, FWW, FUV, FUW, FVW)
         do i = 1, NFAVRE_OUTVARS
            if (trim(tok) .eq. trim(FAVRE_OUT_NAMES(i))) then
               favreVarInclude(i) = .true.
               return
            end if
         end do

         ! Shorthands
         select case (trim(tok))
         case ("Vmean", "V")         ! mean velocity vector
            statsVarInclude(1:3) = .true.
         case ("Rij")                ! all Reynolds stresses
            statsVarInclude(4:9) = .true.
         case ("Fij")                ! all Favre stresses
            favreVarInclude = .true.
         case ("all")                ! everything
            statsVarInclude = .true.
            favreVarInclude = .true.
         end select
         ! Variables not meaningful for stats files (rho, p, Mach, etc.) are silently ignored

      end subroutine applyStatsToken

      subroutine printStatsOutputVariables()
         use Storage, only: NSTAT, statsHasFavre
         implicit none
         integer :: i

         write(STD_OUT,'(/)')
         call Section_Header("Output variables")
         write(STD_OUT,'(/)')
         call SubSection_Header("Selected output variables")

         if (NSTAT .gt. 0) then
            do i = 1, NSTATS_OUTVARS
               if (statsVarInclude(i)) write(STD_OUT,'(30X,A,A)') "* ", trim(STATS_OUT_NAMES(i))
            end do
         end if

         if (statsHasFavre) then
            do i = 1, NFAVRE_OUTVARS
               if (favreVarInclude(i)) write(STD_OUT,'(30X,A,A)') "* ", trim(FAVRE_OUT_NAMES(i))
            end do
         end if

      end subroutine printStatsOutputVariables

      character(len=512) function buildVarsHeader()
         use Storage, only: NSTAT, statsHasFavre
         implicit none
         integer :: i

         buildVarsHeader = 'VARIABLES = "x","y","z"'
         if (NSTAT .gt. 0) then
            do i = 1, NSTATS_OUTVARS
               if (statsVarInclude(i)) then
                  buildVarsHeader = trim(buildVarsHeader) // ',"' // trim(STATS_OUT_NAMES(i)) // '"'
               end if
            end do
         end if
         if (statsHasFavre) then
            do i = 1, NFAVRE_OUTVARS
               if (favreVarInclude(i)) then
                  buildVarsHeader = trim(buildVarsHeader) // ',"' // trim(FAVRE_OUT_NAMES(i)) // '"'
               end if
            end do
         end if

      end function buildVarsHeader

end module Stats2PltModule
