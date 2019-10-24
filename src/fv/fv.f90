!=================================================================================================================================
! Copyright (c) 2010-2016  Prof. Claus-Dieter Munz
! This file is part of FLEXI, a high-order accurate framework for numerically solving PDEs with discontinuous Galerkin methods.
! For more information see https://www.flexi-project.org and https://nrg.iag.uni-stuttgart.de/
!
! FLEXI is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
! as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
!
! FLEXI is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
! of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License v3.0 for more details.
!
! You should have received a copy of the GNU General Public License along with FLEXI. If not, see <http://www.gnu.org/licenses/>.
!=================================================================================================================================
#if FV_ENABLED
#include "flexi.h"

!==================================================================================================================================
!> Module for the Finite Volume sub-cells shock capturing.
!>
!> DG elements, that are detected to contain a shock/high gradients/oscillations/..., can be switched to a Finite Volume scheme.
!> A DG element of polynomial degree N is subdivided into (N+1)^3 sub-cells (to each Gauss Point/DOF one FV sub-cell).
!> The FV sub-cells of such an element are updated using FV method with 2nd order TVD reconstruction (slope limiters).
!==================================================================================================================================
MODULE MOD_FV
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PRIVATE

INTERFACE DefineParametersFV
  MODULE PROCEDURE DefineParametersFV
END INTERFACE

INTERFACE InitFV
  MODULE PROCEDURE InitFV
END INTERFACE

INTERFACE FV_Switch
  MODULE PROCEDURE FV_Switch
END INTERFACE

INTERFACE FV_Info
  MODULE PROCEDURE FV_Info
END INTERFACE

INTERFACE FV_FillIni
  MODULE PROCEDURE FV_FillIni
END INTERFACE

INTERFACE FV_DGtoFV
  MODULE PROCEDURE FV_DGtoFV
END INTERFACE

INTERFACE FinalizeFV
  MODULE PROCEDURE FinalizeFV
END INTERFACE

PUBLIC::DefineParametersFV
PUBLIC::InitFV
PUBLIC::FV_Switch
PUBLIC::FV_Info
PUBLIC::FV_FillIni
PUBLIC::FV_DGtoFV
PUBLIC::FinalizeFV
!==================================================================================================================================

CONTAINS

!==================================================================================================================================
!> Define parameters for FV
!==================================================================================================================================
SUBROUTINE DefineParametersFV()
! MODULES
USE MOD_Globals
USE MOD_ReadInTools ,ONLY: prms,addStrListEntry
#if FV_RECONSTRUCT
USE MOD_FV_Limiter  ,ONLY: DefineParametersFV_Limiter
#endif
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!==================================================================================================================================
CALL prms%SetSection('FV')
CALL prms%CreateRealOption(   'FV_IndUpperThreshold',"Upper threshold: Element is switched from DG to FV if indicator \n"//&
                                                     "rises above this value" )
CALL prms%CreateRealOption(   'FV_IndLowerThreshold',"Lower threshold: Element is switched from FV to DG if indicator \n"//&
                                                     "falls below this value")
CALL prms%CreateLogicalOption('FV_toDG_indicator'   ,"Apply additional Persson indicator to check if DG solution after switch \n"//&
                                                     "from FV to DG is valid.", '.FALSE.')
CALL prms%CreateRealOption   ('FV_toDG_limit'       ,"Threshold for FV_toDG_indicator")
CALL prms%CreateLogicalOption('FV_toDGinRK'         ,"Allow switching of FV elements to DG during Runge Kutta stages. \n"//&
                                                     "This may violated the DG timestep restriction of the element.", '.FALSE.')
CALL prms%CreateLogicalOption('FV_IniSupersample'   ,"Supersample initial solution inside each sub-cell and take mean value \n"// &
                                                     "as average sub-cell value.", '.TRUE.')
CALL prms%CreateLogicalOption('FV_IniSharp'         ,"Maintain a sharp interface in the initial solution in the FV region", '.FALSE.')
#if FV_RECONSTRUCT
CALL DefineParametersFV_Limiter()
#endif
END SUBROUTINE DefineParametersFV

!==================================================================================================================================
!> Read in parameters needed for FV sub-cells (indicator min/max and type of limiter) and allocate several arrays.
!> Build metrics for FV sub-cells and performe initial switch from DG to FV sub-cells for all troubled cells.
!==================================================================================================================================
SUBROUTINE InitFV()
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_FV_Vars
USE MOD_FV_Basis
USE MOD_Basis       ,ONLY: InitializeVandermonde
USE MOD_Indicator   ,ONLY: doCalcIndicator
USE MOD_Mesh_Vars   ,ONLY: nElems,nSides
#if FV_RECONSTRUCT
USE MOD_FV_Limiter
#endif
USE MOD_ReadInTools
USE MOD_IO_HDF5     ,ONLY: AddToElemData,ElementOut
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!==================================================================================================================================
IF(.NOT.FVInitBasisIsDone)THEN
   CALL CollectiveStop(__STAMP__,&
     'InitFV not ready to be called or already called.')
END IF
SWRITE(UNIT_StdOut,'(132("-"))')
SWRITE(UNIT_stdOut,'(A)') ' INIT FV...'

! The indicator value is used to decide where FV sub-cells are needed
doCalcIndicator=.TRUE.

! Read minimal and maximal threshold for the indicator
FV_IndLowerThreshold = GETREAL('FV_IndLowerThreshold','-99.')
FV_IndUpperThreshold = GETREAL('FV_IndUpperThreshold', '99.')

! Read flag indicating, if an additional Persson indicator should check if a FV sub-cells element really contains no oscillations
! anymore.
FV_toDG_indicator = GETLOGICAL('FV_toDG_indicator')
IF (FV_toDG_indicator) FV_toDG_limit = GETREAL('FV_toDG_limit')

! Read flag, which allows switching from FV to DG between the stages of a Runge-Kutta time step
! (this might lead to instabilities, since the time step for a DG element is smaller)
FV_toDGinRK = GETLOGICAL("FV_toDGinRK")

#if FV_RECONSTRUCT
CALL InitFV_Limiter()
#endif

ALLOCATE(FV_Elems(nElems)) ! holds information if element is DG (0) or FV (1)
! All cells are initially DG cells
FV_Elems = 0
CALL AddToElemData(ElementOut,'FV_Elems',IntArray=FV_Elems) ! append this array to HDF5 output files

! The elementwise information of 'FV_Elems' is also needed at the faces and therefore
! is 'prolongated' to the faces into the arrays 'FV_Elems_master/slave'.
! The additional 'FV_Elems_Sum' array sums up these two arrays in the following way:
!     FV_Elems_Sum = FV_Elems_master + 2 * FV_Elems_slave
! This leads to the following information stored in 'FV_Elems_Sum' per face:
!             FV_Elems_Sum  |  0 |  1 |  2 |  3 |
!   master side element is  | DG | FV | DG | FV |
!    slave side element is  | DG | DG | FV | FV |
!ALLOCATE(FV_Elems_master(1:nSides)) ! moved to InitFV_Metrics, since needed there for U_Mortar routine
ALLOCATE(FV_Elems_slave(1:nSides))
ALLOCATE(FV_Elems_Sum(1:nSides))
FV_Elems_master = 0
FV_Elems_slave = 0
FV_Elems_Sum = 0

! arrays for FV/DG statistics
ALLOCATE(FV_Elems_counter(nElems))
ALLOCATE(FV_Elems_Amount(nElems))
FV_Elems_counter  = 0
FV_Switch_counter = 0
FV_Elems_Amount = 0
CALL AddToElemData(ElementOut,'FV_Elems_Amount',RealArray=FV_Elems_Amount)

#if FV_RECONSTRUCT
! Allocate array for multi purposes:
! - For FV elements it stores the slope between the nodes next and second next to the interface.
!    |  x    x    x    x  |                          | = face, x = node
!                  <-->   ^ to this interface
!                    ^ the slope between those nodes
! - For DG elements it stores the solution at the nodes next to the interface.
!    |  x    x    x    x  |                          | = face, x = node
!                         ^ to this interface
!                      ^ the solution at this node
ALLOCATE(FV_multi_master(PP_nVarPrim,0:PP_N,0:PP_NZ,1:nSides))
ALLOCATE(FV_multi_slave (PP_nVarPrim,0:PP_N,0:PP_NZ,1:nSides))
FV_multi_slave = 0.0
FV_multi_master = 0.0

! Allocate array for FD-gradient over faces
!    | x  x  x  x | x  x  x  x |                     | = face, x = node
!               <--->
!                  ^ the slope over the face
ALLOCATE(FV_surf_gradU(PP_nVarPrim,0:PP_N,0:PP_NZ,1:nSides))

! The gradients of the primitive variables are stored at each volume integration point and
! are computed by limiting the slopes to the two adjacent points in the respective direction.
! These are physical gradients, but they are labeled ...xi/eta/zeta, since they are the slopes
! along the xi-/eta-/zeta-lines in physical space. These slopes are required to reconstruct
! the solution at the sub-cell boundaries.
ALLOCATE(gradUxi  (PP_nVarPrim,0:PP_N,0:PP_NZ,0:PP_N,nElems))
ALLOCATE(gradUeta (PP_nVarPrim,0:PP_N,0:PP_NZ,0:PP_N,nElems))
ALLOCATE(gradUzeta(PP_nVarPrim,0:PP_N,0:PP_NZ,0:PP_N,nElems))
gradUxi=0.
gradUeta=0.
gradUzeta=0.
#if PARABOLIC
! Same as gradUxi/eta/zeta, but instead of a TVD-limiter the mean value of the slopes to the
! adjacent points is used. These slopes are used to calculate the physical gradients in
! x-/y-/z-direction, which are required for the parabolic/viscous flux.
! The gradients in x-/y-/z-direction are stored in the gradUx/y/z arrays of the lifting.
ALLOCATE(gradUxi_central  (PP_nVarPrim,0:PP_N,0:PP_N,0:PP_NZ,nElems))
ALLOCATE(gradUeta_central (PP_nVarPrim,0:PP_N,0:PP_N,0:PP_NZ,nElems))
ALLOCATE(gradUzeta_central(PP_nVarPrim,0:PP_N,0:PP_N,0:PP_NZ,nElems))
gradUxi_central  =0.
gradUeta_central =0.
gradUzeta_central=0.
#endif /* PARABOLIC */
#endif /* FV_RECONSTRUCT */

! Options for initial solution
FV_IniSharp       = GETLOGICAL("FV_IniSharp",'.FALSE.')
IF (.NOT.FV_IniSharp) FV_IniSupersample = GETLOGICAL("FV_IniSupersample",'.TRUE.')

FVInitIsDone=.TRUE.
SWRITE(UNIT_stdOut,'(A)')' INIT FV DONE!'
SWRITE(UNIT_StdOut,'(132("-"))')

END SUBROUTINE InitFV

!==================================================================================================================================
!> Performe switching between DG element and FV sub-cells element (and vise versa) depending on the indicator value
!==================================================================================================================================
SUBROUTINE FV_Switch(U,U2,U3,AllowToDG)
! MODULES
USE MOD_PreProc
USE MOD_ChangeBasisByDim,ONLY: ChangeBasisVolume
USE MOD_Indicator_Vars ,ONLY: IndValue
USE MOD_Indicator      ,ONLY: IndPersson
USE MOD_FV_Vars
USE MOD_Analyze
USE MOD_Mesh_Vars      ,ONLY: nElems
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
REAL,INTENT(INOUT)          :: U (PP_nVar,0:PP_N,0:PP_N,0:PP_NZ,1:nElems)
REAL,INTENT(INOUT),OPTIONAL :: U2(PP_nVar,0:PP_N,0:PP_N,0:PP_NZ,1:nElems)
REAL,INTENT(INOUT),OPTIONAL :: U3(PP_nVar,0:PP_N,0:PP_N,0:PP_NZ,1:nElems)
LOGICAL,INTENT(IN)          :: AllowToDG  ! < if .TRUE. FV element is allowed to switch to DG
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL    :: U_DG(PP_nVar,0:PP_N,0:PP_N,0:PP_NZ)
REAL    :: U_FV(PP_nVar,0:PP_N,0:PP_N,0:PP_NZ)
REAL    :: ind
INTEGER :: iElem
!==================================================================================================================================
DO iElem=1,nElems
  IF (FV_Elems(iElem).EQ.0) THEN ! DG Element
    ! Switch DG to FV Element, if Indicator is higher then IndMin
    IF (IndValue(iElem).GT.FV_IndUpperThreshold) THEN
      ! switch Element to FV
      FV_Elems(iElem) = 1
      CALL ChangeBasisVolume(PP_nVar,PP_N,PP_N,FV_Vdm,U(:,:,:,:,iElem),U_FV)
      U(:,:,:,:,iElem) = U_FV
      IF (PRESENT(U2)) THEN
        CALL ChangeBasisVolume(PP_nVar,PP_N,PP_N,FV_Vdm,U2(:,:,:,:,iElem),U_FV)
        U2(:,:,:,:,iElem) = U_FV
      END IF
      IF (PRESENT(U3)) THEN
        CALL ChangeBasisVolume(PP_nVar,PP_N,PP_N,FV_Vdm,U3(:,:,:,:,iElem),U_FV)
        U3(:,:,:,:,iElem) = U_FV
      END IF
    END IF
  ELSE ! FV Element
    ! Switch FV to DG Element, if Indicator is lower then IndMax
    IF ((IndValue(iElem).LT.FV_IndLowerThreshold).AND.AllowToDG) THEN
      CALL ChangeBasisVolume(PP_nVar,PP_N,PP_N,FV_sVdm,U(:,:,:,:,iElem),U_DG)
      IF (FV_toDG_indicator) THEN
        ind = IndPersson(U_DG(:,:,:,:))
        IF (ind.GT.FV_toDG_limit) CYCLE
      END IF
      U(:,:,:,:,iElem) = U_DG
      FV_Elems(iElem) = 0  ! switch Elemnent to DG
      IF (PRESENT(U2)) THEN
        CALL ChangeBasisVolume(PP_nVar,PP_N,PP_N,FV_sVdm,U2(:,:,:,:,iElem),U_DG)
        U2(:,:,:,:,iElem) = U_DG
      END IF
      IF (PRESENT(U3)) THEN
        CALL ChangeBasisVolume(PP_nVar,PP_N,PP_N,FV_sVdm,U3(:,:,:,:,iElem),U_DG)
        U3(:,:,:,:,iElem) = U_DG
      END IF
    END IF
  END IF
END DO !iElem
! collect statistics
FV_Elems_counter  = FV_Elems_counter  + FV_Elems
FV_Switch_counter = FV_Switch_counter + 1
FV_Elems_Amount   = REAL(FV_Elems_Counter)/FV_Switch_counter
END SUBROUTINE FV_Switch

!==================================================================================================================================
!> Print information on the amount of FV subcells
!==================================================================================================================================
SUBROUTINE FV_Info(iter)
! MODULES
USE MOD_Globals
USE MOD_Mesh_Vars    ,ONLY: nGlobalElems
USE MOD_Analyze_Vars ,ONLY: totalFV_nElems
USE MOD_SG_Vars      ,ONLY: nSGElems
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
INTEGER(KIND=8),INTENT(IN) :: iter !< number of iterations
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!==================================================================================================================================
#if USE_MPI
IF(MPIRoot)THEN
  CALL MPI_REDUCE(MPI_IN_PLACE,totalFV_nElems,1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_WORLD,iError)
  ! totalFV_nElems is counted in PrintStatusLine
ELSE
  CALL MPI_REDUCE(totalFV_nElems,0           ,1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_WORLD,iError)
END IF
#endif
SWRITE(UNIT_stdOut,'(A,F8.3,A)')' FV amount %: ', totalFV_nElems / REAL(nGlobalElems) / nSGElems / iter*100
totalFV_nElems = 0
END SUBROUTINE FV_Info


!==================================================================================================================================
!> Initialize all FV elements and overwrite data of DG FillIni. Each subcell is supersampled with PP_N points in each space
!> dimension and the mean value is taken as value for this subcell.
!==================================================================================================================================
SUBROUTINE FV_FillIni()
USE MOD_Globals
USE MOD_PreProc
USE MOD_DG_Vars           ,ONLY: U
USE MOD_Mesh_Vars         ,ONLY: nElems
USE MOD_FV_Vars           ,ONLY: FV_Elems,FV_Vdm,FV_IniSharp,FV_IniSupersample
USE MOD_FV_Basis          ,ONLY: FV_Build_X_w_BdryX
USE MOD_Indicator         ,ONLY: CalcIndicator
USE MOD_Basis             ,ONLY: InitializeVandermonde
USE MOD_Interpolation     ,ONLY: GetNodesAndWeights
USE MOD_ChangeBasis       ,ONLY: ChangeBasis2D_XYZ, ChangeBasis3D_XYZ
USE MOD_ChangeBasisByDim  ,ONLY: ChangeBasisVolume
USE MOD_Mesh_Vars         ,ONLY: Elem_xGP
USE MOD_Equation_Vars     ,ONLY: IniExactFunc
USE MOD_Exactfunc         ,ONLY: ExactFunc
USE MOD_Interpolation_Vars,ONLY: NodeType,NodeTypeVISUInner
USE MOD_ReadInTools       ,ONLY: GETLOGICAL
USE MOD_SG_Vars           ,ONLY: nDimStochPositions,xiDet
USE MOD_SG_Vars           ,ONLY: nQPTotal,xiQP_nDim,SG_Vdm_QuadOrth
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                :: i,iElem, j,k,ii,jj,kk,iVar, iVarDet, iQP, iDimStoch
REAL                   :: FV_w,FV_BdryX(0:PP_N+1)
REAL,DIMENSION(0:PP_N) :: FV_X,xGP,wGP,wBary,SubxGP
REAL                   :: VDM(0:PP_N,0:PP_N,0:PP_N)
REAL,ALLOCATABLE       :: xx(:,:,:,:)
REAL                   :: U_Quad(PP_nVarDet,0:PP_N,0:PP_N,0:PP_NZ,1:nElems,1:nQPTotal)
REAL                   :: tmp_Quad(PP_nVarDet,0:PP_N,0:PP_N,0:PP_NZ,1:nQPTotal)
REAL                   :: tmp(PP_nVar,0:PP_N,0:PP_N,0:PP_NZ)
REAL                   :: Elem_xFV(1:3,0:PP_N,0:PP_N,0:PP_NZ)
REAL                   :: xi(1:SGV_ARRAY_LEN)
!===================================================================================================================================
! initial call of indicator
CALL CalcIndicator(U,0.)
FV_Elems = 0
! Switch DG elements to FV if necessary (converts initial DG solution to FV solution)
CALL FV_Switch(U,AllowToDG=.FALSE.)

IF (.NOT.FV_IniSharp) THEN
  ! Super sample initial solution of all FV elements. Necessary if already initial DG solution contains oscillations, which
  ! may lead to non valid solutions inside a sub-cell.!
  !!! THIS IS EXPENSIVE !!!
  IF (FV_IniSupersample) THEN

    CALL GetNodesAndWeights(PP_N,NodeType,xGP,wGP,wBary)
    CALL FV_Build_X_w_BdryX(PP_N,FV_X,FV_w,FV_BdryX)
    DO i=0,PP_N
      ! compute equidistant supersampling points inside FV sub-cell
      DO j=0,PP_N
        SubxGP(j) = FV_BdryX(i) + (j+0.5)/(PP_N+1)*FV_w
      END DO
      ! build Vandermonde for mapping the whole interval [-1,1] to the i-th FV subcell
      CALL InitializeVandermonde(PP_N,PP_N,wBary,xGP,SubxGP,VDM(:,:,i))
    END DO


    ALLOCATE(xx(1:3,0:PP_N,0:PP_N,0:PP_NZ)) ! coordinates supersampled to FV subcell
    DO iElem=1,nElems
      IF (FV_Elems(iElem).EQ.0) CYCLE ! DG element
      DO k=0,PP_NZ
        DO j=0,PP_N
          DO i=0,PP_N
            ! supersample coordinates to i,j,k-th subcells
#if PP_dim == 3
            CALL ChangeBasis3D_XYZ(3,PP_N,PP_N,Vdm(:,:,i),Vdm(:,:,j),Vdm(:,:,k),Elem_xGP(1:3,:,:,:,iElem),xx)
#else
            CALL ChangeBasis2D_XYZ(3,PP_N,PP_N,Vdm(:,:,i),Vdm(:,:,j),Elem_xGP(1:3,:,:,0,iElem),xx(:,:,:,0))
#endif
            ! evaluate ExactFunc for all supersampled points of subcell (i,j,k)
            DO kk=0,PP_NZ; DO jj=0,PP_N; DO ii=0,PP_N
              DO iQP=1,nQPTotal
                xi(:) = xiDet(:)
                DO iDimStoch=1,PP_nDimStoch
                  xi(nDimStochPositions(iDimStoch)) = xiQP_nDim(iDimStoch,iQP)
                END DO
                CALL ExactFunc(IniExactFunc,0.,xx(1:3,ii,jj,kk),tmp_Quad(:,ii,jj,kk,iQP),xi)
              END DO
              DO iVarDet=1,PP_nVarDet
                tmp(PP_iVarDet,ii,jj,kk)=MATMUL(SG_Vdm_QuadOrth,tmp_Quad(iVarDet,ii,jj,kk,:))
              END DO
            END DO; END DO; END DO
            ! mean value
            DO iVar=1,PP_nVar
              U(iVar,i,j,k,iElem) = SUM(tmp(iVar,:,:,:)) / ((PP_N+1)**2*(PP_NZ+1))
            END DO
          END DO ! i
        END DO ! j
      END DO !k
    END DO ! iElem=1,nElems
    DEALLOCATE(xx)
  END IF
ELSE
  ! maintain a sharp interface in the FV region
  DO iElem=1,nElems
    IF (FV_Elems(iElem).EQ.0) CYCLE ! DG element
    ! get coordinates of the FV elements
    CALL ChangeBasisVolume(3,PP_N,PP_N,FV_Vdm,Elem_xGP(1:3,:,:,:,iElem),Elem_xFV(1:3,:,:,:))
    DO k=0,PP_NZ
      DO j=0,PP_N
        DO i=0,PP_N
          DO iQP=1,nQPTotal
            xi(:) = xiDet(:)
            DO iDimStoch=1,PP_nDimStoch
              xi(nDimStochPositions(iDimStoch)) = xiQP_nDim(iDimStoch,iQP)
            END DO
            CALL ExactFunc(IniExactFunc,0.,Elem_xFV(:,i,j,k),U_Quad(:,i,j,k,iElem,iQP),xi)
          END DO ! iQP
          DO iVarDet=1,PP_nVarDet
            U(PP_iVarDet,i,j,k,iElem)=MATMUL(SG_Vdm_QuadOrth,U_Quad(iVarDet,i,j,k,iElem,:))
          END DO
        END DO ! i
      END DO ! j
    END DO !k
  END DO ! iElem=1,nElems
END IF
END SUBROUTINE FV_FillIni

!==================================================================================================================================
!> Switch DG solution at faces between a DG element and a FV sub-cells element to Finite Volume.
!==================================================================================================================================
SUBROUTINE FV_DGtoFV(nVar,U_master,U_slave)
! MODULES
USE MOD_PreProc
USE MOD_Globals
USE MOD_ChangeBasisByDim ,ONLY: ChangeBasisSurf
USE MOD_FV_Vars
USE MOD_Mesh_Vars   ,ONLY: firstInnerSide,lastMPISide_MINE,nSides
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
INTEGER,INTENT(IN) :: nVar
REAL,INTENT(INOUT) :: U_master(nVar,0:PP_N,0:PP_NZ,1:nSides) !< Solution on master side
REAL,INTENT(INOUT) :: U_slave (nVar,0:PP_N,0:PP_NZ,1:nSides) !< Solution on slave side
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER     :: firstSideID,lastSideID,SideID
!==================================================================================================================================
firstSideID = firstInnerSide
lastSideID  = lastMPISide_MINE

DO SideID=firstSideID,lastSideID
  IF (FV_Elems_Sum(SideID).EQ.2) THEN
    CALL ChangeBasisSurf(nVar,PP_N,PP_N,FV_Vdm,U_master(:,:,:,SideID))
  ELSE IF (FV_Elems_Sum(SideID).EQ.1) THEN
    CALL ChangeBasisSurf(nVar,PP_N,PP_N,FV_Vdm,U_slave (:,:,:,SideID))
  END IF
END DO

END SUBROUTINE FV_DGtoFV


!==================================================================================================================================
!> Finalizes global variables of the module.
!> Deallocate allocatable arrays, nullify pointers, set *InitIsDone = .FALSE.
!==================================================================================================================================
SUBROUTINE FinalizeFV()
! MODULES
USE MOD_FV_Vars
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!==================================================================================================================================
SDEALLOCATE(FV_Elems)
!SDEALLOCATE(FV_Elems_master) ! moved to mesh.f90
SDEALLOCATE(FV_Elems_slave)
SDEALLOCATE(FV_Elems_Counter)
SDEALLOCATE(FV_Elems_Amount)
SDEALLOCATE(FV_Elems_Sum)
#if FV_RECONSTRUCT
SDEALLOCATE(FV_surf_gradU)
SDEALLOCATE(FV_multi_master)
SDEALLOCATE(FV_multi_slave)
SDEALLOCATE(gradUxi)
SDEALLOCATE(gradUeta)
SDEALLOCATE(gradUzeta)
#if PARABOLIC
SDEALLOCATE(gradUxi_central)
SDEALLOCATE(gradUeta_central)
SDEALLOCATE(gradUzeta_central)
#endif
#endif

FVInitIsDone=.FALSE.
END SUBROUTINE FinalizeFV


END MODULE MOD_FV
#endif /* FV_ENABLED */
