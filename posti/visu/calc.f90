!=================================================================================================================================
! Copyright (c) 2016  Prof. Claus-Dieter Munz
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
#include "flexi.h"

!===================================================================================================================================
!> Containes routines that call the equation system specific routines to calculate primitive or derived quantities both for
!> DG and FV.
!===================================================================================================================================
MODULE MOD_Posti_Calc
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PRIVATE

INTERFACE TransformWrapper
  MODULE PROCEDURE TransformWrapper
END INTERFACE
PUBLIC:: TransformWrapper

INTERFACE CalcSurfQuantities_DG
  MODULE PROCEDURE CalcSurfQuantities_DG
END INTERFACE
PUBLIC:: CalcSurfQuantities_DG

#if FV_ENABLED
!INTERFACE CalcQuantities_FV
  !MODULE PROCEDURE CalcQuantities_FV
!END INTERFACE
!PUBLIC:: CalcQuantities_FV

INTERFACE CalcSurfQuantities_FV
  MODULE PROCEDURE CalcSurfQuantities_FV
END INTERFACE
PUBLIC:: CalcSurfQuantities_FV
#endif

INTERFACE FillCopy
  MODULE PROCEDURE FillCopy
END INTERFACE
PUBLIC:: FillCopy


CONTAINS

!===================================================================================================================================
!> 
!===================================================================================================================================
SUBROUTINE TransformWrapper()
USE MOD_Globals
USE MOD_PreProc
USE MOD_Visu_Vars
USE MOD_SG_Vars            ,ONLY: nSGElems,nSGElems_nD
#if PARABOLIC
USE MOD_Lifting_Vars       ,ONLY: gradUx,gradUy,gradUz
USE MOD_Interpolation      ,ONLY: GetVandermonde
USE MOD_Interpolation_Vars ,ONLY: NodeType
USE MOD_ChangeBasisByDim   ,ONLY: ChangeBasisVolume
#endif
USE MOD_SG_Vars            ,ONLY: nQPTotal
USE MOD_EOS_Posti_Vars     ,ONLY: nVarDepEOS
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                :: Psi(0:PP_nCoefM1)
INTEGER             :: iQPInterval,iSGElem,iSGElem_nD(PP_nDimStoch),iDimStoch
!===================================================================================================================================
! calc DG solution
SWRITE(*,*) "Calc quantities"
SWRITE(*,*) "Prepare..."

! some init stuff
!--------------------------------
SDEALLOCATE(UofXiOnThisMesh)
SDEALLOCATE(ThisMeshIsNeeded)
ALLOCATE(UofXiOnThisMesh(1:nSGElems))
ALLOCATE(ThisMeshIsNeeded(1:nSGElems))

nMult=COUNT(doCalcStoch)
nVarPerMesh=COUNT(mapAllVarsToVisuVars(nVarDepEOS+1:nVarAll).GT.0)
nVarMult=nVarVisu-nVarPerMesh
nVarVisu=nVarMult*nMult+nVarPerMesh*nSGElems

nQPUse=MERGE(nQPTotal,0,doCalcMean.OR.doCalcStd)
ThisMeshIsNeeded=(doCalcMean.OR.doCalcStd)

UofXiOnThisMesh=.FALSE.
IF(doCalcUofXi)THEN
  !get index of elem where U is evaluated and get xiEval normalized to local element
  DO iDimStoch=1,PP_nDimStoch
    iSGElem_nD(iDimStoch)=INT((xiEval(iDimStoch)+1.)/2.*nSGElems_nD(iDimStoch))+1
    iSGElem_nD(iDimStoch)=MIN(iSGElem_nD(iDimStoch),nSGElems_nD(iDimStoch))
    iSGElem_nD(iDimStoch)=MAX(iSGElem_nD(iDimStoch),1)
    xiEval(iDimStoch) = (xiEval(iDimStoch)+1.)*nSGElems_nD(iDimStoch) - 2.*iSGElem_nD(iDimStoch) +1
  END DO 
  !get 1d index from nD index
  iSGElem=iSGElem_nD(PP_nDimStoch)
  DO iDimStoch=1,PP_nDimStoch-1
    iSGElem = iSGElem + (iSGElem_nD(iDimStoch)-1)*PRODUCT(nSGElems_nD(iDimStoch+1:PP_nDimStoch))
  END DO 
  UofXiOnThisMesh(iSGElem)=.TRUE.
  ThisMeshIsNeeded(iSGElem)=.TRUE.
  SDEALLOCATE(Psi_UofXi)
  ALLOCATE(Psi_UofXi(0:PP_nCoefM1))
  CALL GetPsiOfXi()
END IF 

! call actual computation routine
!--------------------------------
SWRITE(*,*) "Transform DG..."
SDEALLOCATE(UVisu_DG)
ALLOCATE(UVisu_DG(0:NVisu,0:NVisu,0:ZDIM(NVisu),nElems_DG,1:nVarVisu))
CALL Transform(.FALSE.)

#if FV_ENABLED
IF(nElems_FV.GT.0) THEN 
  SWRITE(*,*) "Transform FV..."
  SDEALLOCATE(UVisu_FV)
  ALLOCATE(UVisu_FV(0:NVisu_FV,0:NVisu_FV,0:ZDIM(NVisu_FV),nElems_FV,1:nVarVisu))
  CALL Transform(.TRUE.)
END IF 
#endif
END SUBROUTINE TransformWrapper

!===================================================================================================================================
!> description
!===================================================================================================================================
SUBROUTINE Transform(isFV) 
! MODULES                                                                                                                          !
USE MOD_Globals
USE MOD_PreProc
USE MOD_Visu_Vars
USE MOD_Interpolation      ,ONLY: GetVandermonde
USE MOD_Interpolation_Vars ,ONLY: NodeType,NodeTypeVISUFVEqui
USE MOD_SG_Vars            ,ONLY: nSGElems
USE MOD_SG_Vars            ,ONLY: wQP_nDim,SG_VdM_OrthQuad
!----------------------------------------------------------------------------------------------------------------------------------!
IMPLICIT NONE
! INPUT / OUTPUT VARIABLES 
LOGICAL,INTENT(IN)     :: isFV
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                :: NLoc,nElemsLoc
INTEGER,POINTER        :: MapLoc(:)
CHARACTER(LEN=255)     :: NodeTypeLoc
REAL(C_DOUBLE),POINTER :: UVisu(:,:,:,:,:)     

REAL,ALLOCATABLE       :: Vdm(:,:)

INTEGER                :: iElemLoc,iElem,iSGElem,iQP
REAL                   :: Psi(0:PP_nCoefM1)
INTEGER                :: iMult,iVarCalc,iVarVisu,iVar
TYPE tUPtr
  REAL,POINTER         :: Ptr(:,:,:,:)
END TYPE tUPtr
TYPE(tUptr)            :: ULoc(3)
!===================================================================================================================================
IF(isFV)THEN
  NLoc=NVisu_FV
  MapLoc => mapFVElemsToAllElems
  nElemsLoc=nElems_FV
  NodeTypeLoc='VISU_FVEQUI'
  UVisu => UVisu_FV
ELSE
  NLoc=NVisu
  MapLoc => mapDGElemsToAllElems
  nElemsLoc=nElems_DG
  NodeTypeLoc=NodeTypeVisuPosti
  UVisu => UVisu_DG
END IF 
ALLOCATE(Vdm(0:NLoc,0:PP_N))

CALL GetVandermonde(PP_N,TRIM(NodeType),NLoc,TRIM(NodeTypeLoc),Vdm,modal=.FALSE.)
SDEALLOCATE(U_NVisu)
ALLOCATE(U_NVisu(PP_nVar,0:NLoc,0:NLoc,0:ZDIM(NLoc)))

IF(doCalcMean.OR.doCalcStd) THEN 
  SDEALLOCATE(UTmpMean)
  ALLOCATE(UTmpMean(0:NLoc,0:NLoc,0:ZDIM(NLoc),nVarCalc))
  ULoc(1)%Ptr => UTmpMean
END IF 
IF(doCalcStd) THEN 
  SDEALLOCATE(UTmpStd)
  ALLOCATE(UTmpStd(0:NLoc,0:NLoc,0:ZDIM(NLoc),nVarCalc))
  ULoc(2)%Ptr => UTmpStd
END IF 
IF(doCalcUofXi) THEN 
  SDEALLOCATE(UTmpUofXi)
  ALLOCATE(UTmpUofXi(0:NLoc,0:NLoc,0:ZDIM(NLoc),nVarCalc))
  ULoc(3)%Ptr => UTmpUofXi
END IF 
SDEALLOCATE(UTmp)
ALLOCATE(UTmp(0:NLoc,0:NLoc,0:ZDIM(NLoc),nVarCalc))

nValN=(/NLoc+1,NLoc+1,ZDIM(NLoc)+1/)
!--------------------------------

DO iElemLoc=1,nElemsLoc
  iElem=MapLoc(iElemLoc)
  IF(doCalcMean.OR.doCalcStd) UTmpMean=0.
  IF(doCalcStd)               UTmpStd=0.
  DO iSGElem=1,nSGElems
    IF(.NOT.ThisMeshIsNeeded(iSGElem)) CYCLE
    CALL InterpolateToVisuCoords(NLoc,Vdm,iElem,iSGElem)
    DO iQP=1,nQPUse
      Psi=SG_VdM_OrthQuad(iQP,:)
      CALL SortAndStochEvalAndCalc(NLoc,Psi)
      UTmpMean=UTmpMean+UTmp*wQP_nDim(iQP)/nSGElems
      IF(doCalcStd) UTmpStd=UTmpStd+UTmp*UTmp*wQP_nDim(iQP)/nSGElems
    END DO 
    IF(UofXiOnThisMesh(iSGElem)) THEN 
      CALL SortAndStochEvalAndCalc(NLoc,Psi_UofXi)
      UTmpUofXi=UTmp
    END IF 
  END DO 
  IF(doCalcStd) UTmpStd=SQRT(MAX(UTmpStd-UTmpMean*UTmpMean,TINY(1.)))
  DO iMult=1,3
    IF(.NOT.doCalcStoch(iMult)) CYCLE
    DO iVar=1,nVarDep
      IF (mapAllVarsToVisuVars(iVar).LE.0) CYCLE
      iVarCalc = mapDepToCalc(iVar)
      iVarVisu = mapAllVarsToVisuVars(iVar)+nVarMult*COUNT(doCalcStoch(1:iMult-1))
      UVisu(:,:,:,iElemLoc,iVarVisu) = ULoc(iMult)%Ptr(:,:,:,iVarCalc)
    END DO
  END DO
END DO 
END SUBROUTINE Transform


!===================================================================================================================================
!> 
!===================================================================================================================================
SUBROUTINE SortAndStochEvalAndCalc(NLoc,Psi)
USE MOD_Globals
USE MOD_PreProc
USE MOD_Visu_Vars
USE MOD_EOS_Posti          ,ONLY: CalcQuantities
USE MOD_StringTools  ,ONLY: STRICMP
IMPLICIT NONE
! INPUT / OUTPUT VARIABLES 
INTEGER,INTENT(IN)  :: NLoc
REAL,INTENT(IN)     :: Psi(0:PP_nCoefM1)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER            :: maskCalc(nVarDep)
INTEGER            :: iVarOut,iVarDet,i,j,k,ind
REAL               :: Tmp(0:PP_nCoefM1)
!#if PARABOLIC
!REAL,ALLOCATABLE,DIMENSION(:,:,:,:,:) :: gradUx_tmp,gradUy_tmp,gradUz_tmp
!#endif
!===================================================================================================================================
!TODO: replace these loops by precomputed mapping
maskCalc = 1
DO iVarOut=1,nVarDep ! iterate over all out variables
  IF (mapDepToCalc(iVarOut).LT.1) CYCLE ! check if variable must be calculated
  DO iVarDet=1,PP_nVarDet ! iterate over all in variables
    IF( STRICMP(VarnamesAll(iVarOut),VarNamesHDF5(iVarDet))) THEN
      DO k=0,ZDIM(NLoc); DO j=0,NLoc; DO i=0,NLoc
        UTmp(i,j,k,mapDepToCalc(iVarOut))=DOT_PRODUCT(Psi,U_NVisu(PP_iVarDet,i,j,k))
      END DO; END DO; END DO! i,j,k=0,PP_N
      maskCalc(iVarOut)=0 ! remove variable from maskCalc, since they now got copied and must not be calculated.
    END IF
  END DO
END DO


IF(TRIM(FileType).EQ.'State')THEN
  IF(withDGOperator.AND.PARABOLIC.EQ.1)THEN
#if PARABOLIC
    !IF(nElems_DG.EQ.nElems)THEN
      !CALL CalcQuantities(nVarCalc,nVal,mapDGElemsToAllElems,mapDepToCalc,UCalc_DG,maskCalc,gradUx,gradUy,gradUz)
    !ELSE
      !ALLOCATE(gradUx_tmp(PP_nVarPrim,0:PP_N,0:PP_N,0:PP_NZ,nElems_DG))
      !ALLOCATE(gradUy_tmp(PP_nVarPrim,0:PP_N,0:PP_N,0:PP_NZ,nElems_DG))
      !ALLOCATE(gradUz_tmp(PP_nVarPrim,0:PP_N,0:PP_N,0:PP_NZ,nElems_DG))
      !! nicer, but only as of gfortran 6+: ALLOCATE(gradUx_tmp,gradUy_tmp,gradUz_tmp,MOLD=gradUx)
      !gradUx_tmp=gradUx(:,:,:,:,mapDGElemsToAllElems)
      !gradUy_tmp=gradUy(:,:,:,:,mapDGElemsToAllElems)
      !gradUz_tmp=gradUz(:,:,:,:,mapDGElemsToAllElems)
      !CALL CalcQuantities(nVarCalc,nVal,mapDGElemsToAllElems,mapDepToCalc,UCalc_DG,maskCalc,gradUx_tmp,gradUy_tmp,gradUz_tmp)
      !DEALLOCATE(gradUx_tmp,gradUy_tmp,gradUz_tmp)
    !END IF
#endif
  ELSE
    CALL CalcQuantities(nVarCalc,nValN,(/1/),mapDepToCalc,UTmp,maskCalc)
  END IF
END IF
END SUBROUTINE SortAndStochEvalAndCalc



!===================================================================================================================================
!> 
!===================================================================================================================================
SUBROUTINE GetPsiOfXi()
USE MOD_Globals
USE MOD_PreProc
USE MOD_Visu_Vars
USE MOD_Basis        ,ONLY: LegendrePolynomialAndDerivative
USE MOD_SG_Quadrature,ONLY: PsiOfXi
USE MOD_SG_Vars      ,ONLY: Distribution,FullOrderToTensorVec
IMPLICIT NONE
! INPUT / OUTPUT VARIABLES 
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER             :: iCoef,iDimStoch
!===================================================================================================================================
Psi_UofXi = 1.
DO iCoef=0,PP_nCoefM1
  DO iDimStoch=1,PP_nDimStoch
    Psi_UofXi(iCoef) = Psi_UofXi(iCoef) * PsiOfXi( FullOrderToTensorVec(iDimStoch,iCoef) , xiEval(iDimStoch) , Distribution )
  END DO
END DO
END SUBROUTINE GetPsiOfXi



!===================================================================================================================================
!> 
!===================================================================================================================================
SUBROUTINE InterpolateToVisuCoords(NLoc,Vdm,iElem,iSGElem)
USE MOD_Globals
USE MOD_PreProc
USE MOD_Visu_Vars
USE MOD_ChangeBasisByDim   ,ONLY: ChangeBasisVolume
IMPLICIT NONE
! INPUT / OUTPUT VARIABLES 
INTEGER,INTENT(IN)  :: NLoc
REAL,INTENT(IN)     :: Vdm(:,:)
INTEGER,INTENT(IN)  :: iElem
INTEGER,INTENT(IN)  :: iSGElem
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER             :: i,j,k,iVar
!===================================================================================================================================
#if FV_ENABLED
IF(FV_Elems_loc(iElem,iSGElem)) THEN 
  DO k=0,PP_NZ; DO j=0,PP_N; DO i=0,PP_N
    DO iVar=1,PP_nVar
      U_NVisu(iVar,i*2:i*2+1, j*2:j*2+1, k*2:k*2+1*(PP_dim-2)) = U(iVar,i,j,k,iElem,iSGElem)
    END DO
  END DO; END DO; END DO
ELSE
#endif
  CALL ChangeBasisVolume(PP_nVar,PP_N,NLoc,Vdm,U(:,:,:,:,iElem,iSGElem),U_NVisu)
#if FV_ENABLED
END IF 
#endif
END SUBROUTINE InterpolateToVisuCoords



!===================================================================================================================================
!> Calc surface quantities for all DG elements.
!> 1. Prolong all independent quantities to the boudary sides that should be visualized, also prolong gradients if they are needed
!> 2. Call CalcQuantities for the surface.
!>
!> This means we only prolong the conservative or additional variables to the boundary. All dependent variables will be calculated
!> on the surface and not prolonged so we don't get interpolation errors.
!===================================================================================================================================
SUBROUTINE CalcSurfQuantities_DG()
USE MOD_Globals
USE MOD_PreProc
USE MOD_Visu_Vars
USE MOD_EOS_Posti          ,ONLY: CalcQuantities
USE MOD_Mesh_Vars          ,ONLY: nBCSides
USE MOD_Mesh_Vars          ,ONLY: NormVec,TangVec1,TangVec2
USE MOD_StringTools        ,ONLY: STRICMP
USE MOD_Interpolation      ,ONLY: GetVandermonde
USE MOD_Interpolation_Vars ,ONLY: NodeType
USE MOD_ChangeBasisByDim   ,ONLY: ChangeBasisSurf
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER            :: maskCalc(nVarDep),nValSide(3)
INTEGER            :: iSide,iSide2
REAL               :: Vdm_N_NCalc(0:NCalc,0:PP_N)
#if PARABOLIC
REAL,ALLOCATABLE   :: gradUxFace(:,:,:,:)
REAL,ALLOCATABLE   :: gradUyFace(:,:,:,:)
REAL,ALLOCATABLE   :: gradUzFace(:,:,:,:)
#endif
REAL,ALLOCATABLE   :: NormVec_loc(:,:,:,:)
REAL,ALLOCATABLE   :: TangVec1_loc(:,:,:,:)
REAL,ALLOCATABLE   :: TangVec2_loc(:,:,:,:)
!===================================================================================================================================
CALL GetVandermonde(PP_N,NodeType,NCalc,NodeType,Vdm_N_NCalc,modal=.FALSE.)

!------ Surface visualization ----------!
nValSide=(/NCalc+1,ZDIM(NCalc)+1,nBCSidesVisu_DG/)

! Allocate array that stores the calculated variables on the visualization boundary.
SDEALLOCATE(USurfCalc_DG)
ALLOCATE(USurfCalc_DG(0:NCalc,0:ZDIM(NCalc),nBCSidesVisu_DG,1:nVarCalc))
#if PARABOLIC
ALLOCATE(gradUxFace(1:PP_nVarPrim,0:NCalc,0:ZDIM(NCalc),nBCSidesVisu_DG))
ALLOCATE(gradUyFace(1:PP_nVarPrim,0:NCalc,0:ZDIM(NCalc),nBCSidesVisu_DG))
ALLOCATE(gradUzFace(1:PP_nVarPrim,0:NCalc,0:ZDIM(NCalc),nBCSidesVisu_DG))
#endif
ALLOCATE(NormVec_loc (1:3,0:NCalc,0:ZDIM(NCalc),nBCSidesVisu_DG))
ALLOCATE(TangVec1_loc(1:3,0:NCalc,0:ZDIM(NCalc),nBCSidesVisu_DG))
ALLOCATE(TangVec2_loc(1:3,0:NCalc,0:ZDIM(NCalc),nBCSidesVisu_DG))

maskCalc=1
CALL ProlongToFace_independent(nVarCalc,nBCSidesVisu_DG,nElems_DG,maskCalc,UCalc_DG,USurfCalc_DG &
#if PARABOLIC
    ,gradUxFace,gradUyFace,gradUzFace &
#endif
    )

IF(TRIM(FileType).EQ.'State')THEN
  DO iSide=1,nBCSides
    iSide2 = mapAllBCSidesToDGVisuBCSides(iSide)
    IF (iSide2.GT.0) THEN
      CALL ChangeBasisSurf(3,PP_N,NCalc,Vdm_N_NCalc,NormVec (:,:,0:PP_NZ,0,iSide),NormVec_loc (:,:,:,iSide2))
      CALL ChangeBasisSurf(3,PP_N,NCalc,Vdm_N_NCalc,TangVec1(:,:,0:PP_NZ,0,iSide),TangVec1_loc(:,:,:,iSide2))
      CALL ChangeBasisSurf(3,PP_N,NCalc,Vdm_N_NCalc,TangVec2(:,:,0:PP_NZ,0,iSide),TangVec2_loc(:,:,:,iSide2))
    END IF
  END DO
  IF(withDGOperator.AND.PARABOLIC.EQ.1)THEN
#if PARABOLIC
    CALL CalcQuantities(nVarCalc,nValSide,mapAllBCSidesToDGVisuBCSides,mapDepToCalc,USurfCalc_DG,maskCalc*(1-DepVolumeOnly),&
        gradUxFace,gradUyFace,gradUzFace,&
        NormVec_loc(:,:,:,:),TangVec1_loc(:,:,:,:),TangVec2_loc(:,:,:,:))
#endif
  ELSE
    CALL CalcQuantities(nVarCalc,nValSide,mapAllBCSidesToDGVisuBCSides,mapDepToCalc,USurfCalc_DG,maskCalc*(1-DepVolumeOnly),&
        NormVec=NormVec_loc(:,:,:,:),TangVec1=TangVec1_loc(:,:,:,:),TangVec2=TangVec2_loc(:,:,:,:))
  END IF
END IF

#if PARABOLIC
DEALLOCATE(gradUxFace)
DEALLOCATE(gradUyFace)
DEALLOCATE(gradUzFace)
#endif
DEALLOCATE(NormVec_loc )
DEALLOCATE(TangVec1_loc)
DEALLOCATE(TangVec2_loc)

END SUBROUTINE CalcSurfQuantities_DG

!===================================================================================================================================
!> Prolong all independent variables (defined as the variables in the state file, e.g. conservative variables for normal
!> Navier Stokes calculations as well as variables that can only be calculated in the volume and must be prolonged anyway)
!> to the faces that should be visualized. If the gradients are needed, they are also prolonged.
!> The solution will be on NCalc, while the gradients are coming from the DG operator and are on PP_N, which means we need to
!> interpolate them to NCalc.
!===================================================================================================================================
SUBROUTINE ProlongToFace_independent(nVar,nSides_calc,nElems_calc,maskCalc,UIn,UBoundary&
#if PARABOLIC
    ,gradUxFace,gradUyFace,gradUzFace &
#endif
    )
USE MOD_PreProc
USE MOD_Globals
USE MOD_Visu_Vars
#if PARABOLIC
USE MOD_Lifting_Vars       ,ONLY: gradUx,gradUy,gradUz
USE MOD_Mesh_Vars          ,ONLY: SideToElem
#endif
USE MOD_Interpolation_Vars ,ONLY: L_Minus,L_Plus,NodeType
USE MOD_Interpolation      ,ONLY: GetVandermonde,GetNodesAndWeights
USE MOD_StringTools        ,ONLY: STRICMP
USE MOD_Mesh_Vars          ,ONLY: nBCSides,S2V2,ElemToSide
USE MOD_ProlongToFace      ,ONLY: EvalElemFace
USE MOD_ChangeBasisByDim   ,ONLY: ChangeBasisSurf
USE MOD_Basis              ,ONLY: LagrangeInterpolationPolys
USE MOD_Mappings           ,ONLY: BuildMappings
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
INTEGER,INTENT(IN)            :: nVar,nSides_calc,nElems_calc
REAL,INTENT(IN)               :: UIn(0:NCalc,0:NCalc,0:ZDIM(NCalc),nElems_calc,1:nVar)
REAL,INTENT(OUT)              :: UBoundary(  0:NCalc,0:ZDIM(NCalc),nSides_calc,1:nVar)
#if PARABOLIC
REAL,INTENT(OUT)              :: gradUxFace(1:PP_nVarPrim,0:NCalc,0:ZDIM(NCalc),nSides_calc)
REAL,INTENT(OUT)              :: gradUyFace(1:PP_nVarPrim,0:NCalc,0:ZDIM(NCalc),nSides_calc)
REAL,INTENT(OUT)              :: gradUzFace(1:PP_nVarPrim,0:NCalc,0:ZDIM(NCalc),nSides_calc)
#endif
INTEGER,INTENT(INOUT)         :: maskCalc(nVarDep)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER             :: iVar,iVarIn,iVarOut,iSide,locSide,iElem,p,q,iElem_DG,iSide_DG
REAL                :: Uface(1,0:NCalc,0:ZDIM(NCalc))
REAL                :: Vdm_N_NCalc(0:NCalc,0:PP_N)
#if PARABOLIC
REAL                :: gradUxFace_tmp( 1:PP_nVarPrim,0:PP_N,0:PP_NZ)
REAL                :: gradUyFace_tmp( 1:PP_nVarPrim,0:PP_N,0:PP_NZ)
REAL                :: gradUxFace_tmp2(1:PP_nVarPrim,0:PP_N,0:PP_NZ)
REAL                :: gradUyFace_tmp2(1:PP_nVarPrim,0:PP_N,0:PP_NZ)
#if PP_dim == 3
REAL                :: gradUzFace_tmp( 1:PP_nVarPrim,0:PP_N,0:PP_NZ)
REAL                :: gradUzFace_tmp2(1:PP_nVarPrim,0:PP_N,0:PP_NZ)
#endif
#endif
REAL,ALLOCATABLE    :: xGP_NCalc(:),wGP_NCalc(:),wBary_NCalc(:),L_Plus_NCalc(:),L_Minus_NCalc(:)
INTEGER,ALLOCATABLE :: S2V2_NCalc(:,:,:,:,:)
!===================================================================================================================================
! We need the evalation of the lagrange polynomials and the mappings on NCalc
IF (PP_NodeType.EQ.1) THEN
  ALLOCATE(xGP_NCalc(0:NCalc), wGP_NCalc(0:NCalc), wBary_NCalc(0:NCalc))
  ALLOCATE(L_Minus_NCalc(0:NCalc), L_Plus_NCalc(0:NCalc))
  CALL GetNodesAndWeights(NCalc,NodeType,xGP_NCalc,wGP_NCalc,wBary_NCalc)
  CALL LagrangeInterpolationPolys(1.,NCalc,xGP_NCalc,wBary_NCalc,L_Plus_NCalc)
  CALL LagrangeInterpolationPolys(-1.,NCalc,xGP_NCalc,wBary_NCalc,L_Minus_NCalc)
END IF

CALL buildMappings(NCalc,S2V2=S2V2_NCalc)

! Loop over all dependent variables
SWRITE (*,*) "ProlongToFace_independent"
DO iVarOut=1,nVarDep ! iterate over all out variables
  IF (mapDepToCalc(iVarOut).LT.1) CYCLE ! check if variable must be calculated
  DO iVarIn=1,nVar_State ! iterate over all out variables
    ! Check if this variable is present in the state file, if so define it as independent
    ! Also define variables as independent that can only be computed in the volume
    IF ((STRICMP(VarnamesAll(iVarOut),VarNamesHDF5(iVarIn))).OR.(DepVolumeOnly(iVarOut).EQ.1)) THEN
      SWRITE(*,*) "  ", TRIM(VarnamesAll(iVarOut))
      iVar=mapDepToCalc(iVarOut)

      DO iElem_DG = 1,nElems_DG                         ! iterate over all DG visu elements
        iElem = mapDGElemsToAllElems(iElem_DG)          ! get global element index
#if PP_dim == 3
        DO locSide=1,6
#else
        DO locSide=2,5
#endif
          iSide = ElemToSide(E2S_SIDE_ID,locSide,iElem) ! get global side index
          IF (iSide.LE.nBCSides) THEN                   ! check if BC side
            iSide_DG = mapAllBCSidesToDGVisuBCSides(iSide)  ! get DG visu side index
            IF (iSide_DG.GT.0) THEN
              IF(PP_NodeType.EQ.1)THEN                  ! prolong solution to face
                CALL EvalElemFace(1,NCalc,UIn(:,:,:,iElem_DG,iVar:iVar),Uface(1:1,:,:),L_Minus_NCalc,L_Plus_NCalc,locSide)
              ELSE
                CALL EvalElemFace(1,NCalc,UIn(:,:,:,iElem_DG,iVar:iVar),Uface(1:1,:,:),locSide)
              END IF
              ! Rotate into coordinate system of master side
              DO q=0,ZDIM(NCalc); DO p=0,NCalc
                UBoundary(p,q,iSide_DG,iVar)=Uface(1,S2V2_NCalc(1,p,q,0,locSide),S2V2_NCalc(2,p,q,0,locSide))
              END DO; END DO
            END IF
          END IF
        END DO
      END DO
      maskCalc(iVarOut) = 0
      EXIT
    END IF
  END DO
END DO


! Also prolong the gradients if parabolic terms are needed and the DG operator is called once.
! The gradients are taken directly from the DG operator and thus live on PP_N. We need to interpolate them to NCalc.
#if PARABOLIC
IF(TRIM(FileType).EQ.'State')THEN
  IF(withDGOperator.AND.PARABOLIC.EQ.1)THEN
    CALL GetVandermonde(PP_N,NodeType,NCalc,NodeType,Vdm_N_NCalc,modal=.FALSE.)
    DO iSide=1,nBCSides
      IF (mapAllBCSidesToDGVisuBCSides(iSide).GT.0) THEN
        iElem = SideToElem(S2E_ELEM_ID,iSide)
        locSide = SideToElem(S2E_LOC_SIDE_ID, iSide)
        IF(PP_NodeType.EQ.1)THEN
          CALL EvalElemFace(PP_nVarPrim,PP_N,gradUx(:,:,:,:,iElem),gradUxFace_tmp,L_Minus,L_Plus,locSide)
          CALL EvalElemFace(PP_nVarPrim,PP_N,gradUy(:,:,:,:,iElem),gradUyFace_tmp,L_Minus,L_Plus,locSide)
#if PP_dim == 3
          CALL EvalElemFace(PP_nVarPrim,PP_N,gradUz(:,:,:,:,iElem),gradUzFace_tmp,L_Minus,L_Plus,locSide)
#endif
        ELSE
          CALL EvalElemFace(PP_nVarPrim,PP_N,gradUx(:,:,:,:,iElem),gradUxFace_tmp,locSide)
          CALL EvalElemFace(PP_nVarPrim,PP_N,gradUy(:,:,:,:,iElem),gradUyFace_tmp,locSide)
#if PP_dim == 3
          CALL EvalElemFace(PP_nVarPrim,PP_N,gradUz(:,:,:,:,iElem),gradUzFace_tmp,locSide)
#endif
        END IF
        iSide_DG = mapAllBCSidesToDGVisuBCSides(iSide)
        ! Rotate into coordinate system of master side
        DO q=0,PP_NZ; DO p=0,PP_N
          gradUxFace_tmp2(:,p,q)=gradUxFace_tmp(:,S2V2(1,p,q,0,locSide),S2V2(2,p,q,0,locSide))
          gradUyFace_tmp2(:,p,q)=gradUyFace_tmp(:,S2V2(1,p,q,0,locSide),S2V2(2,p,q,0,locSide))
#if PP_dim == 3
          gradUzFace_tmp2(:,p,q)=gradUzFace_tmp(:,S2V2(1,p,q,0,locSide),S2V2(2,p,q,0,locSide))
#endif
        END DO; END DO
        ! Interpolate to polynomial degree for calculations
        CALL ChangeBasisSurf(PP_nVarPrim,PP_N,NCalc,Vdm_N_NCalc,gradUxFace_tmp2(:,:,:),gradUxFace(:,:,:,iSide_DG))
        CALL ChangeBasisSurf(PP_nVarPrim,PP_N,NCalc,Vdm_N_NCalc,gradUyFace_tmp2(:,:,:),gradUyFace(:,:,:,iSide_DG))
#if PP_dim == 3
        CALL ChangeBasisSurf(PP_nVarPrim,PP_N,NCalc,Vdm_N_NCalc,gradUzFace_tmp2(:,:,:),gradUzFace(:,:,:,iSide_DG))
#else
        gradUzFace(:,:,:,iSide_DG) = 0.
#endif
      END IF
    END DO
  END IF ! (withDGOperator.AND.PARABOLIC.EQ.1)
END IF
#endif /*PARABOLIC*/

SDEALLOCATE(xGP_NCalc)
SDEALLOCATE(wGP_NCalc)
SDEALLOCATE(wBary_NCalc)
SDEALLOCATE(L_Minus_NCalc)
SDEALLOCATE(L_Plus_NCalc)
END SUBROUTINE ProlongToFace_independent

#if FV_ENABLED
!===================================================================================================================================
!> Calc surface quantities for all FV elements.
!> For this, the already calculated quantities as well as the gradients in the volume will simply be copied to to the side.
!> Also the mesh quantities like normal vectors are prepared and then the CalcQuantities routine is called.
!===================================================================================================================================
SUBROUTINE CalcSurfQuantities_FV()
USE MOD_Globals
USE MOD_PreProc
USE MOD_Visu_Vars
USE MOD_EOS_Posti          ,ONLY: CalcQuantities
USE MOD_Mesh_Vars          ,ONLY: nBCSides,SideToElem,ElemToSide
USE MOD_Mesh_Vars          ,ONLY: NormVec,TangVec1,TangVec2
USE MOD_StringTools        ,ONLY: STRICMP
#if PARABOLIC
USE MOD_Mesh_Vars          ,ONLY: S2V
USE MOD_Lifting_Vars       ,ONLY: gradUx,gradUy,gradUz
#endif
USE MOD_Mappings           ,ONLY: buildMappings
IMPLICIT NONE
! INPUT / OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER            :: maskCalc(nVarDep),nValSide(3)
INTEGER            :: iSide,iSide_FV,iElem,iElem_FV,p,q,dir,ijk(3),locSide
INTEGER,ALLOCATABLE          :: S2V_NVisu(:,:,:,:,:,:)
#if PARABOLIC
INTEGER            :: iVar
REAL,ALLOCATABLE   :: gradUxFace(:,:,:,:)
REAL,ALLOCATABLE   :: gradUyFace(:,:,:,:)
REAL,ALLOCATABLE   :: gradUzFace(:,:,:,:)
#endif
REAL,ALLOCATABLE   :: NormVec_loc(:,:,:,:)
REAL,ALLOCATABLE   :: TangVec1_loc(:,:,:,:)
REAL,ALLOCATABLE   :: TangVec2_loc(:,:,:,:)
!===================================================================================================================================


nValSide=(/NCalc_FV+1,ZDIM(NCalc_FV)+1,nBCSidesVisu_FV/)
CALL buildMappings(NCalc_FV,S2V=S2V_NVisu)
SDEALLOCATE(USurfVisu_FV)
ALLOCATE(USurfVisu_FV(0:NVisu_FV,0:ZDIM(NVisu_FV),0:0,nBCSidesVisu_FV,nVarSurfVisuAll))
! ===  Surface visualization ================================
! copy UCalc_FV to USurfCalc_FV
SDEALLOCATE(USurfCalc_FV)
ALLOCATE(USurfCalc_FV(0:NCalc_FV,0:ZDIM(NCalc_FV),nBCSidesVisu_FV,1:nVarCalc_FV))
DO iElem_FV = 1,nElems_FV                         ! iterate over all FV visu elements
  iElem = mapFVElemsToAllElems(iElem_FV)          ! get global element index
#if PP_dim == 3
        DO locSide=1,6
#else
        DO locSide=2,5
#endif
    iSide = ElemToSide(E2S_SIDE_ID,locSide,iElem) ! get global side index
    IF (iSide.LE.nBCSides) THEN                   ! check if BC side
      iSide_FV = mapAllBCSidesToFVVisuBCSides(iSide)  ! get FV visu side index
      IF (iSide_FV.GT.0) THEN
        DO q=0,ZDIM(NCalc_FV); DO p=0,NCalc_FV          ! map volume solution to surface solution
          ijk = S2V_NVisu(:,0,p,q,0,locSide)
          USurfCalc_FV(p,q,iSide_FV,:) = UCalc_FV(ijk(1),ijk(2),ijk(3),iElem_FV,:)
        END DO; END DO
      END IF
    END IF
  END DO
END DO

ALLOCATE(NormVec_loc (1:3,0:NCalc_FV,0:ZDIM(NCalc_FV),nBCSidesVisu_FV))
ALLOCATE(TangVec1_loc(1:3,0:NCalc_FV,0:ZDIM(NCalc_FV),nBCSidesVisu_FV))
ALLOCATE(TangVec2_loc(1:3,0:NCalc_FV,0:ZDIM(NCalc_FV),nBCSidesVisu_FV))
#if PARABOLIC
ALLOCATE(gradUxFace(1:PP_nVarPrim,0:NCalc_FV,0:ZDIM(NCalc_FV),nBCSidesVisu_FV))
ALLOCATE(gradUyFace(1:PP_nVarPrim,0:NCalc_FV,0:ZDIM(NCalc_FV),nBCSidesVisu_FV))
ALLOCATE(gradUzFace(1:PP_nVarPrim,0:NCalc_FV,0:ZDIM(NCalc_FV),nBCSidesVisu_FV))
#endif
DO iSide=1,nBCSides
  iSide_FV = mapAllBCSidesToFVVisuBCSides(iSide)
  iElem = SideToElem(S2E_ELEM_ID,iSide)
  locSide = SideToElem(S2E_LOC_SIDE_ID, iSide)
  IF (iSide_FV.GT.0) THEN
    DO q=0,PP_NZ; DO p=0,PP_N
#if FV_RECONSTRUCT
      DO dir=1,PP_dim
        NormVec_loc (dir,p*2:p*2+1,q*2:q*2+1*(PP_dim-2),iSide_FV) = NormVec (dir,p,q,1,iSide)
        TangVec1_loc(dir,p*2:p*2+1,q*2:q*2+1*(PP_dim-2),iSide_FV) = TangVec1(dir,p,q,1,iSide)
        TangVec2_loc(dir,p*2:p*2+1,q*2:q*2+1*(PP_dim-2),iSide_FV) = TangVec2(dir,p,q,1,iSide)
      END DO
#endif
#if PARABOLIC
      ijk = S2V(:,0,p,q,0,locSide)
      DO iVar=1,PP_nVarPrim
        gradUxFace(iVar,p*2:p*2+1,q*2:q*2+1*(PP_dim-2),iSide_FV) = gradUx(iVar,ijk(1),ijk(2),ijk(3),iElem)
        gradUyFace(iVar,p*2:p*2+1,q*2:q*2+1*(PP_dim-2),iSide_FV) = gradUy(iVar,ijk(1),ijk(2),ijk(3),iElem)
        gradUzFace(iVar,p*2:p*2+1,q*2:q*2+1*(PP_dim-2),iSide_FV) = gradUz(iVar,ijk(1),ijk(2),ijk(3),iElem)
      END DO
#endif
    END DO; END DO ! p,q=0,PP_N
#if !(FV_RECONSTRUCT)
    NormVec_loc (1:PP_dim,:,:,iSide_FV) = NormVec (1:PP_dim,:,:,1,iSide)
    TangVec1_loc(1:PP_dim,:,:,iSide_FV) = TangVec1(1:PP_dim,:,:,1,iSide)
    TangVec2_loc(1:PP_dim,:,:,iSide_FV) = TangVec2(1:PP_dim,:,:,1,iSide)
#endif
  END IF
END DO


SWRITE(*,*) "[FVRE] CalcSurfQuantities"
maskCalc = DepSurfaceOnly
IF(withDGOperator.AND.PARABOLIC.EQ.1)THEN
#if PARABOLIC
  CALL CalcQuantities(nVarCalc_FV,nValSide,mapAllBCSidesToFVVisuBCSides,mapDepToCalc_FV,USurfCalc_FV,maskCalc*(1-DepVolumeOnly),&
      gradUxFace,gradUyFace,gradUzFace,&
      NormVec_loc(:,:,:,:),TangVec1_loc(:,:,:,:),TangVec2_loc(:,:,:,:))
#endif
ELSE
  CALL CalcQuantities(nVarCalc_FV,nValSide,mapAllBCSidesToFVVisuBCSides,mapDepToCalc_FV,USurfCalc_FV,maskCalc*(1-DepVolumeOnly),&
      NormVec=NormVec_loc(:,:,:,:),TangVec1=TangVec1_loc(:,:,:,:),TangVec2=TangVec2_loc(:,:,:,:))
END IF

END SUBROUTINE CalcSurfQuantities_FV

#endif /* FV_ENABLED */


!==================================================================================================================================
!> Copies variable for given element range from source array to target array using VTK structure
!==================================================================================================================================
SUBROUTINE FillCopy(nVar,NlocIn,NLocOut,nElems,UIn,nElems_calc,indices,UOut,maskCalc)
USE MOD_Visu_Vars
USE MOD_Interpolation_Vars,ONLY:NodeType
USE MOD_StringTools,     ONLY: STRICMP
USE MOD_Interpolation,   ONLY: GetVandermonde
USE MOD_ChangeBasisByDim,ONLY: ChangeBasisVolume
! MODULES
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
INTEGER,INTENT(IN)    :: nVar
INTEGER,INTENT(IN)    :: NlocIn
INTEGER,INTENT(IN)    :: NlocOut
INTEGER,INTENT(IN)    :: nElems
INTEGER,INTENT(IN)    :: nElems_calc
INTEGER,INTENT(IN)    :: indices(nElems_calc)
REAL,INTENT(IN)       :: UIn(nVar,0:NlocIn,0:NlocIn,0:ZDIM(NlocIn),nElems)
REAL,INTENT(OUT)      :: UOut(0:NlocOut,0:NlocOut,0:ZDIM(NlocOut),nElems_calc,nVarCalc)
INTEGER,INTENT(INOUT) :: maskCalc(nVarDep)
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER               :: iVarOut,iVarIn
INTEGER               :: iElem,iElem_calc
REAL                  :: Vdm_NIn_NOut(0:NlocOut,0:NlocIn)
!==================================================================================================================================
IF (NLocOut.NE.NLocIn) CALL GetVandermonde(NlocIn,NodeType,NLocOut,NodeType,Vdm_NIn_NOut,modal=.FALSE.)
! Copy exisiting variables from solution array
DO iVarOut=1,nVarDep ! iterate over all out variables
  IF (mapDepToCalc(iVarOut).LT.1) CYCLE ! check if variable must be calculated
  DO iVarIn=1,nVar! iterate over all in variables
    IF( STRICMP(VarnamesAll(iVarOut),VarNamesHDF5(iVarIn))) THEN
      DO iElem_calc=1,nElems_calc ! copy variable for all elements
        iElem = indices(iElem_calc)
        IF (NLocOut.NE.NLocIn) THEN
          CALL ChangeBasisVolume(NlocIn,NlocOut,Vdm_NIn_NOut,UIn(iVarIn,:,:,:,iElem),UOut(:,:,:,iElem_calc,mapDepToCalc(iVarOut)))
        ELSE
          UOut(:,:,:,iElem_calc,mapDepToCalc(iVarOut)) = UIn(iVarIn,:,:,:,iElem)
        END IF
      END DO ! iElem
      maskCalc(iVarOut)=0 ! remove variable from maskCalc, since they now got copied and must not be calculated.
    END IF
  END DO
END DO
END SUBROUTINE FillCopy


END MODULE MOD_Posti_Calc
