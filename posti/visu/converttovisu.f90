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
!> Contains routines that convert the calculated FV or DG quantities to the visualization grid.
!> The routines are split into surface and volume data. Also there is a routine that handels generic data like additional arrays
!> or data from non-state files.
!===================================================================================================================================
MODULE MOD_Posti_ConvertToVisu
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PRIVATE
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! Private Part ---------------------------------------------------------------------------------------------------------------------
! Public Part ----------------------------------------------------------------------------------------------------------------------

!INTERFACE ConvertToSurfVisu_DG
  !MODULE PROCEDURE ConvertToSurfVisu_DG
!END INTERFACE
!PUBLIC:: ConvertToSurfVisu_DG

INTERFACE ConvertToVisu_GenericData
  MODULE PROCEDURE ConvertToVisu_GenericData
END INTERFACE
PUBLIC:: ConvertToVisu_GenericData

#if FV_ENABLED

!INTERFACE ConvertToSurfVisu_FV
  !MODULE PROCEDURE ConvertToSurfVisu_FV
!END INTERFACE
!PUBLIC:: ConvertToSurfVisu_FV

#if FV_RECONSTRUCT
INTERFACE ConvertToVisu_FV_Reconstruct
  MODULE PROCEDURE ConvertToVisu_FV_Reconstruct
END INTERFACE
PUBLIC:: ConvertToVisu_FV_Reconstruct
#endif /* FV_RECONSTRUCT */
#endif /* FV_ENABLED */

CONTAINS

!!===================================================================================================================================
!!> Perform a ChangeBasis of the calculated surface DG quantities to the visualization grid.
!!===================================================================================================================================
!SUBROUTINE ConvertToSurfVisu_DG()
!USE MOD_Globals
!USE MOD_PreProc
!USE MOD_Visu_Vars
!USE MOD_Interpolation      ,ONLY: GetVandermonde
!USE MOD_ChangeBasisByDim   ,ONLY: ChangeBasisSurf
!USE MOD_Interpolation_Vars ,ONLY: NodeType,NodeTypeVisu
!IMPLICIT NONE
!! INPUT / OUTPUT VARIABLES
!!-----------------------------------------------------------------------------------------------------------------------------------
!! LOCAL VARIABLES
!INTEGER            :: iSide,iVar,iVarVisu,iVarCalc
!REAL,ALLOCATABLE   :: Vdm_N_NVisu(:,:)                  ! Vandermonde from state to visualisation nodes
!!===================================================================================================================================
!SWRITE(*,*) "[DG] convert to surface visu grid"

!! compute UVisu_DG
!ALLOCATE(Vdm_N_NVisu(0:NVisu,0:PP_N))
!CALL GetVandermonde(PP_N,NodeType,NVisu,NodeTypeVisuPosti,Vdm_N_NVisu,modal=.FALSE.)
!! convert DG solution to UVisu_DG
!SDEALLOCATE(USurfVisu_DG)
!ALLOCATE(USurfVisu_DG(0:NVisu,0:ZDIM(NVisu),0:0,nBCSidesVisu_DG,nVarSurfVisuAll))
!DO iVar=1,nVarDep
  !IF (mapAllVarsToSurfVisuVars(iVar).GT.0) THEN
    !iVarCalc = mapDepToCalc(iVar)
    !iVarVisu = mapAllVarsToSurfVisuVars(iVar)
    !DO iSide = 1,nBCSidesVisu_DG
      !CALL ChangeBasisSurf(PP_N,NVisu,Vdm_N_NVisu,USurfCalc_DG(:,:,iSide,iVarCalc),USurfVisu_DG(:,:,0,iSide,iVarVisu))
    !END DO
  !END IF
!END DO
!SDEALLOCATE(Vdm_N_NVisu)
!END SUBROUTINE ConvertToSurfVisu_DG


#if FV_ENABLED
!!===================================================================================================================================
!!> Convert the calculated surface FV quantities to the visualization grid.
!!===================================================================================================================================
!SUBROUTINE ConvertToSurfVisu_FV()
!USE MOD_Globals
!USE MOD_PreProc
!USE MOD_Visu_Vars         ,ONLY: nVarDep,VarnamesAll,mapDepToCalc_FV
!USE MOD_Visu_Vars         ,ONLY: mapAllVarsToSurfVisuVars,USurfVisu_FV,USurfCalc_FV
!#if !(FV_RECONSTRUCT)
!USE MOD_Visu_Vars         ,ONLY: nBCSidesVisu_FV
!#endif
!IMPLICIT NONE
!! INPUT / OUTPUT VARIABLES
!!-----------------------------------------------------------------------------------------------------------------------------------
!! LOCAL VARIABLES
!INTEGER            :: iVar,iVarVisu
!#if !(FV_RECONSTRUCT)
!INTEGER            :: iSide,p,q
!#endif
!!===================================================================================================================================
!SWRITE(*,*) "[FV/FVRE] convert to surface visu grid"
!! compute UVisu_FV
!DO iVar=1,nVarDep
  !iVarVisu = mapAllVarsToSurfVisuVars(iVar)
  !IF (iVarVisu.GT.0) THEN
    !SWRITE(*,*) "    ", TRIM(VarnamesAll(iVar))
!#if FV_RECONSTRUCT
    !USurfVisu_FV(:,:,0,:,iVarVisu) = USurfCalc_FV(:,:,:,mapDepToCalc_FV(iVar))
!#else
    !! No reconstruction: Calculations are done directly on the subcells (PP_N+1), visualization is done on 2*(PP_N+1)-1 points
    !DO iSide = 1,nBCSidesVisu_FV
      !DO q=0,PP_NZ; DO p=0,PP_N
        !USurfVisu_FV(p*2:p*2+1, q*2:q*2+1*(PP_dim-2),0,iSide,iVarVisu) = USurfCalc_FV(p,q,iSide,mapDepToCalc_FV(iVar))
      !END DO; END DO
    !END DO
!#endif
  !END IF
!END DO

!END SUBROUTINE ConvertToSurfVisu_FV



#if FV_RECONSTRUCT
!===================================================================================================================================
!>
!===================================================================================================================================
SUBROUTINE ConvertToVisu_FV_Reconstruct(&
#if PARABOLIC
    gradUx_calc,gradUy_calc,gradUz_calc, &
#endif
    UPrim_Quad,gradUxi_Quad,gradUeta_Quad )
USE MOD_Globals
USE MOD_PreProc
USE MOD_Visu_Vars
USE MOD_Mesh_Vars          ,ONLY: nElems
USE MOD_ReadInTools        ,ONLY: GETINT
USE MOD_Interpolation      ,ONLY: GetVandermonde
USE MOD_Interpolation_Vars ,ONLY: NodeType,NodeTypeVisu
USE MOD_FV_Vars            ,ONLY: FV_dx_XI_L,FV_dx_ETA_L
USE MOD_FV_Vars            ,ONLY: FV_dx_XI_R,FV_dx_ETA_R
#if PP_dim == 3
USE MOD_FV_Vars            ,ONLY: gradUzeta
USE MOD_FV_Vars            ,ONLY: FV_dx_ZETA_L
USE MOD_FV_Vars            ,ONLY: FV_dx_ZETA_R
#endif
USE MOD_EOS_Posti          ,ONLY: GetMaskPrim
#if PARABOLIC
USE MOD_Lifting_Vars       ,ONLY: gradUx, gradUy, gradUz
#endif
IMPLICIT NONE
! INPUT / OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
#if PARABOLIC
REAL,INTENT(OUT),OPTIONAL    :: gradUx_calc  (1:PP_nVarPrim,0:NVisu_FV,0:NVisu_FV,0:ZDIM(NVisu_FV),nElems_FV)
REAL,INTENT(OUT),OPTIONAL    :: gradUy_calc  (1:PP_nVarPrim,0:NVisu_FV,0:NVisu_FV,0:ZDIM(NVisu_FV),nElems_FV)
REAL,INTENT(OUT),OPTIONAL    :: gradUz_calc  (1:PP_nVarPrim,0:NVisu_FV,0:NVisu_FV,0:ZDIM(NVisu_FV),nElems_FV)
#endif
REAL,INTENT(IN),OPTIONAL     :: UPrim_Quad   (PP_nVarPrimDet,0:PP_N,0:PP_N,0:PP_NZ,nElems)
REAL,INTENT(IN),OPTIONAL     :: gradUxi_Quad (PP_nVarPrimDet,0:PP_N,0:PP_NZ,0:PP_N,nElems)
REAL,INTENT(IN),OPTIONAL     :: gradUeta_Quad(PP_nVarPrimDet,0:PP_N,0:PP_NZ,0:PP_N,nElems)
! LOCAL VARIABLES
INTEGER             :: iVar,i,j,k,iElem,iElem_FV
INTEGER             :: iVarCalc
INTEGER             :: nVarPrim,iVarPrim
INTEGER             :: mapUPrim(PP_nVarPrim)
INTEGER             :: mapUCalc(PP_nVarPrim)
INTEGER             :: maskPrim(nVarDep)
!===================================================================================================================================
! Build local maps of maximal size PP_nVarPrim:
! - mapUCalc(1:nVarPrim) = indices of the nVarPrim primitive quantities that should be visualized in the UCalc_FV array
! - mapUPrim(1:nVarPrim) = indices of the nVarPrim primitive quantities in the UPrim array
! Example:
!   If only velocityX and pressure should be visualized then:
!     nVarPrim = 2
!     mapUPrim(1) = 2     mapUCalc(1) = index of velocityX in UCalc_FV
!     mapUPrim(2) = 5     mapUCalc(2) = index of pressure  in UCalc_FV
nVarPrim = 0
iVarPrim = 0
maskPrim = GetMaskPrim()
DO iVar=1,nVarDep
  IF (maskPrim(iVar).GT.0) THEN
    iVarPrim = iVarPrim + 1
    IF (mapDepToCalc_FV(iVar).GT.0) THEN
      nVarPrim = nVarPrim + 1
      mapUPrim(nVarPrim) = iVarPrim
      mapUCalc(nVarPrim) = mapDepToCalc_FV(iVar)
    END IF
  END IF
END DO
! SWRITE(*,*) "  nVarPrim", nVarPrim
! SWRITE(*,*) "  mapUPrim", mapUPrim(1:nVarPrim)
! SWRITE(*,*) "  mapUCalc", mapUCalc(1:nVarPrim)


DO iElem_FV=1,nElems_FV
  iElem = mapFVElemsToAllElems(iElem_FV)
  DO k=0,PP_NZ; DO j=0,PP_N; DO i=0,PP_N
    DO iVar=1,nVarPrim
      iVarPrim = mapUPrim(iVar)
      iVarCalc = mapUCalc(iVar)
      UCalc_FV(i*2  ,j*2  ,k*2  ,iElem_FV,iVarCalc) = UPrim_Quad(iVarPrim,i,j,k,iElem) &
          - gradUxi_Quad  (iVarPrim,j,k,i,iElem) *   FV_dx_XI_L(j,k,i,iElem) &
          - gradUeta_Quad (iVarPrim,i,k,j,iElem) *  FV_dx_ETA_L(i,k,j,iElem)
      UCalc_FV(i*2+1,j*2  ,k*2  ,iElem_FV,iVarCalc) = UPrim_Quad(iVarPrim,i,j,k,iElem) &
          + gradUxi_Quad  (iVarPrim,j,k,i,iElem) *   FV_dx_XI_R(j,k,i,iElem) &
          - gradUeta_Quad (iVarPrim,i,k,j,iElem) *  FV_dx_ETA_L(i,k,j,iElem)
      UCalc_FV(i*2  ,j*2+1,k*2  ,iElem_FV,iVarCalc) = UPrim_Quad(iVarPrim,i,j,k,iElem)  &
          - gradUxi_Quad  (iVarPrim,j,k,i,iElem) *   FV_dx_XI_L(j,k,i,iElem) &
          + gradUeta_Quad (iVarPrim,i,k,j,iElem) *  FV_dx_ETA_R(i,k,j,iElem)
      UCalc_FV(i*2+1,j*2+1,k*2  ,iElem_FV,iVarCalc) = UPrim_Quad(iVarPrim,i,j,k,iElem) &
          + gradUxi_Quad  (iVarPrim,j,k,i,iElem) *   FV_dx_XI_R(j,k,i,iElem) &
          + gradUeta_Quad (iVarPrim,i,k,j,iElem) *  FV_dx_ETA_R(i,k,j,iElem)
#if PP_dim == 3
      UCalc_FV(i*2  ,j*2  ,k*2  ,iElem_FV,iVarCalc) = UCalc_FV(i*2  ,j*2  ,k*2  ,iElem_FV,iVarCalc)  &
          - gradUzeta(iVarPrim,i,j,k,iElem) * FV_dx_ZETA_L(i,j,k,iElem)
      UCalc_FV(i*2+1,j*2  ,k*2  ,iElem_FV,iVarCalc) = UCalc_FV(i*2+1,j*2  ,k*2  ,iElem_FV,iVarCalc)  &
          - gradUzeta(iVarPrim,i,j,k,iElem) * FV_dx_ZETA_L(i,j,k,iElem)
      UCalc_FV(i*2  ,j*2+1,k*2  ,iElem_FV,iVarCalc) = UCalc_FV(i*2  ,j*2+1,k*2  ,iElem_FV,iVarCalc)  &
          - gradUzeta(iVarPrim,i,j,k,iElem) * FV_dx_ZETA_L(i,j,k,iElem)
      UCalc_FV(i*2+1,j*2+1,k*2  ,iElem_FV,iVarCalc) = UCalc_FV(i*2+1,j*2+1,k*2  ,iElem_FV,iVarCalc)  &
          - gradUzeta(iVarPrim,i,j,k,iElem) * FV_dx_ZETA_L(i,j,k,iElem)

      UCalc_FV(i*2  ,j*2  ,k*2+1,iElem_FV,iVarCalc) = UPrim_Quad(iVarPrim,i,j,k,iElem)  &
          - gradUxi  (iVarPrim,j,k,i,iElem) *   FV_dx_XI_L(j,k,i,iElem) &
          - gradUeta (iVarPrim,i,k,j,iElem) *  FV_dx_ETA_L(i,k,j,iElem) &
          + gradUzeta(iVarPrim,i,j,k,iElem) * FV_dx_ZETA_R(i,j,k,iElem)
      UCalc_FV(i*2+1,j*2  ,k*2+1,iElem_FV,iVarCalc) = UPrim_Quad(iVarPrim,i,j,k,iElem) &
          + gradUxi  (iVarPrim,j,k,i,iElem) *   FV_dx_XI_R(j,k,i,iElem) &
          - gradUeta (iVarPrim,i,k,j,iElem) *  FV_dx_ETA_L(i,k,j,iElem) &
          + gradUzeta(iVarPrim,i,j,k,iElem) * FV_dx_ZETA_R(i,j,k,iElem)
      UCalc_FV(i*2  ,j*2+1,k*2+1,iElem_FV,iVarCalc) = UPrim_Quad(iVarPrim,i,j,k,iElem) &
          - gradUxi  (iVarPrim,j,k,i,iElem) *   FV_dx_XI_L(j,k,i,iElem) &
          + gradUeta (iVarPrim,i,k,j,iElem) *  FV_dx_ETA_R(i,k,j,iElem) &
          + gradUzeta(iVarPrim,i,j,k,iElem) * FV_dx_ZETA_R(i,j,k,iElem)
      UCalc_FV(i*2+1,j*2+1,k*2+1,iElem_FV,iVarCalc) = UPrim_Quad(iVarPrim,i,j,k,iElem)  &
          + gradUxi  (iVarPrim,j,k,i,iElem) *   FV_dx_XI_R(j,k,i,iElem) &
          + gradUeta (iVarPrim,i,k,j,iElem) *  FV_dx_ETA_R(i,k,j,iElem) &
          + gradUzeta(iVarPrim,i,j,k,iElem) * FV_dx_ZETA_R(i,j,k,iElem)
#endif
    END DO
  END DO; END DO; END DO
END DO ! iElem_FV


#if PARABOLIC
IF (PRESENT(gradUx_calc).AND.PRESENT(gradUy_calc).AND.PRESENT(gradUz_calc)) THEN
  DO iElem_FV=1,nElems_FV
    iElem = mapFVElemsToAllElems(iElem_FV)
    DO k=0,PP_NZ; DO j=0,PP_N; DO i=0,PP_N
      DO iVar=1,PP_nVarPrim
         gradUx_calc(iVar,i*2:i*2+1, j*2:j*2+1, k*2:k*2+1*(PP_dim-2), iElem_FV) = gradUx(iVar,i,j,k,iElem)
         gradUy_calc(iVar,i*2:i*2+1, j*2:j*2+1, k*2:k*2+1*(PP_dim-2), iElem_FV) = gradUy(iVar,i,j,k,iElem)
         gradUz_calc(iVar,i*2:i*2+1, j*2:j*2+1, k*2:k*2+1*(PP_dim-2), iElem_FV) = gradUz(iVar,i,j,k,iElem)
      END DO
    END DO; END DO; END DO! i,j,k=0,PP_N
  END DO
END IF
#endif
END SUBROUTINE ConvertToVisu_FV_Reconstruct

#endif /* FV_RECONSTRUCT */

#endif /* FV_ENABLED */

!===================================================================================================================================
!> This routine will read all variables that are not conservative or derived quantities and convert the ones that should be
!> visualized to the visu grid.
!> These variables include the additional data from the ElemData and FieldData datasetes as well as other datasets that are
!> present in the HDF5 file. The variables will be named DATASETNAME:VARIABLENAME if a attribute VarNames_DATASETNAME exist
!> where we can read the variable names from. If this  attribute does not exist, the name will be a generic DATASETNAME:1,2... .
!> For each dataset a new Vandermonde matrix is build to convert from the specific polynomial degree to the visu grid,
!> so the datasets are not limited to one polynomial degree. Either elementwise (2 dimensions) or pointwise (5 dimensions) datasets
!> are allowed.
!> The addtional variables will always be sorted AFTER the conservative or derived quantities.
!> If surface visualization is needed, the quantities will simply be prolonged to the surfaces.
!===================================================================================================================================
SUBROUTINE ConvertToVisu_GenericData(statefile)
USE MOD_Globals
USE MOD_PreProc
USE MOD_Visu_Vars
USE MOD_IO_HDF5            ,ONLY: HSize
USE MOD_HDF5_Input         ,ONLY: File_ID,GetVarNames
USE MOD_HDF5_Input         ,ONLY: OpenDataFile,ReadArray,CloseDataFile,DatasetExists,ReadAttribute,GetDataSize
USE MOD_Mesh_Vars          ,ONLY: nElems,offsetElem,nBCSides,ElemToSide
USE MOD_StringTools        ,ONLY: STRICMP,split_string
USE MOD_Interpolation      ,ONLY: GetVandermonde
USE MOD_ChangeBasisByDim   ,ONLY: ChangeBasisVolume,ChangeBasisSurf
USE MOD_Interpolation_Vars ,ONLY: NodeType,NodeTypeVisu
USE MOD_Interpolation_Vars ,ONLY: L_Minus,L_Plus
USE MOD_ProlongToFace      ,ONLY: EvalElemFace
USE MOD_Mappings           ,ONLY: buildMappings
USE MOD_Visu_Avg2D         ,ONLY: Average2D,BuildVandermonds_Avg2D
USE MOD_SG_Vars            ,ONLY: nSGElems
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
CHARACTER(LEN=255),INTENT(IN)  :: statefile   !< HDF5 state file
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                        :: iVarVisu,iElem_DG,iElem_FV,iElem,iVarDataset,iVar,iVar2,ii,jj
INTEGER                        :: substring_count,nDims,nVal,nSize,nSizeZ,p,q
INTEGER                        :: iSide,iSide_DG,iSide_FV,locSide
CHARACTER(LEN=255)             :: substrings(2),DatasetName,VariableName,DataSetOld
LOGICAL                        :: datasetFound,varnamesExist,datasetChanged
REAL,ALLOCATABLE               :: ElemData(:,:,:),FieldData(:,:,:,:,:,:)
REAL,ALLOCATABLE               :: Vdm_DG_Visu(:,:),Vdm_FV_Visu(:,:)
REAL,ALLOCATABLE               :: Uface_tmp(:,:,:),Uface(:,:)
CHARACTER(LEN=255),ALLOCATABLE :: DatasetVarNames(:)
INTEGER,ALLOCATABLE            :: S2V2(:,:,:,:,:)
INTEGER                        :: mapIdentityDG(nElems)
INTEGER                        :: mapIdentityFV(nElems)
INTEGER                        :: mapVarIdentity(1)
REAL,ALLOCATABLE               :: FieldData_DG(:,:,:,:,:),FieldData_FV(:,:,:,:,:)
REAL,ALLOCATABLE               :: FVdouble(:,:)
INTEGER                        :: iVarAdd,iSGElem
!===================================================================================================================================
SWRITE(*,*) "Convert generic datasets to Visu grid"
! Open HDF5 file
CALL OpenDataFile(statefile,create=.FALSE.,single=.FALSE.,readOnly=.TRUE.,communicator=MPI_COMM_WORLD)

DataSetOld = ''  ! Used to decide if arrays and Vandermonde matrix should be re-allocated

mapVarIdentity = 1
mapIdentityDG = 0
mapIdentityFV = 0
DO iElem=1,nElems
  IF (FV_Elems_loc(iElem,1)) THEN !TODO!!!
    mapIdentityFV(iElem)=iElem
  ELSE
    mapIdentityDG(iElem)=iElem
  END IF
END DO

! Loop over all generic variables that should be visualized - sorted after the dependant variables
DO iVar=nVarDep+1,nVarAll
  ! Check if this variable should be visualized
  IF ((mapAllVarsToVisuVars(iVar)).GT.0) THEN
    ! The format of the generic data varnames is DATASETNAME:VARIABLENAME - split into DATASETNAME and VARIABLENAME
    CALL split_string(TRIM(VarnamesAll(iVar)),':',substrings,substring_count)
    ! If we find more than one substring, the variable is additional data
    IF (substring_count.GT.1) THEN
      ! Store dataset and variable name
      DatasetName =  TRIM(substrings(1))
      VariableName = TRIM(substrings(2))

      ! Check if we have a new dataset
      datasetChanged = .NOT.STRICMP(TRIM(DataSetName),TRIM(DataSetOld))

      SWRITE(*,*) "Convert variable ",TRIM(VariableName)," from dataset ", TRIM(DatasetName)

      ! Get metadata if dataset changed
      IF (datasetChanged) THEN
        ! Try to open the dataset
        CALL DatasetExists(File_ID,TRIM(DatasetName),datasetFound)
        ! Abort if the dataset was not found
        IF (.NOT.datasetFound)  THEN
          CALL CloseDataFile()
          CALL Abort(__STAMP__,'Dataset '//TRIM(DatasetName)//' does not exist.')
        END IF
        ! Get dimensions of the dataset and store number of variables as well as size of array
        CALL GetDataSize(File_ID,TRIM(DatasetName),nDims,HSize)
        nVal   = INT(HSize(1))
        nSize  = INT(HSize(2))
#if PP_dim == 3
        nSizeZ = nSize
#else
        nSizeZ = 1
#endif
        SDEALLOCATE(DataSetVarNames)
        CALL GetVarNames("VarNames_"//TRIM(DatasetName),DatasetVarNames,varnamesExist)
      END IF
      WRITE (*,*) "varnamesExist", varnamesExist

      iVarDataset = 0
      ! loop over all varnames
      DO iVar2=nVarDep+1,nVarAll
        CALL split_string(TRIM(VarnamesAll(iVar2)),':',substrings,substring_count)
        IF (substring_count.GT.1) THEN
          ! if dataset is the same increase variable index inside dataset
          IF (STRICMP(substrings(1),DatasetName)) THEN
            iVarDataset = iVarDataset+1
            ! exit if variablename found
            IF (STRICMP(substrings(2),VariableName)) EXIT
          END IF
        END IF
      END DO

      ! Read in the data if we have a new dataset. Also allocate Vandermonde matrix used in conversion to visu grid.
      IF (datasetChanged.AND.(iVarDataset.GT.0)) THEN
        SELECT CASE(nDims)
        CASE(3) ! Elementwise data
          ! Allocate array and read dataset
          SDEALLOCATE(ElemData)
          ALLOCATE(ElemData(nVal,nElems,nSGElems))
          CALL ReadArray(TRIM(DatasetName),3,(/nVal,nElems,nSGElems/),offsetElem,2,RealArray=ElemData)
        CASE(6) ! Pointwise data
          ! Allocate array and read dataset
          SDEALLOCATE(FieldData)
          ALLOCATE(FieldData(nVal,nSize,nSize,nSizeZ,nElems,nSGElems))
          CALL ReadArray(TRIM(DatasetName),6,(/nVal,nSize,nSize,nSizeZ,nElems,nSGElems/),offsetElem,5,RealArray=FieldData)
          ! Get Vandermonde matrix used to convert to the visu grid
          SDEALLOCATE(Vdm_DG_Visu)
          ALLOCATE(Vdm_DG_Visu(0:NVisu,0:nSize-1))
          CALL GetVandermonde(nSize-1,NodeType,NVisu,NodeTypeVisuPosti,Vdm_DG_Visu,modal=.FALSE.)
          SDEALLOCATE(Vdm_FV_Visu)
          ALLOCATE(Vdm_FV_Visu(0:NVisu_FV,0:nSize-1))
          CALL GetVandermonde(nSize-1,NodeType,NVisu_FV,NodeTypeVisuPosti,Vdm_FV_Visu,modal=.FALSE.)
        CASE DEFAULT
          CALL Abort(__STAMP__,'Dataset '//TRIM(DatasetName)//' does not have 3 or 6 dimensions.')
        END SELECT

        ! Store current name of dataset
        DataSetOld = TRIM(DatasetName)
      END IF ! New dataset

      ! Get index of visu array that we should write to
      iVarAdd= mapAllVarsToVisuVars(iVar) - nVarMult
      iVarAdd= nMult*nVarMult+(iVarAdd-1)*nSGElems
      ! Convert the generic data to visu grid
      SELECT CASE(nDims)
      CASE(3) ! Elementwise data
        !IF (Avg2d) THEN
          !UVisu_DG(:,:,:,:,iVarVisu) = 0.
          !UVisu_FV(:,:,:,:,iVarVisu) = 0.
          !DO iElem=1,nElems
            !ii = Elem_IJK(1,iElem)
            !jj = Elem_IJK(2,iElem)
            !IF (FVAmountAvg2D(ii,jj).LE.0.5) THEN
              !iElem_DG = mapElemIJToDGElemAvg2D(ii,jj)
              !UVisu_DG(:,:,:,iElem_DG,iVarVisu) = UVisu_DG(:,:,:,iElem_DG,iVarVisu) + ElemData(iVarDataset,iElem)
            !ELSE
              !iElem_FV = mapElemIJToFVElemAvg2D(ii,jj)
              !UVisu_FV(:,:,:,iElem_FV,iVarVisu) = UVisu_FV(:,:,:,iElem_FV,iVarVisu) + ElemData(iVarDataset,iElem)
            !END IF
          !END DO
          !UVisu_DG(:,:,:,:,iVarVisu) = UVisu_DG(:,:,:,:,iVarVisu) / nElems_IJK(3)
          !UVisu_FV(:,:,:,:,iVarVisu) = UVisu_FV(:,:,:,:,iVarVisu) / nElems_IJK(3)
        !ELSE
          ! Simply write the elementwise data to all visu points
          DO iSGElem=1,nSGElems
            DO iElem_DG=1,nElems_DG
              iElem = mapDGElemsToAllElems(iElem_DG)
              UVisu_DG(:,:,:,iElem_DG,iVarAdd+iSGElem) = ElemData(iVarDataset,iElem,iSGElem)
            END DO
          END DO 
          DO iSGElem=1,nSGElems
            DO iElem_FV=1,nElems_FV
              iElem = mapFVElemsToAllElems(iElem_FV)
              UVisu_FV(:,:,:,iElem_FV,iVarAdd+iSGElem) = ElemData(iVarDataset,iElem,iSGElem)
            END DO
          END DO 
        !END IF
      CASE(6) ! Pointwise data
        !IF (Avg2d) THEN
          !IF (nSize-1.NE.PP_N) THEN
            !CALL CollectiveStop(__STAMP__,&
                !"Avg2D only works for FieldData on PP_N!")
          !END IF
          !ALLOCATE(FieldData_DG(0:nSize-1,0:nSize-1,0:nSizeZ-1,nElems_DG,1))
          !ALLOCATE(FieldData_FV(0:nSize-1,0:nSize-1,0:nSizeZ-1,nElems_FV,1))
          !DO iElem_DG=1,nElems_DG
            !FieldData_DG(:,:,:,iElem_DG,1) = FieldData(iVarDataset,:,:,:,mapDGElemsToAllElems(iElem_DG))
          !END DO
          !DO iElem_FV=1,nElems_FV
            !FieldData_FV(:,:,:,iElem_FV,1) = FieldData(iVarDataset,:,:,:,mapFVElemsToAllElems(iElem_FV))
          !END DO

          !CALL BuildVandermonds_Avg2D(nSize-1,nSize-1)
          !CALL Average2D(1,1,nSize-1,nSize-1,nElems_DG,nElems_FV,NodeType,&
              !FieldData_DG,FieldData_FV, &
              !Vdm_DGToFV,Vdm_FVToDG,Vdm_DGToVisu,FVdouble, &
              !iVar,iVar,mapVarIdentity,UVisu_DG,UVisu_FV)
          !DEALLOCATE(FieldData_DG)
          !DEALLOCATE(FieldData_FV)
        !ELSE
          ! Perform changebasis to visu grid
          DO iSGElem=1,nSGElems
            DO iElem_DG=1,nElems_DG
              iElem = mapDGElemsToAllElems(iElem_DG)
              CALL ChangeBasisVolume(nSize-1,NVisu,Vdm_DG_Visu,FieldData(iVarDataset,:,:,:,iElem,iSGElem),&
                                     UVisu_DG(:,:,:,iElem_DG,iVarAdd+iSGElem))
            END DO
          END DO
          DO iSGElem=1,nSGElems
            DO iElem_FV=1,nElems_FV
              iElem = mapFVElemsToAllElems(iElem_FV)
              CALL ChangeBasisVolume(nSize-1,NVisu_FV,Vdm_FV_Visu,FieldData(iVarDataset,:,:,:,iElem,iSGElem),&
                                     UVisu_FV(:,:,:,iElem_FV,iVarAdd+iSGElem))
            END DO
          END DO
        !END IF
      END SELECT

      !-------------------------------- Surface Visualization ---------------------!
      !IF (doSurfVisu) THEN
        !! Get index of visu array that we should write to
        !iVarVisu= mapAllVarsToSurfVisuVars(iVar)
        !SELECT CASE(nDims)
        !CASE(2) ! Elementwise data
          !! Simply write the elementwise data to all visu points on the surface
          !DO iElem_DG = 1,nElems_DG                         ! iterate over all DG visu elements
            !iElem = mapDGElemsToAllElems(iElem_DG)          ! get global element index
!#if PP_dim == 3
            !DO locSide=1,6
!#else
            !DO locSide=2,5
!#endif
              !iSide = ElemToSide(E2S_SIDE_ID,locSide,iElem) ! get global side index
              !IF (iSide.LE.nBCSides) THEN                   ! check if BC side
                !iSide_DG = mapAllBCSidesToDGVisuBCSides(iSide)  ! get DG visu side index
                !IF (iSide_DG.GT.0) THEN
                  !USurfVisu_DG(:,:,0,iSide_DG,iVarVisu) = ElemData(iVarDataset,iElem)
                !END IF
              !END IF
            !END DO
          !END DO
          !DO iElem_FV = 1,nElems_FV                         ! iterate over all FV visu elements
            !iElem = mapFVElemsToAllElems(iElem_FV)          ! get global element index
!#if PP_dim == 3
            !DO locSide=1,6
!#else
            !DO locSide=2,5
!#endif
              !iSide = ElemToSide(E2S_SIDE_ID,locSide,iElem) ! get global side index
              !IF (iSide.LE.nBCSides) THEN                   ! check if BC side
                !iSide_FV = mapAllBCSidesToFVVisuBCSides(iSide)  ! get FV visu side index
                !IF (iSide_FV.GT.0) THEN
                  !USurfVisu_FV(:,:,0,iSide_FV,iVarVisu) = ElemData(iVarDataset,iElem)
                !END IF
              !END IF
            !END DO
          !END DO
        !CASE(5) ! Pointwise data
          !! If the dataset has changed, reallocate the mapping for rotation to the master side coordinate system
          !! as well as the temporary face data arrays.
          !IF (datasetChanged) THEN
            !SDEALLOCATE(S2V2)
            !CALL buildMappings(nSize-1,S2V2=S2V2) ! Array gets allocated in this routine
            !SDEALLOCATE(Uface_tmp)
            !ALLOCATE(Uface_tmp(1,0:nSize-1,0:nSizeZ-1))
            !SDEALLOCATE(Uface)
            !ALLOCATE(Uface(0:nSize-1,0:nSizeZ-1))
          !END IF
          !! Prolong the pointwise data to the visu face and perform change basis to visu grid
          !DO iElem_DG = 1,nElems_DG                         ! iterate over all DG visu elements
            !iElem = mapDGElemsToAllElems(iElem_DG)          ! get global element index
!#if PP_dim == 3
            !DO locSide=1,6
!#else
            !DO locSide=2,5
!#endif
              !iSide = ElemToSide(E2S_SIDE_ID,locSide,iElem) ! get global side index
              !IF (iSide.LE.nBCSides) THEN                   ! check if BC side
                !iSide_DG = mapAllBCSidesToDGVisuBCSides(iSide)  ! get DG visu side index
                !IF (iSide_DG.GT.0) THEN
                  !IF(PP_NodeType.EQ.1)THEN                  ! prolong solution to face
                    !CALL EvalElemFace(1,nSize-1,FieldData(iVarDataset:iVarDataset,:,:,:,iElem),Uface_tmp(1:1,:,:),&
                                      !L_Minus,L_Plus,locSide)
                  !ELSE
                    !CALL EvalElemFace(1,nSize-1,FieldData(iVarDataset:iVarDataset,:,:,:,iElem),Uface_tmp(1:1,:,:),locSide)
                  !END IF
                  !! Turn into master side coordinate system
                  !DO q=0,nSizeZ-1; DO p=0,nSize-1
                    !Uface(p,q)=Uface_tmp(1,S2V2(1,p,q,0,locSide),S2V2(2,p,q,0,locSide))
                  !END DO; END DO
                  !! Change basis to visu grid
                  !CALL ChangeBasisSurf(nSize-1,NVisu,Vdm_DG_Visu,Uface(:,:),USurfVisu_DG(:,:,0,iSide_DG,iVarVisu))
                !END IF
              !END IF
            !END DO
          !END DO
          !! For FV, we are visualizing on a grid with 2*(PP_N+1) points. This means, to visualize generic datasets on this grid
          !! only makes sense if they are of polynomial degree PP_N! TODO: More general approach?
          !IF (nSize-1.NE.PP_N) THEN
            !SWRITE(*,*) "Can not convert variable ",TRIM(VariableName)," from dataset ", TRIM(DatasetName), "to FV visu grid",&
                        !"since size is not equal to PP_N!"
            !USurfVisu_FV(:,:,0,:,iVarVisu) = 0.
          !ELSE
            !DO iElem_FV = 1,nElems_FV                         ! iterate over all FV visu elements
              !iElem = mapFVElemsToAllElems(iElem_FV)          ! get global element index
!#if PP_dim == 3
            !DO locSide=1,6
!#else
            !DO locSide=2,5
!#endif
                !iSide = ElemToSide(E2S_SIDE_ID,locSide,iElem) ! get global side index
                !IF (iSide.LE.nBCSides) THEN                   ! check if BC side
                  !iSide_FV = mapAllBCSidesToFVVisuBCSides(iSide)  ! get DG visu side index
                  !IF (iSide_FV.GT.0) THEN
                    !CALL EvalElemFace(1,nSize-1,FieldData(iVarDataset:iVarDataset,:,:,:,iElem),Uface_tmp(1:1,:,:),locSide)
                    !! Turn into master side coordinate system
                    !DO q=0,PP_N; DO p=0,PP_N
                      !Uface(p,q)=Uface_tmp(1,S2V2(1,p,q,0,locSide),S2V2(2,p,q,0,locSide))
                    !END DO; END DO
                    !! Change basis to visu grid
                    !CALL ChangeBasisSurf(nSize-1,NVisu_FV,Vdm_FV_Visu,Uface(:,:),USurfVisu_FV(:,:,0,iSide_DG,iVarVisu))
                  !END IF
                !END IF
              !END DO
            !END DO
          !END IF
        !END SELECT
      !END IF ! doSurfVisu

    END IF ! substring_count.GT.1
  END IF ! mapAllVarsToVisuVars(iVar).GT.0

END DO !iVar=1,

! Close HDF5 file
CALL CloseDataFile()

! Cleanup of allocatable arrays
SDEALLOCATE(ElemData)
SDEALLOCATE(FieldData)
SDEALLOCATE(Vdm_DG_Visu)
SDEALLOCATE(Vdm_FV_Visu)
SDEALLOCATE(DataSetVarNames)

SDEALLOCATE(S2V2)
SDEALLOCATE(Uface)
SDEALLOCATE(Uface_tmp)

END SUBROUTINE ConvertToVisu_GenericData

END MODULE MOD_Posti_ConvertToVisu
