#include "flexi.h"

!===================================================================================================================================
!> Module to handle the Recordpoints
!===================================================================================================================================
MODULE MOD_HDF5_OutputRP
! MODULES
USE MOD_IO_HDF5, ONLY: FILE_ID,OpenDataFile,CloseDataFile
IMPLICIT NONE
PRIVATE
!-----------------------------------------------------------------------------------------------------------------------------------
INTERFACE WriteRecordPointstoHDF5
  MODULE PROCEDURE WriteRecordPointstoHDF5
END INTERFACE

PUBLIC :: WriteRecordPointstoHDF5
!===================================================================================================================================

CONTAINS


!===================================================================================================================================
!> Subroutine to write the recordpoints to HDF5 format
!===================================================================================================================================
SUBROUTINE WriteRecordPointstoHDF5(ProjectName,MeshFileName)
! MODULES
USE MOD_Globals
USE MOD_HDF5_Output
USE MOD_Mesh_Vars       ,ONLY:NGeo,nGlobalElems
USE MOD_Interpolation_Vars,ONLY: NodeType
USE MOD_RPSet_Vars
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
CHARACTER(LEN=*),INTENT(IN)    :: ProjectName
CHARACTER(LEN=*),INTENT(IN)    :: MeshFileName
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                        :: iGr,iL,iRP,iPl
INTEGER                        :: i,j
CHARACTER(LEN=255)             :: FileName,FileString,tmp255
CHARACTER(LEN=5)               :: PlaneType
TYPE(tGroup),POINTER           :: Group
TYPE(tLine),POINTER            :: Line
TYPE(tPlane),POINTER           :: Plane
INTEGER,ALLOCATABLE            :: RPset(:),RPset2D(:,:)
REAL,ALLOCATABLE               :: x_dummy(:,:)
INTEGER                        :: s1(1),s2(2)
!===================================================================================================================================
SWRITE(UNIT_stdOut,'(a)',ADVANCE='NO')' WRITE RECORDPOINTS TO HDF5 FILE...'
FileName=TRIM(ProjectName)//'_RPSet'
FileString=TRIM(FileName)//'.h5'
CALL OpenDataFile(TRIM(Filestring),create=.TRUE.,single=.TRUE.,readOnly=.FALSE.)

CALL WriteAttribute(File_ID,'File_Type',1,StrScalar=(/'RecordPoints'/))
CALL WriteAttribute(File_ID,'MeshFile',1,StrScalar=(/TRIM(MeshFileName)/))
CALL WriteAttribute(File_ID,'NGeo',1,IntScalar=NGeo)
CALL WriteAttribute(File_ID,'NodeType',1,StrScalar=(/NodeType/))

! Write RP Set to File -------------------------------------------------------------------------------------------------------------
! 1. Groups
s1=nGroups
CALL WriteArray('GroupNames',1,s1,s1,(/0/),.FALSE.,StrArray=Groups(:)%Name)
DO iGr=1,nGroups
  Group=>Groups(iGr)
  IF(Group%nRP.GT.0) THEN
    ALLOCATE(RPset(1:Group%nRP))
    DO iRP=1,Group%nRP
      RPset(iRP)=Group%RP_ptr(iRP)%RP%ID
    END DO
    s1=Group%nRP
    CALL WriteArray(TRIM(Group%Name),1,s1,s1,(/0/),.FALSE.,IntArray=RPset)
    DEALLOCATE(RPset)
  END IF !Group%nRP.GT.0
END DO !iGr

! 2. Lines
IF(nLines.GT.0) THEN
  s1=nLines
  CALL WriteArray('LineNames',1,s1,s1,(/0/),.FALSE.,StrArray=Lines(:)%Name)
  DO iL=1,nLines
    Line=>Lines(iL)
    ALLOCATE(RPset(1:Line%nRP))
    DO iRP=1,Line%nRP
      RPset(iRP)=Line%RP_ptr(iRP)%RP%ID
    END DO
    s1=Line%nRP
    CALL WriteArray(TRIM(Line%Name),1,s1,s1,(/0/),.FALSE.,IntArray=RPset)
    DEALLOCATE(RPset)
    ! groupID
    CALL WriteAttribute(FILE_ID,'GroupID',1,DataSetname=TRIM(Line%Name),IntScalar=Line%GroupID)
  END DO !iL
END IF

! 3. Points
IF(nPoints.GT.0) THEN
  ALLOCATE(RPset(1:nPoints))
  DO iRP=1,nPoints
    RPset(iRP)=Points(iRP)%RP%ID
  END DO
  s1=nPoints
  CALL WriteArray('Points_IDlist'     ,1,s1,s1,(/0/),.FALSE.,IntArray=RPSet)
  CALL WriteArray('Points_GroupIDlist',1,s1,s1,(/0/),.FALSE.,IntArray=Points(:)%GroupID)
END IF

! 4. Planes
IF(nPlanes.GT.0) THEN
  s1=nPlanes
  CALL WriteArray('PlaneNames',1,s1,s1,(/0/),.FALSE.,StrArray=Planes(:)%Name)
  DO iPl=1,nPlanes
    Plane=>Planes(iPl)
    ALLOCATE(RPset2D(1:Plane%nRP(1),1:Plane%nRP(2)))
    DO j=1,Plane%nRP(2)
      DO i=1,Plane%nRP(1)
        RPset2D(i,j)=Plane%RP_ptr(i,j)%RP%ID
      END DO !i
    END DO !j 
    CALL WriteArray(TRIM(Plane%Name),2,Plane%nRP,Plane%nRP,(/0,0/),.FALSE.,IntArray=RPset2D)
    DEALLOCATE(RPset2D)
    ! groupID
    CALL WriteAttribute(FILE_ID,'GroupID',1,DataSetname=TRIM(Plane%Name),IntScalar=Plane%GroupID)

    ! write normal and tangent vectors in case of the BLPlane
    PlaneType=TRIM(Plane%Name(1:5))
    IF(PlaneType.EQ.TRIM("BLPla")) THEN
      s2=(/3,Plane%nRP(1)/)
      WRITE(tmp255,'(A,A)')TRIM(Plane%Name),'_NormVec'
      CALL WriteArray(tmp255,2,s2,s2,(/0,0/),.FALSE.,RealArray=Plane%NormVec)
      WRITE(tmp255,'(A,A)')TRIM(Plane%Name),'_TangVec'
      CALL WriteArray(tmp255,2,s2,s2,(/0,0/),.FALSE.,RealArray=Plane%TangVec)
    END IF
  END DO !iPl
END IF

! Offset Array
s2=(/2,nGlobalElems/)
CALL WriteArray('OffsetRP',2,s2,s2,(/0,0/),.FALSE.,IntArray=OffsetRP)

! Coordinates
ALLOCATE(x_dummy(1:3,1:nRP_global))
s2=(/3,nRP_global/)
DO iRP=1,nRP_global
  x_dummy(:,iRP)=RPlist(iRP)%RP%xF
END DO
CALL WriteArray('xF_RP',2,s2,s2,(/0,0/),.FALSE.,RealArray=x_dummy)
DO iRP=1,nRP_global
  x_dummy(:,iRP)=RPlist(iRP)%RP%x
END DO
CALL WriteArray('x_RP' ,2,s2,s2,(/0,0/),.FALSE.,RealArray=x_dummy)
DO iRP=1,nRP_global
  x_dummy(:,iRP)=RPlist(iRP)%RP%xi
END DO
CALL WriteArray('xi_RP',2,s2,s2,(/0,0/),.FALSE.,RealArray=x_dummy)
DEALLOCATE(x_dummy)
! Close the file.
CALL CloseDataFile()

SWRITE(UNIT_stdOut,'(a)',ADVANCE='YES')'DONE'
END SUBROUTINE WriteRecordPointstoHDF5




END MODULE MOD_HDF5_OutputRP
