MODULE MOD_VarNameMappingsRP_Vars
!===================================================================================================================================
! VarMappings 
!===================================================================================================================================
! MODULES
IMPLICIT NONE
PUBLIC
! DERIVED QUANTITIES----------------------------------------------------------------------------------------------------------------
TYPE tDerivedQ            
  INTEGER                          :: nVar 
  INTEGER                          :: nVarVisu 
  CHARACTER(LEN=255),ALLOCATABLE   :: VarName(:)
  INTEGER,ALLOCATABLE              :: Ind(:)
  INTEGER,ALLOCATABLE              :: IndGlobal(:)
END TYPE tDerivedQ

INTEGER                            :: max_nVarVisu
TYPE(tDerivedQ)                    :: Cons
TYPE(tDerivedQ)                    :: Prim
!-----------------------------------------------------------------------------------------------------------------------------------
!===================================================================================================================================
END MODULE MOD_VarNameMappingsRP_Vars
