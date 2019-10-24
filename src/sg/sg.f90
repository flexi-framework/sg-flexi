!=================================================================================================================================
! Copyright (c) 2010-2016  Prof. Claus-Dieter Munz
! This file is part of FLEXI, a high-order accurate framework for numerically solving PDEs with Quadontinuous Galerkin methods.
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

!==================================================================================================================================
!> Contains all routines for setup of SG specific operators
!==================================================================================================================================
MODULE MOD_SG
! MODULES
IMPLICIT NONE
PRIVATE
!----------------------------------------------------------------------------------------------------------------------------------
INTERFACE DefineParametersSG
  MODULE PROCEDURE DefineParametersSG
END INTERFACE

INTERFACE InitSG
  MODULE PROCEDURE InitSG
END INTERFACE

INTERFACE FinalizeSG
  MODULE PROCEDURE FinalizeSG
END INTERFACE

PUBLIC :: DefineParametersSG
PUBLIC :: InitSG
PUBLIC :: FinalizeSG
!==================================================================================================================================

CONTAINS

!==================================================================================================================================
!> Define parameters
!==================================================================================================================================
SUBROUTINE DefineParametersSG()
! MODULES
USE MOD_ReadInTools  ,ONLY: prms,addStrListEntry
USE MOD_SG_Operators ,ONLY: DefineParametersSGOperators
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!==================================================================================================================================
CALL prms%SetSection("SG")
CALL prms%CreateIntOption(          'nDimStoch'             ,"Number of stochastic dimensions.")
CALL prms%CreateIntFromStringOption('Distribution'          ,"Distribution function of random variable,"// &
                                                             "uniform or normal",'uniform')
CALL addStrListEntry(               'Distribution','uniform',SG_DIST_UNIFORM)
CALL addStrListEntry(               'Distribution','normal' ,SG_DIST_NORMAL)

CALL prms%CreateIntOption(          'M'                     ,"Stochastic polynomial degree.")
CALL prms%CreateIntOption(          'nQP'                   ,"Number of stochastic integration points.")
CALL prms%CreateIntOption(          'nQPExact'              ,"Number of stochastic integration points for expectation.")

!the following parameters are defined for each stochastic variable
CALL prms%CreateIntFromStringOption('StochVarName'          ,"Name/Usage of random variable,"// &
                                                             "amplitude, wave_speed,wave_length,dmr_angle",multiple=.TRUE.)
CALL addStrListEntry(               'StochVarName','amplitude'     ,SGV_AMPLITUDE)
CALL addStrListEntry(               'StochVarName','wavespeed'     ,SGV_WAVESPEED)
CALL addStrListEntry(               'StochVarName','wavelength'    ,SGV_WAVELENGTH)
CALL addStrListEntry(               'StochVarName','startpos_xy'   ,SGV_STARTPOS_XY)
CALL addStrListEntry(               'StochVarName','refstatevel_xy',SGV_REFSTATEVEL_XY)
CALL addStrListEntry(               'StochVarName','delta_99'      ,SGV_DELTA_99)
CALL addStrListEntry(               'StochVarName','dmr_angle'     ,SGV_DMR_ANGLE)
#if PARABOLIC
CALL addStrListEntry(               'StochVarName','viscosity'     ,SGV_VISCOSITY)
#endif
CALL prms%CreateRealOption(         'StochSigma'            ,"Variance of random input variable"          ,multiple=.TRUE.)
CALL prms%CreateRealOption(         'StochMu'               ,"Expectation of random input variable"       ,multiple=.TRUE.)
CALL prms%CreateRealArrayOption(    'UniformIntervalBounds' ,"",multiple=.TRUE.)
CALL prms%CreateIntOption(          'nSGElems'              ,"number of stochastic intervals"             ,multiple=.TRUE.)

CALL DefineParametersSGOperators()
END SUBROUTINE DefineParametersSG



!==================================================================================================================================
!> Parent routine to initialize SG information and SG operators;
!> Mostly calls other init routines.
!==================================================================================================================================
SUBROUTINE InitSG(postiIn)
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_SG_Vars
USE MOD_Interpolation_Vars ,ONLY: InterpolationInitIsDone
USE MOD_ReadInTools        ,ONLY: GETINT,GETINTFROMSTR
USE MOD_SG_Tensor          ,ONLY: CreateC,GetMultiIndFrom1DInd
USE MOD_SG_Quadrature      ,ONLY: InitQuadrature
USE MOD_SG_Operators       ,ONLY: InitSGOperators
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
LOGICAL,INTENT(IN),OPTIONAL :: postiIn
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
LOGICAL                     :: posti
INTEGER                     :: iDimStoch,MLoc,nDimStochLoc
!==================================================================================================================================
IF(SGInitIsDone.OR.(.NOT.InterpolationInitIsDone))THEN
  CALL CollectiveStop(__STAMP__,'InitSG not ready to be called or already called.')
  RETURN
END IF
SWRITE(UNIT_StdOut,'(132("-"))')
SWRITE(UNIT_stdOut,'(A)') ' INIT SG...'

posti=.FALSE.
IF(PRESENT(postiIn)) posti=postiIn

Distribution = GETINTFROMSTR('Distribution')
#ifdef PP_nDimStoch
nDimStochLoc = GETINT('nDimStoch')
IF(nDimStochLoc.NE.PP_nDimStoch) CALL Abort(__STAMP__,'nDimStoch in ini file is not equal to pre-compiled PP_nDimStoch')
#else
PP_nDimStoch = GETINT('nDimStoch')
#endif

! read in PP_M and get PP_nCoef and PP_nVar from that
#ifdef PP_M
MLoc = GETINT('M')
IF(MLoc.NE.PP_M) CALL Abort(__STAMP__,'M in ini file is not equal to pre-compiled PP_M')
#else
PP_M = GETINT('M')
PP_nCoef = 1
DO iDimStoch=1,PP_nDimStoch
  PP_nCoef = PP_nCoef*(PP_M+iDimStoch)
END DO
DO iDimStoch=1,PP_nDimStoch
  PP_nCoef = PP_nCoef/(PP_nDimStoch-(iDimStoch-1))
END DO
PP_nVar    =PP_nVarDet    *(PP_nCoef)
PP_nVarPrim=PP_nVarPrimDet*(PP_nCoef)
#endif

! Get polynomial degrees in every stoch dimension for every basis function of u
CALL GetMultiIndFrom1DInd(PP_M,PP_nCoefM1,.TRUE.,FullOrderToTensorVec)

! Get priple product tensors
CALL CreateC()

! Multi-Elem SG: Get mySGElem, nSGElems and MPI communicators
CALL StochBuildPartition(posti)

! Read in distribution parameters and names (i.e. which variables are stochastic) of stochastic variables
CALL InitStochVars()

! Init quadrature points and weights, as well as Vdm's to transform between orthogonal basis and quadrature points
! (normal and exact quadrature)
CALL InitQuadrature(posti)

! Init SG Operators (this is limited to setting pointers to the different routines according to
! Stoch1D / Stoch nD and NumFrac / ExactFrac
CALL InitSGOperators()


SGInitIsDone = .TRUE.
SWRITE(UNIT_stdOut,'(A)')' INIT SG DONE!'
SWRITE(UNIT_StdOut,'(132("-"))')
END SUBROUTINE InitSG



!==================================================================================================================================
!> Read in which variables are stochastic, and read in their distribution parameters, such as mean and variance.
!==================================================================================================================================
SUBROUTINE InitStochVars()
! MODULES                                                                                                                          !
USE MOD_Globals
USE MOD_PreProc
USE MOD_SG_Vars
USE MOD_ReadInTools ,ONLY: GETINTFROMSTR,GETREAL,GETREALARRAY,CountOption
!----------------------------------------------------------------------------------------------------------------------------------!
IMPLICIT NONE
! INPUT / OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                     :: iDimStoch,iStochVarName
REAL                        :: tmp(2)
!===================================================================================================================================
SWRITE(UNIT_stdOut,'(A)') ' Initialize stochastic variables ...'
! Define StochVarNamesAll for output
ALLOCATE(StochVarNamesAll(SGV_ARRAY_LEN))
ALLOCATE(DimStochIsSet(SGV_ARRAY_LEN))
StochVarNamesAll(SGV_AMPLITUDE)      = 'Amplitude'
StochVarNamesAll(SGV_WAVESPEED)      = 'Wavespeed'
StochVarNamesAll(SGV_WAVELENGTH)     = 'Wavelength'
StochVarNamesAll(SGV_STARTPOS_XY)    = 'StartPositionXY'
StochVarNamesAll(SGV_REFSTATEVEL_XY) = 'RefStateVelocityXY'
StochVarNamesAll(SGV_DELTA_99)       = 'Delta_99'
StochVarNamesAll(SGV_DMR_ANGLE)      = 'DMR_Angle'
#if PARABOLIC
StochVarNamesAll(SGV_VISCOSITY)      = 'Viscosity'
#endif
! Set Standard values for SGV_Variables
ALLOCATE(xiDet(1:SGV_ARRAY_LEN))
xiDet = 0.
! xiDet is set elsewhere in the code depending on what it replaces, e.g. in InitExactFunc and InitEquation!


! Get information which variables are uncertain and
! Check for multiple definitions of identical input variable
IF(CountOption('StochVarName').GT.PP_nDimStoch) &
  CALL PrintWarning("More StochVars than nDimStoch. Some StochVars will not be considered!")

ALLOCATE(StochVarNames(PP_nDimStoch))
ALLOCATE(nDimStochPositions(1:PP_nDimStoch))
DimStochIsSet(:)  = .FALSE.
DO iDimStoch=1,PP_nDimStoch
  iStochVarName                 = GETINTFROMSTR('StochVarName')
  StochVarNames(iDimStoch)      = StochVarNamesAll(iStochVarName)
  nDimStochPositions(iDimStoch) = iStochVarName
  IF (DimStochIsSet(iStochVarName)) THEN
    CALL Abort(__STAMP__,'ERROR: Stochastic variable is set multiple times in parameter file!')
  ELSE
    DimStochIsSet(iStochVarName) = .TRUE.
  END IF
  IF(iStochVarName.EQ.SGV_VISCOSITY) viscosityDim = iDimStoch
END DO
IF (.NOT.DimStochIsSet(SGV_VISCOSITY)) THEN
  CALL Abort(__STAMP__,'ERROR: Using this branch, please ALWAYS set Viscosity as a stochastic variable. var(mu)=0 is okay.') 
END IF 


! get parameters of variable (mu/sigma or interval bounds)
SELECT CASE(Distribution)
CASE(SG_DIST_NORMAL)
  ALLOCATE(StochSigma(PP_nDimStoch))
  ALLOCATE(StochMu(   PP_nDimStoch))
  StochSigma(:) = 0.
  StochMu(:)    = 0.
  DO iDimStoch=1,PP_nDimStoch
    StochSigma(iDimStoch) = GETREAL('StochSigma')
    StochMu(iDimStoch)    = GETREAL('StochMu')
  END DO
CASE(SG_DIST_UNIFORM)
  ALLOCATE(UniIntBounds(2,PP_nDimStoch))
  UniIntBounds(:,:) = 0.
  DO iDimStoch=1,PP_nDimStoch
    tmp = GETREALARRAY('UniformIntervalBounds',2)
    UniIntBounds(1,iDimStoch)=tmp(1)+(mySGElem_nD(iDimStoch)-1)/REAL(nSGElems_nD(iDimStoch))*(tmp(2)-tmp(1))
    UniIntBounds(2,iDimStoch)=tmp(1)+ mySGElem_nD(iDimStoch)   /REAL(nSGElems_nD(iDimStoch))*(tmp(2)-tmp(1))
  END DO
END SELECT
END SUBROUTINE InitStochVars



!==================================================================================================================================
!> Partition the mesh by numbers of processors. Elements are distributed equally to all processors.
!==================================================================================================================================
SUBROUTINE StochBuildPartition(posti)
! MODULES                                                                                                                          !
USE MOD_Globals
USE MOD_PreProc
USE MOD_SG_Vars
USE MOD_IO_HDF5,     ONLY: GatheredWrite
USE MOD_ReadInTools, ONLY: GETINT
!----------------------------------------------------------------------------------------------------------------------------------!
IMPLICIT NONE
! INPUT / OUTPUT VARIABLES
LOGICAL,INTENT(IN):: posti
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER           :: iSGElem,sumProcs
INTEGER           :: MeshRoot
INTEGER           :: tmp,iDimStoch
!===================================================================================================================================
SWRITE(UNIT_stdOut,'(A)') ' Build partition of SG multi-elems ...'
ALLOCATE(nSGElems_nD(PP_nDimStoch))
ALLOCATE(mySGElem_nD(PP_nDimStoch))
DO iDimStoch=1,PP_nDimStoch
  nSGElems_nD(iDimStoch)=GETINT('nSGElems','1')
END DO
nSGElems=PRODUCT(nSGElems_nD)
IF(nSGElems.GT.1)THEN
  IF(Distribution.NE.SG_DIST_UNIFORM) CALL Abort(__STAMP__,'Multi-element SG only implemented for uniform distribution')
ENDIF

IF(posti) THEN
  mySGElem_nD=1
  RETURN
  !the rest is done in Build_FV_DG_distribution
END IF

IF(nProcessors.LT.nSGElems) &
  CALL CollectiveStop(__STAMP__, "Number of processors smaller than number of Meshfiles")

SDEALLOCATE(nProcsPerSGElem)
ALLOCATE(nProcsPerSGElem(nSGElems))
nProcsPerSGElem = 1

IF(nSGElems.GT.1) gatheredWrite=.FALSE.

nProcsPerSGElem=nProcessors/nSGElems
DO iSGElem=1,MOD(nProcessors,nSGElems)
  nProcsPerSGElem(iSGElem)=nProcsPerSGElem(iSGElem)+1
END DO

! Find the SGElem for each processor
sumProcs = 0
DO iSGElem=1,nSGElems
  IF((sumProcs.LE.myRank).AND.(myRank.LT.sumProcs+nProcsPerSGElem(iSGElem))) THEN
    mySGElem = iSGElem
    EXIT
  END IF
  sumProcs = sumProcs + nProcsPerSGElem(iSGElem)
END DO

!get multi D index from mySGElem
tmp=mySGElem-1
DO iDimStoch=1,PP_nDimStoch-1
  !integer division!
  mySGElem_nD(iDimStoch)=tmp/PRODUCT(nSGElems_nD(iDimStoch+1:PP_nDimStoch))+1
  !rest
  tmp=tmp-(mySGElem_nD(iDimStoch)-1)*PRODUCT(nSGElems_nD(iDimStoch+1:PP_nDimStoch))
END DO
mySGElem_nD(PP_nDimStoch)=tmp+1

! Set root for local SGElem (proc with lowest rank in mesh)
MPI_COMM_SGELEM = MPI_COMM_NULL
CALL MPI_COMM_SPLIT(MPI_COMM_WORLD,mySGElem,myRank,MPI_COMM_SGELEM,iError)
CALL MPI_COMM_RANK(MPI_COMM_SGELEM,mySGElemRank,iError)
MPISGElemRoot = (mySGElemRank.EQ.0)

MPI_COMM_SGELEMROOT = MPI_COMM_NULL
MeshRoot = 0
IF(MPISGElemRoot) MeshRoot = 1
CALL MPI_COMM_SPLIT(MPI_COMM_WORLD,MeshRoot,myRank,MPI_COMM_SGELEMROOT,iError)
END SUBROUTINE StochBuildPartition



!==================================================================================================================================
!> Deallocate SG arrays
!==================================================================================================================================
SUBROUTINE FinalizeSG()
! MODULES
USE MOD_SG_Vars
IMPLICIT NONE
!==================================================================================================================================
SDEALLOCATE(xiDet)
SDEALLOCATE(xiQPRef)
SDEALLOCATE(xiQP)
SDEALLOCATE(xiQP_nDim)
SDEALLOCATE(xiQPRefExact)
SDEALLOCATE(xiQPExact)
SDEALLOCATE(xiQPExact_nDim)
SDEALLOCATE(wQP)
SDEALLOCATE(wQP_nDim)
SDEALLOCATE(wQPExact)
SDEALLOCATE(wQPExact_nDim)
SDEALLOCATE(C_nD)
SDEALLOCATE(C_nD_flux)
SDEALLOCATE(C_nD_idx)
SDEALLOCATE(C_1D_R3)
SDEALLOCATE(C_1D)
SDEALLOCATE(FullOrderToTensorVec)
SDEALLOCATE(SG_Vdm_QuadOrth)
SDEALLOCATE(SG_Vdm_OrthQuad)
SDEALLOCATE(SG_Vdm_QuadOrthExact)
SDEALLOCATE(SG_Vdm_OrthQuadExact)
SDEALLOCATE(StochSigma)
SDEALLOCATE(StochMu)
SDEALLOCATE(UniIntBounds)
SDEALLOCATE(nDimStochPositions)
SDEALLOCATE(StochVarNamesAll)
SDEALLOCATE(StochVarNames)
SDEALLOCATE(nSGElems_nD)
SDEALLOCATE(mySGElem_nD)
SDEALLOCATE(DimStochIsSet)
SDEALLOCATE(SG_Vdm_OrthQuad_1D)
SDEALLOCATE(SG_Vdm_QuadOrth_1D)
SGInitIsDone = .FALSE.
END SUBROUTINE FinalizeSG

END MODULE MOD_SG
