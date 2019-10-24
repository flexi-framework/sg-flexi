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
#include "flexi.h"

!==================================================================================================================================
!> Basic routines performing an analysis of the solution valid for all testcases
!==================================================================================================================================
MODULE MOD_Analyze
! MODULES
IMPLICIT NONE
PRIVATE
!----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
INTERFACE InitAnalyze
  MODULE PROCEDURE InitAnalyze
END INTERFACE

INTERFACE Analyze
  MODULE PROCEDURE Analyze
END INTERFACE

INTERFACE FinalizeAnalyze
  MODULE PROCEDURE FinalizeAnalyze
END INTERFACE


PUBLIC:: Analyze, InitAnalyze, FinalizeAnalyze
!==================================================================================================================================

PUBLIC::DefineParametersAnalyze
CONTAINS

!==================================================================================================================================
!> Define parameters
!==================================================================================================================================
SUBROUTINE DefineParametersAnalyze()
! MODULES
USE MOD_ReadInTools ,ONLY: prms
USE MOD_AnalyzeEquation ,ONLY: DefineParametersAnalyzeEquation
IMPLICIT NONE
!==================================================================================================================================
CALL prms%SetSection("Analyze")
CALL prms%CreateLogicalOption('CalcErrorNorms' , "Set true to compute L2 and LInf error norms at analyze step.",&
                                                 '.TRUE.')
CALL prms%CreateLogicalOption('AnalyzeToFile',   "Set true to output result of error norms to a file (CalcErrorNorms=T)",&
                                                 '.FALSE.')
CALL prms%CreateRealOption(   'Analyze_dt',      "Specifies time intervall at which analysis routines are called.",&
                                                 '0.')
CALL prms%CreateIntOption(    'nWriteData' ,     "Intervall as multiple of Analyze_dt at which HDF5 files "//&
                                                 "(e.g. State,TimeAvg,Fluc) are written.",&
                                                 '1')
CALL prms%CreateIntOption(    'NAnalyze'   ,     "Polynomial degree at which analysis is performed (e.g. for L2 errors). "//&
                                                 "Default: 2*N.")
CALL prms%CreateIntOption(    'AnalyzeExactFunc',"Define exact function used for analyze (e.g. for computing L2 errors). "//&
                                                 "Default: Same as IniExactFunc")
CALL prms%CreateIntOption(    'AnalyzeRefState' ,"Define state used for analyze (e.g. for computing L2 errors). "//&
                                                 "Default: Same as IniRefState")
CALL prms%CreateLogicalOption('doMeasureFlops',  "Set true to measure flop count, if compiled with PAPI.",&
                                                 '.TRUE.')
CALL DefineParametersAnalyzeEquation()
END SUBROUTINE DefineParametersAnalyze

!==================================================================================================================================
!> Initializes variables necessary for analyze subroutines
!> - provides basic quantities like global domain volume, surface area of boundary conditions
!>   or precomputed surface and volume integration weights
!> - initializes other specific analysis and benchmarking routines
!==================================================================================================================================
SUBROUTINE InitAnalyze()
! MODULES
USE MOD_Globals
USE MOD_Preproc
USE MOD_Analyze_Vars
USE MOD_AnalyzeEquation,    ONLY: InitAnalyzeEquation
USE MOD_ReadInTools,        ONLY: GETINT,GETREAL,GETLOGICAL
USE MOD_StringTools,        ONLY: INTTOSTR
USE MOD_Interpolation_Vars, ONLY: xGP,wGP,wBary,InterpolationInitIsDone
USE MOD_Mesh_Vars,          ONLY: nBCs,SurfElem,nSides,AnalyzeSide,sJ,nElems
USE MOD_Equation_Vars,      ONLY: StrVarNamesDet,IniExactFunc,IniRefState
USE MOD_Output,             ONLY: InitOutputToFile
USE MOD_Output_Vars,        ONLY: ProjectName
USE MOD_Benchmarking,       ONLY: InitBenchmarking
USE MOD_Timedisc_Vars,      ONLY: TEnd
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                             :: lastLine(0:4*PP_nVarDet+5)
CHARACTER(MAXVAL(LEN_TRIM(StrVarNamesDet))+10) :: VarNames(4*PP_nVarDet+5)
INTEGER                          :: i,j,k,iSurf,iElem,iSide
!==================================================================================================================================
IF ((.NOT.InterpolationInitIsDone).OR.AnalyzeInitIsDone) THEN
  CALL CollectiveStop(__STAMP__,'InitAnalyse not ready to be called or already called.')
END IF
SWRITE(UNIT_StdOut,'(132("-"))')
SWRITE(UNIT_stdOut,'(A)') ' INIT ANALYZE...'

! Get the various analysis/output variables
doCalcErrorNorms  =GETLOGICAL('CalcErrorNorms' ,'.TRUE.')
doAnalyzeToFile   =GETLOGICAL('AnalyzeToFile','.FALSE.')
AnalyzeExactFunc  =GETINT('AnalyzeExactFunc',INTTOSTR(IniExactFunc))
AnalyzeRefState   =GETINT('AnalyzeRefState' ,INTTOSTR(IniRefState))

Analyze_dt        =GETREAL('Analyze_dt','0.0')
nWriteData        =GETINT('nWriteData' ,'1')
NAnalyze          =GETINT('NAnalyze'   ,INTTOSTR(2*(PP_N+1)))
#if PP_dim == 3
NAnalyzeZ         =NAnalyze
#else
NAnalyzeZ         =0
#endif
! If Analyze_dt is set to 0 (default) or to a negative value, no analyze calls should be performed at all.
! To achieve this, Analyze_dt is set to the final simulation time. This will prevent any calls of the analyze routine
! except at the beginning and the end of the simulation.
IF (Analyze_dt.LE.0.) THEN
  Analyze_dt = TEnd
  nWriteData = 1
END IF

WriteData_dt = Analyze_dt*nWriteData

! precompute integration weights
ALLOCATE(wGPSurf(0:PP_N,0:PP_NZ),wGPVol(0:PP_N,0:PP_N,0:PP_NZ))
#if PP_dim == 3
DO j=0,PP_N; DO i=0,PP_N
  wGPSurf(i,j)  = wGP(i)*wGP(j)
END DO; END DO
DO k=0,PP_N; DO j=0,PP_N; DO i=0,PP_N
  wGPVol(i,j,k) = wGP(i)*wGP(j)*wGP(k)
END DO; END DO; END DO
#else
DO i=0,PP_N
  wGPSurf(i,0)  = wGP(i)
END DO
DO j=0,PP_N; DO i=0,PP_N
  wGPVol(i,j,0) = wGP(i)*wGP(j)
END DO; END DO
#endif

! precompute volume of the domain
ALLOCATE(ElemVol(nElems))
ElemVol=0.
DO iElem=1,nElems
  DO k=0,PP_NZ; DO j=0,PP_N; DO i=0,PP_N
    ElemVol(iElem)=ElemVol(iElem)+wGPVol(i,j,k)/sJ(i,j,k,iElem,0)
  END DO; END DO; END DO !i,j,k
END DO ! iElem
Vol=SUM(ElemVol)


! compute surface of each boundary
ALLOCATE(Surf(nBCs))
Surf=0.
DO iSide=1,nSides
  iSurf=AnalyzeSide(iSide)
  IF(iSurf.EQ.0) CYCLE
  DO j=0,PP_NZ; DO i=0,PP_N
    Surf(iSurf)=Surf(iSurf)+wGPSurf(i,j)*SurfElem(i,j,0,iSide)
  END DO; END DO
END DO
#if USE_MPI
CALL MPI_ALLREDUCE(MPI_IN_PLACE,Vol ,1   ,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,iError)
CALL MPI_ALLREDUCE(MPI_IN_PLACE,Surf,nBCs,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,iError)
#endif /*USE_MPI*/

! Initialize eval routines
CALL InitAnalyzeBasis(PP_N,NAnalyze,xGP,wBary)

IF(doAnalyzeToFile)THEN
  VarNames(1          :  PP_nVarDet)   ='L2_Mean_'//StrVarNamesDet
  VarNames(PP_nVarDet+1  :2*PP_nVarDet)   ='LInf_Mean_'//StrVarNamesDet
  VarNames(2*PP_nVarDet+1         :  3*PP_nVarDet)   ='L2_Std_'//StrVarNamesDet
  VarNames(3*PP_nVarDet+1  :4*PP_nVarDet)   ='LInf_Std_'//StrVarNamesDet
  VarNames(4*PP_nVarDet+1:4*PP_nVarDet+5) = &
    [CHARACTER(9) :: "timesteps","t_CPU","DOF","Ncells","nProcs"] ! gfortran hates mixed length arrays
  FileName_ErrNorm='out.'//TRIM(ProjectName)
  CALL InitOutputToFile(FileName_ErrNorm,'Analyze',4*PP_nVarDet+5,VarNames,lastLine)
  iterRestart    =MAX(lastLine(4*PP_nVarDet+1),0.)
  calcTimeRestart=MAX(lastLine(4*PP_nVarDet+2),0.)
END IF
#ifdef SGdummy
CALL InitAnalyzeEquation()
#endif
CALL InitBenchmarking()

AnalyzeInitIsDone=.TRUE.
SWRITE(UNIT_StdOut,'(A,ES18.9)')' Volume of computational domain : ',Vol
SWRITE(UNIT_stdOut,'(A)')' INIT ANALYZE DONE!'
SWRITE(UNIT_StdOut,'(132("-"))')
END SUBROUTINE InitAnalyze



!==================================================================================================================================
!> Initializes variables necessary for analyse subroutines
!> - Builds Vandermonde to interpolate the solution onto a Gauss-Lobatto mesh at a higher polynomial degree
!> - Precomputes volume interpolation weights
!==================================================================================================================================
SUBROUTINE InitAnalyzeBasis(N_in,Nloc,xGP,wBary)
! MODULES
USE MOD_Globals
USE MOD_Preproc
USE MOD_Analyze_Vars, ONLY: wGPVolAnalyze,Vdm_GaussN_NAnalyze
USE MOD_Basis,        ONLY: InitializeVandermonde
USE MOD_Interpolation,ONLY: GetNodesAndWeights
USE MOD_Interpolation_Vars,ONLY: NodeTypeGL
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER,INTENT(IN)               :: N_in                  !< input polynomial degree
INTEGER,INTENT(IN)               :: Nloc                  !< polynomial degree of analysis polynomial
REAL,INTENT(IN)                  :: xGP(0:N_in)           !< interpolation points
REAL,INTENT(IN)                  :: wBary(0:N_in)         !< barycentric weights
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                             :: XiAnalyze(0:Nloc)
REAL                             :: wAnalyze( 0:Nloc)  ! GL integration weights used for the analyze
INTEGER                          :: i,j
#if PP_dim == 3
INTEGER                          :: k
#endif
!==================================================================================================================================
ALLOCATE(wGPVolAnalyze(0:Nloc,0:Nloc,0:ZDIM(Nloc)),Vdm_GaussN_NAnalyze(0:Nloc,0:N_in))
CALL GetNodesAndWeights(Nloc,NodeTypeGL,XiAnalyze,wAnalyze)
CALL InitializeVandermonde(N_in,Nloc,wBary,xGP,XiAnalyze,Vdm_GaussN_NAnalyze)

#if PP_dim == 3
DO k=0,Nloc; DO j=0,Nloc; DO i=0,Nloc
  wGPVolAnalyze(i,j,k) = wAnalyze(i)*wAnalyze(j)*wAnalyze(k)
END DO; END DO; END DO
#else
DO j=0,Nloc; DO i=0,Nloc
  wGPVolAnalyze(i,j,0) = wAnalyze(i)*wAnalyze(j)
END DO; END DO
#endif

END SUBROUTINE InitAnalyzeBasis



!==================================================================================================================================
!> Controls analysis routines and is called at analyze time levels
!> - calls generic error norm computation
!> - calls equation system specific analysis
!==================================================================================================================================
SUBROUTINE Analyze(Time,iter)
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_Analyze_Vars
USE MOD_AnalyzeEquation,    ONLY: AnalyzeEquation
USE MOD_Output,             ONLY: OutputToFile
USE MOD_Output_Vars,        ONLY: ProjectName
USE MOD_Mesh_Vars,          ONLY: nGlobalElems
USE MOD_Benchmarking,       ONLY: Benchmarking
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
REAL,INTENT(IN)                 :: Time                   !< current simulation time
INTEGER(KIND=8),INTENT(IN)      :: iter                   !< current iteration
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
CHARACTER(LEN=40)               :: formatStr
REAL                            :: CalcTime,RunTime
REAL                            :: L_Inf_Std_Error(PP_nVarDet),L_2_Std_Error(PP_nVarDet)
REAL                            :: L_Inf_Mean_Error(PP_nVarDet),L_2_Mean_Error(PP_nVarDet)
!==================================================================================================================================
! Graphical output
CalcTime=FLEXITIME()
RunTime=CalcTime-StartTime
SWRITE(UNIT_StdOut,'(A14,ES16.7)')' Sim time   : ',Time

! Calculate error norms
IF(doCalcErrorNorms)THEN
  CALL CalcErrorNorms(Time,L_Inf_Std_Error,L_2_Std_Error,L_2_Mean_Error,L_Inf_Mean_Error)
  IF(MPIroot) THEN
    WRITE(formatStr,'(A5,I1,A7)')'(A14,',PP_nVarDet,'ES18.9)'
    WRITE(UNIT_StdOut,formatStr)' L_2_Mean   : ',L_2_Mean_Error
    WRITE(UNIT_StdOut,formatStr)' L_2_Std    : ',L_2_Std_Error
    WRITE(UNIT_StdOut,formatStr)' L_inf_Mean : ',L_Inf_Mean_Error
    WRITE(UNIT_StdOut,formatStr)' L_inf_Std  : ',L_Inf_Std_Error
    IF(doAnalyzeToFile)THEN
      CALL OutputToFile(FileName_ErrNorm,(/Time/),(/4*PP_nVarDet+5,1/),  (/L_2_Mean_Error,L_2_Std_Error, &
                        L_Inf_Mean_Error,L_Inf_Std_Error, REAL(iter)+iterRestart,RunTime+CalcTimeRestart,&
                        REAL(nGlobalElems*(PP_N+1)**3),REAL(nGlobalElems),REAL(nProcessors)/))
    END IF
  END IF !MPIroot
END IF  ! ErrorNorms

CALL AnalyzeEquation(Time)
CALL Benchmarking()

IF(MPIroot .AND. (Time.GT.0.)) THEN
  WRITE(UNIT_StdOut,'(132("."))')
  WRITE(UNIT_stdOut,'(A,A,A,F8.2,A)') ' FLEXI RUNNING ',TRIM(ProjectName),'... [',RunTime,' sec ]'
  WRITE(UNIT_StdOut,'(132("-"))')
  WRITE(UNIT_StdOut,*)
END IF
END SUBROUTINE Analyze


!==================================================================================================================================
!> Calculates L_infinfity and L_2 error norms of state variables using the analyze framework (GL points+weights)
!==================================================================================================================================
SUBROUTINE CalcErrorNorms(Time,L_Inf_Std_Error,L_2_Std_Error,L_2_Mean_Error,L_Inf_Mean_Error)
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_Mesh_Vars,          ONLY: Elem_xGP,sJ,nElems
USE MOD_DG_Vars,            ONLY: U
USE MOD_Exactfunc,          ONLY: ExactFunc
USE MOD_ChangeBasisByDim,   ONLY: ChangeBasisVolume
USE MOD_Analyze_Vars,       ONLY: NAnalyze,NAnalyzeZ,Vdm_GaussN_NAnalyze
USE MOD_Analyze_Vars,       ONLY: wGPVolAnalyze,Vol,AnalyzeExactFunc,AnalyzeRefState
USE MOD_SG_Vars,            ONLY: nQPExactTotal,xiQPExact_nDim,wQPExact_nDim,xiDet
USE MOD_SG_Vars,            ONLY: nDimStochPositions
#if FV_ENABLED
USE MOD_FV_Vars,            ONLY: FV_Elems,FV_Vdm,FV_w
#endif
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
REAL,INTENT(IN)                 :: Time                   !< current simulation time
REAL,INTENT(OUT)                :: L_Inf_Std_Error(  PP_nVarDet)   !< L2 error of the solution
REAL,INTENT(OUT)                :: L_2_Std_Error(PP_nVarDet)   !< LInf error of the solution
REAL,INTENT(OUT)                :: L_Inf_Mean_Error(  PP_nVarDet)   !< L2 error of the solution
REAL,INTENT(OUT)                :: L_2_Mean_Error(PP_nVarDet)   !< LInf error of the solution
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                         :: iElem,i,j,k, iQP, iCoef,iDimStoch
REAL                            :: U_exactQuad(PP_nVarDet),U_exactStd(PP_nVarDet),U_exactMean(PP_nVarDet)
REAL                            :: UStd(PP_nVarDet)
REAL                            :: U_NAnalyze(1:PP_nVar,0:NAnalyze,0:NAnalyze,0:NAnalyzeZ)
REAL                            :: Coords_NAnalyze(3,0:NAnalyze,0:NAnalyze,0:NAnalyzeZ)
REAL                            :: J_NAnalyze(1,0:NAnalyze,0:NAnalyze,0:NAnalyzeZ)
REAL                            :: J_N(1,0:PP_N,0:PP_N,0:PP_NZ)
REAL                            :: wQP_analyze
REAL                            :: xi_analyze(1:SGV_ARRAY_LEN)
#if FV_ENABLED
REAL                            :: FV_w_volume, IntegrationWeight
REAL                            :: U_DG(PP_nVarDet,0:PP_N,0:PP_N,0:PP_NZ)
REAL                            :: U_FV(PP_nVarDet,0:PP_N,0:PP_N,0:PP_NZ)
#endif
!==================================================================================================================================
#if FV_ENABLED
FV_w_volume = FV_w**PP_dim
#endif
! Calculate error norms
L_Inf_Mean_Error(:)=-1.E10
L_Inf_Std_Error(:)=-1.E10
L_2_Mean_Error(:)=0.
L_2_Std_Error(:)=0.
! Interpolate values of Error-Grid from GP's
DO iElem=1,nElems
#if FV_ENABLED
  IF (FV_Elems(iElem).GT.0) THEN ! FV Element
    DO k=0,PP_NZ
      DO j=0,PP_N
        DO i=0,PP_N
          xi_analyze(:) = xiDet(:)
          DO iQP=1,nQPExactTotal
            ! Evaluate xi and integration weight at quadrature point
            DO iDimStoch=1,PP_nDimStoch
              xi_analyze(nDimStochPositions(iDimStoch)) = xiQPExact_nDim(iDimStoch,iQP)
            END DO
            CALL ExactFunc(AnalyzeExactFunc,time,Elem_xGP(1:3,i,j,k,iElem),U_DG(:,i,j,k),RefStateOpt=AnalyzeRefState, xi=xi_analyze)
          END DO
        END DO ! i
      END DO ! j
    END DO ! k
    CALL ChangeBasisVolume(PP_nVar,PP_N,PP_N,FV_Vdm,U_DG(:,:,:,:),U_FV(:,:,:,:))
    DO k=0,PP_NZ
      DO j=0,PP_N
        DO i=0,PP_N
          L_Inf_Mean_Error = MAX(L_Inf_Mean_Error,ABS(U(:,i,j,k,iElem) - U_FV(:,i,j,k)))
          IntegrationWeight = FV_w_volume/sJ(i,j,k,iElem,1)
          ! To sum over the elements, We compute here the square of the L_2 error
          L_2_Mean_Error = L_2_Mean_Error+(U(:,i,j,k,iElem) - U_FV(:,i,j,k))*(U(:,i,j,k,iElem) - U_FV(:,i,j,k))*IntegrationWeight
        END DO ! i
      END DO ! j
    END DO ! k
  ELSE
#endif
   ! Interpolate the physical position Elem_xGP to the analyze position, needed for exact function
   CALL ChangeBasisVolume(3,PP_N,NAnalyze,Vdm_GaussN_NAnalyze,Elem_xGP(1:3,:,:,:,iElem),Coords_NAnalyze(1:3,:,:,:))
   ! Interpolate the Jacobian to the analyze grid: be carefull we interpolate the inverse of the inverse of the jacobian ;-)
   J_N(1,0:PP_N,0:PP_N,0:PP_NZ)=1./sJ(:,:,:,iElem,0)
   CALL ChangeBasisVolume(1,PP_N,NAnalyze,Vdm_GaussN_NAnalyze,J_N,J_NAnalyze)
   ! Interpolate the solution to the analyze grid
   CALL ChangeBasisVolume(PP_nVar,PP_N,NAnalyze,Vdm_GaussN_NAnalyze,U(1:PP_nVar,:,:,:,iElem),U_NAnalyze(1:PP_nVar,:,:,:))
   DO k=0,NAnalyzeZ
     DO j=0,NAnalyze
       DO i=0,NAnalyze
         U_exactMean = 0.
         U_exactStd = 0.

         !DO iQP=1,nQPExact
         !  CALL ExactFunc(AnalyzeExactFunc,time,Coords_NAnalyze(1:3,i,j,k),U_exactQuad,RefStateOpt=AnalyzeRefState, xi=xiQPExact(iQP))
         !  U_exactMean = U_exactMean + wQPExact(iQP)*U_exactQuad
         !  U_exactStd = U_exactStd + wQPExact(iQP)*(U_exactQuad*U_exactQuad)
         !END DO
         !U_exactStd = SQRT(MAX(U_exactStd - U_exactMean*U_exactMean,TINY(1.)))
         !L_Inf_Mean_Error = MAX(L_Inf_Mean_Error,abs(U_NAnalyze(SMODE(0),i,j,k) - U_exactMean))
         !UStd=0.
         !DO iStoch=1,PP_M
         !  UStd=UStd+U_NAnalyze(SMODE(iStoch),i,j,k)**2
         !END DO

         xi_analyze(:) = xiDet(:)
         DO iQP=1,nQPExactTotal
           ! Evaluate xi and integration weight at quadrature point
           DO iDimStoch=1,PP_nDimStoch
             xi_analyze(nDimStochPositions(iDimStoch)) = xiQPExact_nDim(iDimStoch,iQP)
           END DO
           ! Evaluate solution for xi at quadrature points
           CALL ExactFunc(AnalyzeExactFunc,time,Coords_NAnalyze(1:3,i,j,k),U_exactQuad,RefStateOpt=AnalyzeRefState, xi=xi_analyze)
           U_exactMean = U_exactMean + wQPExact_nDim(iQP)* U_exactQuad
           U_exactStd  = U_exactStd  + wQPExact_nDim(iQP)*(U_exactQuad*U_exactQuad)
         END DO

         U_exactStd = SQRT(MAX(U_exactStd - U_exactMean*U_exactMean,TINY(1.)))
         L_Inf_Mean_Error = MAX(L_Inf_Mean_Error,ABS(U_NAnalyze(SMODE(0),i,j,k) - U_exactMean))
         UStd=0.
         DO iCoef=1,PP_nCoefM1
           UStd=UStd+U_NAnalyze(SMODE(iCoef),i,j,k)**2
         END DO

         UStd=SQRT(UStd)
         L_Inf_Std_Error = MAX(L_Inf_Std_Error,abs(UStd - U_exactStd))
         wQP_analyze = wGPVolAnalyze(i,j,k)*J_NAnalyze(1,i,j,k)
         ! To sum over the elements, We compute here the square of the L_2 error
         L_2_Mean_Error = L_2_Mean_Error+(U_NAnalyze(SMODE(0),i,j,k) - U_exactMean)**2*wQP_analyze
         L_2_Std_Error = L_2_Std_Error+(UStd - U_exactStd)**2*wQP_analyze
       END DO ! i
     END DO ! j
   END DO ! k
#if FV_ENABLED
  END IF
#endif
END DO ! iElem=1,nElems

#if USE_MPI
IF(MPIRoot)THEN
  CALL MPI_REDUCE(MPI_IN_PLACE    ,L_2_Mean_Error  ,PP_nVarDet,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,iError)
  CALL MPI_REDUCE(MPI_IN_PLACE    ,L_Inf_Mean_Error,PP_nVarDet,MPI_DOUBLE_PRECISION,MPI_MAX,0,MPI_COMM_WORLD,iError)
  CALL MPI_REDUCE(MPI_IN_PLACE    ,L_2_Std_Error   ,PP_nVarDet,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,iError)
  CALL MPI_REDUCE(MPI_IN_PLACE    ,L_Inf_Std_Error ,PP_nVarDet,MPI_DOUBLE_PRECISION,MPI_MAX,0,MPI_COMM_WORLD,iError)
ELSE
  CALL MPI_REDUCE(L_2_Mean_Error  ,0               ,PP_nVarDet,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,iError)
  CALL MPI_REDUCE(L_Inf_Mean_Error,0               ,PP_nVarDet,MPI_DOUBLE_PRECISION,MPI_MAX,0,MPI_COMM_WORLD,iError)
  CALL MPI_REDUCE(L_2_Std_Error   ,0               ,PP_nVarDet,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,iError)
  CALL MPI_REDUCE(L_Inf_Std_Error ,0               ,PP_nVarDet,MPI_DOUBLE_PRECISION,MPI_MAX,0,MPI_COMM_WORLD,iError)
END IF
#endif

! We normalize the L_2 Error with the Volume of the domain and take into account that we have to use the square root
L_2_Mean_Error = SQRT(L_2_Mean_Error/Vol)
L_2_Std_Error = SQRT(L_2_Std_Error/Vol)

END SUBROUTINE CalcErrorNorms


!==================================================================================================================================
!> Finalizes variables necessary for analyse subroutines
!==================================================================================================================================
SUBROUTINE FinalizeAnalyze()
! MODULES
USE MOD_AnalyzeEquation,   ONLY: FinalizeAnalyzeEquation
USE MOD_Analyze_Vars,      ONLY: AnalyzeInitIsDone,wGPSurf,wGPVol,Surf,wGPVolAnalyze,Vdm_GaussN_NAnalyze,ElemVol
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!==================================================================================================================================
CALL FinalizeAnalyzeEquation()
SDEALLOCATE(Vdm_GaussN_NAnalyze)
SDEALLOCATE(wGPVolAnalyze)
SDEALLOCATE(Surf)
SDEALLOCATE(wGPVol)
SDEALLOCATE(wGPSurf)
SDEALLOCATE(ElemVol)
AnalyzeInitIsDone = .FALSE.
END SUBROUTINE FinalizeAnalyze


END MODULE MOD_Analyze
