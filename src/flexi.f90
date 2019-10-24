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
!> Control program of the Flexi code. Initialization of the computation
!==================================================================================================================================
PROGRAM Flexi
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_Commandline_Arguments
USE MOD_Restart,           ONLY:DefineParametersRestart,InitRestart,Restart,FinalizeRestart
USE MOD_Interpolation,     ONLY:DefineParametersInterpolation,InitInterpolation,FinalizeInterpolation
USE MOD_Mesh,              ONLY:DefineParametersMesh,InitMesh,FinalizeMesh,FinalizePartition
USE MOD_Eos,               ONLY:DefineParametersEos
USE MOD_Exactfunc,         ONLY:DefineParametersExactFunc
USE MOD_Mortar,            ONLY:InitMortar,FinalizeMortar
USE MOD_Equation,          ONLY:DefineParametersEquation,InitEquation,FinalizeEquation
USE MOD_Testcase,          ONLY:DefineParametersTestcase
USE MOD_DG,                ONLY:InitDG,FinalizeDG
#if PARABOLIC
USE MOD_Lifting,           ONLY:DefineParametersLifting,InitLifting,FinalizeLifting
#endif /*PARABOLIC*/
USE MOD_Filter,            ONLY:DefineParametersFilter,InitFilter,FinalizeFilter
USE MOD_SG,                ONLY:DefineParametersSG,InitSG,FinalizeSG
#if SG_HP_LIMITER
USE MOD_HPLimiter,               ONLY:InitHPLimiter,FinalizeHPLimiter
#endif
USE MOD_Overintegration,   ONLY:DefineParametersOverintegration,InitOverintegration,FinalizeOverintegration
USE MOD_IO_HDF5,           ONLY:DefineParametersIO_HDF5,InitIOHDF5
USE MOD_Output,            ONLY:DefineParametersOutput,InitOutput,FinalizeOutput
USE MOD_Analyze,           ONLY:DefineParametersAnalyze,InitAnalyze,FinalizeAnalyze
USE MOD_RecordPoints,      ONLY:DefineParametersRecordPoints,InitRecordPoints,FinalizeRecordPoints
USE MOD_TimeDisc,          ONLY:DefineParametersTimedisc,InitTimeDisc,FinalizeTimeDisc,TimeDisc
USE MOD_MPI,               ONLY:DefineParametersMPI,InitMPI
#if USE_MPI
USE MOD_MPI,               ONLY:InitMPIvars,FinalizeMPI
#endif
USE MOD_Sponge,            ONLY:DefineParametersSponge,InitSponge,FinalizeSponge
#if FV_ENABLED
USE MOD_FV,                ONLY:DefineParametersFV,InitFV,FinalizeFV
USE MOD_FV_Basis,          ONLY:InitFV_Basis,FinalizeFV_Basis
#endif
USE MOD_Indicator,         ONLY:DefineParametersIndicator,InitIndicator,FinalizeIndicator
USE MOD_ReadInTools,       ONLY:prms,IgnoredParameters,PrintDefaultParameterFile,FinalizeParameters,ExtractParameterFile
#ifdef EDDYVISCOSITY
USE MOD_EddyVisc,          ONLY:DefineParametersEddyVisc
USE MOD_EddyVisc_Vars,     ONLY:FinalizeEddyViscosity
#endif
USE MOD_Restart_Vars      ,ONLY:RestartFile
USE MOD_StringTools       ,ONLY:STRICMP, GetFileExtension
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                    :: Time                              !< Used to measure simulation time
LOGICAL                 :: userblockFound
INTEGER                 :: res
!==================================================================================================================================
CALL SetStackSizeUnlimited()
CALL InitMPI()
CALL ParseCommandlineArguments()
! Check if the number of arguments is correct
IF (nArgs.GT.2) THEN
  ! Print out error message containing valid syntax
  CALL CollectiveStop(__STAMP__,'ERROR - Invalid syntax. Please use: flexi parameter.ini [restart.h5] or flexi --help'// &
  '[option/section name] to print help for a single parameter, parameter sections or all parameters.')
END IF
ParameterFile = Args(1)
IF (nArgs.GT.1) THEN
  RestartFile = Args(2)
ELSE IF (STRICMP(GetFileExtension(ParameterFile), "h5")) THEN
  ParameterFile = ".flexi.ini"
  CALL ExtractParameterFile(Args(1), ParameterFile, userblockFound)
  IF (.NOT.userblockFound) THEN
    CALL CollectiveStop(__STAMP__, "No userblock found in state file '"//TRIM(Args(1))//"'")
  END IF
  RestartFile = Args(1)
END IF
CALL DefineParametersMPI()
CALL DefineParametersIO_HDF5()
CALL DefineParametersInterpolation()
CALL DefineParametersSG()
CALL DefineParametersRestart()
CALL DefineParametersOutput()
CALL DefineParametersMesh()
CALL DefineParametersEos()
CALL DefineParametersEquation()
CALL DefineParametersExactFunc()
CALL DefineParametersTestcase()
CALL DefineParametersFilter()
CALL DefineParametersOverintegration()
CALL DefineParametersIndicator()
#if FV_ENABLED
CALL DefineParametersFV()
#endif
#if PARABOLIC
CALL DefineParametersLifting ()
#endif /*PARABOLIC*/
#ifdef EDDYVISCOSITY
CALL DefineParametersEddyVisc()
#endif /*EDDYVISCOSITY*/
CALL DefineParametersSponge()
CALL DefineParametersTimedisc()
CALL DefineParametersAnalyze()
CALL DefineParametersRecordPoints()

! check for command line argument --help or --markdown
IF (doPrintHelp.GT.0) THEN
  CALL PrintDefaultParameterFile(doPrintHelp.EQ.2, Args(1))
  STOP
END IF
CALL prms%read_options(ParameterFile)

CALL InitIOHDF5()
SWRITE(UNIT_stdOut,'(132("="))')
SWRITE(UNIT_stdOut,'(A)') &
"           __________________   _______              __________________   ______      ______   __________________ "
SWRITE(UNIT_stdOut,'(A)') &
"          /                 /) /      /)            /                 /) /      |   _/     /) /                 /)"
SWRITE(UNIT_stdOut,'(A)') &
"         /       __________// /      //            /      ___________// /__     |_/´     _// /____       ______// "
SWRITE(UNIT_stdOut,'(A)') &
"        /      /)__________) /      //            /      /)__________)  (__|          _/´_)  )___/      /)_____)  "
SWRITE(UNIT_stdOut,'(A)') &
"       /      //___         /      //            /      //___              |       _/´_/´       /      //         "
SWRITE(UNIT_stdOut,'(A)') &
"      /           /)       /      //            /           /)             |     /´ /´         /      //          "

SWRITE(UNIT_stdOut,'(A)') &
"     /      _____//       /      //            /      _____//            _/´     |/´          /      //           "
SWRITE(UNIT_stdOut,'(A)') &
"    /      /)____)       /      //            /      /)____)          _/´        |           /      //            "

SWRITE(UNIT_stdOut,'(A)') &
"   /      //            /      //_________   /      //_________   __/´     _     |__   _____/      //____         "
SWRITE(UNIT_stdOut,'(A)') &
"  /      //            /                 /) /                 /) /      _/´ |      /) /                 /)        "
SWRITE(UNIT_stdOut,'(A)') &
" /______//            /_________________// /_________________// /_____/` _/´|_____// /_________________//         "
SWRITE(UNIT_stdOut,'(A)') &
" )______)             )_________________)  )_________________)  )_____)/´   )_____)  )_________________)          "
SWRITE(UNIT_stdOut,'(A)')
SWRITE(UNIT_stdOut,'(132("="))')
! Measure init duration
StartTime=FLEXITIME()

! Initialization
CALL InitInterpolation()
#if FV_ENABLED
CALL InitFV_Basis()
#endif
CALL InitSG()
CALL InitMortar()
CALL InitOutput()
CALL InitMesh(meshMode=2)
#if SG_HP_LIMITER
CALL InitHPLimiter()
#endif
CALL InitRestart()
CALL InitFilter()
CALL InitOverintegration()
CALL InitIndicator()
#if USE_MPI
CALL InitMPIvars()
#endif
CALL InitEquation()
CALL InitDG()
#if FV_ENABLED
CALL InitFV()
#endif
#if PARABOLIC
CALL InitLifting()
#endif /*PARABOLIC*/
CALL InitSponge()
CALL InitTimeDisc()
CALL InitAnalyze()
CALL InitRecordpoints()
CALL IgnoredParameters()
CALL Restart()

! Measure init duration
Time=FLEXITIME()
SWRITE(UNIT_stdOut,'(132("="))')
SWRITE(UNIT_stdOut,'(A,F8.2,A)') ' INITIALIZATION DONE! [',Time-StartTime,' sec ]'
SWRITE(UNIT_stdOut,'(132("="))')

! Run Simulation
CALL TimeDisc()

!Finalize
CALL FinalizeOutput()
CALL FinalizeRecordPoints()
CALL FinalizeAnalyze()
#if PARABOLIC
CALL FinalizeLifting()
#endif /*PARABOLIC*/
#if FV_ENABLED
CALL FinalizeFV()
#endif
CALL FinalizeDG()
CALL FinalizeEquation()
CALL FinalizeInterpolation()
CALL FinalizeTimeDisc()
CALL FinalizeRestart()
CALL FinalizeMesh()
CALL FinalizePartition()
CALL FinalizeMortar()
CALL FinalizeSponge()
CALL FinalizeOverintegration()
CALL FinalizeFilter()
CALL FinalizeSG()
#if SG_HP_LIMITER
CALL FinalizeHPLimiter()
#endif
#if EDDYVISCOSITY
CALL FinalizeEddyViscosity
#endif
#if FV_ENABLED
CALL FinalizeFV_Basis()
#endif
CALL FinalizeIndicator()
! Measure simulation duration
Time=FLEXITIME()
CALL FinalizeParameters()
CALL FinalizeCommandlineArguments()
#if USE_MPI
CALL MPI_FINALIZE(iError)
IF(iError .NE. 0) STOP 'MPI finalize error'
CALL FinalizeMPI()
#endif
SWRITE(UNIT_stdOut,'(132("="))')
SWRITE(UNIT_stdOut,'(A,F8.2,A)') ' FLEXI FINISHED! [',Time-StartTime,' sec ]'
SWRITE(UNIT_stdOut,'(132("="))')
OPEN(unit=99, FILE="walltime.out",STATUS='unknown', ACTION='write')
write(99,*)  Time-StartTime
CLOSE(99,STATUS='KEEP', IOSTAT=res)
END PROGRAM Flexi
