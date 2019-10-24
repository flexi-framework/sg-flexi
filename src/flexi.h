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
!===================================================================================================================================
! Here, preprocessor variables for different equation systems and abbreviations for specific expressions are defined
!===================================================================================================================================

! Abbrevations
#ifdef SUN
#  define __DATE__ '__TIME__ and __DATE__ not'
#  define __TIME__ 'available for SUN COMPILER'
#  define IEEE_ISNAN
#elif SX
#  define __DATE__ '__TIME__ and __DATE__ not'
#  define __TIME__ 'available for SX COMPILER'
#elif PGI
#  define NO_ISNAN
#endif
#ifndef __FILENAME__
#define __FILENAME__ __FILE__
#endif
#define __STAMP__ __FILENAME__,__LINE__,__DATE__,__TIME__

#ifdef GNU
#  define IEEE_IS_NAN ISNAN
#endif

#define SIZEOF_F(x) (STORAGE_SIZE(x)/8)

#ifdef GNU
#define CHECKSAFEINT(x,k)  IF(x>HUGE(1_  k).OR.x<-HUGE(1_  k))       CALL ABORT(__STAMP__,'Integer conversion failed: out of range!')
#define CHECKSAFEREAL(x,k) IF(x>HUGE(1._ k).OR.x<-HUGE(1._ k))       CALL ABORT(__STAMP__,'Real conversion failed: out of range!')
#elif CRAY
#define CHECKSAFEINT(x,k)
#define CHECKSAFEREAL(x,k)
#else
#define CHECKSAFEINT(x,k)  IF(x>HUGE(1_  ## k).OR.x<-HUGE(1_  ## k)) CALL ABORT(__STAMP__,'Integer conversion failed: out of range!')
#define CHECKSAFEREAL(x,k) IF(x>HUGE(1._ ## k).OR.x<-HUGE(1._ ## k)) CALL ABORT(__STAMP__,'Real conversion failed: out of range!')
#endif


! Test for equality: read description in mathtools.f90 for further infos
#define ALMOSTEQUALABSOLUTE(x,y,tol)  (ABS((x)-(y)).LE.(tol))
#define ALMOSTEQUALRELATIVE(x,y,tol)  (ABS((x)-(y)).LE.MAX(ABS(x),ABS(y))*(tol))
#define ALMOSTEQUALABSORREL(x,y,tol)  (ALMOSTEQUALABSOLUTE(x,y,tol) .OR.  ALMOSTEQUALRELATIVE(x,y,tol))
#define ALMOSTEQUALABSANDREL(x,y,tol) (ALMOSTEQUALABSOLUTE(x,y,tol) .AND. ALMOSTEQUALRELATIVE(x,y,tol))

! Define MPI specific write shortcuts
#if USE_MPI
#  define SWRITE IF(MPIRoot) WRITE
#  define IPWRITE(a,b) WRITE(a,b)myRank,
#  define GETTIME(a) a=MPI_WTIME()
#else
#  define SWRITE WRITE
#  define IPWRITE WRITE
#  define GETTIME(a) CALL CPU_TIME(a)
#endif
#define TWRITE WRITE
#define ERRWRITE(a,b) CALL CreateErrFile(); IF(ErrorFiles) WRITE(UNIT_errOut,b)
#define LOGWRITE(a,b) IF(Logging) WRITE(UNIT_logOut,b)
#define SDEALLOCATE(A) IF(ALLOCATED(A))  DEALLOCATE(A)
#define ADEALLOCATE(A) IF(ASSOCIATED(A)) DEALLOCATE(A)

! Loop variables
#define PP_IJK     i,j,k
#define PP_ij      i,j

! Predefined "PARAMETER-like" variables
#define XI_MINUS   5
#define XI_PLUS    3
#define ETA_MINUS  2
#define ETA_PLUS   4
#define ZETA_MINUS 1
#define ZETA_PLUS  6

! Entry position in SideToElem
#define S2E_ELEM_ID        1
#define S2E_NB_ELEM_ID     2
#define S2E_LOC_SIDE_ID    3
#define S2E_NB_LOC_SIDE_ID 4
#define S2E_FLIP           5

! Entry position in ElemToSide
#define E2S_SIDE_ID 1
#define E2S_FLIP    2

! Entry position in BC
#define MI_SIDEID 1
#define MI_FLIP   2

! Entry position in BC
#define BC_TYPE  1
#define BC_STATE 2
#define BC_ALPHA 3

! Entry position in BC
#define SEND 1
#define RECV 2

! Distribution
#define SG_DIST_NORMAL  1
#define SG_DIST_UNIFORM 2

!#define DEBUGMESH

#if !(FV_ENABLED)
#define FV_Elems(x) 0
#define FV_Elems_master(x) 0
#define FV_Elems_slave(x) 0
#define FV_Elems_Sum(x) 0
#endif

#define KILL(x) SWRITE(*,*) __FILE__,__LINE__,x; stop

! overintegration
#define CUTOFF 1
#define CUTOFFCONS 2

#define SMODE(iStoch) 1+iStoch:PP_nVar:PP_nCoef
#define PP_iVarDet  PP_nCoef*(iVarDet-1)+1:PP_nCoef*iVarDet
#define PP_iVarDetQ nQPTotal*(iVarDet-1)+1:nQPTotal*iVarDet
#define PP_iQPCons       iQP:PP_nVarQ:nQPTotal

#define SGI1             1:   PP_nCoef /* density */
#define SGI2    PP_nCoef+1: 2*PP_nCoef /* momentum x */
#define SGI3  2*PP_nCoef+1: 3*PP_nCoef /* momentum y */
#define SGI4  3*PP_nCoef+1: 4*PP_nCoef /* momentum z */
#define SGI5  4*PP_nCoef+1: 5*PP_nCoef /* momentum vector */

#define PP_nCoefM1 PP_nCoef-1

! Entry position in SG Dimensions
#define SGV_AMPLITUDE      1
#define SGV_WAVESPEED      2
#define SGV_WAVELENGTH     3
#define SGV_STARTPOS_XY    4
#define SGV_REFSTATEVEL_XY 5
#define SGV_DELTA_99       6
#define SGV_DMR_ANGLE      7
#define SGV_MACH           8

#define SGV_ARRAY_LEN  8 /*must be equal to maximal number of defined stoch dimensions*/


!2d functionality
#if (PP_dim==2)
#define ZDIM(a) 0
#define PP_NZ 0
#else
#define ZDIM(a) a
#define PP_NZ PP_N
#endif
