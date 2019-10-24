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
!> Contains all routines for setup of SG triple product tensor C.
!==================================================================================================================================
MODULE MOD_SG_Tensor
! MODULES
IMPLICIT NONE
PRIVATE
!----------------------------------------------------------------------------------------------------------------------------------
INTERFACE CreateC
  MODULE PROCEDURE CreateC
END INTERFACE

INTERFACE GetMultiIndFrom1DInd
  MODULE PROCEDURE GetMultiIndFrom1DInd
END INTERFACE

PUBLIC :: CreateC
PUBLIC :: GetMultiIndFrom1DInd
!==================================================================================================================================

CONTAINS

!===================================================================================================================================
!> Fills the Hermite or Legendre polynomial chaos Tensor C_pqr = <P_p*P_q,P_r> for several stochastic dimensions
!> - calls subroutien to build tensor for stoch. 1D
!> - counts nD entries and allocates tensor
!> - fills entries of nD (rank 1) tensor by multiplying entries of (rank 3) stoch. 1D tensor
!===================================================================================================================================
SUBROUTINE CreateC()
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_SG_Vars ,ONLY: C_nD,C_nD_flux,C_nD_idx,SizeC_nD,C_1D_R3,C_1D
USE MOD_SG_Vars ,ONLY: FullOrderToTensorVec
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER,DIMENSION(1:PP_nDimStoch) :: iND,jND,kND
INTEGER                           :: i,j,k,t,s,counter
LOGICAL                           :: nnz_flag
!==================================================================================================================================
SWRITE(UNIT_stdOut,'(A)') ' Create triple product tensor C ...'
! Build 1D C_1D tensor
CALL CreateC_1D()
IF(PP_nDimStoch.EQ.1) RETURN

!----------------------------------------------------------------------------------------------------
! All that follows:
! Get nD Tensor from 1D Tensor

DEALLOCATE(C_1D)

! Calculate number of nnz entries in multi-dim C tensor
!----------------------------------------------------------------------------------------------------
counter = 0
iND = 0
jND = 0
kND = 0
! Calculate number of nnz entries in multi-dim TPM
DO i = 0,PP_nCoefM1  !u
  iND(:) = FullOrderToTensorVec(:,i)
  DO j = 0,i  !v
    jND(:) = FullOrderToTensorVec(:,j)
    DO k = 0,j  !w
      ! Get indices of 1D tensors, that have to be multiplied
      kND(:) = FullOrderToTensorVec(:,k)
#ifdef MANUAL_UNROLL
      ! no need to store C_ijk = 1*1*...*1
      ! C_ijk is 1 in one dimension if i,j, or k is zero
      ! multi dim C_ijk is 1 if 1d C_ijk=1 in all dimensions
      IF (ALL((iND.EQ.0).OR.(jND.EQ.0).OR.(kND.EQ.0))) CYCLE
#endif

      ! Check whether entry is non-zero or symetric
      nnz_flag = .TRUE.
      DO t=1,PP_nDimStoch
        s = iND(t)+jND(t)+kND(t)
        IF ((MODULO(s,2).EQ.1) .OR. (MAX(iND(t),jND(t),kND(t)).GT.s/2)) THEN
          nnz_flag = .FALSE.
          EXIT
        END IF
      END DO

      ! Increase counter if element is non-zero
      IF(nnz_flag) counter = counter+1
    END DO
  END DO
END DO


! Fill nD Tensor entries (as Product of 1D Tensor entries)
!----------------------------------------------------------------------------------------------------

! Allocate sparse multi-dim TPM and corresponding index vector
ALLOCATE(C_nD_idx(3,1:counter))
ALLOCATE(C_nD(1:counter))
ALLOCATE(C_nD_flux(1:counter))
C_nD_idx(:,:) = 0
C_nD_flux(:) = 0.
C_nD(:) = 1.
SizeC_nD = counter

! Init counter variables for sparse representation
counter = 1

! Loop over all (symmetric) full-order elements
DO i = 0,PP_nCoefM1  !u
  iND(:) = FullOrderToTensorVec(:,i)
  ! Update element counter
  DO j=0,i !v
    jND(:) = FullOrderToTensorVec(:,j)
    DO k=0,j !w
      ! Get Indices of 1D tensors, that have to be multiplied
      kND(:) = FullOrderToTensorVec(:,k)
#ifdef MANUAL_UNROLL
      ! no need to store 1*1*1
      ! is 1 in a dimension if i,j, or k is zero
      ! is 1 if 1 in all dimensions
      IF (ALL((iND.EQ.0).OR.(jND.EQ.0).OR.(kND.EQ.0))) CYCLE
#endif

      ! Check whether entry is non-zero or symmetric
      nnz_flag = .TRUE.
      DO t=1,PP_nDimStoch
        s = iND(t)+jND(t)+kND(t)
        IF ((MODULO(s,2).EQ.1) .OR. (MAX(iND(t),jND(t),kND(t)).GT.s/2)) THEN
          nnz_flag = .FALSE.
          EXIT
        END IF
      END DO

      ! Calculate nnz entry and index if element is non-zero
      IF (nnz_flag) THEN
        DO t=1,PP_nDimStoch
          C_nD(counter) = C_nD(counter) * C_1D_R3(iND(t),jND(t),kND(t))
        END DO

#ifndef MANUAL_UNROLL
        ! Build TPM_flux to improve performance in flux calculation
        IF ((i.EQ.j).AND.(j.EQ.k)) THEN ! three identical indices
          C_nD_flux(counter) = C_nD(counter)/6
        ELSEIF ((i.EQ.j).OR.(j.EQ.k).OR.(i.EQ.k)) THEN ! two identical indices
          C_nD_flux(counter) = C_nD(counter)/2
        ELSE
          C_nD_flux(counter) = C_nD(counter)
        END IF
#endif

        ! Save indices corresponding to C_nD entry
        C_nD_idx(1,counter) = i
        C_nD_idx(2,counter) = j
        C_nD_idx(3,counter) = k

        ! Increase Index of sparse representation
        counter = counter+1
      END IF
    END DO
  END DO
END DO

DEALLOCATE(C_1D_R3)

END SUBROUTINE CreateC


!===================================================================================================================================
!> Fills the rank 3 sparse Hermite or Legendre polynomial chaos tensor C_pqr = <P_p*P_q,P_r> for one stochastic dimension.
!> Recursive algorithms are used to keep numbers small (see stochastic Flexi code documentation for details)
!> Assumes normalized polynomials so that C_ppp = 1
!> Also gets rank 1 version with no zero entries by looping over non-zero entries (details in SG flexi code documentation)
!===================================================================================================================================
SUBROUTINE CreateC_1D()
! MODULES
USE MOD_PreProc
USE MOD_SG_Vars ,ONLY: Distribution,C_1D_R3,C_1D
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER :: p,q,r,s,sqr,dqr,i
#ifndef MANUAL_UNROLL
INTEGER :: j,k
#endif
REAL    :: rp,rq,rr
!==================================================================================================================================

!--------------------------------------------------------------------------
! create sparse rank 3 tensor
!--------------------------------------------------------------------------
ALLOCATE(C_1D_R3(0:PP_M,0:PP_M,0:PP_M))

IF (Distribution.EQ.SG_DIST_NORMAL) THEN
  C_1D_R3=0.
  DO p=0,PP_M
    DO dqr=-p,p,2
      q=(p+dqr)/2
      r=(p-dqr)/2
      IF (ABS(dqr).EQ.p) THEN
        C_1D_R3(p,q,r)=1.
      ELSE
        C_1D_R3(p,q,r)=C_1D_R3(p,q-1,r+1)*SQRT(REAL(q)/REAL(r+1))*REAL(p+r-q+2)/REAL(p+q-r)
      END IF
      DO sqr=p+2,2*PP_M-ABS(dqr),2
        q=(sqr+dqr)/2
        r=(sqr-dqr)/2
        C_1D_R3(p,q,r)=C_1D_R3(p,q-1,r-1)*2.*SQRT(REAL(q*r))/REAL(q+r-p)
      END DO
    END DO
  END DO
ELSE IF (Distribution.EQ.SG_DIST_UNIFORM) THEN
  C_1D_R3=0.
  DO p=0,PP_M
    DO q=0,PP_M
      DO r=0,PP_M
        s=(p+q+r)/2
        IF (((MODULO(p+q+r,2).EQ.1) .OR. (ABS(p-q).GT.r) .OR. (r.GT.(p+q)))) THEN
          C_1D_R3(p,q,r) = 0.
        ELSE
          rp=REAL(p)
          rq=REAL(q)
          rr=REAL(r)
          C_1D_R3(p,q,r)=SQRT((2.*rp+1.)*(2.*rq+1.)*(2.*rr+1.))/(rp+rq+rr+1.)&
                    *(AForLege(s-p)*AForLege(s-q)*AForLege(s-r)/AForLege(s))**2
        END IF
      END DO
    END DO
  END DO
END IF !distribution


!--------------------------------------------------------------------------
! create rank one version of tensor from sparse rank 3 tensor
!--------------------------------------------------------------------------

!get number of non-zero entries, then allocate and fill 1D-array with them

#ifdef MANUAL_UNROLL
!--------------------------------------------------------------------------
i=0
DO p=0,PP_M; DO q=0,p; DO r=p-q,q,2
  IF (r.GT.0) i=i+1
END DO; END DO; END DO

ALLOCATE(C_1D(i))

i=0
DO p=0,PP_M
  DO q=0,p
    DO r=p-q,q,2
      IF(r.EQ.0) CYCLE
      i=i+1
      C_1D(i)=C_1D_R3(p,q,r)
    END DO
  END DO
END DO
#else
!--------------------------------------------------------------------------
p=0
DO i=0,PP_M
  DO j=0,i-1
    DO k=i-j,j-1,2
      p=p+1
    END DO
  END DO
  DO j=0,i-1,2
      p=p+1
  END DO
END DO
DO i=0,PP_M,2
  DO j=0,i-1
      p=p+1
  END DO
      p=p+1
END DO

ALLOCATE(C_1D(p))

i=0
DO p=0,PP_M
  DO q=0,p-1
    DO r=p-q,q-1,2
      i=i+1
      C_1D(i)=C_1D_R3(p,q,r)
    END DO
  END DO
  DO q=0,p-1,2
    i=i+1
    C_1D(i)=C_1D_R3(p,p,q)
  END DO
END DO
DO p=0,PP_M,2
  DO q=0,p-1
    i=i+1
    C_1D(i)=C_1D_R3(p,q,q)
  END DO
  i=i+1
  C_1D(i)=C_1D_R3(p,p,p)
END DO
#endif
!--------------------------------------------------------------------------
END SUBROUTINE CreateC_1D



!===================================================================================================================================
!> Evaluates function A(p) for Legendre polynomial chaos calculation
!===================================================================================================================================
RECURSIVE FUNCTION AForLege(p)  RESULT(Resu)
! MODULES
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
REAL                :: Resu
INTEGER, INTENT(IN) :: p
INTEGER             :: q
! LOCAL VARIABLES
!==================================================================================================================================
IF (p.LT.0) THEN
  Resu = 0.
ELSE IF (p.EQ.0) THEN
  Resu = 1.
ELSE
  Resu=1.
  DO q=1,p
    Resu = Resu*SQRT((2.*REAL(q)-1.)/REAL(q))
  END DO
END IF
END FUNCTION AForLege



!===================================================================================================================================
!> Creates an index vector that gives multi D indices for tensor product or full order 1D indices.
!> Full order is used to get order of polynomial in every direction for every full order basis coefficient
!> Tensor product variant is used to get iQP in every direction from global iQP (ranging from 1 to nQPTotal)
!> Convention: full order starts at 0, tensor product starts at 1 (both 1D and nD)
!>
!> Example (PP_nDimStoch=2, PP_M=1, nQP=2):
!>
!>          full order                   tensor product
!>      1D index | nD index            1D index | nD index
!>          0    |   0 0                   1    |   1 1
!>          1    |   0 1                   2    |   1 2
!>          2    |   1 0                   3    |   2 1
!>                                         4    |   2 2
!>
!===================================================================================================================================
SUBROUTINE GetMultiIndFrom1DInd(MaxInd1D,EndInd,isFullOrder,MultiIndArr)
! MODULES
USE MOD_PreProc
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER,INTENT(IN)              :: MaxInd1D
INTEGER,INTENT(IN)              :: EndInd
LOGICAL,INTENT(IN)              :: isFullOrder
INTEGER,ALLOCATABLE,INTENT(OUT) :: MultiIndArr(:,:)
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                         :: i,iDimStoch,StartInd,Ind1D
!==================================================================================================================================
StartInd= MERGE(0,1,isFullOrder)

ALLOCATE(MultiIndArr(PP_nDimStoch,StartInd:EndInd))
MultiIndArr(:,StartInd)= StartInd
DO i=StartInd+1,EndInd
  MultiIndArr(:,i)=MultiIndArr(:,i-1)
  ! increment vector:
  ! increment index of last dim, if this does not trigger break criterion;
  ! else, reset index for last dim and repeat check with next dim index
  DO iDimStoch=PP_nDimStoch,1,-1
    ! Current ind is reset if: sum of ind (full order) or current ind (tensor product) has reached threshold
    Ind1D = MERGE( SUM(MultiIndArr(:,i)),  MultiIndArr(iDimStoch,i) , isFullOrder)
    IF(Ind1D.EQ.MaxInd1D)THEN
      MultiIndArr(iDimStoch,i)=StartInd
    ELSE
      MultiIndArr(iDimStoch,i)=MultiIndArr(iDimStoch,i)+1
      EXIT
    END IF
  END DO
END DO
END SUBROUTINE GetMultiIndFrom1DInd

END MODULE MOD_SG_Tensor
