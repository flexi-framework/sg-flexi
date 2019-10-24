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
!> Contains all routines for setup of stochastic quadrature
!==================================================================================================================================
MODULE MOD_SG_Quadrature
! MODULES
USE MOD_SG_Vars
IMPLICIT NONE
PRIVATE
!----------------------------------------------------------------------------------------------------------------------------------

INTERFACE InitQuadrature
  MODULE PROCEDURE InitQuadrature
END INTERFACE

INTERFACE PsiOfXi
   MODULE PROCEDURE PsiOfXi
END INTERFACE

PUBLIC :: InitQuadrature
PUBLIC :: PsiOfXi
!==================================================================================================================================

CONTAINS



!===================================================================================================================================
!> Build everything related to stochastic quadrature points:
!> - Read in parametes & Allocate arrays
!> - allocate arrays, set pointers and call subroutines to:
!>     - Get quadrature nodes (absolute and reference) and weights, 1D as well as nD
!>     - Get VanDerMonde matrices to switch between orthonormal SG basis and quadrature points in the stochastic space
!===================================================================================================================================
SUBROUTINE InitQuadrature(posti)
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_Basis       ,ONLY: LegendrePolynomialAndDerivative,LegendreGaussNodesAndWeights
USE MOD_ReadInTools ,ONLY: GETINT,CountOption
USE MOD_SG_Tensor   ,ONLY: GetMultiIndFrom1DInd
USE MOD_StringTools ,ONLY: INTTOSTR
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
LOGICAL,INTENT(IN):: posti
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!==================================================================================================================================
SWRITE(UNIT_stdOut,'(A)') ' Initialize quadrature ...'

!standard quadrature
!--------------------------------------------------------------------------------
SWRITE(UNIT_stdOut,'(A)') ' ...standard part'

nQP = GETINT('nQP',INTTOSTR(PP_M+1))
IF (nQP.LE.PP_M) THEN
  nQP=PP_M+1
  CALL PrintWarning("nQP was set lower than or equal to PP_M. Is now set to PP_M+1")
END IF
IF(posti) nQP = MERGE(nQPVisu,nQP,nQPVisu.GT.0)

CALL GetQuadNodesAndWeightsAndVdm(nQP,nQPTotal,&
                                  xiQP, xiQPRef, xiQP_nDim, &
                                  wQP, wQP_nDim, &
                                  SG_Vdm_OrthQuad, SG_Vdm_QuadOrth )

IF(posti) RETURN

!exact quadrature
!--------------------------------------------------------------------------------
SWRITE(UNIT_stdOut,'(A)') ' ...exact part'

nQPExact = GETINT('nQPExact',INTTOSTR(4*nQP))
IF (nQPExact.LT.nQP) THEN
  nQPExact=nQP
  CALL PrintWarning("nQPExact was set lower than nQP. Is now set to nQP")
END IF

CALL GetQuadNodesAndWeightsAndVdm(nQPExact,nQPExactTotal,&
                                  xiQPExact, xiQPRefExact, xiQPExact_nDim, &
                                  wQPExact, wQPExact_nDim, &
                                  SG_Vdm_OrthQuadExact, SG_Vdm_QuadOrthExact )

END SUBROUTINE InitQuadrature



!===================================================================================================================================
!> For normal, or exact quadrature:
!> - Get quadrature nodes (absolute and reference) and weights, 1D as well as nD
!> - Get VanDerMonde matrices to switch between orthonormal SG basis and quadrature points in the stochastic space
!===================================================================================================================================
SUBROUTINE GetQuadNodesAndWeightsAndVdm(nQPLoc, nQPTotalLoc, &
                                        xiQPLoc, xiQPRefLoc, xiQP_nDimLoc,    &
                                        wQPLoc, wQP_nDimLoc, &
                                        Vdm_OrthQuadLoc, Vdm_QuadOrthLoc )
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_Basis       ,ONLY: LegendreGaussNodesAndWeights
USE MOD_ReadInTools ,ONLY: GETINT,CountOption
USE MOD_SG_Tensor   ,ONLY: GetMultiIndFrom1DInd
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER,INTENT(INOUT)          :: nQPLoc
INTEGER,INTENT(INOUT)          :: nQPTotalLoc
REAL,ALLOCATABLE,INTENT(INOUT) :: xiQPLoc(:,:)
REAL,ALLOCATABLE,INTENT(INOUT) :: xiQPRefLoc(:)
REAL,ALLOCATABLE,INTENT(INOUT) :: xiQP_nDimLoc(:,:)
REAL,ALLOCATABLE,INTENT(INOUT) :: wQPLoc(:)
REAL,ALLOCATABLE,INTENT(INOUT) :: wQP_nDimLoc(:)
REAL,ALLOCATABLE,INTENT(INOUT) :: Vdm_OrthQuadLoc(:,:)
REAL,ALLOCATABLE,INTENT(INOUT) :: Vdm_QuadOrthLoc(:,:)
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER              :: iQP,iCoef,iDimStoch,PolyDeg
REAL,ALLOCATABLE     :: Vdm_OrthQuad_1D(:,:)
INTEGER,ALLOCATABLE  :: iQPVec(:,:)
!==================================================================================================================================
nQPTotalLoc = nQPLoc**PP_nDimStoch
!allocate arrays
ALLOCATE(xiQPLoc        (1:PP_nDimStoch,1:nQPLoc     ))
ALLOCATE(xiQPRefLoc     (               1:nQPLoc     ))
ALLOCATE(xiQP_nDimLoc   (1:PP_nDimStoch,1:nQPTotalLoc))
ALLOCATE(wQPLoc         (               1:nQPLoc     ))
ALLOCATE(wQP_nDimLoc    (               1:nQPTotalLoc))
ALLOCATE(Vdm_OrthQuadLoc(1:nQPTotalLoc ,0:PP_nCoefM1 ))
ALLOCATE(Vdm_QuadOrthLoc(0:PP_nCoefM1  ,1:nQPTotalLoc))

!Get 1D nodes and weights in reference space, as well as nD nodes in absolute coordinates
SELECT CASE(Distribution)
CASE(SG_DIST_NORMAL)
  CALL GetHermiteRootsAndWeights(nQPLoc,wQPLoc,xiQPRefLoc)
  DO iDimStoch=1,PP_nDimStoch
    ! Save xi evaluated at quadrature points
    xiQPLoc(iDimStoch,:) = StochMu(iDimStoch)+StochSigma(iDimStoch)*xiQPRefLoc
  END DO
CASE(SG_DIST_UNIFORM)
  CALL LegendreGaussNodesAndWeights(nQPLoc-1,xiQPRefLoc(1:nQPLoc),wQPLoc(1:nQPLoc))
  wQPLoc=0.5*wQPLoc !adaptation for pdf=0.5
  DO iDimStoch=1,PP_nDimStoch
    ! Save xi evaluated at quadrature points
    xiQPLoc(iDimStoch,:) = UniIntBounds(1,iDimStoch)+0.5*(xiQPRefLoc(:)+1.)*(UniIntBounds(2,iDimStoch)-UniIntBounds(1,iDimStoch))
  END DO
END SELECT

! Mapping from 1...nQPTotal to 1..nQP in every dimension (tensor product to 1D)
CALL GetMultiIndFrom1DInd(nQPLoc,nQPTotalLoc,.FALSE.,iQPVec)

! Get MultiD nodes and weight from 1D nodes and weights
wQP_nDimLoc = 1.
Do iQP=1,nQPTotalLoc
  DO iDimStoch=1,PP_nDimStoch
    xiQP_nDimLoc(iDimStoch,iQP) = xiQPLoc(iDimStoch,iQPVec(iDimStoch,iQP))
    wQP_nDimLoc(iQP) = wQP_nDimLoc(iQP) * wQPLoc(iQPVec(iDimStoch,iQP))
  END DO
END DO

!Get VdM's from orthogonal basis to quadrature points and back
ALLOCATE(Vdm_OrthQuad_1D(1:nQPLoc,0:PP_M))

!Get Vdm 1D
DO iCoef=0,PP_M
  DO iQP=1,nQPLoc
    Vdm_OrthQuad_1D(iQP,iCoef)=PsiOfXi(iCoef,xiQPRefLoc(iQP),Distribution)
  END DO
END DO

!Get Vdm at quatrature point as product of orth polynomials in each dimension
Vdm_OrthQuadLoc=1.
DO iCoef=0,PP_nCoefM1
  DO iQP=1,nQPTotalLoc
    DO iDimStoch=1,PP_nDimStoch
      PolyDeg=FullOrderToTensorVec(iDimStoch,iCoef)
      Vdm_OrthQuadLoc(iQP,iCoef)=Vdm_OrthQuadLoc(iQP,iCoef)*Vdm_OrthQuad_1D(iQPVec(iDimStoch,iQP),PolyDeg)
    END DO
  END DO
END DO

DO iQP=1,nQPTotalLoc
  ! Calculate "inverse" Vandermonde matrix
  Vdm_QuadOrthLoc(:,iQP) = Vdm_OrthQuadLoc(iQP,:)*wQP_nDimLoc(iQP)
END DO
END SUBROUTINE GetQuadNodesAndWeightsAndVdm



!===================================================================================================================================
!> Compute Hermite roots (used as quadrature points) and quadrature weights.
!===================================================================================================================================
SUBROUTINE GetHermiteRootsAndWeights(nQPLoc,wQPLoc,xiQPRefLoc)
! MODULES
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER,INTENT(IN)      :: nQPLoc
REAL,INTENT(INOUT)      :: wQPLoc(nQPLoc)
REAL,INTENT(INOUT)      :: xiQPRefLoc(nQPLoc)
!----------------------------------------------------------------------------------------------------------------------------------
! External procedures defined in LAPACK
EXTERNAL DSYEV
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                                              :: h_temp
DOUBLE PRECISION,DIMENSION(1:nQPLoc,1:nQPLoc)     :: ChaosMatrix
INTEGER,PARAMETER                                 :: lwmax=1000
INTEGER                                           :: info
INTEGER                                           :: p
DOUBLE PRECISION,DIMENSION(lwmax)                 :: work
DOUBLE PRECISION,DIMENSION(1:nQPLoc)              :: Roots
!==================================================================================================================================
ChaosMatrix=0.
DO p=1,nQPLoc-1
  ChaosMatrix(p,p+1)=SQRT(REAL(p))
END DO

CALL DSYEV('N','U',nQPLoc,ChaosMatrix,nQPLoc,Roots,work,lwmax,info)
xiQPRefLoc=SNGL(Roots)
DO p=1,nQPLoc
  h_temp=PsiOfXi(nQPLoc-1,xiQPRefLoc(p),SG_DIST_NORMAL)
  wQPLoc(p)=1./(nQPLoc*h_temp*h_temp)
END DO
END SUBROUTINE GetHermiteRootsAndWeights



!===================================================================================================================================
!> evaluate Hermite polynomial of degree pStoch at Xi
!===================================================================================================================================
FUNCTION PsiOfXi(pStoch,Xi,Distribution)
! MODULES
USE MOD_Basis       ,ONLY: LegendrePolynomialAndDerivative
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER,INTENT(IN) :: pStoch  !< (IN)  input polynomial degree
REAL,INTENT(IN)    :: Xi
INTEGER,INTENT(IN) :: Distribution
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER            :: p
REAL               :: PsiOfXi
REAL               :: dummy
REAL,ALLOCATABLE   :: HermiteCoeff(:,:)     !< coefficients of the Hermite polynomials
!==================================================================================================================================
SELECT CASE(Distribution)
CASE(SG_DIST_NORMAL)
  !Get Hermite Coefficients
  !TODO: Do we need the whole matrix?
  ALLOCATE(HermiteCoeff(0:MAX(1,pStoch),0:MAX(1,pStoch)))
  HermiteCoeff(:,:)=0.
  HermiteCoeff(0,0)=1.
  HermiteCoeff(1,1)=1.
  DO p=2,pStoch
    HermiteCoeff(1:p,p)=HermiteCoeff(0:p-1,p-1)/SQRT(REAL(p))
    HermiteCoeff(0:p-2,p)=HermiteCoeff(0:p-2,p)-HermiteCoeff(0:p-2,p-2)*SQRT(REAL(p-1)/REAL(p))
  END DO
  ! P_n(x) = SUM over i of (HermiteCoeff(n,i) * x^i)
  PsiOfXi=0.
  DO p=0,pStoch
    PsiOfXi=PsiOfXi+HermiteCoeff(p,pStoch)*Xi**p
  END DO
CASE(SG_DIST_UNIFORM)
  CALL LegendrePolynomialAndDerivative(pStoch,Xi,PsiOfXi,dummy)
  PsiOfXi=SQRT(2.)*PsiOfXi
END SELECT
END FUNCTION PsiOfXi


END MODULE MOD_SG_Quadrature
