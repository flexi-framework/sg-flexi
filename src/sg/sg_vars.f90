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
!==================================================================================================================================
!> Contains global variables used by the SG module
!==================================================================================================================================
MODULE MOD_SG_Vars
! MODULES
IMPLICIT NONE
PUBLIC
SAVE

!----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
LOGICAL                      :: SGInitIsDone = .FALSE. !< SG routines have been initialized
INTEGER                      :: Distribution           !< distribution: normal or uniform

!Quadrature
!----------------------------------------------------------------------------------------------------

INTEGER                      :: nQP                    !< number of quadrature points in stochastic space per dimension
INTEGER                      :: nQPTotal               !< number of QP in entire spanned stochastic space (nQP**PP_nDimStoch)

REAL,ALLOCATABLE             :: xiQP(:,:)              !< quadrature points in the stoch. space (scaled with mu and sigma)
REAL,ALLOCATABLE             :: xiQP_nDim(:,:)         !<
REAL,ALLOCATABLE             :: xiQPRef(:)             !< quadrature points in the stoch. space (for mu=0 and standard sigma)

REAL,ALLOCATABLE             :: wQP(:)                 !< quadrature weights in the stochastic space
REAL,ALLOCATABLE             :: wQP_nDim(:)            !< quadrature weights in the stochastic space for higher stoch. dimension

REAL,ALLOCATABLE             :: SG_Vdm_OrthQuad(:,:)   !< Vdm for transform from orthonormal basis to quadrature points
REAL,ALLOCATABLE             :: SG_Vdm_QuadOrth(:,:)   !< Vdm for transform from quadrature points to orthonormal basis

! Exact Quadrature

INTEGER                      :: nQPExact               !< nQP for pseudo-exact numer. quadrature (e.g. for reference solution)
INTEGER                      :: nQPExactTotal          !< nQP for pseudo-exact numer. quadratur (nQPExact**PP-nDimStoch)

REAL,ALLOCATABLE             :: xiQPExact(:,:)         !< xiQP for pseudo-exact numer. quadrature for each stochastic dimension
REAL,ALLOCATABLE             :: xiQPExact_nDim(:,:)    !<
REAL,ALLOCATABLE             :: xiQPRefExact(:)        !< xiQPRef for pseudo-exact numer. quadrature

REAL,ALLOCATABLE             :: wQPExact(:)            !< wQP for pseudo-exact numer. quadrature
REAL,ALLOCATABLE             :: wQPExact_nDim(:)       !< wQP for pseudo-exact numer. quadrature for higher stoch. dimension

REAL,ALLOCATABLE             :: SG_Vdm_OrthQuadExact(:,:) !< Vdm for transform from orthonormal basis to quadrature points
REAL,ALLOCATABLE             :: SG_Vdm_QuadOrthExact(:,:) !< Vdm for transform from quadrature points to orthonormal basis

! other 

INTEGER                      :: nQPVisu                !< custom nQP read in from visu prm file

! Orth Basis and Tensor
!----------------------------------------------------------------------------------------------------

INTEGER,ALLOCATABLE          :: FullOrderToTensorVec(:,:) !< Vector to transform from full-order to tensor index
REAL,ALLOCATABLE             :: C_nD(:)                !< polynomial chaos tensor, sparse
REAL,ALLOCATABLE             :: C_nD_flux(:)           !< preconditioned C_nD to speed up sg products
REAL,ALLOCATABLE             :: C_1D_R3(:,:,:)         !< polynomial chaos tensor (rank 3: C_ijk) for 1 stochastic dimension
REAL,ALLOCATABLE             :: C_1D(:)                !< 1D index based version of C_nD (generally no zero entries)

INTEGER                      :: SizeC_nD               !< Size of C_nD tensor saved in variable for speed

INTEGER,ALLOCATABLE          :: C_nD_idx(:,:)          !< index vector for sparse representation of C_nD

! Operators
!----------------------------------------------------------------------------------------------------

LOGICAL                      :: ExactFraction          !< exact computation of fractions instead of numerical (solve Lin.Eq.S.)

REAL                         :: SG_SafetyFactor = 1.0  !< safety factor for imprecise numerical evaluation of min/max vealue of
                                                       !< wave speed over stoch. domain; used for time step and for LF riemann
! StochVars
!----------------------------------------------------------------------------------------------------

CHARACTER(LEN=255),ALLOCATABLE :: StochVarNamesAll(:)  !< array with names of all  stochastic variables
CHARACTER(LEN=255),ALLOCATABLE :: StochVarNames(:)     !< array with names of used stochastic variables
REAL,ALLOCATABLE             :: xiDet(:)               !< deterministic value of variables which could be chosen stochastic
REAL,ALLOCATABLE             :: StochSigma(:)          !< normal distribution: variance of random variable;
REAL,ALLOCATABLE             :: StochMu(:)             !< normal distribution: expectation of random variable
REAL,ALLOCATABLE             :: UniIntBounds(:,:)      !< uniform distribution: lower and upper bound of stoch. var. interval
INTEGER,ALLOCATABLE          :: nDimStochPositions(:)  !<

! Limiter
!----------------------------------------------------------------------------------------------------

REAL,ALLOCATABLE,TARGET      :: t_HPLimiter(:,:,:,:,:)
INTEGER,ALLOCATABLE          :: nLimitedElems(:)

! MultiElem
!----------------------------------------------------------------------------------------------------

INTEGER                      :: nSGElems               !< number of used stoch elements
INTEGER                      :: mySGElem               !< index of own stoch element
INTEGER,ALLOCATABLE          :: nSGElems_nD(:)         !< number of used stoch elements (multi d)
INTEGER,ALLOCATABLE          :: mySGElem_nD(:)         !< index of own stoch element (multi d)

!==================================================================================================================================
END MODULE MOD_SG_Vars
