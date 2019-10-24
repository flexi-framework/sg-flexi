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
!> Contains SG versions of various mathematical operators
!==================================================================================================================================
MODULE MOD_SG_Operators
! MODULES
IMPLICIT NONE
PRIVATE
!----------------------------------------------------------------------------------------------------------------------------------
ABSTRACT INTERFACE
  FUNCTION SG_Fraction_Interface(enum,denom)
    USE MOD_PreProc
    REAL,INTENT(IN),DIMENSION(0:PP_nCoefM1)          :: enum,denom
    REAL           ,DIMENSION(0:PP_nCoefM1)          :: SG_Fraction_Interface
  END FUNCTION
END INTERFACE
PROCEDURE(SG_Fraction_Interface),POINTER :: SG_Fraction !< pointer defining the fraction function

ABSTRACT INTERFACE
  FUNCTION SG_Inv_Interface(v)
    USE MOD_PreProc
    REAL,INTENT(IN),DIMENSION(0:PP_nCoefM1)          :: v
    REAL           ,DIMENSION(0:PP_nCoefM1)          :: SG_Inv_Interface
  END FUNCTION
END INTERFACE
PROCEDURE(SG_Inv_Interface),POINTER :: SG_Inv !< pointer defining the Inv function

ABSTRACT INTERFACE
  PURE FUNCTION  SG_Product_Interface(v1,v2)
    USE MOD_PreProc
    REAL,INTENT(IN),DIMENSION(0:PP_nCoefM1) :: v1,v2
    REAL           ,DIMENSION(0:PP_nCoefM1) :: SG_Product_Interface
  END FUNCTION
END INTERFACE
PROCEDURE(SG_Product_Interface),POINTER :: SG_Product !< pointer defining the product function

INTERFACE SG_DOT_PRODUCT_N
  MODULE PROCEDURE SG_DOT_PRODUCT_N
END INTERFACE

INTERFACE SG_DOT_PRODUCT
  MODULE PROCEDURE SG_DOT_PRODUCT
END INTERFACE

INTERFACE DefineParametersSGOperators
  MODULE PROCEDURE DefineParametersSGOperators
END INTERFACE

INTERFACE InitSGOperators
  MODULE PROCEDURE InitSGOperators
END INTERFACE


PUBLIC :: SG_Product
PUBLIC :: SG_Fraction
PUBLIC :: SG_Root
PUBLIC :: SG_Inv
PUBLIC :: SG_DOT_PRODUCT_N
PUBLIC :: SG_DOT_PRODUCT
PUBLIC :: DefineParametersSGOperators
PUBLIC :: InitSGOperators
!==================================================================================================================================

CONTAINS

!===================================================================================================================================
!> Define Parameters
!===================================================================================================================================
SUBROUTINE DefineParametersSGOperators()
! MODULES
USE MOD_ReadInTools ,ONLY: prms,addStrListEntry
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!==================================================================================================================================
CALL prms%CreateLogicalOption('ExactFraction' , "Exact computation of fractions")
END SUBROUTINE DefineParametersSGOperators



!===================================================================================================================================
!> All initializations related to the different SG operator routines
!===================================================================================================================================
SUBROUTINE InitSGOperators()
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_SG_Vars     ,ONLY: ExactFraction
USE MOD_ReadInTools ,ONLY: GETLOGICAL
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!==================================================================================================================================
SWRITE(UNIT_stdOut,'(A)') ' Initialize SG operators ...'
ExactFraction        = GETLOGICAL('ExactFraction','F')

IF(ExactFraction) THEN
  IF (PP_nDimStoch.EQ.1) THEN
    SG_Fraction => SG_Fraction_Exact_1Dim
    SG_Inv      => SG_Inv_Exact_1Dim
    SG_Product  => SG_Product_1Dim
  ELSE
    SG_Fraction => SG_Fraction_Exact_nDim
    SG_Inv      => SG_Inv_Exact_nDim
    SG_Product  => SG_Product_nDim
  END IF
ELSE
  SG_Fraction => SG_Fraction_Num
  SG_Inv      => SG_Inv_Num
  IF (PP_nDimStoch.EQ.1) THEN
    SG_Product  => SG_Product_1Dim
  ELSE
    SG_Product  => SG_Product_nDim
  END IF
END IF
END SUBROUTINE InitSGOperators



!===================================================================================================================================
!> SG_Inv = 1/v projected onto stochstic domain; exact evaulation
!===================================================================================================================================
FUNCTION SG_Inv_Exact_1Dim(v)
! MODULES
USE MOD_PreProc
USE MOD_SG_Vars , ONLY: C_1D_R3
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
REAL,INTENT(IN),DIMENSION(0:PP_nCoefM1)  :: v
REAL,DIMENSION(0:PP_nCoefM1)             :: SG_Inv_Exact_1Dim
!----------------------------------------------------------------------------------------------------------------------------------
! External procedures defined in LAPACK
EXTERNAL DSYSV
!----------------------------------------------------------------------------------------------------------------------------------

! LOCAL VARIABLES
REAL,DIMENSION(0:PP_nCoefM1,0:PP_nCoefM1)      :: Am
INTEGER                            :: iStoch,info
DOUBLE PRECISION,DIMENSION(100)    :: work
INTEGER,DIMENSION(0:PP_nCoefM1)          :: Pivot_dummy
!==================================================================================================================================
Am=0.
DO iStoch=0,PP_nCoefM1
  Am=Am+v(iStoch)*C_1D_R3(iStoch,:,:)
END DO
SG_Inv_Exact_1Dim(1:PP_nCoefM1)=0.
SG_Inv_Exact_1Dim(0           )=1.
CALL DSYSV('U',PP_nCoef,1,DBLE(Am),PP_nCoef,Pivot_dummy,SG_Inv_Exact_1Dim,PP_nCoef,work,100,info  )
END FUNCTION SG_Inv_Exact_1Dim


!===================================================================================================================================
!> Fraction enum/denom projected onto stochstic domain; exact evaulation
!===================================================================================================================================
FUNCTION SG_Inv_Exact_nDim(v)
! MODULES
USE MOD_PreProc
USE MOD_SG_Vars , ONLY: C_nD_flux,C_nD_idx,SizeC_nD
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
REAL,INTENT(IN),DIMENSION(0:PP_nCoefM1)  :: v
REAL,DIMENSION(0:PP_nCoefM1)             :: SG_Inv_Exact_nDim
!----------------------------------------------------------------------------------------------------------------------------------
! External procedures defined in LAPACK
EXTERNAL DSYSV
!----------------------------------------------------------------------------------------------------------------------------------

! LOCAL VARIABLES
REAL,DIMENSION(0:PP_nCoefM1,0:PP_nCoefM1)   :: Am
INTEGER                                     :: p,info
DOUBLE PRECISION,DIMENSION(100)             :: work
INTEGER,DIMENSION(0:PP_nCoefM1)             :: Pivot_dummy
!==================================================================================================================================
Am=0.

! Calculate upper triangle of symmetric matrix
DO p=1,SizeC_nD
  Am(C_nD_idx(3,p),C_nD_idx(1,p)) = Am(C_nD_idx(3,p),C_nD_idx(1,p)) + v(C_nD_idx(2,p))*C_nD_flux(p)
END DO
DO p=1,SizeC_nD
  Am(C_nD_idx(2,p),C_nD_idx(1,p)) = Am(C_nD_idx(2,p),C_nD_idx(1,p)) + v(C_nD_idx(3,p))*C_nD_flux(p)
END DO
DO p=1,SizeC_nD
  Am(C_nD_idx(3,p),C_nD_idx(2,p)) = Am(C_nD_idx(3,p),C_nD_idx(2,p)) + v(C_nD_idx(1,p))*C_nD_flux(p)
END DO

! Diagonal elements were only added once instead of twice -> double diagonal elements
DO p=0,PP_nCoefM1
  Am(p,p)=Am(p,p)*2.
END DO

SG_Inv_Exact_nDim(1:PP_nCoefM1)=0.
SG_Inv_Exact_nDim(0           )=1.
CALL DSYSV('U',PP_nCoef,1,DBLE(Am),PP_nCoef,Pivot_dummy,SG_Inv_Exact_nDim,PP_nCoef,work,100,info  )

END FUNCTION SG_Inv_Exact_nDim


!===================================================================================================================================
!> SG_Inv = 1/v projected onto stochstic domain; numerical evaulation
!===================================================================================================================================
PURE FUNCTION SG_Inv_Num(v)
! MODULES
USE MOD_PreProc
USE MOD_SG_Vars , ONLY: SG_VdM_OrthQuad,SG_VdM_QuadOrth
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
REAL,INTENT(IN),DIMENSION(0:PP_nCoefM1) :: v
REAL,DIMENSION(0:PP_nCoefM1)            :: SG_Inv_Num
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!==================================================================================================================================
SG_Inv_Num=MATMUL(SG_VdM_QuadOrth, 1./MATMUL(SG_VdM_OrthQuad,v))
END FUNCTION SG_Inv_Num

!===================================================================================================================================
!> SG_Root = sqrt(v) projected onto stochstic domain; numerical evaulation
!===================================================================================================================================
PURE FUNCTION SG_Root(v)
! MODULES
USE MOD_PreProc
USE MOD_SG_Vars , ONLY: SG_VdM_OrthQuad,SG_VdM_QuadOrth
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
REAL,INTENT(IN),DIMENSION(0:PP_M) :: v
REAL,DIMENSION(0:PP_M)            :: SG_Root
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!==================================================================================================================================
SG_Root=MATMUL(SG_VdM_QuadOrth, SQRT(MATMUL(SG_VdM_OrthQuad,v)))
END FUNCTION SG_Root

!===================================================================================================================================
!> Fraction enum/denom projected onto stochstic domain; exact evaulation
!===================================================================================================================================
FUNCTION SG_Fraction_Exact_1Dim(enum,denom)
! MODULES
USE MOD_PreProc
USE MOD_SG_Vars , ONLY: C_1D_R3
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
REAL,INTENT(IN),DIMENSION(0:PP_nCoefM1)  :: enum,denom
REAL,DIMENSION(0:PP_nCoefM1)             :: SG_Fraction_Exact_1Dim
!----------------------------------------------------------------------------------------------------------------------------------
! External procedures defined in LAPACK
EXTERNAL DSYSV
!----------------------------------------------------------------------------------------------------------------------------------

! LOCAL VARIABLES
REAL,DIMENSION(0:PP_nCoefM1,0:PP_nCoefM1)   :: Am
INTEGER                                     :: iCoef,info
DOUBLE PRECISION,DIMENSION(100)             :: work
INTEGER,DIMENSION(0:PP_nCoefM1)             :: Pivot_dummy
!==================================================================================================================================
Am=0.
DO iCoef=0,PP_nCoefM1
  Am=Am+denom(iCoef)*C_1D_R3(iCoef,:,:)
END DO
SG_Fraction_Exact_1Dim=enum
CALL DSYSV('U',PP_nCoef,1,DBLE(Am),PP_nCoef,Pivot_dummy,SG_Fraction_Exact_1Dim,PP_nCoef,work,100,info  )
END FUNCTION SG_Fraction_Exact_1Dim


!===================================================================================================================================
!> Fraction enum/denom projected onto stochstic domain; exact evaulation
!===================================================================================================================================
FUNCTION SG_Fraction_Exact_nDim(enum,denom)
! MODULES
USE MOD_PreProc
USE MOD_SG_Vars , ONLY: C_nD_flux,C_nD_idx,SizeC_nD
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
REAL,INTENT(IN),DIMENSION(0:PP_nCoefM1)  :: enum,denom
REAL,DIMENSION(0:PP_nCoefM1)             :: SG_Fraction_Exact_nDim
!----------------------------------------------------------------------------------------------------------------------------------
! External procedures defined in LAPACK
EXTERNAL DSYSV
!----------------------------------------------------------------------------------------------------------------------------------

! LOCAL VARIABLES
REAL,DIMENSION(0:PP_nCoefM1,0:PP_nCoefM1)   :: Am
INTEGER                                     :: p,info
DOUBLE PRECISION,DIMENSION(100)             :: work
INTEGER,DIMENSION(0:PP_nCoefM1)             :: Pivot_dummy
!==================================================================================================================================
Am=0.

! Calculate upper triangle of symmetric matrix
DO p=1,SizeC_nD
  Am(C_nD_idx(3,p),C_nD_idx(1,p)) = Am(C_nD_idx(3,p),C_nD_idx(1,p)) + denom(C_nD_idx(2,p))*C_nD_flux(p)
END DO
DO p=1,SizeC_nD
  Am(C_nD_idx(2,p),C_nD_idx(1,p)) = Am(C_nD_idx(2,p),C_nD_idx(1,p)) + denom(C_nD_idx(3,p))*C_nD_flux(p)
END DO
DO p=1,SizeC_nD
  Am(C_nD_idx(3,p),C_nD_idx(2,p)) = Am(C_nD_idx(3,p),C_nD_idx(2,p)) + denom(C_nD_idx(1,p))*C_nD_flux(p)
END DO

! Diagonal elements were only added once instead of twice -> double diagonal elements
DO p=0,PP_nCoefM1
  Am(p,p)=Am(p,p)*2.
END DO

SG_Fraction_Exact_nDim=enum
CALL DSYSV('U',PP_nCoef,1,DBLE(Am),PP_nCoef,Pivot_dummy,SG_Fraction_Exact_nDim,PP_nCoef,work,100,info  )

END FUNCTION SG_Fraction_Exact_nDim


!===================================================================================================================================
!> Fraction enum/denom projected onto stochstic domain; numerical evaulation
!===================================================================================================================================
FUNCTION SG_Fraction_Num(enum,denom)
! MODULES
USE MOD_PreProc
USE MOD_SG_Vars , ONLY: SG_VdM_OrthQuad,SG_VdM_QuadOrth
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
REAL,INTENT(IN),DIMENSION(0:PP_nCoefM1) :: enum,denom
REAL,DIMENSION(0:PP_nCoefM1)            :: SG_Fraction_Num
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!==================================================================================================================================
SG_Fraction_Num=MATMUL(SG_VdM_QuadOrth, MATMUL(SG_VdM_OrthQuad,enum)/MATMUL(SG_VdM_OrthQuad,denom))
END FUNCTION SG_Fraction_Num



!===================================================================================================================================
!> Dot product of two stoch. vars
!===================================================================================================================================
PURE FUNCTION SG_DOT_PRODUCT(v1,v2)
! MODULES
USE MOD_PreProc
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
REAL,INTENT(IN),DIMENSION(PP_nCoef*PP_dim) :: v1,v2
REAL,DIMENSION(0:PP_nCoefM1)               :: SG_DOT_PRODUCT
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!==================================================================================================================================
#if (PP_dim == 3)
SG_DOT_PRODUCT=SG_Product(v1(SGI1),v2(SGI1))+SG_Product(v1(SGI2),v2(SGI2))+SG_Product(v1(SGI3),v2(SGI3))
#else
SG_DOT_PRODUCT=SG_Product(v1(SGI1),v2(SGI1))+SG_Product(v1(SGI2),v2(SGI2))
#endif
END FUNCTION SG_DOT_PRODUCT



!===================================================================================================================================
!> Dot product of stoch. var and deterministic var
!===================================================================================================================================
PURE FUNCTION SG_DOT_PRODUCT_N(vS,vN)
! MODULES
USE MOD_PreProc
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
REAL,INTENT(IN),DIMENSION(PP_nCoef*PP_dim) :: vS
REAL,INTENT(IN),DIMENSION(         PP_dim) :: vN
REAL,DIMENSION(0:PP_nCoefM1)               :: SG_DOT_PRODUCT_N
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!==================================================================================================================================
#if (PP_dim == 3)
SG_DOT_PRODUCT_N=vS(SGI1)*vN(1)+vS(SGI2)*vN(2)+vS(SGI3)*vN(3)
#else
SG_DOT_PRODUCT_N=vS(SGI1)*vN(1)+vS(SGI2)*vN(2)
#endif
END FUNCTION SG_DOT_PRODUCT_N


!===================================================================================================================================
!> Product of two stoch vars
!===================================================================================================================================
PURE FUNCTION SG_Product_1Dim(v1,v2)
! MODULES
USE MOD_Preproc
USE MOD_SG_Vars , ONLY: C_1D
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
REAL,INTENT(IN),DIMENSION(0:PP_nCoefM1) :: v1,v2
REAL,DIMENSION(0:PP_nCoefM1)            :: SG_Product_1Dim
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                                    :: p,q,r,s
REAL,DIMENSION(0:PP_nCoefM1,0:PP_nCoefM1)  :: mat
!==================================================================================================================================
DO p=0,PP_nCoefM1
  DO q=0,p-1
    mat(q,p)=v1(p)*v2(q)+v1(q)*v2(p)
  END DO
  mat(p,p)=v1(p)*v2(p)
END DO
SG_Product_1Dim=0.
s=1
DO p=0,PP_nCoefM1
  DO q=0,p-1
    DO r=p-q,q-1,2
      SG_Product_1Dim(r)=SG_Product_1Dim(r)+C_1D(s)*mat(q,p)
      SG_Product_1Dim(q)=SG_Product_1Dim(q)+C_1D(s)*mat(r,p)
      SG_Product_1Dim(p)=SG_Product_1Dim(p)+C_1D(s)*mat(r,q)
      s=s+1
    END DO
  END DO
  DO q=0,p-1,2
    SG_Product_1Dim(p)=SG_Product_1Dim(p)+C_1D(s)*mat(q,p)
    SG_Product_1Dim(q)=SG_Product_1Dim(q)+C_1D(s)*mat(p,p)
    s=s+1
  END DO
END DO
DO p=0,PP_nCoefM1,2
  DO q=0,p-1
    SG_Product_1Dim(p)=SG_Product_1Dim(p)+C_1D(s)*mat(q,q)
    SG_Product_1Dim(q)=SG_Product_1Dim(q)+C_1D(s)*mat(q,p)
    s=s+1
  END DO
  SG_Product_1Dim(p)=SG_Product_1Dim(p)+C_1D(s)*mat(p,p)
  s=s+1
END DO
END FUNCTION SG_Product_1Dim

!===================================================================================================================================
!> Project exact solution u(xIn,t_In,xi) from ExactFunc onto polynomial (Hermite or Legendre) of degree pStoch
!> Two different orders for numerical integration are available (for numerics and for apprximative exact solution)
!===================================================================================================================================
PURE FUNCTION SG_Product_nDim(v1,v2)
! MODULES
USE MOD_Preproc
USE MOD_SG_Vars , ONLY: C_nD_flux,C_nD_idx,SizeC_nD
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
REAL,INTENT(IN),DIMENSION(0:PP_nCoefM1)          :: v1,v2
REAL           ,DIMENSION(0:PP_nCoefM1)          :: SG_Product_nDim
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                                         :: s,i
REAL,DIMENSION(0:PP_nCoefM1,0:PP_nCoefM1)         :: A
!==================================================================================================================================
SG_Product_nDim=0.

! First calculate tensor product of input vectors
DO i=0,PP_nCoefM1
  A(0:i,i) = v1(0:i)*v2(i)+v1(i)*v2(0:i)
END DO

DO s=1,SizeC_nD
  SG_Product_nDim(C_nD_idx(1,s)) = SG_Product_nDim(C_nD_idx(1,s)) + C_nD_flux(s)*A(C_nD_idx(3,s),C_nD_idx(2,s))
END DO
DO s=1,SizeC_nD
  SG_Product_nDim(C_nD_idx(2,s)) = SG_Product_nDim(C_nD_idx(2,s)) + C_nD_flux(s)*A(C_nD_idx(3,s),C_nD_idx(1,s))
END DO
DO s=1,SizeC_nD
  SG_Product_nDim(C_nD_idx(3,s)) = SG_Product_nDim(C_nD_idx(3,s)) + C_nD_flux(s)*A(C_nD_idx(2,s),C_nD_idx(1,s))
END DO

END FUNCTION SG_Product_nDim


END MODULE MOD_SG_Operators
