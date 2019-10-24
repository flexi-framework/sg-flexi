!=================================================================================================================================
! Copyright (c) 2010-2017 Prof. Claus-Dieter Munz
! Copyright (c) 2016-2017 Gregor Gassner (github.com/project-fluxo/fluxo)
! Copyright (c) 2016-2017 Florian Hindenlang (github.com/project-fluxo/fluxo)
! Copyright (c) 2016-2017 Andrew Winters (github.com/project-fluxo/fluxo)
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
#include "eos.h"
!==================================================================================================================================
!> Contains routines to compute the riemann (Advection, Diffusion) for a given Face
!==================================================================================================================================
MODULE MOD_Riemann
! MODULES
IMPLICIT NONE
PRIVATE
!----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
ABSTRACT INTERFACE
  SUBROUTINE RiemannInt(F_L,F_R,U_LL,U_RR,F)
    USE MOD_PreProc
    REAL,DIMENSION(PP_2Var),INTENT(IN) :: U_LL,U_RR
    REAL,DIMENSION(PP_nVar),INTENT(IN) :: F_L,F_R
    REAL,DIMENSION(PP_nVar),INTENT(OUT):: F
  END SUBROUTINE
END INTERFACE

PROCEDURE(RiemannInt),POINTER :: Riemann_pointer    !< pointer defining the standard inner Riemann solver
PROCEDURE(RiemannInt),POINTER :: RiemannBC_pointer  !< pointer defining the standard BC    Riemann solver

INTEGER,PARAMETER      :: PRM_RIEMANN_SAME          = -1
INTEGER,PARAMETER      :: PRM_RIEMANN_LF            = 1
INTEGER,PARAMETER      :: PRM_RIEMANN_HLLE          = 5
#ifdef SGdummy
INTEGER,PARAMETER      :: PRM_RIEMANN_HLLC          = 2
INTEGER,PARAMETER      :: PRM_RIEMANN_ROE           = 3
INTEGER,PARAMETER      :: PRM_RIEMANN_ROEL2         = 32
INTEGER,PARAMETER      :: PRM_RIEMANN_ROEENTROPYFIX = 33
INTEGER,PARAMETER      :: PRM_RIEMANN_HLL           = 4
INTEGER,PARAMETER      :: PRM_RIEMANN_HLLEM         = 6
#endif
#ifdef SPLIT_DG
INTEGER,PARAMETER      :: PRM_RIEMANN_Average       = 0
#endif

INTERFACE InitRiemann
  MODULE PROCEDURE InitRiemann
END INTERFACE

!INTERFACE Riemann
  !MODULE PROCEDURE Riemann
!END INTERFACE

#if PARABOLIC
INTERFACE ViscousFlux
  MODULE PROCEDURE ViscousFlux
END INTERFACE
PUBLIC::ViscousFlux
#endif

INTERFACE FinalizeRiemann
  MODULE PROCEDURE FinalizeRiemann
END INTERFACE


PUBLIC::InitRiemann
PUBLIC::Riemann
PUBLIC::FinalizeRiemann
!==================================================================================================================================

PUBLIC::DefineParametersRiemann
CONTAINS


!==================================================================================================================================
!> Define parameters
!==================================================================================================================================
SUBROUTINE DefineParametersRiemann()
! MODULES
USE MOD_Globals
USE MOD_ReadInTools ,ONLY: prms,addStrListEntry
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!==================================================================================================================================
CALL prms%SetSection("Riemann")
CALL prms%CreateIntFromStringOption('Riemann',   "Riemann solver to be used: LF, HLLC, Roe, RoeEntropyFix, HLL, HLLE, HLLEM", &
                                                 "RoeEntropyFix")
CALL addStrListEntry('Riemann','lf',           PRM_RIEMANN_LF)
CALL addStrListEntry('Riemann','hlle',         PRM_RIEMANN_HLLE)
#ifdef SGdummy
CALL addStrListEntry('Riemann','hllc',         PRM_RIEMANN_HLLC)
CALL addStrListEntry('Riemann','roe',          PRM_RIEMANN_ROE)
CALL addStrListEntry('Riemann','roeentropyfix',PRM_RIEMANN_ROEENTROPYFIX)
CALL addStrListEntry('Riemann','roel2',        PRM_RIEMANN_ROEL2)
CALL addStrListEntry('Riemann','hll',          PRM_RIEMANN_HLL)
CALL addStrListEntry('Riemann','hllem',        PRM_RIEMANN_HLLEM)
#endif
#ifdef SPLIT_DG
CALL addStrListEntry('Riemann','avg',          PRM_RIEMANN_Average)
#endif
CALL prms%CreateIntFromStringOption('RiemannBC', "Riemann solver used for boundary conditions: Same, LF, Roe, RoeEntropyFix, "//&
                                                 "HLL, HLLE, HLLEM",&
                                                 "Same")
CALL addStrListEntry('RiemannBC','lf',           PRM_RIEMANN_LF)
CALL addStrListEntry('RiemannBC','hlle',         PRM_RIEMANN_HLLE)
#ifdef SGdummy
CALL addStrListEntry('RiemannBC','hllc',         PRM_RIEMANN_HLLC)
CALL addStrListEntry('RiemannBC','roe',          PRM_RIEMANN_ROE)
CALL addStrListEntry('RiemannBC','roeentropyfix',PRM_RIEMANN_ROEENTROPYFIX)
CALL addStrListEntry('RiemannBC','roel2',        PRM_RIEMANN_ROEL2)
CALL addStrListEntry('RiemannBC','hll',          PRM_RIEMANN_HLL)
CALL addStrListEntry('RiemannBC','hllem',        PRM_RIEMANN_HLLEM)
#endif
#ifdef SPLIT_DG
CALL addStrListEntry('RiemannBC','avg',          PRM_RIEMANN_Average)
#endif
CALL addStrListEntry('RiemannBC','same',         PRM_RIEMANN_SAME)
END SUBROUTINE DefineParametersRiemann

!==================================================================================================================================!
!> Initialize Riemann solver routines, read inner and BC Riemann solver parameters and set pointers
!==================================================================================================================================!
SUBROUTINE InitRiemann()
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_ReadInTools ,ONLY: GETINTFROMSTR
!----------------------------------------------------------------------------------------------------------------------------------
IMPLICIT NONE
! INPUT / OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                 :: Riemann
#ifdef DEBUG
REAL,DIMENSION(PP_nVar) :: F_L,F_R,F ! dummy variables, only to suppress compiler warnings
REAL,DIMENSION(PP_2Var) :: U_LL,U_RR ! dummy variables, only to suppress compiler warnings
#endif
!==================================================================================================================================
#ifndef SPLIT_DG
Riemann = GETINTFROMSTR('Riemann')
SELECT CASE(Riemann)
CASE(PRM_RIEMANN_LF)
  Riemann_pointer => Riemann_LF
CASE(PRM_RIEMANN_HLLE)
  Riemann_pointer => Riemann_HLLE
#ifdef SGdummy
CASE(PRM_RIEMANN_HLLC)
  Riemann_pointer => Riemann_HLLC
CASE(PRM_RIEMANN_ROE)
  Riemann_pointer => Riemann_Roe
CASE(PRM_RIEMANN_ROEENTROPYFIX)
  Riemann_pointer => Riemann_RoeEntropyFix
CASE(PRM_RIEMANN_ROEL2)
  Riemann_pointer => Riemann_RoeL2
CASE(PRM_RIEMANN_HLL)
  Riemann_pointer => Riemann_HLL
CASE(PRM_RIEMANN_HLLEM)
  Riemann_pointer => Riemann_HLLEM
#endif
CASE DEFAULT
  CALL CollectiveStop(__STAMP__,&
    'Riemann solver not defined!')
END SELECT

Riemann = GETINTFROMSTR('RiemannBC')
SELECT CASE(Riemann)
CASE(PRM_RIEMANN_SAME)
  RiemannBC_pointer => Riemann_pointer
CASE(PRM_RIEMANN_LF)
  RiemannBC_pointer => Riemann_LF
CASE(PRM_RIEMANN_HLLE)
  RiemannBC_pointer => Riemann_HLLE
#ifdef SGdummy
CASE(PRM_RIEMANN_HLLC)
  RiemannBC_pointer => Riemann_HLLC
CASE(PRM_RIEMANN_ROE)
  RiemannBC_pointer => Riemann_Roe
CASE(PRM_RIEMANN_ROEENTROPYFIX)
  RiemannBC_pointer => Riemann_RoeEntropyFix
CASE(PRM_RIEMANN_ROEL2)
  RiemannBC_pointer => Riemann_RoeL2
CASE(PRM_RIEMANN_HLL)
  RiemannBC_pointer => Riemann_HLL
CASE(PRM_RIEMANN_HLLEM)
  RiemannBC_pointer => Riemann_HLLEM
#endif
CASE DEFAULT
  CALL CollectiveStop(__STAMP__,&
    'RiemannBC solver not defined!')
END SELECT

#else
Riemann = GETINTFROMSTR('Riemann')
SELECT CASE(Riemann)
CASE(PRM_RIEMANN_LF)
  Riemann_pointer => Riemann_LF
CASE(PRM_RIEMANN_HLLE)
  Riemann_pointer => Riemann_HLLE
#ifdef SGdummy
CASE(PRM_RIEMANN_ROE)
  Riemann_pointer => Riemann_Roe
CASE(PRM_RIEMANN_ROEENTROPYFIX)
  Riemann_pointer => Riemann_RoeEntropyFix
CASE(PRM_RIEMANN_ROEL2)
  Riemann_pointer => Riemann_RoeL2
CASE(PRM_RIEMANN_Average)
  Riemann_pointer => Riemann_FluxAverage
#endif
CASE DEFAULT
  CALL CollectiveStop(__STAMP__,&
    'Riemann solver not defined!')
END SELECT

Riemann = GETINTFROMSTR('RiemannBC')
SELECT CASE(Riemann)
CASE(PRM_RIEMANN_SAME)
  RiemannBC_pointer => Riemann_pointer
CASE(PRM_RIEMANN_LF)
  RiemannBC_pointer => Riemann_LF
CASE(PRM_RIEMANN_HLLE)
  RiemannBC_pointer => Riemann_HLLE
#ifdef SGdummy
CASE(PRM_RIEMANN_ROE)
  RiemannBC_pointer => Riemann_Roe
CASE(PRM_RIEMANN_ROEENTROPYFIX)
  RiemannBC_pointer => Riemann_RoeEntropyFix
CASE(PRM_RIEMANN_ROEL2)
  RiemannBC_pointer => Riemann_RoeL2
CASE(PRM_RIEMANN_Average)
  RiemannBC_pointer => Riemann_FluxAverage
#endif
CASE DEFAULT
  CALL CollectiveStop(__STAMP__,&
    'RiemannBC solver not defined!')
END SELECT
#endif /*SPLIT_DG*/

#ifdef DEBUG
! ===============================================================================
! Following dummy calls do suppress compiler warnings of unused Riemann-functions
! ===============================================================================
IF (0.EQ.1) THEN
  F_L=1. ;  F_R=1. ;   U_LL=1. ;   U_RR=1.
  CALL Riemann_pointer   (F_L,F_R,U_LL,U_RR,F)
  CALL RiemannBC_pointer (F_L,F_R,U_LL,U_RR,F)
  CALL Riemann_LF   (F_L,F_R,U_LL,U_RR,F)
  CALL Riemann_HLLC (F_L,F_R,U_LL,U_RR,F)
  CALL Riemann_Roe  (F_L,F_R,U_LL,U_RR,F)
  CALL Riemann_RoeEntropyFix(F_L,F_R,U_LL,U_RR,F)
  CALL Riemann_RoeL2(F_L,F_R,U_LL,U_RR,F)
  CALL Riemann_HLL  (F_L,F_R,U_LL,U_RR,F)
  CALL Riemann_HLLE (F_L,F_R,U_LL,U_RR,F)
  CALL Riemann_HLLEM(F_L,F_R,U_LL,U_RR,F)
#ifdef SPLIT_DG
  CALL Riemann_FluxAverage(F_L,F_R,U_LL,U_RR,F)
#endif
END IF
#endif /*DEBUG*/
END SUBROUTINE InitRiemann

!==================================================================================================================================
!> Computes the numerical flux
!> Conservative States are rotated into normal direction in this routine and are NOT backrotated: don't use it after this routine!!
!> Attention 2: numerical flux is backrotated at the end of the routine!!
!==================================================================================================================================
SUBROUTINE Riemann(Nloc,FOut,U_L,U_R,UPrim_L,UPrim_R,nv,t1,t2,doBC)
! MODULES
USE MOD_PreProc
USE MOD_Flux         ,ONLY:EvalEulerFlux1D_fast
USE MOD_SG_Operators, ONLY:SG_Inv,SG_Product,SG_DOT_PRODUCT_N
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
INTEGER,INTENT(IN)                                          :: Nloc       !< local polynomial degree
REAL,DIMENSION(PP_nVar    ,0:Nloc,0:ZDIM(Nloc)),INTENT(IN)  :: U_L        !< conservative solution at left side of the interface
REAL,DIMENSION(PP_nVar    ,0:Nloc,0:ZDIM(Nloc)),INTENT(IN)  :: U_R        !< conservative solution at right side of the interface
REAL,DIMENSION(PP_nVarPrim,0:Nloc,0:ZDIM(Nloc)),INTENT(IN)  :: UPrim_L    !< primitive solution at left side of the interface
REAL,DIMENSION(PP_nVarPrim,0:Nloc,0:ZDIM(Nloc)),INTENT(IN)  :: UPrim_R    !< primitive solution at right side of the interface
REAL,DIMENSION(                      3,0:Nloc,0:ZDIM(Nloc)),INTENT(IN)  :: nv,t1,t2   !< normal vector and tangential vectors at side
LOGICAL,INTENT(IN)                                          :: doBC       !< marker whether side is a BC side
REAL,DIMENSION(PP_nVar    ,0:Nloc,0:ZDIM(Nloc)),INTENT(OUT) :: FOut       !< advective flux
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                 :: i,j,iDim
REAL,DIMENSION(PP_nVar) :: F_L,F_R,F
REAL,DIMENSION(PP_2Var) :: U_LL,U_RR
PROCEDURE(RiemannInt),POINTER :: Riemann_loc !< pointer defining the standard inner Riemann solver
!==================================================================================================================================
IF (doBC) THEN
  Riemann_loc => RiemannBC_pointer
ELSE
  Riemann_loc => Riemann_pointer
END IF

! Momentum has to be rotatet using the normal system individual for each
DO j=0,ZDIM(Nloc); DO i=0,Nloc
  ! left state: U_L

  U_LL(TEMP)=0.
  U_RR(TEMP)=0.
  U_LL(DENS)=U_L(DENS,i,j)
  U_LL(SRHO)=SG_Inv(U_LL(DENS))
  U_LL(ENER)=U_L(ENER,i,j)
  U_LL(PRES)=UPrim_L(PPRES,i,j)


  ! rotate velocity in normal and tangential direction
  U_LL(VEL1)=SG_DOT_PRODUCT_N(UPrim_L(PVELV,i,j),nv(:,i,j))
  U_LL(VEL2)=SG_DOT_PRODUCT_N(UPrim_L(PVELV,i,j),t1(:,i,j))
  U_LL(MOM1)=SG_Product(U_LL(DENS),U_LL(VEL1))
  U_LL(MOM2)=SG_Product(U_LL(DENS),U_LL(VEL2))
#if PP_dim==3
  U_LL(VEL3)=SG_DOT_PRODUCT_N(UPrim_L(PVELV,i,j),t2(:,i,j))
  U_LL(MOM3)=SG_Product(U_LL(DENS),U_LL(VEL3))
#else
  U_LL(VEL3)=0.
  U_LL(MOM3)=0.
#endif
  ! right state: U_R
  U_RR(DENS)=U_R(DENS,i,j)
  U_RR(SRHO)=SG_Inv(U_RR(DENS))
  U_RR(ENER)=U_R(ENER,i,j)
  U_RR(PRES)=UPrim_R(PPRES,i,j)
  ! rotate momentum in normal and tangential direction
  U_RR(VEL1)=SG_DOT_PRODUCT_N(UPRIM_R(PVELV,i,j),nv(:,i,j))
  U_RR(VEL2)=SG_DOT_PRODUCT_N(UPRIM_R(PVELV,i,j),t1(:,i,j))
  U_RR(MOM1)=SG_Product(U_RR(DENS),U_RR(VEL1))
  U_RR(MOM2)=SG_Product(U_RR(DENS),U_RR(VEL2))
#if PP_dim==3
  U_RR(VEL3)=SG_DOT_PRODUCT_N(UPRIM_R(PVELV,i,j),t2(:,i,j))
  U_RR(MOM3)=SG_Product(U_RR(DENS),U_RR(VEL3))
#else
  U_RR(VEL3)=0.
  U_RR(MOM3)=0.
#endif

# ifndef SPLIT_DG
  CALL EvalEulerFlux1D_fast(U_LL,F_L)
  CALL EvalEulerFlux1D_fast(U_RR,F_R)
#endif /*SPLIT_DG*/

  CALL Riemann_loc(F_L,F_R,U_LL,U_RR,F)
  ! Back Rotate the normal flux into Cartesian direction
  Fout(DENS,i,j) = F(DENS)
  DO iDim=1,PP_dim
    Fout( PP_nCoef*iDim+1:PP_nCoef*(1+iDim) ,i,j) = nv(iDim,i,j)*F(MOM1) + t1(iDim,i,j)*F(MOM2) &
#if PP_dim==3
        + t2(iDim,i,j)*F(MOM3)
#else
        + 0.
#endif
  END DO
  Fout(ENER,i,j)=F(ENER)
END DO; END DO

#ifdef DEBUG
! ===============================================================================
! Following dummy calls do suppress compiler warnings of unused Riemann-functions
! ===============================================================================
IF (0.EQ.1) THEN
  WRITE (*,*) t2
END IF
#endif
END SUBROUTINE Riemann



#if PARABOLIC
!==================================================================================================================================
!> Computes the viscous NSE diffusion fluxes in all directions to approximate the numerical flux
!> Actually not a Riemann solver, only here for coding reasons
!==================================================================================================================================
SUBROUTINE ViscousFlux(Nloc,F,UPrim_L,UPrim_R, &
                       gradUx_L,gradUy_L,gradUz_L,gradUx_R,gradUy_R,gradUz_R,nv&
#ifdef EDDYVISCOSITY
                      ,muSGS_L,muSGS_R&
#endif
                      )
! MODULES
USE MOD_Preproc
USE MOD_Flux,ONLY: EvalDiffFlux2D
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
INTEGER,INTENT(IN)                                         :: Nloc     !< local polynomial degree
                                                           !> solution in primitive variables at left/right side of the interface
REAL,DIMENSION(PP_nVarPrim,0:Nloc,0:ZDIM(Nloc)),INTENT(IN)   :: UPrim_L,UPrim_R
                                                           !> solution gradients in x/y/z-direction left/right of the interface
REAL,DIMENSION(PP_nVarPrim,0:Nloc,0:ZDIM(Nloc)),INTENT(IN)   :: gradUx_L,gradUx_R,gradUy_L,gradUy_R,gradUz_L,gradUz_R
REAL,INTENT(IN)                                            :: nv(3,0:Nloc,0:ZDIM(Nloc)) !< normal vector
REAL,INTENT(OUT)                                           :: F(PP_nVar,0:Nloc,0:ZDIM(Nloc)) !< viscous flux
#ifdef EDDYVISCOSITY
                                                           !> eddy viscosity left/right of the interface
REAL,DIMENSION(1,0:Nloc,0:ZDIM(Nloc)),INTENT(IN)             :: muSGS_L,muSGS_R
#endif
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                                              :: p,q
REAL,DIMENSION(PP_nVar,0:Nloc,0:ZDIM(Nloc))            :: diffFluxX_L,diffFluxY_L,diffFluxZ_L
REAL,DIMENSION(PP_nVar,0:Nloc,0:ZDIM(Nloc))            :: diffFluxX_R,diffFluxY_R,diffFluxZ_R
!==================================================================================================================================
! Don't forget the diffusion contribution, my young padawan
! Compute NSE Diffusion flux
  CALL EvalDiffFlux2D(Nloc,diffFluxX_L,diffFluxY_L,diffFluxZ_L,UPrim_L,gradUx_L,gradUy_L,gradUz_L &
#ifdef EDDYVISCOSITY
                    ,muSGS_L&
#endif
      )
  CALL EvalDiffFlux2D(Nloc,diffFluxX_R,diffFluxY_R,diffFluxZ_R,UPrim_R,gradUx_R,gradUy_R,gradUz_R &
#ifdef EDDYVISCOSITY
                    ,muSGS_R&
#endif
      )
! BR1 uses arithmetic mean of the fluxes
DO q=0,ZDIM(Nloc); DO p=0,Nloc
  F(:,p,q)=0.5*(nv(1,p,q)*(diffFluxX_L(1:PP_nVar,p,q)+diffFluxX_R(1:PP_nVar,p,q)) &
               +nv(2,p,q)*(diffFluxY_L(1:PP_nVar,p,q)+diffFluxY_R(1:PP_nVar,p,q)) &
               +nv(3,p,q)*(diffFluxZ_L(1:PP_nVar,p,q)+diffFluxZ_R(1:PP_nVar,p,q)))
END DO; END DO
END SUBROUTINE ViscousFlux
#endif /* PARABOLIC */





!==================================================================================================================================
!> Local Lax-Friedrichs (Rusanov) Riemann solver
!==================================================================================================================================
SUBROUTINE Riemann_LF(F_L,F_R,U_LL,U_RR,F)
! MODULES
USE MOD_PreProc
USE MOD_EOS_Vars      ,ONLY: Kappa
USE MOD_SG_Vars       ,ONLY: SG_Vdm_OrthQuad,SG_SafetyFactor,nQPTotal
#ifdef SPLIT_DG
USE MOD_SplitFlux     ,ONLY: SplitDGSurface_pointer
#endif /*SPLIT_DG*/
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
                                                !> extended solution vector on the left/right side of the interface
REAL,DIMENSION(PP_2Var),INTENT(IN) :: U_LL,U_RR
                                                !> advection fluxes on the left/right side of the interface
REAL,DIMENSION(PP_nVar),INTENT(IN) :: F_L,F_R
REAL,DIMENSION(PP_nVar),INTENT(OUT):: F         !< resulting Riemann flux
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER             :: iVarDet
REAL                :: LambdaMax
REAL                :: c_LL_Quad(1:nQPTotal),c_RR_Quad(1:nQPTotal),LambdaMax_Quad(1:nQPTotal)
REAL                :: U_LL_Quad(PP_2VarQ),U_RR_Quad(PP_2VarQ)
!==================================================================================================================================
! Lax-Friedrichs
DO iVarDet=1,PP_2VarDet-1
  U_RR_Quad(PP_iVarDetQ)=MATMUL(SG_Vdm_OrthQuad,U_RR(PP_iVarDet))
  U_LL_Quad(PP_iVarDetQ)=MATMUL(SG_Vdm_OrthQuad,U_LL(PP_iVarDet))
END DO


c_LL_Quad=SQRT(kappa*U_LL_Quad(QPRES)/U_LL_Quad(QDENS))
c_RR_Quad=SQRT(kappa*U_RR_Quad(QPRES)/U_RR_Quad(QDENS))
LambdaMax_Quad = MAX( ABS(U_RR_Quad(QVEL1)),ABS(U_LL_Quad(QVEL1)) ) + MAX( c_LL_Quad,c_RR_Quad )
LambdaMax = SG_SafetyFactor * MAXVAL(LambdaMax_Quad)
#ifndef SPLIT_DG
F = 0.5*((F_L+F_R) - LambdaMax*(U_RR(CONS) - U_LL(CONS)))
#else
! get split flux
CALL SplitDGSurface_pointer(U_LL,U_RR,F)
! compute surface flux
F = F - 0.5*LambdaMax*(U_RR(CONS) - U_LL(CONS))
#endif /*SPLIT_DG*/

#ifdef DEBUG
! ===============================================================================
! Following dummy calls do suppress compiler warnings of unused Riemann-functions
! ===============================================================================
IF (0.EQ.1) THEN
  WRITE (*,*) F_L,F_R
END IF
#endif /*DEBUG*/
END SUBROUTINE Riemann_LF


#ifdef SGdummy
!=================================================================================================================================
!> Harten-Lax-Van-Leer Riemann solver resolving contact discontinuity
!=================================================================================================================================
PURE SUBROUTINE Riemann_HLLC(F_L,F_R,U_LL,U_RR,F)
! MODULES
USE MOD_PreProc       ,ONLY: PP_nVar
USE MOD_EOS_Vars      ,ONLY: KappaM1
IMPLICIT NONE
!---------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
                                           !> extended solution vector on the left/right side of the interface
REAL,DIMENSION(PP_2Var),INTENT(IN) :: U_LL,U_RR
                                           !> advection fluxes on the left/right side of the interface
REAL,DIMENSION(PP_nVar),INTENT(IN) :: F_L,F_R
REAL,DIMENSION(PP_nVar),INTENT(OUT):: F    !< resulting Riemann flux
!---------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL    :: H_L,H_R
REAL    :: SqrtRho_L,SqrtRho_R,sSqrtRho
REAL    :: RoeVel(3),RoeH,Roec,absVel
REAL    :: Ssl,Ssr,SStar
REAL    :: U_Star(PP_nVar),EStar
REAL    :: sMu_L,sMu_R
!REAL    :: c_L,c_R
!=================================================================================================================================
! HLLC flux

! Version A: Basic Davis estimate for wave speed
!Ssl = U_LL(VEL1) - SPEEDOFSOUND_HE(U_LL)
!Ssr = U_RR(VEL1) + SPEEDOFSOUND_HE(U_RR)

! Version B: Basic Davis estimate for wave speed
!c_L = SPEEDOFSOUND_HE(U_LL)
!c_R = SPEEDOFSOUND_HE(U_RR)
!Ssl = MIN(U_LL(VEL1) - c_L,U_RR(VEL1) - c_R)
!Ssr = MAX(U_LL(VEL1) + c_L,U_RR(VEL1) + c_R)

! Version C: Better Roe estimate for wave speeds Davis, Einfeldt
H_L       = TOTALENTHALPY_HE(U_LL)
H_R       = TOTALENTHALPY_HE(U_RR)
SqrtRho_L = SQRT(U_LL(DENS))
SqrtRho_R = SQRT(U_RR(DENS))
sSqrtRho  = 1./(SqrtRho_L+SqrtRho_R)
! Roe mean values
RoeVel    = (SqrtRho_R*U_RR(VELV) + SqrtRho_L*U_LL(VELV)) * sSqrtRho
RoeH      = (SqrtRho_R*H_R        + SqrtRho_L*H_L       ) * sSqrtRho
absVel    = DOT_PRODUCT(RoeVel,RoeVel)
Roec      = SQRT(KappaM1*(RoeH-0.5*absVel))
Ssl = RoeVel(1) - Roec
Ssr = RoeVel(1) + Roec

! positive supersonic speed
IF(Ssl .GE. 0.)THEN
  F=F_L
! negative supersonic speed
ELSEIF(Ssr .LE. 0.)THEN
  F=F_R
! subsonic case
ELSE
  sMu_L = Ssl - U_LL(VEL1)
  sMu_R = Ssr - U_RR(VEL1)
  SStar = (U_RR(PRES) - U_LL(PRES) + U_LL(MOM1)*sMu_L - U_RR(MOM1)*sMu_R) / (U_LL(DENS)*sMu_L - U_RR(DENS)*sMu_R)
  IF ((Ssl .LE. 0.).AND.(SStar .GE. 0.)) THEN
    EStar  = TOTALENERGY_HE(U_LL) + (SStar-U_LL(VEL1))*(SStar + U_LL(PRES)*U_LL(SRHO)/sMu_L)
    U_Star = U_LL(DENS) * sMu_L/(Ssl-SStar) * (/ 1., SStar, (/U_LL(VEL2),U_LL(VEL3)/), EStar /)
    F=F_L+Ssl*(U_Star-U_LL(CONS))
  ELSE
    EStar  = TOTALENERGY_HE(U_RR) + (SStar-U_RR(VEL1))*(SStar + U_RR(PRES)*U_RR(SRHO)/sMu_R)
    U_Star = U_RR(DENS) * sMu_R/(Ssr-SStar) * (/ 1., SStar, (/U_RR(VEL2),U_RR(VEL3)/), EStar /)
    F=F_R+Ssr*(U_Star-U_RR(CONS))
  END IF
END IF ! subsonic case
END SUBROUTINE Riemann_HLLC


!=================================================================================================================================
!> Roe's approximate Riemann solver
!=================================================================================================================================
SUBROUTINE Riemann_Roe(F_L,F_R,U_LL,U_RR,F)
! MODULES
USE MOD_PreProc   ,ONLY: PP_nVar
USE MOD_EOS_Vars  ,ONLY: kappaM1
#ifdef SPLIT_DG
USE MOD_SplitFlux ,ONLY: SplitDGSurface_pointer
#endif /*SPLIT_DG*/
IMPLICIT NONE
!---------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
                                               !> extended solution vector on the left/right side of the interface
REAL,DIMENSION(PP_2Var),INTENT(IN) :: U_LL,U_RR
                                               !> advection fluxes on the left/right side of the interface
REAL,DIMENSION(PP_nVar),INTENT(IN) :: F_L,F_R
REAL,DIMENSION(PP_nVar),INTENT(OUT):: F        !< resulting Riemann flux
!---------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                    :: H_L,H_R
REAL                    :: SqrtRho_L,SqrtRho_R,sSqrtRho
REAL                    :: RoeVel(3),RoeH,Roec,absVel
REAL,DIMENSION(PP_nVar) :: a,r1,r2,r3,r4,r5  ! Roe eigenvectors
REAL                    :: Alpha1,Alpha2,Alpha3,Alpha4,Alpha5,Delta_U(PP_nVar+1)
!=================================================================================================================================
! Roe flux
H_L       = TOTALENTHALPY_HE(U_LL)
H_R       = TOTALENTHALPY_HE(U_RR)
SqrtRho_L = SQRT(U_LL(DENS))
SqrtRho_R = SQRT(U_RR(DENS))

sSqrtRho  = 1./(SqrtRho_L+SqrtRho_R)
! Roe mean values
RoeVel    = (SqrtRho_R*U_RR(VELV) + SqrtRho_L*U_LL(VELV)) * sSqrtRho
absVel    = DOT_PRODUCT(RoeVel,RoeVel)
RoeH      = (SqrtRho_R*H_R+SqrtRho_L*H_L) * sSqrtRho
Roec      = ROEC_RIEMANN_H(RoeH,RoeVel)

! mean eigenvalues and eigenvectors
a  = (/ RoeVel(1)-Roec, RoeVel(1), RoeVel(1), RoeVel(1), RoeVel(1)+Roec      /)
r1 = (/ 1.,             a(1),      RoeVel(2), RoeVel(3), RoeH-RoeVel(1)*Roec /)
r2 = (/ 1.,             RoeVel(1), RoeVel(2), RoeVel(3), 0.5*absVel          /)
r3 = (/ 0.,             0.,        1.,        0.,        RoeVel(2)           /)
r4 = (/ 0.,             0.,        0.,        1.,        RoeVel(3)           /)
r5 = (/ 1.,             a(5),      RoeVel(2), RoeVel(3), RoeH+RoeVel(1)*Roec /)

! calculate differences
Delta_U(1:5) = U_RR(CONS) - U_LL(CONS)
Delta_U(6)   = Delta_U(5)-(Delta_U(3)-RoeVel(2)*Delta_U(1))*RoeVel(2) - (Delta_U(4)-RoeVel(3)*Delta_U(1))*RoeVel(3)
! calculate factors
Alpha3 = Delta_U(3) - RoeVel(2)*Delta_U(1)
Alpha4 = Delta_U(4) - RoeVel(3)*Delta_U(1)
Alpha2 = ALPHA2_RIEMANN_H(RoeH,RoeVel,Roec,Delta_U)
Alpha1 = 0.5/Roec * (Delta_U(1)*(RoeVel(1)+Roec) - Delta_U(2) - Roec*Alpha2)
Alpha5 = Delta_U(1) - Alpha1 - Alpha2
#ifndef SPLIT_DG
! assemble Roe flux
F=0.5*((F_L+F_R) - &
       Alpha1*ABS(a(1))*r1 - &
       Alpha2*ABS(a(2))*r2 - &
       Alpha3*ABS(a(3))*r3 - &
       Alpha4*ABS(a(4))*r4 - &
       Alpha5*ABS(a(5))*r5)
#else
! get split flux
CALL SplitDGSurface_pointer(U_LL,U_RR,F)
! assemble Roe flux
F = F - 0.5*(Alpha1*ABS(a(1))*r1 + &
             Alpha2*ABS(a(2))*r2 + &
             Alpha3*ABS(a(3))*r3 + &
             Alpha4*ABS(a(4))*r4 + &
             Alpha5*ABS(a(5))*r5)
#endif /*SPLIT_DG*/

#ifdef DEBUG
! ===============================================================================
! Following dummy calls do suppress compiler warnings of unused Riemann-functions
! ===============================================================================
IF (0.EQ.1) THEN
  WRITE (*,*) F_L,F_R
END IF
#endif /*DEBUG*/
END SUBROUTINE Riemann_Roe


!=================================================================================================================================
!> Roe's approximate Riemann solver using the Harten and Hymen II entropy fix, see
!> Pelanti, Marica & Quartapelle, Luigi & Vigevano, L & Vigevano, Luigi. (2018):
!>  A review of entropy fixes as applied to Roe's linearization. 
!=================================================================================================================================
SUBROUTINE Riemann_RoeEntropyFix(F_L,F_R,U_LL,U_RR,F)
! MODULES
USE MOD_PreProc       ,ONLY: PP_nVar
USE MOD_EOS_Vars      ,ONLY: Kappa,KappaM1
#ifdef SPLIT_DG
USE MOD_SplitFlux ,ONLY: SplitDGSurface_pointer
#endif /*SPLIT_DG*/
IMPLICIT NONE
!---------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
                                               !> extended solution vector on the left/right side of the interface
REAL,DIMENSION(PP_2Var),INTENT(IN) :: U_LL,U_RR
                                               !> advection fluxes on the left/right side of the interface
REAL,DIMENSION(PP_nVar),INTENT(IN) :: F_L,F_R
REAL,DIMENSION(PP_nVar),INTENT(OUT):: F        !< resulting Riemann flux
!---------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                 :: iVar
REAL                    :: c_L,c_R
REAL                    :: H_L,H_R
REAL                    :: SqrtRho_L,SqrtRho_R,sSqrtRho,absVel
REAL                    :: RoeVel(3),RoeH,Roec,RoeDens
REAL,DIMENSION(PP_nVar) :: r1,r2,r3,r4,r5,a,al,ar,Delta_U,Alpha  ! Roe eigenvectors
REAL                    :: tmp,da
!=================================================================================================================================
c_L       = SPEEDOFSOUND_HE(U_LL)
c_R       = SPEEDOFSOUND_HE(U_RR)
H_L       = TOTALENTHALPY_HE(U_LL)
H_R       = TOTALENTHALPY_HE(U_RR)
SqrtRho_L = SQRT(U_LL(DENS))
SqrtRho_R = SQRT(U_RR(DENS))

sSqrtRho  = 1./(SqrtRho_L+SqrtRho_R)
! Roe mean values
RoeVel    = (SqrtRho_R*U_RR(VELV) + SqrtRho_L*U_LL(VELV)) * sSqrtRho
RoeH      = (SqrtRho_R*H_R+SqrtRho_L*H_L) * sSqrtRho
absVel    = DOT_PRODUCT(RoeVel,RoeVel)
Roec      = ROEC_RIEMANN_H(RoeH,RoeVel)
RoeDens   = SQRT(U_LL(DENS)*U_RR(DENS))
! Roe+Pike version of Roe Riemann solver

! calculate jump
Delta_U(1)   = U_RR(DENS) - U_LL(DENS)
Delta_U(2:4) = U_RR(VELV) - U_LL(VELV)
Delta_U(5)   = U_RR(PRES) - U_LL(PRES)

! mean eigenvalues and eigenvectors
a  = (/ RoeVel(1)-Roec, RoeVel(1), RoeVel(1), RoeVel(1), RoeVel(1)+Roec      /)
r1 = (/ 1.,             a(1),      RoeVel(2), RoeVel(3), RoeH-RoeVel(1)*Roec /)
r2 = (/ 1.,             RoeVel(1), RoeVel(2), RoeVel(3), 0.5*absVel          /)
r3 = (/ 0.,             0.,        1.,        0.,        RoeVel(2)           /)
r4 = (/ 0.,             0.,        0.,        1.,        RoeVel(3)           /)
r5 = (/ 1.,             a(5),      RoeVel(2), RoeVel(3), RoeH+RoeVel(1)*Roec /)

! calculate wave strenghts
tmp      = 0.5/(Roec*Roec)
Alpha(1) = tmp*(Delta_U(5)-RoeDens*Roec*Delta_U(2))
Alpha(2) = Delta_U(1) - Delta_U(5)*2.*tmp
Alpha(3) = RoeDens*Delta_U(3)
Alpha(4) = RoeDens*Delta_U(4)
Alpha(5) = tmp*(Delta_U(5)+RoeDens*Roec*Delta_U(2))

! Harten+Hyman entropy fix (apply only for acoustic waves, don't fix r)

al(1) = U_LL(VEL1) - c_L
al(2) = U_LL(VEL1)
al(3) = U_LL(VEL1)
al(4) = U_LL(VEL1)
al(5) = U_LL(VEL1) + c_L
ar(1) = U_RR(VEL1) - c_R
ar(2) = U_RR(VEL1)
ar(3) = U_RR(VEL1)
ar(4) = U_RR(VEL1)
ar(5) = U_RR(VEL1) + c_R
! HH1
!IF(ABS(a(1)).LT.da1) a(1)=da1
!IF(ABS(a(5)).LT.da5) a(5)=da5
! HH2
DO iVar=1,5
  da = MAX(0.,a(iVar)-al(iVar),ar(iVar)-a(iVar))

  IF(ABS(a(iVar)).LT.da) THEN
    a(iVar)=0.5*(a(iVar)*a(iVar)/da+da)
  ELSE
    a(iVar) = ABS(a(iVar))
  END IF
END DO

#ifndef SPLIT_DG
! assemble Roe flux
F=0.5*((F_L+F_R)        - &
       Alpha(1)*a(1)*r1 - &
       Alpha(2)*a(2)*r2 - &
       Alpha(3)*a(3)*r3 - &
       Alpha(4)*a(4)*r4 - &
       Alpha(5)*a(5)*r5)
#else
! get split flux
CALL SplitDGSurface_pointer(U_LL,U_RR,F)
! for KG or PI flux eigenvalues have to be altered to ensure consistent KE dissipation
! assemble Roe flux
F= F - 0.5*(Alpha(1)*a(1)*r1 + &
            Alpha(2)*a(2)*r2 + &
            Alpha(3)*a(3)*r3 + &
            Alpha(4)*a(4)*r4 + &
            Alpha(5)*a(5)*r5)
#endif /*SPLIT_DG*/

#ifdef DEBUG
! ===============================================================================
! Following dummy calls do suppress compiler warnings of unused Riemann-functions
! ===============================================================================
IF (0.EQ.1) THEN
  WRITE (*,*) F_L,F_R
END IF
#endif /*DEBUG*/
END SUBROUTINE Riemann_RoeEntropyFix

!=================================================================================================================================
!> low mach number Roe's approximate Riemann solver according to Oßwald(2015)
!=================================================================================================================================
SUBROUTINE Riemann_RoeL2(F_L,F_R,U_LL,U_RR,F)
! MODULES
USE MOD_PreProc   ,ONLY: PP_nVar
USE MOD_EOS_Vars  ,ONLY: kappaM1,kappa
#ifdef SPLIT_DG
USE MOD_SplitFlux ,ONLY: SplitDGSurface_pointer
#endif /*SPLIT_DG*/
IMPLICIT NONE
!---------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
                                               !> extended solution vector on the left/right side of the interface
REAL,DIMENSION(PP_2Var),INTENT(IN) :: U_LL,U_RR
                                               !> advection fluxes on the left/right side of the interface
REAL,DIMENSION(PP_nVar),INTENT(IN) :: F_L,F_R
REAL,DIMENSION(PP_nVar),INTENT(OUT):: F        !< resulting Riemann flux
!---------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                    :: H_L,H_R
REAL                    :: SqrtRho_L,SqrtRho_R,sSqrtRho
REAL                    :: RoeVel(3),RoeH,Roec,absVel
REAL                    :: Ma_loc ! local Mach-Number
REAL,DIMENSION(PP_nVar) :: a,r1,r2,r3,r4,r5  ! Roe eigenvectors
REAL                    :: Alpha1,Alpha2,Alpha3,Alpha4,Alpha5,Delta_U(PP_nVar+1)
!=================================================================================================================================
! Roe flux
H_L       = TOTALENTHALPY_HE(U_LL)
H_R       = TOTALENTHALPY_HE(U_RR)
SqrtRho_L = SQRT(U_LL(DENS))
SqrtRho_R = SQRT(U_RR(DENS))

sSqrtRho  = 1./(SqrtRho_L+SqrtRho_R)
! Roe mean values
RoeVel    = (SqrtRho_R*U_RR(VELV) + SqrtRho_L*U_LL(VELV)) * sSqrtRho
absVel    = DOT_PRODUCT(RoeVel,RoeVel)
RoeH      = (SqrtRho_R*H_R+SqrtRho_L*H_L) * sSqrtRho
Roec      = ROEC_RIEMANN_H(RoeH,RoeVel)

! mean eigenvalues and eigenvectors
a  = (/ RoeVel(1)-Roec, RoeVel(1), RoeVel(1), RoeVel(1), RoeVel(1)+Roec      /)
r1 = (/ 1.,             a(1),      RoeVel(2), RoeVel(3), RoeH-RoeVel(1)*Roec /)
r2 = (/ 1.,             RoeVel(1), RoeVel(2), RoeVel(3), 0.5*absVel          /)
r3 = (/ 0.,             0.,        1.,        0.,        RoeVel(2)           /)
r4 = (/ 0.,             0.,        0.,        1.,        RoeVel(3)           /)
r5 = (/ 1.,             a(5),      RoeVel(2), RoeVel(3), RoeH+RoeVel(1)*Roec /)

! calculate differences
Delta_U(1:5) = U_RR(CONS) - U_LL(CONS)
Delta_U(6)   = Delta_U(5)-(Delta_U(3)-RoeVel(2)*Delta_U(1))*RoeVel(2) - (Delta_U(4)-RoeVel(3)*Delta_U(1))*RoeVel(3)

! low Mach-Number fix
Ma_loc = SQRT(absVel)/(Roec*SQRT(kappa))
Delta_U(2:4) = Delta_U(2:4) * Ma_loc

! calculate factors
Alpha3 = Delta_U(3) - RoeVel(2)*Delta_U(1)
Alpha4 = Delta_U(4) - RoeVel(3)*Delta_U(1)
Alpha2 = ALPHA2_RIEMANN_H(RoeH,RoeVel,Roec,Delta_U)
Alpha1 = 0.5/Roec * (Delta_U(1)*(RoeVel(1)+Roec) - Delta_U(2) - Roec*Alpha2)
Alpha5 = Delta_U(1) - Alpha1 - Alpha2

#ifndef SPLIT_DG
! assemble Roe flux
F=0.5*((F_L+F_R) - &
       Alpha1*ABS(a(1))*r1 - &
       Alpha2*ABS(a(2))*r2 - &
       Alpha3*ABS(a(3))*r3 - &
       Alpha4*ABS(a(4))*r4 - &
       Alpha5*ABS(a(5))*r5)
#else
! get split flux
CALL SplitDGSurface_pointer(U_LL,U_RR,F)
! assemble Roe flux
F = F - 0.5*(Alpha1*ABS(a(1))*r1 + &
             Alpha2*ABS(a(2))*r2 + &
             Alpha3*ABS(a(3))*r3 + &
             Alpha4*ABS(a(4))*r4 + &
             Alpha5*ABS(a(5))*r5)
#endif /*SPLIT_DG*/

#ifdef DEBUG
! ===============================================================================
! Following dummy calls do suppress compiler warnings of unused Riemann-functions
! ===============================================================================
IF (0.EQ.1) THEN
  WRITE (*,*) F_L,F_R
END IF
#endif /*DEBUG*/
END SUBROUTINE Riemann_RoeL2


!=================================================================================================================================
!> Standard Harten-Lax-Van-Leer Riemann solver without contact discontinuity
!=================================================================================================================================
PURE SUBROUTINE Riemann_HLL(F_L,F_R,U_LL,U_RR,F)
! MODULES
USE MOD_PreProc,  ONLY: PP_nVar
USE MOD_EOS_Vars, ONLY: KappaM1
IMPLICIT NONE
!---------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
                                               !> extended solution vector on the left/right side of the interface
REAL,DIMENSION(PP_2Var),INTENT(IN) :: U_LL,U_RR
                                               !> advection fluxes on the left/right side of the interface
REAL,DIMENSION(PP_nVar),INTENT(IN) :: F_L,F_R
REAL,DIMENSION(PP_nVar),INTENT(OUT):: F        !< resulting Riemann flux
!---------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL    :: H_L,H_R
REAL    :: SqrtRho_L,SqrtRho_R,sSqrtRho,absVel
REAL    :: RoeVel(3),RoeH,Roec
REAL    :: Ssl,Ssr
!=================================================================================================================================
H_L       = TOTALENTHALPY_HE(U_LL)
H_R       = TOTALENTHALPY_HE(U_RR)
SqrtRho_L = SQRT(U_LL(DENS))
SqrtRho_R = SQRT(U_RR(DENS))
sSqrtRho  = 1./(SqrtRho_L+SqrtRho_R)
! Roe mean values
RoeVel    = (SqrtRho_R*U_RR(VELV) + SqrtRho_L*U_LL(VELV)) * sSqrtRho
RoeH      = (SqrtRho_R*H_R        + SqrtRho_L*H_L)        * sSqrtRho
absVel    = DOT_PRODUCT(RoeVel,RoeVel)
Roec      = ROEC_RIEMANN_H(RoeH,RoeVel)
! HLL flux
! Basic Davis estimate for wave speed
!Ssl = U_LL(VEL1) - c_L
!Ssr = U_RR(VEL1) + c_R
! Better Roe estimate for wave speeds Davis, Einfeldt
Ssl = RoeVel(1) - Roec
Ssr = RoeVel(1) + Roec
! positive supersonic speed
IF(Ssl .GE. 0.)THEN
  F=F_L
! negative supersonic speed
ELSEIF(Ssr .LE. 0.)THEN
  F=F_R
! subsonic case
ELSE
  F=(Ssr*F_L-Ssl*F_R+Ssl*Ssr*(U_RR(CONS)-U_LL(CONS)))/(Ssr-Ssl)
END IF ! subsonic case
END SUBROUTINE Riemann_HLL

#endif /*SGdummy*/
!=================================================================================================================================
!> Harten-Lax-Van-Leer-Einfeldt Riemann solver
!=================================================================================================================================
SUBROUTINE Riemann_HLLE(F_L,F_R,U_LL,U_RR,F)
!=================================================================================================================================
! MODULES
USE MOD_PreProc
USE MOD_EOS_Vars      ,ONLY: Kappa,KappaM1
USE MOD_SG_Vars       ,ONLY: nQPTotal, SG_Vdm_OrthQuad, SG_Vdm_QuadOrth
IMPLICIT NONE
!---------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
                                               !> extended solution vector on the left/right side of the interface
REAL,DIMENSION(PP_2Var),INTENT(IN) :: U_LL,U_RR
                                               !> advection fluxes on the left/right side of the interface
REAL,DIMENSION(PP_nVar),INTENT(IN) :: F_L,F_R
REAL,DIMENSION(PP_nVar),INTENT(OUT):: F        !< resulting Riemann flux
!---------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL    :: U_LL_Quad(PP_2VarQ),U_RR_Quad(PP_2VarQ)
REAL    :: F_LL_Quad(PP_nVarQ),F_RR_Quad(PP_nVarQ)
REAL    :: c_LL_Quad(nQPTotal),c_RR_Quad(nQPTotal)
REAL    :: F_Quad(PP_nVarQ)
REAL    :: H_L(nQPTotal),H_R(nQPTotal)
REAL    :: SqrtRho_L(nQPTotal) ,SqrtRho_R(nQPTotal) ,sSqrtRho(nQPTotal) ,absVel(nQPTotal)
REAL    :: RoeVel(3*nQPTotal),RoeH(nQPTotal),Roec(nQPTotal)
REAL    :: SsL(nQPTotal),SsR(nQPTotal),beta
INTEGER :: iQP, iVarDet
!=================================================================================================================================
DO iVarDet=1,PP_2VarDet
  U_RR_Quad(PP_iVarDetQ) = MATMUL(SG_Vdm_OrthQuad,U_RR(PP_iVarDet))
  U_LL_Quad(PP_iVarDetQ) = MATMUL(SG_Vdm_OrthQuad,U_LL(PP_iVarDet))
END DO

DO iVarDet=1,PP_nVarDet
  F_RR_Quad(PP_iVarDetQ) = MATMUL(SG_Vdm_OrthQuad,F_R(PP_iVarDet))
  F_LL_Quad(PP_iVarDetQ) = MATMUL(SG_Vdm_OrthQuad,F_L(PP_iVarDet))
END DO

H_L                      = (U_LL_Quad(QENER)+U_LL_Quad(QPRES))*U_LL_Quad(QSRHO)
H_R                      = (U_RR_Quad(QENER)+U_RR_Quad(QPRES))*U_RR_Quad(QSRHO)
SqrtRho_L                = SQRT(U_LL_Quad(QDENS))
SqrtRho_R                = SQRT(U_RR_Quad(QDENS))
sSqrtRho                 = 1./(SqrtRho_L+SqrtRho_R)
! Roe mean values
RoeVel(1:nQPTotal)            = (SqrtRho_R*U_RR_Quad(QVEL1) + SqrtRho_L*U_LL_Quad(QVEL1))* sSqrtRho
RoeVel(nQPTotal+1:2*nQPTotal)      = (SqrtRho_R*U_RR_Quad(QVEL2) + SqrtRho_L*U_LL_Quad(QVEL2))* sSqrtRho
RoeVel(2*nQPTotal+1:3*nQPTotal)    = (SqrtRho_R*U_RR_Quad(QVEL3) + SqrtRho_L*U_LL_Quad(QVEL3))* sSqrtRho
RoeH                     = (SqrtRho_R*H_R        + SqrtRho_L*H_L)                   * sSqrtRho
DO iQP=1,nQPTotal
  absVel(iQP)            = RoeVel(iQP)*RoeVel(iQP) + RoeVel(nQPTotal+iQP)*RoeVel(nQPTotal+iQP) + RoeVel(2*nQPTotal+iQP)+RoeVel(2*nQPTotal+iQP)
END DO !iQP
Roec                     = SQRT(KappaM1*(RoeH-0.5*absVel))
beta                     = BETA_RIEMANN_H()
c_LL_Quad                = SQRT(kappa*U_LL_Quad(QPRES)*U_LL_Quad(QSRHO))
c_RR_Quad                = SQRT(kappa*U_RR_Quad(QPRES)*U_RR_Quad(QSRHO))
SsL                      = MIN(RoeVel(1:nQPTotal)-Roec,U_LL_Quad(QVEL1) - beta*c_LL_Quad, 0.)
SsR                      = MAX(RoeVel(1:nQPTotal)+Roec,U_RR_Quad(QVEL1) + beta*c_RR_Quad, 0.)

DO iQP=1,nQPTotal
  ! positive supersonic speed
  IF(SsL(iQP) .GE. 0.)THEN
    F_Quad(PP_iQPCons)   = F_LL_Quad(PP_iQPCons)
  ! negative supersonic speed
  ELSEIF(SsR(iQP) .LE. 0.)THEN
    F_Quad(PP_iQPCons)   = F_RR_Quad(PP_iQPCons)
  ! subsonic case
  ELSE
    F_Quad(PP_iQPCons)   = (SsR(iQP)*F_LL_Quad(PP_iQPCons)-SsL(iQP)*F_RR_Quad(PP_iQPCons) +SsL(iQP)*SsR(iQP)*(U_RR_Quad(PP_iQPCons)-U_LL_Quad(PP_iQPCons))) &
                           & /(SsR(iQP)-SsL(iQP))
  END IF
END DO !iQP
! transform back
DO iVarDet=1,PP_nVarDet
  F(PP_iVarDet)          = MATMUL(SG_Vdm_QuadOrth,  F_Quad(PP_iVarDetQ))
END DO

END SUBROUTINE Riemann_HLLE

#ifdef SGdummy
!=================================================================================================================================
!> Harten-Lax-Van-Leer-Einfeldt-Munz Riemann solver
!=================================================================================================================================
PURE SUBROUTINE Riemann_HLLEM(F_L,F_R,U_LL,U_RR,F)
!=================================================================================================================================
! MODULES
USE MOD_PreProc       ,ONLY: PP_nVar,PP_nCoef
USE MOD_EOS_Vars      ,ONLY: Kappa,KappaM1
IMPLICIT NONE
!---------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
                                               !> extended solution vector on the left/right side of the interface
REAL,DIMENSION(PP_2Var),INTENT(IN) :: U_LL,U_RR
                                               !> advection fluxes on the left/right side of the interface
REAL,DIMENSION(PP_nVar),INTENT(IN) :: F_L,F_R
REAL,DIMENSION(PP_nVar),INTENT(OUT):: F        !< resulting Riemann flux
!---------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                                   :: H_L,H_R
REAL                                   :: SqrtRho_L,SqrtRho_R,sSqrtRho,absVel
REAL                                   :: RoeVel(3),RoeH,Roec,RoeDens
REAL                                   :: Ssl,Ssr
REAL                                   :: Alpha(2:4),delta,beta
REAL,DIMENSION(PP_nVar)                :: r2,r3,r4  ! Roe eigenvectors + jump in prims
!=================================================================================================================================
H_L       = TOTALENTHALPY_HE(U_LL)
H_R       = TOTALENTHALPY_HE(U_RR)
SqrtRho_L = SQRT(U_LL(DENS))
SqrtRho_R = SQRT(U_RR(DENS))
sSqrtRho  = 1./(SqrtRho_L+SqrtRho_R)
! Roe mean values
RoeVel    = (SqrtRho_R*U_RR(VELV) + SqrtRho_L*U_LL(VELV)) * sSqrtRho
RoeH      = (SqrtRho_R*H_R        + SqrtRho_L*H_L)        * sSqrtRho
absVel    = DOT_PRODUCT(RoeVel,RoeVel)
Roec      = ROEC_RIEMANN_H(RoeH,RoeVel)
RoeDens   = SQRT(U_LL(DENS)*U_RR(DENS))
! HLLEM flux (positively conservative)
beta=BETA_RIEMANN_H()
SsL=MIN(RoeVel(1)-Roec,U_LL(VEL1) - beta*SPEEDOFSOUND_HE(U_LL), 0.)
SsR=MAX(RoeVel(1)+Roec,U_RR(VEL1) + beta*SPEEDOFSOUND_HE(U_RR), 0.)

! positive supersonic speed
IF(Ssl .GE. 0.)THEN
  F=F_L
! negative supersonic speed
ELSEIF(Ssr .LE. 0.)THEN
  F=F_R
! subsonic case
ELSE
  ! delta
  delta = Roec/(Roec+ABS(0.5*(Ssl+Ssr)))

  ! mean eigenvectors
  Alpha(2)   = (U_RR(DENS)-U_LL(DENS))  - (U_RR(PRES)-U_LL(PRES))/(Roec*Roec)
  Alpha(3:4) = RoeDens*((/U_RR(VEL2),U_RR(VEL3)/) - (/U_LL(VEL2),U_LL(VEL3)/))
  r2 = (/ 1., RoeVel(1), RoeVel(2), RoeVel(3), 0.5*absVel /)
  r3 = (/ 0., 0.,        1.,        0.,        RoeVel(2)  /)
  r4 = (/ 0., 0.,        0.,        1.,        RoeVel(3)  /)

  F=(Ssr*F_L-Ssl*F_R + Ssl*Ssr* &
     (U_RR(CONS)-U_LL(CONS) - delta*(r2*Alpha(2)+r3*Alpha(3)+r4*Alpha(4))))/(Ssr-Ssl)
END IF ! subsonic case
END SUBROUTINE Riemann_HLLEM
#endif /*SGdummy*/

#ifdef SPLIT_DG
!==================================================================================================================================
!> Riemann solver using purely the average fluxes
!==================================================================================================================================
SUBROUTINE Riemann_FluxAverage(F_L,F_R,U_LL,U_RR,F)
! MODULES
USE MOD_PreProc
USE MOD_SplitFlux     ,ONLY: SplitDGSurface_pointer
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
                                                !> extended solution vector on the left/right side of the interface
REAL,DIMENSION(PP_2Var),INTENT(IN) :: U_LL,U_RR
                                                !> advection fluxes on the left/right side of the interface
REAL,DIMENSION(PP_nVar),INTENT(IN) :: F_L,F_R
REAL,DIMENSION(PP_nVar),INTENT(OUT):: F         !< resulting Riemann flux
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!==================================================================================================================================
! get split flux
CALL SplitDGSurface_pointer(U_LL,U_RR,F)

#ifdef DEBUG
! ===============================================================================
! Following dummy calls do suppress compiler warnings of unused Riemann-functions
! ===============================================================================
IF (0.EQ.1) THEN
  WRITE (*,*) F_L,F_R
END IF
#endif /*DEBUG*/
END SUBROUTINE Riemann_FluxAverage
#endif /*SPLIT_DG*/

!==================================================================================================================================
!> Finalize Riemann solver routines
!==================================================================================================================================
SUBROUTINE FinalizeRiemann()
! MODULES
!----------------------------------------------------------------------------------------------------------------------------------!
IMPLICIT NONE
! INPUT / OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!==================================================================================================================================
END SUBROUTINE FinalizeRiemann


END MODULE MOD_Riemann
