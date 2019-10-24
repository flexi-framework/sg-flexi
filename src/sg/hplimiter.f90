!=================================================================================================================================
! Copyright (c) 2016  Prof. Claus-Dieter Munz
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
#if SG_HP_LIMITER
#include "flexi.h"
#include "eos.h"
!==================================================================================================================================
!==================================================================================================================================
MODULE MOD_HPLimiter
!----------------------------------------------------------------------------------------------------------------------------------
! MODULES
IMPLICIT NONE
PRIVATE
SAVE
!----------------------------------------------------------------------------------------------------------------------------------
INTERFACE InitHPLimiter
  MODULE PROCEDURE InitHPLimiter
END INTERFACE

INTERFACE HyperbolicityPreservingLimiter
  MODULE PROCEDURE HyperbolicityPreservingLimiter
END INTERFACE

INTERFACE checkAdmissible
  MODULE PROCEDURE checkAdmissible
END INTERFACE

INTERFACE FinalizeHPLimiter
  MODULE PROCEDURE FinalizeHPLimiter
END INTERFACE

PUBLIC:: InitHPLimiter
PUBLIC:: HyperbolicityPreservingLimiter
PUBLIC:: checkAdmissible
PUBLIC:: FinalizeHPLimiter
!==================================================================================================================================

CONTAINS

!==================================================================================================================================
!> Initialize SG information and SG operators
!==================================================================================================================================
SUBROUTINE InitHPLimiter()
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_SG_Vars
USE MOD_IO_HDF5            ,ONLY: AddToFieldData,FieldOut
USE MOD_Mesh_Vars          ,ONLY: nElems
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!==================================================================================================================================
ALLOCATE(t_HPLimiter(1,0:PP_N,0:PP_N,0:PP_NZ,1:nElems))
t_HPLimiter=0.
ALLOCATE(nLimitedElems(1:nSGElems))
CALL AddToFieldData(FieldOut,(/1,PP_N+1,PP_N+1,PP_NZ+1/),'HypPresLim',(/'t_HPLimiter'/),RealArray=t_HPLimiter)
END SUBROUTINE InitHPLimiter

!==================================================================================================================================
!> Hyperbolicity Preserving Limiter, limits polynomial towards admissible cellmean
!==================================================================================================================================
SUBROUTINE HyperbolicityPreservingLimiter()
! MODULES
USE MOD_PreProc
USE MOD_DG_Vars             ,ONLY: U
USE MOD_Mesh_Vars           ,ONLY: nElems
USE MOD_SG_Vars             ,ONLY: nQPTotal,SG_VdM_OrthQuad,t_HPLimiter
USE MOD_Interpolation_Vars  ,ONLY: wGP
#if FV_ENABLED
USE MOD_FV_Vars             ,ONLY: FV_Elems
#if FV_RECONSTRUCT
USE MOD_FV_Vars             ,ONLY: FV_dx_XI_L,FV_dx_ETA_L
USE MOD_FV_Vars             ,ONLY: FV_dx_XI_R,FV_dx_ETA_R
USE MOD_FV_Vars             ,ONLY: gradUxi,gradUeta,FV_Elems
USE MOD_EOS                 ,ONLY: ConsToPrimDet,PrimtoConsDet
#if PP_dim == 3
USE MOD_FV_Vars             ,ONLY: gradUzeta
USE MOD_FV_Vars             ,ONLY: FV_dx_ZETA_L
USE MOD_FV_Vars             ,ONLY: FV_dx_ZETA_R
#endif
#endif
#endif
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                      :: iElem,iQP,iCoef,i,j,k,iVarDet
REAL                         :: UMean(PP_nVarDet)
REAL                         :: UQuad(PP_nVarDet)
REAL                         :: t,t_loc
#if FV_RECONSTRUCT
REAL                         :: UPrim_Quad(PP_nVarPrimDet)
INTEGER                      :: ii,jj,kk
REAL                         :: UPrim_FV_Quad(PP_nVarPrimDet,0:1,0:1,0:ZDIM(1))
REAL                         :: gradUxi_Quad(PP_nVarPrimDet)
REAL                         :: gradUeta_Quad(PP_nVarPrimDet)
#if PP_dim == 3
REAL                         :: gradUzeta_Quad(PP_nVarPrimDet)
#endif
#endif /*FV_RECONSTRUCT*/
!==================================================================================================================================
t_HPLimiter = 0.
! mean value
DO iElem=1,nElems
#if FV_ENABLED
  IF (FV_Elems(iElem).EQ.1) CYCLE ! FV Elem
#endif
  UMean = 0.
  DO k=0,PP_NZ;DO j=0,PP_N;DO i=0,PP_N
    UMean = UMean + U(SMODE(0),i,j,k,iElem)*wGP(i)*wGP(j)*MERGE(wGP(k),2.,PP_dim.EQ.3)
  END DO; END DO; END DO
  UMean = UMean / 8.
  t=0.
  DO k=0,PP_NZ;DO j=0,PP_N;DO i=0,PP_N
    DO iQP=1,nQPTotal
      DO iVarDet=1,PP_nVarDet
        UQuad(iVarDet)=DOT_PRODUCT(SG_VdM_OrthQuad(iQP,:),U(PP_iVarDet,i,j,k,iElem))
      END DO
      CALL GetTheta(UQuad,UMean,t_loc)
      t=MAX(t,t_loc)
    END DO !iQP
  END DO;END DO;END DO
  IF(t.GT.TINY(1.)) THEN
    DO k=0,PP_NZ;DO j=0,PP_N;DO i=0,PP_N
      U(SMODE(0),i,j,k,iElem) = (1.-t)*U(SMODE(0),i,j,k,iElem) + t*UMean
      DO iCoef=1,PP_nCoef-1
        U(SMODE(iCoef),i,j,k,iElem) = (1.-t)*U(SMODE(iCoef),i,j,k,iElem)
      END DO !iCoef
    END DO;END DO;END DO
  END IF
  t_HPLimiter(:,:,:,:,iElem)=t
END DO !iElem


#if FV_ENABLED
DO iElem=1,nElems
  IF(FV_Elems(iElem).EQ.0) CYCLE
  DO k=0,PP_NZ;DO j=0,PP_N;DO i=0,PP_N
    t=0.
    UMean = U(SMODE(0),i,j,k,iElem)
    DO iQP=1,nQPTotal
      DO iVarDet=1,PP_nVarDet
        UQuad(iVarDet)=DOT_PRODUCT(SG_VdM_OrthQuad(iQP,:),U(PP_iVarDet,i,j,k,iElem))
      END DO
#if FV_RECONSTRUCT
      CALL ConsToPrimDet(UPrim_Quad,UQuad)
      DO iVarDet=1,PP_nVarPrimDet
        gradUxi_Quad(  iVarDet)=DOT_PRODUCT(SG_VdM_OrthQuad(iQP,:),  gradUxi(PP_iVarDet,i,k,j,iElem))
        gradUeta_Quad( iVarDet)=DOT_PRODUCT(SG_VdM_OrthQuad(iQP,:), gradUeta(PP_iVarDet,i,k,j,iElem))
#if PP_dim == 3
        gradUzeta_Quad(iVarDet)=DOT_PRODUCT(SG_VdM_OrthQuad(iQP,:),gradUzeta(PP_iVarDet,i,k,j,iElem))
#endif
        UPrim_FV_Quad(iVarDet,:,:,:) = UPrim_Quad(iVarDet)
        UPrim_FV_Quad(iVarDet,0,:,:) = UPrim_FV_Quad(iVarDet,0,:,:) -   gradUxi_Quad(iVarDet) *   FV_dx_XI_L(j,k,i,iElem)
        UPrim_FV_Quad(iVarDet,1,:,:) = UPrim_FV_Quad(iVarDet,1,:,:) +   gradUxi_Quad(iVarDet) *   FV_dx_XI_R(j,k,i,iElem)
        UPrim_FV_Quad(iVarDet,:,0,:) = UPrim_FV_Quad(iVarDet,:,0,:) -  gradUeta_Quad(iVarDet) *  FV_dx_ETA_L(j,k,i,iElem)
        UPrim_FV_Quad(iVarDet,:,1,:) = UPrim_FV_Quad(iVarDet,:,1,:) +  gradUeta_Quad(iVarDet) *  FV_dx_ETA_R(j,k,i,iElem)
#if PP_dim == 3
        UPrim_FV_Quad(iVarDet,:,:,0) = UPrim_FV_Quad(iVarDet,:,:,0) - gradUzeta_Quad(iVarDet) * FV_dx_ZETA_L(j,k,i,iElem)
        UPrim_FV_Quad(iVarDet,:,:,1) = UPrim_FV_Quad(iVarDet,:,:,1) + gradUzeta_Quad(iVarDet) * FV_dx_ZETA_R(j,k,i,iElem)
#endif
      END DO
      DO kk=0,ZDIM(1);DO jj=0,1;DO ii=0,1
        CALL PrimtoConsDet(UPrim_FV_Quad(:,ii,jj,kk),UQuad)
        CALL GetTheta(UQuad,UMean,t_loc)
        t=MAX(t,t_loc)
      END DO;END DO;END DO
#else /*FV_RECONSTRUCT*/
      CALL GetTheta(UQuad,UMean,t_loc)
      t=MAX(t,t_loc)
#endif /*FV_RECONSTRUCT*/
    END DO !iQP
    IF(t.GT.TINY(1.)) THEN
      DO iCoef=1,PP_nCoef-1
        U(SMODE(iCoef),i,j,k,iElem) = (1.-t)*U(SMODE(iCoef),i,j,k,iElem)
      END DO !iCoef
#if FV_RECONSTRUCT
      gradUxi  (:,i,k,j,iElem) = (1.-t) *   gradUxi(:,i,k,j,iElem)
      gradUeta (:,i,k,j,iElem) = (1.-t) *  gradUeta(:,i,k,j,iElem)
      gradUzeta(:,i,k,j,iElem) = (1.-t) * gradUzeta(:,i,k,j,iElem)
#endif /*FV_RECONSTRUCT*/
    END IF
    t_HPLimiter(1,i,j,k,iElem)=t
  END DO;END DO;END DO
END DO
#endif /*FV_ENABLED*/
END SUBROUTINE HyperbolicityPreservingLimiter



!==================================================================================================================================
!> Computes thetha, such that theta*U+(1-theta)*cellmean is admissible
!==================================================================================================================================
SUBROUTINE GetTheta(UQuad,UMean,t_out)
! MODULES
USE MOD_PreProc
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
REAL,INTENT(IN)         :: UQuad(PP_nVarDet)
REAL,INTENT(IN)         :: UMean(PP_nVarDet)
REAL,INTENT(OUT)        :: t_out
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                    :: t1,t2,t3,a,b,c,r
REAL                    :: eps=1.E-8
REAL                    :: UDiff(PP_nVarDet)
REAL                    :: p_sGammaM1
!==================================================================================================================================
p_sGammaM1=UQuad(DENER)-0.5*DOT_PRODUCT(UQuad(DMMV2),UQuad(DMMV2))/UQuad(DDENS)
IF(UQuad(DDENS).GT.eps.AND.p_sGammaM1.GT.eps) THEN  !TODO: Is eps relative or absolute?
  t_out=0.
  RETURN
END IF
UDiff=UQuad-UMean

t1 = (UQuad(DDENS)-eps) / UDiff(DDENS)

a  =                  UDiff(DENER)*UDiff(DDENS)  - 0.5*DOT_PRODUCT(UDiff(DMMV2),UDiff(DMMV2))
b  =      DOT_PRODUCT(UQuad(DMMV2),UDiff(DMMV2)) - UQuad(DENER)*UDiff(DDENS) - UQuad(DDENS)*UDiff(DENER)
c  = -0.5*DOT_PRODUCT(UQuad(DMMV2),UQuad(DMMV2)) + UQuad(DENER)*UQuad(DDENS) - eps

r = -0.5*(b + SIGN(1.,b)*SQRT(b*b-4.*a*c))
IF(a.EQ.0.) THEN
  t3 = 10.
  t2 = -c/b
ELSE
  t2 = r/a
  t3 = c/r
END IF

IF(ISNAN(t1).OR.(t1.GT.1.).OR.(t1.LT.0.)) t1=0.
IF(ISNAN(t2).OR.(t2.GT.1.).OR.(t2.LT.0.)) t2=0.
IF(ISNAN(t3).OR.(t3.GT.1.).OR.(t3.LT.0.)) t3=0.
t_out = MAX(t1,t2,t3)
END SUBROUTINE GetTheta

!==================================================================================================================================
!> Check if cellmean is admissible, i.e. positive density and pressure
!==================================================================================================================================
SUBROUTINE checkAdmissible()
! MODULES
USE MOD_PreProc
USE MOD_Globals
USE MOD_Interpolation_Vars  ,ONLY: wGP
USE MOD_Mesh_Vars           ,ONLY: nElems
USE MOD_DG_Vars             ,ONLY: U
#if FV_ENABLED
USE MOD_FV_Vars             ,ONLY: FV_Elems
#endif
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                      :: iElem,i,j,k
REAL                         :: UMean(PP_nVarDet)
REAL                         :: pressure
!==================================================================================================================================
DO iElem=1,nElems
#if FV_ENABLED
  IF (FV_Elems(iElem).EQ.1) CYCLE
#endif
  DO k=0,PP_NZ;DO j=0,PP_N;DO i=0,PP_N
    UMean = 0.
    UMean = UMean + U(SMODE(0),i,j,k,iElem)*wGP(i)*wGP(j)*MERGE(wGP(k),2.,PP_dim.EQ.3)
  END DO; END DO; END DO
  UMean = UMean / 8.
  pressure= UMean(DENER)-0.5*DOT_PRODUCT(UMean(DMMV2),UMean(DMMV2))/UMean(DDENS)
  IF( (pressure.LT.0.).OR.(UMean(DDENS).LT.0.)) THEN
    ERRWRITE(*,'(A)') 'Cellmean not admissible'
  END IF
END DO
#if FV_ENABLED
DO iElem=1,nElems
  IF (FV_Elems(iElem).EQ.0) CYCLE
  DO k=0,PP_NZ;DO j=0,PP_N;DO i=0,PP_N
    UMean = U(SMODE(0),i,j,k,iElem)! FV Elem
    pressure= UMean(DENER)-0.5*DOT_PRODUCT(UMean(DMMV2),UMean(DMMV2))/UMean(DDENS)
    IF( (pressure.LT.0.).OR.(UMean(DDENS).LT.0.)) THEN
      ERRWRITE(*,'(A)') 'Cellmean not admissible'
    END IF
  END DO; END DO; END DO
END DO
#endif
END SUBROUTINE checkAdmissible

!==================================================================================================================================
!> Initialize SG information and SG operators
!==================================================================================================================================
SUBROUTINE FinalizeHPLimiter()
! MODULES
USE MOD_SG_Vars
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!==================================================================================================================================
SDEALLOCATE(t_HPLimiter)
SDEALLOCATE(nLimitedElems)
END SUBROUTINE FinalizeHPLimiter

END MODULE MOD_HPLimiter
#endif
