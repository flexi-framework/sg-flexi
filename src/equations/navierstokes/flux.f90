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
#include "eos.h"

!==================================================================================================================================
!> Contains the routine EvalFlux3D which computes the complete NSE flux f,g,h for all DOFs in one Element: used in volume integral
!> Contains the routine EvalFlux1D_Adv which computes the Euler flux f for all DOFs of one Element side: used in Riemann_Adv
!==================================================================================================================================
MODULE MOD_Flux
! MODULES
IMPLICIT NONE
PRIVATE
!----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------

!INTERFACE EvalFlux3D
  !MODULE PROCEDURE EvalFlux3D
!END INTERFACE

!INTERFACE EvalEulerFlux1D
  !MODULE PROCEDURE EvalEulerFlux1D
!END INTERFACE
!INTERFACE EvalEulerFlux1D_fast
  !MODULE PROCEDURE EvalEulerFlux1D_fast
!END INTERFACE

#if PARABOLIC
INTERFACE EvalDiffFlux3D
  MODULE PROCEDURE EvalDiffFlux3D
END INTERFACE

INTERFACE EvalDiffFlux2D
  MODULE PROCEDURE EvalDiffFlux2D
END INTERFACE

!INTERFACE EvalDiffFlux1D !only for testing
!  MODULE PROCEDURE EvalDiffFlux1D
!END INTERFACE
#endif /*PARABOLIC*/

PUBLIC::EvalFlux3D,EvalEulerFlux1D_fast
#ifdef SGdummy
PUBLIC::EvalEulerFlux1D
#endif /* SGdummy */
#if PARABOLIC
PUBLIC::EvalDiffFlux3D,EvalDiffFlux2D!,EvalDiffFlux1D
#endif /*PARABOLIC*/
!==================================================================================================================================

CONTAINS

!==================================================================================================================================
!> Compute Navier-Stokes fluxes using the conservative variables and derivatives for every volume Gauss point.
!==================================================================================================================================
SUBROUTINE EvalFlux3D(Nloc,U,UPrim,f,g,h)
! MODULES
USE MOD_PreProc
USE MOD_SG_Operators, ONLY: SG_Product
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
INTEGER,INTENT(IN) :: Nloc
REAL,DIMENSION(PP_nVar    ,0:Nloc,0:Nloc,0:ZDIM(Nloc)),INTENT(IN)  :: U        !< Conservative solution
REAL,DIMENSION(PP_nVarPrim,0:Nloc,0:Nloc,0:ZDIM(Nloc)),INTENT(IN)  :: UPrim    !< Primitive solution
REAL,DIMENSION(PP_nVar    ,0:Nloc,0:Nloc,0:ZDIM(Nloc)),INTENT(OUT) :: f        !< Cartesian flux in x (iVar,i,j,k)
REAL,DIMENSION(PP_nVar    ,0:Nloc,0:Nloc,0:ZDIM(Nloc)),INTENT(OUT) :: g        !< Cartesian flux in y (iVar,i,j,k)
REAL,DIMENSION(PP_nVar    ,0:Nloc,0:Nloc,0:ZDIM(Nloc)),INTENT(OUT) :: h        !< Cartesian flux in z (iVar,i,j,k)
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                :: Ep(1:PP_nCoef)
INTEGER             :: i,j,k
!==================================================================================================================================
DO k=0,ZDIM(Nloc);  DO j=0,Nloc; DO i=0,Nloc
  ! auxiliary variables
  Ep   = U(ENER,i,j,k) + UPrim(PPRES,i,j,k)
#if PP_dim==3
  ! Euler part
  ! Euler fluxes x-direction
  f(SGI1,i,j,k)=U(MOM1,i,j,k)                                                   ! rho*u
  f(SGI2,i,j,k)=SG_Product(U(MOM1,i,j,k),UPrim(PVEL1,i,j,k))+UPrim(PPRES,i,j,k) ! rho*u²+p
  f(SGI3,i,j,k)=SG_Product(U(MOM1,i,j,k),UPrim(PVEL2,i,j,k))                    ! rho*u*v
  f(SGI4,i,j,k)=SG_Product(U(MOM1,i,j,k),UPrim(PVEL3,i,j,k))                    ! rho*u*w
  f(SGI5,i,j,k)=SG_Product(Ep,UPrim(PVEL1,i,j,k))                               ! (rho*e+p)*u
  ! Euler fluxes y-direction
  g(SGI1,i,j,k)=U(MOM2,i,j,k)                                                   ! rho*v
  g(SGI2,i,j,k)=f(SGI3,i,j,k)                                                   ! rho*u*v
  g(SGI3,i,j,k)=SG_Product(U(MOM2,i,j,k),UPrim(PVEL2,i,j,k))+UPrim(PPRES,i,j,k) ! rho*v²+p
  g(SGI4,i,j,k)=SG_Product(U(MOM2,i,j,k),UPrim(PVEL3,i,j,k))                    ! rho*v*w
  g(SGI5,i,j,k)=SG_Product(Ep,UPrim(PVEL2,i,j,k))                               ! (rho*e+p)*v
  ! Euler fluxes z-direction
  h(SGI1,i,j,k)=U(MOM3,i,j,k)                                                   ! rho*v
  h(SGI2,i,j,k)=f(SGI4,i,j,k)                                                   ! rho*u*w
  h(SGI3,i,j,k)=g(SGI4,i,j,k)                                                   ! rho*v*w
  h(SGI4,i,j,k)=SG_Product(U(MOM3,i,j,k),UPrim(PVEL3,i,j,k))+UPrim(PPRES,i,j,k) ! rho*v²+p
  h(SGI5,i,j,k)=SG_Product(Ep,UPrim(PVEL3,i,j,k))                               ! (rho*e+p)*w
#else

  ! Euler part
  ! Euler fluxes x-direction
  f(SGI1,i,j,k)=U(MOM1,i,j,k)                                                   ! rho*u
  f(SGI2,i,j,k)=SG_Product(U(MOM1,i,j,k),UPrim(PVEL1,i,j,k))+UPrim(PPRES,i,j,k) ! rho*u²+p
  f(SGI3,i,j,k)=SG_Product(U(MOM1,i,j,k),UPrim(PVEL2,i,j,k))                    ! rho*u*v
  f(SGI4,i,j,k)=0.
  f(SGI5,i,j,k)=SG_Product(Ep,UPrim(PVEL1,i,j,k))                               ! (rho*e+p)*u
  ! Euler fluxes y-direction
  g(SGI1,i,j,k)=U(MOM2,i,j,k)                                                   ! rho*v
  g(SGI2,i,j,k)=f(SGI3,i,j,k)                                                   ! rho*u*v
  g(SGI3,i,j,k)=SG_Product(U(MOM2,i,j,k),UPrim(PVEL2,i,j,k))+UPrim(PPRES,i,j,k) ! rho*v²+p
  g(SGI4,i,j,k)=0.
  g(SGI5,i,j,k)=SG_Product(Ep,UPrim(PVEL2,i,j,k))                               ! (rho*e+p)*v
  ! Euler fluxes z-direction
  h(:,i,j,k)=0.
#endif
END DO; END DO; END DO ! i,j,k
END SUBROUTINE EvalFlux3D


#if PARABOLIC
!==================================================================================================================================
!> Compute Navier-Stokes diffusive fluxes using the conservative variables and derivatives for every volume Gauss point.
!==================================================================================================================================
SUBROUTINE EvalDiffFlux3D(UPrim,gradUx,gradUy,gradUz,f,g,h,iElem)
! MODULES
USE MOD_PreProc
USE MOD_Equation_Vars,ONLY: s23,s43
USE MOD_EOS_Vars,     ONLY: cp,Pr
USE MOD_Viscosity
#ifdef EDDYVISCOSITY
USE MOD_EddyVisc_Vars,ONLY: muSGS,muSGSmax,PrSGS,eddyViscosity
USE MOD_TimeDisc_Vars,ONLY: CurrentStage
#endif
USE MOD_SG_Operators, ONLY: SG_Product
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
REAL,DIMENSION(PP_nVarPrim,0:PP_N,0:PP_N,0:PP_NZ),INTENT(IN)  :: UPrim                !< Solution vector
REAL,DIMENSION(PP_nVarPrim,0:PP_N,0:PP_N,0:PP_NZ),INTENT(IN)  :: gradUx,gradUy,gradUz !< Gradients in x,y,z directions
REAL,DIMENSION(PP_nVar    ,0:PP_N,0:PP_N,0:PP_NZ),INTENT(OUT) :: f,g,h                !< Cartesian fluxes (iVar,i,j,k)
INTEGER, INTENT(IN)                                           :: iELem                !< element index in global array
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                :: muS(PP_nCoef)
REAL                :: v1(PP_nCoef),v2(PP_nCoef)
REAL                :: tau_xx(PP_nCoef),tau_yy(PP_nCoef),tau_xy(PP_nCoef)
REAL                :: gradT1(PP_nCoef),gradT2(PP_nCoef),lambda(PP_nCoef),prim(PP_nVarPrim)
INTEGER             :: i,j,k
#if PP_dim==3
REAL                :: v3(PP_nCoef)
REAL                :: tau_zz(PP_nCoef),tau_xz(PP_nCoef),tau_yz(PP_nCoef)
REAL                :: gradT3(PP_nCoef)
#endif
!==================================================================================================================================
#if EDDYVISCOSITY
IF(CurrentStage.EQ.1) muSGSmax(iElem)=0.
#endif
DO k=0,PP_NZ;  DO j=0,PP_N; DO i=0,PP_N
  prim = UPrim(   :,i,j,k)
  v1   = UPrim(SGI2,i,j,k)
  v2   = UPrim(SGI3,i,j,k)
  ! Viscous part
  ! ideal gas law
  !TODO
  muS=VISCOSITY_PRIM(prim) !returns mu0
  lambda=THERMAL_CONDUCTIVITY_H(muS) !returns constant mu0*cp/Pr
  !Add turbulent sub grid scale viscosity to mu
#ifdef EDDYVISCOSITY
  IF(CurrentStage.EQ.1) THEN
    CALL eddyViscosity(iElem,i,j,k,muSGS(1,i,j,k,iElem))
  END IF
  muS = muS + muSGS(1,i,j,k,iElem)
  lambda = lambda + muSGS(1,i,j,k,iElem)*cp/PrSGS
#endif

#if PP_dim==3
  v3   = UPrim(SGI4,i,j,k)
  ! gradients of primitive variables are directly available gradU = (/ drho, dv1, dv2, dv3, dT /)

  ! viscous fluxes in x-direction
  tau_xx=SG_Product(muS,( s43*gradUx(SGI2,i,j,k)-s23*gradUy(SGI3,i,j,k)-s23*gradUz(SGI4,i,j,k))) !  4/3*mu*u_x-2/3*mu*v_y -2/3*mu*w*z
  tau_yy=SG_Product(muS,(-s23*gradUx(SGI2,i,j,k)+s43*gradUy(SGI3,i,j,k)-s23*gradUz(SGI4,i,j,k))) ! -2/3*mu*u_x+4/3*mu*v_y -2/3*mu*w*z
  tau_zz=SG_Product(muS,(-s23*gradUx(SGI2,i,j,k)-s23*gradUy(SGI3,i,j,k)+s43*gradUz(SGI4,i,j,k))) ! -2/3*mu*u_x-2/3*mu*v_y +4/3*mu*w*z
  tau_xy=SG_Product(muS,(gradUy(SGI2,i,j,k)+gradUx(SGI3,i,j,k)))                                 ! mu*(u_y+v_x)
  tau_xz=SG_Product(muS,(gradUz(SGI2,i,j,k)+gradUx(SGI4,i,j,k)))                                 ! mu*(u_z+w_x)
  tau_yz=SG_Product(muS,(gradUz(SGI3,i,j,k)+gradUy(SGI4,i,j,k)))                                 ! mu*(y_z+w_y)

  gradT1=gradUx(SGI6,i,j,k)
  gradT2=gradUy(SGI6,i,j,k)
  gradT3=gradUz(SGI6,i,j,k)

  f(SGI1,i,j,k)=0.
  f(SGI2,i,j,k)=-tau_xx                                       ! F_euler-4/3*mu*u_x+2/3*mu*(v_y+w_z)
  f(SGI3,i,j,k)=-tau_xy                                       ! F_euler-mu*(u_y+v_x)
  f(SGI4,i,j,k)=-tau_xz                                       ! F_euler-mu*(u_z+w_x)
  f(SGI5,i,j,k)=-SG_Product(tau_xx,v1)-SG_Product(tau_xy,v2)-SG_Product(tau_xz,v3)-SG_Product(lambda,gradT1)
                                                              ! F_euler-(tau_xx*u+tau_xy*v+tau_xz*w-q_x) q_x=-lambda*T_x
  ! viscous fluxes in y-direction
  g(SGI1,i,j,k)=0.
  g(SGI2,i,j,k)=-tau_xy                                       ! F_euler-mu*(u_y+v_x)
  g(SGI3,i,j,k)=-tau_yy                                       ! F_euler-4/3*mu*v_y+2/3*mu*(u_x+w_z)
  g(SGI4,i,j,k)=-tau_yz                                       ! F_euler-mu*(y_z+w_y)
  g(SGI5,i,j,k)=-SG_Product(tau_xy,v1)-SG_Product(tau_yy,v2)-SG_Product(tau_yz,v3)-SG_Product(lambda,gradT2)
                                                              ! F_euler-(tau_yx*u+tau_yy*v+tau_yz*w-q_y) q_y=-lambda*T_y
  ! viscous fluxes in z-direction
  h(SGI1,i,j,k)=0.
  h(SGI2,i,j,k)=-tau_xz                                       ! F_euler-mu*(u_z+w_x)
  h(SGI3,i,j,k)=-tau_yz                                       ! F_euler-mu*(y_z+w_y)
  h(SGI4,i,j,k)=-tau_zz                                       ! F_euler-4/3*mu*w_z+2/3*mu*(u_x+v_y)
  h(SGI5,i,j,k)=-SG_Product(tau_xz,v1)-SG_Product(tau_yz,v2)-SG_Product(tau_zz,v3)-SG_Product(lambda,gradT3)
                                                              ! F_euler-(tau_zx*u+tau_zy*v+tau_zz*w-q_z) q_z=-lambda*T_z
#else
  ! gradients of primitive variables are directly available gradU = (/ drho, dv1, dv2, dv3, dT /)

  ! viscous fluxes in x-direction
  tau_xx=SG_Product(muS,( s43*gradUx(SGI2,i,j,k)-s23*gradUy(SGI3,i,j,k))) ! 4/3*mu*u_x-2/3*mu*v_y -2/3*mu*w*z
  tau_yy=SG_Product(muS,(-s23*gradUx(SGI2,i,j,k)+s43*gradUy(SGI3,i,j,k))) !-2/3*mu*u_x+4/3*mu*v_y -2/3*mu*w*z
  tau_xy=SG_Product(muS,(gradUy(SGI2,i,j,k)+gradUx(SGI3,i,j,k))       )   !mu*(u_y+v_x)

  gradT1=gradUx(SGI6,i,j,k)
  gradT2=gradUy(SGI6,i,j,k)

  f(SGI1,i,j,k)=0.
  f(SGI2,i,j,k)=-tau_xx                                       ! F_euler-4/3*mu*u_x+2/3*mu*(v_y+w_z)
  f(SGI3,i,j,k)=-tau_xy                                       ! F_euler-mu*(u_y+v_x)
  f(SGI4,i,j,k)=0.
  f(SGI5,i,j,k)=-SG_Product(tau_xx,v1)-SG_Product(tau_xy,v2)-SG_Product(lambda,gradT1)
                                                              ! F_euler-(tau_xx*u+tau_xy*v+tau_xz*w-q_x) q_x=-lambda*T_x
  ! viscous fluxes in y-direction
  g(SGI1,i,j,k)=0.
  g(SGI2,i,j,k)=-tau_xy                                       ! F_euler-mu*(u_y+v_x)
  g(SGI3,i,j,k)=-tau_yy                                       ! F_euler-4/3*mu*v_y+2/3*mu*(u_x+w_z)
  g(SGI4,i,j,k)=0.
  g(SGI5,i,j,k)=-SG_Product(tau_xy,v1)-SG_Product(tau_yy,v2)-SG_Product(lambda,gradT2)
                                                              ! F_euler-(tau_yx*u+tau_yy*v+tau_yz*w-q_y) q_y=-lambda*T_y
  ! viscous fluxes in z-direction
  h(:,i,j,k)=0.
#endif

END DO; END DO; END DO ! i,j,k

#ifdef DEBUG
! ===============================================================================
! Following dummy calls do suppress compiler warnings of unused Riemann-functions
! ===============================================================================
IF (0.EQ.1) THEN
  WRITE (*,*) iElem,gradUz
END IF
#endif
END SUBROUTINE EvalDiffFlux3D
#endif /*PARABOLIC*/


#ifdef SGdummy
!==================================================================================================================================
!> Computes 1D Euler flux using the conservative variables.
!==================================================================================================================================
PURE SUBROUTINE EvalEulerFlux1D(U,F)
! MODULES
USE MOD_PreProc
USE MOD_EOS_Vars ,ONLY:KappaM1
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
REAL,INTENT(IN)     :: U(PP_nVar)   !< vector of conservative variables
REAL,INTENT(OUT)    :: F(PP_nVar)   !< Cartesian flux in "x" direction
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                :: UE(PP_2Var)  !< auxiliary variables
!==================================================================================================================================
! auxiliary variables
! TODO: ATTENTION: Temperature of UE not filled!!!
UE(CONS)=U
UE(SRHO)=1./UE(DENS)
UE(VELV)=VELOCITY_HE(UE)
UE(PRES)=PRESSURE_HE(UE)
! Euler fluxes x-direction
F(SGI1)= U(MOM1)                     ! rho*u
F(SGI2)= U(MOM1)*UE(VEL1)+UE(PRES)   ! rho*u²+p
F(SGI3)= U(MOM1)*UE(VEL2)            ! rho*u*v
#if PP_dim==3
F(SGI4)= U(MOM1)*UE(VEL3)            ! rho*u*w
#else
F(SGI4)=0.
#endif
F(SGI5)=(U(ENER)+UE(PRES))*UE(VEL1)  ! (rho*e+p)*u
END SUBROUTINE EvalEulerFlux1D
#endif /* SGdummy */

!==================================================================================================================================
!> Computes 1D Euler flux using the conservative and primitive variables (for better performance)
!==================================================================================================================================
PURE SUBROUTINE EvalEulerFlux1D_fast(U,F)
! MODULES
USE MOD_PreProc
USE MOD_SG_Operators ,ONLY: SG_Product
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
REAL,INTENT(IN)     :: U(PP_2Var) !< vector of conservative and primitive variables
REAL,INTENT(OUT)    :: F(PP_nVar) !< Cartesian flux in "x" direction
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!==================================================================================================================================
! Euler fluxes x-direction
F(SGI1)= U(MOM1)                              ! rho*u
F(SGI2)= SG_Product(U(MOM1),U(VEL1))+U(PRES)  ! rho*u²+p
F(SGI3)= SG_Product(U(MOM1),U(VEL2))          ! rho*u*v
#if PP_dim==3
F(SGI4)= SG_Product(U(MOM1),U(VEL3))          ! rho*u*w
#else
F(SGI4)= 0.
#endif
F(SGI5)=SG_Product((U(ENER)+U(PRES)),U(VEL1)) ! (:,rho*e+p)*u
END SUBROUTINE EvalEulerFlux1D_fast


#if PARABOLIC
!==================================================================================================================================
!> Compute Navier-Stokes diffusive fluxes using the conservative variables and derivatives for every volume Gauss point.
!==================================================================================================================================
SUBROUTINE EvalDiffFlux2D(Nloc,f,g,h,UPrim_Face,gradUx_Face,gradUy_Face,gradUz_Face&
#ifdef EDDYVISCOSITY
                         ,muSGS_Face&
#endif
                         )
! MODULES
USE MOD_PreProc
USE MOD_Equation_Vars,ONLY: s43,s23
USE MOD_Viscosity
USE MOD_EOS_Vars,     ONLY: cp,Pr
#ifdef EDDYVISCOSITY
USE MOD_EddyVisc_Vars,ONLY: PrSGS
#endif
USE MOD_SG_Operators, ONLY: SG_Product
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
INTEGER,INTENT(IN) :: Nloc                                      !< Polynomial degree of input solution
REAL,INTENT(IN)    :: UPrim_Face( PP_nVarPrim,0:Nloc,0:ZDIM(Nloc))!< U_Face(iVar,i,j,k)
REAL,INTENT(IN)    :: gradUx_Face(PP_nVarPrim,0:Nloc,0:ZDIM(Nloc))!< gradUx_Face(iVar,j,k)
REAL,INTENT(IN)    :: gradUy_Face(PP_nVarPrim,0:Nloc,0:ZDIM(Nloc))!< gradUy_Face(iVar,i,k)
REAL,INTENT(IN)    :: gradUz_Face(PP_nVarPrim,0:Nloc,0:ZDIM(Nloc))!< gradUz_Face(iVar,i,j)
REAL,INTENT(OUT)   :: f(PP_nVar,0:Nloc,0:ZDIM(Nloc))              !< Cartesian fluxes (iVar,i,j)
REAL,INTENT(OUT)   :: g(PP_nVar,0:Nloc,0:ZDIM(Nloc))              !< Cartesian fluxes (iVar,i,j)
REAL,INTENT(OUT)   :: h(PP_nVar,0:Nloc,0:ZDIM(Nloc))              !< Cartesian fluxes (iVar,i,j)
#ifdef EDDYVISCOSITY
REAL,INTENT(IN)    :: muSGS_Face(1,0:Nloc,0:ZDIM(Nloc))
#endif
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL    :: muS(PP_nCoef)
REAL    :: v1(PP_nCoef),v2(PP_nCoef)
REAL    :: tau_xx(PP_nCoef),tau_yy(PP_nCoef),tau_xy(PP_nCoef)
REAL    :: gradT1(PP_nCoef),gradT2(PP_nCoef),lambda(PP_nCoef),prim(PP_nVarPrim)
#if PP_dim==3
REAL    :: tau_zz(PP_nCoef),tau_xz(PP_nCoef),tau_yz(PP_nCoef)
REAL    :: v3(PP_nCoef)
REAL    :: gradT3(PP_nCoef)
#endif
INTEGER :: i,j
!==================================================================================================================================
DO j=0,ZDIM(Nloc) ; DO i=0,Nloc
  prim = UPrim_Face(:,i,j)
  v1=UPrim_Face(SGI2,i,j)
  v2=UPrim_Face(SGI3,i,j)
  ! Viscous part
  ! ideal gas law
  muS=VISCOSITY_PRIM(prim)
  lambda=THERMAL_CONDUCTIVITY_H(muS)

  ! auxiliary variables
  ! In previous versions gradients of conservative variables had been used, see Git commit b984f2895121e236ce24c149ad15615180995b00
  ! gradients of primitive variables are directly available gradU = (/ drho, dv1, dv2, dv3, dp, dT /)
#ifdef EDDYVISCOSITY
  muS = muS + muSGS_Face(1,i,j)
  lambda = lambda + muSGS_Face(1,i,j)*cp/PrSGS
#endif

#if PP_dim==3
  v3=UPrim_Face(SGI4,i,j)
  ! gradients of primitive variables are directly available gradU = (/ drho, dv1, dv2, dv3, dp, dT /)

  tau_xx=SG_Product(muS, s43*gradUx_Face(SGI2,i,j)-s23*gradUy_Face(SGI3,i,j)-s23*gradUz_Face(SGI4,i,j))   !  4/3*mu*u_x-2/3*mu*v_y -2/3*mu*w*z
  tau_yy=SG_Product(muS,(-s23*gradUx_Face(SGI2,i,j)+s43*gradUy_Face(SGI3,i,j)-s23*gradUz_Face(SGI4,i,j))) ! -2/3*mu*u_x+4/3*mu*v_y -2/3*mu*w*z
  tau_zz=SG_Product(muS,(-s23*gradUx_Face(SGI2,i,j)-s23*gradUy_Face(SGI3,i,j)+s43*gradUz_Face(SGI4,i,j))) ! -2/3*mu*u_x-2/3*mu*v_y +4/3*mu*w*z
  tau_xy=SG_Product(muS,(gradUy_Face(SGI2,i,j)+gradUx_Face(SGI3,i,j)))                                    ! mu*(u_y+v_x)
  tau_xz=SG_Product(muS,(gradUz_Face(SGI2,i,j)+gradUx_Face(SGI4,i,j)))                                    ! mu*(u_z+w_x)
  tau_yz=SG_Product(muS,(gradUz_Face(SGI3,i,j)+gradUy_Face(SGI4,i,j)))                                    ! mu*(y_z+w_y)
  gradT3=gradUz_Face(SGI6,i,j)
  gradT1=gradUx_Face(SGI6,i,j)
  gradT2=gradUy_Face(SGI6,i,j)
  ! viscous fluxes in x-direction
  f(SGI1,i,j)=0.
  f(SGI2,i,j)=-tau_xx                                         ! F_euler-4/3*mu*u_x+2/3*mu*(v_y+w_z)
  f(SGI3,i,j)=-tau_xy                                         ! F_euler-mu*(u_y+v_x)
  f(SGI4,i,j)=-tau_xz                                         ! F_euler-mu*(u_z+w_x)
  f(SGI5,i,j)=-SG_Product(tau_xx,v1)-SG_Product(tau_xy,v2)-SG_Product(tau_xz,v3)-SG_Product(lambda,gradT1)
                                                              ! F_euler-(tau_xx*u+tau_xy*v+tau_xz*w-q_x) q_x=-lambda*T_x
  ! viscous fluxes in y-direction
  g(SGI1,i,j)=0.
  g(SGI2,i,j)=-tau_xy                                         ! F_euler-mu*(u_y+v_x)
  g(SGI3,i,j)=-tau_yy                                         ! F_euler-4/3*mu*v_y+2/3*mu*(u_x+w_z)
  g(SGI4,i,j)=-tau_yz                                         ! F_euler-mu*(y_z+w_y)
  g(SGI5,i,j)=-SG_Product(tau_xy,v1)-SG_Product(tau_yy,v2)-SG_Product(tau_yz,v3)-SG_Product(lambda,gradT2)
                                                              ! F_euler-(tau_yx*u+tau_yy*v+tau_yz*w-q_y) q_y=-lambda*T_y
  ! viscous fluxes in z-direction
  h(SGI1,i,j)=0.
  h(SGI2,i,j)=-tau_xz                                         ! F_euler-mu*(u_z+w_x)
  h(SGI3,i,j)=-tau_yz                                         ! F_euler-mu*(y_z+w_y)
  h(SGI4,i,j)=-tau_zz                                         ! F_euler-4/3*mu*w_z+2/3*mu*(u_x+v_y)
  h(SGI5,i,j)=-SG_Product(tau_xz,v1)-SG_Product(tau_yz,v2)-SG_Product(tau_zz,v3)-SG_Product(lambda,gradT3)
                                                              ! F_euler-(tau_zx*u+tau_zy*v+tau_zz*w-q_z) q_z=-lambda*T_z
#else
  ! gradients of primitive variables are directly available gradU = (/ drho, dv1, dv2, dv3, dp, dT /)

  tau_xx=SG_Product(muS,( s43*gradUx_Face(SGI2,i,j)-s23*gradUy_Face(SGI3,i,j))) ! 4/3*mu*u_x-2/3*mu*v_y -2/3*mu*w*z
  tau_yy=SG_Product(muS,(-s23*gradUx_Face(SGI2,i,j)+s43*gradUy_Face(SGI3,i,j))) !-2/3*mu*u_x+4/3*mu*v_y -2/3*mu*w*z
  tau_xy=SG_Product(muS,(gradUy_Face(SGI2,i,j)+gradUx_Face(SGI3,i,j)))          !mu*(u_y+v_x)
  gradT1=gradUx_Face(SGI6,i,j)
  gradT2=gradUy_Face(SGI6,i,j)
  ! viscous fluxes in x-direction
  f(SGI1,i,j)=0.
  f(SGI2,i,j)=-tau_xx                                         ! F_euler-4/3*mu*u_x+2/3*mu*(v_y+w_z)
  f(SGI3,i,j)=-tau_xy                                         ! F_euler-mu*(u_y+v_x)
  f(SGI4,i,j)=0.
  f(SGI5,i,j)=-SG_Product(tau_xx,v1)-SG_Product(tau_xy,v2)-SG_Product(lambda,gradT1)
                                                              ! F_euler-(tau_xx*u+tau_xy*v+tau_xz*w-q_x) q_x=-lambda*T_x
  ! viscous fluxes in y-direction
  g(SGI1,i,j)=0.
  g(SGI2,i,j)=-tau_xy                                         ! F_euler-mu*(u_y+v_x)
  g(SGI3,i,j)=-tau_yy                                         ! F_euler-4/3*mu*v_y+2/3*mu*(u_x+w_z)
  g(SGI4,i,j)=0.
  g(SGI5,i,j)=-SG_Product(tau_xy,v1)-SG_Product(tau_yy,v2)-SG_Product(lambda,gradT2)
                                                              ! F_euler-(tau_yx*u+tau_yy*v+tau_yz*w-q_y) q_y=-lambda*T_y
  ! viscous fluxes in z-direction
  h(:,i,j)=0.
#endif


END DO ; END DO !i,j

#ifdef DEBUG
! ===============================================================================
! Following dummy calls do suppress compiler warnings of unused Riemann-functions
! ===============================================================================
IF (0.EQ.1) THEN
  WRITE (*,*) gradUz_Face
END IF
#endif
END SUBROUTINE EvalDiffFlux2D

#endif /*PARABOLIC*/

END MODULE MOD_Flux
