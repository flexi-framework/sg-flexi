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
!> Subroutines  that provide gas properties and conversion between primitive and conservative variables for a ideal gas.
!==================================================================================================================================
MODULE MOD_EOS
! MODULES
IMPLICIT NONE
PRIVATE
!----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! Private Part ---------------------------------------------------------------------------------------------------------------------
! Public Part ----------------------------------------------------------------------------------------------------------------------

INTERFACE InitEOS
  MODULE PROCEDURE InitEOS
END INTERFACE

INTERFACE ConsToPrim
  MODULE PROCEDURE ConsToPrim
END INTERFACE

INTERFACE PrimToCons
  MODULE PROCEDURE PrimToCons
END INTERFACE

INTERFACE ConsToPrimDet
  MODULE PROCEDURE ConsToPrimDet
  MODULE PROCEDURE ConsToPrimDet_Side
  MODULE PROCEDURE ConsToPrimDet_Volume
END INTERFACE

INTERFACE ConsToPrimSG
  MODULE PROCEDURE ConsToPrimSG
  MODULE PROCEDURE ConsToPrimSG_Side
  MODULE PROCEDURE ConsToPrimSG_Volume
END INTERFACE

INTERFACE PrimToConsDet
  MODULE PROCEDURE PrimToConsDet
  MODULE PROCEDURE PrimToConsDet_Side
  MODULE PROCEDURE PrimToConsDet_Volume
END INTERFACE

INTERFACE PrimToConsSG
  MODULE PROCEDURE PrimToConsSG
  MODULE PROCEDURE PrimToConsSG_Side
  MODULE PROCEDURE PrimToConsSG_Volume
END INTERFACE

INTERFACE PRESSURE_RIEMANN
  MODULE PROCEDURE PRESSURE_RIEMANN
END INTERFACE

PUBLIC::InitEos
PUBLIC::ConsToPrim
PUBLIC::PrimToCons
PUBLIC::ConsToPrimDet
PUBLIC::ConsToPrimSG
PUBLIC::PrimToConsDet
PUBLIC::PrimToConsSG
PUBLIC::PRESSURE_RIEMANN
PUBLIC::DefineParametersEos
!==================================================================================================================================

CONTAINS

!==================================================================================================================================
!> Define parameters for the used eos
!==================================================================================================================================
SUBROUTINE DefineParametersEos()
! MODULES
USE MOD_Globals
USE MOD_ReadInTools ,ONLY: prms
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
!==================================================================================================================================
CALL prms%SetSection("Equation of State")
CALL prms%CreateLogicalOption('UseNonDimensionalEqn',"Set true to compute R and mu from bulk Mach Reynolds (nondimensional form.",&
                                                '.FALSE.')
CALL prms%CreateRealOption(     'BulkMach',     "Bulk Mach     (UseNonDimensionEqn=T)")
CALL prms%CreateRealOption(     'BulkReynolds', "Bulk Reynolds (UseNonDimensionEqn=T)")
CALL prms%CreateRealOption(     'kappa',        "Heat capacity ratio / isentropic exponent", '1.4')
CALL prms%CreateRealOption(     'R',            "Specific gas constant", '287.058')
CALL prms%CreateRealOption(     'Pr',           "Prandtl number", '0.72')
CALL prms%CreateRealOption(     'mu0',          "Dynamic Viscosity", '0.')
CALL prms%CreateRealOption(     'Ts',           "Sutherland's law for variable viscosity: Ts", '110.4')
CALL prms%CreateRealOption(     'Tref',         "Sutherland's law for variable viscosity: Tref ", '280.0')
CALL prms%CreateRealOption(     'ExpoSuth',     "Sutherland's law for variable viscosity: Exponent", '1.5')

END SUBROUTINE DefineParametersEos

!==================================================================================================================================
!> Initialize variables needed by the ideal gas equation of state.
!==================================================================================================================================
SUBROUTINE InitEos()
! MODULES
USE MOD_PreProc
USE MOD_Globals
USE MOD_ReadInTools
USE MOD_EOS_Vars      ,ONLY: Kappa,KappaM1,KappaP1,cp,cv
USE MOD_EOS_Vars      ,ONLY: R,sKappaM1,sKappaP1
#if PARABOLIC
USE MOD_EOS_Vars      ,ONLY: mu0,Pr,KappaSpr
#if PP_VISC == 1
USE MOD_EOS_Vars      ,ONLY: Ts,cSuth
#endif
#if (PP_VISC == 1) || (PP_VISC == 2)
USE MOD_EOS_Vars      ,ONLY: Tref,ExpoSuth
#endif
#endif /*PARABOLIC*/

! IMPLICIT VARIABLE HANDLING
 IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL    :: BulkMach,BulkReynolds
LOGICAL :: UseNonDimensionalEqn=.FALSE.
!==================================================================================================================================
SWRITE(UNIT_StdOut,'(132("-"))')
SWRITE(UNIT_stdOut,'(A)') ' INIT IDEAL GAS...'


UseNonDimensionalEqn=GETLOGICAL('UseNonDimensionalEqn','.FALSE.')
IF(UseNonDimensionalEqn)THEN
  BulkMach    =GETREAL('BulkMach')
  BulkReynolds=GETREAL('BulkReynolds')
END IF

! Gas constants
Kappa    =GETREAL('kappa','1.4')
KappaM1  =Kappa-1.
sKappaM1 =1./KappaM1
KappaP1  =Kappa+1.
sKappaP1 =1./(KappaP1)
IF(.NOT. UseNonDimensionalEqn)THEN
  R=GETREAL('R','287.058')
ELSE
  R=1./(Kappa*BulkMach*BulkMach)
  SWRITE(UNIT_stdOut,'(A,ES16.7)')' |                              R | Set to 1/(Kappa*BulkMach**2)=',R
END IF

cp=R*kappa/(kappa-1.)
cv=R/(kappa-1.)

#if PARABOLIC
Pr       =GETREAL('Pr','0.72')
KappasPr =Kappa/Pr

! Viscosity
#if   PP_VISC == 0
! Constant viscosity
IF(.NOT. UseNonDimensionalEqn)THEN
  mu0=GETREAL('mu0','0.0')
ELSE
  mu0=1./BulkReynolds
  SWRITE(UNIT_stdOut,'(A,ES16.7)')' |                            mu0 | Set to 1/BulkReynolds=',mu0
END IF
#elif PP_VISC == 1
! mu-Sutherland
mu0     =GETREAL('mu0','1.735E-5')
Ts      =GETREAL('Ts','110.4')
Tref    =1./GETREAL('Tref','280.')
ExpoSuth=GETREAL('ExpoSuth','1.5')
Ts      =Ts*Tref
cSuth   =Ts**ExpoSuth*(1+Ts)/(2*Ts*Ts)
#elif PP_VISC == 2
! mu power-law
IF(.NOT. UseNonDimensionalEqn)THEN
  mu0=GETREAL('mu0','0.')
  Tref    =GETREAL('Tref')
ELSE
  mu0=1./BulkReynolds
  SWRITE(UNIT_stdOut,'(A,ES16.7)')' |                            mu0 | Set to 1/BulkReynolds=',mu0
  Tref=1.
  SWRITE(UNIT_stdOut,*)'|                           Tref | Set to 1.'
END IF
ExpoSuth=GETREAL('ExpoSuth')
mu0     =mu0/Tref**ExpoSuth
#endif
#endif /*PARABOLIC*/

SWRITE(UNIT_stdOut,'(A)')' INIT IDEAL-GAS DONE!'
SWRITE(UNIT_StdOut,'(132("-"))')
END SUBROUTINE InitEos

!==================================================================================================================================
!> Transformation from conservative variables to primitive variables for a single state
!==================================================================================================================================
SUBROUTINE ConsToPrim(prim,cons)
! MODULES
USE MOD_Globals ,ONLY: Abort
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
REAL,INTENT(IN)  :: cons(PP_nVarDet)     !< vector of conservative variables
REAL,INTENT(OUT) :: prim(PP_nVarPrimDet) !< vector of primitive variables (density,velocities,temperature,pressure)
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!==================================================================================================================================
CALL abort(__STAMP__,'Cons to Prim must be changed')
#ifdef DEBUG
! ===============================================================================
! Following dummy calls do suppress compiler warnings of unused Riemann-functions
! ===============================================================================
IF (0.EQ.1) THEN
  WRITE (*,*) cons, prim
END IF
#endif /*DEBUG*/
END SUBROUTINE ConsToPrim

!==================================================================================================================================
!> Transformation from conservative variables to primitive variables for a single state
!==================================================================================================================================
SUBROUTINE PrimToCons(prim,cons)
! MODULES
USE MOD_Globals ,ONLY: Abort
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
REAL,INTENT(IN)  :: cons(PP_nVarDet)     !< vector of conservative variables
REAL,INTENT(OUT) :: prim(PP_nVarPrimDet) !< vector of primitive variables (density,velocities,temperature,pressure)
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!==================================================================================================================================
CALL abort(__STAMP__,'Cons to Prim must be changed')
#ifdef DEBUG
! ===============================================================================
! Following dummy calls do suppress compiler warnings of unused Riemann-functions
! ===============================================================================
IF (0.EQ.1) THEN
  WRITE (*,*) cons, prim
END IF
#endif /*DEBUG*/
END SUBROUTINE PrimToCons

!==================================================================================================================================
!> Transformation from conservative variables to primitive variables for a single state
!==================================================================================================================================
PURE SUBROUTINE ConsToPrimDet(prim,cons)
! MODULES
USE MOD_EOS_Vars ,ONLY: KappaM1,R
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
REAL,INTENT(IN)  :: cons(PP_nVarDet)     !< vector of conservative variables
REAL,INTENT(OUT) :: prim(PP_nVarPrimDet) !< vector of primitive variables (density,velocities,temperature,pressure)
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL             :: sRho    ! 1/Rho
!==================================================================================================================================
sRho=1./cons(1)
! density
prim(1)=cons(1)
! velocity
prim(2:3)=cons(2:3)*sRho
#if (PP_dim==3)
prim(4)=cons(4)*sRho
#else
prim(4)=0.
#endif
! pressure
prim(5)=KappaM1*(cons(5)-0.5*SUM(cons(2:4)*prim(2:4)))
! temperature
prim(6) = prim(5)*sRho / R
END SUBROUTINE ConsToPrimDet

!==================================================================================================================================
!> Transformation from conservative variables to primitive variables on a single side
!==================================================================================================================================
PURE SUBROUTINE ConsToPrimDet_Side(Nloc,prim,cons)
! MODULES
USE MOD_Eos_Vars ,ONLY: KappaM1,R
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER,INTENT(IN) :: Nloc
REAL,INTENT(IN)    :: cons(PP_nVarDet    ,0:Nloc,0:ZDIM(Nloc)) !< vector of conservative variables
REAL,INTENT(OUT)   :: prim(PP_nVarPrimDet,0:Nloc,0:ZDIM(Nloc)) !< vector of primitive variables 
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL             :: sRho    ! 1/Rho
INTEGER          :: p,q
!==================================================================================================================================
DO q=0,ZDIM(Nloc); DO p=0,Nloc
  sRho=1./cons(1,p,q)
  ! density
  prim(1,p,q)=cons(1,p,q)
  ! velocity
  prim(2:3,p,q)=cons(2:3,p,q)*sRho
#if (PP_dim==3)
  prim(4,p,q)=cons(4,p,q)*sRho
#else
  prim(4,p,q)=0.
#endif
  ! pressure
  prim(5,p,q)=KappaM1*(cons(5,p,q)-0.5*SUM(cons(2:4,p,q)*prim(2:4,p,q)))
  ! temperature
  prim(6,p,q) = prim(5,p,q)*sRho / R
END DO; END DO
END SUBROUTINE ConsToPrimDet_Side

!==================================================================================================================================
!> Transformation from conservative variables to primitive variables in the whole volume
!==================================================================================================================================
PURE SUBROUTINE ConsToPrimDet_Volume(Nloc,prim,cons)
! MODULES
USE MOD_Eos_Vars  ,ONLY: KappaM1,R
USE MOD_Mesh_Vars ,ONLY: nElems
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER,INTENT(IN) :: Nloc
REAL,INTENT(IN)    :: cons(PP_nVarDet    ,0:Nloc,0:Nloc,0:ZDIM(Nloc),1:nElems) !< vector of conservative variables
REAL,INTENT(OUT)   :: prim(PP_nVarPrimDet,0:Nloc,0:Nloc,0:ZDIM(Nloc),1:nElems) !< vector of primitive variables 
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL             :: sRho    ! 1/Rho
INTEGER          :: i,j,k,iElem
!==================================================================================================================================
DO iElem=1,nElems
  DO k=0,ZDIM(Nloc); DO j=0,Nloc; DO i=0,Nloc
    sRho=1./cons(1,i,j,k,iElem)
    ! density
    prim(1,i,j,k,iElem)=cons(1,i,j,k,iElem)
    ! velocity
    prim(2:3,i,j,k,iElem)=cons(2:3,i,j,k,iElem)*sRho
#if (PP_dim==3)
    prim(4,i,j,k,iElem)=cons(4,i,j,k,iElem)*sRho
#else
    prim(4,i,j,k,iElem)=0.
#endif
    ! pressure
    prim(5,i,j,k,iElem)=KappaM1*(cons(5,i,j,k,iElem)-0.5*SUM(cons(2:4,i,j,k,iElem)*prim(2:4,i,j,k,iElem)))
    ! temperature
    prim(6,i,j,k,iElem) = prim(5,i,j,k,iElem)*sRho / R
  END DO; END DO; END DO! i,j,k=0,Nloc
END DO ! iElem
END SUBROUTINE ConsToPrimDet_Volume

!==================================================================================================================================
!> Transformation from conservative variables to primitive variables for a single state
!==================================================================================================================================
SUBROUTINE ConsToPrimSG(prim,cons)
!PURE SUBROUTINE ConsToPrimSG(prim,cons)
! MODULES
USE MOD_PreProc
USE MOD_EOS_Vars     ,ONLY: KappaM1,R
USE MOD_SG_Operators ,ONLY: SG_Inv,SG_DOT_PRODUCT,SG_Product
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
REAL,INTENT(IN)  :: cons(PP_nVar)     !< vector of conservative variables
REAL,INTENT(OUT) :: prim(PP_nVarPrim) !< vector of primitive variables (density,velocities,temperature,pressure)
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL             :: sRho(1:PP_nCoef)    ! 1/Rho
!==================================================================================================================================
sRho(:)=SG_Inv(cons(DENS))
! density
prim(PDENS)=cons(DENS)
! velocity
prim(PVEL1)=SG_Product(cons(MOM1),sRho)
prim(PVEL2)=SG_Product(cons(MOM2),sRho)
#if (PP_dim==3)
prim(PVEL3)=SG_Product(cons(MOM3),sRho)
#else
prim(PVEL3)=0.
#endif
! pressure
prim(PPRES)=KappaM1*(cons(ENER)-0.5*SG_DOT_PRODUCT(cons(MOMV),prim(PVELV)))
! temperature
prim(PTEMP) = SG_Product(Prim(PPRES),sRho) / R
END SUBROUTINE ConsToPrimSG

!==================================================================================================================================
!> Transformation from conservative variables to primitive variables on a single side
!==================================================================================================================================
SUBROUTINE ConsToPrimSG_Side(Nloc,prim,cons)
! MODULES
USE MOD_PreProc
USE MOD_Eos_Vars     ,ONLY: KappaM1,R
USE MOD_SG_Operators ,ONLY: SG_Inv,SG_DOT_PRODUCT,SG_Product
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER,INTENT(IN) :: Nloc
REAL,INTENT(IN)    :: cons(PP_nVar    ,0:Nloc,0:ZDIM(Nloc)) !< vector of conservative variables
REAL,INTENT(OUT)   :: prim(PP_nVarPrim,0:Nloc,0:ZDIM(Nloc)) !< vector of primitive variables 
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL             :: sRho(1:PP_nCoef)    ! 1/Rho
INTEGER          :: p,q
!==================================================================================================================================
DO q=0,ZDIM(Nloc); DO p=0,Nloc
  sRho(:)=SG_Inv(cons(DENS,p,q))
  ! density
  prim(PDENS,p,q)=cons(DENS,p,q)
  ! velocity
  prim(PVEL1,p,q)=SG_Product(cons(MOM1,p,q),sRho)
  prim(PVEL2,p,q)=SG_Product(cons(MOM2,p,q),sRho)
#if (PP_dim==3)
  prim(PVEL3,p,q)=SG_Product(cons(MOM3,p,q),sRho)
#else
  prim(PVEL3,p,q)=0.
#endif
  ! pressure
  prim(PPRES,p,q)=KappaM1*(cons(ENER,p,q)-0.5*SG_DOT_PRODUCT(cons(MOMV,p,q),prim(PVELV,p,q)))
  ! temperature
  prim(PTEMP,p,q) = SG_Product(Prim(PPRES,p,q),sRho) / R
END DO; END DO
END SUBROUTINE ConsToPrimSG_Side

!==================================================================================================================================
!> Transformation from conservative variables to primitive variables in the whole volume
!==================================================================================================================================
SUBROUTINE ConsToPrimSG_Volume(Nloc,prim,cons)
! MODULES
USE MOD_PreProc
USE MOD_Eos_Vars     ,ONLY: KappaM1,R
USE MOD_Mesh_Vars    ,ONLY: nElems
USE MOD_SG_Operators ,ONLY: SG_Inv,SG_DOT_PRODUCT,SG_Product
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER,INTENT(IN) :: Nloc
REAL,INTENT(IN)    :: cons(PP_nVar    ,0:Nloc,0:Nloc,0:ZDIM(Nloc),1:nElems) !< vector of conservative variables
REAL,INTENT(OUT)   :: prim(PP_nVarPrim,0:Nloc,0:Nloc,0:ZDIM(Nloc),1:nElems) !< vector of primitive variables 
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL             :: sRho(1:PP_nCoef)    ! 1/Rho
INTEGER          :: i,j,k,iElem
!==================================================================================================================================
DO iElem=1,nElems
  DO k=0,ZDIM(Nloc); DO j=0,Nloc; DO i=0,Nloc
    sRho(:)=SG_Inv(cons(DENS,i,j,k,iElem))
    ! density
    prim(PDENS,i,j,k,iElem)=cons(DENS,i,j,k,iElem)
    ! velocity
    prim(PVEL1,i,j,k,iElem)=SG_Product(cons(MOM1,i,j,k,iElem),sRho)
    prim(PVEL2,i,j,k,iElem)=SG_Product(cons(MOM2,i,j,k,iElem),sRho)
#if (PP_dim==3)
    prim(PVEL3,i,j,k,iElem)=SG_Product(cons(MOM3,i,j,k,iElem),sRho)
#else
    prim(PVEL3,i,j,k,iElem)=0.
#endif
    ! pressure
    prim(PPRES,i,j,k,iElem)=KappaM1*(cons(ENER,i,j,k,iElem)-0.5*SG_DOT_PRODUCT(cons(MOMV,i,j,k,iElem),prim(PVELV,i,j,k,iElem)))
    ! temperature
    prim(PTEMP,i,j,k,iElem) = SG_Product(Prim(PPRES,i,j,k,iElem),sRho) / R
  END DO; END DO; END DO! i,j,k=0,Nloc
END DO ! iElem
END SUBROUTINE ConsToPrimSG_Volume

!==================================================================================================================================
!> Transformation from primitive to conservative variables
!==================================================================================================================================
PURE SUBROUTINE PrimToConsDet(prim,cons)
! MODULES
USE MOD_Eos_Vars ,ONLY: sKappaM1
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
REAL,INTENT(IN)  :: prim(PP_nVarPrimDet) !< vector of primitive variables
REAL,INTENT(OUT) :: cons(PP_nVarDet)     !< vector of conservative variables
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!==================================================================================================================================
! density
cons(1)=prim(1)
! momentum
cons(2:3)=prim(2:3)*prim(1)
#if (PP_dim==3)
cons(4)=prim(4)*prim(1)
#else
cons(4)=0.
#endif
! energy
cons(5)=sKappaM1*prim(5)+0.5*SUM(cons(2:4)*prim(2:4))
END SUBROUTINE PrimToConsDet

!==================================================================================================================================
!> Transformation from primitive to conservative variables on a single side
!==================================================================================================================================
PURE SUBROUTINE PrimToConsDet_Side(Nloc,prim,cons)
! MODULES
USE MOD_Eos_Vars ,ONLY: sKappaM1
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER,INTENT(IN):: Nloc
REAL,INTENT(IN)   :: prim(PP_nVarPrimDet,0:Nloc,0:ZDIM(Nloc)) !< vector of primitive variables
REAL,INTENT(OUT)  :: cons(PP_nVarDet    ,0:Nloc,0:ZDIM(Nloc))     !< vector of conservative variables
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER           :: p,q
!==================================================================================================================================
DO q=0,ZDIM(Nloc); DO p=0,Nloc
  ! density
  cons(1,p,q)=prim(1,p,q)
  ! momentum
  cons(2:3,p,q)=prim(2:3,p,q)*prim(1,p,q)
#if (PP_dim==3)
  cons(4,p,q)=prim(4,p,q)*prim(1,p,q)
#else
  cons(4,p,q)=0.
#endif
  ! energy
  cons(5,p,q)=sKappaM1*prim(5,p,q)+0.5*SUM(cons(2:4,p,q)*prim(2:4,p,q))
END DO; END DO ! p,q=0,Nloc
END SUBROUTINE PrimToConsDet_Side

!==================================================================================================================================
!> Transformation from primitive to conservative variables in the whole volume
!==================================================================================================================================
PURE SUBROUTINE PrimToConsDet_Volume(Nloc,prim,cons)
! MODULES
USE MOD_Eos_Vars  ,ONLY: sKappaM1
USE MOD_Mesh_Vars ,ONLY: nElems
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER,INTENT(IN):: Nloc
REAL,INTENT(IN)   :: prim(PP_nVarPrimDet,0:Nloc,0:Nloc,0:ZDIM(Nloc),1:nElems)     !< vector of primitive variables
REAL,INTENT(OUT)  :: cons(PP_nVarDet    ,0:Nloc,0:Nloc,0:ZDIM(Nloc),1:nElems)     !< vector of conservative variables
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER           :: p,q,r,iElem
!==================================================================================================================================
DO iElem=1,nElems
  DO r=0,ZDIM(Nloc); DO q=0,Nloc; DO p=0,Nloc
    ! density
    cons(1,p,q,r,iElem)=prim(1,p,q,r,iElem)
    ! momentum
    cons(2:3,p,q,r,iElem)=prim(2:3,p,q,r,iElem)*prim(1,p,q,r,iElem)
#if (PP_dim==3)
    cons(4,p,q,r,iElem)=prim(4,p,q,r,iElem)*prim(1,p,q,r,iElem)
#else
    cons(4,p,q,r,iElem)=0.
#endif
    ! energy
    cons(5,p,q,r,iElem)=sKappaM1*prim(5,p,q,r,iElem)+0.5*SUM(cons(2:4,p,q,r,iElem)*prim(2:4,p,q,r,iElem))
  END DO; END DO; END DO ! p,q=0,Nloc
END DO
END SUBROUTINE PrimToConsDet_Volume

!==================================================================================================================================
!> Transformation from primitive to conservative variables
!==================================================================================================================================
SUBROUTINE PrimToConsSG(prim,cons)
! MODULES
USE MOD_PreProc
USE MOD_Eos_Vars     ,ONLY: sKappaM1
USE MOD_SG_Operators ,ONLY: SG_Inv,SG_DOT_PRODUCT,SG_Product
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
REAL,INTENT(IN)  :: prim(PP_nVarPrim) !< vector of primitive variables
REAL,INTENT(OUT) :: cons(PP_nVar)     !< vector of conservative variables
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!==================================================================================================================================
! density
cons(DENS)=prim(PDENS)
! momentum
cons(MOM1)=SG_Product(prim(PVEL1),prim(PDENS))
cons(MOM2)=SG_Product(prim(PVEL2),prim(PDENS))
#if (PP_dim==3)
cons(MOM3)=SG_Product(prim(PVEL3),prim(PDENS))
#else
cons(MOM3)=0.
#endif
! energy
cons(ENER)=sKappaM1*prim(PPRES)+0.5*SG_DOT_PRODUCT(cons(MOMV),prim(PVELV))
END SUBROUTINE PrimToConsSG



!==================================================================================================================================
!> Transformation from primitive to conservative variables on a single side
!==================================================================================================================================
SUBROUTINE PrimToConsSG_Side(Nloc,prim,cons)
! MODULES
USE MOD_PreProc
USE MOD_Eos_Vars     ,ONLY: sKappaM1
USE MOD_SG_Operators ,ONLY: SG_Inv,SG_DOT_PRODUCT,SG_Product
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER,INTENT(IN):: Nloc
REAL,INTENT(IN)   :: prim(PP_nVarPrim,0:Nloc,0:ZDIM(Nloc)) !< vector of primitive variables
REAL,INTENT(OUT)  :: cons(PP_nVar    ,0:Nloc,0:ZDIM(Nloc))     !< vector of conservative variables
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER           :: p,q
!==================================================================================================================================
DO q=0,ZDIM(Nloc); DO p=0,Nloc
  ! density
  cons(DENS,p,q)=prim(PDENS,p,q)
  ! momentum
  cons(MOM1,p,q)=SG_Product(prim(PVEL1,p,q),prim(PDENS,p,q))
  cons(MOM2,p,q)=SG_Product(prim(PVEL2,p,q),prim(PDENS,p,q))
#if (PP_dim==3)
  cons(MOM3,p,q)=SG_Product(prim(PVEL3,p,q),prim(PDENS,p,q))
#else
  cons(MOM3,p,q)=0.
#endif
  ! energy
  cons(ENER,p,q)=sKappaM1*prim(PPRES,p,q)+0.5*SG_DOT_PRODUCT(cons(MOMV,p,q),prim(PVELV,p,q))
END DO; END DO ! p,q=0,Nloc
END SUBROUTINE PrimToConsSG_Side



!==================================================================================================================================
!> Transformation from primitive to conservative variables in the whole volume
!==================================================================================================================================
SUBROUTINE PrimToConsSG_Volume(Nloc,prim,cons)
! MODULES
USE MOD_PreProc
USE MOD_Eos_Vars     ,ONLY: sKappaM1
USE MOD_Mesh_Vars    ,ONLY: nElems
USE MOD_SG_Operators ,ONLY: SG_Inv,SG_DOT_PRODUCT,SG_Product
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER,INTENT(IN):: Nloc
REAL,INTENT(IN)   :: prim(PP_nVarPrim,0:Nloc,0:Nloc,0:ZDIM(Nloc),1:nElems)     !< vector of primitive variables
REAL,INTENT(OUT)  :: cons(PP_nVar    ,0:Nloc,0:Nloc,0:ZDIM(Nloc),1:nElems)     !< vector of conservative variables
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER           :: p,q,r,iElem
!==================================================================================================================================
DO iElem=1,nElems
  DO r=0,ZDIM(Nloc); DO q=0,Nloc; DO p=0,Nloc
    ! density
    cons(DENS,p,q,r,iElem)=prim(PDENS,p,q,r,iElem)
    ! momentum
    cons(MOM1,p,q,r,iElem)=SG_Product(prim(PVEL1,p,q,r,iElem),prim(PDENS,p,q,r,iElem))
    cons(MOM2,p,q,r,iElem)=SG_Product(prim(PVEL2,p,q,r,iElem),prim(PDENS,p,q,r,iElem))
#if (PP_dim==3)
    cons(MOM3,p,q,r,iElem)=SG_Product(prim(PVEL3,p,q,r,iElem),prim(PDENS,p,q,r,iElem))
#else
    cons(MOM3,p,q,r,iElem)=0.
#endif
    ! energy
    cons(ENER,p,q,r,iElem)=sKappaM1*prim(PPRES,p,q,r,iElem)+0.5*SG_DOT_PRODUCT(cons(MOMV,p,q,r,iElem),prim(PVELV,p,q,r,iElem))
  END DO; END DO; END DO ! p,q=0,Nloc
END DO
END SUBROUTINE PrimToConsSG_Volume

!==================================================================================================================================
!> Riemann solver function to get pressure at BCs
!==================================================================================================================================
PURE FUNCTION PRESSURE_RIEMANN(U_Prim)
!==================================================================================================================================
! MODULES
USE MOD_Eos_Vars ,ONLY: Kappa,KappaM1,sKappaM1,sKappaP1
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE 
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
REAL,INTENT(IN) :: U_Prim(PP_nVarPrimDet)
REAL            :: PRESSURE_RIEMANN
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
REAL     :: kappaFac,ar,br,P_RP
!==================================================================================================================================
kappaFac=2.*Kappa*sKappaM1
IF(U_Prim(2) .LE. 0.)THEN ! rarefaction
  P_RP=U_Prim(5) * MAX(0.0001,(1.+0.5*KappaM1*U_Prim(2)/SQRT(Kappa*U_Prim(5)/U_Prim(1))))**kappaFac
ELSE ! shock
  ar=2.*sKappaP1/U_Prim(1)
  br=KappaM1*sKappaP1*U_Prim(5)
  P_RP=U_Prim(5)+U_Prim(2)/ar*0.5*(U_Prim(2)+SQRT(U_Prim(2)*U_Prim(2)+4.*ar*(U_Prim(5)+br)))
END IF
PRESSURE_RIEMANN=P_RP
END FUNCTION PRESSURE_RIEMANN


END MODULE MOD_EOS
