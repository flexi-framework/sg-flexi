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

!==================================================================================================================================
!> Contains the routines for computing the fluxes for the Split-DG algorithm
!==================================================================================================================================
#include "flexi.h"
#include "eos.h"

MODULE MOD_SplitFlux
! MODULES
IMPLICIT NONE
PRIVATE
!----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
ABSTRACT INTERFACE
  SUBROUTINE VolumeFlux(URef,UPrimRef,U,UPrim,MRef,Met,Flux)
    USE MOD_PreProc
    REAL,DIMENSION(PP_nVar    ),INTENT(IN)  :: URef,U
    REAL,DIMENSION(PP_nVarPrim),INTENT(IN)  :: UPrimRef,UPrim
    REAL,DIMENSION(1:3        ),INTENT(IN)  :: MRef,Met
    REAL,DIMENSION(PP_nVar    ),INTENT(OUT) :: Flux
  END SUBROUTINE
END INTERFACE

ABSTRACT INTERFACE
  SUBROUTINE SurfaceFlux(U_LL,U_RR,F)
    USE MOD_PreProc
    REAL,DIMENSION(PP_2Var),INTENT(IN)  :: U_LL,U_RR
    REAL,DIMENSION(PP_nVar),INTENT(OUT) :: F
  END SUBROUTINE
END INTERFACE

PROCEDURE(VolumeFlux),POINTER    :: SplitDGVolume_pointer    !< pointer defining the SpliDG formulation beeing used
PROCEDURE(SurfaceFlux),POINTER   :: SplitDGSurface_pointer   !< pointer defining the SpliDG formulation beeing used

INTEGER,PARAMETER      :: PRM_SPLITDG_SD          = 0
INTEGER,PARAMETER      :: PRM_SPLITDG_MO          = 1
INTEGER,PARAMETER      :: PRM_SPLITDG_DU          = 2
INTEGER,PARAMETER      :: PRM_SPLITDG_KG          = 3
INTEGER,PARAMETER      :: PRM_SPLITDG_PI          = 4

INTERFACE InitSplitDG
  MODULE PROCEDURE InitSplitDG
END INTERFACE

PUBLIC::InitSplitDG,DefineParametersSplitDG
PUBLIC::SplitDGSurface_pointer,SplitDGVolume_pointer
!==================================================================================================================================

CONTAINS

!==================================================================================================================================
!> Define parameters
!==================================================================================================================================
SUBROUTINE DefineParametersSplitDG()
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
CALL prms%SetSection("SplitDG")
CALL prms%CreateIntFromStringOption('SplitDG',"SplitDG formulation to be used: SD, MO, DU, KG, PI")
CALL addStrListEntry('SplitDG','sd',           PRM_SPLITDG_SD)
CALL addStrListEntry('SplitDG','mo',           PRM_SPLITDG_MO)
CALL addStrListEntry('SplitDG','du',           PRM_SPLITDG_DU)
CALL addStrListEntry('SplitDG','kg',           PRM_SPLITDG_KG)
CALL addStrListEntry('SplitDG','pi',           PRM_SPLITDG_PI)

END SUBROUTINE DefineParametersSplitDG

!==================================================================================================================================!
!> Initialize function pointers for the specific split version in use 
!==================================================================================================================================!
SUBROUTINE InitSplitDG()
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_ReadInTools ,ONLY: GETINTFROMSTR
!----------------------------------------------------------------------------------------------------------------------------------
IMPLICIT NONE
! INPUT / OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                     :: SplitDG
#ifdef DEBUG
REAL,DIMENSION(PP_nVar    ) :: U         ! dummy variables, only to suppress compiler warnings
REAL,DIMENSION(PP_nVarPrim) :: UPrim     ! dummy variables, only to suppress compiler warnings
REAL,DIMENSION(PP_nVar    ) :: f,g,h     ! dummy variables, only to suppress compiler warnings
REAL,DIMENSION(PP_2Var    ) :: U_LL,U_RR ! dummy variables, only to suppress compiler warnings
#endif
!==================================================================================================================================
! check if Gauss-Lobatto-Pointset is beeing used
#if (PP_NodeType==1)
CALL CollectiveStop(__STAMP__,&
  'Wrong Pointset: Gauss-Lobatto-Points are mandatory for using SplitDG !')
#endif
! set pointers
SplitDG = GETINTFROMSTR('SplitDG')
SELECT CASE(SplitDG)
#ifdef SGdummy
CASE(PRM_SPLITDG_SD)
  SplitDGVolume_pointer  => SplitVolumeFluxSD
  SplitDGSurface_pointer => SplitSurfaceFluxSD
CASE(PRM_SPLITDG_MO)
  SplitDGVolume_pointer  => SplitVolumeFluxMO
  SplitDGSurface_pointer => SplitSurfaceFluxMO
CASE(PRM_SPLITDG_DU)
  SplitDGVolume_pointer  => SplitVolumeFluxDU
  SplitDGSurface_pointer => SplitSurfaceFluxDU
CASE(PRM_SPLITDG_KG)
  SplitDGVolume_pointer  => SplitVolumeFluxKG
  SplitDGSurface_pointer => SplitSurfaceFluxKG
#endif /*SGdummy*/
CASE(PRM_SPLITDG_PI)
  SplitDGVolume_pointer  => SplitVolumeFluxPI
  SplitDGSurface_pointer => SplitSurfaceFluxPI
CASE DEFAULT
  CALL CollectiveStop(__STAMP__,&
    'SplitDG formulation not defined!')
END SELECT

#ifdef DEBUG
! ===============================================================================
! Following dummy calls do suppress compiler warnings of unused Riemann-functions
! ===============================================================================
IF (0.EQ.1) THEN
  U=1. ;  UPrim=1. ;   U_LL=1. ;   U_RR=1.
  CALL SplitDGVolume_pointer  (U,UPrim,U,UPrim,f,g,h)
  CALL SplitDGSurface_pointer (U_LL,U_RR,F)
#ifdef SGdummy
  CALL SplitVolumeFluxSD     (U,UPrim,U,UPrim,f,g,h)
  CALL SplitSurfaceFluxSD    (U_LL,U_RR,F)
  CALL SplitVolumeFluxMO     (U,UPrim,U,UPrim,f,g,h)
  CALL SplitSurfaceFluxMO    (U_LL,U_RR,F)
  CALL SplitVolumeFluxDU     (U,UPrim,U,UPrim,f,g,h)
  CALL SplitSurfaceFluxDU    (U_LL,U_RR,F)
  CALL SplitVolumeFluxKG     (U,UPrim,U,UPrim,f,g,h)
  CALL SplitSurfaceFluxKG    (U_LL,U_RR,F)
#endif /*SGdummy*/
  CALL SplitVolumeFluxPI     (U,UPrim,U,UPrim,f,g,h)
  CALL SplitSurfaceFluxPI    (U_LL,U_RR,F)
END IF
#endif
END SUBROUTINE InitSplitDG

#ifdef SGdummy
!==================================================================================================================================
!> Computes the Split-Flux retaining the standart NS-Equations
!> Attention 1: Factor 2 from differentiation matrix is already been considered
!==================================================================================================================================
SUBROUTINE SplitVolumeFluxSD(URef,UPrimRef,U,UPrim,MRef,Met,Flux)
! MODULES
USE MOD_PreProc
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
REAL,DIMENSION(PP_nVar    ),INTENT(IN)  :: URef,U ! conserved variables
REAL,DIMENSION(PP_nVarPrim),INTENT(IN)  :: UPrimRef,UPrim ! primitive variables
REAL,DIMENSION(1:3        ),INTENT(IN)  :: MRef,Met ! metric terms
REAL,DIMENSION(PP_nVar    ),INTENT(OUT) :: Flux ! flux in reverence space
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                                    :: EpRef,Ep ! auxilery variable for (rho*e+p)
REAL,DIMENSION(PP_nVar    )             :: fTilde,gTilde ! flux in physical space
#if PP_dim == 3
REAL,DIMENSION(PP_nVar    )             :: hTilde        ! flux in physical space
#endif
!==================================================================================================================================
! compute auxilery variables
EpRef = URef(SGI5) + UPrimRef(SGI5)
Ep    = U(SGI5) + UPrim(SGI5)

! local Euler fluxes x-direction
  fTilde(SGI1) = (URef(SGI2) + U(SGI2))                                           ! {rho*u}
  fTilde(SGI2) = (URef(SGI2)*UPrimRef(SGI2)+UPrimRef(SGI5) + U(SGI2)*UPrim(SGI2)+UPrim(SGI5)) ! {rho*u²}+{p}
  fTilde(SGI3) = (URef(SGI2)*UPrimRef(SGI3) + U(SGI2)*UPrim(SGI3))                      ! {rho*u*v}
#if PP_dim == 3
  fTilde(SGI4) = (URef(SGI2)*UPrimRef(SGI4) + U(SGI2)*UPrim(SGI4))                      ! {rho*u*w}
#else
  fTilde(SGI4) = 0.
#endif
  fTilde(SGI5) = (EpRef*UPrimRef(SGI2) + Ep*UPrim(SGI2))                          ! {(rho*e+p)*u}
! local Euler fluxes y-direction
  gTilde(SGI1) = (URef(SGI3) + U(SGI3))                                           ! {rho*v}
  gTilde(SGI2) = (URef(SGI2)*UPrimRef(SGI3) + U(SGI2)*UPrim(SGI3))                      ! {rho*u*v}
  gTilde(SGI3) = (URef(SGI3)*UPrimRef(SGI3)+UPrimRef(SGI5) + U(SGI3)*UPrim(SGI3)+UPrim(SGI5)) ! {rho*v²}+{p}
#if PP_dim == 3
  gTilde(SGI4) = (URef(SGI3)*UPrimRef(SGI4) + U(SGI3)*UPrim(SGI4))                      ! {rho*v*w}
#else
  gTilde(SGI4) = 0.
#endif
  gTilde(SGI5) = (EpRef*UPrimRef(SGI3) + Ep*UPrim(SGI3))                          ! {(rho*e+p)*v}
#if PP_dim == 3
! local Euler fluxes z-direction
  hTilde(SGI1) = (URef(SGI4) + U(SGI4))                                           ! {rho*w}
  hTilde(SGI2) = (URef(SGI2)*UPrimRef(SGI4) + U(SGI2)*UPrim(SGI4))                      ! {rho*u*w}
  hTilde(SGI3) = (URef(SGI3)*UPrimRef(SGI4) + U(SGI3)*UPrim(SGI4))                      ! {rho*v*w}
  hTilde(SGI4) = (URef(SGI4)*UPrimRef(SGI4)+UPrimRef(SGI5) + U(SGI4)*UPrim(SGI4)+UPrim(SGI5)) ! {rho*v²+p}
  hTilde(SGI5) = (EpRef*UPrimRef(SGI4) + Ep*UPrim(SGI4))                          ! {(rho*e+p)*w}
#endif

! transform into reference space
Flux(:) = 0.5*(MRef(1)+Met(1))*fTilde(:) + &
#if PP_dim == 3
          0.5*(MRef(3)+Met(3))*hTilde(:) + &
#endif
          0.5*(MRef(2)+Met(2))*gTilde(:)

END SUBROUTINE SplitVolumeFluxSD

!==================================================================================================================================
!> Computes the surface flux for the split formulation retaining the standart NS-Equations
!==================================================================================================================================
SUBROUTINE SplitSurfaceFluxSD(U_LL,U_RR,F)
! MODULES
USE MOD_PreProc
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
REAL,DIMENSION(PP_2Var),INTENT(IN)  :: U_LL,U_RR ! variables at the left-/right-Surfaces
REAL,DIMENSION(PP_nVar),INTENT(OUT) :: F ! resulting flux
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!==================================================================================================================================

F(SGI1)= 0.5*(U_LL(MOM1)+U_RR(MOM1))                                                ! {rho*u}
F(SGI2)= 0.5*(U_LL(MOM1)*U_LL(VEL1)+U_LL(PRES)+U_RR(MOM1)*U_RR(VEL1)+U_RR(PRES))    ! {rho*u²}+{p}
F(SGI3)= 0.5*(U_LL(MOM1)*U_LL(VEL2)+U_RR(MOM1)*U_RR(VEL2))                          ! {rho*u*v}
#if PP_dim == 3
F(SGI4)= 0.5*(U_LL(MOM1)*U_LL(VEL3)+U_RR(MOM1)*U_RR(VEL3))                          ! {rho*u*w}
#else
F(SGI4)= 0.
#endif
F(SGI5)= 0.5*((U_LL(ENER)+U_LL(PRES))*U_LL(VEL1)+(U_RR(ENER)+U_RR(PRES))*U_RR(VEL1))! {(rho*e+p)*u}

END SUBROUTINE SplitSurfaceFluxSD

!==================================================================================================================================
!> Computes the Split-Flux retaining the formulation of Ducros
!> Attention 1: Factor 2 from differentiation matrix is already been considered
!==================================================================================================================================
SUBROUTINE SplitVolumeFluxDU(URef,UPrimRef,U,UPrim,MRef,Met,Flux)
! MODULES
USE MOD_PreProc
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
REAL,DIMENSION(PP_nVar    ),INTENT(IN)  :: URef,U ! conserved variables
REAL,DIMENSION(PP_nVarPrim),INTENT(IN)  :: UPrimRef,UPrim ! primitive variables
REAL,DIMENSION(1:3        ),INTENT(IN)  :: MRef,Met ! metric terms
REAL,DIMENSION(PP_nVar    ),INTENT(OUT) :: Flux ! flux in reverence space
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL,DIMENSION(PP_nVar    )             :: fTilde,gTilde ! flux in physical space
#if PP_dim == 3
REAL,DIMENSION(PP_nVar    )             :: hTilde        ! flux in physical space
#endif
!==================================================================================================================================

! local Euler fluxes x-direction
  fTilde(SGI1) = 0.5*(URef(SGI1)+U(SGI1))*(UPrimRef(SGI2)+UPrim(SGI2))                          ! {rho}*{u}
  fTilde(SGI2) = 0.5*(URef(SGI2)+U(SGI2))*(UPrimRef(SGI2)+UPrim(SGI2)) + (UPrimRef(SGI5)+UPrim(SGI5)) ! {rho*u}*{u}+{p}
  fTilde(SGI3) = 0.5*(URef(SGI3)+U(SGI3))*(UPrimRef(SGI2)+UPrim(SGI2))                          ! {rho*v}*{u}
#if PP_dim == 3
  fTilde(SGI4) = 0.5*(URef(SGI4)+U(SGI4))*(UPrimRef(SGI2)+UPrim(SGI2))                          ! {rho*w}*{u}
#else
  fTilde(SGI4) = 0.
#endif
  fTilde(SGI5) = 0.5*(URef(SGI5)+U(SGI5)+UPrimRef(SGI5)+UPrim(SGI5))*(UPrimRef(SGI2)+UPrim(SGI2))     ! ({rho*e}+{p})*{u}
! local Euler fluxes y-direction
  gTilde(SGI1) = 0.5*(URef(SGI1)+U(SGI1))*(UPrimRef(SGI3)+UPrim(SGI3))                          ! {rho}*{v}
  gTilde(SGI2) = 0.5*(URef(SGI2)+U(SGI2))*(UPrimRef(SGI3)+UPrim(SGI3))                          ! {rho*u}*{v}
  gTilde(SGI3) = 0.5*(URef(SGI3)+U(SGI3))*(UPrimRef(SGI3)+UPrim(SGI3)) + (UPrimRef(SGI5)+UPrim(SGI5)) ! {rho*v}*{v}+{p}
#if PP_dim == 3
  gTilde(SGI4) = 0.5*(URef(SGI4)+U(SGI4))*(UPrimRef(SGI3)+UPrim(SGI3))                          ! {rho*w}*{v}
#else
  gTilde(SGI4) = 0.
#endif
  gTilde(SGI5) = 0.5*(URef(SGI5)+U(SGI5)+UPrimRef(SGI5)+UPrim(SGI5))*(UPrimRef(SGI3)+UPrim(SGI3))     ! ({rho*e}+{p})*{v}
#if PP_dim == 3
! local Euler fluxes z-direction
  hTilde(SGI1) = 0.5*(URef(SGI1)+U(SGI1))*(UPrimRef(SGI4)+UPrim(SGI4))                          ! {rho}*{w}
  hTilde(SGI2) = 0.5*(URef(SGI2)+U(SGI2))*(UPrimRef(SGI4)+UPrim(SGI4))                          ! {rho*u}*{w}
  hTilde(SGI3) = 0.5*(URef(SGI3)+U(SGI3))*(UPrimRef(SGI4)+UPrim(SGI4))                          ! {rho*v}*{w}
  hTilde(SGI4) = 0.5*(URef(SGI4)+U(SGI4))*(UPrimRef(SGI4)+UPrim(SGI4)) + (UPrimRef(SGI5)+UPrim(SGI5)) ! {rho*w}*{w}+{p}
  hTilde(SGI5) = 0.5*(URef(SGI5)+U(SGI5)+UPrimRef(SGI5)+UPrim(SGI5))*(UPrimRef(SGI4)+UPrim(SGI4))     ! ({rho*e}+{p})*{w}
#endif

! transform into reference space
Flux(:) = 0.5*(MRef(1)+Met(1))*fTilde(:) + &
#if PP_dim == 3
          0.5*(MRef(3)+Met(3))*hTilde(:) + &
#endif
          0.5*(MRef(2)+Met(2))*gTilde(:)

END SUBROUTINE SplitVolumeFluxDU

!==================================================================================================================================
!> Computes the surface flux for the split formulation of Ducros
!==================================================================================================================================
SUBROUTINE SplitSurfaceFluxDU(U_LL,U_RR,F)
! MODULES
USE MOD_PreProc
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
REAL,DIMENSION(PP_2Var),INTENT(IN)  :: U_LL,U_RR ! variables at the left-/right-Surfaces
REAL,DIMENSION(PP_nVar),INTENT(OUT) :: F ! resulting flux
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!==================================================================================================================================
F(SGI1)= 0.25*(U_LL(DENS)+U_RR(DENS))*(U_LL(VEL1)+U_RR(VEL1))                               ! {rho}*{u}
F(SGI2)= 0.25*(U_LL(MOM1)+U_RR(MOM1))*(U_LL(VEL1)+U_RR(VEL1)) + 0.5*(U_LL(PRES)+U_RR(PRES)) ! {rho*u}*{u}+{p}
F(SGI3)= 0.25*(U_LL(MOM2)+U_RR(MOM2))*(U_LL(VEL1)+U_RR(VEL1))                               ! {rho*v}*{u}
#if PP_dim == 3
F(SGI4)= 0.25*(U_LL(MOM3)+U_RR(MOM3))*(U_LL(VEL1)+U_RR(VEL1))                               ! {rho*w}*{u}
#else
F(SGI4)= 0.
#endif
F(SGI5)= 0.25*(U_LL(ENER)+U_RR(ENER)+U_LL(PRES)+U_RR(PRES))*(U_LL(VEL1)+U_RR(VEL1))         ! ({rho*e}+{p})*{u}

END SUBROUTINE SplitSurfaceFluxDU

!==================================================================================================================================
!> Computes the Split-Flux retaining the formulation of Kennedy and Gruber
!> Attention 1: Factor 2 from differentiation matrix is already been considered
!==================================================================================================================================
SUBROUTINE SplitVolumeFluxKG(URef,UPrimRef,U,UPrim,MRef,Met,Flux)
! MODULES
USE MOD_PreProc
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
REAL,DIMENSION(PP_nVar    ),INTENT(IN)  :: URef,U ! conserved variables
REAL,DIMENSION(PP_nVarPrim),INTENT(IN)  :: UPrimRef,UPrim ! primitive variables
REAL,DIMENSION(1:3        ),INTENT(IN)  :: MRef,Met ! metric terms
REAL,DIMENSION(PP_nVar    ),INTENT(OUT) :: Flux ! flux in reverence space
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                                    :: e,eRef ! auxilery variables for the specific energy
REAL,DIMENSION(PP_nVar    )             :: fTilde,gTilde ! flux in physical space
#if PP_dim == 3
REAL,DIMENSION(PP_nVar    )             :: hTilde        ! flux in physical space
#endif
!==================================================================================================================================

! specific energy
eRef = URef(SGI5)/URef(SGI1)
e    = U(SGI5)/U(SGI1)

! local Euler fluxes x-direction
  fTilde(SGI1) = 0.5* (URef(SGI1)+U(SGI1))*(UPrimRef(SGI2)+UPrim(SGI2))                             ! {rho}*{u}
  fTilde(SGI2) = 0.25*(URef(SGI1)+U(SGI1))*(UPrimRef(SGI2)+UPrim(SGI2))**2 + (UPrimRef(SGI5)+UPrim(SGI5)) ! {rho}*{u}²+{p}
  fTilde(SGI3) = 0.25*(URef(SGI1)+U(SGI1))*(UPrimRef(SGI2)+UPrim(SGI2))*(UPrimRef(SGI3)+UPrim(SGI3))      ! {rho}*{u}*{v}
#if PP_dim == 3
  fTilde(SGI4) = 0.25*(URef(SGI1)+U(SGI1))*(UPrimRef(SGI2)+UPrim(SGI2))*(UPrimRef(SGI4)+UPrim(SGI4))      ! {rho}*{u}*{w}
#else
  fTilde(SGI4) = 0.
#endif
  fTilde(SGI5) = 0.25*(URef(SGI1)+U(SGI1))*(UPrimRef(SGI2)+UPrim(SGI2))*(eRef+e) + &
              0.5* (UPrimRef(SGI5)+UPrim(SGI5))*(UPrimRef(SGI2)+UPrim(SGI2))                     ! {rho}*{e}*{u}+{p}*{u}
! local Euler fluxes y-direction
  gTilde(SGI1) = 0.5 *(URef(SGI1)+U(SGI1))*(UPrimRef(SGI3)+UPrim(SGI3))                             ! {rho}*{v}
  gTilde(SGI2) = fTilde(SGI3)                                                              ! {rho}*{v}*{u}
  gTilde(SGI3) = 0.25*(URef(SGI1)+U(SGI1))*(UPrimRef(SGI3)+UPrim(SGI3))**2 + (UPrimRef(SGI5)+UPrim(SGI5)) ! {rho}*{v}²+{p}
#if PP_dim == 3
  gTilde(SGI4) = 0.25*(URef(SGI1)+U(SGI1))*(UPrimRef(SGI3)+UPrim(SGI3))*(UPrimRef(SGI4)+UPrim(SGI4))      ! {rho}*{v}*{w}
#else
  gTilde(SGI4) = 0.
#endif
  gTilde(SGI5) = 0.25*(URef(SGI1)+U(SGI1))*(UPrimRef(SGI3)+UPrim(SGI3))*(eRef+e) + &
              0.5* (UPrimRef(SGI5)+UPrim(SGI5))*(UPrimRef(SGI3)+UPrim(SGI3))                     ! {rho}*{e}*{v}+{p}*{v}
#if PP_dim == 3
! local Euler fluxes z-direction
  hTilde(SGI1) = 0.5 *(URef(SGI1)+U(SGI1))*(UPrimRef(SGI4)+UPrim(SGI4))                             ! {rho}*{w}
  hTilde(SGI2) = fTilde(SGI4)                                                              ! {rho}*{w}*{u}
  hTilde(SGI3) = gTilde(SGI4)                                                              ! {rho}*{w}*{v}
  hTilde(SGI4) = 0.25*(URef(SGI1)+U(SGI1))*(UPrimRef(SGI4)+UPrim(SGI4))**2 + (UPrimRef(SGI5)+UPrim(SGI5)) ! {rho}*{w}²+{p}
  hTilde(SGI5) = 0.25*(URef(SGI1)+U(SGI1))*(UPrimRef(SGI4)+UPrim(SGI4))*(eRef+e) + &
              0.5 *(UPrimRef(SGI5)+UPrim(SGI5))*(UPrimRef(SGI4)+UPrim(SGI4))                     ! {rho}*{e}*{w}+{p}*{w}
#endif

! transform into reference space
Flux(:) = 0.5*(MRef(1)+Met(1))*fTilde(:) + &
#if PP_dim == 3
          0.5*(MRef(3)+Met(3))*hTilde(:) + &
#endif
          0.5*(MRef(2)+Met(2))*gTilde(:)

END SUBROUTINE SplitVolumeFluxKG

!==================================================================================================================================
!> Computes the surface flux for the split formulation of Kennedy and Gruber
!==================================================================================================================================
SUBROUTINE SplitSurfaceFluxKG(U_LL,U_RR,F)
! MODULES
USE MOD_PreProc
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
REAL,DIMENSION(PP_2Var),INTENT(IN)  :: U_LL,U_RR ! variables at the left-/right-Surfaces
REAL,DIMENSION(PP_nVar),INTENT(OUT) :: F ! resulting flux
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                                :: e_LL,e_RR ! auxilery variables for the specific energy
!==================================================================================================================================
! specific energy
e_LL = U_LL(ENER)/U_LL(DENS)
e_RR = U_RR(ENER)/U_RR(DENS)
!compute flux
F(SGI1)= 0.25* (U_LL(DENS)+U_RR(DENS))*(U_LL(VEL1)+U_RR(VEL1))                                  ! {rho}*{u}
F(SGI2)= 0.125*(U_LL(DENS)+U_RR(DENS))*(U_LL(VEL1)+U_RR(VEL1))**2 + 0.5*(U_LL(PRES)+U_RR(PRES)) ! {rho}*{u}²+{p}
F(SGI3)= 0.125*(U_LL(DENS)+U_RR(DENS))*(U_LL(VEL1)+U_RR(VEL1))*(U_LL(VEL2)+U_RR(VEL2))          ! {rho}*{u}*{v}
#if PP_dim == 3
F(SGI4)= 0.125*(U_LL(DENS)+U_RR(DENS))*(U_LL(VEL1)+U_RR(VEL1))*(U_LL(VEL3)+U_RR(VEL3))          ! {rho}*{u}*{w}
#else
F(SGI4)= 0.
#endif
F(SGI5)= 0.125*(U_LL(DENS)+U_RR(DENS))*(e_LL+e_RR)*(U_LL(VEL1)+U_RR(VEL1)) + &
      0.25 *(U_LL(PRES)+U_RR(PRES))*(U_LL(VEL1)+U_RR(VEL1))                                  ! {rho}*{e}*{u}+{p}*{u}

END SUBROUTINE SplitSurfaceFluxKG

!==================================================================================================================================
!> Computes the Split-Flux retaining the formulation of Morinishi
!> Attention 1: Factor 2 from differentiation matrix is already been considered
!==================================================================================================================================
SUBROUTINE SplitVolumeFluxMO(URef,UPrimRef,U,UPrim,MRef,Met,Flux)
! MODULES
USE MOD_PreProc
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
REAL,DIMENSION(PP_nVar    ),INTENT(IN)  :: URef,U ! conserved variables
REAL,DIMENSION(PP_nVarPrim),INTENT(IN)  :: UPrimRef,UPrim ! primitive variables
REAL,DIMENSION(1:3        ),INTENT(IN)  :: MRef,Met ! metric terms
REAL,DIMENSION(PP_nVar    ),INTENT(OUT) :: Flux ! flux in reverence space
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                                    :: EpRef,Ep ! auxilery variable for inner energy + pressure
REAL,DIMENSION(PP_nVar    )             :: fTilde,gTilde ! flux in physical space
#if PP_dim == 3
REAL,DIMENSION(PP_nVar    )             :: hTilde        ! flux in physical space
#endif
!==================================================================================================================================

! internal energy + pressure
EpRef = URef(SGI5)-0.5*URef(SGI1)*(UPrimRef(SGI2)**2+UPrimRef(SGI3)**2+UPrimRef(SGI4)**2)+UPrimRef(SGI5)
Ep    = U(SGI5)-0.5*U(SGI1)*(UPrim(SGI2)**2+UPrim(SGI3)**2+UPrim(SGI4)**2)+UPrim(SGI5)

! local Euler fluxes x-direction
  fTilde(SGI1) =     (URef(SGI2)+U(SGI2))                                                 ! {rho*u}
  fTilde(SGI2) = 0.5*(URef(SGI2)+U(SGI2))*(UPrimRef(SGI2)+UPrim(SGI2)) + (UPrimRef(SGI5)+UPrim(SGI5)) ! {rho*u}*{u}+{p}
  fTilde(SGI3) = 0.5*(Uref(SGI2)+U(SGI2))*(UPrimRef(SGI3)+UPrim(SGI3))                          ! {rho*u}*{v}
#if PP_dim == 3
  fTilde(SGI4) = 0.5*(URef(SGI2)+U(SGI2))*(UPrimRef(SGI4)+UPrim(SGI4))                          ! {rho*u}*{w}
#else
  fTilde(SGI4) = 0.
#endif
  fTilde(SGI5) = (EpRef*UPrimRef(SGI2)+Ep*UPrim(SGI2)) + &                                !{(rho*e_int+p)*u} +
              0.5*(URef(SGI2)*UPrimRef(SGI2)+U(SGI2)*UPrim(SGI2))*(UPrimRef(SGI2)+UPrim(SGI2)) + & !{rho*u²}*{u} +
              0.5*(URef(SGI2)*UPrimRef(SGI3)+U(SGI2)*UPrim(SGI3))*(UPrimRef(SGI3)+UPrim(SGI3)) + & !{rho*u*v}*{v} +
              0.5*(URef(SGI2)*UPrimRef(SGI4)+U(SGI2)*UPrim(SGI4))*(UPrimRef(SGI4)+UPrim(SGI4)) - & !{rho*u*w}*{w} -
              0.5*(URef(SGI2)*UPrimRef(SGI2)*UPrimRef(SGI2)+U(SGI2)*UPrim(SGI2)*UPrim(SGI2))   - & !1/2*({rho*u³} +
              0.5*(URef(SGI2)*UPrimRef(SGI3)*UPrimRef(SGI3)+U(SGI2)*UPrim(SGI3)*UPrim(SGI3))   - & !{rho*u*v²} +
              0.5*(URef(SGI2)*UPrimRef(SGI4)*UPrimRef(SGI4)+U(SGI2)*UPrim(SGI4)*UPrim(SGI4))       !{rho*u*w²})
! local Euler fluxes y-direction
  gTilde(SGI1) =     (URef(SGI3)+U(SGI3))                                                 ! {rho*v}
  gTilde(SGI2) = 0.5*(Uref(SGI3)+U(SGI3))*(UPrimRef(SGI2)+UPrim(SGI2))                          ! {rho*v}*{u}
  gTilde(SGI3) = 0.5*(URef(SGI3)+U(SGI3))*(UPrimRef(SGI3)+UPrim(SGI3)) + (UPrimRef(SGI5)+UPrim(SGI5)) ! {rho*v}*{v}+{p}
#if PP_dim == 3
  gTilde(SGI4) = 0.5*(URef(SGI3)+U(SGI3))*(UPrimRef(SGI4)+UPrim(SGI4))                          ! {rho*v}*{w}
#else
  gTilde(SGI4) = 0.
#endif
  gTilde(SGI5) = (EpRef*UPrimRef(SGI3)+Ep*UPrim(SGI3)) + &                                !{(rho*e_int+p)*v} +
              0.5*(URef(SGI3)*UPrimRef(SGI2)+U(SGI3)*UPrim(SGI2))*(UPrimRef(SGI2)+UPrim(SGI2)) + & !{rho*v*u}*{u} +
              0.5*(URef(SGI3)*UPrimRef(SGI3)+U(SGI3)*UPrim(SGI3))*(UPrimRef(SGI3)+UPrim(SGI3)) + & !{rho*v²}*{v} +
              0.5*(URef(SGI3)*UPrimRef(SGI4)+U(SGI3)*UPrim(SGI4))*(UPrimRef(SGI4)+UPrim(SGI4)) - & !{rho*v*w}*{w} -
              0.5*(URef(SGI3)*UPrimRef(SGI2)*UPrimRef(SGI2)+U(SGI3)*UPrim(SGI2)*UPrim(SGI2))   - & !1/2*({rho*v*u²} +
              0.5*(URef(SGI3)*UPrimRef(SGI3)*UPrimRef(SGI3)+U(SGI3)*UPrim(SGI3)*UPrim(SGI3))   - & !{rho*v³} +
              0.5*(URef(SGI3)*UPrimRef(SGI4)*UPrimRef(SGI4)+U(SGI3)*UPrim(SGI4)*UPrim(SGI4))       !{rho*v*w²})
#if PP_dim == 3
! local Euler fluxes z-direction
  hTilde(SGI1) =     (URef(SGI4)+U(SGI4))                                                 ! {rho*w}
  hTilde(SGI2) = 0.5*(Uref(SGI4)+U(SGI4))*(UPrimRef(SGI2)+UPrim(SGI2))                          ! {rho*w}*{u}
  hTilde(SGI3) = 0.5*(URef(SGI4)+U(SGI4))*(UPrimRef(SGI3)+UPrim(SGI3))                          ! {rho*w}*{v}
  hTilde(SGI4) = 0.5*(URef(SGI4)+U(SGI4))*(UPrimRef(SGI4)+UPrim(SGI4)) + (UPrimRef(SGI5)+UPrim(SGI5)) ! {rho*w}*{w}+{p}
  hTilde(SGI5) = (EpRef*UPrimRef(SGI4)+Ep*UPrim(SGI4)) + &                                !{(rho*e_int+p)*w} +
              0.5*(URef(SGI4)*UPrimRef(SGI2)+U(SGI4)*UPrim(SGI2))*(UPrimRef(SGI2)+UPrim(SGI2)) + & !{rho*w*u}*{u} +
              0.5*(URef(SGI4)*UPrimRef(SGI3)+U(SGI4)*UPrim(SGI3))*(UPrimRef(SGI3)+UPrim(SGI3)) + & !{rho*w*v}*{v} +
              0.5*(URef(SGI4)*UPrimRef(SGI4)+U(SGI4)*UPrim(SGI4))*(UPrimRef(SGI4)+UPrim(SGI4)) - & !{rho*w²}*{w} -
              0.5*(URef(SGI4)*UPrimRef(SGI2)*UPrimRef(SGI2)+U(SGI4)*UPrim(SGI2)*UPrim(SGI2))   - & !1/2*({rho*w*u²} +
              0.5*(URef(SGI4)*UPrimRef(SGI3)*UPrimRef(SGI3)+U(SGI4)*UPrim(SGI3)*UPrim(SGI3))   - & !{rho*w*v²} +
              0.5*(URef(SGI4)*UPrimRef(SGI4)*UPrimRef(SGI4)+U(SGI4)*UPrim(SGI4)*UPrim(SGI4))       !{rho*w³})
#endif

! transform into reference space
Flux(:) = 0.5*(MRef(1)+Met(1))*fTilde(:) + &
#if PP_dim == 3
          0.5*(MRef(3)+Met(3))*hTilde(:) + &
#endif
          0.5*(MRef(2)+Met(2))*gTilde(:)

END SUBROUTINE SplitVolumeFluxMO

!==================================================================================================================================
!> Computes the surface flux for the split formulation of Morinishi
!==================================================================================================================================
SUBROUTINE SplitSurfaceFluxMO(U_LL,U_RR,F)
! MODULES
USE MOD_PreProc
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
REAL,DIMENSION(PP_2Var),INTENT(IN)  :: U_LL,U_RR ! variables at the left-/right-Surfaces
REAL,DIMENSION(PP_nVar),INTENT(OUT) :: F ! resulting flux
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                                :: Ep_LL,Ep_RR
!==================================================================================================================================
! internal energy + pressure
Ep_LL = U_LL(ENER)-0.5*U_LL(DENS)*(U_LL(VEL1)**2+U_LL(VEL2)**2+U_LL(VEL3)**2)+U_LL(PRES)
Ep_RR = U_RR(ENER)-0.5*U_RR(DENS)*(U_RR(VEL1)**2+U_RR(VEL2)**2+U_RR(VEL3)**2)+U_RR(PRES)

! compute flux
F(SGI1)= 0.5 *(U_LL(MOM1)+U_RR(MOM1))                                                       ! {rho*u}
F(SGI2)= 0.25*(U_LL(MOM1)+U_RR(MOM1))*(U_LL(VEL1)+U_RR(VEL1)) + 0.5*(U_LL(PRES)+U_RR(PRES)) ! {rho*u}*{u}+{p}
F(SGI3)= 0.25*(U_LL(MOM1)+U_RR(MOM1))*(U_LL(VEL2)+U_RR(VEL2))                               ! {rho*u}*{v}
#if PP_dim == 3
F(SGI4)= 0.25*(U_LL(MOM1)+U_RR(MOM1))*(U_LL(VEL3)+U_RR(VEL3))                               ! {rho*u}*{w}
#else
F(SGI4)= 0.
#endif
F(5)= 0.5 *(Ep_LL*U_LL(VEL1)+Ep_RR*U_RR(VEL1)) +  &                                      !{(rho*e_int+p)*u} +
      0.25*(U_LL(MOM1)*U_LL(VEL1)+U_RR(MOM1)*U_RR(VEL1))*(U_LL(VEL1)+U_RR(VEL1)) + &     !{rho*u²}*{u} +
      0.25*(U_LL(MOM1)*U_LL(VEL2)+U_RR(MOM1)*U_RR(VEL2))*(U_LL(VEL2)+U_RR(VEL2)) + &     !{rho*u*v}*{v} +
      0.25*(U_LL(MOM1)*U_LL(VEL3)+U_RR(MOM1)*U_RR(VEL3))*(U_LL(VEL3)+U_RR(VEL3)) - &     !{rho*u*w}*{w} -
      0.25*(U_LL(MOM1)*U_LL(VEL1)*U_LL(VEL1)+U_RR(MOM1)*U_RR(VEL1)*U_RR(VEL1)) - &       !1/2*({rho*u³} +
      0.25*(U_LL(MOM1)*U_LL(VEL2)*U_LL(VEL2)+U_RR(MOM1)*U_RR(VEL2)*U_RR(VEL2)) - &       !{rho*u*v²} +
      0.25*(U_LL(MOM1)*U_LL(VEL3)*U_LL(VEL3)+U_RR(MOM1)*U_RR(VEL3)*U_RR(VEL3))           !{rho*u*w²})

END SUBROUTINE SplitSurfaceFluxMO
#endif /*SGdummy*/

!==================================================================================================================================
!> Computes the Split-Flux retaining the formulation of Pirozzoli
!> Attention 1: Factor 2 from differentiation matrix is already been considered
!==================================================================================================================================
SUBROUTINE SplitVolumeFluxPI(URef,UPrimRef,U,UPrim,MRef,Met,Flux)
! MODULES
USE MOD_PreProc
USE MOD_SG_Operators  ,ONLY:SG_Product,SG_Fraction
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
REAL,DIMENSION(PP_nVar    ),INTENT(IN)  :: URef,U ! conserved variables
REAL,DIMENSION(PP_nVarPrim),INTENT(IN)  :: UPrimRef,UPrim ! primitive variables
REAL,DIMENSION(1:3        ),INTENT(IN)  :: MRef,Met ! metric terms
REAL,DIMENSION(PP_nVar    ),INTENT(OUT) :: Flux ! flux in reverence space
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL,DIMENSION(0:PP_nCoefM1)                        :: e,eRef ! auxilery variables for the specific enthalpy
REAL,DIMENSION(PP_nVar    )             :: fTilde,gTilde ! flux in physical space
#if PP_dim == 3
REAL,DIMENSION(PP_nVar    )             :: hTilde        ! flux in physical space
#endif
!==================================================================================================================================

! specific enthalpy
eRef = SG_Fraction((URef(SGI5)+UPrimRef(SGI5)),URef(SGI1))
e    = SG_Fraction((U(SGI5)+UPrim(SGI5)),U(SGI1))

! local Euler fluxes x-direction
  fTilde(SGI1) = 0.5 *SG_Product((URef(SGI1)+U(SGI1)),(UPrimRef(SGI2)+UPrim(SGI2)))                             ! {rho}*{u}
  fTilde(SGI2) = SG_Product((URef(SGI1)+U(SGI1)),(UPrimRef(SGI2)+UPrim(SGI2)))  ! {rho}*{u}
  fTilde(SGI2) = 0.25*(SG_Product(fTilde(SGI2),(UPrimRef(SGI2)+UPrim(SGI2)))) + (UPrimRef(SGI5)+UPrim(SGI5)) ! ({rho}*{u})*{u}+{p}
  fTilde(SGI3) = SG_Product((URef(SGI1)+U(SGI1)),(UPrimRef(SGI2)+UPrim(SGI2)))   ! {rho}*{u}
  fTilde(SGI3) = 0.25*SG_Product(fTilde(SGI3),(UPrimRef(SGI3)+UPrim(SGI3)))      ! ({rho}*{u})*{v}
#if PP_dim == 3
  fTilde(SGI4) = SG_Product((URef(SGI1)+U(SGI1)),(UPrimRef(SGI2)+UPrim(SGI2)))                   ! {rho}*{u}
  fTilde(SGI4) = 0.25*SG_Product(fTilde(SGI4),(UPrimRef(SGI4)+UPrim(SGI4)))                      ! ({rho}*{u})*{w}
#else
  fTilde(SGI4) = 0.
#endif
  fTilde(SGI5) = SG_Product((URef(SGI1)+U(SGI1)),(UPrimRef(SGI2))+UPrim(SGI2))                    ! {rho}*{h}
  fTilde(SGI5) = 0.25*SG_Product(fTilde(SGI5),(eRef+e))                                           ! ({rho}*{h})*{u}
! local Euler fluxes y-direction
  gTilde(SGI1) = 0.5 *SG_Product((URef(SGI1)+U(SGI1)),(UPrimRef(SGI3)+UPrim(SGI3)))                             ! {rho}*{v}
  gTilde(SGI2) = fTilde(SGI3)                                                              ! {rho}*{v}*{u}
  gTilde(SGI3) = SG_Product((URef(SGI1)+U(SGI1)),(UPrimRef(SGI3)+UPrim(SGI3))) ! {rho}*{v}²+{p}
  gTilde(SGI3) = 0.25*SG_Product(gTilde(SGI3),(UPrimRef(SGI3)+UPrim(SGI3))) + (UPrimRef(SGI5)+UPrim(SGI5)) ! {rho}*{v}²+{p}
#if PP_dim == 3
  gTilde(SGI4) = SG_Product((URef(SGI1)+U(SGI1)),(UPrimRef(SGI3)+UPrim(SGI3)))             ! {rho}*{v}*{w}
  gTilde(SGI4) = 0.25*SG_Product(gTilde(SGI4),(UPrimRef(SGI4)+UPrim(SGI4)))      ! {rho}*{v}*{w}
#else
  gTilde(SGI4) = 0.
#endif
  gTilde(SGI5) = SG_Product((URef(SGI1)+U(SGI1)),(UPrimRef(SGI3)+UPrim(SGI3)))                    ! {rho}*{h}*{v}
  gTilde(SGI5) = 0.25*SG_Product(gTilde(SGI5),(eRef+e))                    ! {rho}*{h}*{v}
#if PP_dim == 3
! local Euler fluxes z-direction
  hTilde(SGI1) = 0.5 *SG_Product((URef(SGI1)+U(SGI1)),(UPrimRef(SGI4)+UPrim(SGI4)))                             ! {rho}*{w}
  hTilde(SGI2) = fTilde(SGI4)                                                              ! {rho}*{w}*{u}
  hTilde(SGI3) = gTilde(SGI4)                                                              ! {rho}*{w}*{v}
  hTilde(SGI4) = SG_Product((URef(SGI1)+U(SGI1)),(UPrimRef(SGI4)+UPrim(SGI4))) ! {rho}*{w}²+{p}
  hTilde(SGI4) = 0.25*SG_Product(hTilde(SGI4),(UPrimRef(SGI4)+UPrim(SGI4))) + (UPrimRef(SGI5)+UPrim(SGI5)) ! {rho}*{w}²+{p}
  hTilde(SGI5) = SG_Product((URef(SGI1)+U(SGI1)),(UPrimRef(SGI4)+UPrim(SGI4)))                    ! {rho}*{h}*{w}
  hTilde(SGI5) = 0.25*SG_Product(hTilde(SGI5),(eRef+e))                    ! {rho}*{h}*{w}
#endif

! transform into reference space
Flux(:) = 0.5*(MRef(1)+Met(1))*fTilde(:) + &
#if PP_dim == 3
          0.5*(MRef(3)+Met(3))*hTilde(:) + &
#endif
          0.5*(MRef(2)+Met(2))*gTilde(:)

END SUBROUTINE SplitVolumeFluxPI

!==================================================================================================================================
!> Computes the surface flux for the split formulation of Pirozzoli
!==================================================================================================================================
SUBROUTINE SplitSurfaceFluxPI(U_LL,U_RR,F)
! MODULES
USE MOD_PreProc
USE MOD_SG_Operators  ,ONLY:SG_Product,SG_Fraction
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
REAL,DIMENSION(PP_2Var),INTENT(IN)  :: U_LL,U_RR ! variables at the left-/right-Surfaces
REAL,DIMENSION(PP_nVar),INTENT(OUT) :: F ! resulting flux
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL ,DIMENSION(0:PP_nCoefM1)            :: e_LL,e_RR ! auxilery variables for the specific energy
!==================================================================================================================================
! specific energy
e_LL = SG_Fraction((U_LL(ENER)+U_LL(PRES)),U_LL(DENS))
e_RR = SG_Fraction((U_RR(ENER)+U_RR(PRES)),U_RR(DENS))
!compute flux
F(SGI1)= 0.25* SG_Product((U_LL(DENS)+U_RR(DENS)),(U_LL(VEL1)+U_RR(VEL1)))                                  ! {rho}*{u}
F(SGI2)= SG_Product((U_LL(DENS)+U_RR(DENS)),(U_LL(VEL1)+U_RR(VEL1))) ! {rho}*{u}²+{p}
F(SGI2)= 0.125*SG_Product(F(SGI2),(U_LL(VEL1)+U_RR(VEL1))) + 0.5*(U_LL(PRES)+U_RR(PRES)) ! {rho}*{u}²+{p}
F(SGI3)= SG_Product((U_LL(DENS)+U_RR(DENS)),(U_LL(VEL1)+U_RR(VEL1)))          ! {rho}*{u}*{v}
F(SGI3)= 0.125*SG_Product(F(SGI3),(U_LL(VEL2)+U_RR(VEL2)))          ! {rho}*{u}*{v}
#if PP_dim == 3
F(SGI4)= SG_Product((U_LL(DENS)+U_RR(DENS)),(U_LL(VEL1)+U_RR(VEL1)))          ! {rho}*{u}*{w}
F(SGI4)= 0.125*SG_Product(F(SGI4),(U_LL(VEL3)+U_RR(VEL3)))          ! {rho}*{u}*{w}
#else
F(SGI4)= 0.
#endif
F(SGI5)= SG_Product((U_LL(DENS)+U_RR(DENS)),(e_LL+e_RR))                      ! {rho}*{h}*{u}
F(SGI5)= 0.125*SG_Product(F(SGI5),(U_LL(VEL1)+U_RR(VEL1)))                      ! {rho}*{h}*{u}

END SUBROUTINE SplitSurfaceFluxPI
END MODULE MOD_SplitFlux
