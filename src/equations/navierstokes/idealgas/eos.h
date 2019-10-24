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
! Define variables for normal and extended state vector
! Normal   U(1:5)  with conservative variables
! Extended U(1:11) with conservative and primitive variables

#define PP_2VarDet 11
#define PP_2Var    11*PP_nCoef
#define PP_2VarQ   11*nQPTotal
#define PP_nVarQ   6*nQPTotal

#define CONS              1: 5*PP_nCoef /* conservative variables */
#define PRIM   5*PP_nCoef+1:11*PP_nCoef /* primitive variables */

! conservative variables
#define DENS              1:   PP_nCoef /* density */
#define MOM1     PP_nCoef+1: 2*PP_nCoef /* momentum x */
#define MOM2   2*PP_nCoef+1: 3*PP_nCoef /* momentum y */
#define MOM3   3*PP_nCoef+1: 4*PP_nCoef /* momentum z */
#define MOMV     PP_nCoef+1: 4*PP_nCoef /* momentum vector */
#define MMV2     PP_nCoef+1: (1+PP_dim)*PP_nCoef /* momentum vector */
#define ENER   4*PP_nCoef+1: 5*PP_nCoef /* energy */

! general SG index vectors for det variables 1-6
#define SGI1              1:   PP_nCoef
#define SGI2     PP_nCoef+1: 2*PP_nCoef
#define SGI3   2*PP_nCoef+1: 3*PP_nCoef
#define SGI4   3*PP_nCoef+1: 4*PP_nCoef
#define SGI5   4*PP_nCoef+1: 5*PP_nCoef
#define SGI6   5*PP_nCoef+1: 6*PP_nCoef

! primitive (extended) variables
#define SRHO   5*PP_nCoef+1: 6*PP_nCoef /* specific volume (1./density) */
#define VEL1   6*PP_nCoef+1: 7*PP_nCoef /* velocity x */
#define VEL2   7*PP_nCoef+1: 8*PP_nCoef /* velocity y */
#define VEL3   8*PP_nCoef+1: 9*PP_nCoef /* velocity z */
#define VELV   6*PP_nCoef+1: 9*PP_nCoef /* velocity range */
#define VLV2   6*PP_nCoef+1: (6+PP_dim)*PP_nCoef /* velocity range */
#define PRES   9*PP_nCoef+1:10*PP_nCoef /* pressure */
#define TEMP  10*PP_nCoef+1:11*PP_nCoef /* temperature */

! primitive (non-extended) variables
#define PDENS             1: 1*PP_nCoef /* specific volume (1./density) */
#define PVEL1    PP_nCoef+1: 2*PP_nCoef /* velocity x */
#define PVEL2  2*PP_nCoef+1: 3*PP_nCoef /* velocity y */
#define PVEL3  3*PP_nCoef+1: 4*PP_nCoef /* velocity z */
#define PVELV    PP_nCoef+1: 4*PP_nCoef /* velocity range */
#define PVLV2    PP_nCoef+1: (1+PP_dim)*PP_nCoef /* velocity range */
#define PPRES  4*PP_nCoef+1: 5*PP_nCoef /* pressure */
#define PTEMP  5*PP_nCoef+1: 6*PP_nCoef /* temperature */

! conservative variables quadrature points
#define QDENS              1:   nQPTotal /* density */
#define QMOM1     nQPTotal+1: 2*nQPTotal /* momentum x */
#define QMOM2   2*nQPTotal+1: 3*nQPTotal /* momentum y */
#define QMOM3   3*nQPTotal+1: 4*nQPTotal /* momentum z */
#define QMOMV     nQPTotal+1: 4*nQPTotal /* momentum vector */
#define QMMV2     nQPTotal+1: (1+PP_dim)*nQPTotal /* momentum vector */
#define QENER   4*nQPTotal+1: 5*nQPTotal /* energy */

! primitive (non-extended) variables quadrature points
#define QPSRHO             1: 1*nQPTotal /* specific volume (1./density) */
#define QPVEL1    nQPTotal+1: 2*nQPTotal /* velocity x */
#define QPVEL2  2*nQPTotal+1: 3*nQPTotal /* velocity y */
#define QPVEL3  3*nQPTotal+1: 4*nQPTotal /* velocity z */
#define QPVELV    nQPTotal+1: 4*nQPTotal /* velocity range */
#define QPVLV2    nQPTotal+1: (1+PP_dim)*nQPTotal /* velocity range */
#define QPPRES  4*nQPTotal+1: 5*nQPTotal /* pressure */
#define QPTEMP  5*nQPTotal+1: 6*nQPTotal /* temperature */

! conservative variables quadrature points
#define QDENS              1:   nQPTotal /* density */
#define QMOM1     nQPTotal+1: 2*nQPTotal /* momentum x */
#define QMOM2   2*nQPTotal+1: 3*nQPTotal /* momentum y */
#define QMOM3   3*nQPTotal+1: 4*nQPTotal /* momentum z */
#define QMOMV     nQPTotal+1: 4*nQPTotal /* momentum vector */
#define QMMV2     nQPTotal+1: (1+PP_dim)*nQPTotal /* momentum vector */
#define QENER   4*nQPTotal+1: 5*nQPTotal /* energy */
#define QSRHO   5*nQPTotal+1: 6*nQPTotal /* specific volume (1./density) */
#define QVEL1   6*nQPTotal+1: 7*nQPTotal /* velocity x */
#define QVEL2   7*nQPTotal+1: 8*nQPTotal /* velocity y */
#define QVEL3   8*nQPTotal+1: 9*nQPTotal /* velocity z */
#define QVELV   6*nQPTotal+1: 9*nQPTotal /* velocity range */
#define QVLV2   6*nQPTotal+1: (6+PP_dim)*nQPTotal /* velocity range */
#define QPRES   9*nQPTotal+1:10*nQPTotal /* pressure */
#define QTEMP  10*nQPTotal+1:11*nQPTotal /* temperature */

! deterministic conservative and primitive variables
#define DCONS 1:5  /* conservative variables */
#define DPRIM 6:11 /* primitive variables */

! deterministic conservative variables
#define DDENS  1   /* density */
#define DMOM1  2   /* momentum x */
#define DMOM2  3   /* momentum y */
#define DMOM3  4   /* momentum z */
#define DMOMV  2:4 /* momentum vector */
#define DMMV2  2:1+PP_dim /* momentum vector */
#define DENER  5   /* energy */

! deterministic primitive (extended) variables
#define DSRHO  6   /* specific volume (1./density) */
#define DVEL1  7   /* velocity x */
#define DVEL2  8   /* velocity y */
#define DVEL3  9   /* velocity z */
#define DVELV  7:9 /* velocity range */
#define DVLV2  7:6+PP_dim /* velocity range */
#define DPRES  10  /* pressure */
#define DTEMP  11  /* temperature */


! routines to compute physical quantities
#define KAPPASPR_MAX_TIMESTEP_H()      (MAX(4./3.,KappasPr))
#define THERMAL_CONDUCTIVITY_H(mu)     (mu*cp/Pr)
#define TOTAL_TEMPERATURE_H(T,Mach)    (T*(1+0.5*(kappa-1)*Mach**2))
#define TOTAL_PRESSURE_H(p,Mach)       (p/((1+0.5*(kappa-1)*Mach**2)**(-kappa/(kappa-1.))))
#define BETA_RIEMANN_H()               (SQRT(0.5*kappaM1/kappa))
#define ROEC_RIEMANN_H(RoeH,RoeVel)    (SQRT(kappaM1*(RoeH-0.5*DOT_PRODUCT(RoeVel,RoeVel))))
#define ALPHA2_RIEMANN_H(RoeH,RoeVel,Roec,Delta_U)     (kappaM1/(Roec*Roec) * (Delta_U(1)*(RoeH-RoeVel(1)*RoeVel(1)) - Delta_U(6) + RoeVel(1)*Delta_U(2)))

! routines to compute physical quantities from conservative variables or extended variables
! conservative
#define VELOCITY_H(U,sRho)             (U(DMOMV)*sRho)
#define SPEEDOFSOUND_H(p,sRho)         (SQRT(Kappa*p*sRho))
#define TOTALENERGY_H(U,sRho,Vel)      (U(DENER)/U(DDENS))
#define TOTALENTHALPY_H(U,p,sRho)      ((U(DENER)+p)*sRho)
#define ENTROPY_H(U,T)                 (R*(sKappaM1*LOG(T)-LOG(U(DDENS))))
#define TEMPERATURE_H(U)               ((U(DENER)-0.5*DOT_PRODUCT(U(DMOMV),U(DMOMV))/U(DDENS))/(U(DDENS)*cv))

! extended (NOTE: compute from cons. When computing derived (neither prim or cons) variables
! assume that both prim and cons vars are filled
#define VELOCITY_HE(UE)                (UE(DMOMV)*UE(DSRHO))
#define PRESSURE_HE(UE)                (KappaM1*(UE(DENER)-0.5*DOT_PRODUCT(UE(DVELV),UE(DMOMV))))
#define SPEEDOFSOUND_HE(UE)            (SQRT(Kappa*UE(DPRES)*UE(DSRHO)))
#define TOTALENERGY_HE(UE)             (UE(DENER)*UE(DSRHO))
#define TOTALENTHALPY_HE(UE)           ((UE(DENER)+UE(DPRES))*UE(DSRHO))
#define TEMPERATURE_HE(UE)             (UE(DPRES)*UE(DSRHO)/R)
#define ENERGY_HE(UE)                  (sKappaM1*UE(DPRES)+0.5*DOT_PRODUCT(UE(DMOMV),UE(DVELV)))

#if PP_VISC == 0
#define VISCOSITY_PRIM(U)              mu0
#elif PP_VISC == 1
#define VISCOSITY_PRIM(U)              muSuth(U(6))
#elif PP_VISC == 2
#define VISCOSITY_PRIM(U)              mu0*U(6)**ExpoSuth
#endif

#if PP_VISC == 0
#define VISCOSITY_TEMPERATURE(T)       mu0
#elif PP_VISC == 1
#define VISCOSITY_TEMPERATURE(T)       muSuth(T)
#elif PP_VISC == 2
#define VISCOSITY_TEMPERATURE(T)       mu0*T**ExpoSuth
#endif
