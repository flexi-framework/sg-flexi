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

!==================================================================================================================================
!> Contains the different Surface integral formulations
!> Computes the Surface integral for all faces using U and updates Ut
!> Computes only inner surface integrals!
!> Surface integrals are separated for each direction
!==================================================================================================================================
MODULE MOD_SurfInt
IMPLICIT NONE
PRIVATE

#define WITHnVar 1

INTERFACE SurfInt
  MODULE PROCEDURE SurfInt
END INTERFACE

INTERFACE DoSurfInt
  MODULE PROCEDURE DoSurfInt
END INTERFACE

PUBLIC::SurfInt,DoSurfInt

CONTAINS
#include "surfint.t90"
END MODULE MOD_SurfInt

!==================================================================================================================================

MODULE MOD_SurfIntCons
IMPLICIT NONE
PRIVATE

#ifndef PP_M
#define WITHnVar 1
#else
#undef WITHnVar
INTEGER,PARAMETER :: TP_nVar = PP_nVar
#endif

INTERFACE SurfIntCons
  MODULE PROCEDURE SurfInt
END INTERFACE

INTERFACE DoSurfIntCons
  MODULE PROCEDURE DoSurfInt
END INTERFACE

PUBLIC::SurfIntCons,DoSurfIntCons

CONTAINS
#include "surfint.t90"
END MODULE MOD_SurfIntCons

!==================================================================================================================================

MODULE MOD_SurfIntPrim
IMPLICIT NONE
PRIVATE

#ifndef PP_M
#define WITHnVar 1
#else
#undef WITHnVar
INTEGER,PARAMETER :: TP_nVar = PP_nVarPrim
#endif

INTERFACE SurfIntPrim
  MODULE PROCEDURE SurfInt
END INTERFACE

INTERFACE DoSurfIntPrim
  MODULE PROCEDURE DoSurfInt
END INTERFACE

PUBLIC::SurfIntPrim,DoSurfIntPrim

CONTAINS
#include "surfint.t90"
END MODULE MOD_SurfIntPrim
