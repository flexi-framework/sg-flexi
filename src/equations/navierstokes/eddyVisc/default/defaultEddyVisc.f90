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
!> Subroutines needed by the eddy viscosity model. This is the default  model, which means no eddy viscosity.
!> We need to initialize some arrays that are used in the interfaces to the rounties. The eddy viscosity model itself will
!> return 0 as turbulent viscosity.
!==================================================================================================================================
MODULE MOD_DefaultEddyVisc
! MODULES
IMPLICIT NONE
PRIVATE

PUBLIC::DefaultEddyVisc,FinalizeDefaultEddyViscosity
!===================================================================================================================================

CONTAINS

!===================================================================================================================================
!> Dummy for default eddy viscosity (meaning no eddy viscosity), do nothing since the muSGS arrays will be passed here and they
!> are zero all the time.
!===================================================================================================================================
SUBROUTINE DefaultEddyVisc(iElem,i,j,k,muSGS)
! MODULES
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER,INTENT(IN)                        :: iElem             !< index of current element
!> indices of the current volume point
INTEGER,INTENT(IN)                        :: i,j,k
!> gradients of the velocities w.r.t. all directions
REAL,INTENT(INOUT)                        :: muSGS             !< local SGS viscosity
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
#ifdef DEBUG
! ===============================================================================
! Following dummy calls do suppress compiler warnings of unused Riemann-functions
! ===============================================================================
IF (0.EQ.1) THEN
  muSGS = i+j+k+iElem 
END IF
#endif
END SUBROUTINE DefaultEddyVisc

!===============================================================================================================================
!> Deallocate arrays and finalize variables used by Smagorinsky SGS model
!===============================================================================================================================
SUBROUTINE FinalizeDefaultEddyviscosity()
! MODULES
USE MOD_EddyVisc_Vars
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!===============================================================================================================================
END SUBROUTINE FinalizeDefaultEddyViscosity

END MODULE MOD_DefaultEddyVisc
