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
#include "flexi.h"
#include "eos.h"

!==================================================================================================================================
!==================================================================================================================================
MODULE MOD_EOS_Posti_Vars
! MODULES
USE MOD_PreProc
IMPLICIT NONE
PUBLIC
SAVE

#ifdef SGdummy
#if PARABOLIC
#define CUT(x)
INTEGER,PARAMETER :: nVarDepEOS=36
#else
#define CUT(x) x!
INTEGER,PARAMETER :: nVarDepEOS=19
#endif
#else/*SGdummy*/
#define CUT(x) x!
INTEGER,PARAMETER :: nVarDepEOS=19
#endif /*SGdummy*/

! ATTENTION: The first     5 variables must be the conservative ones
!            The following 5 variables must be the primitive ones
!           E
!           n
!           e                                                               W
!           r                                                               a
!           g                                                               l
!           y                    E                        V N               l
!           S            V       n       P                o o               F
!           t            e     E t   T   r                r r               r W
!           a            l     n h   o   e                t m               i a
! W         g            o     e a   t   s                i a               c l
! i         n            c V   r l   a T s                c l         W W W t l
! t         a            i e   g p   l o u                i i         a a a i H
! h         t         T  t l   y y   T t r                t z         l l l o e
! D         i         e  y o   S S   e a e          V V V y e   D Q   l l l n a
! G   M M M o V V V   m  M c   t t   m l T          o o o M d   i C S F F F M t
! O   o o o n e e e P p  a i   a a   p P i          r r r a H   l r c r r r a T
! p D m m m D l l l r e  g t   g g E e r m          t t t g e L a i h i i i g r
! e e e e e e o o o e r  n y   n n n r e e          i i i n l a t t l c c c n a
! r n n n n n c c c s a  i S   a a t a s D          c c c i i m a e i t t t i n
! a s t t t s i i i s t  t o M t t r t s e          i i i t c b t r e i i i t s
! t i u u u i t t t u u  u u a i i o u u r          t t t u i d i i r o o o u f
! o t m m m t y y y r r  d n c o o p r r i          y y y d t a o o e n n n d e x y z
! r y X Y Z y X Y Z e e  e d h n n y e e v          X Y Z e y 2 n n n X Y Z e r + + +
INTEGER,DIMENSION(1:nVarDepEOS,0:nVarDepEOS),PARAMETER :: DepTableEOS = TRANSPOSE(RESHAPE(&
(/&
  0,1,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0, CUT(&)  0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 ,& !1  Density
  0,0,1,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0, CUT(&)  0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 ,& !2  MomentumX
  0,0,0,1,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0, CUT(&)  0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 ,& !3  MomentumY
  0,0,0,0,1,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0, CUT(&)  0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 ,& !4  MomentumZ
  0,0,0,0,0,1,0,0,0,0,0, 0,0,0,0,0,0,0,0,0, CUT(&)  0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 ,& !5  EnergyStagnationDensity
  0,1,1,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0, CUT(&)  0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 ,& !6  VelocityX
  0,1,0,1,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0, CUT(&)  0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 ,& !7  VelocityY
  0,1,0,0,1,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0, CUT(&)  0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 ,& !8  VelocityZ
  0,1,1,1,1,1,0,0,0,0,0, 0,0,0,0,0,0,0,0,0, CUT(&)  0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 ,& !9  Pressure
  0,1,0,0,0,0,0,0,0,1,0, 0,0,0,0,0,0,0,0,0, CUT(&)  0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 ,& !10 Temperature
  0,1,1,1,1,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0, CUT(&)  0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 ,& !11 VelocityMagnitude
  0,1,0,0,0,0,0,0,0,1,0, 0,0,0,0,0,0,0,0,0, CUT(&)  0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 ,& !12 VelocitySound
  0,0,0,0,0,0,0,0,0,0,0, 1,1,0,0,0,0,0,0,0, CUT(&)  0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 ,& !13 Mach
  0,1,0,0,0,1,0,0,0,0,0, 0,0,0,0,0,0,0,0,0, CUT(&)  0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 ,& !14 EnergyStagnation
  0,1,0,0,0,0,0,0,0,1,0, 0,0,0,1,0,0,0,0,0, CUT(&)  0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 ,& !15 EnthalpyStagnation
  0,1,0,0,0,0,0,0,0,0,1, 0,0,0,0,0,0,0,0,0, CUT(&)  0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 ,& !16 Entropy
  0,0,0,0,0,0,0,0,0,0,1, 1,0,0,0,0,0,0,0,0, CUT(&)  0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 ,& !17 TotalTemperature
  0,1,0,0,0,0,0,0,0,1,0, 1,0,0,0,0,0,0,0,0, CUT(&)  0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 ,& !18 TotalPressure
  1,1,1,1,1,1,0,0,0,0,0, 0,0,0,0,0,0,0,0,0  CUT(&) ,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 ,& !19 PressureTimeDeriv
#ifdef SGdummy
#if PARABOLIC
  1,0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0         ,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 ,& !20 VorticityX
  1,0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0         ,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 ,& !21 VorticityY
  1,0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0         ,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 ,& !22 VorticityZ
  1,0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0         ,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 ,& !23 VorticityMagnitude
  1,1,1,1,1,0,0,0,0,0,0, 1,0,0,0,0,0,0,0,0         ,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0 ,& !24 NormalizedHelicity
  1,0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0         ,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 ,& !25 Lambda2
  1,0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0         ,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 ,& !26 Dilatation
  1,0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0         ,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 ,& !27 QCriterion
  1,0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0         ,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 ,& !28 Schlieren
  1,0,0,0,0,0,0,0,0,0,1, 0,0,0,0,0,0,0,0,0         ,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 ,& !29 WallFrictionX
  1,0,0,0,0,0,0,0,0,0,1, 0,0,0,0,0,0,0,0,0         ,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 ,& !30 WallFrictionY
  1,0,0,0,0,0,0,0,0,0,1, 0,0,0,0,0,0,0,0,0         ,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 ,& !31 WallFrictionZ
  1,0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0         ,0,0,0,0,0,0,0,0,0,1,1,1,0,0,0,0,0 ,& !32 WallFrictionMagnitude
  1,0,0,0,0,0,0,0,0,0,1, 0,0,0,0,0,0,0,0,0         ,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 ,& !33 WallHeatTransfer
  1,0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0         ,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0 ,& !34 x+
  1,0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0         ,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0 ,& !35 y+
  1,0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0         ,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0  & !36 z+
#endif
#endif /*SGdummy*/
/),(/nVarDepEOS+1,nVarDepEOS/)))

! Mark all quantities that can be calculated exclusively on the surface
INTEGER,DIMENSION(1:nVarDepEOS),PARAMETER :: DepSurfaceOnlyEOS = &
(/  0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0  CUT(&) ,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1 &
/)

! Mark all quantities that can be calculated exclusively in the volume and must be prolonged to the surface from the volume
INTEGER,DIMENSION(1:nVarDepEOS),PARAMETER :: DepVolumeOnlyEOS = &
(/  0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,1  CUT(&) ,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 &
/)

#if FV_ENABLED && FV_RECONSTRUCT
!           E
!           n
!           e
!           r
!           g
!           y
!           S
!           t
!           a
! W         g
! i         n
! t         a
! h         t         T
! D         i         e
! G   M M M o V V V   m
! O   o o o n e e e P p
! p D m m m D l l l r e
! e e e e e e o o o e r
! r n n n n n c c c s a
! a s t t t s i i i s t
! t i u u u i t t t u u
! o t m m m t y y y r r
! r y X Y Z y X Y Z e e
INTEGER,DIMENSION(1:PP_nVarDet,0:nVarDepEOS),PARAMETER :: DepTablePrimToCons =TRANSPOSE(RESHAPE(&
(/&
  0,0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0, CUT(&)  0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 ,& !1 Density
  0,1,0,0,0,0,1,0,0,0,0, 0,0,0,0,0,0,0,0,0, CUT(&)  0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 ,& !2 MomentumX
  0,1,0,0,0,0,0,1,0,0,0, 0,0,0,0,0,0,0,0,0, CUT(&)  0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 ,& !3 MomentumY
  0,1,0,0,0,0,0,0,1,0,0, 0,0,0,0,0,0,0,0,0, CUT(&)  0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 ,& !4 MomentumZ
  0,1,0,0,0,0,1,1,1,1,0, 0,0,0,0,0,0,0,0,0  CUT(&) ,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0  & !5 EnergyStagnationDensity
/),(/nVarDepEOS+1,PP_nVarDet/)))
#endif
#undef CUT

CHARACTER(LEN=255),DIMENSION(nVarDepEOS),PARAMETER :: DepNames = &
(/ CHARACTER(LEN=255) ::    &
"Density"                  ,& !1
"MomentumX"                ,& !2
"MomentumY"                ,& !3
"MomentumZ"                ,& !4
"EnergyStagnationDensity"  ,& !5
"VelocityX"                ,& !6
"VelocityY"                ,& !7
"VelocityZ"                ,& !8
"Pressure"                 ,& !9
"Temperature"              ,& !10
"VelocityMagnitude"        ,& !11
"VelocitySound"            ,& !12
"Mach"                     ,& !13
"EnergyStagnation"         ,& !14
"EnthalpyStagnation"       ,& !15
"Entropy"                  ,& !16
"TotalTemperature"         ,& !17
"TotalPressure"            ,& !18
"PressureTimeDeriv"         & !19
#ifdef SGdummy
#if PARABOLIC
,"VorticityX"              ,& !20
"VorticityY"               ,& !21
"VorticityZ"               ,& !22
"VorticityMagnitude"       ,& !23
"NormalizedHelicity"       ,& !24
"Lambda2"                  ,& !25
"Dilatation"               ,& !26
"QCriterion"               ,& !27
"Schlieren"                ,& !28
"WallFrictionX"            ,& !29
"WallFrictionY"            ,& !30
"WallFrictionZ"            ,& !31
"WallFrictionMagnitude"    ,& !32
"WallHeatTransfer"         ,& !33
"x+"                       ,& !34
"y+"                       ,& !35
"z+"                        & !36
#endif /*PARABOLIC*/
#endif /*SGdummy*/
/)

END MODULE MOD_EOS_Posti_Vars
