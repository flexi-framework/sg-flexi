#!/bin/bash

# ====================================================================================================
# EXCEPTION FOR DYNAMIC POLYNOMIAL DEGREE (REST OF MAIN PROGRAM IS AT END OF FILE)
# ====================================================================================================

# get parameters
# ----------------------------------------------------------------------------------------------------
SrcDir=$1
TargetDir=$2
M=$3
nDim=$4
nCoef=$5

# exception
# ----------------------------------------------------------------------------------------------------
if [ "$M" = "M" ] ; then
  if [ -d "$TargetDir" ]; then
      rm -r $TargetDir
  fi
  mkdir $TargetDir
  cp $SrcDir/src/sg/*.f90 $TargetDir
  exit 0
fi

# ====================================================================================================
# FUNCTIONS
# ====================================================================================================

# append string to string array
# ----------------------------------------------------------------------------------------------------
Append()
{
local -n arr=$1
arr="$arr
$2"
}

# insert string array to file
# ----------------------------------------------------------------------------------------------------
InsertString()
{
local -n arr=$1
IFS=$'\n'
for line in ${arr[*]}; do
   if [ "$line" = "!" ]; then
      sed -i -e "${nLineStart}G" $TargetFile
   else
      sed -i -e "${nLineStart}i$line" $TargetFile
   fi
   nLineStart=$(($nLineStart+1))
done
}

# append  " res($1) += C($nC) * mat($2,$3) "  to ProdStr
# ----------------------------------------------------------------------------------------------------
AppendToStrings()
{
MapPQto1DInd $2 $3
if [ ${MapTMtoMat[$nTM]} -eq 0 ] ; then
   nMat=$(($nMat+1))
   MapTMtoMat[$nTM]=$nMat
   if [ $2 -eq $3 ] ; then
      Append MatStr "mat(${MapTMtoMat[$nTM]})  =  v1($2) \* v2($2)"
   else
      Append MatStr "mat(${MapTMtoMat[$nTM]})  =  v1($3) \* v2($2)  +  v1($2) \* v2($3)"
   fi
fi
if ${FirstOccurence[$1]} ; then
   Append ProdStr "SG_Product_${DimStr}Dim($1)  =                         $CStr  mat(${MapTMtoMat[$nTM]})"
   FirstOccurence[$1]=false
else
   Append ProdStr "SG_Product_${DimStr}Dim($1)  =  SG_Product_${DimStr}Dim($1)  +  $CStr  mat(${MapTMtoMat[$nTM]})"
fi
}

# get 1D index of triangular matrix (incl. diagonal) 
# ----------------------------------------------------------------------------------------------------
MapPQto1DInd()
{
nTM=0
for ((i=1; i <= $2 ; i++)); do
   nTM=$(($nTM+$i))
done
nTM=$(($nTM+$1))
}

# stoch 1D: append all lines of the form SG_Product(...) = ...
# ----------------------------------------------------------------------------------------------------
#GenerateOperationStrings_1Dim()
#{
## p=0 only has one entry (0,0,0), which is covered separately:
#Append ProdStr "SG_Product_${DimStr}Dim(0)  =  v1(0) \* v2(0)"
#for ((p=1; p <= $M ; p++)); do
   #for ((q=0; q <= $p ; q++)); do
      #for ((r=$(($p-$q)); r <= $q ; r+=2)); do
         #isOne=false; [ "$r" == 0 ] && isOne=true
         #PermutatePQR
      #done
   #done
#done
#}

#wrapper for all AppendToStrings funcs for one p/q/r combination
# ----------------------------------------------------------------------------------------------------
PermutatePQR()
{
if $isOne ; then
   CStr="          "
else
   nC=$(($nC+1))
   CStr="C_${DimStr}D($nC)  \*"
fi
if [  $p -gt $q -a $q -gt $r ] ; then
   AppendToStrings $r $q $p
   AppendToStrings $q $r $p
   AppendToStrings $p $r $q
elif [ $p -eq $q -a $q -gt $r ]; then
   AppendToStrings $p $r $p
   AppendToStrings $r $p $p
elif [ $p -gt $q -a $q -eq $r ]; then
   AppendToStrings $r $r $p
   AppendToStrings $p $r $r
elif [ $p -eq $q -a $q -eq $r ]; then
   AppendToStrings $p $p $p
fi
}

#sum entries of array "arr"
# ----------------------------------------------------------------------------------------------------
GetSum()
{
sum=0
for t in ${arr[@]}; do
   sum=$(($sum+$t))
done
}

# reset full order vector to [0,0, ... ,-1]
# ----------------------------------------------------------------------------------------------------
ResetVec()
{
local -n arr=$1
for ((t=1; t < $nDim ; t++)); do
   arr[$t]=0
done 
arr[$nDim]=-1
}

# get new full order vector from last one if 1d index is incremented
# ----------------------------------------------------------------------------------------------------
IncrementVec()
{
local -n arr=$1
GetSum arr
for ((tt=0; tt < $nDim ; tt++)); do
   t=$(($nDim-$tt))
   if [ $sum -eq $M ] ; then
      arr[$t]=0
      GetSum arr
   else
      arr[$t]=$((${arr[$t]}+1))
      break
   fi
done 
}

#check from pVec etc. whether C_pqr = 1 (is the case if C1D_pqr in every dimension is 1)
# ----------------------------------------------------------------------------------------------------
CheckIfOne()
{
isOne=true
for ((t=1; t <= $nDim ; t++)); do
   # if index of at least one dimension contains no zeros, then C_pqr != 1
   if [ ${pVec[$t]} -gt 0 -a ${qVec[$t]} -gt 0 -a ${rVec[$t]} -gt 0 ] ; then
      isOne=false
      break
   fi
done 
}

# stoch nD: append all lines of the form SG_Product(...) = ...
# ----------------------------------------------------------------------------------------------------
GenerateOperationStrings_nDim()
{
# initialize vecs
for ((t=1; t <= $nDim ; t++)); do
   pVec[$t]=0; qVec[$t]=0; rVec[$t]=0
done 
pVec[$nDim]=-1; qVec[$nDim]=-1; rVec[$nDim]=-1
# p=0 only has one entry (0,0,0), which is covered separately:
IncrementVec pVec
Append ProdStr "SG_Product_${DimStr}Dim(0)  =  v1(0) \* v2(0)"
for ((p=1; p < $nCoef ; p++)); do
   IncrementVec pVec
   ResetVec qVec
   for ((q=0; q <= $p ; q++)); do
      IncrementVec qVec
      ResetVec rVec
      for ((r=0; r <= $q ; r++)); do
         IncrementVec rVec
         
         # Check whether entry is non-zero
         NonZero=true
         for ((t=1; t <= $nDim ; t++)); do
            s=$((${pVec[$t]}+${qVec[$t]}+${rVec[$t]}))
            if [ $(( $s % 2)) -eq 1 ] ; then
               NonZero=false
               break
            fi
            s=$(($s/2))
            if [ ${pVec[$t]} -gt $s -o ${qVec[$t]} -gt $s -o  ${rVec[$t]} -gt $s ] ; then
               NonZero=false
               break
            fi
         done

         # Calculate nnz entry and index if element is non-zero
         if $NonZero ; then
            # check whether C_pqr = 1
            CheckIfOne
            PermutatePQR
         fi
      done
   done
done
}

# ====================================================================================================
# PROGRAM START
# ====================================================================================================

# initializations
# ----------------------------------------------------------------------------------------------------
nC=0          # index counter for C_nD
nMat=0        # index counter for mat
MatStr="!"    # string with all precomputations mat(.)=v1(.)*v2(.)+v1(.)*v2(.)
ProdStr="!"   # string with all main operations Prod += C * mat
Footer="!"
# initialize array for check whether to write Prod = mat (first occurence) or Prod += mat (else)
FirstOccurence[0]=false
for ((i=1; i < $nCoef ; i++)); do
   FirstOccurence[$i]=true
done
# initialize mapping from triangular matrix to only needed matrix entries
MapPQto1DInd $(($nCoef-1)) $(($nCoef-1)) # get size of triangular matrix
for ((t=1; t <= $nTM ; t++)); do
   MapTMtoMat[$t]=0
done 
# names in files are either 1Dim or nDim
if [ $nDim -eq 1 ] ; then
  DimStr="1"
else
  DimStr="n"
fi
FunctionStartString="PURE FUNCTION SG_Product_${DimStr}Dim(v1,v2)"
FunctionEndString="END FUNCTION SG_Product_${DimStr}Dim"

# write function header
# ----------------------------------------------------------------------------------------------------
Header="$FunctionStartString
! MODULES
USE MOD_Preproc
USE MOD_SG_Vars , ONLY : C_${DimStr}D
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
REAL,INTENT(IN),DIMENSION(0:PP_nCoefM1) :: v1,v2
REAL,DIMENSION(0:PP_nCoefM1)            :: SG_Product_${DimStr}Dim
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES"

# write actual sg prod steps
# ----------------------------------------------------------------------------------------------------
GenerateOperationStrings_nDim

# complete strings
# ----------------------------------------------------------------------------------------------------
if [ $nMat -gt 0 ] ; then
   Append Header "REAL,DIMENSION($nMat)  :: mat"
fi
Append Header "!=================================================================================================================================="
Append MatStr "!"
Append Footer "$FunctionEndString"

# copy source files and modify them 
# ----------------------------------------------------------------------------------------------------

if [ -d "$TargetDir" ]; then
    rm -r $TargetDir
fi
mkdir $TargetDir
# copy everything into preproc TargetDir (which is where the compiler looks for source files)
cp $SrcDir/src/sg/*.f90 $TargetDir
# get line number of start and end of function
TargetFile="$TargetDir/sg_operators.f90"
nLineStart=`grep -n "$FunctionStartString" $TargetFile |cut -f1 -d:`
nLineEnd=`grep -n "$FunctionEndString" $TargetFile |cut -f1 -d:`
# cut function SG_Product_1Dim from file 
sed -i -e "$nLineStart"\,"$nLineEnd"d $TargetFile
# insert strings where function SG_Product_1Dim was
InsertString Header
InsertString MatStr
InsertString ProdStr
InsertString Footer

# ====================================================================================================
