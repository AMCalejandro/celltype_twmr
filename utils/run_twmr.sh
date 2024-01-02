#!/usr/bin/env bash


#WORKDIR="$@"
WORKDIR=$1
NGWAS=$2
NEQTL=$3
CHRFILES=$WORKDIR/chr*

#CHRFILES=/home/rstudio/celltyping_wd/celltype_twmr/test_epqtl/workdir/eqtlgen_height/chr2
# echo -e "$NGWAS \n"
# echo -e "$NEQTL \n"
# echo -e "$CHRFILES \n"

for f in $CHRFILES
do

#chr=${f#*/chr}
#echo -e "Working on chr$chr... \n\n"
allgenes=$(ls $f/)

for gene in $allgenes
do
#echo -e "$gene \n"
if test -f "$f/$gene/$gene.matrix"; then
  ~/celltyping_wd/celltype_twmr/test_epqtl/utils/MR.R --args $gene $f/$gene/$gene.matrix $f/$gene/$gene.ld $NGWAS $NEQTL
else
  echo "$f/$gene/$gene.matrix does not exist."
fi
done

done

echo -e "\n Finished \n"



# TESTS

## Tests in 1 given chromosome only
# f=/home/rstudio/celltyping_wd/celltype_twmr/test_epqtl/workdir/eqtlgen_pdrisk/chr1
# 
# allgenes=$(ls $f/)
# 
# for gene in $allgenes
# do
# #echo -e "$gene \n"
# if test -f "$f/$gene/$gene.matrix"; then
#   ~/celltyping_wd/celltype_twmr/test_epqtl/utils/MR.R --args $gene $f/$gene/$gene.matrix $f/$gene/$gene.ld
# else 
#   echo "$f/$gene/$gene.matrix does not exist."
# fi
# done