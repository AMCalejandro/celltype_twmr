#!/usr/bin/env bash


#inputs=~/celltyping_wd/celltype_twmr/test_epqtl/data/gwasqtl_qc/
#abswd=~/celltyping_wd/celltype_twmr/test_epqtl/workdir/
#inputs=/home/rstudio/celltyping_wd/celltype_twmr/test_epqtl/data/tmp/sbatch_test/
inputs=/mnt/acarrasco/RDS/acarrasco/ANALYSES_WORKSPACE/CELLTYPING_PROJECT/CELLTYPING_V2/celltype_twmr/test_epqtl/data/tmp/sbatch_test/
abswd=/mnt/acarrasco/RDS/acarrasco/ANALYSES_WORKSPACE/CELLTYPING_PROJECT/CELLTYPING_V2/celltype_twmr/test_epqtl/workdir/
files=$(ls $inputs)

for file in $files
do
eqtlname=$(echo $file | cut -d_ -f1)
gwasname=$(echo $file | cut -d_ -f2)
out="${eqtlname}_${gwasname}"
wd="${abswd}${eqtlname}_${gwasname}"
echo $inputs$file
echo $wd
echo $eqtlname
echo -e "\n"
sbatch -t 24:00:00 -n 1 -o log_${out}.txt --wrap="Rscript ~/celltyping_wd/celltype_twmr/test_epqtl/utils/prepare_twmrinput.R --args $wd $inputs$file $eqtlname"
done
