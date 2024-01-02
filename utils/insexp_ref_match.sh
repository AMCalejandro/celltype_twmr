#!/usr/bin/env sh

# Take args from command line prompt
while getopts r:g:o: flag
do
    case "${flag}" in
        r) snps=${OPTARG};;
        g) gwaseqtl=${OPTARG};;
        o) out=${OPTARG};;
    esac
done

echo "Finding GWAS / eQTL significant SNPS that match the reference panel"
echo "Writing to: $out"
#echo $(head -n 1 $gwaseqtl) > $out
sed '1q' $gwaseqtl > $out
awk 'NR==FNR{A[$1];next}$1 in A' $snps $gwaseqtl >> $out

