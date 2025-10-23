#!/bin/bash

set -e

vcf_dir=$1;
vcf_suffix=$2;
sample=$3;
out_dir=$4;
out_suffix=$5;

vcf_file=$vcf_dir$sample$vcf_suffix;
tmp1=$out_dir$sample".tmp1.txt";
tmp2=$out_dir$sample".tmp2.txt";
tmp3=$out_dir$sample".tmp3.txt";
het_pos_file=$out_dir$sample$out_suffix;

if [ -f $vcf_file ]
then

if [ ! -f $het_pos_file ]
then

bcftools-1.9 view $vcf_file | awk -F'\t' '{ print $8}' | awk -F';' '{for (i=1;i<=NF;i++) {if ($i ~ /^AF=/) {print $i}}}' | sed 's/AF=//g' > $tmp1
bcftools-1.9 view $vcf_file | awk -F'\t' '{ print $2}' | grep -v -e '^$' | grep -v "^POS" > $tmp2
paste -d'\t' $tmp2 $tmp1 > $tmp3
cat $tmp3 | awk -F'\t' '{ if($2>0.2 && $2<0.8) print $1}' > $het_pos_file
rm $tmp1 $tmp2 $tmp3

else
	echo $het_pos_file" already found" 
fi
else
	echo $vcf_file" could not be found" 
fi