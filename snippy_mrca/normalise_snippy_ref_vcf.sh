#!/bin/bash

set -e

# This bash script is used to normalise variants produced by Snippy 

vcf_in=$1
vcf_out=$2
ref_fasta=$3


### File names

vcf_tmp1=`echo $vcf_in | sed 's/.vcf.gz$/.sorted.vcf/g'`;


### Commands

# Sorting and indexing VCF

bcftools-1.9 sort -o $vcf_tmp1 -O v $vcf_in

bgzip-1.9 $vcf_tmp1

tabix-1.9 -p vcf $vcf_tmp1".gz"


# Forth, normalising variants

bcftools-1.9 norm -f $ref_fasta $vcf_tmp1".gz" > $vcf_out

# Finally, removed temporary files

rm $vcf_tmp1".gz" $vcf_tmp1".gz.tbi"