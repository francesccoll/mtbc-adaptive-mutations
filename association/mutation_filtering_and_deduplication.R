
# FC: note that patient_id in this script corresponds to strain_id in efm_wh_extended.v1 dataset
library(gtsummary)

install.packages("splitstackshape")
install.packages("stringr")
install.packages("dplyr")
library(splitstackshape)
library(stringr)
library(dplyr)
library(tidyr)

###############################################################################################################
####                                         INPUT/OUTPUT FILES                                            ####
###############################################################################################################

dataset = "mtb_wh_extended.v1";

# Change ref and ref_file variables to change reference genome
ref = "h37rv"; ref_file = "H37Rv";

setwd("./");
dir = "./"
combined_table = paste(dir, "data/",dataset,".mrca.snippy_snpEff_annotated_variants.",ref,".in_ref_VCF.kept.txt.gz", sep = ""); 
gff_table_file = paste(dir,"annotation/Mycobacterium_tuberculosis_",ref_file,"_gff_v4.gff.table.txt", sep = "");
gff_table_updated = paste(dir,"annotation/Mycobacterium_tuberculosis_",ref_file,"_txt_v4.txt", sep = "");
low_complexity_file = paste(dir,"annotation/Mycobacterium_tuberculosis_ancestral_reference",".dustmasker.chr_pos.txt.gz", sep = "");
repetitive_regions_file = paste(dir,"annotation/Mycobacterium_tuberculosis_ancestral_reference",".repetitive_regions.tab", sep = "");
holt2018_marin2022_file = paste(dir, "annotation/holt2018_and_marin2022.masked_regions.tab", sep = "");

filtered_mutations_table = paste(dir,dataset,".mrca.snippy_snpEff_annotated_variants.",ref,".in_ref_VCF.filtered.txt", sep = "");

###############################################################################################################
####                                              MAIN CODE                                                ####
###############################################################################################################

mutations = read.delim(gzfile(combined_table), sep = "\t")
dim(mutations)
# [1] 43584    31 


###### Mutations at regions of low complexity need to be removed
lc_table = read.delim(gzfile(low_complexity_file) , sep = "\t", header = T)
dim(lc_table)
# [1] 4411531       3


low_complexity = lc_table$low_complexity[mutations$pos]
mutations = cbind(mutations, low_complexity)

mutations = subset(mutations, low_complexity == "no")
dim(mutations)
# [1] 41419    32


###### Mutations at repetitive regions need to be removed
rep_table = read.table(repetitive_regions_file, sep = "\t", header = F)
rep_table = subset(rep_table, select = -c(V2))
dim(rep_table)

# [1] 1078   2
rep_chr_pos = vector()
for(r in 1:nrow(rep_table))
{
  rep_chr_pos = c(rep_chr_pos, seq(rep_table[r,1], rep_table[r,2], 1))
}
length(unique(rep_chr_pos))
# [1] 607478 

tmp = which(!is.na(match(mutations$pos, rep_chr_pos)))
if(length(tmp)>0){ mutations = mutations[-tmp,]; }
dim(mutations)
# [1] 36601    32

#### additional mutations to be removed 
holt2018_and_marin2022 = read.table(holt2018_marin2022_file, sep = "\t", header = T)

# write.table(holt2018_and_marin2022, file= "/Users/helenzhang/Desktop/LSHTM/thesis/association/holt2018_and_marin2022.csv", sep = '\t', col.names = T, row.names = F, quote = F)
holt2018_and_marin2022 = read.csv("/Users/helenzhang/Desktop/LSHTM/thesis/association/holt2018_and_marin2022.csv", sep = ",")
holt2018_and_marin2022 =separate(data = holt2018_and_marin2022,col = region, into = c("start", "end"), sep = ",")
dim(holt2018_and_marin2022)
# [1] 957   4

masked_chr_pos = vector()
for(x in 1:nrow(holt2018_and_marin2022))
{
  masked_chr_pos = c(masked_chr_pos, seq(holt2018_and_marin2022[x,3], holt2018_and_marin2022[x,4], 1))
}
length(unique(masked_chr_pos))
# [1] 467531

tmp1 = which(!is.na(match(mutations$pos, masked_chr_pos)))
if(length(tmp1)>0){ mutations = mutations[-tmp1,]; }
dim(mutations)
# [1] 35305    32


#### Number of unique samples (i.e. pairs) with mutatations
laSamples = unique(as.vector(mutations$sample))
length(laSamples)
# [1] 2975 

# mutations found in multiple isolates of the same patient need to be removed (duplicated)

# NOTE: the function get_patient_id actually get the strain_id (the last part of sample ID)
get_patient_id = function(x)
{ 
  tmp = unlist(strsplit(as.character(x),"[.]"))
  y = paste(tmp[-(1:3)],collapse=".");
  if(grepl("^GCA_",tmp[2])){ y = paste(tmp[-(1:4)],collapse="."); }
  return(y);
}

mutations$patient_id = sapply(mutations$sample, get_patient_id)
patients = unique(as.vector(mutations$patient_id))
length(patients)
# [1] 503 

# NOTE: here the vector 'patients' actually contains strain ids > to do: make sure all these strain ids are also found in metadata file efm_wh_extended.v1.metadata.mi.qc.csv

mutations_to_remove = vector();
for(p in 1:length(patients))
{
  tmp = which(mutations$patient_id == patients[p])
  tmp2 = which(duplicated(mutations$pos[tmp]))
  if(length(tmp2)>0){ mutations_to_remove = c(mutations_to_remove, tmp[tmp2]); }
}
mutations = mutations[-mutations_to_remove,];
dim(mutations)
# [1] 3296   32 

# NOTE: a total of 3296 mutations can be attributed to point mutations (as opposed to recombination)


write.table(mutations, file=filtered_mutations_table, sep = '\t', col.names = T, row.names = F, quote = F)






