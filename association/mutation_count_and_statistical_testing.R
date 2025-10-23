
library(dplyr)
library(tidyr)
###############################################################################################################
####                                         INPUT/OUTPUT FILES                                            ####
###############################################################################################################

dataset = "mtb_wh_extended.v1."

# Change ref and ref_file variables to change reference genome
ref = "h37rv.CDS_only"; ref_file = "H37Rv.CDS_only";
# ref = "h37rv.TU_CDS"; ref_file = "H37Rv.TU_CDS";
# ref = "h37rv.TU_promoter"; ref_file = "H37Rv.TU_promoter";

# set working directory
dir = "./"

# Mutations table
filtered_mutations_table = paste(dir,"data/","mtb_wh_extended.v1.mrca.snippy_snpEff_annotated_variants.h37rv.in_ref_VCF.filtered.txt", sep = "");

# Metadata table
metadata_table = paste(dir,"data/","mtb_wh_extended.v1.basic_metadata.mi.qc.csv", sep = "");
tb_profiler_file = paste(dir,"data/","mtb_wh_extended.v1.tb-profiler.strain.csv", sep = "");

# Annotation files
gff_txt_file = paste(dir,"annotation/Mycobacterium_tuberculosis_H37Rv_txt_v4.txt", sep = "");
all_transcription_unit_file = paste(dir, "annotation/All-transcription-units-of-M.-tuberculosis-H37Rv.txt", sep = "")
all_genes_file = paste(dir, "annotation/All-instances-of-Genes-in-Mycobacterium-tuberculosis-H37Rv.txt", sep = "")

# Change unit of association
feature_type = "CDS_only_functional"; # only mutations with HIGH or MODERATE functional annotation at CDS are included
# feature_type = "TU_CDS_only_functional"; # only mutations with HIGH or MODERATE functional annotation at CDS are included
# feature_type = "TU_promoter"; # putative TU promoter region 100 bp upstream

# Other parameteres
min_UA_length = 300; # minimum length of unit of association (in base pairs)
mutations_count_table = paste(dataset,".",feature_type,".",ref,".in_ref_VCF.mutations_count.txt",sep=""); 
mutations_count_table_assoc = paste(dataset,".",feature_type,".",ref,".in_ref_VCF.mutations_count.assoc.txt",sep="");


###############################################################################################################
####                                              FUNCTIONS                                                ####
###############################################################################################################

# Taken from https://github.com/johnlees/paired-samples/blob/master/poisson_test.R

poisson_test <- function(mutations, length, total_mutations, genome_length)
{
  print(length)
  length <- as.numeric(length)
  mutations <- as.numeric(mutations)
  p.value<-poisson.test(mutations,r=total_mutations*(length/genome_length), alternative = "greater")['p.value']
  return(as.numeric(p.value))
}


###############################################################################################################
####                                              FEATURE LINKING                                          ####
###############################################################################################################
mutations = read.delim(gzfile(filtered_mutations_table), sep = "\t")
dim(mutations)
# [1] 3296   32

# #load and edit TU annotation file 
all_transcription_unit = read.delim(all_transcription_unit_file, sep = "\t", header = T)
all_transcription_unit <- separate(all_transcription_unit, Sequence...coordinates.of.DNA.region, c("start", "end"), sep = ";")
all_transcription_unit$start <- gsub("\\[TU185E-*", "", all_transcription_unit$start)
all_transcription_unit$start <- gsub("\\&rarr", "", all_transcription_unit$start)
all_transcription_unit$end <- gsub("\\]", "", all_transcription_unit$end)
all_transcription_unit <- separate(all_transcription_unit, start, c("name", "start"), sep = " ")
all_transcription_unit <- subset(all_transcription_unit, select = c("Transcription.Units", "start", "end", "Genes.of.transcription.unit", "Direction"))
all_transcription_unit$start <- gsub(",", "", all_transcription_unit$start)
all_transcription_unit$end <- gsub(",", "", all_transcription_unit$end)
all_transcription_unit$start <- as.numeric(all_transcription_unit$start)
all_transcription_unit$end <- as.numeric(all_transcription_unit$end)

dim(all_transcription_unit)
# [1] 2573    5

# ##add TU ID to mutations table 
for (x in mutations$pos) {
  for (y in 1:nrow(all_transcription_unit)) {
    if (x >= all_transcription_unit[y,2] && x <= all_transcription_unit[y,3]) {
      mutations$id[which(mutations$pos == x)] = all_transcription_unit[y,1]
    }
  }
}
mutations = rename(mutations, TU_id = id)

#all genes file 
all_genes = read.delim(all_genes_file, sep = "\t", header = T)
dim(all_genes)
# [1] 4072    7

## make a TU table (similar to gene_table)
TU_table = subset(all_transcription_unit, select = c("Transcription.Units", "start", "end", "Direction"))
colnames(TU_table) = c("TU_id", "start", "end", "Direction")

#add all genes under a TU
TU_table$TU_genes = "-"
for (i in 1:nrow(TU_table)) {
  id = unlist(strsplit(as.character(TU_table$TU_id[i]) , split = "//"))
  tmp = match(id, all_transcription_unit$Transcription.Units)
  if (length(tmp)>0)
  {
    TU_table$TU_genes[i] = as.character(paste(unique(unlist(strsplit(as.character(all_transcription_unit$Genes.of.transcription.unit[tmp]) , split = "//"))),collapse = ";"))
    }
}

## make a TU promoter table 
TU_promoter_table = subset(TU_table, select = c("TU_id","start", "end", "Direction"))
TU_promoter_table$promoter_start = "-"
TU_promoter_table$promoter_end = "-"

for (p in 1:nrow(TU_promoter_table)) {
  direction = TU_promoter_table$Direction[p]
  if (direction == "+") 
  {
    TU_promoter_table$promoter_end[p] = TU_promoter_table$start[p]
    TU_promoter_table$promoter_start[p] = TU_promoter_table$start[p] - 100
  }
}
for (p in 1:nrow(TU_promoter_table)) {
  direction = TU_promoter_table$Direction[p]
  if (direction == "-") 
  {
    TU_promoter_table$promoter_start[p] = TU_promoter_table$end[p]
    TU_promoter_table$promoter_end[p] = TU_promoter_table$end[p] + 100
  }
}
TU_promoter_table$promoter_id = paste(TU_promoter_table$TU_id, "promoter", sep ="_")
TU_promoter_table = subset(TU_promoter_table, select = c("promoter_id", "Direction", "promoter_start", "promoter_end"))
TU_promoter_table = TU_promoter_table[which(!duplicated(TU_promoter_table$promoter_id)),]
TU_promoter_table = TU_promoter_table[which(TU_promoter_table$promoter_start != -99),]
TU_promoter_table = arrange(TU_promoter_table, promoter_id)
TU_promoter_table$promoter_start = as.numeric(TU_promoter_table$promoter_start)
TU_promoter_table$promoter_end = as.numeric(TU_promoter_table$promoter_end)
  
#add TU genes to mutations table 
mutations$gene = "-"
for (i in 1:nrow(mutations)) {
  id = unlist(strsplit(as.character(mutations$TU_id[i]) , split = "//"))
  id = gsub(" ", "", id)
  tmp = match(id, TU_table$TU_id)
  if (length(tmp)>0)
  {
    mutations$gene[i] = as.character(paste(unique(unlist(strsplit(as.character(TU_table$TU_genes[tmp]) , split = ";"))),collapse = ";"))
  }
}

#add promoter ID to mutations table 
mutations$TU_promoter_id = "-"
for (x in mutations$pos) {
  for (y in 1:nrow(TU_promoter_table)) {
    if (x >= TU_promoter_table[y,3] && x <= TU_promoter_table[y,4]) {
      mutations$TU_promoter_id[which(mutations$pos == x)] = TU_promoter_table[y,1]
    }
  }
}

# ##add locus tag to mutations table 
gene_table = read.delim(gff_txt_file, sep = "\t", header = T)
dim(gene_table)
# [1] 4031    8
# [1] 4187   34 --> updated table 

# only keep the columns needed from the gene table 
gene_table_slim <- subset(gene_table, select = c(Locus, Name, Product, Feature, Start, Stop, Strand))

mutations$locus_tag = "-"
mutations$locus_name = "-"
mutations$feature = "-"
mutations$product = "-"
mutations$feature = "-"
for (x in mutations$pos) {
  for (y in 1:nrow(gene_table_slim)) {
    if (x >= gene_table_slim[y,5] && x <= gene_table_slim[y,6]) {
      mutations$locus_tag[which(mutations$pos == x)] = gene_table_slim[y,1]
      mutations$locus_name[which(mutations$pos == x)] = gene_table_slim[y,2]
      mutations$feature[which(mutations$pos == x)] = gene_table_slim[y,4]
      mutations$product[which(mutations$pos == x)] = gene_table_slim[y,3]
    }
  }
}

# save a version of the mutation table with added columns 
write.table(mutations, paste(dir, "filtered_mutations_added_columns.csv", sep = ""), sep = '\t', col.names = T, row.names = F, quote = F)

# upload the saved version
mutations = read.delim(paste(dir, "filtered_mutations_added_columns.csv", sep = ""), sep = "\t", header = T, quote = "")
dim(mutations)
# [1] 3296   35

###############################################################################################################
####                                              MAIN CODE                                                ####
###############################################################################################################


metadata = read.csv(metadata_table, sep = "\t", header = T)
dim(metadata)
# [1] 11995    10
tb_profiler = read.csv(tb_profiler_file, sep = "\t")
metadata = merge(metadata, tb_profiler, by.x = 'isolate_id', by.y = 'sample')
dim(metadata)
# [1] 11995    12

#### If feature_type contains "CDS_only", then keep only mutations on CDS (intergenic are removed)
if(grepl("CDS_only",feature_type))
{
  mutations = subset(mutations, feature_type=="transcript")
}
dim(mutations)
# [1] 2907   35 > CDS_only
# [1] 2907   35 > TU_CDS_only  
# [1] 2907   35 > submodule_CDS_only_functional


if(grepl("_promoter",feature_type))
{
  mutations = subset(mutations, annotation_impact=="MODIFIER")
}
dim(mutations)
# [1] 420   35 > TU_promoter


##### Creating a table with mutations counts per loci

## Loading CDS annotation file
if(grepl("^CDS_only",feature_type))
{
  cds_table_file = paste(dir,"annotation/Mycobacterium_tuberculosis_H37Rv_gff_v4.cds_200bp_downstream.txt.gz", sep = "");
  cds_table = read.delim(gzfile(cds_table_file) , sep = "\t", header = T)
  dim(cds_table)
  # [1] 4411531       6
}


## Getting list of unique loci Ids depending on unit of association selected
if(grepl("^CDS_only",feature_type))
{
  loci = unique(as.vector(gene_table_slim$Locus))
  length(loci)
  # [1] 4031 
  # [1] 4173 --> updated 
}

if(grepl("^TU",feature_type) & !grepl("^TU_cluster",feature_type) & !grepl("^TU_promoter",feature_type))
{
  loci = unique(unlist(strsplit(as.vector(TU_table$TU_id),";"))); # number of unique TU
  tmp = which(loci=="-"); if(length(tmp)>0){ loci = loci[-tmp]; }
  length(loci)
  # [1] 2573 > TU_CDS_only_functional
}

if(grepl("^TU_promoter",feature_type))
{
  loci = unique(unlist(strsplit(as.vector(TU_promoter_table$promoter_id),";"))); # number of unique TU
  tmp = which(loci=="-"); if(length(tmp)>0){ loci = loci[-tmp]; }
  length(loci)
  # [1] 2572 > TU_promoter
}


## Counting mutations per loci
loci_count_table = mat.or.vec(length(loci), 8); # table to keep locus name, mutations count, locus length, samples count

if(grepl("^CDS_only",feature_type))
{
  colnames(loci_count_table) = c('locus_id','feature_type','mutations_count','patient_count','patient_id','samples_count','samples_id','isolate_id')
}

if(grepl("^TU_CDS_only_functional",feature_type))
{
  colnames(loci_count_table) = c('TU_id','feature_type','mutations_count','patient_count','patient_id','samples_count','samples_id','isolate_id')
}

if(grepl("^TU_promoter",feature_type))
{
  colnames(loci_count_table) = c('promoter_id','feature_type','mutations_count','patient_count','patient_id','samples_count','samples_id','isolate_id')
}

for(c in 1:length(loci))
{
  print(c)
  locus_id = as.character(loci[c]); locus_id_grep = paste(locus_id,sep="");
  if(feature_type=="CDS"){ mutations_idx = which(grepl(locus_id_grep,mutations$locus_tag)); }
  if(grepl("^CDS_only",feature_type)){ mutations_idx = which(mutations$locus_tag == locus_id); }
  if(grepl("^TU",feature_type) & !grepl("^TU_cluster",feature_type) & !grepl("^TU_promoter",feature_type)){ mutations_idx = which(grepl(locus_id_grep,mutations$TU_id)); }
  if(grepl("^TU_promoter",feature_type)){ mutations_idx = which(grepl(locus_id_grep,mutations$TU_promoter_id)); }
  mutations_count = length(mutations_idx)
  samples_count = length(unique(mutations[mutations_idx,"sample"])); samples_id = paste(unique(mutations[mutations_idx,"sample"]),collapse = ";");
  patient_count = length(unique(mutations[mutations_idx,"patient_id"])); patient_id = paste(unique(mutations[mutations_idx,"patient_id"]),collapse = ";");
  isolate_id = paste(unique(mutations[mutations_idx,"isolate_id"]),collapse = ";");
  loci_count_table[c,] = c(locus_id,feature_type, mutations_count, patient_count, patient_id, samples_count, samples_id, isolate_id);
}
# at the end of this loops, mutations counts are CHAR type 

dim(loci_count_table)
# [1] 4173    8 > updated CDS only 
# [1] 2573    8 > updated TU_CDS_only_functional
# [1] 2572    8 > updated TU_promoter
# [1] 289     8 > submodule_CDS_only_functional


write.table(loci_count_table, file = mutations_count_table, sep = '\t', col.names = T, row.names = F, quote = F);


###################################################################################################################
########                                          ADD LOCI ANNOTATION                                          ####
###################################################################################################################

if(grepl("^CDS_only",feature_type)){ loci_ann_table = read.delim(gff_txt_file, sep = "\t", header = T); }
if(grepl("^TU",feature_type) & !grepl("^TU_cluster",feature_type)){ loci_ann_table = TU_table; }
if(grepl("^TU_promoter",feature_type) & !grepl("^TU_cluster",feature_type)){ loci_ann_table = TU_promoter_table; }

dim(loci_ann_table)
# [1] 4187    8 > CDS only  
# [1] 2573    7 > TU_CDS_only_functional
# [1] 2572    4 > TU_promoter

loci_count_table = read.delim(mutations_count_table, sep = '\t', header = T)
dim(loci_count_table)
# [1] 4173    8 > CDS_only
# [1] 2573    8 > TU_CDS_only_functional
# [1] 2572    8 > TU_promoter

if(grepl("^CDS_only",feature_type)){ loci_count_table = merge(loci_count_table, loci_ann_table, by.x = "locus_id", by.y = "Locus"); }
if(grepl("^TU",feature_type) & !grepl("^TU_cluster",feature_type) & !grepl("^TU_promoter",feature_type)){ loci_count_table = merge(loci_count_table, loci_ann_table, by.x = "TU_id", by.y = "TU_id"); }
if(grepl("^TU_promoter",feature_type)){ loci_count_table = merge(loci_count_table, loci_ann_table, by.x = "promoter_id", by.y = "promoter_id"); }

dim(loci_count_table)
# [1] 4187   41 > CDS_only
# [1] 2573    14 > TU_CDS_only_functional
# [1] 2572    11 > TU_promoter

if(grepl("^CDS_only",feature_type)){ loci_count_table$gene_length = (loci_count_table$Stop - loci_count_table$Start)+1; }
if(grepl("^TU",feature_type) & !grepl("^TU_cluster",feature_type) & !grepl("^TU_promoter",feature_type) & grepl("CDS_only",feature_type)){ loci_count_table$gene_length = (loci_count_table$end - loci_count_table$start)+1; }
if(grepl("^TU_promoter",feature_type)){loci_count_table$gene_length = 100; }

all_transcription_unit$TU_length = (all_transcription_unit$end - all_transcription_unit$start)+1

###################################################################################################################
########                                          ASSOCIATION ANALYSIS                                         ####
###################################################################################################################

total_number_mutations = nrow(mutations)
chromosome_num_loci = nrow(loci_count_table)
chromosome_length = 4411532

#### Removing loci with length shorter than min_UA_length (expect for TU_promoter, which are all 200bp)
if(grepl("^CDS_only",feature_type))
{
  loci_count_table = loci_count_table[-which(loci_count_table[,"gene_length"]<min_UA_length),]
}

if(grepl("^TU",feature_type) & !grepl("^TU_cluster",feature_type) & !grepl("^TU_promoter",feature_type) & !grepl("CDS_only",feature_type))
{
  loci_count_table = loci_count_table[-which(loci_count_table[,"TU_length_mader2016"]<min_UA_length),]
}

if(grepl("^TU",feature_type) & !grepl("^TU_cluster",feature_type) & !grepl("^TU_promoter",feature_type) & grepl("CDS_only",feature_type))
{
  loci_count_table = loci_count_table[-which(loci_count_table[,"gene_length"]<min_UA_length),]
}

dim(loci_count_table)
# [1] 3728   42 > CDS_only 
# [1] 2406   15 > TU_CDS_only_functional
# [1] 2572   12 > TU promoter 

if(grepl("^CDS_only",feature_type))
{
  p_values <- apply(loci_count_table,1,function(x) {poisson_test(x[["mutations_count"]], x[["gene_length"]], total_number_mutations, chromosome_length)})
}
if(grepl("^TU",feature_type) & !grepl("^TU_cluster",feature_type) & !grepl("^TU_promoter",feature_type) & !grepl("CDS_only",feature_type))
{
  p_values <- apply(loci_count_table,1,function(x) {poisson_test(x[["mutations_count"]], x[["TU_length_mader2016"]], total_number_mutations, chromosome_length)})
}

if(grepl("^TU",feature_type) & !grepl("^TU_cluster",feature_type) & !grepl("^TU_promoter",feature_type) & grepl("CDS_only",feature_type))
{
  p_values <- apply(loci_count_table,1,function(x) {poisson_test(x[["mutations_count"]], x[["gene_length"]], total_number_mutations, chromosome_length)})
}

if(grepl("^TU_promoter",feature_type))
{
  p_values <- apply(loci_count_table,1,function(x) {poisson_test(x[["mutations_count"]], x[["gene_length"]], total_number_mutations, chromosome_length)})
}

loci_count_table = cbind(loci_count_table, p_values)
p_values_adj <- p.adjust(p_values, method="BH", n=chromosome_num_loci)
loci_count_table = cbind(loci_count_table, p_values_adj)
p_value_cutoff = rep(as.character(0.05/chromosome_num_loci), nrow(loci_count_table))
loci_count_table = cbind(loci_count_table, p_value_cutoff)

# for CDS_only
loci_count_table_slim = subset(loci_count_table, select = c(locus_id, feature_type, mutations_count, patient_count, patient_id, samples_count, samples_id, isolate_id, 
                                                            Feature, Start, Stop, Strand, Name, Function, Product, Comments, gene_length, p_values, p_values_adj, p_value_cutoff))


# add locus IDs and names to TU CDS only loci count table 
loci_count_table$locus_id = "-"
loci_count_table$Name = "-"
for(r in 1:nrow(loci_count_table))
{
  samples = unlist(strsplit(loci_count_table$TU_id[r] , split = ";"))
  tmp = grep(samples, mutations$TU_id)
  if(length(tmp)>0)
  {
    loci_count_table$locus_id[r] = as.character(paste(unique(mutations$locus_tag[tmp]),collapse = ";"))
    loci_count_table$Name[r] = as.character(paste(unique(mutations$locus_name[tmp]),collapse = ";"))
  }
}

# add locus IDs and names to TU promoter loci count table 
# used gene name instead of locus id
loci_count_table$locus_id = "-"
loci_count_table$Name = "-"
for(r in 1:nrow(loci_count_table))
{
  samples = unlist(strsplit(loci_count_table$promoter_id[r] , split = ";"))
  tmp = grep(samples, mutations$TU_promoter_id)
  if(length(tmp)>0)
  {
    loci_count_table$locus_id[r] = as.character(paste(unique(mutations$locus_tag[tmp]),collapse = ";"))
    loci_count_table$Name[r] = as.character(paste(unique(mutations$gene_name[tmp]),collapse = ";"))
  }
}

# for CDS only
write.table(loci_count_table_slim, file=mutations_count_table_assoc, sep = '\t', col.names = T, row.names = F, quote = F)

# for everything else 
write.table(loci_count_table, file=mutations_count_table_assoc, sep = '\t', col.names = T, row.names = F, quote = F)




