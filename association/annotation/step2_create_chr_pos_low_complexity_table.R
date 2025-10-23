

dustmasker_output_file="/Users/francesccoll/fellowship/2.projects/67.mtb_within_host/reference/Mycobacterium_tuberculosis_ancestral_reference.dustmasker.txt"
chr_table_file="/Users/francesccoll/fellowship/2.projects/67.mtb_within_host/association_mrca/annotation/Mycobacterium_tuberculosis_H37Rv_gff_v4.cds_200bp_downstream.txt.gz"
output_file="/Users/francesccoll/fellowship/2.projects/67.mtb_within_host/association_mrca/annotation/Mycobacterium_tuberculosis_ancestral_reference.dustmasker.chr_pos.txt"

lmD = read.delim(dustmasker_output_file, sep = " ", header = T)
dim(lmD)
# [1] 765   3

window_lc = 10; # window around low complexity region to label as low complexity region as well

low_complexity_positions = vector()
for(r in 1:nrow(lmD))
{
  from = as.numeric(lmD[r,1]); to = as.numeric(lmD[r,3]); from = from - window_lc; to = to + window_lc;
  low_complexity_positions = c(low_complexity_positions, seq(from, to, 1))
}
low_complexity_positions = unique(sort(low_complexity_positions))
length(low_complexity_positions)
# [1] 94755


chr_table = read.delim(gzfile(chr_table_file) , sep = "\t", header = T)
chr_table = chr_table[,c(1:2)]
dim(chr_table)
# [1] 4411531       2


# Adding low complexity position to chromosome table
low_complexity = rep("no", nrow(chr_table))
low_complexity[low_complexity_positions] = "yes"
chr_table = cbind(chr_table, low_complexity)

write.table(chr_table, file=output_file, sep = '\t', col.names = T, row.names = F, quote = F)





