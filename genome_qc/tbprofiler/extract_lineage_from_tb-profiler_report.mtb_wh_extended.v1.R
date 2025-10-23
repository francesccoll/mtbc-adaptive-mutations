

table0 = as.matrix(read.delim("../datasets/mtb_wh_extended.v1.run_accessions.txt", sep = "\t", header = F))
dim(table0)
# [1] 19582     1
# [1] 19681     1
colnames(table0) = "sample"

table1 = read.delim("mtb_wh_extended.v1.tb-profiler.strain.tmp.csv", sep = "\t", header = T)
dim(table1)
# [1] 19250     2
# [1] 19863     2

table2 = merge(table0, table1, by.x = "sample", by.y = "sample", all.x = T)
dim(table2)
# [1] 19582     2
# [1] 19681     2

# to label mixed strains
table2$tb_profiler_strain = gsub(" ", "", table2$tb_profiler_strain)
mixed_strain = rep(NA, nrow(table2))
tmp = which(grepl("lineage", table2$tb_profiler_strain)==TRUE)
mixed_strain[tmp] = "no"
tmp = which(grepl(";", table2$tb_profiler_strain)==TRUE)
mixed_strain[tmp] = "yes"
table2 = cbind(table2, mixed_strain)

# samples with no lineage assigned > the few investigated had very low coverage
table2$sample[which(is.na(table2$tb_profiler_strain))]

lineages = sort(table(table2$tb_profiler_strain[which(table2$mixed_strain=="no")]))
out_table = mat.or.vec(length(lineages), 2)
out_table[,1] = names(lineages)
out_table[,2] = lineages
colnames(out_table) = c("lineage", "num_samples")

write.table(table2, file = "mtb_wh_extended.v1.tb-profiler.strain.csv", sep = "\t", col.names = T, row.names = F, quote = F)
write.table(out_table, file = "mtb_wh_extended.v1.tb-profiler.num_strain.csv", sep = "\t", col.names = T, row.names = F, quote = F)

## NOTE: file mtb_wh_extended.v1.tb-profiler.num_strain.csv was manually extended to add the field 'alignment', that is, how samples will be grouped 
# into multiple alignments for phylogenetic analyses


