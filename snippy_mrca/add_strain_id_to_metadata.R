

metadata_file = "/Users/francesccoll/fellowship/2.projects/67.mtb_within_host/metadata/mtb_wh_extended.v1.basic_metadata.csv"
strain_file = "/Users/francesccoll/fellowship/2.projects/67.mtb_within_host/identify_ancestral/isolate_ids.all_lineages.csv"
metadata_file2 = "/Users/francesccoll/fellowship/2.projects/67.mtb_within_host/snippy_mrca/mtb_wh_extended.v1.basic_metadata.mi.qc.csv"

metadata = read.delim(metadata_file, sep = "\t", header = T)
dim(metadata)
# [1] 22209     9

strain = read.delim(strain_file, sep = "\t", header = T)
dim(strain)
# [1] 11666     2

metadata2 = merge(strain, metadata, by.x = "isolate_id", by.y = "run_accession", all.x = T)
dim(metadata2)
# [1] 11995    10

# NOTE: some run accessions may be duplicated across different strain_ids/studies???

write.table(metadata2, file = metadata_file2, sep = "\t", col.names = T, row.names = F, quote = F)



