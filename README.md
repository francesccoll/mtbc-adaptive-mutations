# *In host* mutational adaptation of *Mycobacterium tuberculosis* complex strains

This GitHub project contains the data and code necessary to reproduce the findings of the study '*In host* mutational adaptation of *Mycobacterium tuberculosis* complex strains', and includes the following directories:
* genome_qc: scripts applied to individual MTBC genomes to derive genome QC metrics.
* genome_analyses: genome analyses applied to individual MTBC genomes (short-read de novo assembly and mapping).
* population_analyses: code to compute phylogenies and pairwise SNP distances (analyses conducted per lineage or grups of related sub-lineages).
* identify_ancestral: pipeline to reconstruct the nucleotide sequence of the most recent common ancestor (MRCA) of all isolates from the same strain.
* snippy_mrca: code to call genetic variants against the MRCA sequence of MTBC strains: detection of de novo acquired genetic variants during infection.
* association: R code used to filter within-host de novo acquired mutations, and to identify loci enriched by these mutations.
* dr_acquisition: files and code to determine drug-resistance acquisition rates.


# Citation
Zhang H, Medina-Jaudes N, Forcada-Nadal A, *et al.* *In host* mutational adaptation of *Mycobacterium tuberculosis* complex strains. BioRxiv.
