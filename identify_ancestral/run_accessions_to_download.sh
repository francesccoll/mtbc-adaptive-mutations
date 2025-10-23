# NOTE: the command below is used to select the isolates to download, which include multiple isolates per host + selected outgroup
cat *.output_table.mp.csv | awk -F'\t' '{ print $2";"$6}' | grep -v "^host_isolates" | tr ';' '\n' | sort | uniq > run_accessions_to_download.txt
