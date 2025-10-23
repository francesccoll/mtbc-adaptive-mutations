cat mtb_wh_extended.v1.all_lineages.output_table.mp.csv | grep -v "host_id" | awk -F'\t' '{ print $1}' | sort | uniq > mtb_wh_extended.v1.all_lineages.output_table.mp.strain_ids.txt
