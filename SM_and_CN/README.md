## SM_and_CN

## Copy Number data
ASCAT data was made by ascat_script (https://github.com/Kazuki526/ascat_scripts)

- extract_all_PASS_maf.pl  
Extract mutations which are FILTER==PASS in all four mutation caller in TCGA.

- merge_ascat_and_allPASS_maf.pl  
connecting copy number infromation from ASCAT to somatic mutation

- merge_ascat_and_mutect_maf.pl  
Extract the read count information of all mutations extracted by mutect to calculate Beta, and add ASCAT information.

- calculate_beta_alpha.pl  
Calculate alpha and beta using the information above.

- extract_LOH_mutation.R  
Extracting LOH mutation 
