# General
My repository for random code. Someday it will be organized!  

## Extracting genes from GWAS summary statistics
First, download gene positions using UCSC table browser. Make a space-delimited table called gene_positions.tab that has GENE_NAME, CHR, BP_START, BP_END.  

This loop goes line-by-line in the gene_positions.tab file, turns each line (gene) into an array, and searches for these positions in your GWAS summary stats file. In this case, column 1 is CHR and column 2 is BP in gwas_summary_stats.tab.  

```
cat gene_positions.tab | while read line
do
IFS=' ' read -a arr <<< "$line"
awk -v VAR="${arr[*]}" '{split(VAR,a," "); \
if ($1==a[2] && $2>=a[3] && $2<=a[4]) {$(NF+1)=a[1]; print}}' gwas_summary_stats.tab > ${arr[0]}_stats.txt
done
