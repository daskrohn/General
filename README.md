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
````

Compile all files into one, take only unique positions, and add header:
````
cat *_stats.txt > all_sumStats.txt

awk '!_[$3]++' all_sumStats.txt > unique_all_sumStats.txt
mv unique_all_sumStats.txt all_sumStats.txt

sed -i '1i CHR BP SNP A1 A2 beta se Zscore MAF P N' all_sumStats.txt
````
Use pruning to avoid over-correction:
````
cut -f3 -d' ' all_sumStats.txt > extract_variants.txt
plink --bfile ~/runs/krohn/krohn/Projects/GWAS/PRS/imptd_RBD_small-ctl --extract extract_variants.txt --indep-pairwise 250 5 0.5 --out pruning
wc -l pruning.prune.in # this number is number of tests
# in this case, tests = 358899
````
Sort by p-value and (optional) select only those with p < 0.05:
````
sort -gk10 all_sumStats.txt > p.all_sumStats.txt
awk '{if ($10<=0.05) print $0}' p.all_sumStats.txt > filtered_sumStats.txt
````
(Optional) clump summary stats for top variants:
````
plink --bfile ~/runs/krohn/krohn/Projects/GWAS/PRS/imptd_RBD_small-ctl --clump p.all_sumStats.txt --clump-kb 250 --clump-p1 0.05 --clump-p2 1 --clump-r2 0.5 --out p.all_sumStats.txt
````
Adjust p-values (Benjamini-Hochberg FDR):
```R
require(data.table)
data = fread("filtered_sumStats.txt", header = T)

data$p_FDR = p.adjust(data$P, "BH", 358899)

write.table(data, file = "FDR_filtered_sumStats.txt", sep = "\t", row.names = F, quote = F)
```
