# PCA Generation for 450K (and previous WES releases with GL) and loading of variants

Overview of the PCA generation pipeline for the 450K.

[Backman, J.D., Li, A.H., Marcketta, A. et al. Exome sequencing and analysis of 454,787 UK Biobank participants. Nature 599, 628–634 (2021).](https://doi.org/10.1038/s41586-021-04103-z)

## 1 Setup directories 

Creates required output directories

## 2 Generate Subject Lists

Steps:

- 2.1: Genetate Subject Lists
    - Loads source subjects from *UKB_Freeze_450.GL.splitmulti.fam* 
    - Adds subject genders: reported and sequenced
    - Removes discordant genders (turns out these are already removed from the source data)
    - Adds related flag using *UKB_Freeze_450.NF.pVCF.genome.3rd_degree_unrelated*
    - Bridges sample_id to eid using *WES_26041_450k_bridge.csv*
    - Saves bridged samples: *reference/subjects.tsv*, *reference/subject_ids.tsv*
    - Saves samples not-bridged (removed): *reference/removed_subjects.tsv*, *reference/removed_subject_ids.tsv*
    - Removes subjects > 1% missingness
    - Builds populations (self-reported ethnicities) using [field 21000](https://biobank.ctsu.ox.ac.uk/crystal/field.cgi?id=21000)
- 2.2: Set removed (non-bridged) IDs to -1
- 2.3: Move original BED/BIM to processed directory

## 3 Prepare Plink

Steps:

- 3.1: Subsetting source PLINK data
- 3.2: Linking source PLINK data
    - Fix CPRAs (IDs): 2nd column has `1:69159:I:1` or `1:69159:D:1` instead of `1:69159:G:GT` or `1:69159:GT:G`

## 4 PCA

Generates PCs using [Eigensoft's](https://github.com/DReichLab/EIG) [convertf](https://github.com/DReichLab/EIG/tree/master/CONVERTF) and [smartpca](https://github.com/DReichLab/EIG/tree/master/POPGEN).

Steps:

- 4.1: PLINK tasks for the `$E` population (initial pass)
    - Apply basic filters: maf &gt; 0.001, hwe p &gt; 1x10-12, genotype missingness &lt; 2%
    - Remove regions of [high LD](https://genome.sph.umich.edu/wiki/Regions_of_high_linkage_disequilibrium_(LD))
    - LD SNP pruning “indep-pairwise 1000 100 0.15” (black pop: 1000 100 0.1)
    - Remove INDELS (“snps-only")
- 4.2: Modify FAM file for the `$E` population (initial pass)
    - Adds *unrelated|related* to generated FAM file from 4.1
- 4.3: Run Eigensoft tasks for the `$E` population (initial pass)
    - Generates **SIX** eigenvectors with no outlier removal using *unrelated* population and fastmode
- 4.4: Find outliers for the `$E` population
    - Identifying subjects that are outliers (+/- 3SD) for any of the first 6 inital PCs
    - Generates the final eids for the population
- 4.5: PLINK tasks for the `$E` population (final pass)
    - Filter subjects (using final eids); results in final
    - Creates a backup of FAM file before modifying
- 4.6: Modify FAM file for the `$E` population (final pass)
    - Adds *unrelated|related* to generated FAM file from 4.5
- 4.7: Running Eigensoft tasks for `$E` population (final pass)
    - Generates **30** eigenvectors with no outlier removal using *unrelated* population and fastmode

Where `$E` are the populations *asian black chinese white*

Detail steps given in [Detailed PCA Overview](##4-pca-overview)

## 5 Covariate Sets Generation

Uses the generated PCs, filtered EIDs and additional UK Biobank fields to generate various covariate sets bases on the populations. 

For the 450K, an additional EID (4931830) was identified as an outlier and removed.

Fields created:

- PC_*: The PCs generated in step 4.7 
- Is_Relative: based on the string from the `related` column from step 4.7
    - Description says computed by P4, didn't we just used the supplied file (TODO DOUBLE CHECK) 
- Age: Age at recruitment [field 21022](https://biobank.ndph.ox.ac.uk/ukb/field.cgi?id=21022)
- Sex: Self-reported sex [Field 31](https://biobank.ndph.ox.ac.uk/ukb/field.cgi?id=31); set to 0 if female, 1 if male.
    - Genetic sex is [Field 22001](https://biobank.ctsu.ox.ac.uk/crystal/field.cgi?id=22001). Individuals with gender (sex) discordance were removed from the primary data (not in UKB_Freeze_450.GL.splitmulti.fam)
- Scottish_Recruited: recruited from Edinburgh or Glasgow, based on Recruitment Assessment Centre [field 54](https://biobank.ndph.ox.ac.uk/ukb/field.cgi?id=54); set to 1 if recruited, 0 otherwise.
- Welsh_Recruited: recruited from Cardiff, Swansea or Wrexham [field 54](https://biobank.ndph.ox.ac.uk/ukb/field.cgi?id=54); set to 1 if recruited, 0 otherwise.
- Has_GP: Uses the BioBank function `get_gp_clinical` which returns clinical data on eids; set to 1 if has GP, 0 otherwise

Sets created:

- GL_450K_`POP`_COVARIATES: Covariates for `POP` individuals, excluding relatives
- GL_450K_`POP`_COVARIATES_RELATED: Covariates for `POP` individuals, including relatives

## 6 Loading Variants

In *load/load_450k/load_450K.Rmd*:

    - load variants found in *data/pVCF/all_variants/processed/ref.bim* (This was manually moved in step 2.3) as chrom, pos, ref, alt. 
    - load genotypes found in *data/pVCF/all_variants/processed/bgen_per_chromosome/chr%i.bgen* using scidb's `load_bgen`

# Details

## 2 Generate Subject Lists

### 2.1 gen_subject_lists

R script


### 2.2 Remove IDs 

Run notebook *fix_fam_file.Rmd* which sets removed IDs to negative numbers. 

This allows skipping another plink run on the very large `.bed` file, creating another version without those genotypes. In Step 3, PLINK should then just skip these IDs.

Inputs:
    - */data/pVCF/all_variants/UKB_Freeze_450.GL.splitmulti.fam* (source data)
    - *reference/removed_subject_ids.tsv* (generated in 2.1)
Outputs:
    - */data/pVCF/all_variants/processed/UKB_Freeze_450.GL.splitmulti.fam*

### 2.3 Move BIM/BAM files


Inputs:
    - */data/pVCF/all_variants/UKB_Freeze_450.GL.splitmulti.[bim|bam]* (moves source data)
Outputs:
    - */data/pVCF/all_variants/processed/UKB_Freeze_450.GL.splitmulti.[bim|bam]*

Have
    ref.bed  ref.bim  ref.fam  UKB_Freeze_450.GL.splitmulti.bim.BAK  UKB_Freeze_450.GL.splitmulti.fam.BAK

So what gives?

## 3 Prepare PLINK

### 3.1 subset_plink_data


### 3.2 link_plink_data


## 4 PCA Overview

### 4.1 run_plink_tasks_initial

Sets working directory: `cd output/$E/initial`

Task 1: Apply basic filters

```bash
        $PLINK2 \
            --bfile $DATA_PATH/$BFILE \
            --geno 0.02 \
            --hwe 1e-12 \
            --keep eids.tsv \
            --maf 0.001 \
            --make-bed \
            --out plink_task_1 \
            --output-chr 26
```

Input: *data/pVCF/all_variants/processed/UKB_Freeze_450.GL.splitmulti.[bed|bim|fam]*
Output: *plink_task_1.[bed|bim|fam]*

Task 2: Remove regions of long LD

```bash
        $PLINK2 \
            --bfile plink_task_1 \
            --exclude range ../../../reference/excluded_ranges.bed \
            --make-bed \
            --out plink_task_2 \
            --output-chr 26
```

The `--output-chr 26` option converts "chr1" to "1" in the IDs

Inputs:
	- *plink_task_1.[bed|bim|fam]*
	- *reference/excluded_ranges.bed*  (original is *pca/200k/reference/excluded_ranges.bed* via Meg)
Output: *plink_task_2.[bed|bim|fam]*

Task 3: LD pruning

**NOTE**, uses different R^2 threshold for `black` population

```bash
        [ $E = 'black' ] && T=0.10 || T=0.15
        $PLINK2 \
            --bfile plink_task_2 \
            --indep-pairwise 1000 100 $T \
            --out plink_task_3
```

Input: *plink_task_2.[bed|bim|fam]*
Outputs: *plink_task_3.prune.[in|out]*

Task 4: Export filtered and LD-pruned results

```bash
        $PLINK2 \
            --bfile plink_task_2 \
            --extract plink_task_3.prune.in \
            --make-bed \
            --out plink_task_4 \
            --output-chr 26
```

Input: *plink_task_2.[bed|bim|fam]*
Output:  *plink_task_4.[bed|bim|fam]*

Task 5: Remove indels so SNP names are not too long

Proabably could combine with Task 4 (check for other usage of plink_task_4)

```bash
        $PLINK2 \
            --bfile plink_task_4 \
            --make-bed \
            --out plink_task_5 \
            --output-chr 26 \
            --snps-only
```

Input: *plink_task_4.[bed|bim|fam]*
Output:  *plink_task_5.[bed|bim|fam]*

### 4.2 modify_fam_initial

Adds **unrelated/related** field as 6th column for use in `convertf` and `smartpca`

```R
E <- args[1]
setwd(paste0('output/', E, '/initial'))

f <- data.table::fread(FAM_PATH, header = F)
s <- data.table::fread(SUBJECTS_PATH, header = T)
related <- s[s$related]$eid
f$V7 <- 'unrelated'
f$V7[f$V1 %in% related] <- 'related' 
data.table::fwrite(f[, c('V1', 'V2', 'V3', 'V4', 'V5', 'V7')], FAM_PATH, col.names = F, sep = '\t')
```

Inputs:
    - *reference/subjects.tsv* (generated in `2_1_gen_subject_lists.R`)
    - *output/$E/intitial/plink_task_5.fam* (generated in `4_1_run_plink_tasks_initial.sh` task 5 )
Outputs:
    - *output/$E/intitial/plink_task_5.fam.BAK* (**copy of original**)
    - *output/$E/intitial/plink_task_5.fam* (**replaces original**)


### 4.3 run_eig_tasks_initial

Sets working directory: `cd output/$E/initial`

Task 1: Convert BED genotypes to EIGENSTRAT format

`$CONVERTF -p ../../../reference/convertf_initial.par`
	
In *convertf_initial.par*

```
genotypename: plink_task_5.bed
snpname: plink_task_5.bim
indivname: plink_task_5.fam
outputformat: EIGENSTRAT
genotypeoutname: eig_task_1.eigenstratgeno
snpoutname: eig_task_1.snp
indivoutname: eig_task_1.ind
```

Inputs:
	- *output/$E/initial/plink_task_5.[bed|bim|fam]* (from 4.1 Task 4 and 4.2)
	
Outputs:
	- *output/$E/initial/eig_task_1.eigenstratgeno*
	- *output/$E/initial/eig_task_1.snp*
	- *output/$E/initial/eig_task_1.ind*

Task 2: Generate SIX eigenvectors (PCs)

`$SMARTPCA -p ../../../reference/smartpca_initial.par`

In *smartpca_initial.par*:

```
genotypename: eig_task_1.eigenstratgeno
snpname: eig_task_1.snp
indivname: eig_task_1.ind
evecoutname: eig_task_2.evec
evaloutname: eig_task_2.eval
altnormstyle: NO
numoutevec: 6
numoutlieriter: 0
numoutlierevec: 0
outliersigmathresh: 0
snpweightoutname: eig_task_2.snpweight.txt
poplistname: ../../../reference/poplistname.txt
fastmode: YES
```

Inputs:
	- *output/$E/initial/eig_task_1.eigenstratgeno*
	- *output/$E/initial/eig_task_1.snp*
	- *output/$E/initial/eig_task_1.ind*
	- *../../../reference/poplistname.txt* (uses only **unrelated**)
	
Outputs:
    - *output/$E/initial/eig_task_2.evec* (used in 4.4)
    - *output/$E/initial/eig_task_2.eval*
    - *output/$E/initial/eig_task_2.snpweight.txt*
	
Defaults are:
- `numoutevec`: 10
- `numoutlieriter`: 5. To turn off outlier removal, set this parameter to 0 (which it is)
- `snpweightoutname`: output file containing SNP weightings of each
  principal component.  Note that this output file does not contain entries
  for monomorphic SNPs from the input .snp file.
- `poplistname`:  If wishing to infer eigenvectors using only individuals from a 
  subset of populations, and then project individuals from all populations 
  onto those eigenvectors, this input file contains a list of population names,
  one population name per line, which will be used to infer eigenvectors.  
  It is assumed that the population of each individual is specified in the 
  indiv file.  Default is to use individuals from all populations.

### 4.4 find_outliers

Identifying subjects that are outliers (+/- 3SD) for any of the first 6 inital PCs

```R
setwd(paste0('output/', E, '/initial'))

print('Reading in evec file...')
evec <- data.table::fread('eig_task_2.evec', header = F)
evec$V1 <- as.numeric(sapply(strsplit(evec$V1, ':'), function(x) { x[1] }))
names(evec) <- c('eid', 'pc1', 'pc2', 'pc3', 'pc4', 'pc5', 'pc6', 'poplistname')

print('Calculating PC1 cutoffs...')
mean_pc1 <- mean(evec$pc1)
cut_pc1 <- 3 * sd(evec$pc1)
lower_bound_pc1 <- mean_pc1 - cut_pc1
upper_bound_pc1 <- mean_pc1 + cut_pc1
evec$out1 = ifelse(evec$pc1 > lower_bound_pc1 & evec$pc1 < upper_bound_pc1, 0, 1)

rint('Identifying subjects that are outliers for any of the first 6 inital PCs...')
evec$out = ifelse(evec$out1 == 0 & evec$out2 == 0 & evec$out3 == 0 & evec$out4 == 0 & evec$out5 == 0 & evec$out6 == 0, 0, 1)

print('Writing analysis results...')
data.table::fwrite(evec, file = 'pc_cuts.tsv', sep = '\t')

print('Writing EIDs for final pass...')
eids_in <- subset(evec, evec$out == 0)[, c('eid', 'eid')]
data.table::fwrite(eids_in, file = '../final/eids.tsv', col.names = F, sep = '\t')
```

Inputs: 
    - *output/$E/initial/eig_task_2.evec* (from 4.3)
Outputs:
    - *output/$E/initial/pc_cuts.tsv* (eig_task_2.evec with 'out' columns added)
    - *output/$E/final/eids.tsv*  (used in 4.5)

### 4.5 run_plink_tasks_final

Sets working directory: `cd output/$E/final`

Filter subjects

```bash
        $PLINK2 \
            --bfile ../initial/plink_task_5 \
            --keep eids.tsv \
            --make-bed \
            --out plink_task_1 \
            --output-chr 26
```

Inputs: 
    - *output/$E/final/eids.tsv* (generated in 4.4)
    - *output/$E/initial/plink_task_5.[bed|bim|fam]* (generated in `4_1_run_plink_tasks_initial.sh` task 5 and `4_2_modify_fam_initial.R`)
Outputs:
    - *output/$E/final/plink_task_1.[bed|bim|fam]*

### 4.6 modify_fam_final

```R
E <- args[1]

setwd(paste0('output/', E, '/final'))
```

Inputs:
    - *output/$E/final/plink_task_1.[bed|bim|fam]* (4.5's version of 4.1 Task 1's output) 
Outputs:
    - *output/$E/final/plink_task_1.fam.BAK* (**copy of original**)
        - attempts to restore if run again
    - *output/$E/final/plink_task_1.fam* (**replaces original**) 

### 4.7 run_eig_tasks_final

Sets working directory: `cd output/$E/final`

Task 1: Convert BED genotypes to EIGENSTRAT format

`$CONVERTF -p ../../../reference/convertf_final.par`

In *convertf_final.par*

```
genotypename: plink_task_1.bed
snpname: plink_task_1.bim
indivname: plink_task_1.fam
outputformat: EIGENSTRAT
genotypeoutname: eig_task_1.eigenstratgeno
snpoutname: eig_task_1.snp
indivoutname: eig_task_1.ind
```

Inputs:
    - *output/$E/final/plink_task_1.[bed|bim|fam]* (from 4.5 and 4.6)
Outputs:
	- *output/$E/final/eig_task_1.eigenstratgeno*
	- *output/$E/final/eig_task_1.snp*
	- *output/$E/final/eig_task_1.ind*


Task 2: Generate <> eigenvectors (PCs)

`$SMARTPCA -p ../../../reference/smartpca_final.par`

In *smartpca_final.par*:

```
genotypename: eig_task_1.eigenstratgeno
snpname: eig_task_1.snp
indivname: eig_task_1.ind
evecoutname: eig_task_2.evec
evaloutname: eig_task_2.eval
altnormstyle: NO
numoutevec: 30
numoutlieriter: 0
numoutlierevec: 0
outliersigmathresh: 0
snpweightoutname: eig_task_2.snpweight.txt
poplistname: ../../../reference/poplistname.txt
fastmode: YES
```

Inputs:
	- *output/$E/final/eig_task_1.eigenstratgeno*
	- *output/$E/final/eig_task_1.snp*
	- *output/$E/final/eig_task_1.ind*
	- *../../../reference/poplistname.txt* (uses only **unrelated**)

Outputs:
    - *output/$E/final/eig_task_2.evec* (used by `load_pca_to_scidb.Rmd`)
    - *output/$E/final/eig_task_2.eval*
    - *output/$E/final/eig_task_2.snpweight.txt*

