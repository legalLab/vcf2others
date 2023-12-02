# vcf2others
Convert VCF to Other Formats and More

``` r
library(vcf2others)
```

# Introduction

This package has two main functionalities. First, it is used for
filtering, subsetting, merging and otherwise wrangling VCF files.
Second, it is used for converting the VCF file to other population
genetic and phylogenetic formats for downstream analyses.

# Installation

This package needs to be installed from GitHub.
`devtools::install_github("legalLab/vcf2others")`

# How to use the functions of this package

Following are examples of the usage of the functions of this package.

## Load VCF and associated files of individuals and groups

First we need to load our files. In addition to the VCF, we also need to
have a white/black list of individuals (read from a file of one column
with an ‘id’ header), a strata file associating individuals with groups
(read from a file of two tab-separated columns with an ‘id’ and ‘pop’
headers), and a file of groups (read from a file of one column with a
‘pop’ header). These auxilary files are used for keeping/excluding
individuals and groups of individuals during filtering and during
conversion to other data formats. The strata file is also necessary
since many of the population genetic formats require information on
grouping of individuals.

``` r
library(vcfR)
#> 
#>    *****       ***   vcfR   ***       *****
#>    This is vcfR 1.14.0 
#>      browseVignettes('vcfR') # Documentation
#>      citation('vcfR') # Citation
#>    *****       *****      *****       *****
library(dplyr)
library(vcf2others)

# path to example files
fpath <- paste0(system.file("extdata", package="vcf2others"), "/")
# set path and project name
path <- "/home/tomas/git/legal_public/packages/data/"
project <- "trigonatus_all_sub"

# load vcf
vcf <- read.vcfR(paste0(fpath, project, ".vcf.gz"))
#> Scanning file to determine attributes.
#> File attributes:
#>   meta lines: 22
#>   header_line: 23
#>   variant count: 10000
#>   column count: 21
#> Meta line 22 read in.
#> All meta lines processed.
#> gt matrix initialized.
#> Character matrix gt created.
#>   Character matrix gt rows: 10000
#>   Character matrix gt cols: 21
#>   skip: 0
#>   nrows: 10000
#>   row_num: 0
#> Processed variant 1000Processed variant 2000Processed variant 3000Processed variant 4000Processed variant 5000Processed variant 6000Processed variant 7000Processed variant 8000Processed variant 9000Processed variant 10000Processed variant: 10000
#> All variants processed

# read individuals to include
indivs <- read.table(paste0(fpath, "indivs_b"), header = TRUE)$id %>%
  as.character()

# check if all indivs are in vcf_names
if (any(!(indivs %in% colnames(vcf@gt)[-1]))) stop(paste("Some individuals in list not in VCF"))

# read all sample names in vcf
vcf_names <- colnames(vcf@gt)[-1] %>%
  as_tibble() %>%
  rename(id = 1) %>%
  mutate(id, id = as.character(id))

# read sample to group assignment
strata <- read.table(paste0(fpath, "strata"), header = TRUE) %>%
  as_tibble() %>%
  mutate(id, id = as.character(id))

# assign samples in vcf to groups
strata <- vcf_names %>%
  left_join(strata)
#> Joining with `by = join_by(id)`

# check if all samples in vcf are assigned to groups
if (any(is.na(strata$pop) == TRUE)) {
  print("vcf has individuals not assigned to a group")
  print(strata[is.na(strata$pop) == TRUE,])
} else {
  indiv_group <- as.factor(strata$pop)
}

# read groups to be used; filter only those that are actually present
groups <- read.table(paste0(fpath, "groups"), header = TRUE)[,1] %>%
  .[. %in% unique(strata$pop)] %>%
  as.factor()
```

## Filtering, subsetting, merging and otherwise wrangling VCF files

All functions can be used directly on a VCF, but many of them are called
by the other functions of this package. For example
`vcf_filter_invariant()` will remove invariant SNPs from a VCF, however,
VCF by definition should not have invariant SNPs. However, if a VCF is
subsetted or certain individuals are removed, SNPs can become invariant
and need to be removed.

### Subsetting

The functions `vcf_extract_indivs()`, `vcf_extract_pops()` and
`vcf_sub_loci()` are used to subset a VCF by specific individuals, a
group (all individuals of a group/population), or a random subset of
loci, respectively. Subsetting individuals and groups can be done either
as a whitelist (keep the listed individuals/groups) or blacklist (remove
the listed individuals/groups).

``` r
# subset VCF by individuals (default whitelist = TRUE, remove invariants = TRUE)
indivs1 <- c("CTGA_H4683", "CTGA_H4689")
vcf1 <- vcf_extract_indivs(vcf, indivs1)

# subset VCF by groups (default whitelist = TRUE, remove invariants = TRUE)
groups1 <- c("bob", "joe") %>%
  as.factor()
vcf1 <- vcf_extract_pops(vcf, indiv_group, groups1)

# subset VCF by n random loci (default n_loci = 1000)
vcf1 <- vcf_sub_loci(vcf, n_loci = 5000)
```

### Filtering

Filtering consists of removing individuals and loci based on specific
properties. The function `vcf_filter_missing_indivs()` removes
individuals above a % missing SNP threshold. The other functions filter
on properties of individuals SNPs across individuals. The function
`vcf_filter_missingness()` removes SNPs above a % missing threshold,
`vcf_filter_quality()` removes SNPs below a quality threshold,
`vcf_filter_maf()` removes SNPs with minor allele frequency below a %
threshold, `vcf_filter_hets()` removes SNPs with heterozigosity above a
% threshold, `vcf_filter_rank()` removes SNPs with rank below a %
threshold (rank is specific to VCFs generated by DiscoSnpRad and
measures the likelihood of observed SNP variation being due to
paralogs), and `vcf_filter_coverage()` removes SNPs below a specific
coverage (makes them missing in the individual). Finally the
`vcf_filter_invariant()` removes invariant SNPs (SNPs may become
invariant after subsetting and filtering). The functions
`vcf_filter_oneSNP()` and `vcf_filter_multiSNP()` filter a VCF to remove
linked SNPs (retain only unliked SNPs) or retain only linked SNPs,
respectively. In general, all analyses assume that SNPs are unlinked,
however, some analyses such as FineRadStructure work with linked SNPs.

``` r
# filter VCF for analyses (unlinked SNPs)
vcf_oneSNP <- vcf_extract_indivs(vcf, indivs, whitelist = FALSE) %>%
  vcf_filter_missing_indivs(.9) %>%
  vcf_filter_rank(.5) %>%
  vcf_filter_maf(.03) %>%
  vcf_filter_coverage(6) %>%
  vcf_filter_oneSNP() %>%
  vcf_filter_missingness(.4) %>%
  vcf_filter_missing_indivs(.2)
#> removed samples are: CTGA_H3327 
#>  removed samples are: CTGA_H3600 
#> [1] "final % missing data in VCF is 57.17 %"
#> [1] "final % missing data in VCF is 25.42 %"
#> removed samples are: CTGA_H3646 
#>  removed samples are: CTGA_H4598 
#>  removed samples are: CTGA_H4613 
#>  removed samples are: CTGA_H4683 
#>  removed samples are: CTGA_H4689 
#> [1] "final % missing data in VCF is 2.92 %"

# filter VCF for analyses (linked SNPs)
vcf_multiSNP <- vcf_extract_indivs(vcf, indivs, whitelist = FALSE) %>%
  vcf_filter_missing_indivs(.9) %>%
  vcf_filter_rank(.5) %>%
  vcf_filter_maf(.03) %>%
  vcf_filter_coverage(6) %>%
  vcf_filter_multiSNP() %>%
  vcf_filter_missingness(.4) %>%
  vcf_filter_missing_indivs(.2)
#> removed samples are: CTGA_H3327 
#>  removed samples are: CTGA_H3600 
#> [1] "final % missing data in VCF is 57.17 %"
#> [1] "final % missing data in VCF is 25.2 %"
#> removed samples are: CTGA_H3646 
#>  removed samples are: CTGA_H4613 
#>  removed samples are: CTGA_H4683 
#>  removed samples are: CTGA_H4689 
#> [1] "final % missing data in VCF is 5.15 %"
```

### Extracting, merging and adding

The functions `vcf_merge()`, `vcf_add_indivs()`, `vcf_extract_indivs()`
and `vcf_extract_pops()` are used to merge two VCF files, add
individuals from a second VCF into a first VCF, and extract
individuals/groups from a VCF saving it into a separate VCF. The
`vcf_add_indivs()`, `vcf_extract_indivs()` and `vcf_extract_pops()`
functions were designed in mind to work with outgroups. Because of
phylogenetic effects, outgroups have relatively more missing data than
ingroups, and so often are removed during filtering. The solution is to
extract the outgroup taxa with `vcf_extract_indivs()` or
`vcf_extract_pops()`, filter the ingroup VCF, and add the outgroup taxa
with `vcf_add_indivs()`. It is important to NOT filter invariant loci
during the filtering of the ingroup, and only filter them
`vcf_filter_invariant()` after adding back the outgroup. The functions
`vcf_extract_indivs()`, `vcf_extract_pops()`,
`vcf_filter_missing_indivs()` and `vcf_filter_coverage()` by default
filter invariants, but this can be turned off using the “f_invar =
FALSE” flag.

``` r
# extract individuals from VCF (keep all loci)
indivs1 <- c("CTGA_H4644", "CTGA_H4661", "CTGA_H4683", "CTGA_H4689")
vcf_outgrp <- vcf_extract_indivs(vcf, indivs1, f_invar = FALSE)
vcf_ingrp <- vcf_extract_indivs(vcf, indivs1, whitelist = FALSE, f_invar = FALSE)

# extract groups of individuals from VCF (keep all loci)
groups1 <- as.factor("outgrp")
vcf_outgrp <- vcf_extract_pops(vcf, indiv_group, groups1, f_invar = FALSE)
vcf_ingrp <- vcf_extract_pops(vcf, indiv_group, groups1, whitelist = FALSE, f_invar = FALSE)

# merge vcf_outgrp into vcf_ingrp
vcf1 <- vcf_merge(vcf_ingrp, vcf_outgrp)
#> Joining with `by = join_by(FORMAT, id)`

# add individuals from vcf_outgrp into vcf_ingrp (whitelist = TRUE is default)
indivs1 <- c("CTGA_H4683", "CTGA_H4689")
vcf1 <- vcf_add_indivs(vcf_ingrp, vcf_outgrp, indivs1)
#> Joining with `by = join_by(FORMAT, id)`
```

### Filtering with outgroup

Outgroups, which tend to have many fewer individuals than ingroups, has
more missing data than ingroups due to phylogenetic effects. Filtering
the entire dataset often results in removal of the outgroup taxa when
`vcf_filter_invariant()` is used during filtering to remove individuals
with % missing data above some threshold. Therefore, a better strategy
may be to extract the outgroup, filter the ingroup, add the outgroup,
and remove any invariant SNPs.

``` r
# extract outgroup from VCF (keep all loci)
groups1 <- as.factor("outgrp")
vcf_outgrp <- vcf_extract_pops(vcf, indiv_group, groups1, f_invar = FALSE)

# filter ingroup VCF for analyses then add outgroups
vcf1 <- vcf_extract_pops(vcf, indiv_group, groups1, whitelist = FALSE, f_invar = FALSE) %>%
  vcf_filter_missing_indivs(.9, f_invar = TRUE) %>%
  vcf_filter_rank(.5) %>%
  vcf_filter_maf(.03) %>%
  vcf_filter_coverage(6) %>%
  vcf_filter_oneSNP() %>%
  vcf_filter_missingness(.4) %>%
  vcf_filter_missing_indivs(.2, f_invar = TRUE) %>%
  vcf_merge(vcf_outgrp) %>%
  vcf_filter_invariant()
#> removed samples are: CTGA_H3327 
#>  removed samples are: CTGA_H3600 
#> [1] "final % missing data in VCF is 41.24 %"
#> [1] "final % missing data in VCF is 11.45 %"
#> removed samples are: CTGA_H3646 
#> [1] "final % missing data in VCF is 6.32 %"
#> Joining with `by = join_by(FORMAT, id)`
```

## Converting a VCF file to other population genetic and phylogenetic formats

All the following functions will take a VCF, groups and indiv_group
(individuals to groups relationship), and convert the VCF to other
populations genetic and phylogenetic formats, exporting/writing a file
of this format. By default all functions will include missing data. If
that is not desired, only 100% complete data matrices can be exported
using the “inc_missing = FALSE” flag. The function `vcf2genlight()`
automatically returns a genlight object and optionally can also
export/write a file of this format to the working directory. Generally
the `vcf2genlight()` function is called within a script using functions
of the `adegenet` package rather than importing the genlight object.

``` r
project <- "trigonatus"

##########
# datasets for analyses with unlinked SNPs
vcf <- vcf_oneSNP

# read all sample names in vcf
vcf_names <- colnames(vcf@gt)[-1] %>%
  as_tibble() %>%
  rename(id = 1) %>%
  mutate(id, id = as.character(id))

# read sample to group assignment
strata <- read.table(paste0(fpath, "strata"), header = TRUE) %>%
  as_tibble() %>%
  mutate(id, id = as.character(id))

# assign samples in vcf to groups
strata <- vcf_names %>%
  left_join(strata)
#> Joining with `by = join_by(id)`

# check if all samples in vcf are assigned to groups
if (any(is.na(strata$pop) == TRUE)) {
  print("vcf has individuals not assigned to a group")
  print(strata[is.na(strata$pop) == TRUE,])
} else {
  indiv_group <- as.factor(strata$pop)
}

# read groups to be used; filter only those that are actually present
groups <- read.table(paste0(fpath, "groups"), header = TRUE)[,1] %>%
  .[. %in% unique(strata$pop)] %>%
  as.factor()

##########
# export data formats
# migrate-n https://peterbeerli.com/migrate-html5/
vcf2migrate(vcf, indiv_group, groups, out_file = paste0(path, project, '_migrateN.txt'))
# migrate-n - new more compact heterozygosity format
vcf2migrate(vcf, indiv_group, groups, out_file = paste0(path, project, '_migrateH.txt'), method = "H")
# arlequin 
vcf2arlequin(vcf, indiv_group, groups, out_file = paste0(path, project, '.arp'))
# structure 
vcf2structure(vcf, indiv_group, groups, out_file = paste0(path, project, '.str'))
# faststucture 
vcf2structure(vcf, indiv_group, groups, out_file = paste0(path, project, '.fstr'), method = "F")
# genepop 
vcf2genepop(vcf, indiv_group, groups, out_file = paste0(path, project, '.gen'))
# smartsnp 
vcf2smartsnp(vcf, indiv_group, groups, out_file = paste0(path, project, '.smartsnp'))
# eigenstrat 
vcf2eigenstrat(vcf, indiv_group, groups, out_file = paste0(path, project, '_eigenstrat'))
# bayescan 
vcf2bayescan(vcf, indiv_group, groups, out_file = paste0(path, project, '.bayescan'))
# treemix 
vcf2treemix(vcf, indiv_group, groups, out_file = paste0(path, project, '.treemix'))
# apparent 
vcf2apparent(vcf, indiv_group, groups, out_file = paste0(path, project, '.apparent'))
# related 
vcf2related(vcf, indiv_group, groups, out_file = paste0(path, project, '.related'))
# dataframe with population membership in the last column
vcf2df(vcf, indiv_group, groups, TRUE, out_file = paste0(path, project, '.df'))
# snapp 
vcf2snapp(vcf, indiv_group, groups, out_file = paste0(path, project, '_snapp.nex'))
# nexus - only SNPs, meant for SDVq analyses
vcf2nexus(vcf, indiv_group, groups, out_file = paste0(path, project, '_sdvq.nex'))
# fasta 
vcf2fasta(vcf, indiv_group, groups, TRUE, out_file = paste0(path, project, '.fas'))

##########
# datasets for analyses with linked SNPs
vcf <- vcf_multiSNP

# read all sample names in vcf
vcf_names <- colnames(vcf@gt)[-1] %>%
  as_tibble() %>%
  rename(id = 1) %>%
  mutate(id, id = as.character(id))

# read sample to group assignment
strata <- read.table(paste0(fpath, "strata"), header = TRUE) %>%
  as_tibble() %>%
  mutate(id, id = as.character(id))

# assign samples in vcf to groups
strata <- vcf_names %>%
  left_join(strata)
#> Joining with `by = join_by(id)`

# check if all samples in vcf are assigned to groups
if (any(is.na(strata$pop) == TRUE)) {
  print("vcf has individuals not assigned to a group")
  print(strata[is.na(strata$pop) == TRUE,])
} else {
  indiv_group <- as.factor(strata$pop)
}

# read groups to be used; filter only those that are actually present
groups <- read.table(paste0(fpath, "groups"), header = TRUE)[,1] %>%
  .[. %in% unique(strata$pop)] %>%
  as.factor()

##########
# export data formats
# fineRadStructure - expects VCF of linked SNPs
vcf2fineRadStructure(vcf, indiv_group, groups, out_file = paste0(path, project, '.finerad'))
```
