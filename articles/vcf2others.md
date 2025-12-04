# vcf2others

## Installation

This package needs to be installed from GitHub.
`devtools::install_github("legalLab/vcf2others")`

## Introduction

This package has two main functionalities. First, it is used for
filtering, subsetting, merging and otherwise wrangling VCF files.
Second, it is used for converting the VCF file to other population
genetic and phylogenetic formats for downstream analyses.

## How to use the functions of this package

Following are examples of the usage of the functions of this package.

### Load VCF and associated files of individuals and groups

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
#>    This is vcfR 1.15.0 
#>      browseVignettes('vcfR') # Documentation
#>      citation('vcfR') # Citation
#>    *****       *****      *****       *****
library(dplyr)
#> 
#> Attaching package: 'dplyr'
#> The following objects are masked from 'package:stats':
#> 
#>     filter, lag
#> The following objects are masked from 'package:base':
#> 
#>     intersect, setdiff, setequal, union
library(vcf2others)

# set path to example files and project name
data_path <- paste0(system.file("extdata", package="vcf2others"), "/")
project <- "crocs_"
postfix <- "discosnp_sub"

# load vcf
vcf <- read.vcfR(paste0(data_path, project, postfix, ".vcf.gz"))
#> Scanning file to determine attributes.
#> File attributes:
#>   meta lines: 22
#>   header_line: 23
#>   variant count: 10000
#>   column count: 26
#> Meta line 22 read in.
#> All meta lines processed.
#> gt matrix initialized.
#> Character matrix gt created.
#>   Character matrix gt rows: 10000
#>   Character matrix gt cols: 26
#>   skip: 0
#>   nrows: 10000
#>   row_num: 0
#> Processed variant 1000Processed variant 2000Processed variant 3000Processed variant 4000Processed variant 5000Processed variant 6000Processed variant 7000Processed variant 8000Processed variant 9000Processed variant 10000Processed variant: 10000
#> All variants processed

# read individuals to include
indivs <- read.table(paste0(data_path, "indivs_b"), header = TRUE)$id %>%
  as.character()

# check if all indivs are in vcf_names
if (any(!(indivs %in% colnames(vcf@gt)[-1]))) stop(paste("Some individuals in list not in VCF"))

# get groups and group info on samples in vcf
# default file names for strata is "strata" and for groups "groups"
project_info <- get_vcf_group_info(vcf, data_path)
#> Joining with `by = join_by(id)`
```

### Assessment of Missing Data

The function
[`assess_vcf_missing_data()`](https://legallab.github.io/vcf2others/reference/assess_vcf_missing_data.md)
is used to generate a table and graph of missing data for each sample in
a VCF. It is useful for visualizing data before and after filtering, and
evaluate the effect of filtering. The function accepts parameters which
form part of the name of the input file. The idea is that the VCF name
contains information on the name of the project (project), how the VCF
was extracted (postfix), and how it was filtered (fltr). The “postfix”
and “fltr” can be left blank if the VCF name does not contain this
information. Samples are ordered and colored by group assignment.

``` r
# define results path
res_path <- data_path
# species for plot title
species <- "Paleosuchus/Caiman"
# filter - filter parameters in file name
fltr <- ""

assess_vcf_missing_data(vcf, data_path, res_path, project, postfix, fltr, species)
```

### Basic stats for individuals in VCF

The function
[`vcf_stats()`](https://legallab.github.io/vcf2others/reference/vcf_stats.md)
is used to generate a table with basic stats for each sample in the VCF.
These stats include average read depth, heterozygosity, number of
heterozygotes, homozygotes, reference homozygotes and alternate
homogygotes, percent missing SNPs, total missing SNPs and total SNPs.

``` r
# get a table of basic sample stats

vcf_stats(vcf, res_path, project)
```

### Filtering, subsetting, merging and otherwise wrangling VCF files

All functions can be used directly on a VCF, but many of them are called
by the other functions of this package. For example
[`vcf_filter_invariant()`](https://legallab.github.io/vcf2others/reference/vcf_filter_invariant.md)
will remove invariant SNPs from a VCF, however, VCF by definition should
not have invariant SNPs. However, if a VCF is subsetted or certain
individuals are removed, SNPs can become invariant and need to be
removed.

#### Subsetting

The functions
[`vcf_extract_indivs()`](https://legallab.github.io/vcf2others/reference/vcf_extract_indivs.md),
[`vcf_extract_pops()`](https://legallab.github.io/vcf2others/reference/vcf_extract_pops.md)
and
[`vcf_sub_loci()`](https://legallab.github.io/vcf2others/reference/vcf_sub_loci.md)
are used to subset a VCF by specific individuals, a group (all
individuals of a group/population), or a random subset of loci,
respectively. Subsetting individuals and groups can be done either as a
whitelist (keep the listed individuals/groups) or blacklist (remove the
listed individuals/groups).

``` r
# subset VCF by individuals (default whitelist = TRUE, remove invariants = TRUE)
indivs1 <- c("CTGA_H4635", "CTGA_H4667")
vcf1 <- vcf_extract_indivs(vcf, indivs1)

# subset VCF by groups (default whitelist = TRUE, remove invariants = TRUE)
groups1 <- c("paleosuchus_1", "paleosuchus_2") %>%
  as.factor()
vcf1 <- vcf_extract_pops(vcf, project_info$indiv_group, groups1)

# subset VCF by n random loci (default n_loci = 1000)
vcf1 <- vcf_sub_loci(vcf, n_loci = 5000)
```

#### Filtering

Filtering consists of removing individuals and loci based on specific
properties. The function
[`vcf_filter_missing_indivs()`](https://legallab.github.io/vcf2others/reference/vcf_filter_missing_indivs.md)
removes individuals above a % missing SNP threshold. The other functions
filter on properties of individuals SNPs across individuals. The
function
[`vcf_filter_missingness()`](https://legallab.github.io/vcf2others/reference/vcf_filter_missingness.md)
removes SNPs above a % missing threshold,
[`vcf_filter_quality()`](https://legallab.github.io/vcf2others/reference/vcf_filter_quality.md)
removes SNPs below a quality threshold,
[`vcf_filter_maf()`](https://legallab.github.io/vcf2others/reference/vcf_filter_maf.md)
removes SNPs with minor allele frequency below a % threshold,
[`vcf_filter_hets()`](https://legallab.github.io/vcf2others/reference/vcf_filter_hets.md)
removes SNPs with heterozigosity above a % threshold,
[`vcf_filter_rank()`](https://legallab.github.io/vcf2others/reference/vcf_filter_rank.md)
removes SNPs with rank below a % threshold (rank is specific to VCFs
generated by DiscoSnpRad and measures the likelihood of observed SNP
variation being due to paralogs), and
[`vcf_filter_coverage()`](https://legallab.github.io/vcf2others/reference/vcf_filter_coverage.md)
removes SNPs below a specific coverage (makes them missing in the
individual). Finally the
[`vcf_filter_invariant()`](https://legallab.github.io/vcf2others/reference/vcf_filter_invariant.md)
removes invariant SNPs (SNPs may become invariant after subsetting and
filtering). The functions
[`vcf_filter_oneSNP()`](https://legallab.github.io/vcf2others/reference/vcf_filter_oneSNP.md)
and
[`vcf_filter_multiSNP()`](https://legallab.github.io/vcf2others/reference/vcf_filter_multiSNP.md)
filter a VCF to remove linked SNPs (retain only unliked SNPs) or retain
only linked SNPs, respectively. In general, all analyses assume that
SNPs are unlinked, however, some analyses such as FineRadStructure work
with linked SNPs.

``` r
# filter VCF for analyses (unlinked SNPs)
vcf_oneSNP <- vcf_extract_indivs(vcf, indivs, whitelist = FALSE) %>%
  vcf_filter_missing_indivs(.8) %>%
  vcf_filter_rank(.5) %>%
  vcf_filter_maf(.03) %>%
  vcf_filter_coverage(6) %>%
  vcf_filter_oneSNP() %>%
  vcf_filter_missingness(.2) %>%
  vcf_filter_missing_indivs(.3)
#> Removed samples are:
#> CTGA_H4663
#> CTGA_H4668
#> Final % missing data in VCF is 55.12%
#> Final % missing data in VCF is 13.92%
#> Removed samples are:
#> CTGA_H4666
#> CTGA_H4667
#> CTGA_H4669
#> Final % missing data in VCF is 7.1%

# filter VCF for analyses (linked SNPs)
vcf_multiSNP <- vcf_extract_indivs(vcf, indivs, whitelist = FALSE) %>%
  vcf_filter_missing_indivs(.8) %>%
  vcf_filter_rank(.5) %>%
  vcf_filter_maf(.03) %>%
  vcf_filter_coverage(6) %>%
  vcf_filter_multiSNP() %>%
  vcf_filter_missingness(.2) %>%
  vcf_filter_missing_indivs(.3)
#> Removed samples are:
#> CTGA_H4663
#> CTGA_H4668
#> Final % missing data in VCF is 55.12%
#> Warning: In dplyr::row_number(): 
#> ℹ Expression not supported in Arrow
#> → Pulling data into R
#> Final % missing data in VCF is 12.89%
#> Removed samples are:
#> CTGA_H4667
#> CTGA_H4669
#> Final % missing data in VCF is 8.27%
```

#### Extracting, merging and adding

The functions
[`vcf_merge()`](https://legallab.github.io/vcf2others/reference/vcf_merge.md),
[`vcf_add_indivs()`](https://legallab.github.io/vcf2others/reference/vcf_add_indivs.md),
[`vcf_extract_indivs()`](https://legallab.github.io/vcf2others/reference/vcf_extract_indivs.md)
and
[`vcf_extract_pops()`](https://legallab.github.io/vcf2others/reference/vcf_extract_pops.md)
are used to merge two VCF files, add individuals from a second VCF into
a first VCF, and extract individuals/groups from a VCF saving it into a
separate VCF. The
[`vcf_add_indivs()`](https://legallab.github.io/vcf2others/reference/vcf_add_indivs.md),
[`vcf_extract_indivs()`](https://legallab.github.io/vcf2others/reference/vcf_extract_indivs.md)
and
[`vcf_extract_pops()`](https://legallab.github.io/vcf2others/reference/vcf_extract_pops.md)
functions were designed in mind to work with outgroups. Because of
phylogenetic effects, outgroups have relatively more missing data than
ingroups, and so often are removed during filtering. The solution is to
extract the outgroup taxa with
[`vcf_extract_indivs()`](https://legallab.github.io/vcf2others/reference/vcf_extract_indivs.md)
or
[`vcf_extract_pops()`](https://legallab.github.io/vcf2others/reference/vcf_extract_pops.md),
filter the ingroup VCF, and add the outgroup taxa with
[`vcf_add_indivs()`](https://legallab.github.io/vcf2others/reference/vcf_add_indivs.md).
It is important to NOT filter invariant loci during the filtering of the
ingroup, and only filter them
[`vcf_filter_invariant()`](https://legallab.github.io/vcf2others/reference/vcf_filter_invariant.md)
after adding back the outgroup. The functions
[`vcf_extract_indivs()`](https://legallab.github.io/vcf2others/reference/vcf_extract_indivs.md),
[`vcf_extract_pops()`](https://legallab.github.io/vcf2others/reference/vcf_extract_pops.md),
[`vcf_filter_missing_indivs()`](https://legallab.github.io/vcf2others/reference/vcf_filter_missing_indivs.md)
and
[`vcf_filter_coverage()`](https://legallab.github.io/vcf2others/reference/vcf_filter_coverage.md)
by default filter invariants, but this can be turned off using the
“f_invar = FALSE” flag.

``` r
# extract individuals from VCF (keep all loci)
indivs1 <- c("CTGA_H4635", "CTGA_H4637", "CTGA_H4643", "CTGA_H4645", "CTGA_H4647")
vcf_outgrp <- vcf_extract_indivs(vcf, indivs1, f_invar = FALSE)
vcf_ingrp <- vcf_extract_indivs(vcf, indivs1, whitelist = FALSE, f_invar = FALSE)

# extract groups of individuals from VCF (keep all loci)
groups1 <- as.factor("caiman")
vcf_outgrp <- vcf_extract_pops(vcf, project_info$indiv_group, groups1, f_invar = FALSE)
vcf_ingrp <- vcf_extract_pops(vcf, project_info$indiv_group, groups1, whitelist = FALSE, f_invar = FALSE)

# merge vcf_outgrp into vcf_ingrp
vcf1 <- vcf_merge(vcf_ingrp, vcf_outgrp)
#> Joining with `by = join_by(FORMAT, id)`

# add individuals from vcf_outgrp into vcf_ingrp (whitelist = TRUE is default)
indivs1 <- c("CTGA_H4637", "CTGA_H4643", "CTGA_H4645")
vcf1 <- vcf_add_indivs(vcf_ingrp, vcf_outgrp, indivs1)
#> Joining with `by = join_by(FORMAT, id)`
```

#### Filtering with outgroup

Outgroups, which tend to have many fewer individuals than ingroups, has
more missing data than ingroups due to phylogenetic effects. Filtering
the entire dataset often results in removal of the outgroup taxa when
[`vcf_filter_invariant()`](https://legallab.github.io/vcf2others/reference/vcf_filter_invariant.md)
is used during filtering to remove individuals with % missing data above
some threshold. Therefore, a better strategy may be to extract the
outgroup, filter the ingroup, add the outgroup, and remove any invariant
SNPs.

``` r
# extract outgroup from VCF (keep all loci)
groups1 <- as.factor("caiman")
vcf_outgrp <- vcf_extract_pops(vcf, project_info$indiv_group, groups1, f_invar = FALSE)

# filter ingroup VCF for analyses then add outgroups
vcf1 <- vcf_extract_pops(vcf, project_info$indiv_group, groups1, whitelist = FALSE, f_invar = FALSE) %>%
  vcf_filter_missing_indivs(.8, f_invar = FALSE) %>%
  vcf_filter_rank(.5) %>%
  vcf_filter_maf(.03) %>%
  vcf_filter_coverage(6) %>%
  vcf_filter_oneSNP() %>%
  vcf_filter_missingness(.2) %>%
  vcf_filter_missing_indivs(.3, f_invar = FALSE) %>%
  vcf_merge(vcf_outgrp) %>%
  vcf_filter_invariant()
#> Removed samples are:
#> CTGA_H4663
#> CTGA_H4668
#> Final % missing data in VCF is 57.65%
#> Final % missing data in VCF is 7.5%
#> Removed samples are:
#> CTGA_H4667
#> Final % missing data in VCF is 2.78%
#> Joining with `by = join_by(FORMAT, id)`
```

### Converting a VCF file to other population genetic and phylogenetic formats

All the following functions will take a VCF, groups and indiv_group
(individuals to groups relationship), and convert the VCF to other
populations genetic and phylogenetic formats, exporting/writing a file
of this format. By default all functions will include missing data. If
that is not desired, only 100% complete data matrices can be exported
using the “inc_missing = FALSE” flag. The function
[`vcf2genlight()`](https://legallab.github.io/vcf2others/reference/vcf2genlight.md)
automatically returns a genlight object and optionally can also
export/write a file of this format to the working directory. Generally
the
[`vcf2genlight()`](https://legallab.github.io/vcf2others/reference/vcf2genlight.md)
function is called within a script using functions of the `adegenet`
package rather than importing the genlight object.

``` r
res_path <- data_path
project <- "trigonatus_"
postfix <- "discosnp_sub"

##########
# datasets for analyses with unlinked SNPs
vcf <- vcf_oneSNP
# get groups and group info on samples in vcf
project_info <- get_vcf_group_info(vcf, data_path)
#> Joining with `by = join_by(id)`

##########
# export data formats
# migrate-n https://peterbeerli.com/migrate-html5/
vcf2migrate(vcf, project_info$indiv_group, project_info$groups, out_file = paste0(res_path, project, postfix, '_migrate.txt'))
# arlequin http://cmpg.unibe.ch/software/arlequin35/
vcf2arlequin(vcf, project_info$indiv_group, project_info$groups, out_file = paste0(res_path, project, postfix, '.arp'))
# structure https://web.stanford.edu/group/pritchardlab/structure.html
vcf2structure(vcf, project_info$indiv_group, project_info$groups, out_file = paste0(res_path, project, postfix, '.str'))
# faststucture http://rajanil.github.io/fastStructure/
vcf2structure(vcf, project_info$indiv_group, project_info$groups, out_file = paste0(res_path, project, postfix, '.fstr'), method = "F")
# genepop https://gitlab.mbb.univ-montp2.fr/francois/genepop
vcf2genepop(vcf, project_info$indiv_group, project_info$groups, out_file = paste0(res_path, project, postfix, '.gen'))
# smartsnp https://github.com/ChristianHuber/smartsnp
vcf2smartsnp(vcf, project_info$indiv_group, project_info$groups, out_file = paste0(res_path, project, postfix, '.smartsnp'))
# eigenstrat https://github.com/DReichLab/EIG/tree/master
vcf2eigenstrat(vcf, project_info$indiv_group, project_info$groups, out_file = paste0(res_path, project, postfix, '_eigenstrat'))
# bayescan https://github.com/mfoll/BayeScan
vcf2bayescan(vcf, project_info$indiv_group, project_info$groups, out_file = paste0(res_path, project, postfix, '.bayescan'))
# treemix https://bitbucket.org/nygcresearch/treemix/wiki/Home
vcf2treemix(vcf, project_info$indiv_group, project_info$groups, out_file = paste0(res_path, project, postfix, '.treemix'))
# apparent https://github.com/halelab/apparent/tree/master
vcf2apparent(vcf, project_info$indiv_group, project_info$groups, out_file = paste0(res_path, project, postfix, '.apparent'))
# related https://github.com/timothyfrasier/related
vcf2related(vcf, project_info$indiv_group, project_info$groups, out_file = paste0(res_path, project, postfix, '.related'))
# dataframe with population membership in the last column
vcf2df(vcf, project_info$indiv_group, project_info$groups, out_file = paste0(res_path, project, postfix, '.df'))
# snapp https://www.beast2.org/snapp/
vcf2snapp(vcf, project_info$indiv_group, project_info$groups, out_file = paste0(res_path, project, postfix, '_snapp.nex'))
# nexus - only SNPs, meant for SDVq analyses https://www.asc.ohio-state.edu/kubatko.2/software/SVDquartets/
vcf2nexus(vcf, project_info$indiv_group, project_info$groups, out_file = paste0(res_path, project, postfix, '_sdvq.nex'))
# fasta https://www.ncbi.nlm.nih.gov/genbank/fastaformat/
vcf2fasta(vcf, project_info$indiv_group, project_info$groups, out_file = paste0(res_path, project, postfix, '.fna'))

##########
# datasets for analyses with linked SNPs
vcf <- vcf_multiSNP
# get groups and group info on samples in vcf
project_info <- get_vcf_group_info(vcf, data_path)
#> Joining with `by = join_by(id)`

##########
# export data formats
# fineRadStructure - expects VCF of linked SNPs
vcf2fineRadStructure(vcf, project_info$indiv_group, project_info$groups, out_file = paste0(res_path, project, postfix, '.finerad'))
```
