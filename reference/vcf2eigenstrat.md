# vcf2eigenstrat

converts vcfR format data to Eigenstrat infiles

in part based on vcfR2migrate function (vcfR package)

## Usage

``` r
vcf2eigenstrat(
  vcf,
  ind_pop,
  keep_pop,
  sex = "U",
  rel_pos = 0,
  inc_missing = TRUE,
  out_file = "eigenstrat_infile"
)
```

## Arguments

- vcf:

  -\> vcfR object

- ind_pop:

  -\> population assignment of individuals in vcf (factor)

- keep_pop:

  -\> population(s) of interest to include in Eigenstrat infiles
  (factor)

- inc_missing:

  -\> include missing data (logical)

- out_file:

  -\> name of file to output (Eigenstrat infiles)

## Value

nothing

## Details

This function converts the vcfR object to a Eigenstrat formatted input
files When list of sexes is not provided, lists all individuals as
unknown When relative position on chromosome (cM distance or similar) is
not provides, list as 0 The function will remove indels, and
multiallelic loci, and optionally loci with missing data

## Author

Tomas Hrbek December 2022

## Examples

``` r
vcf2eigenstrat(vcf = my_vcf, ind_pop = ind_pop, keep_pop = keepers, sex = list_of_sex, rel_pos = marker_cM_map, inc_missing = TRUE, out_file = "Eigenstrat")
#> Error: object 'my_vcf' not found
vcf2eigenstrat(my_vcf, ind_pop, keepers, out_file = "Eigenstrat")
#> Error: object 'my_vcf' not found
vcf2eigenstrat(my_vcf, ind_pop, keepers)
#> Error: object 'my_vcf' not found
```
