# vcf2genlight

converts vcfR format data to Genlight infile

in part based on vcfR2genlight function (vcfR package)

## Usage

``` r
vcf2genlight(
  vcf,
  ind_pop,
  keep_pop,
  ploidy = 2,
  inc_missing = TRUE,
  save = FALSE,
  out_file = "genlight_infile"
)
```

## Arguments

- vcf:

  -\> vcfR object

- ind_pop:

  -\> population assignment of individuals in vcf (factor)

- keep_pop:

  -\> population(s) of interest to include in Genlight infile (factor)

- ploidy:

  -\> ploidy level (default 2)

- inc_missing:

  -\> include missing data (logical)

- save:

  -\> save to file (logical - default TRUE) or just return Genlight
  object

## Value

Genlight object

## Details

This function converts the vcfR object to a Genlight formatted input
file This function labels populations. The function will remove indels,
and multiallelic loci, and optionally loci with missing data

## Author

Tomas Hrbek December 2022

## Examples

``` r
vcf2genlight(vcf = my_vcf, ind_pop = ind_pop, keep_pop = keepers, ploidy = 2, inc_missing = TRUE, out_file = "Genlight_infile")
#> Error: object 'my_vcf' not found
vcf2genlight(my_vcf, ind_pop, keepers, out_file = "Genlight_infile")
#> Error: object 'my_vcf' not found
vcf2genlight(my_vcf, ind_pop, keepers)
#> Error: object 'my_vcf' not found
```
