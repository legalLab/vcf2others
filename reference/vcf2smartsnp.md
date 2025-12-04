# vcf2smartsnp

converts vcfR format data to smartsnp tabular format infile

in part based on vcfR2genlight function (vcfR package)

## Usage

``` r
vcf2smartsnp(
  vcf,
  ind_pop,
  keep_pop,
  inc_missing = TRUE,
  out_file = "smartsnp_infile.txt"
)
```

## Arguments

- vcf:

  -\> vcfR object

- ind_pop:

  -\> population assignment of individuals in vcf (factor)

- keep_pop:

  -\> population(s) of interest to include in genotype table infile
  (factor)

- inc_missing:

  -\> include missing data (logical)

- out_file:

  -\> name of file to output (genotype table infile)

## Value

nothing

## Details

This function converts the vcfR object to a genotype table formatted
input file For use in smartsnp, but keeps sample names as column names
The function will remove indels, and multiallelic loci, and optionally
loci with missing data

## Author

Tomas Hrbek October 2023

## Examples

``` r
vcf2smartsnp(vcf = my_vcf, ind_pop = ind_pop, keep_pop = keepers, inc_missing = TRUE, out_file = "SmartSNP_infile.txt")
#> Error: object 'my_vcf' not found
vcf2smartsnp(my_vcf, ind_pop, keepers, out_file = "SmartSNP_infile.txt")
#> Error: object 'my_vcf' not found
vcf2smartsnp(my_vcf, ind_pop, keepers)
#> Error: object 'my_vcf' not found
```
