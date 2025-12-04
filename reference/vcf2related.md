# vcf2related

converts vcfR format data to Related infile

in part based on vcfR2genepop function (vcfR package)

## Usage

``` r
vcf2related(
  vcf,
  ind_pop,
  keep_pop,
  inc_missing = TRUE,
  out_file = "related_infile.txt"
)
```

## Arguments

- vcf:

  -\> vcfR object

- ind_pop:

  -\> population assignment of individuals in vcf (factor)

- keep_pop:

  -\> population(s) of interest to include in Related infile (factor)

- inc_missing:

  -\> include missing data (logical)

- out_file:

  -\> name of file to output (Related infile)

## Value

nothing

## Details

This function converts the vcfR object to a Related formatted input file
The function will remove indels, and multiallelic loci, and optionally
loci with missing data 'A', 'C', 'G', 'T' is 01, 02, 03, 04

## Author

Tomas Hrbek April 2023

## Examples

``` r
vcf2related(vcf = my_vcf, ind_pop = ind_pop, keep_pop = keepers, inc_missing = TRUE, out_file = "Genepop_infile.gen")
#> Error: object 'my_vcf' not found
vcf2related(my_vcf, ind_pop, keepers, out_file = "Related_infile.txt")
#> Error: object 'my_vcf' not found
vcf2related(my_vcf, ind_pop, keepers)
#> Error: object 'my_vcf' not found
```
