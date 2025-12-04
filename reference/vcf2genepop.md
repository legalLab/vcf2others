# vcf2genepop

converts vcfR format data to Genepop infile

includes population names that can be ready by read_genepop

in part based on vcfR2migrate function (vcfR package)

## Usage

``` r
vcf2genepop(
  vcf,
  ind_pop,
  keep_pop,
  inc_missing = TRUE,
  out_file = "genepop_infile.txt"
)
```

## Arguments

- vcf:

  -\> vcfR object

- ind_pop:

  -\> population assignment of individuals in vcf (factor)

- keep_pop:

  -\> population(s) of interest to include in Genepop infile (factor)

- inc_missing:

  -\> include missing data (logical)

- out_file:

  -\> name of file to output (Genepop infile)

## Value

nothing

## Details

This function converts the vcfR object to a Genepop formatted input file
This function labels populations. To read labeled populations use
"read_genepop" function The function will remove indels, and
multiallelic loci, and optionally loci with missing data 01, 02, 03, 04
is 'A', 'C', 'G', 'T'

## Author

Tomas Hrbek January 2021

## Examples

``` r
vcf2genepop(vcf = my_vcf, ind_pop = ind_pop, keep_pop = keepers, inc_missing = TRUE, out_file = "Genepop_infile.gen")
#> Error: object 'my_vcf' not found
vcf2genepop(my_vcf, ind_pop, keepers, out_file = "Genepop_infile.gen")
#> Error: object 'my_vcf' not found
vcf2genepop(my_vcf, ind_pop, keepers)
#> Error: object 'my_vcf' not found
```
