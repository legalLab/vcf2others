# vcf2apparent

converts vcfR format data to Apparent infile

## Usage

``` r
vcf2apparent(
  vcf,
  ind_pop,
  keep_pop,
  key = "All",
  inc_missing = TRUE,
  out_file = "apparent_infile.txt"
)
```

## Arguments

- vcf:

  -\> vcfR object

- ind_pop:

  -\> population assignment of individuals in vcf (factor)

- keep_pop:

  -\> population(s) of interest to include in Apparent infile (factor)

- key:

  -\> relationship type (All, Pa, Mo, Fa, Off); default All (character)

- inc_missing:

  -\> include missing data (logical)

- out_file:

  -\> name of file to output (Apparent infile)

## Value

nothing

## Details

This function converts the vcfR object to a Apparent formatted input
file The function will remove indels, and multiallelic loci, and
optionally loci with missing data

## Author

Tomas Hrbek February 2023

## Examples

``` r
vcf2apparent(vcf = my_vcf, ind_pop = ind_pop, keep_pop = keepers, key = key, inc_missing = TRUE, out_file = "Genepop_infile.gen")
#> Error: object 'my_vcf' not found
vcf2apparent(my_vcf, ind_pop, keepers, out_file = "Apparent_infile.txt")
#> Error: object 'my_vcf' not found
vcf2apparent(my_vcf, ind_pop, keepers)
#> Error: object 'my_vcf' not found
```
