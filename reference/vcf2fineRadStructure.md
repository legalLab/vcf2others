# vcf2fineRadStructure

converts vcfR format data to fineRadStructure infile

in part based on vcfR2genepop function

## Usage

``` r
vcf2fineRadStructure(
  vcf,
  ind_pop,
  keep_pop,
  inc_missing = TRUE,
  out_file = "fineRadStructure_infile.txt"
)
```

## Arguments

- vcf:

  -\> vcfR object

- ind_pop:

  -\> population assignment of individuals in vcf (factor)

- keep_pop:

  -\> population(s) of interest to include in fineRadStructure infile
  (factor)

- inc_missing:

  -\> include missing data (logical)

- out_file:

  -\> name of file to output (fineRadStructure infile)

## Value

nothing

## Details

This function converts the vcfR object to a fineRadStructure formatted
input file The function expects multiallelic loci

## Author

Tomas Hrbek July 2022

## Examples

``` r
vcf2fineRadStructure(vcf = my_vcf, ind_pop = ind_pop, keep_pop = keepers, inc_missing = TRUE, out_file = "fineRadStructure_infile.txt")
#> Error: object 'my_vcf' not found
vcf2fineRadStructure(my_vcf, ind_pop, keepers, out_file = "fineRadStructure_infile.txt")
#> Error: object 'my_vcf' not found
vcf2fineRadStructure(my_vcf, ind_pop, keepers)
#> Error: object 'my_vcf' not found
```
