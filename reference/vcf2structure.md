# vcf2structure

converts vcfR format data to Structure or FastStructure infile

in part based on vcfR2migrate function (vcfR package)

## Usage

``` r
vcf2structure(
  vcf,
  ind_pop,
  keep_pop,
  inc_missing = TRUE,
  out_file = "structure.str",
  method = "S"
)
```

## Arguments

- vcf:

  -\> vcfR object

- ind_pop:

  -\> population assignment of individuals in vcf (factor)

- keep_pop:

  -\> population(s) of interest to include in Structure infile (factor)

- inc_missing:

  -\> include missing data (logical)

- out_file:

  -\> name of file to output (Structure infile)

- method:

  -\> Structure or FastStructure format

## Value

nothing

## Details

This function converts the vcfR object to a Structure or FastStructure
formatted input file The function will remove indels, and multiallelic
loci, and optionally loci with missing data

## Author

Tomas Hrbek December 2020

## Examples

``` r
vcf2structure(vcf = my_vcf, ind_pop = ind_pop, keep_pop = keepers, inc_missing = TRUE, out_file = "Structure_infile.str", method = "S")
#> Error: object 'my_vcf' not found
vcf2structure(my_vcf, ind_pop, keepers, out_file = "Structure_infile.str")
#> Error: object 'my_vcf' not found
vcf2structure(my_vcf, ind_pop, keepers)
#> Error: object 'my_vcf' not found
```
