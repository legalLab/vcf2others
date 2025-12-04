# vcf2snapp

converts vcfR format data to SNAPP (Nexus) format infile

## Usage

``` r
vcf2snapp(vcf, ind_pop, keep_pop, inc_missing = TRUE, out_file = "snapp.nex")
```

## Arguments

- vcf:

  -\> vcfR object

- ind_pop:

  -\> population assignment of individuals in vcf (factor)

- keep_pop:

  -\> population(s) of interest to include in SNAPP infile (factor)

- inc_missing:

  -\> include missing data (logical)

- out_file:

  -\> name of file to output (SNAPP infile)

## Value

nothing

## Details

This function converts the vcfR object to a SNAPP (Nexus) formatted
input file The function will remove indels, and multiallelic loci, and
optionally loci with missing data

## Author

Tomas Hrbek December 2020

## Examples

``` r
vcf2snapp(vcf = my_vcf, ind_pop = ind_pop, keep_pop = keepers, inc_missing = TRUE, out_file = "SNAPP_infile.nex")
#> Error: object 'my_vcf' not found
vcf2snapp(my_vcf, ind_pop, keepers, out_file = "SNAPP_infile.nex")
#> Error: object 'my_vcf' not found
vcf2snapp(my_vcf, ind_pop, keepers)
#> Error: object 'my_vcf' not found
```
