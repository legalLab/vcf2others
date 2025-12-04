# vcf2nexus

converts vcfR format data to Nexus format infile

## Usage

``` r
vcf2nexus(vcf, ind_pop, keep_pop, inc_missing = TRUE, out_file = "nexus.nex")
```

## Arguments

- vcf:

  -\> vcfR object

- ind_pop:

  -\> population assignment of individuals in vcf (factor)

- keep_pop:

  -\> population(s) of interest to include in Nexus infile (factor)

- inc_missing:

  -\> include missing data (logical)

- out_file:

  -\> name of file to output (Nexus infile)

## Value

nothing

## Details

This function converts the vcfR object to a Nexus formatted input file
The function will remove indels, and multiallelic loci, and optionally
loci with missing data

## Author

Tomas Hrbek February 2021

## Examples

``` r
vcf2nexus(vcf = my_vcf, ind_pop = ind_pop, keep_pop = keepers, inc_missing = TRUE, out_file = "Nexus_infile.nex")
#> Error: object 'my_vcf' not found
vcf2nexus(my_vcf, ind_pop, keepers, out_file = "Nexus_infile.nex")
#> Error: object 'my_vcf' not found
vcf2nexus(my_vcf, ind_pop, keepers)
#> Error: object 'my_vcf' not found
```
