# vcf_extract_pops

extract vcfR format data by group

## Usage

``` r
vcf_extract_pops(vcf, ind_pop, keep_pop, whitelist = TRUE, f_invar = TRUE)
```

## Arguments

- vcf:

  -\> vcfR object

- ind_pop:

  -\> group assignment of individuals in vcf (factor)

- keep_pop:

  -\> group(s) of interest to include/exclude in vcf infile (factor)

## Value

extracted vcfR object

## Details

This function extracts the vcfR object by group(s) (individuals assigned
to one or more group), returning new vcfR object

## Author

Tomas Hrbek December 2020

## Examples

``` r
vcf_extract_pops(vcf = my_vcf, ind_pop = ind_pop, keep_pop = keepers, whitelist = TRUE)
#> Error: object 'ind_pop' not found
vcf_extract_pops(vcf = my_vcf, ind_pop = ind_pop, keep_pop = non_keepers, whitelist = FALSE)
#> Error: object 'ind_pop' not found
vcf_extract_pops(my_vcf, ind_pop, keepers)
#> Error: object 'ind_pop' not found
```
