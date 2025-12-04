# vcf_extract_indivs

extract vcfR format data by individuals

## Usage

``` r
vcf_extract_indivs(vcf, indiv, whitelist = TRUE, f_invar = TRUE)
```

## Arguments

- vcf:

  -\> vcfR object

- indiv:

  -\> individuals to retain/drop in vcf (factor)

## Value

extracted vcfR object

## Details

This function extracts the vcfR object by individuals, returning new
vcfR object

## Author

Tomas Hrbek December 2020

## Examples

``` r
vcf_extract_indivs(vcf = my_vcf, indiv = indivs_to_keep, whitelist = TRUE)
#> Error: object 'my_vcf' not found
vcf_extract_indivs(vcf = my_vcf, indiv = indivs_to_drop, whitelist = FALSE)
#> Error: object 'my_vcf' not found
vcf_extract_indivs(my_vcf, indivs_to_keep)
#> Error: object 'my_vcf' not found
```
