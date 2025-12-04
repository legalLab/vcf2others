# vcf_add_indivs

add individuals to vcf

## Usage

``` r
vcf_add_indivs(vcf, vcf1, indiv, whitelist = TRUE)
```

## Arguments

- vcf:

  -\> vcfR object

- vcf1:

  -\> vcfR object

- indiv:

  -\> individuals to add into vcf from vcf1 (factor)

## Value

augmented vcfR object

## Details

This function adds individuals from one vcfR object to another vcfR
object, returning new vcfR object

## Author

Tomas Hrbek October 2023

## Examples

``` r
vcf_add_indivs(vcf = my_vcf, vcf1 = other_vcf, indiv = indivs_to_add, whitelist = TRUE)
#> Error: object 'other_vcf' not found
vcf_add_indivs(vcf = my_vcf, vcf1 = other_vcf, indiv = indivs_to_not_add, whitelist = FALSE)
#> Error: object 'other_vcf' not found
vcf_add_indivs(my_vcf, other_vcf, indivs_to_add)
#> Error: object 'other_vcf' not found
```
