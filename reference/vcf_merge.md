# vcf_merge

merge two vcf; merge second vcf into first one

## Usage

``` r
vcf_merge(vcf, vcf1)
```

## Arguments

- vcf:

  -\> vcfR object

- vcf1:

  -\> vcfR object

## Value

augmented vcfR object

## Details

This function adds all individuals from a second vcfR object into a
first vcfR object, returning new vcfR object

## Author

Tomas Hrbek October 2023

## Examples

``` r
vcf_merge(vcf = my_vcf, vcf1 = other_vcf)
#> Error: object 'my_vcf' not found
vcf_merge(my_vcf, other_vcf)
#> Error: object 'my_vcf' not found
```
