# vcf_filter_coverage

remove genotypes below coverage threshold from vcfR format data

## Usage

``` r
vcf_filter_coverage(vcf, cover = 10, f_invar = TRUE)
```

## Arguments

- vcf:

  -\> vcfR object

## Value

subsetted vcfR object

## Details

This function removes genotypes below coverage threshold from the vcfR
object

## Author

Tomas Hrbek November 2022

## Examples

``` r
vcf_filter_coverage(vcf = my_vcf, cover = 10)
#> Error: object 'my_vcf' not found
vcf_filter_coverage(my_vcf)
#> Error: object 'my_vcf' not found
```
