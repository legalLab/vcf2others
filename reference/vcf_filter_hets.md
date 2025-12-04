# vcf_filter_hets

remove loci above het threshold from vcfR format data

high hets are indicative of paralogs

## Usage

``` r
vcf_filter_hets(vcf, hets = 0.4)
```

## Arguments

- vcf:

  -\> vcfR object

## Value

subsetted vcfR object

## Details

This function removes loci above het threshold from the vcfR object

## Author

Tomas Hrbek October 2022

## Examples

``` r
vcf_filter_hets(vcf = my_vcf, hets = .4)
#> Error: object 'my_vcf' not found
vcf_filter_hets(my_vcf)
#> Error: object 'my_vcf' not found
```
