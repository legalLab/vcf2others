# vcf_filter_maf

remove loci below MAF threshold from vcfR format data

## Usage

``` r
vcf_filter_maf(vcf, maf = 0.05)
```

## Arguments

- vcf:

  -\> vcfR object

## Value

subsetted vcfR object

## Details

This function removes loci below MAF threshold from the vcfR object

## Author

Tomas Hrbek October 2022

## Examples

``` r
vcf_filter_maf(vcf = my_vcf, maf = .05)
#> Error: object 'my_vcf' not found
vcf_filter_maf(my_vcf)
#> Error: object 'my_vcf' not found
```
