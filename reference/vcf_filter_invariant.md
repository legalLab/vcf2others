# vcf_filter_invariant

remove invariant loci from vcfR format data

## Usage

``` r
vcf_filter_invariant(vcf)
```

## Arguments

- vcf:

  -\> vcfR object

## Value

subsetted vcfR object

## Details

This function removes invariant loci from the vcfR object This might be
desirable after subsetting the vcf by individuals

## Author

Tomas Hrbek February 2022

## Examples

``` r
vcf_filter_invariant(vcf = my_vcf)
#> Error: object 'my_vcf' not found
vcf_filter_invariant(my_vcf)
#> Error: object 'my_vcf' not found
```
