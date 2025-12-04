# vcf_filter_quality

remove loci below quality threshold from vcfR format data

## Usage

``` r
vcf_filter_quality(vcf, qual = 20)
```

## Arguments

- vcf:

  -\> vcfR object

## Value

subsetted vcfR object

## Details

This function removes loci below quality threshold from the vcfR object

## Author

Tomas Hrbek February 2022

## Examples

``` r
vcf_filter_quality(vcf = my_vcf, qual = 20)
#> Error in h(simpleError(msg, call)): error in evaluating the argument 'x' in selecting a method for function 'getQUAL': object 'my_vcf' not found
vcf_filter_quality(my_vcf)
#> Error in h(simpleError(msg, call)): error in evaluating the argument 'x' in selecting a method for function 'getQUAL': object 'my_vcf' not found
```
