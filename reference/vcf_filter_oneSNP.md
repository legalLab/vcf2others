# vcf_filter_oneSNP

subsets vcfR format data keeping only 1 SNP per locus

## Usage

``` r
vcf_filter_oneSNP(vcf)
```

## Arguments

- vcf:

  -\> vcfR object

## Value

subsetted vcfR object

## Details

This function subsets the vcfR object keeping only 1 SNP per locus,
returning new vcfR object The first SNP independent of quality is taken
(may mofify this in future)

## Author

Tomas Hrbek February 2022

## Examples

``` r
vcf_filter_oneSNP(vcf = my_vcf)
#> Error in h(simpleError(msg, call)): error in evaluating the argument 'x' in selecting a method for function 'getCHROM': object 'my_vcf' not found
vcf_filter_oneSNP(my_vcf)
#> Error in h(simpleError(msg, call)): error in evaluating the argument 'x' in selecting a method for function 'getCHROM': object 'my_vcf' not found
```
