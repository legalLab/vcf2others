# vcf_filter_multiSNP

subsets vcfR format data keeping only loci with 2+ SNPs

## Usage

``` r
vcf_filter_multiSNP(vcf, minS = 2, maxS = 5)
```

## Arguments

- vcf:

  -\> vcfR object

## Value

subsetted vcfR object

## Details

This function subsets the vcfR object keeping only loci with between min
and max \# of SNPs per locus, returning new vcfR object default min = 2
and max = 5 SNPs per locus Recommended as input for fineRADstructure
analyses

## Author

Tomas Hrbek February 2022

## Examples

``` r
vcf_filter_multiSNP(vcf = my_vcf, minS = 2, maxS = 5)
#> Error in h(simpleError(msg, call)): error in evaluating the argument 'x' in selecting a method for function 'getCHROM': object 'my_vcf' not found
vcf_filter_multiSNP(vcf = my_vcf)
#> Error in h(simpleError(msg, call)): error in evaluating the argument 'x' in selecting a method for function 'getCHROM': object 'my_vcf' not found
vcf_filter_multiSNP(my_vcf)
#> Error in h(simpleError(msg, call)): error in evaluating the argument 'x' in selecting a method for function 'getCHROM': object 'my_vcf' not found
```
