# vcf_filter_rank

remove loci below rank threshold from vcfR format data

rank is calculated as sqrt(chi-sqr/n) of allele read counts

used for paralog detection - very low rank values (\<0.4)

rank is calculated in DiscoSNP, registered as Pk in INFO

## Usage

``` r
vcf_filter_rank(vcf, rank = 0.4)
```

## Arguments

- vcf:

  -\> vcfR object

## Value

subsetted vcfR object

## Details

This function removes loci below rank threshold from the vcfR object

## Author

Tomas Hrbek February 2022

## Examples

``` r
vcf_filter_quality(vcf = my_vcf, rank = .4)
#> Error in vcf_filter_quality(vcf = my_vcf, rank = 0.4): unused argument (rank = 0.4)
vcf_filter_quality(my_vcf)
#> Error in h(simpleError(msg, call)): error in evaluating the argument 'x' in selecting a method for function 'getQUAL': object 'my_vcf' not found
```
