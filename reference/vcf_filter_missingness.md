# vcf_filter_missingness

subsets vcfR format data by % missing data

## Usage

``` r
vcf_filter_missingness(vcf, p_miss)
```

## Arguments

- vcf:

  -\> vcfR object

- p_miss:

  -\> max missing data per locus as decimal (numeric)

## Value

subsetted vcfR object

## Details

This function subsets the vcfR object by % missing data, returning new
vcfR object

## Author

Tomas Hrbek February 2022

## Examples

``` r
vcf_filter_missingness(vcf = my_vcf, p_miss = p_miss)
#> Error: object 'my_vcf' not found
vcf_filter_missingness(my_vcf, p_miss)
#> Error: object 'my_vcf' not found
```
