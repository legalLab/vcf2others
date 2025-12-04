# vcf_filter_missing_indivs

remove from vcfR format data indivs with \>% missing data

## Usage

``` r
vcf_filter_missing_indivs(vcf, p_miss, f_invar = TRUE)
```

## Arguments

- vcf:

  -\> vcfR object

- p_miss:

  -\> max missing data per locus as decimal (numeric)

## Value

subsetted vcfR object

## Details

This function subsets the vcfR object by % missing data within and
individual, returning new vcfR object

## Author

Tomas Hrbek October 2022

## Examples

``` r
vcf_filter_missing_indivs(vcf = my_vcf, p_miss = p_miss)
#> Error: object 'my_vcf' not found
vcf_filter_missing_indivs(my_vcf, p_miss)
#> Error: object 'my_vcf' not found
```
