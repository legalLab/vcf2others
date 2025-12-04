# vcf_stats

calculates basic stats of each samples from vcfR format data

high heterozygosity is indicative of potential contamination

average read depth per individual

missing data per individual

## Usage

``` r
vcf_stats(vcf, res_path, project)
```

## Arguments

- vcf:

  -\> vcfR object

## Value

table of statistics

## Details

This function calculates average read depth, heterozygosity number of
heterozygotes, number of reference and alternative homozygotes, missing
data and total number SNPs of each sample in an vcfR object

## Author

Tomas Hrbek August 2024

## Examples

``` r
vcf_stats(vcf = my_vcf)
#> Error: object 'my_vcf' not found
vcf_stats(my_vcf)
#> Error: object 'my_vcf' not found
```
