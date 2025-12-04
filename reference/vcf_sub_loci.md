# vcf_sub_loci

randomly subsets vcfR format data to \# of loci

## Usage

``` r
vcf_sub_loci(vcf, n_loci = 1000)
```

## Arguments

- vcf:

  -\> vcfR object

- n_loci:

  -\> number of loci to subset (numeric)

## Value

subsetted vcfR object

## Details

This function subsets the vcfR object to specific \# of loci, returning
new vcfR object

## Author

Tomas Hrbek January 2023

## Examples

``` r
vcf_sub_loci(vcf = my_vcf, n_loci = n_loci)
#> Error: object 'my_vcf' not found
vcf_sub_loci(my_vcf, n_loci)
#> Error: object 'my_vcf' not found
```
