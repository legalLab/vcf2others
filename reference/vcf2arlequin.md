# vcf2arlequin

converts vcfR format data to Arlequin infile

in part based on vcfR2migrate function (vcfR package)

## Usage

``` r
vcf2arlequin(
  vcf,
  ind_pop,
  keep_pop,
  inc_missing = TRUE,
  out_file = "arlequin.arp"
)
```

## Arguments

- vcf:

  -\> vcfR object

- ind_pop:

  -\> population assignment of individuals in vcf (factor)

- keep_pop:

  -\> population(s) of interest to include in Arlequin infile (factor)

- inc_missing:

  -\> include missing data (logical)

- out_file:

  -\> name of file to output (Arlequin infile)

## Value

nothing

## Details

This function converts the vcfR object to a Arlequin formatted input
file The function will remove indels, and multiallelic loci, and
optionally loci with missing data

## Author

Tomas Hrbek December 2020

## Examples

``` r
vcf2arlequin(vcf = my_vcf, ind_pop = ind_pop, keep_pop = keepers, inc_missing = TRUE, out_file = "Arlequin_infile.arp")
#> Error: object 'my_vcf' not found
vcf2arlequin(my_vcf, ind_pop, keepers, out_file = "Arlequin_infile.arp")
#> Error: object 'my_vcf' not found
vcf2arlequin(my_vcf, ind_pop, keepers)
#> Error: object 'my_vcf' not found
```
