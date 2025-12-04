# vcf2treemix

converts vcfR format data to TreeMix infile

in part based on vcfR2migrate function (vcfR package)

## Usage

``` r
vcf2treemix(
  vcf,
  ind_pop,
  keep_pop,
  inc_missing = TRUE,
  out_file = "treemix_infile.txt"
)
```

## Arguments

- vcf:

  -\> vcfR object

- ind_pop:

  -\> population assignment of individuals in vcf (factor)

- keep_pop:

  -\> population(s) of interest to include in TreeMix infile (factor)

- inc_missing:

  -\> include missing data (logical)

- out_file:

  -\> name of file to output (TreeMix infile)

## Value

nothing

## Details

This function converts the vcfR object to a TreeMix formatted input file
The function will remove indels, and multiallelic loci, and optionally
loci with missing data

## Author

Tomas Hrbek November 2022

## Examples

``` r
vcf2treemix(vcf = my_vcf, ind_pop = ind_pop, keep_pop = keepers, inc_missing = TRUE, out_file = "TreeMix_infile.txt")
#> Error: object 'my_vcf' not found
vcf2treemix(my_vcf, ind_pop, keepers, out_file = "TreeMix_infile.txt")
#> Error: object 'my_vcf' not found
vcf2treemix(my_vcf, ind_pop, keepers)
#> Error: object 'my_vcf' not found
```
