# vcf2df

converts vcfR format data to DataFrame infile

adds last column of an individual's group origin

motivated to provide a dataframe for plotting PCA results

for use by modified PCA function from adegenet

## Usage

``` r
vcf2df(
  vcf,
  ind_pop,
  keep_pop,
  alleles = TRUE,
  inc_missing = TRUE,
  out_file = "dataframe_infile.txt"
)
```

## Arguments

- vcf:

  -\> vcfR object

- ind_pop:

  -\> population assignment of individuals in vcf (factor)

- keep_pop:

  -\> population(s) of interest to include in DataFrame infile (factor)

- alleles:

  -\> record data as alleles (GATC) or SNPs (012) (logical)

- inc_missing:

  -\> include missing data (logical)

## Value

nothing

## Details

This function converts the vcfR object to a DataFrame formatted input
file This function adds group (population) membership as last column The
function will remove indels, and multiallelic loci, and optionally loci
with missing data

## Author

Tomas Hrbek December 2022

## Examples

``` r
vcf2df(vcf = my_vcf, ind_pop = ind_pop, keep_pop = keepers, alleles = TRUE, inc_missing = TRUE, out_file = "dataframe_infile.txt")
#> Error: object 'my_vcf' not found
vcf2df(my_vcf, ind_pop, keepers, out_file = "dataframe_infile.txt")
#> Error: object 'my_vcf' not found
vcf2df(my_vcf, ind_pop, keepers)
#> Error: object 'my_vcf' not found
```
