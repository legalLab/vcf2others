# vcf2migrate

converts vcfR format data to MigrateN infile

adapted vcfR2migrate function (vcfR package) to allow for inclusion of
missing data in Migrate output

converts vcfR format data to MigrateN infile

adapted vcfR2migrate function (vcfR package) to allow for inclusion of
missing data in Migrate output

## Usage

``` r
vcf2migrate(
  vcf,
  ind_pop,
  keep_pop,
  inc_missing = TRUE,
  out_file = "migrateN_infile.txt",
  method = "N"
)

vcf2migrate(
  vcf,
  ind_pop,
  keep_pop,
  inc_missing = TRUE,
  out_file = "migrateN_infile.txt",
  method = "N"
)
```

## Arguments

- vcf:

  -\> vcfR object

- ind_pop:

  -\> population assignment of individuals in vcf (factor)

- keep_pop:

  -\> population(s) of interest to include in MigrateN infile (factor)

- inc_missing:

  -\> include missing data (logical)

- out_file:

  -\> name of file to output (MigrateN infile)

- method:

  -\> classic or het format

## Value

nothing

nothing

## Details

This function converts the vcfR object to a MigrateN formatted input
file The function will remove indels, and multiallelic loci, and
optionally loci with missing data

This function converts the vcfR object to a MigrateN formatted input
file The function will remove indels, and multiallelic loci, and
optionally loci with missing data

## Author

Tomas Hrbek August 2020

## Examples

``` r
vcf2migrate(vcf = my_vcf, ind_pop = ind_pop, keep_pop = keepers, inc_missing = TRUE, out_file = "MigrateN_infile.txt", method = "N")
#> Error: object 'my_vcf' not found
vcf2migrate(my_vcf, ind_pop, keepers, out_file = "MigrateN_infile.txt")
#> Error: object 'my_vcf' not found
vcf2migrate(my_vcf, ind_pop, keepers)
#> Error: object 'my_vcf' not found

vcf2migrate(vcf = my_vcf, ind_pop = ind_pop, keep_pop = keepers, inc_missing = TRUE, out_file = "MigrateN_infile.txt")
#> Error: object 'my_vcf' not found
vcf2migrate(my_vcf, ind_pop, keepers, out_file = "MigrateN_infile.txt")
#> Error: object 'my_vcf' not found
vcf2migrate(my_vcf, ind_pop, keepers)
#> Error: object 'my_vcf' not found
```
