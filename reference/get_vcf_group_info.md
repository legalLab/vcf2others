# get_vcf_group_info

associate samples with group info based on "strata" file

"strata" is a tsv file with 2+ columns, with 2 columns name id and pop

## Usage

``` r
get_vcf_group_info(vcf, data_path, strt = "strata", grps = "groups")
```

## Arguments

- vcf:

  -\> vcfR object

- data_path:

  -\> path to data (where strata and groups are located)

- strt:

  -\> file with 2+ columns, one id column and one pop column (tsv file,
  header = TRUE)

- grps:

  -\> one column file listing groups to include, column name pop (header
  = TRUE)

## Value

list of factors of sample-to-group assignments and groups

## Details

This function creates list of factors of sample-to-group assignments and
groups based on a strata and group file

## Author

Tomas Hrbek July 2022

## Examples

``` r
get_vcf_group_info(vcf = my_vcf, data_path = data_path, strt = "strata", grps = "groups")
#> Error: object 'my_vcf' not found
get_vcf_group_info(my_vcf, data_path, strt = "strata", grps = "groups")
#> Error: object 'my_vcf' not found
get_vcf_group_info(my_vcf, data_path)
#> Error: object 'my_vcf' not found
```
