# assess_vcf_missing_data

quantifying missing data of all samples in VCF

inspired by
https://grunwaldlab.github.io/Population_Genetics_in_R/qc.html

## Usage

``` r
assess_vcf_missing_data(
  vcf,
  data_path,
  res_path,
  project,
  postfix,
  fltr,
  species,
  strt = "strata"
)
```

## Arguments

- vcf:

  -\> vcfR object

- data_path:

  -\> path to data (where strata is located)

- res_path:

  -\> path to results (directory for output dataframe and plot)

- project:

  -\> project base name (character)

- postfix:

  -\> VCF extraction method (character)

- fltr:

  -\> VCF filtration method (character)

- species:

  -\> species name for plot (character)

- strt:

  -\> file with 2+ columns, one id column and one pop column (tsv file,
  header = TRUE)

## Value

dataframe and plot of missing data for each individual in VCF

## Details

This function generates a dataframe of absolute and relative missing
data per sample This function plots relative missing data per
individual, with individuals sorted by group The postfix and fltr
parameters can have empty values ("")

## Author

Tomas Hrbek July 2022

## Examples

``` r
assess_vcf_missing_data(vcf = my_vcf, data_path = data_path, res_path = res_path, project = project_name, postfix = extraction_info, fltr = filter_parms, species = species_name, strt = "strata")
#> Error: object 'data_path' not found
assess_vcf_missing_data(my_vcf, data_path, res_path, project_name, extraction_info, filter_parms, species_name, strt = "strata")
#> Error: object 'data_path' not found
assess_vcf_missing_data(my_vcf, data_path, res_path, project_name, postfix, fltr, species)
#> Error: object 'data_path' not found
```
