# vcf2fasta

converts vcfR format data to Fasta format infile

## Usage

``` r
vcf2fasta(
  vcf,
  ind_pop,
  keep_pop,
  interleaved = FALSE,
  inc_missing = TRUE,
  out_file = "fasta.fas"
)
```

## Arguments

- vcf:

  -\> vcfR object

- ind_pop:

  -\> population assignment of individuals in vcf (factor)

- keep_pop:

  -\> population(s) of interest to include in Fasta infile (factor)

- inc_missing:

  -\> include missing data (logical)

- out_file:

  -\> name of file to output (Fasta infile)

## Value

nothing

## Details

This function converts the vcfR object to a Fasta formatted input file
The function will remove indels, and multiallelic loci, and optionally
loci with missing data

## Author

Tomas Hrbek August 2022

## Examples

``` r
vcf2fasta(vcf = my_vcf, ind_pop = ind_pop, keep_pop = keepers, interleaved = FALSE, inc_missing = TRUE, out_file = "Fasta_infile.fas")
#> Error: object 'my_vcf' not found
vcf2fasta(my_vcf, ind_pop, keepers, out_file = "Fasta_infile.fas")
#> Error: object 'my_vcf' not found
vcf2fasta(my_vcf, ind_pop, keepers)
#> Error: object 'my_vcf' not found
```
