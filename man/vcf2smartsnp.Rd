% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/vcf2smartsnp.R
\name{vcf2smartsnp}
\alias{vcf2smartsnp}
\title{vcf2smartsnp}
\usage{
vcf2smartsnp(
  vcf,
  ind_pop,
  keep_pop,
  inc_missing = TRUE,
  out_file = "smartsnp_infile.txt"
)
}
\arguments{
\item{vcf}{-> vcfR object}

\item{ind_pop}{-> population assignment of individuals in vcf (factor)}

\item{keep_pop}{-> population(s) of interest to include in genotype table infile (factor)}

\item{inc_missing}{-> include missing data (logical)}

\item{out_file}{-> name of file to output (genotype table infile)}
}
\value{
nothing
}
\description{
converts vcfR format data to smartsnp tabular format infile

in part based on vcfR2genlight function (vcfR package)
}
\details{
This function converts the vcfR object to a genotype table formatted input file
For use in smartsnp, but keeps sample names as column names
The function will remove indels, and multiallelic loci, and optionally loci with missing data
}
\examples{
vcf2smartsnp(vcf = my_vcf, ind_pop = ind_pop, keep_pop = keepers, inc_missing = TRUE, out_file = "SmartSNP_infile.txt")
vcf2smartsnp(my_vcf, ind_pop, keepers, out_file = "SmartSNP_infile.txt")
vcf2smartsnp(my_vcf, ind_pop, keepers)
}
\author{
Tomas Hrbek October 2023
}
