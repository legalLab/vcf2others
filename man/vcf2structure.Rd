% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/vcf2structure.R
\name{vcf2structure}
\alias{vcf2structure}
\title{vcf2structure}
\usage{
vcf2structure(
  vcf,
  ind_pop,
  keep_pop,
  inc_missing = TRUE,
  out_file = "structure.str",
  method = "S"
)
}
\arguments{
\item{vcf}{-> vcfR object}

\item{ind_pop}{-> population assignment of individuals in vcf (factor)}

\item{keep_pop}{-> population(s) of interest to include in Structure infile (factor)}

\item{inc_missing}{-> include missing data (logical)}

\item{out_file}{-> name of file to output (Structure infile)}

\item{method}{-> Structure or FastStructure format}
}
\value{
nothing
}
\description{
converts vcfR format data to Structure or FastStructure infile

in part based on vcfR2migrate function (vcfR package)
}
\details{
This function converts the vcfR object to a Structure or FastStructure formatted input file
The function will remove indels, and multiallelic loci, and optionally loci with missing data
}
\examples{
vcf2structure(vcf = my_vcf, ind_pop = ind_pop, keep_pop = keepers, inc_missing = TRUE, out_file = "Structure_infile.str", method = "S")
vcf2structure(my_vcf, ind_pop, keepers, out_file = "Structure_infile.str")
vcf2structure(my_vcf, ind_pop, keepers)

}
\author{
Tomas Hrbek December 2020
}
