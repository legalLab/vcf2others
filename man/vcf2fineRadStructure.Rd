% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/vcf2fineRadStructure.R
\name{vcf2fineRadStructure}
\alias{vcf2fineRadStructure}
\title{vcf2fineRadStructure}
\usage{
vcf2fineRadStructure(
  vcf,
  ind_pop,
  keep_pop,
  inc_missing = TRUE,
  out_file = "fineRadStructure_infile.txt"
)
}
\arguments{
\item{vcf}{-> vcfR object}

\item{ind_pop}{-> population assignment of individuals in vcf (factor)}

\item{keep_pop}{-> population(s) of interest to include in fineRadStructure infile (factor)}

\item{inc_missing}{-> include missing data (logical)}

\item{out_file}{-> name of file to output (fineRadStructure infile)}
}
\value{
nothing
}
\description{
converts vcfR format data to fineRadStructure infile

in part based on vcfR2genepop function
}
\details{
This function converts the vcfR object to a fineRadStructure formatted input file
The function expects multiallelic loci
}
\examples{
vcf2fineRadStructure(vcf = my_vcf, ind_pop = ind_pop, keep_pop = keepers, inc_missing = TRUE, out_file = "fineRadStructure_infile.txt")
vcf2fineRadStructure(my_vcf, ind_pop, keepers, out_file = "fineRadStructure_infile.txt")
vcf2fineRadStructure(my_vcf, ind_pop, keepers)

}
\author{
Tomas Hrbek July 2022
}
