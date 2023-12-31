% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/vcf_sub_loci.R
\name{vcf_sub_loci}
\alias{vcf_sub_loci}
\title{vcf_sub_loci}
\usage{
vcf_sub_loci(vcf, n_loci = 1000)
}
\arguments{
\item{vcf}{-> vcfR object}

\item{n_loci}{-> number of loci to subset (numeric)}
}
\value{
subsetted vcfR object
}
\description{
randomly subsets vcfR format data to # of loci
}
\details{
This function subsets the vcfR object to specific # of loci, returning new vcfR object
}
\examples{
vcf_sub_loci(vcf = my_vcf, n_loci = n_loci)
vcf_sub_loci(my_vcf, n_loci)

}
\author{
Tomas Hrbek January 2023
}
