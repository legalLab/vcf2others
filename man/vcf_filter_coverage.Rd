% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/vcf_filter_coverage.R
\name{vcf_filter_coverage}
\alias{vcf_filter_coverage}
\title{vcf_filter_coverage}
\usage{
vcf_filter_coverage(vcf, cover = 10, f_invar = TRUE)
}
\arguments{
\item{vcf}{-> vcfR object}
}
\value{
subsetted vcfR object
}
\description{
remove genotypes below coverage threshold from vcfR format data
}
\details{
This function removes genotypes below coverage threshold from the vcfR object
}
\examples{
vcf_filter_coverage(vcf = my_vcf, cover = 10)
vcf_filter_coverage(my_vcf)

}
\author{
Tomas Hrbek November 2022
}
