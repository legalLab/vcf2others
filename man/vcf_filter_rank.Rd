% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/vcf_filter_rank.R
\name{vcf_filter_rank}
\alias{vcf_filter_rank}
\title{vcf_filter_rank}
\usage{
vcf_filter_rank(vcf, rank = 0.4)
}
\arguments{
\item{vcf}{-> vcfR object}
}
\value{
subsetted vcfR object
}
\description{
remove loci below rank threshold from vcfR format data

rank is calculated as sqrt(chi-sqr/n) of allele read counts

used for paralog detection - very low rank values (<0.4)

rank is calculated in DiscoSNP, registered as Pk in INFO
}
\details{
This function removes loci below rank threshold from the vcfR object
}
\examples{
vcf_filter_quality(vcf = my_vcf, rank = .4)
vcf_filter_quality(my_vcf)

}
\author{
Tomas Hrbek February 2022
}
