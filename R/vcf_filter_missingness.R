#' @title vcf_filter_missingness
#' @description subsets vcfR format data by % missing data
#' @author Tomas Hrbek February 2022
#'
#' @param vcf -> vcfR object
#' @param p_miss -> max missing data per locus as decimal (numeric)
#' @export
#' @return subsetted vcfR object
#'
#' @details
#' This function subsets the vcfR object by % missing data, returning new vcfR object
#'
#' @examples
#' vcf_filter_missingness(vcf = my_vcf, p_miss = p_miss)
#' vcf_filter_missingness(my_vcf, p_miss)
#'

vcf_filter_missingness <- function(vcf, p_miss) {
  gt <- vcfR::extract.gt(vcf, convertNA = T)
  # get number of samples in vcf
  n_samples <- ncol(gt)

  # keep only those loci with < % missing data
  vcf <- vcf[rowSums(is.na(gt)) < floor(n_samples*p_miss),]

  # print VCF matrix completeness
  vcf1 <- vcf_filter_oneSNP(vcf)
  gt <- vcfR::extract.gt(vcf1, convertNA = T)
  p_missing <- sum(is.na(gt)) / length(gt)
  print(paste("final % missing data in VCF is", round(p_missing*100, 2), "%", sep = " "))

  return(vcf)
}
