#' @title vcf_filter_missing_indivs
#' @description remove from vcfR format data indivs with >% missing data
#' @author Tomas Hrbek October 2022
#'
#' @param vcf -> vcfR object
#' @param p_miss -> max missing data per locus as decimal (numeric)
#' @export
#' @return subsetted vcfR object
#'
#' @details
#' This function subsets the vcfR object by % missing data within
#' and individual, returning new vcfR object
#'
#' @examples
#' vcf_filter_missing_indivs(vcf = my_vcf, p_miss = p_miss)
#' vcf_filter_missing_indivs(my_vcf, p_miss)
#'

vcf_filter_missing_indivs <- function(vcf, p_miss, f_invar = TRUE) {
  gt <- vcfR::extract.gt(vcf, convertNA = T)
  # get number of snps in vcf
  n_snps <- nrow(gt)

  # report which samples above threshold
  samples_removed <- colnames(gt)[colSums(is.na(gt)) > floor(n_snps*p_miss)]
  cat(paste("removed samples are:", samples_removed, "\n", sep = " "))

  # keep only those loci with < % missing data
  vcf <- vcf[, c(TRUE, colSums(is.na(gt)) < floor(n_snps*p_miss))] %>%
    {if (f_invar == TRUE) vcf_filter_invariant(.) else .}

  # print VCF matrix completeness
  vcf1 <- vcf_filter_oneSNP(vcf)
  gt <- vcfR::extract.gt(vcf1, convertNA = T)
  p_missing <- sum(is.na(gt)) / length(gt)
  print(paste("final % missing data in VCF is", round(p_missing*100, 2), "%", sep = " "))

  return(vcf)
}
