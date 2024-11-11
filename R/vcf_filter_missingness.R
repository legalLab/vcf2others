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
  gt <- vcfR::extract.gt(vcf, convertNA = TRUE) %>%
    t() %>%
    tibble::as_tibble()
  # get number of samples in vcf
  n_samples <- ncol(gt)

  rows_to_keep <- arrow::as_arrow_table(gt) %>%
    dplyr::summarise(across(everything(), ~ sum(is.na(.)) / n_samples < p_miss)) %>%
    dplyr::collect() %>%
    unlist()
  vcf1 <- vcf[rows_to_keep, ]

  # print VCF matrix completeness
  vcf1 <- vcf_filter_oneSNP(vcf)
  gt <- vcfR::extract.gt(vcf1, convertNA = TRUE) %>%
    tibble::as_tibble()
  p_missing <- arrow::as_arrow_table(gt) %>%
    dplyr::summarise(across(everything(), ~ sum(is.na(.)))) %>%
    dplyr::collect() %>%
    unlist() %>%
    {sum(.) / (ncol(gt) * nrow(gt))}

  print(paste("final % missing data in VCF is", round(p_missing*100, 2), "%", sep = " "))

  return(vcf)
}
