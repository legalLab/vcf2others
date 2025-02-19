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
  gt <- vcfR::extract.gt(vcf, convertNA = TRUE)
  # get number of samples in vcf
  n_samples <- ncol(gt)

  # keep only those loci with <= % missing data
  rows_to_keep <- rowSums(is.na(gt)) / n_samples <= p_miss
  vcf <- vcf[rows_to_keep, ]

  # print VCF matrix completeness
  gt <- vcfR::extract.gt(vcf, convertNA = TRUE) %>%
    tibble::as_tibble()
  p_missing <- arrow::as_arrow_table(gt) %>%
    dplyr::summarise(across(everything(), ~ sum(is.na(.)))) %>%
    dplyr::collect() %>%
    unlist() %>%
    {sum(.) / (ncol(gt) * nrow(gt))}

  cat(paste0("Final % missing data in VCF is ", round(p_missing*100, 2), "%\n"))

  return(vcf)
}
