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
  gt <- vcfR::extract.gt(vcf, convertNA = TRUE) %>%
    tibble::as_tibble()
  # get number of snps in vcf
  n_snps <- nrow(gt)

  cols_to_keep <- arrow::as_arrow_table(gt) %>%
    dplyr::summarise(across(everything(), ~ mean(is.na(.))  < p_miss)) %>%
    dplyr::collect() %>%
    unlist()
  vcf <- vcf[, c(TRUE, cols_to_keep)] %>%
    {if (f_invar == TRUE) vcf_filter_invariant(.) else .}

  removed_samples <- colnames(gt)[!cols_to_keep]
  cat(paste0("Removed samples are: ", paste(removed_samples, collapse = "\n"), "\n"))

  # print VCF matrix completeness
  vcf1 <- vcf_filter_oneSNP(vcf)
  gt <- vcfR::extract.gt(vcf1, convertNA = TRUE) %>%
    tibble::as_tibble()
  p_missing <- arrow::as_arrow_table(gt) %>%
    dplyr::summarise(across(everything(), ~ sum(is.na(.)))) %>%
    dplyr::collect() %>%
    unlist() %>%
    {sum(.) / (ncol(gt) * nrow(gt))}

  print(paste0("final % missing data in VCF is ", round(p_missing*100, 2), "%"))

  return(vcf)
}
