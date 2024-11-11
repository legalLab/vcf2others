#' @title vcf_filter_invariant
#' @description remove invariant loci from vcfR format data
#' @author Tomas Hrbek February 2022
#'
#' @param vcf -> vcfR object
#' @export
#' @return subsetted vcfR object
#'
#' @details
#' This function removes invariant loci from the vcfR object
#' This might be desirable after subsetting the vcf by individuals
#'
#' @examples
#' vcf_filter_invariant(vcf = my_vcf)
#' vcf_filter_invariant(my_vcf)
#'

vcf_filter_invariant <- function(vcf) {
  gt <- vcfR::extract.gt(vcf, convertNA = TRUE) %>%
    t() %>%
    tibble::as_tibble()

  # remove invariant loci
  rows_to_keep <- arrow::as_arrow_table(gt) %>%
    dplyr::summarise(across(everything(), ~ length(unique(na.omit(.))) > 1)) %>%
    dplyr::collect() %>%
    unlist()
  vcf <- vcf[rows_to_keep, ]

  return(vcf)
}
