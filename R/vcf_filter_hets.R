#' @title vcf_filter_hets
#' @description remove loci above het threshold from vcfR format data
#' @description high hets are indicative of paralogs
#' @author Tomas Hrbek October 2022
#'
#' @param vcf -> vcfR object
#' @export
#' @return subsetted vcfR object
#'
#' @details
#' This function removes loci above het threshold from the vcfR object
#'
#' @examples
#' vcf_filter_hets(vcf = my_vcf, hets = .4)
#' vcf_filter_hets(my_vcf)
#'

vcf_filter_hets <- function(vcf, hets = .4) {
  vcf <- vcfR::extract.indels(vcf, return.indels = FALSE)
  vcf <- vcf[vcfR::is.biallelic(vcf), ]
  gt <- vcfR::extract.gt(vcf, return.alleles = FALSE, convertNA = TRUE) %>%
    t() %>%
    tibble::as_tibble()
  rows_to_keep <- arrow::as_arrow_table(gt) %>%
    dplyr::mutate(across(everything(), ~ case_when(
      . == "0/0" | . == "0|0" ~ 0L,
      . == "1/1" | . == "1|1" ~ 1L,
      . == "0/1" | . == "0|1" | . == "1/0" | . == "1|0" ~ 2L,
      TRUE ~ NA_integer_
    ))) %>%
    dplyr::summarise(across(everything(), ~ sum(. == 2, na.rm = TRUE) / sum(!is.na(.)) < hets)) %>%
    dplyr::collect() %>%
    unlist()

  vcf <- vcf[rows_to_keep, ]

  return(vcf)
}
