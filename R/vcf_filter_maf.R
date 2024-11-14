#' @title vcf_filter_maf
#' @description remove loci below MAF threshold from vcfR format data
#' @author Tomas Hrbek October 2022
#'
#' @param vcf -> vcfR object
#' @export
#' @return subsetted vcfR object
#'
#' @details
#' This function removes loci below MAF threshold from the vcfR object
#'
#' @examples
#' vcf_filter_maf(vcf = my_vcf, maf = .05)
#' vcf_filter_maf(my_vcf)
#'

vcf_filter_maf <- function(vcf, maf = .05) {
  vcf <- vcfR::extract.indels(vcf, return.indels = FALSE)
  vcf <- vcf[vcfR::is.biallelic(vcf), ]
  gt <- vcfR::extract.gt(vcf, return.alleles = FALSE, convertNA = TRUE) %>%
    tibble::as_tibble()

  gt_table <- arrow::as_arrow_table(gt) %>%
    dplyr::mutate(across(everything(), ~ case_when(
      . == "0/0" | . == "0|0" ~ 0L,
      . == "1/1" | . == "1|1" ~ 2L,
      . == "0/1" | . == "0|1" | . == "1/0" | . == "1|0" ~ 1L,
      TRUE ~ NA_integer_
    ))) %>%
    dplyr::collect()

  rows_to_keep <- apply(gt_table, MARGIN = 1,
                  function(x){sum(as.numeric(x), na.rm = TRUE) /
                      (2*sum(!is.na(x), na.rm = TRUE))}) >= maf
  vcf <- vcf[rows_to_keep, ]

  return(vcf)
}
