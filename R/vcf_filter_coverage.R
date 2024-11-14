#' @title vcf_filter_coverage
#' @description remove genotypes below coverage threshold from vcfR format data
#' @author Tomas Hrbek November 2022
#'
#' @param vcf -> vcfR object
#' @export
#' @return subsetted vcfR object
#'
#' @details
#' This function removes genotypes below coverage threshold from the vcfR object
#'
#' @examples
#' vcf_filter_coverage(vcf = my_vcf, cover = 10)
#' vcf_filter_coverage(my_vcf)
#'

vcf_filter_coverage <- function(vcf, cover = 10, f_invar = TRUE) {
  dp <- vcfR::extract.gt(vcf, element = "DP", as.numeric = TRUE) %>%
    tibble::as_tibble()

  # create a matrix indicating positions below the coverage threshold
  mask <- arrow::as_arrow_table(dp) %>%
    dplyr::mutate(across(everything(), ~ . < cover)) %>%
    dplyr::collect() %>%
    as.matrix()

  # replace loci below threshold in genotype data
  gt_matrix <- as.matrix(vcf@gt[,-1])
  gt_matrix[mask] <- "./.:0:.,.,.:0,0:0,0"
  vcf@gt[,-1] <- gt_matrix

  # remove missing data
  if (f_invar == TRUE) vcf <- vcf_filter_invariant(vcf)

  return(vcf)
}
