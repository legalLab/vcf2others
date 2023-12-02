#' @title vcf_merge
#' @description merge two vcf; merge second vcf into first one
#' @author Tomas Hrbek October 2023
#'
#' @param vcf -> vcfR object
#' @param vcf1 -> vcfR object
#' @export
#' @return augmented vcfR object
#'
#' @importFrom magrittr %>%
#'
#' @details
#' This function adds all individuals from a second vcfR object into a first vcfR object, returning new vcfR object
#'
#' @examples
#' vcf_merge(vcf = my_vcf, vcf1 = other_vcf)
#' vcf_merge(my_vcf, other_vcf)
#'

vcf_merge <- function(vcf, vcf1) {
  vcf_gt <- vcf@gt %>%
    as.data.frame() %>%
    dplyr::mutate(id = vcfR::getID(vcf))
  vcf1_gt <- vcf1@gt %>%
    as.data.frame() %>%
    dplyr::mutate(id = vcfR::getID(vcf1))

  vcf_merged_gt <-dplyr::left_join(vcf_gt, vcf1_gt) %>%
    dplyr::select(-id) %>%
    as.matrix()

  vcf@gt <- vcf_merged_gt

  return(vcf)
}
