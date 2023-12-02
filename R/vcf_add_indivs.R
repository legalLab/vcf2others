#' @title vcf_add_indivs
#' @description add individuals to vcf
#' @author Tomas Hrbek October 2023
#'
#' @param vcf -> vcfR object
#' @param vcf1 -> vcfR object
#' @param indiv -> individuals to add into vcf from vcf1 (factor)
#' @export
#' @return augmented vcfR object
#'
#' @importFrom magrittr %>%
#'
#' @details
#' This function adds individuals from one vcfR object to another vcfR object, returning new vcfR object
#'
#' @examples
#' vcf_add_indivs(vcf = my_vcf, vcf1 = other_vcf, indiv = indivs_to_add, whitelist = TRUE)
#' vcf_add_indivs(vcf = my_vcf, vcf1 = other_vcf, indiv = indivs_to_not_add, whitelist = FALSE)
#' vcf_add_indivs(my_vcf, other_vcf, indivs_to_add)
#'

vcf_add_indivs <- function(vcf, vcf1, indiv, whitelist = TRUE) {
  # read all sample names in vcf
  vcf1_names <- colnames(vcf1@gt)[-1]

  # allow empty indiv list - keep all individuals
  if (length(indiv) == 0){
    indiv <- vcf1_names
    whitelist <- TRUE
  }

  vcf_gt <- vcf@gt %>%
    as.data.frame() %>%
    dplyr::mutate(id = vcfR::getID(vcf))
  vcf1_gt <- vcf1@gt %>%
    as.data.frame() %>%
    dplyr::mutate(id = vcfR::getID(vcf1))

  ids <- which(vcf1_names %in% indiv)
  ids <- ids+1
  if (whitelist == TRUE){
    vcf_merged_gt <-dplyr::left_join(vcf_gt, vcf1_gt[, c(1,ids,ncol(vcf1_gt))]) %>%
      dplyr::select(-id) %>%
      as.matrix()
  }
  else {
    vcf_merged_gt <-dplyr::left_join(vcf_gt, vcf1_gt[, -c(ids)]) %>%
      dplyr::select(-id) %>%
      as.matrix()
  }
  vcf@gt <- vcf_merged_gt

  return(vcf)
}
