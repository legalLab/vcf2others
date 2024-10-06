#' @title vcf_extract_indivs
#' @description extract vcfR format data by individuals
#' @author Tomas Hrbek December 2020
#'
#' @param vcf -> vcfR object
#' @param indiv -> individuals to retain/drop in vcf (factor)
#' @export
#' @return extracted vcfR object
#'
#' @details
#' This function extracts the vcfR object by individuals, returning new vcfR object
#'
#' @examples
#' vcf_extract_indivs(vcf = my_vcf, indiv = indivs_to_keep, whitelist = TRUE)
#' vcf_extract_indivs(vcf = my_vcf, indiv = indivs_to_drop, whitelist = FALSE)
#' vcf_extract_indivs(my_vcf, indivs_to_keep)
#'

vcf_extract_indivs <- function(vcf, indiv, whitelist = TRUE, f_invar = TRUE) {
  # read all sample names in vcf
  vcf_names <- colnames(vcf@gt)[-1]

  # allow empty indiv list - keep all individuals
  if (length(indiv) == 0){
    indiv <- vcf_names
    whitelist <- TRUE
  }

  ids <- which(vcf_names %in% indiv)
  ids <- ids+1
  if (whitelist == TRUE){
    vcf <- vcf[, c(1,ids)] %>%
      {if (f_invar == TRUE & length(ids) > 1) vcf_filter_invariant(.) else .}
  }
  else {
    vcf <- vcf[, -c(ids)] %>%
      {if (f_invar == TRUE & length(ids) > 1) vcf_filter_invariant(.) else .}
  }

  return(vcf)
}
