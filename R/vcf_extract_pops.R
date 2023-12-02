#' @title vcf_extract_pops
#' @description extract vcfR format data by group
#' @author Tomas Hrbek December 2020
#'
#' @param vcf -> vcfR object
#' @param ind_pop -> group assignment of individuals in vcf (factor)
#' @param keep_pop -> group(s) of interest to include/exclude in vcf infile (factor)
#' @export
#' @return extracted vcfR object
#'
#' @details
#' This function extracts the vcfR object by group(s) (individuals
#' assigned to one or more group), returning new vcfR object
#'
#' @examples
#' vcf_extract_pops(vcf = my_vcf, ind_pop = ind_pop, keep_pop = keepers, whitelist = TRUE)
#' vcf_extract_pops(vcf = my_vcf, ind_pop = ind_pop, keep_pop = non_keepers, whitelist = FALSE)
#' vcf_extract_pops(my_vcf, ind_pop, keepers)
#'

vcf_extract_pops <- function(vcf, ind_pop, keep_pop, whitelist = TRUE, f_invar = TRUE) {
  ids <- which(ind_pop %in% keep_pop)
  ids <- ids+1
  if (whitelist == TRUE){
    vcf <- vcf[, c(1,ids)] %>%
      {if (f_invar == TRUE) vcf_filter_invariant(.) else .}
  }
  else {
    vcf <- vcf[, -c(ids)] %>%
      {if (f_invar == TRUE) vcf_filter_invariant(.) else .}
  }

  return(vcf)
}
