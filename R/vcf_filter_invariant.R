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
  gt <- vcfR::extract.gt(vcf, convertNA = TRUE)

  # remove invariant loci
  unique_counts <- apply(gt, 1, function(x) length(unique(na.omit(x))))
  rows_to_keep <- unique_counts > 1
  vcf <- vcf[rows_to_keep, ]

  return(vcf)
}

