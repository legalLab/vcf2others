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
  gt <- vcfR::extract.gt(vcf, convertNA = T)
  # remove invariant loci
  y <- vector(length = nrow(gt))
  for(i in 1:length(y)) {
    # keep if num unique loci > 1
    y[i] <- length(unique(na.omit(gt[i,]))) > 1
  }
  vcf <- vcf[y,]

  return(vcf)
}
