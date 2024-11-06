#' @title vcf_sub_loci
#' @description randomly subsets vcfR format data to # of loci
#' @author Tomas Hrbek January 2023
#'
#' @param vcf -> vcfR object
#' @param n_loci -> number of loci to subset (numeric)
#' @export
#' @return subsetted vcfR object
#'
#' @details
#' This function subsets the vcfR object to specific # of loci, returning new vcfR object
#'
#' @examples
#' vcf_sub_loci(vcf = my_vcf, n_loci = n_loci)
#' vcf_sub_loci(my_vcf, n_loci)
#'

vcf_sub_loci <- function(vcf, n_loci = 1000) {
  gt <- vcfR::extract.gt(vcf, convertNA = T)
  # get number of snps in vcf
  n_snps <- nrow(gt)
  
  # subsample and keep the same order of snps
  x <- sort(sample(1:n_snps, n_loci))

  vcf <- vcf[x,]

  return(vcf)
}
