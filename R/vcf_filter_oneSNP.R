#' @title vcf_filter_oneSNP
#' @description subsets vcfR format data keeping only 1 SNP per locus
#' @author Tomas Hrbek February 2022
#'
#' @param vcf -> vcfR object
#' @export
#' @return subsetted vcfR object
#'
#' @details
#' This function subsets the vcfR object keeping only 1 SNP per locus, returning new vcfR object
#' The first SNP independent of quality is taken (may mofify this in future)
#'
#' @examples
#' vcf_filter_oneSNP(vcf = my_vcf)
#' vcf_filter_oneSNP(my_vcf)
#'

vcf_filter_oneSNP <- function(vcf) {
  # read all loci names in vcf
  chrom <- vcfR::getCHROM(vcf)

  # keep only those loci with 1 SNP
  vcf <- vcf[!duplicated(chrom),]

  return(vcf)
}
