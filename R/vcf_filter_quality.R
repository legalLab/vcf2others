#' @title vcf_filter_quality
#' @description remove loci below quality threshold from vcfR format data
#' @author Tomas Hrbek February 2022
#'
#' @param vcf -> vcfR object
#' @export
#' @return subsetted vcfR object
#'
#' @details
#' This function removes loci below quality threshold from the vcfR object
#'
#' @examples
#' vcf_filter_quality(vcf = my_vcf, qual = 20)
#' vcf_filter_quality(my_vcf)
#'

vcf_filter_quality <- function(vcf, qual = 20) {
  if(any(is.na(vcfR::getQUAL(vcf)))){
    print("No quality information in VCF; keeping VCF as is")
  } else {
    # keep only those loci with minimum quality
    vcf <- vcf[vcfR::getQUAL(vcf) >= qual, ]
  }

  return(vcf)
}
