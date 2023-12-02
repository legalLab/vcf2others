#' @title vcf_filter_rank
#' @description remove loci below rank threshold from vcfR format data
#' @description rank is calculated as sqrt(chi-sqr/n) of allele read counts
#' @description used for paralog detection - very low rank values (<0.4)
#' @description rank is calculated in DiscoSNP, registered as Pk in INFO
#' @author Tomas Hrbek February 2022
#'
#' @param vcf -> vcfR object
#' @export
#' @return subsetted vcfR object
#'
#' @importFrom magrittr %>%
#'
#' @details
#' This function removes loci below rank threshold from the vcfR object
#'
#' @examples
#' vcf_filter_quality(vcf = my_vcf, rank = .4)
#' vcf_filter_quality(my_vcf)
#'

vcf_filter_rank <- function(vcf, rank = .4) {
  snp_rank <- stringr::str_extract(vcf@fix[,8], "Rk=[0-9|.]*") %>%
    stringr::str_extract("[^Rk=]+")
  # keep only those loci with minimum rank
  vcf <- vcf[snp_rank >= rank,]

  return(vcf)
}
