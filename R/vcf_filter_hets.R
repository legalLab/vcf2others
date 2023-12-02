#' @title vcf_filter_hets
#' @description remove loci above het threshold from vcfR format data
#' @description high hets are indicative of paralogs
#' @author Tomas Hrbek October 2022
#'
#' @param vcf -> vcfR object
#' @export
#' @return subsetted vcfR object
#'
#' @details
#' This function removes loci above het threshold from the vcfR object
#'
#' @examples
#' vcf_filter_hets(vcf = my_vcf, hets = .4)
#' vcf_filter_hets(my_vcf)
#'

vcf_filter_hets <- function(vcf, hets = .4) {
  vcf <- vcfR::extract.indels(vcf, return.indels = F)
  vcf <- vcf[vcfR::is.biallelic(vcf), ]
  gt <- vcfR::extract.gt(vcf, return.alleles = F, convertNA = T)
  gt[is.na(gt)] <- "NA"
  gt[gt == "0/0" | gt == "0|0"] <- "0"
  gt[gt == "1/1" | gt == "1|1"] <- "1"
  gt[gt == "0/1" | gt == "0|1" | gt == "1/0" | gt == "1|0"] <- "2"

  homo <- apply(gt, MARGIN = 1, function(x){sum(x == 0 | x == 1, na.rm = TRUE)})
  hetero <- apply(gt, MARGIN = 1, function(x){sum(x == 2, na.rm = TRUE)})

  vcf <- vcf[hetero/(homo+hetero) <= hets,]

  return(vcf)
}
