#' @title vcf_filter_maf
#' @description remove loci below MAF threshold from vcfR format data
#' @author Tomas Hrbek October 2022
#'
#' @param vcf -> vcfR object
#' @export
#' @return subsetted vcfR object
#'
#' @details
#' This function removes loci below MAF threshold from the vcfR object
#'
#' @examples
#' vcf_filter_maf(vcf = my_vcf, maf = .05)
#' vcf_filter_maf(my_vcf)
#'

vcf_filter_maf <- function(vcf, maf = .05) {
  # remove NAs from a vector
  rm_na <- function(x) {
    y <- x[!is.na(x)]
    return(y)
  }

  vcf <- vcfR::extract.indels(vcf, return.indels = F)
  vcf <- vcf[vcfR::is.biallelic(vcf), ]
  gt <- vcfR::extract.gt(vcf, return.alleles = F, convertNA = T)
  gt[is.na(gt)] <- "NA"
  gt[gt == "0/0" | gt == "0|0"] <- "0"
  gt[gt == "1/1" | gt == "1|1"] <- "1"
  gt[gt == "0/1" | gt == "0|1" | gt == "1/0" | gt == "1|0"] <- "2"

  alleles_all <- apply(gt, MARGIN = 1, function(x){2*length(rm_na(x))})
  allele_rare <- apply(gt, MARGIN = 1, function(x){sum(x == 2, na.rm = TRUE) + 2*sum(x == 1, na.rm = TRUE)})

  vcf <- vcf[allele_rare/alleles_all >= maf,]

  return(vcf)
}
