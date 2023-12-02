#' @title vcf_filter_multiSNP
#' @description subsets vcfR format data keeping only loci with 2+ SNPs
#' @author Tomas Hrbek February 2022
#'
#' @param vcf -> vcfR object
#' @export
#' @return subsetted vcfR object
#'
#' @importFrom magrittr %>%
#'
#' @details
#' This function subsets the vcfR object keeping only loci with between min and max # of SNPs
#' per locus, returning new vcfR object
#' default min = 2 and max = 5 SNPs per locus
#' Recommended as input for fineRADstructure analyses
#'
#' @examples
#' vcf_filter_multiSNP(vcf = my_vcf, minS = 2, maxS = 5)
#' vcf_filter_multiSNP(vcf = my_vcf)
#' vcf_filter_multiSNP(my_vcf)
#'

vcf_filter_multiSNP <- function(vcf, minS = 2, maxS = 5) {
  # read all loci names in vcf
  chrom <- vcfR::getCHROM(vcf)
  id <-  vcfR::getID(vcf)
  chrom_id <- cbind(chrom, id) %>%
    dplyr::as_tibble()

  keeper <- chrom_id %>%
    dplyr::count(chrom) %>%
    dplyr::filter(n >= minS & n <= maxS) %>%
    dplyr::select(-n) %>%
    as.matrix() %>%
    as.character()

  # keep only those loci with between minS and maxS SNPs
  vcf <- vcf[chrom %in% keeper, ]

  return(vcf)
}
