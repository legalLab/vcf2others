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
  chrom_id_table <- arrow::Table$create(
    chrom = vcfR::getCHROM(vcf),
    id = vcfR::getID(vcf)
  )

  chrom = vcfR::getCHROM(vcf)
  keeper <- chrom_id_table %>%
    dplyr::group_by(chrom) %>%
    dplyr::summarise(n = dplyr::n()) %>%
    dplyr::filter(n >= minS & n <= maxS) %>%
    dplyr::select(chrom) %>%
    dplyr::collect() %>%
    dplyr::pull(chrom)

  # keep only those loci with between minS and maxS SNPs
  vcf <- vcf[chrom %in% keeper, ]

  return(vcf)
}
