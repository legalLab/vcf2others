#' @title vcf_filter_oneSNP
#' @description subsets vcfR format data keeping only 1 SNP per locus
#' @author Tomas Hrbek February 2022
#'
#' @param vcf -> vcfR object
#' @param block_size -> size of linked SNP blocks (integer)
#' @export
#' @return subsetted vcfR object
#'
#' @details
#' This function subsets the vcfR object keeping only 1 SNP per locus, returning new vcfR object
#' Locus is defined as a different chromosome or a block within a chromosome
#' The first SNP independent of quality is taken (may mofify this in the future)
#'
#' @examples
#' vcf_filter_oneSNP(vcf = my_vcf)
#' vcf_filter_oneSNP(my_vcf, block_size = 1000)
#' vcf_filter_oneSNP(my_vcf)
#'

vcf_filter_oneSNP <- function(vcf, block_size = 10000) {
  # read all loci names and positions of their SNPs in VCF
  chrom_pos_table <- arrow::Table$create(
    chrom = vcfR::getCHROM(vcf),
    pos = vcfR::getPOS(vcf),
    index = 1:length(vcfR::getCHROM(vcf))
  )

  # compute blocks within each chromosome
  selected_indices <- chrom_pos_table %>%
    dplyr::group_by(chrom) %>%
    dplyr::mutate(first_pos = min(pos)) %>%
    dplyr::mutate(block = ((pos - first_pos) %/% block_size) + 1) %>%
    dplyr::ungroup() %>%
    dplyr::arrange(chrom, block, pos) %>%
    # select the first SNP from each chromosome/block combination
    dplyr::distinct(chrom, block, .keep_all = TRUE) %>%
    dplyr::collect() %>%
    dplyr::pull(index)

  # filter VCF to keep only selected SNPs
  vcf <- vcf[selected_indices, ]

  return(vcf)
}
