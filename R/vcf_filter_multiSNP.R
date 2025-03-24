#' @title vcf_filter_multiSNP
#' @description subsets vcfR format data keeping only loci with 2+ SNPs, selecting upper limit of SNPs per locus
#' @author Tomas Hrbek February 2022
#'
#' @param vcf -> vcfR object
#' @param block_size -> size of linked SNP blocks (integer)
#' @param minS -> minimum linked block size (integer)
#' @param maxS -> maximum number of selected linked SNPs per block (integer)
#' @export
#' @return subsetted vcfR object
#'
#' @importFrom magrittr %>%
#'
#' @details
#' This function subsets the vcfR object keeping only loci with between min and max # of SNPs
#' per locus, returning new vcfR object
#' default min = 2 and max = 5 SNPs per locus
#' (recommended as input for fineRADstructure analyses)
#'
#' @examples
#' vcf_filter_multiSNP(vcf = my_vcf, block_size = 10000, minS = 2, maxS = 5)
#' vcf_filter_multiSNP(vcf = my_vcf)
#' vcf_filter_multiSNP(my_vcf)
#'

vcf_filter_multiSNP <- function(vcf, block_size = 10000, minS = 2, maxS = 5) {
  chrom_pos_table <- arrow::Table$create(
    chrom = vcfR::getCHROM(vcf),
    pos = vcfR::getPOS(vcf),
    index = 1:length(vcfR::getCHROM(vcf))
  )

  # identify blocks with at least minS SNPs
  blocks_with_enough_snps <- chrom_pos_table %>%
    dplyr::group_by(chrom) %>%
    dplyr::mutate(first_pos = min(pos)) %>%
    dplyr::mutate(block = ((pos - first_pos) %/% block_size) + 1) %>%
    dplyr::ungroup() %>%
    # Count SNPs in each block
    dplyr::group_by(chrom, block) %>%
    dplyr::summarize(snps_in_block = n(), .groups = 'drop') %>%
    dplyr::filter(snps_in_block >= minS) %>%
    dplyr::select(chrom, block) %>%
    dplyr::collect()

  # collect all data and do the final filtering in R
  selected_indices <- chrom_pos_table %>%
    dplyr::group_by(chrom) %>%
    dplyr::mutate(first_pos = min(pos)) %>%
    dplyr::mutate(block = ((pos - first_pos) %/% block_size) + 1) %>%
    dplyr::ungroup() %>%
    # join with the blocks that have enough SNPs
    dplyr::inner_join(blocks_with_enough_snps, by = c('chrom', 'block')) %>%
    # group by chrom and block
    dplyr::group_by(chrom, block) %>%
    dplyr::arrange(chrom, block, pos) %>%
    # create a rank within each group
    dplyr::mutate(rank = dplyr::row_number()) %>% # row_number() is not arrow compatible and is pulled into dplyr
    # keep SNPs with rank <= maxS
    dplyr::filter(rank <= maxS) %>%
    dplyr::ungroup() %>%
    dplyr::collect() %>%
    dplyr::pull(index)

  # Filter VCF to keep only selected SNPs
  vcf <- vcf[selected_indices, ]

  return(vcf)
}
