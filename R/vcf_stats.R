#' @title vcf_stats
#' @description calculates basic stats of each samples from vcfR format data
#' @description high heterozygosity is indicative of potential contamination
#' @description average read depth per individual
#' @description missing data per individual
#'
#' @author Tomas Hrbek August 2024
#'
#' @param vcf -> vcfR object
#' @export
#' @return table of statistics
#'
#' @details
#' This function calculates average read depth, heterozygosity
#' number of heterozygotes, number of reference and alternative homozygotes,
#' missing data and total number SNPs of each sample in an vcfR object
#'
#' @examples
#' vcf_stats(vcf = my_vcf)
#' vcf_stats(my_vcf)
#'

vcf_stats <- function(vcf, res_path, project) {
  vcf <- vcfR::extract.indels(vcf, return.indels = FALSE)
  vcf <- vcf[vcfR::is.biallelic(vcf), ]
  dp_table <- vcfR::extract.gt(vcf, element = "DP", as.numeric = TRUE) %>%
    tibble::as_tibble() %>%
    arrow::as_arrow_table()
  gt <- vcfR::extract.gt(vcf, return.alleles = FALSE, convertNA = TRUE) %>%
    tibble::as_tibble()
  gt_table <- arrow::as_arrow_table(gt) %>%
    dplyr::mutate(across(everything(), ~ dplyr::case_when(
      . == "0/0" | . == "0|0" ~ 0L,
      . == "1/1" | . == "1|1" ~ 2L,
      . == "0/1" | . == "0|1" | . == "1/0" | . == "1|0" ~ 1L,
      TRUE ~ NA_integer_
    )))
  samples <- colnames(gt)

  homo_ref <- gt_table %>%
    dplyr::summarise(across(everything(), ~ sum(. == 0, na.rm = TRUE))) %>%
    dplyr::collect() %>%
    unlist()
  homo_alt <- gt_table %>%
    dplyr::summarise(across(everything(), ~ sum(. == 2, na.rm = TRUE))) %>%
    dplyr::collect() %>%
    unlist()
  hetero <- gt_table %>%
    dplyr::summarise(across(everything(), ~ sum(. == 1, na.rm = TRUE))) %>%
    dplyr::collect() %>%
    unlist()
  miss <- gt_table %>%
    dplyr::summarise(across(everything(), ~ sum(is.na(.)))) %>%
    dplyr::collect() %>%
    unlist()
  total <- homo_ref + homo_alt + hetero
  loc_n <- nrow(gt)
  depth <- dp_table %>%
    dplyr::summarise(across(everything(), ~ mean(., na.rm = TRUE))) %>%
    dplyr::collect() %>%
    unlist()

  stats <- tibble::tibble(
    id = samples,
    read_depth = depth,
    heterozygosity = hetero / total,
    heterozygotes = hetero,
    homozygotes = homo_ref + homo_alt,
    homozygotes_ref = homo_ref,
    homozygotes_alt = homo_alt,
    missing_p = miss / loc_n,
    missing = miss,
    non_missing = total,
    total = loc_n
  )

  write.table(stats, file = paste0(res_path, project, "_stats.csv"), row.names = FALSE, quote = FALSE, sep = ",")

  # invisible return of the first argument so function can be used in a pipe
  invisible(vcf)
}
