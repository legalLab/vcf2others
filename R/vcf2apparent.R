#' @title vcf2apparent
#' @description converts vcfR format data to Apparent infile
#' @author Tomas Hrbek February 2023
#'
#' @param vcf -> vcfR object
#' @param ind_pop -> population assignment of individuals in vcf (factor)
#' @param keep_pop -> population(s) of interest to include in Apparent infile (factor)
#' @param key -> relationship type (All, Pa, Mo, Fa, Off); default All (character)
#' @param inc_missing -> include missing data (logical)
#' @param out_file -> name of file to output (Apparent infile)
#' @export
#' @return nothing
#'
#' @details
#' This function converts the vcfR object to a Apparent formatted input file
#' The function will remove indels, and multiallelic loci, and optionally loci with missing data
#'
#' @examples
#' vcf2apparent(vcf = my_vcf, ind_pop = ind_pop, keep_pop = keepers, key = key, inc_missing = TRUE, out_file = "Genepop_infile.gen")
#' vcf2apparent(my_vcf, ind_pop, keepers, out_file = "Apparent_infile.txt")
#' vcf2apparent(my_vcf, ind_pop, keepers)
#'

vcf2apparent <- function (vcf, ind_pop, keep_pop, key = "All", inc_missing = TRUE,
                          out_file = "apparent_infile.txt")
{
  if (class(vcf) != "vcfR") {
    stop(paste("Expecting an object of class vcfR, received a",
               class(vcf), "instead"))
  }
  if (class(ind_pop) != "factor" | class(keep_pop) != "factor") {
    stop(paste("Expecting population vector, received a",
               class(ind_pop), "and", class(keep_pop), "instead"))
  }
  vcf <- vcfR::extract.indels(vcf, return.indels = FALSE)
  vcf <- vcf[vcfR::is.biallelic(vcf), ]
  if (inc_missing == FALSE) {
    gt <- vcfR::extract.gt(vcf, convertNA = TRUE)
    vcf <- vcf[!rowSums(is.na(gt)), ]
  }
  vcf2 <- vcf_extract_pops(vcf, ind_pop, keep_pop)
  gt <- vcfR::extract.gt(vcf2, return.alleles = TRUE, convertNA = TRUE) %>%
    tibble::as_tibble()
  gt_table <- arrow::as_arrow_table(gt) %>%
    dplyr::mutate(across(everything(), ~ dplyr::case_when(
      . == "." ~ "-/-",
      . == "A/A" | . == "A|A" ~ "A/A",
      . == "A/C" | . == "C/A" | . == "A|C" | . == "C|A" ~ "A/C",
      . == "A/G" | . == "G/A" | . == "A|G" | . == "G|A" ~ "A/G",
      . == "A/T" | . == "T/A" | . == "A|T" | . == "T|A" ~ "A/T",
      . == "C/C" | . == "C|C" ~ "C/C",
      . == "C/G" | . == "G/C" | . == "C|G" | . == "G|C" ~ "C/G",
      . == "C/T" | . == "T/C" | . == "C|T" | . == "T|C" ~ "C/T",
      . == "G/G" | . == "G|G" ~ "G/G",
      . == "G/T" | . == "T/G" | . == "G|T" | . == "T|G" ~ "G/T",
      . == "T/T" | . == "T|T" ~ "T/T",
      TRUE ~ .
    ))) %>%
    dplyr::collect() %>%
    t() %>%
    as.data.frame() %>%
    dplyr::mutate(key = key) %>%
    dplyr::relocate(key, .before = everything())

  #suppressWarnings(utils::write.table(gt_table, file = out_file, quote = FALSE, sep = "\t", col.names = TRUE, append = TRUE))

  utils::write.table(gt_table, file = out_file, quote = FALSE, sep = "\t", col.names = FALSE, append = FALSE)

  invisible(vcf)
}
