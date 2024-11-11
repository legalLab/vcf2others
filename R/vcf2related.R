#' @title vcf2related
#' @description converts vcfR format data to Related infile
#' @description in part based on vcfR2genepop function (vcfR package)
#' @author Tomas Hrbek April 2023
#'
#' @param vcf -> vcfR object
#' @param ind_pop -> population assignment of individuals in vcf (factor)
#' @param keep_pop -> population(s) of interest to include in Related infile (factor)
#' @param inc_missing -> include missing data (logical)
#' @param out_file -> name of file to output (Related infile)
#' @export
#' @return nothing
#'
#' @details
#' This function converts the vcfR object to a Related formatted input file
#' The function will remove indels, and multiallelic loci, and optionally loci with missing data
#' 'A', 'C', 'G', 'T' is 01, 02, 03, 04
#'
#' @examples
#' vcf2related(vcf = my_vcf, ind_pop = ind_pop, keep_pop = keepers, inc_missing = TRUE, out_file = "Genepop_infile.gen")
#' vcf2related(my_vcf, ind_pop, keepers, out_file = "Related_infile.txt")
#' vcf2related(my_vcf, ind_pop, keepers)
#'

vcf2related <- function (vcf, ind_pop, keep_pop, inc_missing = TRUE,
                         out_file = "related_infile.txt")
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
    vcf <- vcf[!rowSums((is.na(gt))), ]
  }
  vcf2 <- vcf_extract_pops(vcf, ind_pop, keep_pop)

  gt <- vcfR::extract.gt(vcf2, return.alleles = TRUE, convertNA = TRUE) %>%
    tibble::as_tibble()
  gt_table <- arrow::as_arrow_table(gt) %>%
    dplyr::mutate(across(everything(), ~ dplyr::case_when(
      . == "." ~ "NA\tNA",
      . == "A/A" | . == "A|A" ~ "01\t01",
      . == "A/C" | . == "C/A" | . == "A|C" | . == "C|A" ~ "01\t02",
      . == "A/G" | . == "G/A" | . == "A|G" | . == "G|A" ~ "01\t03",
      . == "A/T" | . == "T/A" | . == "A|T" | . == "T|A" ~ "01\t04",
      . == "C/C" | . == "C|C" ~ "02\t02",
      . == "C/G" | . == "G/C" | . == "C|G" | . == "G|C" ~ "02\t03",
      . == "C/T" | . == "T/C" | . == "C|T" | . == "T|C" ~ "02\t04",
      . == "G/G" | . == "G|G" ~ "03\t03",
      . == "G/T" | . == "T/G" | . == "G|T" | . == "T|G" ~ "03\t04",
      . == "T/T" | . == "T|T" ~ "04\t04",
      TRUE ~ .
    ))) %>%
    dplyr::collect() %>%
    as.data.frame()

  utils::write.table(gt_table, file = out_file, quote = FALSE, sep = "\t", col.names = FALSE, row.names = TRUE)

  invisible(vcf)
}
