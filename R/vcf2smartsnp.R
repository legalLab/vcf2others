#' @title vcf2smartsnp
#' @description converts vcfR format data to smartsnp tabular format infile
#' @description in part based on vcfR2genlight function (vcfR package)
#' @author Tomas Hrbek October 2023
#'
#' @param vcf -> vcfR object
#' @param ind_pop -> population assignment of individuals in vcf (factor)
#' @param keep_pop -> population(s) of interest to include in genotype table infile (factor)
#' @param inc_missing -> include missing data (logical)
#' @param out_file -> name of file to output (genotype table infile)
#' @export
#' @return nothing
#'
#' @details
#' This function converts the vcfR object to a genotype table formatted input file
#' For use in smartsnp, but keeps sample names as column names
#' The function will remove indels, and multiallelic loci, and optionally loci with missing data
#'
#' @examples
#' vcf2smartsnp(vcf = my_vcf, ind_pop = ind_pop, keep_pop = keepers, inc_missing = TRUE, out_file = "SmartSNP_infile.txt")
#' vcf2smartsnp(my_vcf, ind_pop, keepers, out_file = "SmartSNP_infile.txt")
#' vcf2smartsnp(my_vcf, ind_pop, keepers)

vcf2smartsnp <- function (vcf, ind_pop, keep_pop, inc_missing = TRUE,
                          out_file = "smartsnp_infile.txt")
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
  gt <- vcfR::extract.gt(vcf2, return.alleles = FALSE, convertNA = TRUE) %>%
    tibble::as_tibble()
  gt_table <- arrow::as_arrow_table(gt) %>%
    dplyr::mutate(across(everything(), ~ dplyr::case_when(
      is.na(.) ~ "9",
      . == "0/0" | . == "0|0" ~ "0",
      . == "1/1" | . == "1|1" ~ "2",
      . == "0/1" | . == "0|1" | . == "1/0" | . == "1|0" ~ "1",
      TRUE ~ as.character(.)
    ))) %>%
    dplyr::collect() %>%
    as.data.frame()

  utils::write.table(gt_table, file = out_file, quote = FALSE, sep = " ", row.names = FALSE, col.names = TRUE, append = FALSE)

  invisible(vcf)
}
