#' @title vcf2nexus
#' @description converts vcfR format data to Nexus format infile
#' @author Tomas Hrbek February 2021
#'
#' @param vcf -> vcfR object
#' @param ind_pop -> population assignment of individuals in vcf (factor)
#' @param keep_pop -> population(s) of interest to include in Nexus infile (factor)
#' @param inc_missing -> include missing data (logical)
#' @param out_file -> name of file to output (Nexus infile)
#' @export
#' @return nothing
#'
#' @details
#' This function converts the vcfR object to a Nexus formatted input file
#' The function will remove indels, and multiallelic loci, and optionally loci with missing data
#'
#' @examples
#' vcf2nexus(vcf = my_vcf, ind_pop = ind_pop, keep_pop = keepers, inc_missing = TRUE, out_file = "Nexus_infile.nex")
#' vcf2nexus(my_vcf, ind_pop, keepers, out_file = "Nexus_infile.nex")
#' vcf2nexus(my_vcf, ind_pop, keepers)
#'

vcf2nexus <-function (vcf, ind_pop, keep_pop, inc_missing = TRUE, out_file = "nexus.nex")
{
  #    if (!require(ape)) {
  #        install.packages("ape")
  #    }
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
      . == "." ~ "?",
      . == "A/A" | . == "A|A" ~ "A",
      . == "A/C" | . == "C/A" | . == "A|C" | . == "C|A" ~ "M",
      . == "A/G" | . == "G/A" | . == "A|G" | . == "G|A" ~ "R",
      . == "A/T" | . == "T/A" | . == "A|T" | . == "T|A" ~ "W",
      . == "C/C" | . == "C|C" ~ "C",
      . == "C/G" | . == "G/C" | . == "C|G" | . == "G|C" ~ "S",
      . == "C/T" | . == "T/C" | . == "C|T" | . == "T|C" ~ "Y",
      . == "G/G" | . == "G|G" ~ "G",
      . == "G/T" | . == "T/G" | . == "G|T" | . == "T|G" ~ "K",
      . == "T/T" | . == "T|T" ~ "T",
      TRUE ~ .
    ))) %>%
    dplyr::collect() %>%
    t() %>%
    as.matrix()

  ape::write.nexus.data(gt_table, out_file, format = "DNA", interleaved = FALSE)

  invisible(vcf)
}
