#' @title vcf2snapp
#' @description converts vcfR format data to SNAPP (Nexus) format infile
#' @author Tomas Hrbek December 2020
#'
#' @param vcf -> vcfR object
#' @param ind_pop -> population assignment of individuals in vcf (factor)
#' @param keep_pop -> population(s) of interest to include in SNAPP infile (factor)
#' @param inc_missing -> include missing data (logical)
#' @param out_file -> name of file to output (SNAPP infile)
#' @export
#' @return nothing
#'
#' @details
#' This function converts the vcfR object to a SNAPP (Nexus) formatted input file
#' The function will remove indels, and multiallelic loci, and optionally loci with missing data
#'
#' @examples
#' vcf2snapp(vcf = my_vcf, ind_pop = ind_pop, keep_pop = keepers, inc_missing = TRUE, out_file = "SNAPP_infile.nex")
#' vcf2snapp(my_vcf, ind_pop, keepers, out_file = "SNAPP_infile.nex")
#' vcf2snapp(my_vcf, ind_pop, keepers)
#'

vcf2snapp <-function (vcf, ind_pop, keep_pop, inc_missing = TRUE, out_file = "snapp.nex")
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
  gt <- vcfR::extract.gt(vcf2, return.alleles = FALSE, convertNA = TRUE) %>%
    tibble::as_tibble()
  gt_table <- arrow::as_arrow_table(gt) %>%
    dplyr::mutate(across(everything(), ~ dplyr::case_when(
      is.na(.) ~ "?",
      . == "0/0" | . == "0|0" ~ "0",
      . == "1/1" | . == "1|1" ~ "2",
      . == "0/1" | . == "0|1" | . == "1/0" | . == "1|0" ~ "1",
      TRUE ~ .
    ))) %>%
    dplyr::collect() %>%
    t() %>%
    as.matrix()

  ape::write.nexus.data(gt_table, out_file, format = "standard", interleaved = FALSE)

  # fix symbols in nexus so that snapp.xml file is correctly formated
  fix_symbols <- paste(c("sed -i 's/symbols=\"0123456789\"/symbols=\"012\"/'", out_file), collapse = " ")
  system(fix_symbols)

  invisible(vcf)
}
