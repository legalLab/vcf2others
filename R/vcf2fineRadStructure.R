#' @title vcf2fineRadStructure
#' @description converts vcfR format data to fineRadStructure infile
#' @description in part based on vcfR2genepop function
#' @author Tomas Hrbek July 2022
#'
#' @param vcf -> vcfR object
#' @param ind_pop -> population assignment of individuals in vcf (factor)
#' @param keep_pop -> population(s) of interest to include in fineRadStructure infile (factor)
#' @param inc_missing -> include missing data (logical)
#' @param out_file -> name of file to output (fineRadStructure infile)
#' @export
#' @return nothing
#'
#' @importFrom magrittr %>%
#'
#' @details
#' This function converts the vcfR object to a fineRadStructure formatted input file
#' The function expects multiallelic loci
#'
#' @examples
#' vcf2fineRadStructure(vcf = my_vcf, ind_pop = ind_pop, keep_pop = keepers, inc_missing = TRUE, out_file = "fineRadStructure_infile.txt")
#' vcf2fineRadStructure(my_vcf, ind_pop, keepers, out_file = "fineRadStructure_infile.txt")
#' vcf2fineRadStructure(my_vcf, ind_pop, keepers)
#'

vcf2fineRadStructure <- function (vcf, ind_pop, keep_pop, inc_missing = TRUE,
                                  out_file = "fineRadStructure_infile.txt")
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
  if (inc_missing == FALSE) {
    gt <- vcfR::extract.gt(vcf, convertNA = TRUE)
    vcf <- vcf[!rowSums((is.na(gt))), ]
  }
  vcf2 <- vcf_extract_pops(vcf, ind_pop, keep_pop)
  # extract chromosome names
  chrom <- stringr::str_extract(vcfR::getID(vcf2), "^[^_]*") %>%
    unique()
  # get list of SNPs by chromosomes
  snp_list <- lapply(chrom, function(x){vcf2[x == stringr::str_extract(vcfR::getID(vcf2), "^[^_]*"), ]})

  gt <- vcfR::extract.gt(vcf2, return.alleles = FALSE, convertNA = TRUE)
  # write IDs to file
  utils::write.table(gt[0,], quote = FALSE, sep = "\t", file = out_file)

  for (i in seq_along(snp_list)) {
    # Extract genotype data
    gt <- vcfR::extract.gt(snp_list[[i]], return.alleles = TRUE, convertNA = TRUE) %>%
      tibble::as_tibble()

    allele1 <- arrow::as_arrow_table(gt) %>%
      dplyr::mutate(across(everything(), ~ if_else(. == ".", NA_character_, .))) %>%
      dplyr::mutate(across(everything(), ~ substr(., 1, 1))) %>%
      dplyr::collect() %>%
      apply(2, function(x) glue::glue_collapse(x, sep = ""))

    allele2 <- arrow::as_arrow_table(gt) %>%
      dplyr::mutate(across(everything(), ~ if_else(. == ".", NA_character_, .))) %>%
      dplyr::mutate(across(everything(), ~ substr(., 3, 3))) %>%
      dplyr::collect() %>%
      apply(2, function(x) glue::glue_collapse(x, sep = ""))

    genotype <- apply(data.frame(allele1, allele2), 1, function(x) glue::glue_collapse(x, sep = "/")) %>%
      ifelse(is.na(.), "", .)

    utils::write.table(genotype, quote = FALSE, col.names = FALSE, row.names = FALSE,
                       sep = "\t", file = out_file, append = TRUE)
  }

  invisible(vcf)
}
