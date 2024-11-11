#' @title vcf2eigenstrat
#' @description converts vcfR format data to Eigenstrat infiles
#' @description in part based on vcfR2migrate function (vcfR package)
#' @author Tomas Hrbek December 2022
#'
#' @param vcf -> vcfR object
#' @param ind_pop -> population assignment of individuals in vcf (factor)
#' @param keep_pop -> population(s) of interest to include in Eigenstrat infiles (factor)
#' @param inc_missing -> include missing data (logical)
#' @param out_file -> name of file to output (Eigenstrat infiles)
#' @export
#' @return nothing
#'
#' @details
#' This function converts the vcfR object to a Eigenstrat formatted input files
#' When list of sexes is not provided, lists all individuals as unknown
#' When relative position on chromosome (cM distance or similar) is not provides, list as 0
#' The function will remove indels, and multiallelic loci, and optionally loci with missing data
#'
#' @examples
#' vcf2eigenstrat(vcf = my_vcf, ind_pop = ind_pop, keep_pop = keepers, sex = list_of_sex, rel_pos = marker_cM_map, inc_missing = TRUE, out_file = "Eigenstrat")
#' vcf2eigenstrat(my_vcf, ind_pop, keepers, out_file = "Eigenstrat")
#' vcf2eigenstrat(my_vcf, ind_pop, keepers)
#'

vcf2eigenstrat <- function (vcf, ind_pop, keep_pop, sex = "U", rel_pos = 0, inc_missing = TRUE,
                            out_file = "eigenstrat_infile")
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
      TRUE ~ .
    ))) %>%
    dplyr::collect()

  df_geno <- apply(as.matrix(gt_table), 1, function(x) glue::glue_collapse(x))
  utils::write.table(df_geno, file = paste(out_file, ".geno", sep = ""),
                     quote = FALSE, sep = " ", col.names = FALSE, row.names = FALSE)

  df_ind <- cbind(colnames(gt), sex, as.character(ind_pop))
  utils::write.table(df_ind, file = paste(out_file, ".ind", sep = ""),
                     quote = FALSE, sep = "\t", col.names = FALSE, row.names = FALSE)

  #rel_pos <- vcfR::getPOS(vcf2) / 1000  # Calculate relative positions - fix
  df_snp <- cbind(vcfR::getID(vcf2), vcfR::getCHROM(vcf2), rel_pos,
                  vcfR::getPOS(vcf2), vcfR::getREF(vcf2), vcfR::getALT(vcf2))
  utils::write.table(df_snp, file = paste(out_file, ".snp", sep = ""),
                     quote = FALSE, sep = "\t", col.names = FALSE, row.names = FALSE)

  invisible(vcf)
}
