#' @title vcf2df
#' @description converts vcfR format data to DataFrame infile
#' @description adds last column of an individual's group origin
#' @description motivated to provide a dataframe for plotting PCA results
#' @description for use by modified PCA function from adegenet
#' @author Tomas Hrbek December 2022
#'
#' @param vcf -> vcfR object
#' @param ind_pop -> population assignment of individuals in vcf (factor)
#' @param keep_pop -> population(s) of interest to include in DataFrame infile (factor)
#' @param alleles -> record data as alleles (GATC) or SNPs (012) (logical)
#' @param inc_missing -> include missing data (logical)
#' @export
#' @return nothing
#'
#' @details
#' This function converts the vcfR object to a DataFrame formatted input file
#' This function adds group (population) membership as last column
#' The function will remove indels, and multiallelic loci, and optionally loci with missing data
#'
#' @examples
#' vcf2df(vcf = my_vcf, ind_pop = ind_pop, keep_pop = keepers, alleles = TRUE, inc_missing = TRUE, out_file = "dataframe_infile.txt")
#' vcf2df(my_vcf, ind_pop, keepers, out_file = "dataframe_infile.txt")
#' vcf2df(my_vcf, ind_pop, keepers)

vcf2df <- function (vcf, ind_pop, keep_pop, alleles = TRUE, inc_missing = TRUE,
                    out_file = "dataframe_infile.txt")
{
  if (class(vcf) != "vcfR") {
    stop(paste("Expecting an object of class vcfR, received a",
               class(vcf), "instead"))
  }
  if (class(ind_pop) != "factor" | class(keep_pop) != "factor") {
    stop(paste("Expecting population vector, received a",
               class(ind_pop), "and", class(keep_pop), "instead"))
  }
  vcf <- vcfR::extract.indels(vcf, return.indels = F)
  vcf <- vcf[vcfR::is.biallelic(vcf), ]
  if (inc_missing == FALSE) {
    gt <- vcfR::extract.gt(vcf, convertNA = T)
    vcf <- vcf[!rowSums(is.na(gt)), ]
  }
  vcf2 <- vcf_extract_pops(vcf, ind_pop, keep_pop)
  if (alleles == TRUE) {
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
      as.data.frame()
  }
  else {
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
      as.data.frame()
  }

  utils::write.table(gt_table, file = out_file, quote = FALSE, sep = "\t", col.names = TRUE, row.names = TRUE)

  df <- cbind(gt_table, as.character(ind_pop))
  colnames(df)[ncol(df)] <- "group"
  df <- df[,-(2:(ncol(df)-1))]

  utils::write.table(df, file = paste(out_file, "groups", sep = "."), quote = FALSE, sep = "\t", col.names = TRUE, row.names = TRUE)

  invisible(vcf)
}
