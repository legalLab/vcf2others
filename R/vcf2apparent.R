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
  vcf <- vcfR::extract.indels(vcf, return.indels = F)
  vcf <- vcf[vcfR::is.biallelic(vcf), ]
  if (inc_missing == FALSE) {
    gt <- vcfR::extract.gt(vcf, convertNA = T)
    vcf <- vcf[!rowSums(is.na(gt)), ]
  }
  vcf2 <- vcf_extract_pops(vcf, ind_pop, keep_pop)
  gt <- vcfR::extract.gt(vcf2, return.alleles = T, convertNA = T)
  gt[gt == "."] <- "-/-"
  gt[gt == "A/A" | gt == "A|A"] <- "A/A"
  gt[gt == "A/C" | gt == "C/A" | gt == "A|C" | gt == "C|A"] <- "A/C"
  gt[gt == "A/G" | gt == "G/A" | gt == "A|G" | gt == "G|A"] <- "A/G"
  gt[gt == "A/T" | gt == "T/A" | gt == "A|T" | gt == "T|A"] <- "A/T"
  gt[gt == "C/C" | gt == "C|C"] <- "C/C"
  gt[gt == "C/G" | gt == "G/C" | gt == "C|G" | gt == "G|C"] <- "C/G"
  gt[gt == "C/T" | gt == "T/C" | gt == "C|T" | gt == "T|C"] <- "C/T"
  gt[gt == "G/G" | gt == "G|G"] <- "G/G"
  gt[gt == "G/T" | gt == "T/G" | gt == "G|T" | gt == "T|G"] <- "G/T"
  gt[gt == "T/T" | gt == "T|T"] <- "T/T"
  gt <- t(gt)
  gt <- cbind(key, gt)

  #suppressWarnings(utils::write.table(gt[0,], file = out_file, quote = FALSE, sep = "\t", col.names = TRUE, append = TRUE))

  utils::write.table(gt, file = out_file, quote = FALSE, sep = "\t", col.names = FALSE, append = FALSE)

  return(invisible(NULL))
}
