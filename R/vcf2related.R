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
  vcf <- vcfR::extract.indels(vcf, return.indels = F)
  vcf <- vcf[vcfR::is.biallelic(vcf), ]
  if (inc_missing == FALSE) {
    gt <- vcfR::extract.gt(vcf, convertNA = T)
    vcf <- vcf[!rowSums((is.na(gt))), ]
  }
  vcf2 <- vcf_extract_pops(vcf, ind_pop, keep_pop)

  gt <- vcfR::extract.gt(vcf2, return.alleles = T, convertNA = T) #convertNA not working here
  gt[gt == "."] <- "NA\tNA"
  gt[gt == "A/A" | gt == "A|A"] <- "01\t01"
  gt[gt == "A/C" | gt == "C/A" | gt == "A|C" | gt == "C|A"] <- "01\t02"
  gt[gt == "A/G" | gt == "G/A" | gt == "A|G" | gt == "G|A"] <- "01\t03"
  gt[gt == "A/T" | gt == "T/A" | gt == "A|T" | gt == "T|A"] <- "01\t04"
  gt[gt == "C/C" | gt == "C|C"] <- "02\t02"
  gt[gt == "C/G" | gt == "G/C" | gt == "C|G" | gt == "G|C"] <- "02\t03"
  gt[gt == "C/T" | gt == "T/C" | gt == "C|T" | gt == "T|C"] <- "02\t04"
  gt[gt == "G/G" | gt == "G|G"] <- "03\t03"
  gt[gt == "G/T" | gt == "T/G" | gt == "G|T" | gt == "T|G"] <- "03\t04"
  gt[gt == "T/T" | gt == "T|T"] <- "04\t04"
  gt <- t(gt)

  utils::write.table(gt, file = out_file, quote = FALSE, sep = "\t", col.names = FALSE, row.names = TRUE)

  return(invisible(NULL))
}
