#' @title vcf2fasta
#' @description converts vcfR format data to Fasta format infile
#' @author Tomas Hrbek August 2022
#'
#' @param vcf -> vcfR object
#' @param ind_pop -> population assignment of individuals in vcf (factor)
#' @param keep_pop -> population(s) of interest to include in Fasta infile (factor)
#' @param inc_missing -> include missing data (logical)
#' @param out_file -> name of file to output (Fasta infile)
#' @export
#' @return nothing
#'
#' @details
#' This function converts the vcfR object to a Fasta formatted input file
#' The function will remove indels, and multiallelic loci, and optionally loci with missing data
#'
#' @examples
#' vcf2fasta(vcf = my_vcf, ind_pop = ind_pop, keep_pop = keepers, interleaved = FALSE, inc_missing = TRUE, out_file = "Fasta_infile.fas")
#' vcf2fasta(my_vcf, ind_pop, keepers, out_file = "Fasta_infile.fas")
#' vcf2fasta(my_vcf, ind_pop, keepers)
#'

vcf2fasta <-function (vcf, ind_pop, keep_pop, interleaved = FALSE, inc_missing = TRUE, out_file = "fasta.fas")
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
  vcf <- vcfR::extract.indels(vcf, return.indels = F)
  vcf <- vcf[vcfR::is.biallelic(vcf), ]
  if (inc_missing == FALSE) {
    gt <- vcfR::extract.gt(vcf, convertNA = T)
    vcf <- vcf[!rowSums(is.na(gt)), ]
  }
  vcf2 <- vcf_extract_pops(vcf, ind_pop, keep_pop)
  gt <- vcfR::extract.gt(vcf2, return.alleles = T, convertNA = T)
  gt[gt == "."] <- "?"
  gt[gt == "A/A" | gt == "A|A"] <- "A"
  gt[gt == "A/C" | gt == "C/A" | gt == "A|C" | gt == "C|A"] <- "M"
  gt[gt == "A/G" | gt == "G/A" | gt == "A|G" | gt == "G|A"] <- "R"
  gt[gt == "A/T" | gt == "T/A" | gt == "A|T" | gt == "T|A"] <- "W"
  gt[gt == "C/C" | gt == "C|C"] <- "C"
  gt[gt == "C/G" | gt == "G/C" | gt == "C|G" | gt == "G|C"] <- "S"
  gt[gt == "C/T" | gt == "T/C" | gt == "C|T" | gt == "T|C"] <- "Y"
  gt[gt == "G/G" | gt == "G|G"] <- "G"
  gt[gt == "G/T" | gt == "T/G" | gt == "G|T" | gt == "T|G"] <- "K"
  gt[gt == "T/T" | gt == "T|T"] <- "T"
  gt <- t(gt)

  ifelse(interleaved == FALSE, nbcol <- -1, nbcol <- 6)

  ape::write.dna(gt, out_file, format = "fasta", nbcol = nbcol, colsep = "")

  return(invisible(NULL))
}
