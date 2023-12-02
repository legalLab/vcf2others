#' @title vcf2genlight
#' @description converts vcfR format data to Genlight infile
#' @description in part based on vcfR2genlight function (vcfR package)
#' @author Tomas Hrbek December 2022
#'
#' @param vcf -> vcfR object
#' @param ind_pop -> population assignment of individuals in vcf (factor)
#' @param keep_pop -> population(s) of interest to include in Genlight infile (factor)
#' @param inc_missing -> include missing data (logical)
#' @param ploidy -> ploidy level (default 2)
#' @param save -> save to file (logical - default TRUE) or just return Genlight object
#' @export
#' @return Genlight object
#'
#' @details
#' This function converts the vcfR object to a Genlight formatted input file
#' This function labels populations.
#' The function will remove indels, and multiallelic loci, and optionally loci with missing data
#'
#' @examples
#' vcf2genlight(vcf = my_vcf, ind_pop = ind_pop, keep_pop = keepers, ploidy = 2, inc_missing = TRUE, out_file = "Genlight_infile")
#' vcf2genlight(my_vcf, ind_pop, keepers, out_file = "Genlight_infile")
#' vcf2genlight(my_vcf, ind_pop, keepers)

vcf2genlight <- function (vcf, ind_pop, keep_pop, ploidy = 2, inc_missing = TRUE,
                          save = FALSE, out_file = "genlight_infile")
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
  if (any(is.na(vcfR::getID(vcf2)))) {
    vcf2 <- vcfR::addID(vcf2)
  }
  gt <- vcfR::extract.gt(vcf2, return.alleles = F, convertNA = T)
  gt[is.na(gt)] <- "NA"
  gt[gt == "0/0" | gt == "0|0"] <- "0"
  gt[gt == "1/1" | gt == "1|1"] <- "2"
  gt[gt == "0/1" | gt == "0|1" | gt == "1/0" | gt == "1|0"] <- "1"

  # create genelight object
  x <- suppressWarnings(new("genlight", t(gt)))
  adegenet::chromosome(x) <- vcfR::getCHROM(vcf2)
  adegenet::position(x) <- vcfR::getPOS(vcf2)
  adegenet::pop(x) <- ind_pop[ind_pop %in% keep_pop]
  adegenet::ploidy(x) <- ploidy
  adegenet::strata(x) <- data.frame(groups = adegenet::pop(x))

  if (save == TRUE) save(x, file = out_file)

  return(x)
}
