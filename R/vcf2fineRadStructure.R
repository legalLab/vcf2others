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
  vcf <- vcfR::extract.indels(vcf, return.indels = F)
  if (inc_missing == FALSE) {
    gt <- vcfR::extract.gt(vcf, convertNA = T)
    vcf <- vcf[!rowSums((is.na(gt))), ]
  }
  vcf2 <- vcf_extract_pops(vcf, ind_pop, keep_pop)
  # extract chromosome names
  chrom <- stringr::str_extract(vcfR::getID(vcf2), "^[^_]*") %>%
    unique()
  # get list of SNPs by chromosomes
  snp_list <- lapply(chrom, function(x){vcf2[x == stringr::str_extract(vcfR::getID(vcf2), "^[^_]*"), ]})

  gt <- vcfR::extract.gt(vcf2, return.alleles = F, convertNA = T)
  # write IDs to file
  utils::write.table(gt[0,], quote = FALSE, sep = "\t", file = out_file)

  for (i in 1:length(snp_list)) {
    gt <- vcfR::extract.gt(snp_list[[i]], return.alleles = T, convertNA = T) #convertNA not working here

    allele1 <- apply(gt, MARGIN = 2, function(x){substr(x, 1, 1)})
    rownames(allele1) <- NULL
    allele1 <- t(allele1)
    allele1[allele1 == "."] <- NA
    allele1a <- apply(allele1, MARGIN = 1, function(x){glue::glue_collapse(x)})
    allele2 <- apply(gt, MARGIN = 2, function(x){substr(x, 3, 3)})
    rownames(allele2) <- NULL
    allele2 <- t(allele2)
    allele2[allele2 == "."] <- NA
    allele2a <- apply(allele2, MARGIN = 1, function(x){glue::glue_collapse(x)})
    allele <- data.frame(allele1a, allele2a)
    allele$frs <- if_else(allele$allele1a == allele$allele2a,
                          allele$allele1a, paste0(allele$allele1a, "/", allele$allele2a))
    allele$frs[is.na(allele$frs)] <- ""

    utils::write.table(t(allele$frs), quote = FALSE, col.names = FALSE, row.names = FALSE,
                sep = "\t", file = out_file, append = TRUE)
  }

  return(invisible(NULL))
}
