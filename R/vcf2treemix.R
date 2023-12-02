#' @title vcf2treemix
#' @description converts vcfR format data to TreeMix infile
#' @description in part based on vcfR2migrate function (vcfR package)
#' @author Tomas Hrbek November 2022
#'
#' @param vcf -> vcfR object
#' @param ind_pop -> population assignment of individuals in vcf (factor)
#' @param keep_pop -> population(s) of interest to include in TreeMix infile (factor)
#' @param inc_missing -> include missing data (logical)
#' @param out_file -> name of file to output (TreeMix infile)
#' @export
#' @return nothing
#'
#' @details
#' This function converts the vcfR object to a TreeMix formatted input file
#' The function will remove indels, and multiallelic loci, and optionally loci with missing data
#'
#' @examples
#' vcf2treemix(vcf = my_vcf, ind_pop = ind_pop, keep_pop = keepers, inc_missing = TRUE, out_file = "TreeMix_infile.txt")
#' vcf2treemix(my_vcf, ind_pop, keepers, out_file = "TreeMix_infile.txt")
#' vcf2treemix(my_vcf, ind_pop, keepers)
#'

vcf2treemix <- function (vcf, ind_pop, keep_pop, inc_missing = TRUE,
                         out_file = "treemix_infile.txt")
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
  vcf_list <- lapply(keep_pop, function(x) {
    vcf[, c(TRUE, x == ind_pop)]
  })
  names(vcf_list) <- keep_pop

  for (i in 1:length(vcf_list)) {
    gt <- vcfR::extract.gt(vcf_list[[i]], return.alleles = F, convertNA = T) #convertNA not working here
    allele1 <- apply(gt, MARGIN = 2, function(x){substr(x, 1, 1)})
    allele2 <- apply(gt, MARGIN = 2, function(x){substr(x, 3, 3)})

    REF <- apply(cbind(allele1, allele2), 1, function(x){sum(x == 0, na.rm = TRUE)})
    ALT <- apply(cbind(allele1, allele2), 1, function(x){sum(x == 1, na.rm = TRUE)})
    GEN <- apply(cbind(REF, ALT), 1, function(x){glue::glue_collapse(x, sep = ",")})

    if (i == 1) {
      df <- GEN
    } else {
      df <- cbind(df, GEN)
    }
  }

  colnames(df) <- keep_pop
  utils::write.table(df, file = out_file, quote = FALSE, sep = " ", col.names = TRUE, row.names = FALSE)

  return(invisible(NULL))
}
