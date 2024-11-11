#' @title vcf2bayescan
#' @description converts vcfR format data to Bayescan infile
#' @description in part based on vcfR2migrate function (vcfR package)
#' @author Tomas Hrbek November 2022
#'
#' @param vcf -> vcfR object
#' @param ind_pop -> population assignment of individuals in vcf (factor)
#' @param keep_pop -> population(s) of interest to include in Bayescan infile (factor)
#' @param inc_missing -> include missing data (logical)
#' @param out_file -> name of file to output (Bayescan infile)
#' @export
#' @return nothing
#'
#' @details
#' This function converts the vcfR object to a Bayescan formatted input file
#' The function will remove indels, and multiallelic loci, and optionally loci with missing data
#'
#' @examples
#' vcf2bayescan(vcf = my_vcf, ind_pop = ind_pop, keep_pop = keepers, inc_missing = TRUE, out_file = "Bayescan_infile.txt")
#' vcf2bayescan(my_vcf, ind_pop, keepers, out_file = "Bayescan_infile.txt")
#' vcf2bayescan(my_vcf, ind_pop, keepers)
#'

vcf2bayescan <- function (vcf, ind_pop, keep_pop, inc_missing = TRUE,
                          out_file = "bayescan_infile.txt")
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
    vcf <- vcf[!rowSums((is.na(gt))), ]
  }
  vcf_list <- lapply(keep_pop, function(x) {
    vcf[, c(TRUE, x == ind_pop)]
  })
  names(vcf_list) <- keep_pop
  pop_list <- vector(mode = "list", length = length(vcf_list))
  names(pop_list) <- names(vcf_list)

  for (i in seq_along(vcf_list)) {
    gt <- vcfR::extract.gt(vcf_list[[i]], return.alleles = FALSE, convertNA = TRUE) %>%
      tibble::as_tibble()

    allele1 <- arrow::as_arrow_table(gt) %>%
      dplyr::mutate(across(everything(), ~ substr(., 1, 1))) %>%
      dplyr::collect()

    allele2 <- arrow::as_arrow_table(gt) %>%
      dplyr::mutate(across(everything(), ~ substr(., 3, 3))) %>%
      dplyr::collect()

    pop_list[[i]][[3]] <- rowSums(cbind(allele1, allele2) == "0", na.rm = TRUE)
    pop_list[[i]][[4]] <- rowSums(cbind(allele1, allele2) == "1", na.rm = TRUE)

    pop_list[[i]][[1]] <- pop_list[[i]][[3]] + pop_list[[i]][[4]]
    pop_list[[i]][[2]] <- rep(2, length(pop_list[[i]][[1]]))
  }

  write(paste("[loci]=", nrow(vcf), sep = ""), file = out_file)
  write("", file = out_file, append = TRUE)
  write(paste("[populations]=", length(pop_list), sep = ""), file = out_file, append = TRUE)

  for (i in seq_along(pop_list)) {
    write("", file = out_file, append = TRUE)
    write(paste("[pop]=", i, sep = ""), file = out_file, append = TRUE)
    utils::write.table(pop_list[[i]], file = out_file, quote = FALSE, sep = "\t",
                       col.names = FALSE, row.names = TRUE, append = TRUE)
  }

  invisible(vcf)
}
