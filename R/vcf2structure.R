#' @title vcf2structure
#' @description converts vcfR format data to Structure or FastStructure infile
#' @description in part based on vcfR2migrate function (vcfR package)
#' @author Tomas Hrbek December 2020
#'
#' @param vcf -> vcfR object
#' @param ind_pop -> population assignment of individuals in vcf (factor)
#' @param keep_pop -> population(s) of interest to include in Structure infile (factor)
#' @param inc_missing -> include missing data (logical)
#' @param out_file -> name of file to output (Structure infile)
#' @param method -> Structure or FastStructure format
#' @export
#' @return nothing
#'
#' @details
#' This function converts the vcfR object to a Structure or FastStructure formatted input file
#' The function will remove indels, and multiallelic loci, and optionally loci with missing data
#'
#' @examples
#' vcf2structure(vcf = my_vcf, ind_pop = ind_pop, keep_pop = keepers, inc_missing = TRUE, out_file = "Structure_infile.str", method = "S")
#' vcf2structure(my_vcf, ind_pop, keepers, out_file = "Structure_infile.str")
#' vcf2structure(my_vcf, ind_pop, keepers)
#'

vcf2structure <-function (vcf, ind_pop, keep_pop, inc_missing = TRUE, out_file = "structure.str", method = "S")
{
  method <- match.arg(method, c("S", "F"), several.ok = FALSE)
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
  vcf_list <- lapply(keep_pop, function(x) {
    vcf[, c(TRUE, x == ind_pop)]
  })
  names(vcf_list) <- keep_pop
  pop_list <- vector(mode = "list", length = length(vcf_list))
  names(pop_list) <- names(vcf_list)

  # process genotype data
  for (i in seq_along(vcf_list)) {
    gt <- vcfR::extract.gt(vcf_list[[i]], return.alleles = FALSE, convertNA = TRUE) %>%
      tibble::as_tibble()

    # split alleles
    allele1 <- arrow::as_arrow_table(gt) %>%
      dplyr::mutate(across(everything(), ~ substr(., 1, 1))) %>%
      dplyr::mutate(across(everything(), ~ if_else(is.na(.), "-9", .))) %>%
      dplyr::collect() %>%
      t() %>%
      as.data.frame()
    rownames(allele1) <- paste0(rownames(allele1), "_1")

    allele2 <- arrow::as_arrow_table(gt) %>%
      dplyr::mutate(across(everything(), ~ substr(., 3, 3))) %>%
      dplyr::mutate(across(everything(), ~ if_else(is.na(.), "-9", .))) %>%
      dplyr::collect() %>%
      t() %>%
      as.data.frame()
    rownames(allele2) <- paste0(rownames(allele2), "_2")

    pop_list[[i]][[1]] <- allele1
    pop_list[[i]][[2]] <- allele2
  }

  # remove existing output file
  if (file.exists(out_file)) {
    file.remove(out_file)
  }

  # write to output file
  if (method == "S") {
    for (i in seq_along(pop_list)) {
      for (j in 1:nrow(pop_list[[i]][[1]])) {
        utils::write.table(t(c(names(pop_list[[i]][[1]][j, 1]), i, pop_list[[i]][[1]][j, ])), file = out_file,
                           append = TRUE, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
        utils::write.table(t(c(names(pop_list[[i]][[2]][j, 1]), i, pop_list[[i]][[2]][j, ])), file = out_file,
                           append = TRUE, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
      }
    }
  }
  else if (method == "F") {
    fill <- rep(0, 4)
    for (i in seq_along(pop_list)) {
      for (j in 1:nrow(pop_list[[i]][[1]])) {
        utils::write.table(t(c(names(pop_list[[i]][[1]][j, 1]), i, fill, pop_list[[i]][[1]][j, ])), file = out_file,
                           append = TRUE, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
        utils::write.table(t(c(names(pop_list[[i]][[2]][j, 1]), i, fill, pop_list[[i]][[2]][j, ])), file = out_file,
                           append = TRUE, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
      }
    }
  }
  else {
    stop("You should never get here!\nMethods are S or F.")
  }

  invisible(vcf)
}
