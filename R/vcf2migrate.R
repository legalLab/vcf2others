#' @title vcf2migrate
#' @description converts vcfR format data to MigrateN infile
#' @description adapted vcfR2migrate function (vcfR package) to allow for inclusion of missing data in Migrate output
#' @author Tomas Hrbek August 2020
#'
#' @param vcf -> vcfR object
#' @param ind_pop -> population assignment of individuals in vcf (factor)
#' @param keep_pop -> population(s) of interest to include in MigrateN infile (factor)
#' @param inc_missing -> include missing data (logical)
#' @param out_file -> name of file to output (MigrateN infile)
#' @param method -> classic or het format
#' @export
#' @return nothing
#'
#' @details
#' This function converts the vcfR object to a MigrateN formatted input file
#' The function will remove indels, and multiallelic loci, and optionally loci with missing data
#'
#' @examples
#' vcf2migrate(vcf = my_vcf, ind_pop = ind_pop, keep_pop = keepers, inc_missing = TRUE, out_file = "MigrateN_infile.txt")
#' vcf2migrate(my_vcf, ind_pop, keepers, out_file = "MigrateN_infile.txt")
#' vcf2migrate(my_vcf, ind_pop, keepers)
#'

vcf2migrate <- function (vcf, ind_pop, keep_pop, inc_missing = TRUE,
                         out_file = "migrateN_infile.txt", method = "N")
{
  method <- match.arg(method, c("N", "H"), several.ok = FALSE)
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

  if (method == "N") {
    myHeader <- c("N", length(vcf_list), nrow(vcf_list[[1]]))
    pop_list <- vector(mode = "list", length = length(vcf_list))
    names(pop_list) <- names(vcf_list)

    for (i in seq_along(vcf_list)) {
      gt <- vcfR::extract.gt(vcf_list[[i]], return.alleles = TRUE, convertNA = TRUE) %>%
        tibble::as_tibble()

      allele1 <- arrow::as_arrow_table(gt) %>%
        dplyr::mutate(across(everything(), ~ if_else(. == ".", "?/?", .))) %>%
        dplyr::mutate(across(everything(), ~ substr(., 1, 1))) %>%
        dplyr::collect() %>%
        t() %>%
        as.data.frame()
      rownames(allele1) <- paste0(rownames(allele1), "_1")

      allele2 <- arrow::as_arrow_table(gt) %>%
        dplyr::mutate(across(everything(), ~ if_else(. == ".", "?/?", .))) %>%
        dplyr::mutate(across(everything(), ~ substr(., 3, 3))) %>%
        dplyr::collect() %>%
        t() %>%
        as.data.frame()
      rownames(allele2) <- paste0(rownames(allele2), "_2")

      pop_list[[i]][[1]] <- allele1
      pop_list[[i]][[2]] <- allele2
    }

    # write header and metadata
    write(myHeader, file = out_file, ncolumns = length(myHeader), sep = "\t")
    write(rep(1, times = ncol(pop_list[[1]][[1]])), file = out_file,
          ncolumns = ncol(pop_list[[1]][[1]]), append = TRUE, sep = "\t")

    # write allele data
    for (i in seq_along(pop_list)) {
      popName <- c(2 * nrow(pop_list[[i]][[1]]), names(pop_list)[i])
      write(popName, file = out_file, ncolumns = length(popName), append = TRUE, sep = "\t")

      for (j in seq_len(ncol(pop_list[[i]][[1]]))) {
        utils::write.table(pop_list[[i]][[1]][, j], file = out_file,
                           append = TRUE, quote = FALSE, sep = "\t", row.names = TRUE, col.names = FALSE)
        utils::write.table(pop_list[[i]][[2]][, j], file = out_file,
                           append = TRUE, quote = FALSE, sep = "\t", row.names = TRUE, col.names = FALSE)
      }
    }
  }
  else if (method == "H") {
    myHeader <- c("H", length(vcf_list), nrow(vcf_list[[1]]))
    pop_list <- vector(mode = "list", length = length(vcf_list))
    names(pop_list) <- names(vcf_list)

    for (i in seq_along(vcf_list)) {
      myMat <- matrix(nrow = nrow(vcf_list[[i]]), ncol = 6)
      var_info <- tibble::as_tibble(vcf_list[[i]]@fix[, 1:2, drop = FALSE])
      var_info$mask <- TRUE

      gt <- vcfR::extract.gt(vcf_list[[i]], return.alleles = FALSE, convertNA = TRUE) %>%
        tibble::as_tibble()
      popSum <- arrow::as_arrow_table(gt) %>%
        dplyr::mutate(across(everything(), ~ dplyr::case_when(
          . == "0/0" | . == "0|0" ~ "REF",
          . == "1/1" | . == "1|1" ~ "ALT",
          . == "0/1" | . == "1/0" | . == "0|1" | . == "1|0" ~ "HET",
          TRUE ~ .
        ))) %>%
        dplyr::summarise(
          allele_counts = dplyr::across(everything(), ~ paste0(
            sum(. == "REF", na.rm = TRUE), ",",
            sum(. == "ALT", na.rm = TRUE)
          ))
        ) %>%
        dplyr::collect() %>%
        as.data.frame()

      myMat[, 1] <- paste(vcf_list[[i]]@fix[, "CHROM"], vcf_list[[i]]@fix[, "POS"], sep = "_")
      myMat[, 2] <- vcf_list[[i]]@fix[, "REF"]
      myMat[, 4] <- vcf_list[[i]]@fix[, "ALT"]
      allele_counts <- strsplit(as.character(popSum$allele_counts), split = ",", fixed = TRUE)
      myMat[, 3] <- sapply(allele_counts, function(x) ifelse(length(x) >= 1, x[1], 0))
      myMat[, 5] <- sapply(allele_counts, function(x) ifelse(length(x) >= 2, x[2], 0))
      myMat[, 3][is.na(myMat[, 3])] <- 0
      myMat[, 5][is.na(myMat[, 5])] <- 0
      myMat[, 6] <- as.numeric(myMat[, 3]) + as.numeric(myMat[, 5])

      pop_list[[i]] <- myMat
    }

    write(myHeader, file = out_file, ncolumns = length(myHeader), sep = "\t")

    for (i in seq_along(pop_list)) {
      popName <- c(pop_list[[i]][1, 6], names(pop_list[i]))
      write(popName, file = out_file, ncolumns = length(popName), append = TRUE, sep = "\t")
      utils::write.table(pop_list[[i]], file = out_file, append = TRUE, quote = FALSE,
                         sep = "\t", row.names = FALSE, col.names = FALSE)
    }
  }
  else {
    stop("You should never get here!\nMethods are N or H.")
  }

  invisible(vcf)
}
