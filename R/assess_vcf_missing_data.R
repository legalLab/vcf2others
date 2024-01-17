#' @title assess_vcf_missing_data
#' @description quantifying missing data of all samples in VCF
#' @description inspired by https://grunwaldlab.github.io/Population_Genetics_in_R/qc.html
#' @author Tomas Hrbek July 2022
#'
#' @param vcf -> vcfR object
#' @param data_path -> path to data (where strata is located)
#' @param res_path -> path to results (directory for output dataframe and plot)
#' @param project -> project base name (character)
#' @param postfix -> VCF extraction method (character)
#' @param fltr -> VCF filtration method (character)
#' @param species -> species name for plot (character)
#' @param strt -> file with 2+ columns, one id column and one pop column (tsv file, header = TRUE)
#' @export
#' @return dataframe and plot of missing data for each individual in VCF
#'
#' @importFrom magrittr %>%
#'
#' @details
#' This function generates a dataframe of absolute and relative missing data per sample
#' This function plots relative missing data per individual, with individuals sorted by group
#' The postfix and fltr parameters can have empty values ("")
#'
#' @examples
#' assess_vcf_missing_data(vcf = my_vcf, data_path = data_path, res_path = res_path, project = project_name, postfix = extraction_info, fltr = filter_parms, species = species_name, strt = "strata")
#' assess_vcf_missing_data(my_vcf, data_path, res_path, project_name, extraction_info, filter_parms, species_name, strt = "strata")
#' assess_vcf_missing_data(my_vcf, data_path, res_path, project_name, postfix, fltr, species)
#'

assess_vcf_missing_data <- function(vcf, data_path, res_path, project, postfix, fltr, species, strt = "strata") {
  # read sample to group assignment
  strata <- read.table(paste0(data_path, strt), header = TRUE) %>%
    tibble::as_tibble() %>%
    dplyr::mutate(id, id = as.character(id))

  # extract genotypes
  gt <- vcfR::extract.gt(vcf, element = "GT", as.numeric = TRUE)

  # get missing number per individual and % missing
  n_miss <- apply(gt, MARGIN = 2, function(x){sum(is.na(x))})
  p_miss <- n_miss/nrow(gt)

  # organize as dataframe and save
  df <- cbind(p_miss, n_miss) %>%
    tibble::as_tibble() %>%
    dplyr::mutate(n_total = nrow(gt)) %>%
    dplyr::mutate(id = colnames(gt)) %>%
    dplyr::left_join(strata) %>%
    dplyr::relocate(c(id, group = pop))

  write.table(df, file = paste0(res_path, project, postfix, fltr, "_missingness"), row.names = FALSE, quote = FALSE)

  # plot missingness
  ordr <- df %>%
    dplyr::arrange(group) %>%
    dplyr::pull(id)
  p <- ggplot2::ggplot(df, ggplot2::aes(x=factor(id, level=ordr), y=p_miss, fill=group)) +
    ggplot2::geom_bar(stat = "identity") +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = .5, hjust = 1, size = ggplot2::rel(.8))) +
    ggplot2::labs(x = "sample", y = "% missing", title = bquote(paste(italic(.(species)), " ", "(", .(postfix), ")")))

  # output missingness to pdf, svg and png
  ggplot2::ggsave(plot=p, filename=paste0(res_path, project, postfix, fltr, "_missingness.pdf"), width=6, height=4, bg="transparent", limitsize=FALSE)
  ggplot2::ggsave(plot=p, filename=paste0(res_path, project, postfix, fltr, "_missingness.svg"), device="svg", width=6, height=4, bg="transparent", limitsize=FALSE)
  ggplot2::ggsave(plot=p, filename=paste0(res_path, project, postfix, fltr, "_missingness.png"), device="png", width=6, height=4, bg="transparent", limitsize=FALSE)
}
