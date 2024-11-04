#' @title get_vcf_group_info
#' @description associate samples with group info based on "strata" file
#' @description "strata" is a tsv file with 2+ columns, with 2 columns name id and pop
#' @author Tomas Hrbek July 2022
#'
#' @param vcf -> vcfR object
#' @param data_path -> path to data (where strata and groups are located)
#' @param strt -> file with 2+ columns, one id column and one pop column (tsv file, header = TRUE)
#' @param grps -> one column file listing groups to include, column name pop (header = TRUE)
#' @export
#' @return list of factors of sample-to-group assignments and groups
#'
#' @importFrom magrittr %>%
#'
#' @details
#' This function creates list of factors of sample-to-group assignments and groups based on a strata and group file
#'
#' @examples
#' get_vcf_group_info(vcf = my_vcf, data_path = data_path, strt = "strata", grps = "groups")
#' get_vcf_group_info(my_vcf, data_path, strt = "strata", grps = "groups")
#' get_vcf_group_info(my_vcf, data_path)
#'

get_vcf_group_info <- function(vcf, data_path, strt = "strata", grps = "groups") {
  # read all sample names in vcf
  vcf_names <- colnames(vcf@gt)[-1] %>%
    tibble::as_tibble() %>%
    dplyr::rename(id = 1) %>%
    dplyr::mutate(id = as.character(id))

  # read sample to group assignment
  strata <- read.table(paste0(data_path, strt), header = TRUE) %>%
    tibble::as_tibble() %>%
    dplyr::mutate(id = as.character(id))

  # assign samples in vcf to groups
  strata <- vcf_names %>%
    dplyr::left_join(strata)

  # check if all samples in vcf are assigned to groups
  if (any(is.na(strata$pop) == TRUE)) {
    print("vcf has individuals not assigned to a group")
    print(strata[is.na(strata$pop) == TRUE,])
  } else {
    indiv_group <- as.factor(strata$pop)
  }

  # read groups to be used; filter only those that are actually present
  groups <- read.table(paste0(data_path, grps), header = TRUE)$pop %>%
    .[. %in% unique(strata$pop)] %>%
    as.factor()

  return(list(indiv_group = indiv_group, groups = groups))
}
