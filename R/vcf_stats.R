#' @title vcf_stats
#' @description calculates basic stats of each samples from vcfR format data
#' @description high heterozygosity is indicative of potential contamination
#' @description average read depth per individual
#' @description missing data per individual
#' 
#' @author Tomas Hrbek August 2024
#'
#' @param vcf -> vcfR object
#' @export
#' @return table of statistics
#'
#' @details
#' This function calculates average read depth, heterozygosity
#' number of heterozygotes, number of reference and alternative homozygotes, 
#' missing data and total number SNPs of each sample in an vcfR object
#'
#' @examples
#' vcf_stats(vcf = my_vcf)
#' vcf_stats(my_vcf)
#'

vcf_stats <- function(vcf, res_path, project) {
  vcf <- vcfR::extract.indels(vcf, return.indels = F)
  vcf <- vcf[vcfR::is.biallelic(vcf), ]
  dp <- vcfR::extract.gt(vcf, element = "DP", as.numeric = TRUE)
  gt <- vcfR::extract.gt(vcf, return.alleles = F, convertNA = T)
  samples <- colnames(gt)
  gt[is.na(gt)] <- "NA"
  gt[gt == "0/0" | gt == "0|0"] <- "0"
  gt[gt == "1/1" | gt == "1|1"] <- "1"
  gt[gt == "0/1" | gt == "0|1" | gt == "1/0" | gt == "1|0"] <- "2"

  homo_ref <- apply(gt, MARGIN = 2, function(x){sum(x == 0 , na.rm = TRUE)})
  homo_alt <- apply(gt, MARGIN = 2, function(x){sum(x == 1, na.rm = TRUE)})
  hetero <- apply(gt, MARGIN = 2, function(x){sum(x == 2, na.rm = TRUE)})
  miss <- apply(gt, MARGIN = 2, function(x){sum(x == "NA")})
  total <- homo_ref+homo_alt+hetero
  loc_n <- apply(gt, MARGIN = 2, function(x){length(x)})
  depth <- apply(dp, MARGIN = 2, function(x){mean(x)})
  
  stats <- rbind(samples, depth, hetero/total, hetero, (homo_ref+homo_alt), 
                 homo_ref, homo_alt, miss/loc_n, miss, total, loc_n) %>%
    t() %>%
    as.data.frame()
  colnames(stats) <- c('id', 'read_depth', 'heterozygosity', 'heterozygotes', 'homozygotes',
                       'homozygotes_ref', 'homozygotes_alt', 'missing_p', 'missing', 
                      'non_missing', 'total')
  
  write.table(stats, file = paste0(res_path, project, "stats.csv"), row.names = FALSE, quote = FALSE, sep = ",")
  
  # invisible return of the first argument so function can be used in a pipe
  invisible(vcf)
}
