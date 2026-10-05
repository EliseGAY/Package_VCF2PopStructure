#' filter_vcf_by_dp : Filter Inds and Pos by DP Thresholds
#'
#' This function filters a VCF file based on the depth of coverage (DP) for both positions (pos) and individuals (ind). 
#' It removes positions and individuals whose DP is outside the specified minimum and maximum thresholds.
#'
#' @param vcf A VCF object (vcfR format) to be filtered.
#' @param min_dp_pos (double) The minimum allowed DP for positions.
#' @param max_dp_pos (double) The maximum allowed DP for positions.
#' @param min_dp_ind (double) The minimum allowed DP for individuals.
#' @param max_dp_ind (double) The maximum allowed DP for individuals.
#'
#' @return A filtered VCF object where positions and individuals are removed based on the DP thresholds.
#' @import vcfR
#' @export
#'
#' @examples
#' # Filter a VCF file based on DP thresholds for both pos and ind
#' filtered_vcf <- filter_vcf_by_dp(vcf_file, 10, 100, 20, 200)

filter_vcf_by_dp <- function(vcf, min_dp_pos, max_dp_pos, min_dp_ind, max_dp_ind) {
  
  # Extract the depth of coverage (DP) for all positions
  dp_data <- extract.gt(vcf, element = 'DP', as.numeric = TRUE)
  
  # Calculate the mean DP for each ind (columns)
  mean_dp_ind <- colMeans(dp_data, na.rm = TRUE)
  
  # Filter ind based on the DP thresholds
  valid_inds <- names(mean_dp_ind[mean_dp_ind >= min_dp_ind & 
                                    mean_dp_ind <= max_dp_ind])
  
  # Filter the VCF by ind
  filtered_vcf_ind <- vcf[, c("FORMAT", valid_inds)]
  
  
  # Calculate the mean DP for each position (rows)
  mean_dp_pos <- rowMeans(dp_data, na.rm = TRUE)
  
  # Filter positions based on the DP thresholds
  filtered_vcf_ind_pos <- filtered_vcf_ind[mean_dp_pos >= min_dp_pos & 
                                             mean_dp_pos <= max_dp_pos, ]
  
  
  return(filtered_vcf_ind_pos)
}

#' Filters_na_Ind_Pos : Filters VCF based on missing data rate for individuals and positions
#'
#' This function filters a VCF object by the missing data rate for both individuals (samples)
#' and positions (variants). The user can specify the ploidy level (1 or 2), and the function
#' will compute the missing data rate based on either the `.` (for ploidy 1) or `./.` (for ploidy 2)
#' genotype representation.
#'
#' @param vcf A VCF object (vcfR object) to be filtered.
#' @param rate_na_max_POS Maximum allowed percentage of missing data per position.
#' @param rate_na_max_Ind Maximum allowed percentage of missing data per individual.
#' @param ploidy Ploidy level (1 or 2) to determine the missing data representation (`.` or `./.`).
#' 
#' @return A filtered VCF object.
#' @export

Filters_na_Ind_Pos <- function(vcf, rate_na_max_POS, rate_na_max_Ind, ploidy) {
  
  # Check for valid ploidy value
  if (!ploidy %in% c(1, 2)) {
    stop("Ploidy must be either 1 or 2.")
  }
  
  # Get genotype data from the VCF
  gt <- extract.gt(vcf, element = "GT", mask = FALSE, as.numeric = FALSE, return.alleles = FALSE, convertNA = FALSE)
  
  # Set the missing genotype string based on ploidy
  missing_str <- if (ploidy == 1) "." else "./."
  
  # Calculate missing data rate for positions (rows)
  missing_positions <- rowSums(gt == missing_str)
  rate_na_pos <- (missing_positions / ncol(gt)) * 100
  positions_kept <- which(rate_na_pos <= rate_na_max_POS)
  
  # Filter the VCF by positions
  vcf_filtered_na <- vcf[positions_kept, ]
  
  # Re-extract genotype data for the filtered individuals
  gt_filtered <- extract.gt(vcf_filtered_na, element = "GT", mask = FALSE, as.numeric = FALSE, return.alleles = FALSE, convertNA = FALSE)
  
  # Calculate missing data rate for individuals (columns)
  missing_individuals <- colSums(gt_filtered == missing_str)
  rate_na_ind <- (missing_individuals / nrow(gt)) * 100
  samples_kept <- names(rate_na_ind[rate_na_ind <= rate_na_max_Ind])
  
  # Filter the VCF by individuals
  vcf_filtered_ind <- vcf_filtered_na[, c("FORMAT", samples_kept)]
  
  
  # Return the filtered VCF
  return(vcf_filtered_ind)
}

#' Filters_Het_Ind_Pos : Filter Individuals and Positions by Alternative Rate
#'
#' This function filters a VCF file based on the alternative allele rate both by individuals and by positions. 
#' The function calculates the percentage of heterozygous (alternative) genotypes for each individual and for each position,
#' and filters the VCF file according to the provided thresholds. The function also handles ploidy (1 or 2) by adjusting
#' the genotypes accordingly ("0/1" for ploidy 2, "1" for ploidy 1).
#'
#' @param vfc A VCF object (vcfR format) to be filtered.
#' @param ploidy (integer) The ploidy level of the dataset. Either 1 or 2.
#' @param rate_het_max_POS (double) The maximum allowed alternative allele rate (percent) for positions.
#' @param rate_het_max_Ind (double) The maximum allowed alternative allele rate (percent) for individuals.
#' @param rate_het_min_pos (integer) The minimum number of alternative alleles allowed per position.
#'
#' @return A filtered VCF object where positions and individuals are removed based on the alternative allele rate filters.
#' @import vcfR
#' @export
#'
#' @examples
#' # Filter a VCF file based on alternative rate thresholds for ploidy 2
#' filtered_vcf <- Filters_Het_Ind_Pos(vcf_file, 2, 30.0, 50.0, 0.0)
#'
#' # Filter a VCF file based on alternative rate thresholds for ploidy 1
#' filtered_vcf <- Filters_Het_Ind_Pos(vcf_file, 1,  30.0, 50.0, 0.0)

Filters_Het_Ind_Pos <- function(vcf, ploidy, rate_het_max_POS, rate_het_max_Ind) {
  
  # Validate ploidy value
  if (!ploidy %in% c(1, 2)) {
    stop("Ploidy must be either 1 or 2.")
  }
  
  # Define genotype pattern based on ploidy
  het_pattern <- if (ploidy == 1) "1" else "0/1"
  
  # Extract genotypes from VCF
  gt <- extract.gt(VCF_P1, element = "GT", mask = FALSE, as.numeric = FALSE, return.alleles = FALSE, convertNA = FALSE)
  
  # Calculate alternative allele rate for individuals (columns)
  somme_het_ind <- colSums(gt == het_pattern)
  rate_het_ind <- (somme_het_ind / nrow(gt)) * 100
  samples_kept <- names(rate_het_ind[rate_het_ind <= rate_het_max_Ind])
  
  # Filter VCF by individuals
  vcf_filtered_ind <- VCF_P1[, c("FORMAT", samples_kept)]
  
  # Re-extract genotypes for the filtered VCF
  gt_1 <- extract.gt(vcf_filtered_ind, element = "GT", mask = FALSE, as.numeric = FALSE, return.alleles = FALSE, convertNA = FALSE)
  
  # Calculate alternative allele count for positions (rows)
  somme_het_pos <- rowSums(gt_1 == het_pattern)
  rate_het_pos <- (somme_het_pos / ncol(gt)) * 100
  # Apply minimum threshold for alternative alleles at positions
  pos_kept_min <- which(rate_het_pos >= rate_het_max_POS)
  vcf_filtered_het <- vcf_filtered_ind[pos_kept_min, ]
  
  # Return the filtered VCF
  return(vcf_filtered_het)
}
