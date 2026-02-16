#' @param : genotype table, make sure that 'na' are encoded in R readable. In format :
#' 
#' loci_table = data.frame(I=c("0/1" , "0/1" , "1/1", "0/0" , "0/0" , NA, "0/1" , "0/1" , "1/1"),
#' 
#'                        J=c("0/1" , "0/1" , "0/1", NA , "0/0" , NA, "0/1" , "0/1", NA ),
#' 
#'                        K=c("0/1" , "1/1" , "1/1", "0/0" , "0/0" , "1/1", "1/1" , "1/1", "0/0" ))
#'
#' @param vec_pop : list of pop in order of genotype table with the following encoding : 
#' 
#' pop_vec=c("A","A","A", "B", "B", "B", "C", "C", "C")
#' 
#' samples=c("s1","s2","s3", "s4", "s5", "s6", "s7", "s8", "s9")
#' 
#' pop_table_test = data.frame(cbind(pop_vec, samples))

#' @examples
#' loci_table_T = t(loci_table)
#' colnames(loci_table_T) = samples

#' # add info on GT table in 0,1,2 genotype codes
#' loci_table_T_CV = Convert_GT(loci_table_T)
#'   s1 s2 s3 s4 s5 s6 s7 s8 s9
#' I  1  1  2  0  0 NA  1  1  2
#' J  1  1  1 NA  0 NA  1  1 NA
#' K  1  2  2  0  0  2  2  2  0
Convert_GT=function(GT_table){
  if(all(GT_table[!is.na(GT_table)] %in% c("0/1", "0/0", "1/1"))){
    GT_num <- matrix(rep(NA, length(GT_table)), dimnames = list(row.names(GT_table), colnames(GT_table)), 
                     nrow = dim(GT_table)[1], ncol = dim(GT_table)[2])
    
    GT_num[GT_table == "0/0"] <- 0L
    GT_num[GT_table == "0/1"] <- 1L
    GT_num[GT_table == "1/1"] <- 2L
  }
  else{print("only genotype 0/0 , 0/1 , 1/1, NA are allowed")}
  GT_num = as.data.frame(GT_num)
  return(GT_num)
}

FormatLociTable<-function(lociTable){
  lociTable = as.data.frame(lociTable)
  
}

#' PurgeRefInPops : detect and discard homozygous loci (usefull for op subset table)
#' 
#' @param loci_table : genotype table, make sure that 'na' are encoded in R readable. In format :
#' 
#' loci_table (df) : 
#' 
#'     s1 s2 s3 s4 s5 s6 s7 s8 s9
#' 
#' I  1  1  2  0  0 NA  1  1  2
#' 
#' J  1  1  1 NA  0 NA  1  1 NA
#' 
#' K  1  2  2  0  0  2  2  2  0
#' 
#' L  0  0  0  0  0  0  0  0  0
#' 
#' @param pop_table : list of pop in order of genotype table with the following encoding : 
#' 
#' pop_vec=c("A","A","A", "B", "B", "B", "C", "C", "C")
#' 
#' samples=c("s1","s2","s3", "s4", "s5", "s6", "s7", "s8", "s9")
#' 
#' pop_table_test = data.frame(cbind(pop_vec, samples))
#' 
#' @returns subset of loci table
#' 
#' @examples
#' loci_table_T = t(loci_table)
#' colnames(loci_table_T) = samples
#' loci_table_T_CV = Convert_GT(loci_table_T)
#' SubSnp = PurgeRefInPops(loci_table = loci_table_T_CV_ex, pop_table = pop_table_test)
#'   s1 s2 s3 s4 s5 s6 s7 s8 s9
#' I  1  1  2  0  0 NA  1  1  2
#' J  1  1  1 NA  0 NA  1  1 NA
#' K  1  2  2  0  0  2  2  2  0
#' 
PurgeRefInPops <- function(loci_table, pop_table) {
  # check for dataframe loci_table
  # get pop list : col2 = Samples, col1 = populations
  # keep only metadata samples present in loci table
  pop_meta_filtered <- pop_table[pop_table[,2] %in% colnames(loci_table),]
  print(pop_meta_filtered)
  
  # build pop list (population -> samples)
  pop_list <- split(pop_meta_filtered[[2]], pop_meta_filtered[[1]])
  
  # loop on each locus (line)
  is_snp_any_pop <- apply(loci_table, 1, function(gt){
    
    # initiate the presence of SNPs
    pop_with_snps = 0 
    
    # count the nb of genotype element in each pop
    # if = 1 (only one element either "1" or "0") if > 1 there is two different gt
    
    for(popi_sample in pop_list){
      locus_pop_i = gt[popi_sample]
      locus_pop_i_noNA = locus_pop_i[!is.na(locus_pop_i)]
      if (length(unique(locus_pop_i_noNA)) > 1) {
        pop_with_snps = 1
        break   # STOP checking other pops
      }
    }
    # return a bolean to this keep/remove locus 
    pop_with_snps >= 1
  })
  # subset loci table
  pruged_Loci = loci_table[is_snp_any_pop, , drop = FALSE]
  return(pruged_Loci)
}

#' NaByPop : detect and discard loci if they have a NA rate > threshold
#' 
#' @param loci_table : genotype table, make sure that 'na' are encoded in R readable. In format :
#' loci_table (df) : 
#'     s1 s2 s3 s4 s5 s6 s7 s8 s9
#' I  1  1  2  0  0 NA  1  1  2
#' J  1  1  1 NA  0 NA  1  1 NA
#' K  1  2  2  0  0  2  2  2  0
#' L  0  0  0  NA  NA  NA  0  1  0
#' @param pop_table : list of pop in order of genotype table with the following encoding : 
#' pop_vec=c("A","A","A", "B", "B", "B", "C", "C", "C")
#' samples=c("s1","s2","s3", "s4", "s5", "s6", "s7", "s8", "s9")
#' pop_table_test = data.frame(cbind(pop_vec, samples))
#' @returns subset of loci table
#' @examples
#' loci_table_T = t(loci_table)
#' colnames(loci_table_T) = samples
#' loci_table_T_CV = Convert_GT(loci_table_T)
#' SubSnp = NaByPop(loci_table = loci_table_T, pop_table = pop_table_test)
#'   s1 s2 s3 s4 s5 s6 s7 s8 s9
#' I  1  1  2  0  0 NA  1  1  2
#' K  1  2  2  0  0  2  2  2  0
#' 
NaByPop <- function(loci_table, table_pop, threshold = 0.2) {
  
  # get pop list : col2 = Samples, col1 = populations
  # keep only metadata samples present in loci table
  pop_meta_filtered <- table_pop[table_pop[,2] %in% colnames(loci_table),]
  
  # build pop list (population -> samples)
  pop_list <- split(pop_meta_filtered[, 2], pop_meta_filtered[, 1])  
  
  keep_locus <- apply(loci_table, 1, function(gt) {
    all(sapply(pop_list, function(samples) {
      mean(is.na(gt[samples])) <= 0.2
    }))
  })
  
  loci_table[keep_locus, , drop = FALSE]
}

#' getSamplePop : Extract one of pop in a metadata table
#' @param pop_name (char)
#' @param pop_table (dataframe) : contains samples ID in col1 and Pop names in col2
#' @returns list of samples named by pop 
#' @examples 
#' > pop_table_test
#' pop_vec samples
#' 1       A      s1
#' 2       A      s2
#' 3       A      s3
#' 4       B      s4
#' 5       B      s5
#' 6       B      s6
#' 7       C      s7
#' 8       C      s8
#' 9       C      s9
#' > getSamplePop(pop_name = "A", pop_table_test)
#' [1] "s1" "s2" "s3"
#' attr(,"levels")
#' [1] "A"
getSamplePop<-function(pop_name, pop_table){
  pop1 = pop_table[which(pop_table[,1] == pop_name[1]),][,2]
  levels(pop1) = pop_name[1]
  return(pop1)
}

# getSamplePairPop : Extract pair of pop in a metadata table
#' @param Pair_i vector of 2 pop names (char)
#' @param pop_table (dataframe) : contains samples ID in col1 and Pop names in col2
#' @returns list of samples named by pop 
#' @examples 
#'  > pop_table_test
#' pop_vec samples
#' 1       A      s1
#' 2       A      s2
#' 3       A      s3
#' 4       B      s4
#' 5       B      s5
#' 6       B      s6
#' 7       C      s7
#' 8       C      s8
#' 9       C      s9
#' > getSamplePairPop(Pair_i = c("A", "B"), pop_table_test)
#' $A
#' [1] "s1" "s2" "s3"
#'
#' $B
#' [1] "s4" "s5" "s6"
getSamplePairPop<-function(Pair_i, pop_table){
  pop1 = pop_table[which(pop_table[,1]== Pair_i[1]),][,2]
  pop2 = pop_table[which(pop_table[,1]== Pair_i[2]),][,2]
  list_pair = list(pop1, pop2)
  names(list_pair) = c(Pair_i[1], Pair_i[2])
  return(list_pair)
}

#' getSamplePop : Extract one of pop in a metadata table
#' @param pop_vec (char  vec)
#' @returns df of pairwise pop couple (by row)
#' @examples 
#' pairs_pop_sample = getPairwisePop(unique(pop_table_test[,1]))
getPairwisePop<-function(pop_vec){
  pops = unique(pop_vec)
  pairs = list()
  k = 1
  for (i in seq_len(length(pops) - 1)) {
    for (j in (i + 1):length(pops)) {
      pairs[[k]] = c(pops[i], pops[j])
      k = k + 1
    }
  }
  
  do.call(rbind, pairs) |> as.data.frame() # eq. to  rbind(pairs[[1]], pairs[[2]], pairs[[3]])
}

#' getAlleleFreq : get the Alt and Ref allele frequencies in genotype table (0,1,2 encoded)
#' 
#' @param : locus table, make sure that 'NA' are encoded in R readable. 
#' 
#' In format : Obtain by convert_GT()
#' 
#'   s1 s2 s3 s4 s5 s6 s7 s8 s9
#' 
#' I  1  1  2  0  0 NA  1  1  2
#' 
#' @return List of alt and ref allele frequencies in the pop. Na are removed from the sample number
#' 
#' @examples getAlleleFreq(locus)
#' $alt_freq
#' [1] 0.4
#' 
#' $ref_freq
#' [1] 0.6

getAlleleFreq = function(locus){
  ## TODO : Managment and verif of LociPairs table (Na and df)
  
  if(all(locus[!is.na(locus)] %in% c(0, 1, 2, NA))){
    
    tot_alt_freq = as.numeric(rowSums(locus, na.rm = T)) / as.numeric(2*ncol(locus))
    tot_ref_freq = 1-tot_alt_freq
    freq_list = list(tot_alt_freq, tot_ref_freq)
    names(freq_list) = c("alt_freq", "ref_freq")
  }else{stop("Invalid genotype detected")}
  return(freq_list)
}


#' getAltCount : get the Alt and Ref allele count in genotype table (0,1,2 encoded)
#' 
#' @param : locus table, make sure that 'NA' are encoded in R readable. 
#' 
#' @return List of alt allele count in the pop.
#' 
#' @examples getAlleleFreq(locus)
#' $alt_freq
#' [1] 0.4
#' 
#' $ref_freq
#' [1] 0.6

getAlleleCount= function(locus){
  
  if(all(locus[!is.na(locus)] %in% c(0, 1, 2, NA))){
    count_list=list()
    tot_alt_count = as.numeric(rowSums(locus, na.rm = T)) 
    count_list[["alt_count"]] = tot_alt_count
	count_list[["ref_count"]] = (2*ncol(locus)) - count_list$alt_count
  }else{stop("Invalid genotype detected")}
  return(count_list)
}

# getAlleleFreqByPop : Extract pair of pop in a metadata table
#' @param loci_table dataframe of allele count in loci (rows) and in samples (col) from one or several pop
#' @param pop_table (dataframe) : contains samples ID in col1 and Pop names in col2
#' @returns list of allele frequencies for each pop
#' @examples
#' > locus_table
#' s1 s2 s3 s4 s5 s6 s7 s8 s9
#' I  1  1  2  0  0 NA  1  1  2
#' > pop_table_test
#' pop_vec samples
#' 1       A      s1
#' 2       A      s2
#' 3       A      s3
#' 4       B      s4
#' 5       B      s5
#' 6       B      s6
#' 7       C      s7
#' 8       C      s8
#' 9       C      s9
#' > summary(res)
#'         Length Class  Mode
#'pop1 2      -none- list
#'pop2   2      -none- list
#'pop3 2      -none- list
#'> summary(res$pop1)
#'         Length Class  Mode   
#'alt_freq 181043 -none- numeric
#'ref_freq 181043 -none- numeric


getAlleleFreqByPop=function(loci_table, pop_table){
  # sanity check
  stopifnot(is.data.frame(loci_table))
  
  # check that genotypes are valid
  if (!all(loci_table[!is.na(loci_table)] %in% c(0, 1, 2))) {
    stop("Invalid genotype detected in loci_table. Only 0, 1, 2, or NA allowed.")
  }
  ## TODO : Managment and verif of LociPairs table (Na and df)
  samples_2pops <- colnames(loci_table)
  
  # check all samples exist in pop_table
  missing_in_pop_table <- setdiff(samples_2pops, pop_table[,2])
  if (length(missing_in_pop_table) > 0) {
    stop(
      "ERROR: loci_table contains samples not present in pop_table:\n",
      paste(missing_in_pop_table, collapse = ", ")
    )
  }
  
  # get populations present in this subset, sorted for stable order
  pop_level <- sort(unique(pop_table[,1][pop_table[,2] %in% samples_2pops]))
  
  # initiate result list
  res <- list()
   # loop over populations
  for (pop_i in pop_level) {
    # extract sample IDs for this population
    samples_pop_i <- pop_table[pop_table[,1] == pop_i & pop_table[,2] %in% samples_2pops, 2]
    
    # subset loci_table for these samples
    loci_pop_i <- loci_table[, samples_pop_i, drop = FALSE]
    res[[pop_i]] = getAlleleFreq(loci_pop_i)
    }
  return(res)
}



# getAlleleCountByPop : Extract allele count from genotype table by pop use getAlleleCount function 
#' @param loci_table dataframe of allele count in loci (rows) and in samples (col) from one or several pop
#' @param pop_table (dataframe) : contains samples ID in col1 and Pop names in col2
#' @returns list of allele count for each pop
#' @examples
#' > locus_table
#' A tibble: 181,043 × 31
#'   `1622W1_S254` `1623W1_S276` `1624W1_S392` `1625Q_S431` `1625W1_S384` `1626W1_S364` `1627W1_S385` `1629W1_S379` `1681W1_S374` `1682W1_S380` `1684W1_S376` `1686W1_S402`
#' 1             0             0             0            0             0             0             0             0             0             0             0             0
#' 2             0             0             0            0             0             0             0             0             0             0             0             0
#' 3             0             0             0            0             0             0             0             0             0             0             0             0
# ℹ 181,033 more rows
# ℹ 19 more variables 

#' > pop_table_test
#' pop_vec samples
#' 1       A      s1
#' 2       A      s2
#' 3       A      s3
#' 4       B      s4
#' 5       B      s5
#' 6       B      s6
#' 7       C      s7
#' 8       C      s8
#' 9       C      s9

#' by_popCount = getAlleleCountByPop(loci_table = loci_table_T_CV, pop_table = metadata)
#' by_popCount[[1]]$alt_count
#' by_popCount$Pop1$ref_count

getAlleleCountByPop=function(loci_table, pop_table){
   # sanity check
  stopifnot(is.data.frame(loci_table))
  
  # check that genotypes are valid
  if (!all(loci_table[!is.na(loci_table)] %in% c(0, 1, 2))) {
    stop("Invalid genotype detected in loci_table. Only 0, 1, 2, or NA allowed.")
  }
  samples_2pops <- colnames(loci_table)
  
  # check all samples exist in pop_table
  missing_in_pop_table <- setdiff(samples_2pops, pop_table[,2])
  if (length(missing_in_pop_table) > 0) {
    stop(
      "ERROR: loci_table contains samples not present in pop_table:\n",
      paste(missing_in_pop_table, collapse = ", ")
    )
  }
  
  # get populations present in this subset, sorted for stable order
  pop_level <- sort(unique(pop_table[,1][pop_table[,2] %in% samples_2pops]))
  
  # initiate result list
  res <- list()
  # loop over populations
  for (pop_i in pop_level) {
    # extract sample IDs for this population
    samples_pop_i <- pop_table[pop_table[,1] == pop_i & pop_table[,2] %in% samples_2pops, 2]
    print("get count:")
	print(pop_i)
	print(samples_pop_i)
    # subset loci_table for these samples
    loci_pop_i <- loci_table[, samples_pop_i, drop = FALSE]
    
    # compute allele counts
    res[[pop_i]] <- getAlleleCount(loci_pop_i)
  }
  
  return(res)

}
