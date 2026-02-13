#' Compute Pi on loci
#' @param loci : dataframe format with 0,1,2 encoded genotype, samples in column and ind in lines
#' @param pop_table : data.frame with col1 : pop (char) and col2 : ind (car)
#' @return vec_pi_pop : vector of Pi within and Pi between for the loci
#' @examples : run Pi computation
#' 
#' # input : genotype table, make sure that na are encoded in upper case, R readable. In forme :
#' loci_table = data.frame(I=c("0/1" , "0/1" , "1/1", "0/0" , "0/0" , NA, "0/1" , "0/1" , "1/1"),
#'                       J=c("0/1" , "0/1" , "0/1", NA , "0/0" , NA, "0/1" , "0/1", NA ),
#'                       K=c("0/1" , "1/1" , "1/1", "0/0" , "0/0" , "1/1", "1/1" , "1/1", "0/0" ))
#'
#' # vec_pop : list of pop in order of genotype table with the following encoding : 
#' pop_vec=c("A","A","A", "B", "B", "B", "C", "C", "C")
#' samples=c("s1","s2","s3", "s4", "s5", "s6", "s7", "s8", "s9")
#' pop_table_test = data.frame(cbind(pop_vec, samples))
#'
#' loci_table_T = t(loci_table)
#' colnames(loci_table_T) = samples
#'
#' # add info on GT table in 0,1,2 genotype codes
#' loci_table_T_CV = Convert_GT(loci_table_T)
#' 
#' # compute Pi :
#' getPiWithinBetween(loci = loci_test, pop_table = pop_table_test) 

#' [1] 0.3333333 0.6666667 
#' attr(,"levels") 
#' [1] "Pi_within" "Pi_between"

getPiWithinBetween<-function(loci_pairs, pop_table){
  # TODO : Check for each loci to have at least one SNPs # 
  # get pop names 
  samples_2pops = colnames(loci_pairs)
  pop_level = unique(pop_table[,1][pop_table[,2] %in% samples_2pops])
  
  stopifnot(length(pop_level) == 2)

  # put a warnings if pop_level > 2, because it will not trigger an error
  
  # compute alt and ref freq in the two pops 
  table_freq=getAlleleFreqByPop(loci_pairs, pop_table)

  # compute Pi within : TODO ### NA managment ?? 2pq in popA * 2pq in popB / 2 --> expected het
  Pi_within = 0.5 * (2 * table_freq[[1]]$alt_freq*(1-table_freq[[1]]$alt_freq) + (2 * table_freq[[2]]$alt_freq*(1-table_freq[[2]]$alt_freq)))

  # compute Pi between :
  Pi_between = table_freq[[1]]$alt_freq*(1-table_freq[[2]]$alt_freq) + table_freq[[2]]$alt_freq*(1-table_freq[[1]]$alt_freq)

  # output vector
  vec_pi_pop = list(Pi_within, Pi_between)
  names(vec_pi_pop) = c("Pi_within","Pi_between")
  
  return(vec_pi_pop)
}

#' Compute Hudson Fst per SNP
#'
#' @param loci_table_T Transposed genotype table with 0,1,2 encoding
#' @param pop_table Data.frame with pop_vec and sample names
#' @param Na_rate = 0.2
#' @param MAF_threshold = 0.05
#' @return List of Fst values per SNP and population pair
#' 
#' @examples : run Hudson Fst 
#' 
#' # input : genotype table, make sure that na are encoded in upper case, R readable. In forme :
#' loci_table = data.frame(I=c("0/1" , "0/1" , "1/1", "0/0" , "0/0" , NA, "0/1" , "0/1" , "1/1"),
#'                       J=c("0/1" , "0/1" , "0/1", NA , "0/0" , NA, "0/1" , "0/1", NA ),
#'                       K=c("0/1" , "1/1" , "1/1", "0/0" , "0/0" , "1/1", "1/1" , "1/1", "0/0" ))
#'
#' # vec_pop : list of pop in order of genotype table with the following encoding : 
#' pop_vec=c("A","A","A", "B", "B", "B", "C", "C", "C")
#' samples=c("s1","s2","s3", "s4", "s5", "s6", "s7", "s8", "s9")
#' pop_table_test = data.frame(cbind(pop_vec, samples))
#'
#' loci_table_T = t(loci_table)
#' colnames(loci_table_T) = samples
#'
#' # add info on GT table in 0,1,2 genotype codes
#' loci_table_T_CV = Convert_GT(loci_table_T)
#'
#' # get hudson Fst by pop pairs written in pop_table_test
#' getFstBySNP(loci_table_T = loci_table_T_CV, pop_table = pop_table_test)
#' # test on a subset of table :
#' pop_table_test_sub = pop_table_test[c(1:6),]
#' getFstBySNP(loci_table_T = loci_table_T_CV, pop_table = pop_table_test_sub)
#' 
getFstBySNP<-function(loci_table_T, pop_table,Na_rate = 0.2, MAF_threshold = 0.05){
  Fst_Pairs=list()
  # get the rigth pop in pop table
  samples_2pops = colnames(loci_table_T)
  pop_level = unique(pop_table[,1][pop_table[,2] %in% samples_2pops])
  print(pop_level)
  # Get the pairwise pop couple
  pop_tbl = getPairwisePop(pop_level)
  print(pop_tbl)
  # but pair_name is not define 
  
  for (pair_ni in seq_len(nrow(pop_tbl))) {
    # get samples from the two pop : # rather do a more general function : get sample from pop ID and run it twice
    Sample_pop = getSamplePairPop(Pair_i = as.vector(unlist(pop_tbl[c(pair_ni), c(1,2)])), pop_table)
    print(Sample_pop)
    # get  table of genotype of the two pop
    sample_vec = c(Sample_pop[[1]], Sample_pop[[2]])
    loci_pair <- loci_table_T[, sample_vec, drop = FALSE]
	
	# intialize the Fst by snp vector to the corresponding pairs
	pair_name = paste(pop_tbl[pair_ni, 1],
                      pop_tbl[pair_ni, 2],
                      sep = "_")
    Fst_Pairs[[pair_name]] <- numeric() # carefull if the size is the nrow , and fst are lower than nrow the value by default is 0

    # Sanity check :
    #---------------------#
    # NA rate check
	na_pop1 = rowMeans(is.na(loci_pair[, Sample_pop[[1]]])) > Na_rate
	na_pop2 = rowMeans(is.na(loci_pair[, Sample_pop[[2]]])) > Na_rate

	loci_pair_na = loci_pair[!(na_pop1 | na_pop2), ,drop = F]     
    
	# monomorphic in both populations
	rs <- rowSums(loci_pair_na)
	loci_pair_na_bi = loci_pair_na[!(rs == 0 | rs == 2 * ncol(loci_pair_na)), ,drop = F]

    # MAF 
    freqs <- getAlleleFreq(loci_pair_na_bi)
    maf <- pmin(freqs$alt_freq, freqs$ref_freq)
	loci_pair_na_bi_maf = loci_pair_na_bi[maf > MAF_threshold, ,drop = F]
	
    # compute Pis : 
    Pi_W_B = getPiWithinBetween(loci_pair_na_bi_maf, pop_table) ### warning here
    Fst_SNP = (Pi_W_B$Pi_between - Pi_W_B$Pi_within) / Pi_W_B$Pi_between

    # return wrong list, have to contains all the Fst_SNP
    Fst_Pairs[[pair_name]] = Fst_SNP
}
  return(Fst_Pairs)
}

#' Compute pairwise Hudson Fst on whole loci
#'
#' @param loci_table_T Transposed genotype table with 0,1,2 encoding
#' @param pop_table Data.frame with pop_vec and sample names
#' @param Na_rate = 0.2
#' @param MAF_threshold = 0.05
#' @return List of Fst values and population pair
#' 
#' @examples
#'  : run Global Hudson Fst 
#' 
#' # input : genotype table, make sure that na are encoded in upper case, R readable. In forme :
#' loci_table = data.frame(I=c("0/1" , "0/1" , "1/1", "0/0" , "0/0" , NA, "0/1" , "0/1" , "1/1"),
#'                       J=c("0/1" , "0/1" , "0/1", NA , "0/0" , NA, "0/1" , "0/1", NA ),
#'                       K=c("0/1" , "1/1" , "1/1", "0/0" , "0/0" , "1/1", "1/1" , "1/1", "0/0" ))
#'
#' # vec_pop : list of pop in order of genotype table with the following encoding : 
#' pop_vec=c("A","A","A", "B", "B", "B", "C", "C", "C")
#' samples=c("s1","s2","s3", "s4", "s5", "s6", "s7", "s8", "s9")
#' pop_table_test = data.frame(cbind(pop_vec, samples))
#'
#' loci_table_T = t(loci_table)
#' colnames(loci_table_T) = samples
#'
#' # add info on GT table in 0,1,2 genotype codes
#' loci_table_T_CV = Convert_GT(loci_table_T)
#'
#' # get hudson Fst by pop pairs written in pop_table_test
#' getGlobalFst(loci_table_T = loci_table_T_CV, pop_table = pop_table_test)
#'> Fst_Global
#'$A_B
#'[1] 0.266881
#'$B_C
#'[1] 0.3564457
#'
#'$A_C
#'[1] 0.3770475

getGlobalFst <- function(loci_table_T, pop_table, Na_rate = 0.2, MAF_threshold = 0.05) {
  
  # get the rigth pop in pop table
  samples_2pops = colnames(loci_table_T)
  pop_level = unique(pop_table[,1][pop_table[,2] %in% samples_2pops])
  print(pop_level)
  # Get the pairwise pop couple
  pop_tbl = getPairwisePop(pop_level)
  
  FstGlobal_Pairs <- list()
  
  for (pair_ni in seq_len(nrow(pop_tbl))) {
    
    # populations in the pair
    pair_pops <- as.character(pop_tbl[pair_ni, ])
    print(pair_pops)
    # samples for the two populations
    Sample_pop <- getSamplePairPop(Pair_i = pair_pops, pop_table = pop_table)
    
    sum_pi_within  <- 0
    sum_pi_between <- 0
    tot_loci = 0
 
    # add a function to convert locus in table 
    loci_pair <- loci_table_T[, c(Sample_pop[[1]], Sample_pop[[2]]), drop = FALSE]
    loci_pair <- as.data.frame(loci_pair)
    
    # Sanity check :
    #---------------------#
    # NA rate check
	na_pop1 = rowMeans(is.na(loci_pair[, Sample_pop[[1]]])) > Na_rate
	na_pop2 = rowMeans(is.na(loci_pair[, Sample_pop[[2]]])) > Na_rate

	loci_pair_na = loci_pair[!(na_pop1 | na_pop2), ,drop = F]     
    
	# monomorphic in both populations
	rs <- rowSums(loci_pair_na)
	loci_pair_na_bi = loci_pair_na[!(rs == 0 | rs == 2 * ncol(loci_pair_na)), ,drop = F]

    # MAF 
    freqs <- getAlleleFreq(loci_pair_na_bi)
    maf <- pmin(freqs$alt_freq, freqs$ref_freq)
	loci_pair_na_bi_maf = loci_pair_na_bi[maf > MAF_threshold, ,drop = F]
	
	
    print("nb of loci : ")
    print(dim(loci_pair_na_bi_maf)[1])
    # Hudson global Fst
	
	# compute Pis : 
    Pi_W_B = getPiWithinBetween(loci_pair_na_bi_maf, pop_table) ### warning here
    
	Fst_SNP=sum(Pi_W_B$Pi_between - Pi_W_B$Pi_within, na.rm = T) / sum(Pi_W_B$Pi_between, na.rm=T)
	print(Fst_SNP)
    pair_name=paste(pair_pops[1], pair_pops[2], sep = "_")
    FstGlobal_Pairs[[pair_name]]=Fst_SNP
  }
  return(FstGlobal_Pairs)
}
