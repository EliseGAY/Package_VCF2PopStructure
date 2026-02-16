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

getPiWithinBetween_Freq<-function(loci_pairs, pop_table){
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




################### TODO #####################################
getPiWithinBetween_Count<-function(loci_pairs, pop_table){

  # Samples present in the locus table
  samples_2pops <- colnames(loci_pairs)

  # Subset metadata
  pop_sub <- pop_table[pop_table[,2] %in% samples_2pops, , drop = FALSE]

  # Get the two populations in a stable order
  pop_level <- sort(unique(pop_sub[,1]))
  stopifnot(length(pop_level) == 2)

  # Build ordered population list (NO split!)
  pop_list <- lapply(pop_level, function(p){
    pop_sub[pop_sub[,1] == p, 2]
  })
  names(pop_list) <- pop_level
  
  # compute alt and ref freq in the two pops 
  table_counts=getAlleleCountByPop(loci_pairs, pop_table)
  Pi_within=list()
  nb_chr_pop=list()
  
  for(i in c(1,2)){
	  # set parameters for Pi :
	  nb_chr = 2*length(pop_list[[i]])
	  n_pairs = (nb_chr * (nb_chr - 1)) / 2
	  n_pairs_ref = ((nb_chr - table_counts[[i]]$alt_count) * (nb_chr - table_counts[[i]]$alt_count - 1)) / 2
	  n_paris_alt = (table_counts[[i]]$alt_count * (table_counts[[i]]$alt_count - 1)) / 2
	  npairs_diff = n_pairs - (n_pairs_ref + n_paris_alt)
	  
	  # mean pairwise difference whithin pop
	  Pi_within[[i]] = npairs_diff / n_pairs
	  nb_chr_pop[[i]] = nb_chr
	  }
	  
  # carefull about the nb of Pi within in each pop.
  Pi_within_mean =(nb_chr_pop[[1]] * Pi_within[[1]] + nb_chr_pop[[2]] * Pi_within[[2]]) / (nb_chr_pop[[1]] + nb_chr_pop[[2]])
  # compute Pi between :
  # ALT–ALT + REF–REF
  cross_pairs = nb_chr_pop[[1]] * nb_chr_pop[[2]] 
  same_pairs = (table_counts[[1]]$alt_count * table_counts[[2]]$alt_count) + (table_counts[[1]]$ref_count * table_counts[[2]]$ref_count)
  # ALT - REF pairs count
  tot_cross_pair_diff = cross_pairs - same_pairs
  Pi_between = tot_cross_pair_diff / cross_pairs 
  
  # output vector
  vec_pi_pop = list(Pi_within_mean, Pi_between)
  names(vec_pi_pop) = c("Pi_within","Pi_between")
  
  return(vec_pi_pop)
}

################### TODO #####################################


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
getFstBySNP_Freq<-function(loci_table_T, pop_table,Na_rate = 0.2, MAF_threshold = 0.05){
  Fst_Pairs=list()
  # get the rigth pop in pop table
  samples_2pops = colnames(loci_table_T)
  pop_level = unique(pop_table[,1][pop_table[,2] %in% samples_2pops])
  # Get the pairwise pop couple
  pop_tbl = getPairwisePop(pop_level)
  # but pair_name is not define 
  
  for (pair_ni in seq_len(nrow(pop_tbl))) {
    # get samples from the two pop : # rather do a more general function : get sample from pop ID and run it twice
    Sample_pop = getSamplePairPop(Pair_i = as.vector(unlist(pop_tbl[c(pair_ni), c(1,2)])), pop_table)
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
    Pi_W_B = getPiWithinBetween_Freq(loci_pair_na_bi_maf, pop_table) ### warning here
    Fst_SNP = (Pi_W_B$Pi_between - Pi_W_B$Pi_within) / Pi_W_B$Pi_between

    # return wrong list, have to contains all the Fst_SNP
    Fst_Pairs[[pair_name]] = Fst_SNP
}
  return(Fst_Pairs)
}
################### TODO AND DEBUG FUNCTION #####################################

getFstBySNP_Count<-function(loci_table_T, pop_table,Na_rate = 0.2, MAF_threshold = 0.05){
  Fst_Pairs=list()
  # get the rigth pop in pop table
  samples_2pops = colnames(loci_table_T)
  pop_level = unique(pop_table[,1][pop_table[,2] %in% samples_2pops])
  # Get the pairwise pop couple
  pop_tbl = getPairwisePop(pop_level)
  # but pair_name is not define 
  
  for (pair_ni in seq_len(nrow(pop_tbl))) {
	print(pair_ni)
	
    # get samples from the two pop : # rather do a more general function : get sample from pop ID and run it twice
    Sample_pop = getSamplePairPop(Pair_i = as.vector(unlist(pop_tbl[c(pair_ni), c(1,2)])), pop_table)
	# get  table of genotype of the two pop
    sample_vec = c(Sample_pop[[1]], Sample_pop[[2]])
    loci_pair <- loci_table_T[, sample_vec, drop = FALSE]
	print(colnames(loci_pair))
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
    Pi_W_B = getPiWithinBetween_Count(loci_pair_na_bi_maf, pop_table)
	# All NA in PiWithin or Between gonna return NA FST
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
  # Get the pairwise pop couple
  pop_tbl = getPairwisePop(pop_level)
  
  FstGlobal_Pairs <- list()
  
  for (pair_ni in seq_len(nrow(pop_tbl))) {
    
    # populations in the pair
    pair_pops <- as.character(pop_tbl[pair_ni, ])
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
    Pi_W_B = getPiWithinBetween_Freq(loci_pair_na_bi_maf, pop_table) ### warning here
    
	Fst_SNP=sum(Pi_W_B$Pi_between - Pi_W_B$Pi_within, na.rm=T) / sum(Pi_W_B$Pi_between, na.rm=T)
	print(Fst_SNP)
    pair_name=paste(pair_pops[1], pair_pops[2], sep = "_")
    FstGlobal_Pairs[[pair_name]]=Fst_SNP
  }
  return(FstGlobal_Pairs)
}

#' Compute pairwise Hudson Fst on ALT count 
#'
#' @param loci_table_T Transposed genotype table with 0,1,2 encoding
#' @param pop_table Data.frame with pop_vec and sample names
#' @param Na_rate = 0.2
#' @param MAF_threshold = 0.05
#' @return List of Fst values and population pair
#' 
#' @examples
#' 
#' `1622W1_S254` `1623W1_S276` `1624W1_S392` `1625Q_S431` `1625W1_S384` `1626W1_S364` `1627W1_S385` `1629W1_S379` `1681W1_S374` `1682W1_S380` `1684W1_S376` `1686W1_S402`
#' 1             0             0             0            0             0             0             0             0             0             0             0             0
#' 2             0             0             0            0             0             0             0             0             0             0             0             0
#' 3             0             0             0            0             0             0             0             0             0             0             0             0
#'# ℹ 19 more variables: `
#'> metadata
#'   Population   GT_sample
#'1    pop1 1551W1_S395
#'2    pop1 1553W1_S314
#'3    pop1 1554W1_S341
#'4    pop1 1555W1_S397
#' 5    pop1 1556W1_S337

#' Fst_GlobalCount = getGlobalFst_Count(loci_table_T = loci_table_T_CV, 
#'                          pop_table = meta_sub, Na_rate = 0.2, MAF_threshold = 0)
#'						  tibble(loci_table_T_CV)
#' [1] "pop1" "pop2"  
#' [1] "get count:"
#' [1] "pop1"
#' [1] "1551W1_S395" "1553W1_S314" "1554W1_S341" "1555W1_S397" "1556W1_S337" "1557Q_S432"  "1557W1_S275" "1558W1_S354" "1559W1_S357" "1560W1_S300" "1562W1_S242" "1563W1_S281"
#' [1] "get count:"
#' [1] "pop2"
#' [1] "1622W1_S254" "1623W1_S276" "1624W1_S392" "1625Q_S431"  "1625W1_S384" "1626W1_S364" "1627W1_S385" "1629W1_S379"
#' [1] "nb of loci : "
#' [1] 54200
#' [1] 0.2316629

getGlobalFst_Count <- function(loci_table_T, pop_table, Na_rate = 0.2, MAF_threshold = 0.05) {
  
  # get the rigth pop in pop table
  samples_2pops = colnames(loci_table_T)
  pop_level = unique(pop_table[,1][pop_table[,2] %in% samples_2pops])
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
	
    #Hudson global Fst
	
	# compute Pis : 
    Pi_W_B = getPiWithinBetween_Count(loci_pair_na_bi_maf, pop_table)
  	# All NA in PiWithin or Between gonna return NA FST.
	Fst_SNP=sum((Pi_W_B$Pi_between - Pi_W_B$Pi_within), na.rm = T) / sum(Pi_W_B$Pi_between, na.rm=T)
	print("nb of loci : ")
    print(dim(loci_pair_na_bi_maf)[1])
	print(Fst_SNP)
    pair_name=paste(pair_pops[1], pair_pops[2], sep = "_")
    FstGlobal_Pairs[[pair_name]]=Fst_SNP
  }
  return(FstGlobal_Pairs)
}
