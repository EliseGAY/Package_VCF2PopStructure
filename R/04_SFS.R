# Make SFS ##########################################
#' Get2DSFS: Compute pairwise Hudson Fst on whole loci
#' @import comprehenr
#' @param loci_table_T Transposed genotype table with 0,1,2 encoding
#' @return vector of counts of minor allele frequencies
#' 
#' @examples :
#' Compute folded SFS as following : 
#'   count_snp=as.numeric(sum(current_pos))
#'   count_ref=as.numeric(nb_samples*2 - count_snp)
#'   minor_count = pmin(count_ref,count_snp)
#' 
#' loci_table_Test=loci_table_convert[c(1:50000),]
#' dim(loci_table_Test)
#' SFS_folded=Get_FoldedSFS(loci_table_Test)
#' [1] "nb_class:"
#' [1] 31
#' SFS_folded
#' [1]     1     2     3     4     5     6     7     8     9    10    11    12    13    14    15    16    17    18    19    20    21    22    23    24    25    26    27    28    29 
#' #     25810 12912  8908  6501  5006  3969  3148  2772  2340  2079  1934  1724  1473  1450  1400  1325  1221  1339  1203  1115  1137  1225  1090  1054  1039  1038  1082   971   962 
#  #    30    31 
#  #    993   443 

Get_FoldedSFS = function(loci_table_T){
  # set the classes
  nb_samples=dim(loci_table_T)[2]
  count_class=seq_len(nb_samples)
  print("nb_class:")
  print(length(count_class))
  
  # remove loci if NA or non-variant loci
  Purged_loci <- loci_table_T[rowSums(is.na(loci_table_T)) == 0 & rowSums(loci_table_T) != 0 & rowSums(loci_table_T) != nb_samples*2,, drop=FALSE]
  
  # to continue to vectorization
  count_snp=as.numeric(rowSums(Purged_loci))
  count_ref=as.numeric(nb_samples*2 - count_snp)
  minor_count = pmin(count_ref,count_snp)
  SFS = table(minor_count)
  
  return(SFS)
}


# Make 2D SFS ##########################################
#' Get2DSFS: Compute pairwise Hudson Fst on whole loci
#' @import comprehenr
#' @param loci_table_T Transposed genotype table with 0,1,2 encoding
#' @param pop_table Data.frame with pop_vec and sample names
#' @param pop_pairs Pair of pop names
#' @return count matrix of 10 class allele frequencies
#' 
#'  @description
#'  Discard all position with NA or monomorphic in both pops
#'  Get the minor allele frequencies by pop as following : 
#'      # count of minor and major allele
#'      count_snp_i=as.numeric(rowSums(current_pop_i))
#'      count_ref_i=as.numeric(size_pop_i*2 - rowSums(current_pop_i))
#'      # get minor allele count
#'      minor_count_i = pmin(count_ref_i,count_snp_i)
#'      # get minor allele freq
#'      minor_freq_i = minor_count_i / size_pop_i
#'      # assign class to freq
#'      class_i_vec <- findInterval(minor_freq_i, breaks)
#'      freq_pop[[pop]] = class_i_vec

#' 
#' @examples :
#' 
#' loci_table_test=loci_table_convert[c(1:50000),, drop = FALSE]
#' 
#' metadata table : 
#'      Population   GT_sample
#' 2     pop1 s1
#' 3     pop1 s2
#' 4     pop2 s5
#' 5     pop2 s6 ....
#'       .......
#'      
#' SFS_2D = Get2DSFS(pop_table = metadata, pop_pair = c("pop1", "pop2"), loci_table_T = loci_table_convert)
#' 
#' [1] "pop sizes"
#' [1]  12 11
#' 
#' > SFS_2D
  
#'          1      2      3      4      5      6      7      8      9     10     11
#'  1  118068   1741   1327    986    977    872    871    850    795    788    388
#'  2   10934    197    156    141    147    165    139    122    144    111     51
#'  3    7515    151    169    118    119    130    142    115    122    130     69
#'  4    5203    169    153    115    127    121    130    116    124    102     44
#'  5    3931    160    134    123    145    150    132    148    107    134     67
#'  6    5432    319    292    245    267    250    238    256    229    214    120
#'  7    2026    155    140    108    128    119    125    116    118    119     49
#'  8    1808    114    139    125    123     91    123    102     97    112     56
#'  9    1623    157    153    110    136    106    102     93     93    130     51
#'  10   1525    129    135     96    130    118    115    103     96     94     73
#'  11    771     66     68     72     67     55     59     46     45     46     40

Get2DSFS = function(pop_table, pop_pair, loci_table_T){
  
  # samples for the two populations
  samples_names = colnames(loci_table_T)
  # create a list of samples and their pop names
  meta_sub = pop_table[pop_table[,2] %in% samples_names,]
  samples_list=split(meta_sub$GT_sample, meta_sub$Population)
  # Get the user-chosen pop pair
  Sample_2pop <- samples_list[pop_pair]
  print("pop sizes")
  print(c(length(Sample_2pop[[1]]), length(Sample_2pop[[2]])))
  
  # initialization of the 2D-SFS :
  n_class=10
  interval = to_list(for(i in seq(1:n_class)) seq(0,1,by=0.1)[c(i, i+1)])
  SFS_2D <- matrix(0, n_class, n_class)
  colnames(SFS_2D)=as.character(interval)
  rownames(SFS_2D)=as.character(interval)
  
  # initalize the class storage
  breaks <- seq(0, 1,length.out = 11)
  freq_pop=list()
  
  #remove na and monomorphic site
  Purged_loci <- loci_table_T[rowSums(is.na(loci_table_T)) == 0 &  rowSums(loci_table_T) != 0, c(Sample_2pop[[1]], Sample_2pop[[2]]), drop=FALSE]
  
    # Loop over the 2 pop to store the current minor allele frequencies
    for(pop in c(1,2)){
      # get only one pop 
      current_pop_i = loci_table_T[,Sample_2pop[[pop]]]
      
      # samples size in pop i
      size_pop_i=length(colnames(current_pop_i))
      # count of minor and major allele
      count_snp_i=as.numeric(rowSums(current_pop_i))
      count_ref_i=as.numeric(size_pop_i*2 - rowSums(current_pop_i))
      # get minor allele count
      minor_count_i = pmin(count_ref_i,count_snp_i)
      # get minor allele freq
      minor_freq_i = minor_count_i / size_pop_i
      # assign class to freq
      class_i_vec <- findInterval(minor_freq_i, breaks)
      freq_pop[[pop]] = class_i_vec
    }
  
    SFS_2D <- table(freq_pop[[1]], freq_pop[[2]])
    return(SFS_2D)
  }

