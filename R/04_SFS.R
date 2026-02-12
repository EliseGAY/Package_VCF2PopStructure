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
#'   minor_count = min(count_ref,count_snp)
#' 
#' loci_table_Test=loci_table_convert[c(1:50000),]
#' dim(loci_table_Test)
#' SFS_folded=Get_FoldedSFS(loci_table_Test)
#' [1] "nb_class:"
#' [1] 31
#' unlist(SFS_folded)
#' [1] 6181 3751 2238 1773 1373 1034  789  765  674  537  444  399  393  350  298  298  285  265  246  211  212  214  179  192  181  171  167  138  149  152  109

Get_FoldedSFS = function(loci_table_T){
  # set the classes
  nb_samples=dim(loci_table_T)[2]
  count_class=seq_len(nb_samples)
  print("nb_class:")
  print(length(count_class))
  
  list_class=list()  
  for(i in c(seq_len(nb_samples))){
    list_class[i]=0
  }
  # Loop over loci to compute frequencies
  for(i in seq_len(nrow(loci_table_T))){
    current_pos = loci_table_T[i,]
    # remove loci if NA or non-variant loci
    if(NA %in% current_pos || sum(current_pos, na.rm = T) == 0){next}
    count_snp=as.numeric(sum(current_pos))
    count_ref=as.numeric(nb_samples*2 - count_snp)
    minor_count = min(count_ref,count_snp)
    if (count_snp %in% count_class){
      list_class[minor_count] = as.numeric(list_class[minor_count]) + 1
    }
  }
  return(list_class)
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
#'    count_snp_i=as.numeric(sum(current_pop_i))
#'    count_ref_i=as.numeric(size_pop_i*2 - count_snp_i)
#'    # get minor allele count
#'    minor_count_i = min(count_ref_i,count_snp_i)
#'    # get minor allele freq
#'    minor_freq_i = minor_count_i / size_pop_i
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
#' [1]  8 11
#' 
#' > SFS_2D
#' c(0, 0.1) c(0.1, 0.2) c(0.2, 0.3) c(0.3, 0.4) c(0.4, 0.5) c(0.5, 0.6) c(0.6, 0.7) c(0.7, 0.8) c(0.8, 0.9) c(0.9, 1)
#' c(0, 0.1)        4742        1567        1150         826         825         627         658         585         535       585
#' c(0.1, 0.2)      7010         350         307         245         259         257         277         263         257       218
#' c(0.2, 0.3)      4168         260         213         202         240         250         213         239         216       202
#' c(0.3, 0.4)      3110         223         219         170         198         186         190         194         202       204
#' c(0.4, 0.5)         0           0           0           0           0           0           0           0           0         0
#' c(0.5, 0.6)      2739         228         226         203         195         198         191         192         186       182
#' c(0.6, 0.7)      2367         205         201         157         160         186         184         191         164       166
#' c(0.7, 0.8)      2371         215         218         162         208         199         210         153         148       187
#' c(0.8, 0.9)      1923         225         224         167         185         176         166         165         171       163
#' c(0.9, 1)           0           0           0           0           0           0           0           0           0         0

Get2DSFS = function(pop_table, pop_pair, loci_table_T){
  
  # samples for the two populations
  samples_names = colnames(loci_table_T)
  
  # create a list of samples and their pop names
  meta_sub = pop_table[pop_table[,2] %in% samples_names,]
  samples_list=split(meta_sub$GT_sample, meta_sub$Population)
  # Get the user-chosen pop pair
  Sample_2pop <- samples_list[pop_pairs]
  print("pop sizes")
  print(c(length(Sample_2pop[[1]]), length(Sample_2pop[[2]])))
  # initialization of the 2D-SFS :
  n_class=10
  interval = to_list(for(i in seq(1:n_class)) seq(0,1,by=0.1)[c(i, i+1)])
  SFS_2D <- matrix(0, n_class, n_class)
  colnames(SFS_2D)=as.character(interval)
  rownames(SFS_2D)=as.character(interval)
  
  for (i in seq_len(nrow(loci_table_T))) {
    # subset the genotype for each locus with the 2 pops to compute 2D sfs
    current_pos <- loci_table_T[i,c(Sample_2pop[[1]], Sample_2pop[[2]]), drop = FALSE]
    
    if(NA %in% current_pos || sum(current_pos, na.rm = T) == 0){next}
    class_loci_i=list()
    
    # Loop over the 2 pop to store the current minor allale frequencies
    for(pop in c(1,2)){
      # get only one pop 
      current_pop_i = current_pos[,Sample_2pop[[pop]]]
      
      # samples size in pop i
      size_pop_i=length(Sample_2pop[[pop]])
      # count of minor and major allele
      count_snp_i=as.numeric(sum(current_pop_i))
      count_ref_i=as.numeric(size_pop_i*2 - count_snp_i)
      # get minor allele count
      minor_count_i = min(count_ref_i,count_snp_i)
      # get minor allele freq
      minor_freq_i = minor_count_i / size_pop_i
      # assign class to freq
      class_i=to_vec(for(i in 1:length(interval)) if(minor_freq_i >= as.numeric(unlist(interval[i])[1]) && minor_freq_i < as.numeric(unlist(interval[i])[2])) i)
      class_loci_i[pop]=class_i
    }
    
    SFS_2D[unlist(class_loci_i)[1], unlist(class_loci_i)[2]] = SFS_2D[unlist(class_loci_i)[1], unlist(class_loci_i)[2]] + 1
    
  }
  return(SFS_2D)
}
