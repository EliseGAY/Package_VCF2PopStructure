#'GetGenoFreqByPop
#' @param geno_table : dataframe format with 0/1 diploid genotype format or 0 ,1 haplpoid genotype, samples in column and position in lines
#' @param pop_list : list of samples (values) in each pop (names)
#' @return geno_freq : table with  col : "chr", "position", "POP1_het_freq", "POP1__hom_alt_freq", "POP1_hom_ref_freq" 
#' @description
#' Na are removed from the genotype count but not the corresponding sample. The function should be run on table with no missing for more accuracy
#' 
#' @examples : 
#' loci_table = extract.gt(VCFR_data, element = "GT")
#' loci_table = as.data.frame(loci_table)
#' loci_table[c(1:10),c(1:10)]
#'             1622W1_S254 1623W1_S276 1624W1_S392 1625Q_S431 1625W1_S384 1626W1_S364 1627W1_S385 1629W1_S379 1681W1_S374 1682W1_S380
#' ptg000007l_9661          0/0         0/0         0/0        0/0         0/0         0/0         0/0         0/0         0/0         0/0
#' ptg000007l_9668          0/0         0/0         0/0        0/0         0/0         0/0         0/0         0/0         0/0         0/0
#' ptg000007l_9670          0/0         0/0         0/0        0/0         0/0         0/0         0/0         0/0         0/0         0/0
#' ptg000007l_9678          0/0         0/0         0/0        0/0         0/0         0/0         0/0         0/0         1/1         1/1
#' 
#' pop_list = pop_list
#'$POP1
#' [1] "1681W1_S374" "1682W1_S380" "1684W1_S376" "1686W1_S402" 
#' $POP2
#' [1] "1551W1_S395" "1553W1_S314" "1554W1_S341" "1555W1_S397" "1556W1_S337" 
#' $POP3
#' [1] "1622W1_S254" "1623W1_S276" "1624W1_S392" "1625Q_S431"  "1625W1_S384"
#' 
#' table = GetGenoFreqByPop(geno_table = loci_table, chr_name = chr, ploidy = "diploid", pop_list = pop_list)
#' 
#' chr position POP1_het_freq POP1_hom_alt_freq POP1_hom_ref_freq POP1_hom_tot_freq POP1_het_freq ... same with POP2 and POP3
#' 1  ptg000007l     9661         0.0000000             0.0000000             1.0000000             1.0000000        0.00000000
#' 2  ptg000007l     9668         0.0000000             0.0000000             1.0000000             1.0000000        0.00000000
#' 3  ptg000007l     9670         0.0000000             0.0000000             1.0000000             1.0000000        0.08333333
#' 4  ptg000007l     9678         0.2727273             0.4545455             0.2727273             0.7272727        0.00000000
#' 5  ptg000007l     9687         0.0000000             0.0000000             1.0000000             1.0000000        0.00000000

# compute genotype count in all ind
GetGenoFreqByPop = function(geno_table, chr_name=NULL, ploidy=NULL, pop_list){
  
  # check for pop list and genot table match
  if (FALSE %in% c(as.vector(unlist(pop_list))[c(1,2,3)] %in% colnames(geno_table))){
    stop("pop_list has much samples than genotype table")
  }
  # subset the geno table according to pop list
  geno_all=as.data.frame(loci_table[,as.vector(unlist(pop_list))])
  geno_freq=data.frame("chr"=rep(chr_name,dim(geno_all)[1]))
  # add position column
  position_vector=as.numeric(str_remove(rownames(geno_all), paste(chr,"_",sep="")))
  # Add vector of position in the sequencing depth table 
  geno_freq[,"position"]=position_vector
  
  # check and set poidy
  if (ploidy == "haploid") {
    genotype=list("ref" = "0", "alt"= "1")
  } else if (ploidy == "diploid") {
    genotype=list("ref" = "0/0", "het"="0/1", alt = "1/1")
  } else {
    stop("ERROR: ploidy values are haploid or diploid with 0/ and /1 units \n")}
  
  for (i in seq(1:length(pop_list))){
    i_name = names(pop_list[i])
    # get temporary table corresponding to pop i
    geno_temp=as.data.frame(geno_all[,as.vector(pop_list[[i]])])
    
    # set colnames with pop
    het = paste(i_name, "het_freq", sep = "_")
    hom_alt = paste(i_name, "hom_alt_freq", sep = "_")
    hom_ref = paste(i_name, "hom_ref_freq", sep = "_")
    hom_all = paste(i_name, "hom_tot_freq", sep = "_")
    
    # compute all frequencies
    geno_freq[,het] = rowSums(geno_temp==genotype$het, na.rm = T)/dim(geno_temp)[2]
    geno_freq[,hom_alt] = rowSums(geno_temp==genotype$alt, na.rm = T)/dim(geno_temp)[2]
    geno_freq[,hom_ref]  = rowSums(geno_temp==genotype$ref, na.rm = T)/dim(geno_temp)[2]
    geno_freq[,hom_all]  = (rowSums(geno_temp == genotype$ref, na.rm = T) +  rowSums(geno_temp == genotype$alt, na.rm = T))/dim(geno_temp)[2]
  }
  return(geno_freq) 
}
