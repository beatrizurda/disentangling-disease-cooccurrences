#####################################################################################
# COMPARING WITH THE GENOME MEDICINE EPIDEMIOLOGICAL NETWORK
#
#  
# 
# 
# Beatriz Urda García 2025, Aug-Sept
######################################################################################


setwd('~/Desktop/ANALYSIS/Genome_medicine/')
source('ukb_comparison_tools.R')

# To test the script
network_type = "new_dsn"
split_duplicates = FALSE    # This is the default option
set.seed(1)


compute_overlap_Dong_et_al <- function(network_type,split_duplicates, distance="spearman", alpha = 0.05){

  # Network from the Genome Medicine Dong et al.
  dong_epidem =  read.csv2('ukb_multimorbidity.csv',stringsAsFactors = F,sep=",",header=T)
  colnames(dong_epidem)[1:2] = c('Dis1', 'Dis2')

  # OPEN THE NETWORK
  if(network_type == "dsn"){
    # ICD9-level DSN
    dsn = read.csv2('../Network_building/Defined_networks/Final_networks_RNAseq_paper/icd9_pairwise_union_spearman_distance_sDEGs_pos_network.txt',
                    stringsAsFactors = F,sep="\t",header=T)
    dim(dsn) # 347

    # Transforming DSN into ICD10 codes
    dsn = network_from_icd9_to_icd10(dsn, icd_dic); dim(dsn) # 417
    length(unique(union(dsn$Dis1, dsn$Dis2))) # 41
    
    distance=""

  }else if(network_type == "ssn"){
    # SSN
    ssn_filename = "../Network_building/Defined_networks/metapatients_and_disease/metap_dis_pairwise_union_spearman_distance_sDEGs_pos_network.txt"
    
    # Selecting D-M network with disease names (BreastCancer_1 --> BreastCancer)
    ssn = get_dis_representation_of_ssn(ssn_filename); nrow(ssn)
    
    # Transforming SSN into ICD10 codes
    ssn = network_from_dis_to_icd10(ssn); nrow(ssn)
    dsn = ssn
    
    distance=""
    
  }else if(network_type == "new_dsn"){
    #  DSN
    if(distance == "spearman"){
      netfilename = '../../GEVariability/Network_building/Defined_networks/all_diseases/DEGs/pairwise_union_spearman_distance_sDEGs_pos_network.txt'
    }else if(distance == "cosine"){
      netfilename = '../../GEVariability/Network_building/Defined_networks/all_diseases/DEGs/pairwise_union_cosine_distance_sDEGs_pos_network.txt'
    }else if(distance == "fisher"){
      netfilename = '../../GEVariability/Network_building/Defined_networks/all_diseases/DEGs/pairwise_union_spearman_cosine_distance_sDEGs_merged_pos_network.txt'
    }
    dsn = read.csv2(netfilename,
                    stringsAsFactors = F,sep="\t",header=T)
    nrow(dsn)
    
    if(alpha != 0.05){
      dsn$adj_pvalue = as.numeric(as.character(dsn$adj_pvalue))
      dsn = dsn[dsn$adj_pvalue <= alpha, ]
    }
    nrow(dsn)
    str(dsn)

    # Transform into icd10 codes
    dsn = network_from_dis_to_icd10(dsn)

    # Split obesity_t2b interactions in two
    # dsn = split_obesity_t2b_in_network(dsn); nrow(dsn)
    # length(unique(union(dsn$Dis1, dsn$Dis2))) # 88
    
    # Split duplicated ICD10 codes in two
    dsn = split_duplicates_in_dong_et_al(dsn)
    
  }else if(network_type == "new_ssn"){
    # SSN
    if(distance == "spearman"){
      netfilename = '../../GEVariability/Network_building/Defined_networks/all_diseases/metapatients_and_disease/metap_dis_pairwise_union_spearman_distance_sDEGs_pos_network.txt'
    }else if(distance == "cosine"){
      netfilename = '../../GEVariability/Network_building/Defined_networks/all_diseases/metapatients_and_disease/metap_dis_pairwise_union_cosine_distance_sDEGs_pos_network.txt'
    }else if(distance == "fisher"){
      netfilename = '../../GEVariability/Network_building/Defined_networks/all_diseases/metapatients_and_disease/metap_dis_pairwise_union_spearman_cosine_distance_sDEGs_merged_pos_network.txt'
    }
    
    # Selecting D-M network with disease names (BreastCancer_1 --> BreastCancer)
    ssn = get_dis_representation_of_ssn(netfilename); nrow(ssn)
    
    if(alpha != 0.05){
      ssn$adj_pvalue = as.numeric(as.character(ssn$adj_pvalue))
      ssn = ssn[ssn$adj_pvalue <= alpha, ]
    }
    
    # Transforming SSN into ICD10 codes
    ssn = network_from_dis_to_icd10(ssn); nrow(ssn)
    dsn = ssn
    
    # Split duplicated ICD10 codes in two
    dsn = split_duplicates_in_dong_et_al(dsn)
    
  }


  # Handling ICD10 synonyms in Dong et al. i.e.,Dong interactions with multiple ICD10 codes in different lines
  if(split_duplicates == TRUE){    # Split Dong interactions with multiple ICD10 codes in different lines
    dong_epidem = split_duplicates_in_dong_et_al(dong_epidem)
  }else{                           # Add synonyms to the DSN network
    syn_groups = identity_synonyms_in_Dong(dong_epidem)
    syn_dis = get_diseases_with_synonyms(syn_groups); length(syn_dis) # 54

    # Put diseases synonyms in the DSN
    original_dsn = dsn
    dsn = add_dis_synonyms_to_network(dsn, syn_dis, syn_groups)
  }

  original_dim_dsn = nrow(dsn)

  # Sort the interactions
  dsn = sort_icd_interactions(dsn)
  dong_epidem = sort_icd_interactions(dong_epidem)

  # Remove duplicated interactions
  dsn = remove_duplicated_interactions(dsn) # 354 unique interactions
  dong_epidem = remove_duplicated_interactions(dong_epidem) # 11285 - No duplicated interactions

  # Remove interactions where icd Dis1 == icd Dis2
  dsn = dsn[dsn$Dis1 != dsn$Dis2, ]; dim(dsn) # 353 unique interactions
  dong_epidem = dong_epidem[dong_epidem$Dis1 != dong_epidem$Dis2, ]; dim(dong_epidem) # 11285 - the same as before

  unique_int_dsn = nrow(dsn)
  unique_int_dong = nrow(dong_epidem)

  # Common ICD10s
  icds_dsn = unique(union(dsn$Dis1, dsn$Dis2)); length(icds_dsn) # 41
  icds_dong = unique(union(dong_epidem$Dis1, dong_epidem$Dis2)); length(icds_dong) # 467 --> 438

  common_icds = intersect(icds_dsn, icds_dong); length(common_icds) # 30

  # Keep networks entailing common icds
  dsn = dsn[(dsn$Dis1 %in% common_icds) & (dsn$Dis2 %in% common_icds), ]; dim(dsn)  # 205
  dong_epidem = dong_epidem[(dong_epidem$Dis1 %in% common_icds) & (dong_epidem$Dis2 %in% common_icds), ]; dim(dong_epidem) # 117
  length(unique(union(dong_epidem$Dis1, dong_epidem$Dis2))) # 28 -- "C71" "C73" are the ones wo links
  length(unique(union(dsn$Dis1, dsn$Dis2))) # 30

  # Compute the overlap
  dsn_graph <- graph_from_data_frame(dsn, directed=FALSE)
  dong_epidem_graph <- graph_from_data_frame(dong_epidem, directed=FALSE)
  intersect <- intersection(dsn_graph,dong_epidem_graph)
  overlap = gsize(intersect) # Number of edges # 55

  precision = (overlap/(dim(dsn)[1]))*100; precision      # 26.82927 PRECISION
  recall = (overlap/(dim(dong_epidem)[1]))*100; recall    # 47.00855 RECALL


  ## Generate a network with all my possible interactions between common diseases
  N <- length(common_icds) # 30

  connected_net <- data.frame(Dis1=c(),Dis2=c(),Distance=c())
  set.seed(5)
  for (k1 in 1:(N-1)){
    for (k2 in (k1+1):N){
      if(is.na(common_icds[k2])== T){
        print(k1)
        print(k2)
        break
      }
      connected_net <- rbind(connected_net,data.frame(common_icds[k1],common_icds[k2],1))
    }
  }

  ## Sort the connected network
  colnames(connected_net) <- c('Dis1','Dis2','Distance')
  connected_net$Dis1 = as.character(connected_net$Dis1); connected_net$Dis2 = as.character(connected_net$Dis2)
  connected_net <- sort_icd_interactions(connected_net)
  dim(connected_net) # ((30*30) - 30)/2
  # 435 possible interactions

  ## Generate 10.000 iterations of: select random network with my nº of interactions and compute overlap with Dong et al.
  ## Shuffling labels and keeping the degree distribution of the graph
  n_iters <- 10000
  random_overlaps <- c()
  random_overlaps_second <- c()
  set.seed(5) ### ESTABLECER SEMILLA

  for(k in 1:n_iters){
    # p-value for the comparison with Dong et al. (RECALL)
    cgraph <- rewire(dsn_graph, with=keeping_degseq(loops=FALSE, niter = ecount(dsn_graph)*10))
    # str(cgraph)
    cintersect <- intersection(cgraph,dong_epidem_graph)
    # cintersect
    gsize(cintersect)
    random_overlaps <- append(random_overlaps,gsize(cintersect))
    # print_all(rewire(mygraph, with = keeping_degseq(loops=FALSE, niter = vcount(mygraph) * 10)))
  }

  for(k in 1:n_iters){ # it has to go into another loop to keep consistent p-values (reproducibility)
    # p-value for the comparison with the DSN (PRECISION)
    cgraph_second <- rewire(dong_epidem_graph, with=keeping_degseq(loops=FALSE, niter = ecount(dong_epidem_graph)*10))
    # str(cgraph)
    cintersect_second <- intersection(cgraph_second,dsn_graph)
    # cintersect
    gsize(cintersect_second)
    random_overlaps_second <- append(random_overlaps_second,gsize(cintersect_second))
    # print_all(rewire(mygraph, with = keeping_degseq(loops=FALSE, niter = vcount(mygraph) * 10)))
  }

  # Comparison with Dong et al. --> Recall
  pval <- length(which(random_overlaps>=gsize(intersect)))/length(random_overlaps)
  pval
  # 0.326

  # Comparison with DSN --> Precision
  pval_second <- length(which(random_overlaps_second>=gsize(intersect)))/length(random_overlaps_second)
  pval_second
  # 0.6076

  outdf = c(net, split_duplicates, original_dim_dsn, unique_int_dsn, unique_int_dong, length(icds_dsn), length(icds_dong), length(common_icds),
            nrow(dsn), nrow(dong_epidem), overlap, recall, precision, pval, pval_second, distance)
  names(outdf) = c("net","split_duplicates","original_dim_net", "unique_int_net", "unique_int_dong", "icds_net", "icds_dong", "common_icds",
                   "unique_int_net_common_icds", "unique_int_dong_common_icds", "overlap","recall", "precision","pval_recall", "pval_second_precision", "distance")
  return(outdf)
}

# dsn_output = compute_overlap_Dong_et_al(network_type = "dsn", split_duplicates=FALSE)


# To test the script
network_type = "dsn"
split_duplicates = FALSE

output = c()
for(option in c(FALSE,TRUE)){
  for(net in c("dsn","ssn")){
    print(option); print(net)
    coutput = compute_overlap_Dong_et_al(network_type = net, split_duplicates = option)
    output = rbind(output, coutput)
  }
}

output = c()
for(option in c(FALSE,TRUE)){
  for(net in c("new_dsn","new_ssn")){
    net = "new_ssn"
    for(distance in c("spearman", "cosine", "fisher")){
      print(option); print(net); print(distance)
      coutput = compute_overlap_Dong_et_al(network_type = net, split_duplicates = option, distance)
      output = rbind(output, coutput)
    }
  }
}

output = c()
for(option in c(FALSE,TRUE)){
  for(net in c("new_dsn","new_ssn")){
    for(distance in c("spearman", "cosine", "fisher")){
      print(option); print(net); print(distance)
      coutput = compute_overlap_Dong_et_al(network_type = net, split_duplicates = option, distance)
      output = rbind(output, coutput)
    }
  }
}
output = as.data.frame(output)
output1 = output

# RUN alpha 0.1:
# for(net in c("new_dsn","new_ssn")){
#   for(distance in c("spearman", "cosine", "fisher")){
#     # coutput = compute_overlap_Dong_et_al(network_type = "new_dsn", split_duplicates = FALSE, "spearman", alpha=0.01)
#     coutput = compute_overlap_Dong_et_al(network_type = net, split_duplicates = FALSE, distance, alpha=0.01)
#     print(coutput)
#     output = rbind(output, coutput)
#   }
# }


write.table(output, file="Overlaps_with_Dong_epidemiology_small_alpha.txt", sep="\t", row.names = FALSE, col.names = TRUE, quote = FALSE)



#### RUN A GIVEN SET OF PARAMETERS
# Network from the Genome Medicine Dong et al. 
dong_epidem =  read.csv2('ukb_multimorbidity.csv',stringsAsFactors = F,sep=",",header=T)
colnames(dong_epidem)[1:2] = c('Dis1', 'Dis2')


# OPEN THE NETWORK
if(network_type == "dsn"){
  # ICD9-level DSN 
  dsn = read.csv2('../Network_building/Defined_networks/Final_networks_RNAseq_paper/icd9_pairwise_union_spearman_distance_sDEGs_pos_network.txt',
                  stringsAsFactors = F,sep="\t",header=T)
  dim(dsn) # 347
  
  # Transforming DSN into ICD10 codes
  dsn = network_from_icd9_to_icd10(dsn, icd_dic); dim(dsn) # 417
  length(unique(union(dsn$Dis1, dsn$Dis2))) # 41

}else if(network_type == "new_dsn"){
  # ICD9-level DSN 
  dsn = read.csv2('../../GEVariability/Network_building/Defined_networks/all_diseases/DEGs/pairwise_union_spearman_distance_sDEGs_pos_network.txt',
                  stringsAsFactors = F,sep="\t",header=T)
  nrow(dsn)
  
  # Transform into icd10 codes
  dsn = network_from_dis_to_icd10(dsn)
  
  # Split obesity_t2b interactions in two
  dsn = split_obesity_t2b_in_network(dsn); nrow(dsn)
  length(unique(union(dsn$Dis1, dsn$Dis2))) # 88
  
}else{
  # SSN
  ssn_filename = "../Network_building/Defined_networks/metapatients_and_disease/metap_dis_pairwise_union_spearman_distance_sDEGs_pos_network.txt"
  
  # Selecting D-M network with disease names (BreastCancer_1 --> BreastCancer)
  ssn = get_dis_representation_of_ssn(ssn_filename); nrow(ssn)
  
  # Transforming SSN into ICD10 codes
  ssn = network_from_dis_to_icd10(ssn); nrow(ssn)
  dsn = ssn
}


# Handling ICD10 synonyms in Dong et al. i.e.,Dong interactions with multiple ICD10 codes in different lines 
if(split_duplicates == TRUE){    # Split Dong interactions with multiple ICD10 codes in different lines
  dong_epidem = split_duplicates_in_dong_et_al(dong_epidem)
}else{                           # Add synonyms to the DSN network
  syn_groups = identity_synonyms_in_Dong(dong_epidem)
  syn_dis = get_diseases_with_synonyms(syn_groups); length(syn_dis) # 54
  
  # Put diseases synonyms in the DSN
  original_dsn = dsn
  dsn = add_dis_synonyms_to_network(dsn, syn_dis, syn_groups)
}

# Sort the interactions
dsn = sort_icd_interactions(dsn)
dong_epidem = sort_icd_interactions(dong_epidem)

# Remove duplicated interactions
dsn = remove_duplicated_interactions(dsn) # 354 unique interactions
dong_epidem = remove_duplicated_interactions(dong_epidem) # 11285 - No duplicated interactions

# Remove interactions where icd Dis1 == icd Dis2
dsn = dsn[dsn$Dis1 != dsn$Dis2, ]; dim(dsn) # 353 unique interactions
dong_epidem = dong_epidem[dong_epidem$Dis1 != dong_epidem$Dis2, ]; dim(dong_epidem) # 11285 - the same as before

# Common ICD10s
icds_dsn = unique(union(dsn$Dis1, dsn$Dis2)); length(icds_dsn) # 41
icds_dong = unique(union(dong_epidem$Dis1, dong_epidem$Dis2)); length(icds_dong) # 467 --> 438

common_icds = intersect(icds_dsn, icds_dong); length(common_icds) # 30

# Keep networks entailing common icds
dsn = dsn[(dsn$Dis1 %in% common_icds) & (dsn$Dis2 %in% common_icds), ]; dim(dsn)  # 205
dong_epidem = dong_epidem[(dong_epidem$Dis1 %in% common_icds) & (dong_epidem$Dis2 %in% common_icds), ]; dim(dong_epidem) # 117
length(unique(union(dong_epidem$Dis1, dong_epidem$Dis2))) # 28 -- "C71" "C73" are the ones wo links
length(unique(union(dsn$Dis1, dsn$Dis2))) # 30

# Compute the overlap
dsn_graph <- graph_from_data_frame(dsn, directed=FALSE)
dong_epidem_graph <- graph_from_data_frame(dong_epidem, directed=FALSE)
intersect <- intersection(dsn_graph,dong_epidem_graph)
overlap = gsize(intersect) # Number of edges # 55

precision = (overlap/(dim(dsn)[1]))*100; precision      # 26.82927 PRECISION
recall = (overlap/(dim(dong_epidem)[1]))*100; recall    # 47.00855 RECALL


## Generate a network with all my possible interactions between common diseases
N <- length(common_icds) # 30

connected_net <- data.frame(Dis1=c(),Dis2=c(),Distance=c())
set.seed(5)
for (k1 in 1:(N-1)){
  for (k2 in (k1+1):N){
    if(is.na(common_icds[k2])== T){
      print(k1)
      print(k2)
      break
    }
    connected_net <- rbind(connected_net,data.frame(common_icds[k1],common_icds[k2],1))
  }
}

## Sort the connected network
colnames(connected_net) <- c('Dis1','Dis2','Distance')
connected_net$Dis1 = as.character(connected_net$Dis1); connected_net$Dis2 = as.character(connected_net$Dis2)
connected_net <- sort_icd_interactions(connected_net)
dim(connected_net) # ((30*30) - 30)/2
# 435 possible interactions

## Generate 10.000 iterations of: select random network with my nº of interactions and compute overlap with Dong et al. 
## Shuffling labels and keeping the degree distribution of the graph
n_iters <- 10000 
random_overlaps <- c()
random_overlaps_second <- c()
set.seed(5) ### ESTABLECER SEMILLA

for(k in 1:n_iters){
  # p-value for the comparison with Dong et al. (RECALL)
  cgraph <- rewire(dsn_graph, with=keeping_degseq(loops=FALSE, niter = ecount(dsn_graph)*10))
  # str(cgraph)
  cintersect <- intersection(cgraph,dong_epidem_graph)
  # cintersect
  gsize(cintersect)
  random_overlaps <- append(random_overlaps,gsize(cintersect))
  # print_all(rewire(mygraph, with = keeping_degseq(loops=FALSE, niter = vcount(mygraph) * 10)))
}

for(k in 1:n_iters){ # it has to go into another loop to keep consistent p-values (reproducibility)
  # p-value for the comparison with the DSN (PRECISION)
  cgraph_second <- rewire(dong_epidem_graph, with=keeping_degseq(loops=FALSE, niter = ecount(dong_epidem_graph)*10))
  # str(cgraph)
  cintersect_second <- intersection(cgraph_second,dsn_graph)
  # cintersect
  gsize(cintersect_second)
  random_overlaps_second <- append(random_overlaps_second,gsize(cintersect_second))
  # print_all(rewire(mygraph, with = keeping_degseq(loops=FALSE, niter = vcount(mygraph) * 10)))
}

# Comparison with Dong et al. --> Recall
pval <- length(which(random_overlaps>=gsize(intersect)))/length(random_overlaps)
pval
# 0.326

# Comparison with DSN --> Precision
pval_second <- length(which(random_overlaps_second>=gsize(intersect)))/length(random_overlaps_second)
pval_second
# 0.6076

