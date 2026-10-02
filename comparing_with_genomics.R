#####################################################################################
# COMPARING WITH THE GENOME MEDICINE EPIDEMIOLOGICAL NETWORK
# ANALYZING GENETIC VS. NON-GENETIC ORIGIN OF DISEASE CO-OCCURRENCES
#  
# 
# Beatriz Urda García 2021, 2026
######################################################################################
  

setwd('~/Desktop/ANALYSIS/Genome_medicine/')
source('ukb_comparison_tools.R')

args = commandArgs(trailingOnly = TRUE)

if(length(args) != 5){
  stop("Please provide the arguments")
}else{
  network_type = args[1] # dsn / ssn
  split_duplicates = args[2] # TRUE/FALSE
  only_interpretable_diseases = args[3] # TRUE/FALSE
  filt=args[4] # To test Leave-one-out. By default it is ""
  sign_threshold = as.numeric(args[5]) # By default it is 0.01 (1% FDR)
  
}

carguments = paste("",network_type, filt, sign_threshold, "split_dup", split_duplicates, "only_interp_dis", only_interpretable_diseases, sep="_")

# Write results in a file
outputname = paste0("comparisons/comparison_",network_type,"_","splitduplicates_",split_duplicates,"_only_interprt_diseases_",only_interpretable_diseases,"_",filt,"_",sign_threshold,".txt")
write("COMPARISON WITH THE GENOME MEDICINE\n", file=outputname, append=FALSE)
write(paste("Split duplicated:",split_duplicates), file=outputname, append=TRUE)


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
  dsn = read.csv2('../../GEVariability/Network_building/Defined_networks/all_diseases/DEGs/pairwise_union_cosine_distance_sDEGs_pos_network.txt',
                  stringsAsFactors = F,sep="\t",header=T)
  nrow(dsn)
  
  if(filt != ""){
    int_to_remove = union(which(dsn$Dis1 %in% filt), which(dsn$Dis2 %in% filt))
    if(length(int_to_remove) > 0){
      dsn = dsn[-int_to_remove, ]
    }
  }

  if(sign_threshold != 0.05){
    print("Interactions network - before / after filtering: ")
    print(nrow(dsn))
    dsn$adj_pvalue = as.numeric(as.character(dsn$adj_pvalue))
    dsn = dsn[dsn$adj_pvalue <= 0.01, ]
  }
  print(nrow(dsn))
  
  # Transform into icd10 codes
  dsn = network_from_dis_to_icd10(dsn)
  
  dups_index = grep(";",newmeta$icd10)
  dups = trimws(newmeta[dups_index, ]$icd10)
  for(dup in dups){
    print(dup)
    dsn = split_dup_in_network(dsn, dup)
  }
  
  # Split obesity_t2b interactions in two
  # dsn = split_obesity_t2b_in_network(dsn); nrow(dsn)
  length(unique(union(dsn$Dis1, dsn$Dis2))) # 88
  
}else if(network_type == "ssn"){
  # SSN
  ssn_filename = "../Network_building/Defined_networks/metapatients_and_disease/metap_dis_pairwise_union_spearman_distance_sDEGs_pos_network.txt"
  
  # Selecting D-M network with disease names (BreastCancer_1 --> BreastCancer)
  ssn = get_dis_representation_of_ssn(ssn_filename); nrow(ssn)
  
  # Transforming SSN into ICD10 codes
  ssn = network_from_dis_to_icd10(ssn); nrow(ssn)
  dsn = ssn
}else{  # new_ssn   # Check that it works!
  print("new_ssn")
  # SSN
  ssn_filename = "../../GEVariability/Network_building/Defined_networks/all_diseases/metapatients_and_disease/metap_dis_pairwise_union_cosine_distance_sDEGs_pos_network.txt"
  
  # Selecting D-M network with disease names (BreastCancer_1 --> BreastCancer)
  ssn = get_dis_representation_of_ssn(ssn_filename); nrow(ssn)
  
  if(filt != ""){
    int_to_remove = union(which(ssn$Dis1 %in% filt), which(ssn$Dis2 %in% filt))
    if(length(int_to_remove) > 0){
      ssn = ssn[-int_to_remove, ]
    }
  }
  
  
  if(sign_threshold != 0.05){
    print("Interactions network - before / after filtering: ")
    print(nrow(ssn))
    ssn$adj_pvalue = as.numeric(as.character(ssn$adj_pvalue))
    ssn = ssn[ssn$adj_pvalue <= 0.01, ]
  }
  print(nrow(ssn))
  
  # Transforming SSN into ICD10 codes
  ssn = network_from_dis_to_icd10(ssn); nrow(ssn)
  dsn = ssn
  
  dups_index = grep(";",newmeta$icd10)
  dups = trimws(newmeta[dups_index, ]$icd10)
  for(dup in dups){
    print(dup)
    dsn = split_dup_in_network(dsn, dup)
  }
  
}


write(paste("\nDim DSN:",nrow(dsn)), file=outputname, append=TRUE)
write(paste("Dim Dong et al:",nrow(dong_epidem)), file=outputname, append=TRUE)

check_if_dong_network_is_correct(dong_epidem)

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


nrow(dsn)
nrow(dong_epidem)
# Sort the interactions
dsn = sort_icd_interactions(dsn, sort_names = TRUE); nrow(dsn)
dong_epidem = sort_icd_interactions(dong_epidem, sort_names = TRUE, disease_names = c("Description1","Description2")); nrow(dong_epidem)
check_if_dong_network_is_correct(dong_epidem) # Perfect: 
# [1] "F00" "F01" "F03" "G30"
# [1] 4

head(dsn)
head(dong_epidem)
length(unique(union(dong_epidem$Dis1, dong_epidem$Dis2))) # 438

dim(dsn) # 417
dim(dong_epidem) # 11285

# Remove duplicated interactions
dsn = remove_duplicated_interactions(dsn) # 354 unique interactions
dong_epidem = remove_duplicated_interactions(dong_epidem) # 11285 - No duplicated interactions
# check_if_dong_network_is_correct(dong_epidem)

write("\nRemoving duplicated interactions:", file=outputname, append=TRUE)
write(paste("Dim DSN:",dim(dsn)[1]), file=outputname, append=TRUE)
write(paste("Dim Dong et al:",dim(dong_epidem)[1]), file=outputname, append=TRUE)

# Remove interactions where icd Dis1 == icd Dis2
dsn = dsn[dsn$Dis1 != dsn$Dis2, ]; dim(dsn) # 353 unique interactions
dong_epidem = dong_epidem[dong_epidem$Dis1 != dong_epidem$Dis2, ]; dim(dong_epidem) # 11285 - the same as before

write("\nRemoving interactions where icd Dis1 == icd Dis2:", file=outputname, append=TRUE)
write(paste("Dim DSN:",dim(dsn)[1]), file=outputname, append=TRUE)
write(paste("Dim Dong et al:",dim(dong_epidem)[1]), file=outputname, append=TRUE)

# Saving the Dong et al network with all the interactions
dong_epidem_all = dong_epidem

# Common ICD10s
icds_dsn = unique(union(dsn$Dis1, dsn$Dis2)); length(icds_dsn) # 41
icds_dong = unique(union(dong_epidem$Dis1, dong_epidem$Dis2)); length(icds_dong) # 467 --> 438

common_icds = intersect(icds_dsn, icds_dong); length(common_icds) # 30 --> 72

# Save the icds that are only in Dong
# icds_only_dong = setdiff(icds_dong, icds_dsn)
# write.table(icds_only_dong, file="icds_only_dong.txt", sep="\t", quote = FALSE)

write("\nNumber of unique ICD10s:", file=outputname, append=TRUE)
write(paste("DSN:",length(icds_dsn)), file=outputname, append=TRUE)
write(paste("Dong et al:", length(icds_dong)), file=outputname, append=TRUE)
write(paste("In common:", length(common_icds) ), file=outputname, append=TRUE)
write("They are:\n", file=outputname, append=TRUE)
write(paste(common_icds, collapse = "\t"), file=outputname, append=TRUE)

# Keep networks entailing common icds
dsn = dsn[(dsn$Dis1 %in% common_icds) & (dsn$Dis2 %in% common_icds), ]; dim(dsn)  # 205
dong_epidem = dong_epidem[(dong_epidem$Dis1 %in% common_icds) & (dong_epidem$Dis2 %in% common_icds), ]; dim(dong_epidem) # 117
length(unique(union(dong_epidem$Dis1, dong_epidem$Dis2))) # 28 -- "C71" "C73" are the ones wo links
length(unique(union(dsn$Dis1, dsn$Dis2))) # 30

write("\nKeep networks entailing common icds:", file=outputname, append=TRUE)
write(paste("Dim DSN:",nrow(dsn)), file=outputname, append=TRUE)
write(paste("Dim Dong et al:",nrow(dong_epidem)), file=outputname, append=TRUE)

# Compute the overlap
dsn_graph <- graph_from_data_frame(dsn, directed=FALSE)
dong_epidem_graph <- graph_from_data_frame(dong_epidem, directed=FALSE)
intersect <- intersection(dsn_graph,dong_epidem_graph)
overlap = gsize(intersect) # Number of edges # 55

precision = (overlap/(dim(dsn)[1]))*100; precision      # 26.82927 PRECISION
recall = (overlap/(dim(dong_epidem)[1]))*100; recall    # 47.00855 RECALL

write("\nComputing the overlap:", file=outputname, append=TRUE)
write(paste("Precision:",precision), file=outputname, append=TRUE)
write(paste("Recall:",recall), file=outputname, append=TRUE)

#### Assign a significance to the overlap

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

try(dev.off())
outputname_plots = paste0("Plots/Upset_plots_and_euler_diagrams","_",network_type,"_splitduplicates_",split_duplicates,"_only_interprt_diseases_",only_interpretable_diseases,"_",filt,"_",sign_threshold,".pdf")
pdf(file=outputname_plots, width = 9, height = 8)

# Fisher's exact test
dim_dsn = dim(dsn)[1] # 205
all_possible_int = dim(connected_net)[1] # 435
all_possible_graph = graph_from_data_frame(connected_net, directed=FALSE)
not_dsn_graph = difference(all_possible_graph, dsn_graph); gsize(not_dsn_graph)   # 230
not_dsn_in_dong = intersection(not_dsn_graph,dong_epidem_graph); gsize(not_dsn_in_dong)   # 62

## Contingency table
int_dsn = c(overlap, (dim_dsn - overlap))
int_not_dsn = c(gsize(not_dsn_in_dong), (gsize(not_dsn_graph)-gsize(not_dsn_in_dong)))

fsh_test = fisher.test(as.matrix(cbind(int_dsn, int_not_dsn)), alternative = "greater") 

write("\nFisher's exact test:", file=outputname, append=TRUE)
write(paste("pv:",fsh_test$p.value), file=outputname, append=TRUE)
write(paste("Odds's ratio:",fsh_test$estimate), file=outputname, append=TRUE)

####### Of the disease co-occurrences not recapitulated by Dong et al. (i.e., not genetically explained), how many are captured in our network? #######
library(readxl)
library(UpSetR) # Upset plot
library(eulerr) # Euler diagram

gene_df = as.data.frame(read_excel("Paper_and_data/gene interpretable multimorbidity.xlsx"))
geneticcorr_df = as.data.frame(read_excel("Paper_and_data/genetic correlation interpretable multimorbidity.xlsx"))
pathways_df = as.data.frame(read_excel("Paper_and_data/pathway interpretable multimorbidity.xlsx"))
ppi_df = as.data.frame(read_excel("Paper_and_data/ppi interpretable multimorbidity.xlsx"))
snp_df = as.data.frame(read_excel("Paper_and_data/snp interpretable multimorbidity.xlsx"))

# Remove SNP entries with missing Disease2 information
dim(snp_df)
snp_df = snp_df[!is.na(snp_df$Disease2), ]
dim(snp_df)

# Standardize disease column names across dfs to ensure compatibility with downstream functions
colnames(gene_df)[1:2] = c("Dis1", "Dis2")
colnames(geneticcorr_df)[1:2] = c("Dis1", "Dis2")
colnames(pathways_df)[1:2] = c("Dis1", "Dis2")
colnames(ppi_df)[1:2] = c("Dis1", "Dis2")
colnames(snp_df)[1:2] = c("Dis1", "Dis2")

# Split duplicates
if(split_duplicates == TRUE){
  dim(gene_df) # 1463
  gene_df = split_duplicates_in_dong_et_al(gene_df); dim(gene_df) # 1700
  dim(geneticcorr_df) # 1970
  geneticcorr_df = split_duplicates_in_dong_et_al(geneticcorr_df); dim(geneticcorr_df) # 2322
  dim(pathways_df) # 1959
  pathways_df = split_duplicates_in_dong_et_al(pathways_df); dim(pathways_df) # 2218
  dim(ppi_df) # 1959
  ppi_df = split_duplicates_in_dong_et_al(ppi_df); dim(ppi_df) # 2218
  dim(snp_df) # 147
  snp_df = split_duplicates_in_dong_et_al(snp_df); dim(snp_df) # 158
}else{
  dim(gene_df) # 1463
  dim(geneticcorr_df) # 1970
  dim(pathways_df) # 1959
  dim(ppi_df) # 1959
  dim(snp_df) # 147
}


write("\nMOLECULAR LAYERS:", file=outputname, append=TRUE)
write("\nDealing with duplicated interactions accordingly:", file=outputname, append=TRUE)
write(paste("Dim Genes:", dim(gene_df)[1]), file=outputname, append=TRUE)
write(paste("Dim Genetic Corr:", dim(geneticcorr_df)[1]), file=outputname, append=TRUE)
write(paste("Dim Pathways:", dim(pathways_df)[1]), file=outputname, append=TRUE)
write(paste("Dim PPIs:", dim(ppi_df)[1]), file=outputname, append=TRUE)
write(paste("Dim SNPs:", dim(snp_df)[1]), file=outputname, append=TRUE)

##### Get the diseases in the comorbidity networks explained by each molecular layer ####
# See if they match the diseases with described molecular information. 

dis_snps = unique(union(snp_df$Dis1, snp_df$Dis2)); length(dis_snps) # 63
dis_genes = unique(union(gene_df$Dis1, gene_df$Dis2)); length(dis_genes) # 245
dis_pathways = unique(union(pathways_df$Dis1, pathways_df$Dis2)); length(dis_pathways) # 256
dis_ppi = unique(union(ppi_df$Dis1, ppi_df$Dis2)); length(dis_ppi) # 264
dis_genetic_corr = unique(union(geneticcorr_df$Dis1, geneticcorr_df$Dis2)); length(dis_genetic_corr) # 176

interp1 = unique(Reduce(union,list(dis_snps, dis_genes, dis_pathways, dis_ppi, dis_genetic_corr))) # 310
interp2 = get_interpretable_diseases_from_Dong()
length(interp1) # 310
length(interp2) # 330


setdiff(interp1, interp2)
#  "N99" ; 
# "N99" %in% dis_genetic_corr # TRUE
setdiff(interp2, interp1)
# [1] "N63"     "K07"     "I60"     "L05"     "F33"     "N76"     "N90"     "N89"     "N02"    
# [10] "M62"     "N97"     "N86"     "N64"     "L90;L91" "H90"     "N70"     "L89"     "K13"    
# [19] "M89"     "N88"     "L81"

interpretable_diseases = union(interp1, interp2); length(interpretable_diseases) # 331 --> splitting diseases: 375
common_interpretable_diseases = length(intersect(interpretable_diseases, common_icds)) # 57 --> splitting diseases: 57

dsn_icds = union(dsn$Dis1, dsn$Dis2)
eu = euler(list(dsn=dsn_icds, interpretable=interpretable_diseases), shape="ellipse")
plot(eu, quantities=TRUE, main="Number of diseases") # 6 of the ICD10s are not interpretable diseases!!!!
#####

# Sort the dataframes
gene_df = sort_icd_interactions(gene_df)
geneticcorr_df = sort_icd_interactions(geneticcorr_df)
pathways_df = sort_icd_interactions(pathways_df)
ppi_df = sort_icd_interactions(ppi_df)
snp_df = sort_icd_interactions(snp_df)

# Create a column with the explainable pairs
gene_df$pairs = paste(gene_df$Dis1, gene_df$Dis2, sep="_")
geneticcorr_df$pairs = paste(geneticcorr_df$Dis1, geneticcorr_df$Dis2, sep="_")
pathways_df$pairs = paste(pathways_df$Dis1, pathways_df$Dis2, sep="_")
ppi_df$pairs = paste(ppi_df$Dis1, ppi_df$Dis2, sep="_")
snp_df$pairs = paste(snp_df$Dis1, snp_df$Dis2, sep="_")

# Create a Venn Diagram of them
dim(dong_epidem_all); dim(dong_epidem) # 13002 # 117
dong_epidem_all$pairs = paste(dong_epidem_all$Dis1, dong_epidem_all$Dis2, sep="_")
dsn$pairs = paste(dsn$Dis1, dsn$Dis2, sep="_")

# Create Dong's dictionary of ICD10s
dong_dis_dic = as.matrix(dong_epidem_all[, c("Dis1", "Description1", "Category1")]); nrow(dong_dis_dic)
dong_dis_dic = rbind(dong_dis_dic, as.matrix(dong_epidem_all[, c("Dis2", "Description2", "Category2")])); nrow(dong_dis_dic)
dong_dis_dic = unique(dong_dis_dic); nrow(dong_dis_dic)
dong_dis_dic = as.data.frame(dong_dis_dic)
colnames(dong_dis_dic) = c("icd10", "Disease_name_dong", "Disease_category")
length(unique(dong_dis_dic$Disease_category)) # 24
dong_dis_dic_cat = dong_dis_dic; nrow(dong_dis_dic_cat) # Contains duplicated ICD10 (the ones that correspond to several categories. E.g. F17)
dong_dis_dic = dong_dis_dic[!duplicated(dong_dis_dic$icd10), c("icd10", "Disease_name_dong")]; nrow(dong_dis_dic_cat)

# Checks
# length(intersect(common_icds, dong_dis_dic$icd10)) # 72
# length(common_icds) # 72

# try(dev.off())
# outputname_plots = paste0("Plots/Upset_plots_and_euler_diagrams","_",network_type,"_splitduplicates_",split_duplicates,"_only_interprt_diseases_",only_interpretable_diseases,".pdf")
# pdf(file=outputname_plots, width = 9, height = 9)

if(only_interpretable_diseases == TRUE){
  dong_epidem_all = dong_epidem_all[(dong_epidem_all$Dis1 %in% interpretable_diseases) & (dong_epidem_all$Dis2 %in% interpretable_diseases), ]
  dong_epidem = dong_epidem[(dong_epidem$Dis1 %in% interpretable_diseases) & (dong_epidem$Dis2 %in% interpretable_diseases), ]
  dsn = dsn[(dsn$Dis1 %in% interpretable_diseases) & (dsn$Dis2 %in% interpretable_diseases), ]
  
  snp_df = snp_df[(snp_df$Dis1 %in% interpretable_diseases) & (snp_df$Dis2 %in% interpretable_diseases), ]
  gene_df = gene_df[(gene_df$Dis1 %in% interpretable_diseases) & (gene_df$Dis2 %in% interpretable_diseases), ]
  pathways_df = pathways_df[(pathways_df$Dis1 %in% interpretable_diseases) & (pathways_df$Dis2 %in% interpretable_diseases), ]
  ppi_df = ppi_df[(ppi_df$Dis1 %in% interpretable_diseases) & (ppi_df$Dis2 %in% interpretable_diseases), ]
  geneticcorr_df = geneticcorr_df[(geneticcorr_df$Dis1 %in% interpretable_diseases) & (geneticcorr_df$Dis2 %in% interpretable_diseases), ]
  
  write("\nKeeping only Interpretable Diseases for Dong, net and molecular Dong:", file=outputname, append=TRUE)
  write("\nNetworks' dimensions::", file=outputname, append=TRUE)
  write(paste("Dim dong_epidem_all:", dim(dong_epidem_all)[1]), file=outputname, append=TRUE)
  write(paste("Dim dong_epidem:", dim(dong_epidem)[1]), file=outputname, append=TRUE)
  write(paste("Dim dsn:", dim(dsn)[1]), file=outputname, append=TRUE)
  
  write(paste("Dim snp_df:", dim(snp_df)[1]), file=outputname, append=TRUE)
  write(paste("Dim gene_df:", dim(gene_df)[1]), file=outputname, append=TRUE)
  write(paste("Dim pathways_df:", dim(pathways_df)[1]), file=outputname, append=TRUE)
  write(paste("Dim ppi_df:", dim(ppi_df)[1]), file=outputname, append=TRUE)
  write(paste("Dim geneticcorr_df:", dim(geneticcorr_df)[1]), file=outputname, append=TRUE)
}else{
  write("\nWithout filtering for Interpretable Diseases:", file=outputname, append=TRUE)
  write("\nNetworks' dimensions::", file=outputname, append=TRUE)
  write(paste("Dim dong_epidem_all:", dim(dong_epidem_all)[1]), file=outputname, append=TRUE)
  write(paste("Dim dong_epidem:", dim(dong_epidem)[1]), file=outputname, append=TRUE)
  write(paste("Dim dsn:", dim(dsn)[1]), file=outputname, append=TRUE)
  
  write(paste("Dim snp_df:", dim(snp_df)[1]), file=outputname, append=TRUE)
  write(paste("Dim gene_df:", dim(gene_df)[1]), file=outputname, append=TRUE)
  write(paste("Dim pathways_df:", dim(pathways_df)[1]), file=outputname, append=TRUE)
  write(paste("Dim ppi_df:", dim(ppi_df)[1]), file=outputname, append=TRUE)
  write(paste("Dim geneticcorr_df:", dim(geneticcorr_df)[1]), file=outputname, append=TRUE)
  
}

#### MUTUAL INFORMATION
df = dong_epidem_all[, c("Dis1", "Dis2", "pairs")]; head(df)
df$snps = NA; df$genes = NA; df$ppis = NA; df$pathways = NA; df$snps = NA; df$geneticcorr = NA; 
df$geneexpr = NA; df$Distance = NA; df$adj_pvalue = NA
for(k in 1:nrow(df)){    # De las interacciones de Dong, are they in the molecular layers?
  cpair = df$pairs[k]
  if(cpair %in% snp_df$pairs){df$snps[k] = 1}else{df$snps[k] = 0}
  if(cpair %in% gene_df$pairs){df$genes[k] = 1}else{df$genes[k] = 0}
  if(cpair %in% pathways_df$pairs){df$pathways[k] = 1}else{df$pathways[k] = 0}
  if(cpair %in% ppi_df$pairs){df$ppis[k] = 1}else{df$ppis[k] = 0}
  if(cpair %in% geneticcorr_df$pairs){df$geneticcorr[k] = 1}else{df$geneticcorr[k] = 0}
  if((df$Dis1[k] %in% common_icds) &  (df$Dis2[k] %in% common_icds)){
    if(cpair %in% dsn$pairs){
      df$geneexpr[k] = 1
      df$Distance[k] = dsn[dsn$pairs == cpair, ]$Distance
      df$adj_pvalue[k] = dsn[dsn$pairs == cpair, ]$adj_pvalue
    }else{df$geneexpr[k] = 0}
  }else{
    df$geneexpr[k] = "-"
  }
}

dong_epidem$pairs = paste(dong_epidem$Dis1, dong_epidem$Dis2, sep="_")

# All the interactions (with the selected parameters in args)
pairs_list = list(snp_df$pairs, gene_df$pairs, pathways_df$pairs, ppi_df$pairs, geneticcorr_df$pairs, intersect(dsn$pairs, dong_epidem_all$pairs))
# All the UKB interactions (with the selected parameters in args)
pairs_list_common = list(intersect(snp_df$pairs, dong_epidem$pairs), intersect(gene_df$pairs, dong_epidem$pairs), 
                         intersect(pathways_df$pairs, dong_epidem$pairs), intersect(ppi_df$pairs, dong_epidem$pairs), 
                         intersect(geneticcorr_df$pairs, dong_epidem$pairs), 
                         intersect(dsn$pairs, dong_epidem$pairs))

nrow(df)
dfcommon = df[df$geneexpr != "-", ]; nrow(dfcommon)
dfcommon$geneexpr = as.numeric(dfcommon$geneexpr)

length(union(df$Dis1, df$Dis2))
length(union(dfcommon$Dis1, dfcommon$Dis2))
length(intersect(union(dfcommon$Dis1, dfcommon$Dis2), dong_dis_dic$icd10)) # 53
setdiff(union(dfcommon$Dis1, dfcommon$Dis2), dong_dis_dic$icd10)

# SAVE DF DATAFRAME FOR MUTUAL INFORMATION LATER ON
df_mutual_info = df

libmeta = read.csv("../../GEVariability/disease_metadata.csv", sep="\t", stringsAsFactors = FALSE)

# Add disease name to dfcommon
dfcommon2 = merge(dfcommon, dong_dis_dic, all.x=TRUE, all.y = FALSE, by.x = "Dis2", by.y="icd10")
dfcommon2 = merge(dfcommon2, dong_dis_dic, all.x=TRUE, all.y = FALSE, by.x = "Dis1", by.y="icd10")
head(dfcommon2)
cncol = ncol(dfcommon2)
dfcommon2 <- dfcommon2[,c(1:(ncol(dfcommon2)-2), ncol(dfcommon2), ncol(dfcommon2)-1)]
colnames(dfcommon2)[(cncol-1):cncol] = c("Dis1_name", "Dis2_name")
write.table(dfcommon2,file=paste0("Results/dong_interactions_in_molecular_networks",carguments,".txt"),sep="\t",row.names=F, col.names=T, quote=FALSE)

# Add disease name to df
df2 = merge(df, dong_dis_dic, all.x=TRUE, all.y = FALSE, by.x = "Dis2", by.y="icd10")
df2 = merge(df2, dong_dis_dic, all.x=TRUE, all.y = FALSE, by.x = "Dis1", by.y="icd10")
head(df2)
cncol = ncol(df2)
df2 <- df2[,c(1:(ncol(df2)-2), ncol(df2), ncol(df2)-1)]
colnames(df2)[(cncol-1):cncol] = c("Dis1_name", "Dis2_name")
write.table(df2,file=paste0("Results/dong_interactions_in_molecular_networks_all",carguments,".txt"),sep="\t",row.names=F, col.names=T, quote=FALSE)

in_mol = dong_dis_dic$icd10 %in% unique(union(dfcommon2$Dis1, dfcommon2$Dis2))
cdislist = dong_dis_dic[in_mol, ]; nrow(cdislist) # 56 --> 61

# CREATE DF WITH THE %EIs EXPLAINED BY EACH MOLECULAR LAYER *BY DISEASE*: diseases_df
diseases_df = data.frame()
for(i in 1:nrow(cdislist)){
  dis = cdislist$icd10[i]; dis  
  dis_name = cdislist$Disease_name_dong[i]
  
  cdf = dfcommon2[dfcommon2$Dis1 == dis | dfcommon2$Dis2 == dis, ]; head(cdf)
  cnint = nrow(cdf)
  counts = colSums(cdf == 1)
  sum_rows_genetic = rowSums(cdf[, c("snps","genes","ppis","pathways","geneticcorr")])
  sum_genetic = length(sum_rows_genetic[sum_rows_genetic >= 1])
  sum_rows_molecular = rowSums(cdf[, c("snps","genes","ppis","pathways","geneticcorr","geneexpr")])
  sum_molecular = length(sum_rows_molecular[sum_rows_molecular >= 1])
  
  diseases_df = rbind(diseases_df, c(dis, dis_name, cnint, (counts["snps"]/cnint)*100, (counts["genes"]/cnint)*100,
                                     (counts["ppis"]/cnint)*100, (counts["pathways"]/cnint)*100, (counts["geneticcorr"]/cnint)*100,
                                     (sum_genetic/cnint)*100,
                                     (counts["geneexpr"]/cnint)*100,
                                     (sum_molecular/cnint)*100))
  
}
colnames(diseases_df) = c("icd10", "disease_name", "N_interactions", "Perc_snps", "Perc_genes", "Perc_ppis", "Perc_pathways", "Perc_geneticcorr", "Perc_genetic","Perc_geneexpr", "Perc_molecular")

# identify the columns to convert to numeric
cols_to_convert = c("N_interactions", "Perc_snps", "Perc_genes", "Perc_ppis", "Perc_pathways", "Perc_geneticcorr", "Perc_genetic","Perc_geneexpr", "Perc_molecular")

# use apply() and as.numeric() to convert the selected columns to numeric
diseases_df[cols_to_convert] <- apply(diseases_df[cols_to_convert], 2, as.numeric)

# check the data types of the columns after conversion
str(diseases_df)
diseases_df$expr_vs_genetic = diseases_df$Perc_geneexpr / diseases_df$Perc_genetic

# Add the disease category: diseases_df_cat is the same as diseases_df but ICD10s with multiple disease categories are splitted in two rows, one per category
# ICDs with multiple disease cats: F17 (tobacco) and F32 (depression)
nrow(diseases_df) # 56
diseases_df_cat = merge(diseases_df, dong_dis_dic_cat[, c("icd10", "Disease_category")], by="icd10", all.x=TRUE, all.y=FALSE)
nrow(diseases_df_cat) 

library(ggrepel)
library(plotly)
library(scales)
library(scales) # for label_number()

# PLOT %EIs BY DISEASE

dev.off()
pdf("Plots/examples.pdf", width = 11, height = 9)
max_ratio = max(diseases_df$expr_vs_genetic[is.finite(diseases_df$expr_vs_genetic)])
g = ggplot(diseases_df[is.finite(diseases_df$expr_vs_genetic), ], aes(x=Perc_genetic, y=Perc_geneexpr, label=disease_name, size=N_interactions, colour=expr_vs_genetic)) +
  theme_classic() +
  scale_x_continuous(limits=c(0, 100)) +
  scale_y_continuous(limits=c(0, 100)) +
  geom_point() + 
  geom_text_repel(size=2.5, box.padding = 0.5) +
  geom_point(data = diseases_df[!is.finite(diseases_df$expr_vs_genetic), ], aes(x=Perc_genetic, y=Perc_geneexpr, size=N_interactions), colour = "#000CB5") +   # muted("#0111FF") -- #16B5AA green color
  geom_text_repel(data=diseases_df[!is.finite(diseases_df$expr_vs_genetic), ], aes(x=Perc_genetic, y=Perc_geneexpr, label=disease_name), size=2.5, colour="#000CB5", box.padding=0.5) +
  scale_colour_gradient2(low="#F2A90A", high="#0111FF", mid="#969696", midpoint=1, na.value = "#969696", space = "Lab", 
                         transform = transform_log(), labels = label_number(accuracy = 0.1))+
  geom_abline(intercept = 0, slope = 1, color = "black", linetype = "dotted")+
  labs(x = "% EIs captured by genetics", y = "% EIs captured by transcriptomics",
       size = "#EIs", colour="Ratio") +
  ggtitle("Percentage of EIs explained by transcriptomics vs. genetics")+ 
  theme(plot.title = element_text(hjust = 0.5))
  
g
ggsave("Plots/map_EIs.pdf", plot = g, width = 10.5, height = 9) 

### Save the colors of the ratio
# Extract the color scale
gg_build <- ggplot_build(g)
color_scale1 <- gg_build$data[[1]]
color_scale2 <- gg_build$data[[3]]
color_scale = rbind(color_scale1, color_scale2); nrow(color_scale) == nrow(diseases_df)
head(color_scale); # 56 diseases
diseases_df = merge(diseases_df, color_scale[, c("label", "colour")], by.x="disease_name", by.y="label", all.x=TRUE, all.y=FALSE)
colnames(diseases_df)[ncol(diseases_df)] = "color_ratio"

# N interactions
# Function to transform T and G to a -1 to 1 scale
transform_ratio_to_linear_scale <- function(t, g) {
  scaled_value <- (g - t) / (g + t)
  return(scaled_value)
}
transform_ratio_to_sigmoid_scale <- function(ratio) {
  log_ratio <- log(ratio)
  scaled_ratio <- (2 / (1 + exp(-log_ratio))) - 1
  return(scaled_ratio)
}
diseases_df$transformed_ratio_linear = transform_ratio_to_linear_scale(diseases_df$Perc_geneexpr, diseases_df$Perc_genetic)
diseases_df$transformed_ratio_sigmoid = transform_ratio_to_sigmoid_scale(diseases_df$expr_vs_genetic)
g = ggplot(diseases_df, aes(x=transformed_ratio_linear, y=N_interactions, label=disease_name, size=N_interactions, colour=expr_vs_genetic)) +
  theme_classic() +
  scale_x_continuous(limits=c(-1, 1)) +
  geom_point() + 
  scale_colour_gradient2(low="#F2A90A", high="#0111FF", mid="#969696", midpoint=1, na.value = "#969696", space = "Lab", 
                         transform = transform_log(), labels = label_number(accuracy = 0.1))+
  labs(x = "Ratio", y = "Number of interactions",
       size = "#EIs", colour="Ratio") +
  ggtitle("Percentage of EIs explained by transcriptomics vs. genetics")+ 
  theme(plot.title = element_text(hjust = 0.5))

g

g = ggplot(diseases_df, aes(x=transformed_ratio_linear, y=N_interactions, label=disease_name, size=Perc_molecular, colour=expr_vs_genetic)) +
  theme_classic() +
  scale_x_continuous(limits=c(-1, 1)) +
  geom_point() + 
  geom_text_repel(size=2.5, box.padding = 0.5) +
  scale_colour_gradient2(low="#F2A90A", high="#0111FF", mid="#969696", midpoint=1, na.value = "#969696", space = "Lab",
                         transform = transform_log(), labels = label_number(accuracy = 0.1))+
  labs(x = "Transcriptomic and genetic contributions", y = "Number of interactions",
       size = "% EIs molecular", colour="Ratio") +
  ggtitle("Percentage of EIs explained by transcriptomics vs. genetics")+ 
  theme(plot.title = element_text(hjust = 0.5))
g


diseases_df$distance_from_zero = abs(diseases_df$transformed_ratio_linear) - 0

# Calculate correlation
correlation_result <- cor.test(abs(diseases_df$distance_from_zero), diseases_df$N_interactions, alternative = "less")
print(correlation_result)
correlation_result <- cor.test(abs(diseases_df$distance_from_zero), diseases_df$N_interactions, method = "spearman", alternative = "less")
print(correlation_result)
correlation_result <- cor.test(abs(diseases_df$distance_from_zero), diseases_df$N_interactions, method = "kendall", alternative = "less")
print(correlation_result)
kendall_tau = correlation_result$estimate; kendall_p_value = correlation_result$p.value
# Kendall's Tau (ττ): -0.2085 --- p-value: 0.02571

# The good one of the # EIs
g = ggplot(diseases_df[is.finite(diseases_df$expr_vs_genetic), ], aes(x=transformed_ratio_linear, y=N_interactions, label=disease_name, size = N_interactions, colour=expr_vs_genetic)) +
  theme_classic() +
  scale_x_continuous(limits=c(-1, 1)) +
  geom_point(data = diseases_df[!is.finite(diseases_df$expr_vs_genetic), ], aes(x=transformed_ratio_linear, y=N_interactions, size=N_interactions), colour = "#000CB5") +
  geom_point() + 
  scale_colour_gradient2(low="#F2A90A", high="#0111FF", mid="#969696", midpoint=1, na.value = "#969696", space = "Lab",transform = transform_log(),
                         labels = label_number(accuracy = 0.01)) + # Rounds to 2 decimals
  
  geom_vline(xintercept = 0, color = "black", linetype = "dotted")+
  labs(x = "Relative transcriptomic and genetic contributions", y = "Number of interactions",
       size = "#EIs", colour="Ratio") +
  theme(plot.title = element_text(hjust = 0.5)) + 
  # Annotate Kendall's Tau with tau symbol and p-value
  annotate("text", x = 0.4, y = max(diseases_df$N_interactions, na.rm = TRUE) - 1, 
           label = bquote(R^2 == 0.11 ~ ", p-value = 0.024"), 
           hjust = 0, size = 3.5, color = "black")+
  annotate("text", x = 0.4, y = max(diseases_df$N_interactions, na.rm = TRUE) - 3, 
           label = bquote(tau == .(round(kendall_tau,3)) * "," ~ "p-value =" ~ .(round(kendall_p_value,3))), 
           hjust = 0, size = 3.5, color = "black")
g
ggsave("Plots/Number_EIs_1.pdf", plot = g, width = 8, height = 7) 

g = ggplot(diseases_df, aes(x=distance_from_zero, y=N_interactions, label=disease_name)) +
  theme_classic() +
  scale_x_continuous(limits=c(0, 1)) +
  geom_point() + 
  labs(x = "Distance to zero", y = "Number of interactions",
       size = "#EIs", colour="Ratio") ++ 
  theme(plot.title = element_text(hjust = 0.5)) + 
  # Annotate Kendall's Tau with tau symbol and p-value
  annotate("text", x = 0.65, y = max(diseases_df$N_interactions, na.rm = TRUE) - 1, 
           label = bquote(tau == .(round(kendall_tau,3)) * "," ~ "p-value =" ~ .(round(kendall_p_value,3))), 
           hjust = 0, size = 4, color = "black")
g
ggsave("Plots/Number_EIs_2.pdf", plot = g, width = 8, height = 7)

g = ggplot(diseases_df, aes(x=Perc_geneticcorr, y=N_interactions, label=disease_name)) +
  theme_classic() +
  geom_point() + 
  labs(x = "Perc_geneticcorr", y = "Number of interactions",
       size = "#EIs", colour="Ratio") +
  theme(plot.title = element_text(hjust = 0.5)) + 
  # Annotate Kendall's Tau with tau symbol and p-value
  annotate("text", x = 0.65, y = max(diseases_df$N_interactions, na.rm = TRUE) - 1, 
           label = bquote(tau == .(round(kendall_tau,3)) * "," ~ "p-value =" ~ .(round(kendall_p_value,3))), 
           hjust = 0, size = 4, color = "black")
g

# Extended version of the previous plot
g = ggplot(diseases_df[is.finite(diseases_df$expr_vs_genetic), ], aes(x=distance_from_zero, y=N_interactions, label=disease_name, size = N_interactions, colour=expr_vs_genetic)) +
  theme_classic() +
  scale_x_continuous(limits=c(-1, 1)) +
  geom_point(data = diseases_df[!is.finite(diseases_df$expr_vs_genetic), ], aes(x=distance_from_zero, y=N_interactions, size=N_interactions), colour = "#000CB5") +
  geom_point() + 
  scale_colour_gradient2(low="#F2A90A", high="#0111FF", mid="#969696", midpoint=1, na.value = "#969696", space = "Lab", transform = transform_log(),
                         labels = label_number(accuracy = 0.01)) + # Rounds to 2 decimals
  
  geom_vline(xintercept = 0, color = "black", linetype = "dotted")+
  labs(x = "Transcriptomic and genetic contributions", y = "Number of interactions",
       size = "#EIs", colour="Ratio") +
  theme(plot.title = element_text(hjust = 0.5)) + 
  # Annotate Kendall's Tau with tau symbol and p-value
  annotate("text", x = 0.4, y = max(diseases_df$N_interactions, na.rm = TRUE) - 3, 
           label = bquote(R^2 == 0.11 ~ ", p-value = 0.024"), 
           hjust = 0, size = 3.5, color = "black")+
  annotate("text", x = 0.4, y = max(diseases_df$N_interactions, na.rm = TRUE) - 1, 
           label = bquote(tau == .(round(kendall_tau,3)) * "," ~ "p-value =" ~ .(round(kendall_p_value,3))), 
           hjust = 0, size = 3.5, color = "black")
g


## Fitting models
# Fit a Quadratic Regression Model
# Create a squared term for 'Ratio'
diseases_df$Ratio_squared <- diseases_df$transformed_ratio_linear^2

# Fit a quadratic regression model
model_quadratic <- lm(N_interactions ~ transformed_ratio_linear + Ratio_squared, data = diseases_df)

# Check the summary of the model
summary(model_quadratic)

# Plot data points
ggplot(diseases_df, aes(x = transformed_ratio_linear, y = N_interactions)) +
  geom_point() +
  # Add the quadratic regression line
  stat_smooth(method = "lm", formula = y ~ x + I(x^2), color = "blue", se = TRUE) +
  labs(x = "Ratio", y = "Number of Interactions",
       title = "Quadratic Relationship between Ratio and Number of Interactions") +
  theme_minimal()


library(mgcv)

# Fit a GAM model
model_gam <- gam(N_interactions ~ s(transformed_ratio_linear), data = diseases_df)
summary(model_gam)

# Visualize the GAM fit
ggplot(diseases_df, aes(x = transformed_ratio_linear, y = N_interactions)) +
  geom_point() +
  stat_smooth(method = "gam", formula = y ~ s(x), color = "green", se = TRUE) +
  labs(x = "Ratio", y = "Number of Interactions",
       title = "GAM Fit for Relationship between Ratio and Number of Interactions") +
  theme_minimal()

# Fit a GAM model with smooth terms
model_gam <- gam(N_interactions ~ s(Perc_geneexpr) + s(Perc_genetic) + s(Perc_molecular), 
                 data = diseases_df)

# Summarize the model
summary(model_gam)

# Fit a GAM model with smooth terms
model_gam <- gam(N_interactions ~ s(Perc_geneexpr) + s(Perc_genetic), 
                 data = diseases_df)

# Summarize the model
summary(model_gam)

# Fit a GAM model with smooth terms
model_gam <- gam(N_interactions ~ s(Perc_geneexpr) + s(Perc_genetic) + s(Perc_geneexpr*Perc_genetic) + s(Perc_molecular), 
                 data = diseases_df)

# Summarize the model
summary(model_gam)

# Fit a GAM model with smooth terms
model_gam <- gam(N_interactions ~ s(Perc_geneexpr*Perc_genetic) + s(Perc_molecular), 
                 data = diseases_df)

# Summarize the model
summary(model_gam)

# Fit a GAM model with smooth terms
model_gam <- gam(N_interactions ~ s(Perc_molecular), 
                 data = diseases_df)

# Summarize the model
summary(model_gam)

# Fit a GAM model with smooth terms
model_gam <- gam(N_interactions ~ s(Perc_molecular) + s(Perc_geneexpr*Perc_genetic), 
                 data = diseases_df)

# Summarize the model
summary(model_gam)

# Fit a GAM model with smooth terms
model_gam <- gam(N_interactions ~ s(Perc_molecular) + s(Perc_geneticcorr), 
                 data = diseases_df)

# Summarize the model
summary(model_gam)

# Fit a GAM model with smooth terms
model_gam <- gam(N_interactions ~ s(Perc_molecular) + s(Perc_geneticcorr) + s(Perc_geneexpr*Perc_genetic), 
                 data = diseases_df)

# Summarize the model
summary(model_gam)

# Fit a GAM model with smooth terms
model_gam <- gam(N_interactions ~ s(Perc_snps)+s(Perc_genes)+s(Perc_geneticcorr)+s(Perc_genetic)+s(Perc_geneexpr), 
                 data = diseases_df)

# Summarize the model
summary(model_gam)

# Fit a GAM model with smooth terms
model_gam <- gam(N_interactions ~ s(Perc_geneticcorr)+Perc_molecular, 
                 data = diseases_df)

# Summarize the model
summary(model_gam)

# Fit a GAM model with smooth terms
model_gam <- gam(N_interactions ~ s(Perc_molecular) + s(transformed_ratio_linear), 
                 data = diseases_df)

# Summarize the model
summary(model_gam)

# Fit a GAM model with smooth terms
model_gam <- gam(N_interactions ~ s(distance_from_zero), 
                 data = diseases_df)

# Summarize the model
summary(model_gam)

# Correlation between Perc_<label> and N_interactions
# SNPs
correlation_result <- cor.test(abs(diseases_df$Perc_snps), diseases_df$N_interactions, method = "spearman", alternative = "greater")
print(correlation_result)
correlation_result <- cor.test(abs(diseases_df$Perc_snps), diseases_df$N_interactions, method = "kendall", alternative = "greater")
print(correlation_result)
# kendall_tau = correlation_result$estimate; kendall_p_value = correlation_result$p.value

# SNPs
correlation_result <- cor.test(abs(diseases_df$Perc_snps), diseases_df$N_interactions, method = "spearman", alternative = "greater")
print(correlation_result)
correlation_result <- cor.test(abs(diseases_df$Perc_snps), diseases_df$N_interactions, method = "kendall", alternative = "greater")
print(correlation_result)

# Perc_genes
correlation_result <- cor.test(abs(diseases_df$Perc_genes), diseases_df$N_interactions, method = "spearman", alternative = "greater")
print(correlation_result)
correlation_result <- cor.test(abs(diseases_df$Perc_genes), diseases_df$N_interactions, method = "kendall", alternative = "greater")
print(correlation_result)

# Perc_ppis
correlation_result <- cor.test(abs(diseases_df$Perc_ppis), diseases_df$N_interactions, method = "spearman", alternative = "greater")
print(correlation_result)
correlation_result <- cor.test(abs(diseases_df$Perc_ppis), diseases_df$N_interactions, method = "kendall", alternative = "greater")
print(correlation_result)

# Perc_pathways
correlation_result <- cor.test(abs(diseases_df$Perc_pathways), diseases_df$N_interactions, method = "spearman", alternative = "greater")
print(correlation_result)
correlation_result <- cor.test(abs(diseases_df$Perc_pathways), diseases_df$N_interactions, method = "kendall", alternative = "greater")
print(correlation_result)

# Perc_geneticcorr
correlation_result <- cor.test(abs(diseases_df$Perc_geneticcorr), diseases_df$N_interactions, method = "spearman", alternative = "greater")
print(correlation_result)
correlation_result <- cor.test(abs(diseases_df$Perc_geneticcorr), diseases_df$N_interactions, method = "kendall", alternative = "greater")
print(correlation_result)

# Perc_genetic
correlation_result <- cor.test(abs(diseases_df$Perc_genetic), diseases_df$N_interactions, method = "spearman", alternative = "greater")
print(correlation_result)
correlation_result <- cor.test(abs(diseases_df$Perc_genetic), diseases_df$N_interactions, method = "kendall", alternative = "greater")
print(correlation_result)

# Perc_geneexpr
correlation_result <- cor.test(abs(diseases_df$Perc_geneexpr), diseases_df$N_interactions, method = "spearman", alternative = "greater")
print(correlation_result)
correlation_result <- cor.test(abs(diseases_df$Perc_geneexpr), diseases_df$N_interactions, method = "kendall", alternative = "greater")
print(correlation_result)

# Perc_molecular
correlation_result <- cor.test(abs(diseases_df$Perc_molecular), diseases_df$N_interactions, method = "spearman", alternative = "greater")
print(correlation_result)
correlation_result <- cor.test(abs(diseases_df$Perc_molecular), diseases_df$N_interactions, method = "kendall", alternative = "greater")
print(correlation_result)

library(scatterpie)
# Add a new column for the complementary value (100 - pie_value)
diseases_df$pie_value_complement <- 100 - diseases_df$Perc_molecular

g <- ggplot(diseases_df, aes(x = transformed_ratio_linear, y = N_interactions, label = disease_name)) +
  theme_classic() +
  scale_x_continuous(limits = c(-1, 1)) +
  geom_scatterpie(aes(x = transformed_ratio_linear, y = N_interactions), 
                  # data = diseases_df, 
                  cols = c("Perc_molecular","pie_value_complement"), 
                  # long_format=TRUE,
                  donut_radius=.5,
                  color = NA) +  # `color = NA` removes the pie outline
  # scale_colour_gradient2(low = "#F2A90A", high = "#0111FF", mid = "#969696", midpoint = 1, na.value = "#969696", space = "Lab", transform = transform_log()) +
  labs(x = "Ratio", y = "Number of interactions",
       size = "#EIs", colour = "Ratio") +
  ggtitle("Percentage of EIs explained by transcriptomics vs. genetics") +
  theme(plot.title = element_text(hjust = 0.5))

print(g)

ggplot() + 
  scale_x_continuous(limits = c(-1, 1)) +
  # scale_y_continuous(limits = c(0, 50)) +
  geom_scatterpie(aes(x=transformed_ratio_linear, y=N_interactions), data=diseases_df,
                               cols=c("Perc_molecular","pie_value_complement")) + 
                              coord_fixed(ratio = 0.2)   # Adjust this ratio to balance x and y scales
                              # + coord_equal() 
                           
diseases_df$transformed_ratio_linear_scaled = diseases_df$transformed_ratio_linear *20
head(diseases_df$transformed_ratio_linear_scaled)

ggplot() + geom_scatterpie(aes(x=transformed_ratio_linear_scaled, y=N_interactions), data=diseases_df,
                           cols=c("Perc_molecular","pie_value_complement")) + coord_equal()  


library(tidyr)
d <- data.frame(x=rnorm(5), y=rnorm(5))
d$A <- abs(rnorm(5, sd=1))
d$B <- abs(rnorm(5, sd=2))
d$C <- abs(rnorm(5, sd=3))
d <- tidyr::gather(d, key="letters", value="value", -x:-y)
ggplot() + geom_scatterpie(aes(x=x, y=y), data=d, cols="letters", long_format=TRUE) + coord_fixed()
p1 <- ggplot() + 
  geom_scatterpie(
    mapping = aes(x=x, y=y), data=d, cols="letters", 
    long_format=TRUE, 
    donut_radius=.5
  ) + 
  coord_fixed()
p1


g = ggplot(diseases_df, aes(x=transformed_ratio_sigmoid, y=Perc_molecular, label=disease_name, size=N_interactions, colour=Perc_molecular)) +
  theme_classic() +
  scale_x_continuous(limits=c(-1, 1)) +
  geom_point() + 
  labs(x = "Ratio", y = "Number of interactions",
       size = "#EIs", colour="Ratio") +
  ggtitle("Percentage of EIs explained by transcriptomics vs. genetics")+ 
  theme(plot.title = element_text(hjust = 0.5))

g

if(network_type == "new_ssn"){
  g = ggplot(diseases_df[diseases_df$transformed_ratio_linear, ], aes(x=transformed_ratio_linear, y=N_interactions, label=disease_name, size=N_interactions, colour=Perc_molecular)) +
    theme_classic() +
    scale_x_continuous(limits=c(-1, 1)) +
    geom_point() + 
    labs(x = "Ratio", y = "Number of interactions",
         size = "#EIs", colour="Ratio") +
    ggtitle("Percentage of EIs explained by transcriptomics vs. genetics")+ 
    theme(plot.title = element_text(hjust = 0.5))
  
  g
  
  g = ggplot(diseases_df[diseases_df$transformed_ratio_linear, ], aes(x=transformed_ratio_linear, y=Perc_molecular, label=disease_name, size=N_interactions, colour=Perc_molecular)) +
    theme_classic() +
    scale_x_continuous(limits=c(-1, 1)) +
    geom_point() + 
    labs(x = "Ratio", y = "Number of interactions",
         size = "#EIs", colour="Ratio") +
    ggtitle("Percentage of EIs explained by transcriptomics vs. genetics")+ 
    theme(plot.title = element_text(hjust = 0.5))
  
  g
  
}


# Categorize into groups
group <- cut(diseases_df$transformed_ratio_linear, breaks = c(-1, -0.5, 0.5, 1), labels = c("Low", "Mid", "High"))

# Check if Y differs across groups
kruskal.test(diseases_df$N_interactions ~ group)

model <- lm(diseases_df$N_interactions ~ diseases_df$transformed_ratio_linear + I(diseases_df$transformed_ratio_linear^2))
summary(model)

###### HERITABILITY
# Add the heritability by Westergaard et al. 
west <- read.delim("Paper_and_data/website_icd10_lvl3.tsv", sep = "\t", header = TRUE)
head(west)
colnames(west) = c("icd10", "term","disease_name", "chapter","chapter_name","h2_west", "C1_west", "C2_west")

westA = west[west$term == "A", ]
westSp = west[west$term == "Sp", ]
westA = westA[, c("icd10","h2_west")]
westSp = westSp[, c("icd10","h2_west")]
colnames(westSp) = c("icd10", "h2_west_Sp")

dim(diseases_df)
diseases_df = merge(diseases_df, westA, by='icd10', all.x = TRUE, all.y=FALSE)
dim(diseases_df)
diseases_df = merge(diseases_df, westSp, by='icd10', all.x = TRUE, all.y=FALSE)
dim(diseases_df)

g = ggplot(diseases_df, aes(x=Perc_genetic, y=Perc_geneexpr, label=disease_name, size=N_interactions, colour=h2_west)) +
  theme_classic() +
  scale_x_continuous(limits=c(0, 100)) +
  scale_y_continuous(limits=c(0, 100)) +
  geom_point() + 
  geom_text_repel(size=2.5, box.padding = 0.5) +
  geom_abline(intercept = 0, slope = 1, color = "black", linetype = "dotted")+
  labs(x = "% EIs captured by genetics", y = "% EIs captured by transcriptomics",
       size = "#EIs", colour="h2") +
  ggtitle("Percentage of EIs explained by transcriptomics vs. genetics")+ 
  theme(plot.title = element_text(hjust = 0.5))
g

g = ggplot(diseases_df, aes(x=Perc_genetic, y=Perc_geneexpr, label=disease_name, size=N_interactions, colour=h2_west_Sp)) +
  theme_classic() +
  scale_x_continuous(limits=c(0, 100)) +
  scale_y_continuous(limits=c(0, 100)) +
  geom_point() + 
  geom_text_repel(size=2.5, box.padding = 0.5) +
  geom_abline(intercept = 0, slope = 1, color = "black", linetype = "dotted")+
  labs(x = "% EIs captured by genetics", y = "% EIs captured by transcriptomics",
       size = "#EIs", colour="h2 Spousal") +
  ggtitle("Percentage of EIs explained by transcriptomics vs. genetics")+ 
  theme(plot.title = element_text(hjust = 0.5))

g

other_possible_colors = c("#F2AC83", "#8CDADB", "#FFA46B", "#6AFDFF", "#FFB109")

dev.off()

pdf("Plots/barplots.pdf", width = 13, height = 9)

# Keep first 50 characters of the disease names
nchars = 50
diseases_df$short_disease_name <- ifelse(nchar(diseases_df$disease_name) > nchars, paste0(substr(diseases_df$disease_name, 1, nchars), "..."), diseases_df$disease_name)

g = ggplot(diseases_df, aes(x=reorder(short_disease_name, -expr_vs_genetic), y=expr_vs_genetic, fill=N_interactions)) + 
  geom_bar(stat = "identity") +
  theme_classic() +
  labs(x = "Diseases", y = "Ratio of the % EIs captured by transcriptomics vs. genetics",
       fill="#EIs") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  geom_abline(intercept = 1, slope = 0, color = "black", linetype = "dotted") +
  theme(plot.margin = unit(c(1, 1, 1, 4), "cm"))
g

g = ggplot(diseases_df, aes(x=reorder(short_disease_name, -expr_vs_genetic), y=expr_vs_genetic, fill=Perc_geneexpr)) + 
  geom_bar(stat = "identity") +
  theme_classic() +
  labs(x = "Diseases", y = "Ratio of the % EIs captured by transcriptomics vs. genetics",
       fill="%EIs captured by transcriptomics") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  geom_abline(intercept = 1, slope = 0, color = "black", linetype = "dotted") +
  theme(plot.margin = unit(c(1, 1, 1, 4), "cm"))
  # theme(legend.position = "bottom")
g
dev.off()

pdf("Plots/summary_individual_diseases.pdf", width = 11, height = 14)
# Geom_points - originak
g = ggplot(diseases_df, aes(x=Perc_geneexpr, y=reorder(short_disease_name, Perc_geneexpr)), alpha=0.4) + 
  geom_point(aes(color = "Gene Expression"), alpha=0.6) +
  geom_point(aes(x=Perc_snps, y=reorder(short_disease_name, Perc_geneexpr), color="SNPs"), alpha=0.4, position = position_jitter(width = 0.5)) +
  geom_point(aes(x=Perc_genes, y=reorder(short_disease_name, Perc_geneexpr), color="Genes"), alpha=0.4, position = position_jitter(width = 0.5)) +
  geom_point(aes(x=Perc_ppis, y=reorder(short_disease_name, Perc_geneexpr), color="PPIs"), alpha=0.4, position = position_jitter(width = 0.5)) +
  geom_point(aes(x=Perc_pathways, y=reorder(short_disease_name, Perc_geneexpr), color="Pathways"), alpha=0.4, position = position_jitter(width = 0.5)) +
  geom_point(aes(x=Perc_geneticcorr, y=reorder(short_disease_name, Perc_geneexpr), color="Genetic Corr."), alpha=0.4, position = position_jitter(width = 0.5)) +
  scale_color_manual(values=c("black", "red", "pink", "orange", "green", "blue"), 
                     labels = c("Gene Expression", "SNPs", "Genes", "PPIs", "Pathways", "Genetic Corr.")) +
  theme_classic() +
  labs(x = "% EIs captured by the molecular layers", y = "Diseases", color="Molecular layers") +
  # theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(plot.margin = unit(c(1, 1, 1, 1), "cm")) +
  theme(legend.position = "bottom")
g

# LOLIPOP CHART - the final one (w/ correct jitter and legend)
set.seed(1)
g <- ggplot(diseases_df, aes(y = reorder(short_disease_name, Perc_geneexpr))) +
  geom_point(aes(x = Perc_geneexpr, color = "Gene Expression"), alpha = 0.65, position = position_jitter(width = 0.5)) +
  geom_segment(aes(x = 0, xend = Perc_snps, yend = reorder(short_disease_name, Perc_geneexpr), color = "SNPs"), alpha = 0.3, position = position_jitter(width = 0.5)) +
  geom_segment(aes(x = 0, xend = Perc_genes, yend = reorder(short_disease_name, Perc_geneexpr), color = "Genes"), alpha = 0.3, position = position_jitter(width = 0.5)) +
  geom_segment(aes(x = 0, xend = Perc_ppis, yend = reorder(short_disease_name, Perc_geneexpr), color = "PPIs"), alpha = 0.3, position = position_jitter(width = 0.5)) +
  geom_segment(aes(x = 0, xend = Perc_pathways, yend = reorder(short_disease_name, Perc_geneexpr), color = "Pathways"), alpha = 0.3, position = position_jitter(width = 0.5)) +
  geom_segment(aes(x = 0, xend = Perc_geneticcorr, yend = reorder(short_disease_name, Perc_geneexpr), color = "Genetic Corr."), alpha = 0.3, position = position_jitter(width = 0.5)) +
  scale_color_manual(values = c("Gene Expression" = "black", #"#00BFBF",
                                "SNPs" = "#FF5733",
                                "Genes" = "#FF851B",
                                "PPIs" = "#FFC300",
                                "Pathways" = "#0074D9",
                                "Genetic Corr." = "#9A33FF"),
                     limits = c("Gene Expression", "PPIs", "SNPs", "Pathways", "Genes","Genetic Corr.")) +
  theme_classic() +
  labs(x = "% EIs captured by the molecular layers", y = "Diseases", color = "Molecular layers") +
  theme(plot.margin = unit(c(1, 1, 1, 1), "cm")) +
  theme(legend.position = "bottom")

# Build the plot and retrieve the actual positions of the segments
plot_build <- ggplot_build(g)
segments <- plot_build$data[2:7]

# Combine the list of data frames into a single data frame
combined_segments <- do.call(rbind, segments)
combined_segments$xend <- as.numeric(as.character(combined_segments$xend))
combined_segments$yend <- as.numeric(as.character(combined_segments$yend))

set.seed(1)
# Add a point at the exact end position of each segment
p <- g +
  geom_point(data = combined_segments, aes(x = xend, y = yend), alpha = 0.6, color = combined_segments$colour)
p


dev.off()
set.seed(1)

library(fmsb)
library(tidyverse)
library(dplyr)
library(ggradar)
cdis = "F33"

to_plot_all = diseases_df[, c("icd10","Perc_geneexpr", "Perc_snps", "Perc_genes", "Perc_ppis", "Perc_pathways", "Perc_geneticcorr")]
head(to_plot_all)
# Divide columns by 100
scaled_all = to_plot_all
scaled_all[, 2:ncol(scaled_all)] <- to_plot_all[, 2:ncol(to_plot_all)] / 100
# colnames(to_plot_all) <- c("Group","Var1","Var2","Var3","Var4","Var5","Var6")
# GGRADAR (remove to remove the two last rows, which are min and max values)
ggradar(scaled_all[1:20, ])
col_mean <- apply(scaled_all[, 2:ncol(scaled_all)], 2, mean)

# Do this for disease categories

# Get average profile
head(scaled_all)
col_mean <- apply(scaled_all[, 2:ncol(scaled_all)], 2, mean)
# Put together the summary of columns
col_summary <- t(data.frame(Max = c(0, rep(1,ncol(scaled_all)-1)), 
                            Min = c(0, rep(0,ncol(scaled_all)-1)), 
                            Average = c(0,col_mean)))
colnames(col_summary)[1] = "icd10"
# Bind variables summary to the data
to_plot_all2 <- as.data.frame(rbind(col_summary, scaled_all))
head(to_plot_all2)
ggradar(to_plot_all2[3:5, 1:ncol(to_plot_all2)])


try(dev.off())
pdf("Plots/radar_charts.pdf", width = 11, height = 9)

# Restore the standard par() settings
try(par <- par(opar))
# Sort by Ratio of %EIs transcriptomics / %EIs genetics
diseases_df <- diseases_df[order(desc(diseases_df$expr_vs_genetic)), ]
for(dis in diseases_df$icd10){
  create_disease_radar_chart(dis, diseases_df, color=NULL, average_profile=FALSE)
}

dev.off()

pdf("Plots/radar_charts_disease_color.pdf", width = 7, height = 5)

# Restore the standard par() settings
try(par <- par(opar))
# Sort by Ratio of %EIs transcriptomics / %EIs genetics
diseases_df <- diseases_df[order(desc(diseases_df$expr_vs_genetic)), ]
for(dis in diseases_df$icd10){
  ccolor = diseases_df[diseases_df$icd10 == dis, ]$color_ratio
  create_disease_radar_chart(dis, diseases_df, color=ccolor, average_profile=FALSE)
}

dev.off()

pdf("Plots/radar_charts_disease_cat.pdf", width = 11, height = 9)
# Sort by Ratio of %EIs transcriptomics / %EIs genetics
# diseases_df <- diseases_df[order(desc(diseases_df$expr_vs_genetic)), ]
Ncat = as.data.frame(table(diseases_df_cat$Disease_category))
Ncat$Var1 = as.character(Ncat$Var1)

# Calcular la media de 'expr_vs_genetic' por categoría de 'Disease_category'
cat_df <- aggregate(icd10 ~ Disease_category, diseases_df_cat, length)
cat_df$cat_mean_ratio <- aggregate(expr_vs_genetic ~ Disease_category, diseases_df_cat, mean)$expr_vs_genetic
cat_df$cat_Nint <- aggregate(N_interactions ~ Disease_category, diseases_df_cat, sum)$N_interactions
cat_df$cat_mean_Nint <- aggregate(N_interactions ~ Disease_category, diseases_df_cat, mean)$N_interactions
cat_df
colnames(cat_df)[2] = "Ndis"
# cat_df = cat_df[order(cat_df$cat_mean_Nint, decreasing = TRUE), ]
cat_df = cat_df[order(cat_df$cat_mean_ratio, decreasing = TRUE), ]

try(par <- par(opar))

pfcol=NA  # options: NA / "default"
pdf(paste0("Plots/radar_charts_disease_cat_everything_",pfcol,".pdf"), width = 11, height = 10)
# Everything
for(cat in cat_df$Disease_category){
  # cat = diseases_df_cat$Disease_category[2] # to test
  cdislist = diseases_df_cat[diseases_df_cat$Disease_category == cat, ]$icd10
  cNcat = Ncat[Ncat$Var1 == cat, ]$Freq
  count1 = create_diseases_radar_chart(cdislist, diseases_df, average_profile=TRUE, selected_dis_average=TRUE, individual_dis=TRUE, title=cat, subtitle=cNcat, pfcol=pfcol)
  head(count1)
}
dev.off()
try(dev.off())
pdf(paste0("Plots/radar_charts_disease_cat_diseases_w_total_average_",pfcol,".pdf"), width = 11, height = 10)

# Individual diseases to total average
for(cat in cat_df$Disease_category){
  # cat = diseases_df_cat$Disease_category[2] # to test
  cdislist = diseases_df_cat[diseases_df_cat$Disease_category == cat, ]$icd10
  cNcat = Ncat[Ncat$Var1 == cat, ]$Freq
  count3 = create_diseases_radar_chart(cdislist, diseases_df, average_profile=TRUE, selected_dis_average=FALSE, individual_dis=TRUE, title=cat, subtitle=cNcat, pfcol=pfcol)
  head(count3)
}
dev.off()

pdf("Plots/radar_charts_disease_cataverages_w_total_average_one_page.pdf", width = 8.27, height = 11.69) # dimensions of an A4
# Averages with reference
# to save plot
op <- par(mar = c(1, 1, 1, 1))
par(mfrow = c(6,3))
par(mai = c(0.5, 0.5, 0.5, 0.5))
for(cat in cat_df$Disease_category){
  # cat = diseases_df_cat$Disease_category[2] # to test
  cdislist = diseases_df_cat[diseases_df_cat$Disease_category == cat, ]$icd10
  cNcat = Ncat[Ncat$Var1 == cat, ]$Freq
  title = paste0(cat," (",cat_df[cat_df$Disease_category == cat, ]$cat_Nint,")")
  create_diseases_radar_chart(cdislist, diseases_df, average_profile=TRUE, selected_dis_average=TRUE, individual_dis=FALSE, title=title, subtitle=cNcat, legend=FALSE)
}
par(op)
dev.off()

pdf("Plots/radar_charts_disease_cataverages_w_total_average.pdf", width = 11, height = 9)
# to save text
category_averages = data.frame()
k = 0
for(cat in cat_df$Disease_category){
  print(cat)
  # cat = diseases_df_cat$Disease_category[2] # to test
  cdislist = diseases_df_cat[diseases_df_cat$Disease_category == cat, ]$icd10
  cNcat = Ncat[Ncat$Var1 == cat, ]$Freq
  title = paste0(cat," (NEIs=",cat_df[cat_df$Disease_category == cat, ]$cat_Nint,")")
  count2 = create_diseases_radar_chart(cdislist, diseases_df, average_profile=TRUE, selected_dis_average=TRUE, individual_dis=FALSE, title=title, subtitle=cNcat)
  
  # Create table
  head(count2)
  if (k == 0){
    category_averages = count2
    category_averages$category = c("Max", "Min", "TAverage", cat)
    k = k+1
  }else{
    tobind = count2[rownames(count2) == "Average", ]; tobind$category = cat
    category_averages = rbind(category_averages, tobind)
  }
}
dev.off()

pdf("Plots/Dis_category_description.pdf", width = 8, height = 10)

g1 = ggplot(cat_df, aes(x=reorder(Disease_category, cat_mean_Nint), y=Ndis, fill=Ndis)) + 
  geom_bar(stat="identity") + 
  geom_text(aes(label = Ndis, y = Ndis + 0.1), size = 3) +
  scale_fill_gradient(low = "#FBB4AE", high = "#B3CDE3") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(legend.position = "none")+
  xlab("Disease Category") +
  ylab("Number of diseases")+
  labs(fill="N diseases")
# g1

g2 = ggplot(cat_df, aes(x=reorder(Disease_category, cat_mean_Nint), y=cat_Nint, fill=cat_Nint)) + 
  geom_bar(stat="identity") + 
  geom_text(aes(label = cat_Nint, y = cat_Nint + 0.1), size = 3) +
  scale_fill_gradient(low = "#FBB4AE", high = "#B3CDE3") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(legend.position = "none")+
  xlab("") +
  ylab("Total number of interactions")+
  labs(fill="N diseases")
# g2

g3 = ggplot(diseases_df_cat, aes(x=reorder(Disease_category, N_interactions, FUN=mean), y=N_interactions, fill=N_interactions)) + 
  # geom_point(position = position_jitter(width = 0.2), size = 2)+
  # geom_violin() + 
  geom_boxplot(alpha=0.5) + 
  geom_point(position=position_jitterdodge(jitter.width=0, jitter.height=0.6), alpha=0.7)+
  # geom_text(aes(label = Freq, y = Freq + 0.1), size = 3) +
  # scale_fill_gradient(low = "#FBB4AE", high = "#B3CDE3") +
  # scale_color_gradient(low = "#FBB4AE", high = "#B3CDE3") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none") +
  xlab("") +
  ylab("Number of interactions")
# g3


# Combine the two plots vertically using grid.arrange from gridExtra package
library(gridExtra)
grid.arrange(g2,g3, g1, nrow = 3, heights = c(1, 0.9, 0.6))

dev.off()

write.table(diseases_df_cat, file = "Results/icd10_intersect_dong_df.txt", sep = "\t", row.names = FALSE, quote=FALSE)


# create a sample data.frame with 5 features and a variable
set.seed(123)
my_df <- data.frame(
  feature1 = sample(LETTERS[1:5], 50, replace = TRUE),
  feature2 = sample(LETTERS[1:5], 50, replace = TRUE),
  feature3 = sample(LETTERS[1:5], 50, replace = TRUE),
  feature4 = sample(LETTERS[1:5], 50, replace = TRUE),
  feature5 = sample(LETTERS[1:5], 50, replace = TRUE),
  value = rnorm(50)
)


library(tidyr)
my_df_long <- gather(diseases_df, key = "feature", value = "icd10", -icd10)


#### MUTUAL INFORMATION
library(infotheo) # mutinformation
library(aricode) # NMI
mutinformation(df$snps, df$snps, method="emp")
mutinformation(df$snps, df$genes, method="emp")

cicds = unique(union(df$Dis1, df$Dis2))
length(cicds)
length(common_icds)
df = df[(df$Dis1 %in% common_icds) & (df$Dis2 %in% common_icds), ]
df$geneexpr = as.numeric(df$geneexpr)

MI_all = mutinformation(df[, c(4:9)])
MI_common = mutinformation(dfcommon[, c(4:9)])

H_all = diag(MI_all); H_all
H_common = diag(MI_common); H_common

nMI_all = MI_all
nMI_common = MI_common
for(k in 1:length(H_all)){
  print(k)
  nMI_all[k, ] = MI_all[k, ] / H_all[k]
  nMI_common[k, ] = MI_common[k, ] / H_common[k]
}

nMI_all*100
nMI_common*100

library(lsa)
cosine(as.matrix(df[, c(4:9)]))
cosine(as.matrix(dfcommon[, c(4:9)]))

cosine(df$snps, df$snps)
cosine(df$snps, df$genes)

head(df)
minf = matrix(data=NA, nrow=6, ncol=6)
mjaccard =  matrix(data=NA, nrow=6, ncol=6)
mintersect =  matrix(data=NA, nrow=6, ncol=6)
for(k1 in 4:9){
  for(k2 in 4:9){
    # print(paste(k1, k2))
    minf[k1-3, k2-3] = NMI(df[,k1], df[,k2])
    mjaccard[k1-3, k2-3] = jaccard_similarity(pairs_list[[k1-3]], pairs_list[[k2-3]])
    mintersect[k1-3, k2-3] = length(intersect(pairs_list[[k1-3]], pairs_list[[k2-3]]))
  }
}

minf2 = matrix(data=NA, nrow=6, ncol=6)
mjaccard2 =  matrix(data=NA, nrow=6, ncol=6)
mintersect2 =  matrix(data=NA, nrow=6, ncol=6)
for(k1 in 4:9){
  for(k2 in 4:9){
    # print(paste(k1, k2))
    minf2[k1-3, k2-3] = NMI(dfcommon[,k1], dfcommon[,k2])
    mjaccard2[k1-3, k2-3] = jaccard_similarity(pairs_list_common[[k1-3]], pairs_list_common[[k2-3]])
    # print(jaccard_similarity(pairs_list_common[[k1-3]], pairs_list_common[[k2-3]]))
    mintersect2[k1-3, k2-3] = length(intersect(pairs_list_common[[k1-3]], pairs_list_common[[k2-3]]))
  }
}

colnames(minf) = colnames(df)[4:9]
rownames(minf) = colnames(df)[4:9]

colnames(minf2) = colnames(df)[4:9]
rownames(minf2) = colnames(df)[4:9]

colnames(mjaccard) = colnames(df)[4:9]
rownames(mjaccard) = colnames(df)[4:9]

colnames(mjaccard2) = colnames(df)[4:9]
rownames(mjaccard2) = colnames(df)[4:9]

colnames(mintersect) = colnames(df)[4:9]
rownames(mintersect) = colnames(df)[4:9]

colnames(mintersect2) = colnames(df)[4:9]
rownames(mintersect2) = colnames(df)[4:9]

minf*100
minf2*100
mjaccard
mjaccard2

# PLOTTING
require(reshape)

# MI
plot_comparison_molecular_layers <- function(df1, df2, title, legend_name, by100=TRUE, remove_last=TRUE, upper_triang = FALSE){
  # To test
  # df1 = mjaccard
  # df2 = mjaccard
  # legend_name = "test"
  
  
  if(by100){
    df1 = df1*100
    df2 = df2*100
  }
  
  if(remove_last){
    df1 = df1[1:5, 1:5]
    names1 = c("snps", "genes", "ppis", "pathways", "geneticcorr")
  }else{
    names1 = c("snps", "genes", "ppis", "pathways", "geneticcorr", "geneexpr")
  }
  
  if(upper_triang){
    df1[lower.tri(df1)] = NA
    df2[lower.tri(df2)] = NA
  }
  
  df1m = na.omit(melt(t(df1))); head(df1m)
  df2m = na.omit(melt(t(df2))); head(df2m)
  
  df1m$X1 = factor(df1m$X1, ordered=TRUE, levels = names1) 
  df1m$X2 = factor(df1m$X2, ordered=TRUE, levels = rev(names1))
  
  df2m$X1 = factor(df2m$X1, ordered=TRUE, levels = c("snps", "genes", "ppis", "pathways", "geneticcorr", "geneexpr")) 
  df2m$X2 = factor(df2m$X2, ordered=TRUE, levels = c("geneexpr", "geneticcorr", "pathways", "ppis", "genes", "snps"))
  
  max_color =  "#400A80"  #"#730F73" "#781078"  "#4F0A4F" # steelblue "darkblue"
  
  g = ggplot(df1m, aes(X1, X2, fill=value)) + geom_tile() +
    theme_minimal() + theme(panel.grid.major = element_line(colour = "white"))+  # axis.ticks = element_line(colour="grey50")
    geom_text(aes(label = round(value,2)), color="black") +
    scale_fill_gradient(low = "white", high = max_color) +
    labs(fill=legend_name, title = title) + ylab("") +
    xlab("Y") + ylab("X") +
  theme(plot.title = element_text(hjust = 0.5), axis.text=element_text(size=11))
  print(g)
  
  g = ggplot(df2m, aes(X1, X2, fill=value)) + geom_tile() +
    theme_minimal() + theme(panel.grid.major = element_line(colour = "white"))+
    geom_text(aes(label = round(value,2)), color="black") +
    scale_fill_gradient(low = "white", high = max_color) +
    labs(fill=legend_name, title = title) + ylab("") +
    xlab("Y") + ylab("X") +
  theme(plot.title = element_text(hjust = 0.5), axis.text=element_text(size=11))
  print(g)
  
}

plot_comparison_molecular_layers(nMI_all, nMI_common, "I(X;Y) / H(X)", "I(X;Y) / H(X)")
plot_comparison_molecular_layers(minf, minf2, "Normalized mutual information", "Normalized I", upper_triang = TRUE)
plot_comparison_molecular_layers(mjaccard, mjaccard2, "Jaccard Index", "Jaccard Index", remove_last = FALSE, upper_triang = TRUE)
plot_comparison_molecular_layers(mjaccard, mjaccard2, "Jaccard Index", "Jaccard Index", upper_triang = TRUE)
plot_comparison_molecular_layers(mintersect, mintersect2, "Common Epidemiological Interactions (EI)", "Common EI", by100 = FALSE, remove_last = FALSE, upper_triang = TRUE)
plot_comparison_molecular_layers(mintersect, mintersect2, "Common Epidemiological Interactions (EI)", "Common EI", by100 = FALSE, upper_triang = TRUE)
mrecall = (mintersect / diag(mintersect))*100
mrecall2 = (mintersect2 / diag(mintersect2))*100
plot_comparison_molecular_layers(mrecall, mrecall2, "Recall of Y from X", "Recall", by100 = FALSE, remove_last = FALSE)


# Plot of Dong + 5 levels in Dong
listInput = list(Dong=dong_epidem_all$pairs, 
                 SNPs=snp_df$pairs,
                 Genes=gene_df$pairs,
                 Pathways=pathways_df$pairs,
                 PPIs=ppi_df$pairs,
                 Genetic_Correlation=geneticcorr_df$pairs)

# Euler diagram
eu = euler(listInput, shape="ellipse")
plot(eu, quantities=TRUE, main="Entire Dong et al. and molecular levels")

# Upset plot
upset(fromList(listInput), nsets=6, order.by = "freq")

# Another Eurler
eu = euler(list(Pathways=pathways_df$pairs,
                PPIs=ppi_df$pairs,
                Dong=dong_epidem_all$pairs), 
           shape="ellipse")
plot(eu, quantities=TRUE,
     main="Notice that the previous Euler diagram was wrong")


# Plot of Dong + 5 levels in Dong + DSN
listInput = list(Dong=dong_epidem_all$pairs, 
                 DSN=dsn$pairs,
                 SNPs=snp_df$pairs,
                 Genes=gene_df$pairs,
                 Pathways=pathways_df$pairs,
                 PPIs=ppi_df$pairs,
                 Genetic_Correlation=geneticcorr_df$pairs)

# Euler
eu = euler(listInput, shape="ellipse")
plot(eu, quantities=TRUE, main="Entire Dong, molecular levels and DSN") # no sale bien el plot

# Upset
upset(fromList(listInput), nsets=7, order.by = "freq")
# upset(fromList(listInput), nsets=7, order.by = "freq", keep.order = TRUE, empty.intersections = "on")

# Euler Dong & DSN only
eu = euler(list(dsn=dsn$pairs, Dong=dong_epidem_all$pairs), shape="ellipse")
plot(eu, quantities=TRUE, main="Dong and DSN overlap")

### UpSet plot with Dong et al entailing common diseases
dong_epidem$pairs = paste(dong_epidem$Dis1, dong_epidem$Dis2, sep="_")

## Filtering for common diseases in the rest
cgene_df = gene_df[(gene_df$Dis1 %in% common_icds) & (gene_df$Dis2 %in% common_icds), ]
cgeneticcorr_df = geneticcorr_df[(geneticcorr_df$Dis1 %in% common_icds) & (geneticcorr_df$Dis2 %in% common_icds), ]
cpathways_df = pathways_df[(pathways_df$Dis1 %in% common_icds) & (pathways_df$Dis2 %in% common_icds), ]
cppi_df = ppi_df[(ppi_df$Dis1 %in% common_icds) & (ppi_df$Dis2 %in% common_icds), ]
csnp_df = snp_df[(snp_df$Dis1 %in% common_icds) & (snp_df$Dis2 %in% common_icds), ]

length(unique(cgene_df$pairs)) 
length(unique(cgeneticcorr_df$pairs))
length(unique(cpathways_df$pairs))
length(unique(cppi_df$pairs))
length(unique(csnp_df$pairs))

####### Common icds: Genetic, Dong and Expression
listInput = list(Dong=dong_epidem$pairs, 
                 Genetic=unique(c(csnp_df$pairs,cgene_df$pairs,cpathways_df$pairs, cppi_df$pairs, cgeneticcorr_df$pairs)),
                 Expression=intersect(dsn$pairs, dong_epidem$pairs))

# Euler diagram
eu = euler(listInput, shape="ellipse")
plot(eu, quantities=TRUE, main="Entire Dong et al. and molecular levels")

# Upset plot
upset(fromList(listInput), nsets=6, order.by = "freq")

listInput = list(Genetic=unique(c(csnp_df$pairs,cgene_df$pairs,cpathways_df$pairs, cppi_df$pairs, cgeneticcorr_df$pairs)),
                 Expression=intersect(dsn$pairs, dong_epidem$pairs))

# Euler diagram
eu = euler(listInput, shape="ellipse")
plot(eu, quantities=TRUE, main="Entire Dong et al. and molecular levels")

# Upset plot
upset(fromList(listInput), nsets=6, order.by = "freq")

### BARPLOTS WITH THE RECALLS
df_numbers = c(length(intersect(unique(cgene_df$pairs), dong_epidem$pairs)),
              length(intersect(unique(cgeneticcorr_df$pairs), dong_epidem$pairs)),
              length(intersect(unique(cpathways_df$pairs), dong_epidem$pairs)),
              length(intersect(unique(cppi_df$pairs), dong_epidem$pairs)),
              length(intersect(unique(csnp_df$pairs), dong_epidem$pairs)),
              length(intersect(unique(dsn$pairs), dong_epidem$pairs))
              )
df_names <- c("Genes","Genetic Correlation","Pathways",
                         "PPIs", "SNPs", "Gene Expression")
barplotdf = data.frame("overlap" = df_numbers, "layer" = as.character(df_names))
barplotdf$perc_overlap = (barplotdf$overlap / length(unique(dong_epidem$pairs)))*100

barp1 = ggplot(barplotdf, aes(reorder(layer, overlap), overlap, fill=overlap)) + 
  geom_bar(stat="identity") +
  geom_text(aes(label=overlap), vjust=-0.5)+
  theme_classic()+
  theme(plot.title = element_text(hjust = 0.5), axis.text=element_text(size=11))+
  ggtitle("Number of interactions explained by each molecular information")+
  xlab("Molecular information")+
  ylab("Number of interactions")
# barp1

barp2 = ggplot(barplotdf, aes(reorder(layer, perc_overlap), perc_overlap, fill=overlap)) + 
  geom_bar(stat="identity") +
  geom_text(aes(label=round(perc_overlap, digits = 2)), vjust=-0.5)+
  scale_fill_gradient(low="#174445", high="#8CDADB")+
  ggtitle("Percentage of interactions from Dong et al. explained by each molecular information")+
  xlab("Molecular information")+
  ylab("Percentage of interactions")+
  theme_classic()+
  theme(plot.title = element_text(hjust = 0.5), axis.text=element_text(size=11))
# barp2

barp2B = ggplot(barplotdf, aes(reorder(layer, perc_overlap), perc_overlap, fill=overlap)) + 
  geom_bar(stat="identity") +
  geom_text(aes(label=round(perc_overlap, digits = 2)), vjust=-0.5, size=4)+
  scale_x_discrete(
    labels = c(
      "SNPs" = "SNPs",
      "PPIs" = "PPIs",
      "Genes" = "Genes",
      "Genetic Correlation" = "Genetic\ncorrelation",
      "Pathways" = "Pathways",
      "Gene Expression" = "Gene\nexpression"
    )
  )+
  scale_fill_gradient(low="#174445", high="#8CDADB")+
  # ggtitle("Percentage of interactions from Dong et al. explained by each molecular information")+
  xlab("Molecular information")+
  ylab("Percentage of interactions")+
  theme_classic()+
  theme(plot.title = element_text(hjust = 0.5), 
        # Axis titles
        axis.title.x = element_text(
          size = 15,
          margin = margin(t = 7)
        ),
        
        axis.title.y = element_text(
          size = 15,
          margin = margin(r = 7)
        ),
        
        # Category names
        axis.text.x = element_text(
          size = 14,
          lineheight = 0.85
        ),
        
        # Y-axis numbers
        axis.text.y = element_text(
          size = 12
        ),
        legend.position = "none")


write("\nKeeping only common ICD10s:", file=outputname, append=TRUE)
write(paste("Dim Genes:", dim(cgene_df)[1]), file=outputname, append=TRUE)
write(paste("Dim Genetic Corr:", dim(cgeneticcorr_df)[1]), file=outputname, append=TRUE)
write(paste("Dim Pathways:", dim(cpathways_df)[1]), file=outputname, append=TRUE)
write(paste("Dim PPIs:", dim(cppi_df)[1]), file=outputname, append=TRUE)
write(paste("Dim SNPs:", dim(csnp_df)[1]), file=outputname, append=TRUE)

listInput = list(Dong=dong_epidem$pairs, 
                 DSN=dsn$pairs,
                 SNPs=csnp_df$pairs,
                 Genes=cgene_df$pairs,
                 Pathways=cpathways_df$pairs,
                 PPIs=cppi_df$pairs,
                 Genetic_Correlation=cgeneticcorr_df$pairs)

# Euler
eu = euler(listInput, shape="ellipse")
plot(eu, quantities=TRUE, main="Dong, molecular levels and DSN over common ICD10s")

# Upset
upset(fromList(listInput), nsets=7, order.by = "freq")
# upset(fromList(listInput), nsets=7, order.by = "freq", keep.order = TRUE, empty.intersections = "on")
upset(fromList(listInput), nsets=7, order.by = "freq", group.by = "sets")



## From the interactions in Dong et al. 
listInput = list(Dong=dong_epidem$pairs, 
                  DSN=intersect(dsn$pairs, dong_epidem$pairs),
                  SNPs=csnp_df$pairs,
                  Genes=cgene_df$pairs,
                  Pathways=cpathways_df$pairs,
                  PPIs=cppi_df$pairs,
                  'Genetic Correlation'=cgeneticcorr_df$pairs)

# Euler
eu = euler(listInput, shape="ellipse")
plot(eu, quantities=TRUE, main="Dong interactions over the common ICD10s")

# Upset
upset(fromList(listInput), nsets=7, order.by = "freq")
# upset(fromList(listInput), nsets=7, order.by = "freq", keep.order = TRUE, empty.intersections = "on")
upset(fromList(listInput), nsets=7, order.by = "freq", group.by = "sets")

# Colores Upset plot
#low="#174445", high="#8CDADB"
###### FIGURE 1D
try(dev.off())
outputname_plots = paste0("Plots/Final_Upset_plot","_",network_type,"_splitduplicates_",split_duplicates,"_only_interprt_diseases_",only_interpretable_diseases,"_",filt,"_",sign_threshold,"_2",".pdf")
pdf(file=outputname_plots, width = 8, height = 7)
upset(fromList(listInput), nsets=7, order.by = "freq", text.scale = c(1.5, 1.4, 1.5, 1.4, 1.6, 1.1),
      # text.scale = 1.5,
      queries = list(list(query = intersects, params = list("Dong"), color = "#F2AC83", active = T), 
                     list(query = intersects, params = list("DSN", "Dong"), color = "#8CDADB", active = T))
)
upset(fromList(listInput), nsets=7, order.by = "freq", text.scale = c(1.5, 1.4, 1.5, 1.4, 1.6, 1.1), show.numbers = "no",
      # text.scale = 1.5,
      queries = list(list(query = intersects, params = list("Dong"), color = "#F2AC83", active = T), 
                     list(query = intersects, params = list("DSN", "Dong"), color = "#8CDADB", active = T))
      )
try(dev.off())


listInput1 = listInput

union_dong_molecular = unique(c(csnp_df$pairs, cgene_df$pairs, cpathways_df$pairs, cppi_df$pairs, cgeneticcorr_df$pairs))
length(union_dong_molecular)
union_dong_molecular_in_dong = intersect(union_dong_molecular, dong_epidem$pairs); length(union_dong_molecular_in_dong) # All in Dong
write(paste("\nInteractions in Dong explained by Dong (at least one of the layers):", length(union_dong_molecular_in_dong)), file=outputname, append=TRUE)

length(dong_epidem$pairs) # 117
length(unique(dong_epidem$pairs)) # 117
write(paste("\nUnique interactions in Dong:", length(unique(dong_epidem$pairs))), file=outputname, append=TRUE)
# Recapitulamos 55/117 = 47% de las interacciones, de las cuales 21 únicamente se explican en la DSN y 34 tienen alguna otra explicación molecular

# OVERLAP between Molecular Dong and DSN
# Dong molecular
union_dong_molecular_dsn =  intersect(union_dong_molecular, dsn$pairs); length(union_dong_molecular_dsn)
write(paste("\nInteractions in Dong molecular:", length(union_dong_molecular)), file=outputname, append=TRUE)
write(paste("Interactions in Dong molecular in DSN (overlap):", length(union_dong_molecular_dsn)), file=outputname, append=TRUE)
write(paste("Interactions in Dong molecular in DSN (precision):", length(union_dong_molecular_dsn)/length(dsn$pairs)), file=outputname, append=TRUE)
write(paste("Interactions in Dong molecular in DSN (recall):", length(union_dong_molecular_dsn)/length(union_dong_molecular)), file=outputname, append=TRUE)
# Dong not molecular
dong_not_molecular = setdiff(dong_epidem$pairs, union_dong_molecular)
dong_not_molecular_dsn =  intersect(dong_not_molecular, dsn$pairs); length(dong_not_molecular_dsn)
write(paste("\nInteractions in Dong NOT molecular:", length(dong_not_molecular)), file=outputname, append=TRUE)
write(paste("Interactions in Dong NOT molecular in DSN (overlap):", length(dong_not_molecular_dsn)), file=outputname, append=TRUE)
write(paste("Interactions in Dong NOT molecular in DSN (precision):", length(dong_not_molecular_dsn)/length(dsn$pairs)), file=outputname, append=TRUE)
write(paste("Interactions in Dong NOT molecular in DSN (recall):", length(dong_not_molecular_dsn)/length(dong_not_molecular)), file=outputname, append=TRUE)

# See how they overlap with the DSN (include it in the Venn Diagram / Euler diagram)
common_int = intersect(dong_epidem$pairs, dsn$pairs); length(common_int) # 55
write(paste("\nInteractions in Dong explained by the DSN:", length(common_int)), file=outputname, append=TRUE)
listInput = list(DSN=common_int,
                 SNPs=intersect(csnp_df$pairs, common_int),
                 Genes=intersect(cgene_df$pairs, common_int),
                 Pathways=intersect(cpathways_df$pairs, common_int),
                 PPIs=intersect(cppi_df$pairs, common_int),
                 'Genetic Correlation'=intersect(cgeneticcorr_df$pairs, common_int))

# Euler
eu = euler(listInput, shape="ellipse")
plot(eu, quantities=TRUE, main="Over the overlap between Dong and DSN")

# Upset
upset(fromList(listInput), nsets=7, order.by = "freq")
listInput2 = listInput

upset(fromList(listInput), nsets=7, order.by = "freq",
      queries = list(list(query = intersects, params = list("DSN"), color = "#8CDADB", active = T))
)

# SAVE BARPLOTS
try(dev.off())
outputname_plots = paste0("Plots/Upset_plots_and_euler_diagrams","_",network_type,"_splitduplicates_",split_duplicates,"_only_interprt_diseases_",only_interpretable_diseases,"_",filt,"_",sign_threshold,"_2",".pdf")
pdf(file=outputname_plots, width = 8, height = 9)  # Old  width = 8, height = 9

barp1
barp2
barp2B

try(dev.off())
