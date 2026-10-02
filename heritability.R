#################################################################################################
# Analysis on heritability and the genetic, non-genetic origin of disease co-occurrences
#
# Author: Beatriz Urda García, 2024
##################################################################################################

library(ggrepel)
library(ggpubr)
library(lsa)
library(gridExtra)
library(ggExtra)
source("ukb_comparison_tools.R")

setwd('~/Desktop/ANALYSIS/Genome_medicine/')

remove_dis=FALSE  # FALSE, I consider all diseases
split_dup = FALSE # FALSE, 
# The results are available for all combinations of the following parameters: Min_n_int and with_icds_synonyms.
min_n_int = 0    # 3 and 0. 0 by defect and 0 to compare with
with_icds_synonyms = FALSE # FALSE, very consistent results than with TRUE 

# OPEN FILE WITH EIs
if(split_dup){
  EIs_path = 'Results/dong_interactions_in_molecular_networks_new_ssn_0.01_split_dup_TRUE_only_interp_dis_TRUE.txt' # OLD
}else{
  EIs_path = 'Results/dong_interactions_in_molecular_networks_new_ssn_0.01_split_dup_FALSE_only_interp_dis_TRUE.txt'  # NEW
  EIs_path_disease_level = 'Results/dong_interactions_in_molecular_networks_new_dsn_0.01_split_dup_FALSE_only_interp_dis_TRUE.txt'  # NEW
}

# OUTPUT FILENAMES
filename_sufix = paste0("_Split_dup_",substr(split_dup, 1, 1),"__","Min_n_int_",min_n_int,"__Remove_dis_",substr(remove_dis, 1, 1))
if(with_icds_synonyms){
  outdir="Results/with_icds_synonyms/"
  filename_sufix = paste0(filename_sufix,paste0("__with_ICD10_synonyms"))
}else{
  outdir="Results/wo_icds_synonyms/"
}

# considering SSN
eis = read.csv2(EIs_path,stringsAsFactors = F,sep="\t",header=T); dim(eis) # 476  13
# Add a new column 'genetic' that is 1 if at least one of the specified columns is 1
eis$genetic <- apply(eis[, c("snps", "genes", "ppis", "pathways", "geneticcorr")], 1, function(row) ifelse(any(row == 1), 1, 0))
eis$molecular <- apply(eis[, c("snps", "genes", "ppis", "pathways", "geneticcorr", "geneexpr")], 1, function(row) ifelse(any(row == 1), 1, 0))

# considering DSN
eis_disease_level = read.csv2(EIs_path_disease_level,stringsAsFactors = F,sep="\t",header=T); dim(eis) # 
# Add a new column 'genetic' that is 1 if at least one of the specified columns is 1
eis_disease_level$genetic <- apply(eis_disease_level[, c("snps", "genes", "ppis", "pathways", "geneticcorr")], 1, function(row) ifelse(any(row == 1), 1, 0))
eis_disease_level$molecular <- apply(eis_disease_level[, c("snps", "genes", "ppis", "pathways", "geneticcorr", "geneexpr")], 1, function(row) ifelse(any(row == 1), 1, 0))

icd_of_interest = "J45"
cdf = eis_disease_level[eis_disease_level$Dis1 == icd_of_interest | eis_disease_level$Dis2 == icd_of_interest, ]
c_interactions = nrow(cdf); c_interactions
tt = table(cdf$geneexpr)
head(tt)
(tt/c_interactions)*100


# SSN
ssn_filename = "../../GEVariability/Network_building/Defined_networks/all_diseases/metapatients_and_disease/metap_dis_pairwise_union_cosine_distance_sDEGs_network_1FDR.txt"
# ssn = read.delim(ssn_filename, sep = "\t", header = TRUE)

# Get dis representation of SSN with DM and DD.
raw_ssn = read.csv2(ssn_filename,stringsAsFactors = F,sep="\t",header=T)
ssn = get_dis_representation_of_ssn(ssn_filename, int="DM-DD"); dim(ssn) #  35100     5
head(ssn)
nrow(ssn)

ssn$Distance = as.numeric(ssn$Distance); ssn$pvalue = as.numeric(ssn$pvalue); ssn$adj_pvalue = as.numeric(ssn$adj_pvalue);
str(ssn)

ssn = sort_icd_interactions(ssn)
head(ssn)
nrow(ssn)

ssn$pairs = paste(ssn$Dis1, ssn$Dis2, sep="_")

# Remove M-M interactions between the same diseases
to_remove = which(ssn$Dis1 == ssn$Dis2); length(to_remove)
ssn = ssn[-to_remove, ]; nrow(ssn)

# Keep the positive interactions
ssn = ssn[ssn$Distance >= 0, ]

top_ssn = data.frame()
for(pair in unique(ssn$pairs)){
  # pair = ssn$pairs[2] # to test
  
  cdf = ssn[ssn$pairs == pair, ]
  if(nrow(cdf) == 1){
    top_ssn = rbind(top_ssn, cdf)
  }else{
    to_append = cdf[which.max(cdf$Distance), ]
    top_ssn = rbind(top_ssn, to_append)
  }
}

head(top_ssn)

fssn = network_from_dis_to_icd10(top_ssn); dim(fssn) # 7157    8


#### HERITABILITY
# Read the tab-delimited text file
df <- read.delim("Results/icd10_intersect_dong_df.txt", sep = "\t", header = TRUE); head(df)
icds = unique(df$icd10); length(icds)
ini_hdf<- read.delim("Paper_and_data/heritability.csv", sep = "\t", header = TRUE)
colnames(ini_hdf) = c("icd10", "disease_name", "h2_Agg","h2_Agg_catch","h2_AE", "rango", "h2_ACE", "match", "Literature", "Papers", "Notes")


# write.table(ini_hdf, "Paper_and_data/heritability_v0.csv", row.names = FALSE, sep="\t")
hdf = ini_hdf
hdf$disease_name = NULL

nrow(df) # 58
nrow(hdf) # 56
length(unique(hdf$icd10))
# 35

df = df[!duplicated(df$icd10),]

head(df)

mdf = merge(df, hdf, by="icd10", all.x=TRUE, all.y=FALSE)
nrow(mdf)
# mdf = 

# Add the heritability by Westergaard et al. 
west <- read.delim("Paper_and_data/website_icd10_lvl3.tsv", sep = "\t", header = TRUE)
head(west)
colnames(west) = c("icd10", "term","disease_name", "chapter","chapter_name","h2_west", "C1_west", "C2_west")
west0 = west

###############MEJORAR PARA EL TOTAL DEL UKB DE DONG""
if(with_icds_synonyms){ # Considering heritabilities of ICD10s with synonyms in Dong et al. 
  icds_with_synonyms = find_icds_with_synonyms(mdf)
  print(icds_with_synonyms)
  dim(west)
  west = add_average_heritability_to_west_df(west, icds_with_synonyms)
  dim(west) # 15 more than before, perfect!
}

westA = west[west$term == "A", ]
westSp = west[west$term == "Sp", ]
westSp = westSp[, c("icd10","h2_west")]
colnames(westSp) = c("icd10", "h2_west_Sp")
nrow(mdf)

mdf = merge(mdf, westA[, c("icd10","h2_west", "C1_west", "C2_west","chapter","chapter_name")], by="icd10", all.x=TRUE, all.y=FALSE)
nrow(mdf)
mdf = merge(mdf, westSp[, c("icd10","h2_west_Sp")], by="icd10", all.x=TRUE, all.y=FALSE)
nrow(mdf)
mdf$h2_catch_west = rowMeans(mdf[, c("h2_Agg_catch", "h2_west")], na.rm=TRUE)
mdf$match = as.numeric(mdf$match)

# ALL DONG (all original ICD10s in the UKB)
head(ceis_like)
# Add a column with 'genetic'
ceis_like <- ceis_like %>%
  mutate(genetic = as.integer(if_any(c(snps, genes, ppis, pathways, geneticcorr), ~ . == 1)))

all_diseases = union(ceis_like$Dis1, ceis_like$Dis2); length(all_diseases); # 310
vertices_df2 = data.frame(icd10 = all_diseases, 
                          stringsAsFactors = FALSE)
vertices_df2 = merge(vertices_df2, westA[, c("icd10","h2_west", "C1_west", "C2_west","disease_name","chapter","chapter_name")], by="icd10", all.x=TRUE, all.y=FALSE)
nrow(vertices_df2)
vertices_df2 = merge(vertices_df2, westSp[, c("icd10","h2_west_Sp")], by="icd10", all.x=TRUE, all.y=FALSE)
nrow(vertices_df2) # 310

###### ASSORTATIVITY AND HERITABILITY ######
# Compute assortativity of EIs in the molecular layers based on heritability.
# Get ICD10s in ssn and west
icd10_in_ssn = unique(union(fssn$Dis1, fssn$Dis2)); length(icd10_in_ssn)
icd10_in_ssn_west = intersect(icd10_in_ssn, westA$icd10); length(icd10_in_ssn_west)
cfssn = fssn[(fssn$Dis1 %in% icd10_in_ssn_west) & (fssn$Dis2 %in% icd10_in_ssn_west), ]
nrow(eis)
ceis = eis[(eis$Dis1 %in% icd10_in_ssn_west) & (eis$Dis2 %in% icd10_in_ssn_west), ]; nrow(ceis)
nrow(fssn)
nrow(cfssn)

pdf(paste0(outdir,"assortativity_h2_",filename_sufix,".pdf"))

# Initialize an empty data frame to store results
assortativity_h2 <- data.frame(network = character(), layer = character(), assortativity = numeric(), stringsAsFactors = FALSE)

vertices_df = westA[westA$icd10 %in% icd10_in_ssn_west, ]
vertices_df$h2_log = log10(vertices_df$h2_west)
# vertices_df2$h2_log = log10(vertices_df2$h2_west)

c_asortativity = compute_heritability_assortativity(cfssn, vertices_df, attribute_name = "h2_west") # -0.01485
c_asortativity_log10 = compute_heritability_assortativity(cfssn, vertices_df, attribute_name = "h2_log") # -0.01587975


vertices_df = westA[westA$icd10 %in% union(eis$Dis1, eis$Dis2), ]
assortativity_h2 <- rbind(assortativity_h2, data.frame(network = "cfssn", layer = "", assortativity = c_asortativity))

# For the entire Dong network.
c_asortativity = compute_heritability_assortativity(edges_df=ceis[, c("Dis1", "Dis2")],
                                                    vertices_df, 
                                                    attribute_name = "h2_west") # -0.01485
# Save the result in the results_df
assortativity_h2 <- rbind(assortativity_h2, data.frame(network = "ceis", layer = "UKB EIs", assortativity = c_asortativity))

for(layer in c("snps", "genes", "ppis", "pathways", "geneticcorr", "geneexpr", "genetic", "molecular")){
  # layer = "snps"
  print(layer)
  edges_df = ceis[, c("Dis1", "Dis2", layer)]; edges_df = edges_df[edges_df[[layer]] == 1, ]
  c_asortativity = compute_heritability_assortativity(edges_df, vertices_df, attribute_name = "h2_west") # -0.01485
  # Save the result in the results_df
  assortativity_h2 <- rbind(assortativity_h2, data.frame(network = "ceis", layer = layer, assortativity = c_asortativity))
}
write.csv(assortativity_h2, file = paste0(outdir,"assortativity_h2_",filename_sufix,".csv"), row.names = FALSE)


# REPEAT THE SAME BUT FOR ALL DONG (ALL DISEASES IN THE UKB)
# Initialize an empty data frame to store results
assortativity_h2_2 <- data.frame(network = character(), layer = character(), assortativity = numeric(), stringsAsFactors = FALSE)

c_asortativity2 = compute_heritability_assortativity2(edges_df=ceis_like[, c("Dis1", "Dis2")],
                                                    vertices_df=vertices_df2, 
                                                    attribute_name = "h2_west") # -0.03485832; considering dups: -0.03464347

for(layer in c("snps", "genes", "ppis", "pathways", "geneticcorr", "genetic")){
  # layer = "snps"
  print(layer)
  edges_df = ceis_like[, c("Dis1", "Dis2", layer)]; edges_df = edges_df[edges_df[[layer]] == 1, ]
  c_asortativity = compute_heritability_assortativity2(edges_df, vertices_df2, attribute_name = "h2_west") # -0.01485
  # Save the result in the results_df
  assortativity_h2_2 <- rbind(assortativity_h2_2, data.frame(network = "ceis_like", layer = layer, assortativity = c_asortativity))
}
head(assortativity_h2_2)
# network       layer assortativity
# 1 ceis_like        snps  0.1495815807
# 2 ceis_like       genes -0.0221568148
# 3 ceis_like        ppis -0.0286061709
# 4 ceis_like    pathways -0.0451617948
# 5 ceis_like geneticcorr  0.0008307966
# 6 ceis_like     genetic -0.0348583214

# Considering dups:
# network       layer assortativity
# 1 ceis_like        snps   0.154351100
# 2 ceis_like       genes  -0.031153907
# 3 ceis_like        ppis  -0.032612793
# 4 ceis_like    pathways  -0.046118048
# 5 ceis_like geneticcorr   0.004083564
# 6 ceis_like     genetic  -0.034643465

### NEW
assort_test_vs_zero <- function(edges_df, vertices_df, attr = "h2_west",
                                nperm = 10000, seed = 42) {
  set.seed(seed)
  # 1) normalize edge cols
  if (all(c("Dis1","Dis2") %in% names(edges_df))) {
    edges_df <- setNames(edges_df[, c("Dis1","Dis2")], c("from","to"))
  }
  stopifnot(all(c("from","to") %in% names(edges_df)))
  
  # 2) vertices: expect `icd10` and the attribute
  verts <- vertices_df[, c("icd10", attr)]
  colnames(verts) <- c("name","val")
  
  # 3) build graph (and simplify just in case)
  g <- graph_from_data_frame(edges_df, directed = FALSE, vertices = verts)
  # try(g <- simplify(g, remove.multiple = TRUE, remove.loops = TRUE))
  
  # 4) drop nodes without the attribute
  g <- delete_vertices(g, which(is.na(V(g)$val)))
  
  # 5) drop isolated vertices (degree == 0)
  if (vcount(g) > 0) {
    g <- delete_vertices(g, which(degree(g) == 0))
  }
  plot(g)
  
  # bail out if nothing left
  if (ecount(g) == 0 || vcount(g) < 2) {
    return(list(r_obs = NA, p = NA, ci_low = NA, ci_high = NA,
                n_vertices = vcount(g), n_edges = ecount(g)))
  }
  
  # 6) observed assortativity
  r_obs <- assortativity(g, values = V(g)$val, directed = FALSE)
  
  # 7) permutation null: shuffle h2 among remaining vertices
  vals <- V(g)$val; n <- length(vals)
  r_null <- replicate(nperm, {
    assortativity(g, values = sample(vals, n, replace = FALSE), directed = FALSE)
  })
  
  # 8) two-sided p vs 0 and 95% CI from permutation distribution
  p  <- (1 + sum(abs(r_null) >= abs(r_obs))) / (nperm + 1)
  ci <- quantile(r_null, c(.025, .975), na.rm = TRUE)
  p_pos <- (1 + sum(r_null >= r_obs)) / (length(r_null) + 1)
  
  p
  ci
  p_pos
  
  list(r_obs = r_obs, p = p, ci_low = unname(ci[1]), ci_high = unname(ci[2]), p_pos = p_pos,
       n_vertices = vcount(g), n_edges = ecount(g))
}

edges_snps <- subset(ceis, snps == 1, select = c(Dis1, Dis2))
res_snps <- assort_test_vs_zero(edges_snps, vertices_df, attr = "h2_west",
                                nperm = 10000, seed = 42)
res_snps
# $r_obs  -> your 0.36
# $p      -> permutation p-value vs 0
# $ci_low, $ci_high -> null CI
# $p_pos
# [1] 0.00369963

edges_snps <- subset(ceis_like, snps == 1, select = c(Dis1, Dis2))
res_snps <- assort_test_vs_zero(edges_snps, vertices_df2, attr = "h2_west",
                                nperm = 10000, seed = 42)
res_snps
# $r_obs
# [1] 0.1543511
# 
# $p
# [1] 0.05429457
# 
# $ci_low
# [1] -0.1672513
# 
# $ci_high
# [1] 0.125544
# 
# $p_pos
# [1] 0.01259874 --> under the directional hypothesis is significant, yet small modest effect. 
# 
# $n_vertices
# [1] 57
# 
# $n_edges
# [1] 131



###
# Create the plot

cassortativity_h2 <- assortativity_h2[-1, ] # Remove the first row
cassortativity_h2$layer_name = c("UKB EIs", "SNPs", "Genes", "PPIs", "Pathways", "Genetic Corr.", "Gene Expr", "Genetic", "Molecular")
gg = ggplot(cassortativity_h2, aes(x = reorder(layer_name, -assortativity), y = assortativity, fill=assortativity)) +
  geom_bar(stat = "identity", alpha=0.79) +
  theme_classic() +
  geom_text(aes(label = round(assortativity, 3), vjust = ifelse(assortativity >= 0, -0.5, 1.5)), size=3.5) +
  scale_fill_gradient(low = "#E789F0", high = "#4F0A4F") +
  # geom_text(aes(label = round(assortativity, 4)), vjust = -0.5) +
  labs(title = "",
       x = "Layer",
       y = "h2 Assortativity",
       fill="Assortativity") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
gg
ggsave(paste0(outdir,"h2_assortativity",filename_sufix,".pdf"), plot = gg, width = 7, height = 7)


dev.off()






# Calculate correlation and p-value
correlation <- cor(mdf$Perc_genetic, mdf$h2_Agg, use = "complete.obs");correlation
p_value <- cor.test(mdf$Perc_genetic, mdf$h2_Agg)$p.value; p_value

cor.test(mdf$Perc_genetic, mdf$h2_Agg)
cor.test(mdf$Perc_genetic, mdf$h2_Agg, method="spearman")

if(remove_dis){
  diseases_to_remove = c("J33 Nasal polyp",
                         # "A04 Other bacterial intestinal infections",
                         # "K63 Other diseases of intestine",
                         # "F43 Reaction to severe stress, and adjustment disorders",
                         "I71 Aortic aneurysm and dissection")
  fmdf = mdf[!(mdf$disease_name %in% diseases_to_remove), ]
}else{
  fmdf = mdf
}

# FILENAMES - original place
if(min_n_int != 0){
  
  print("Removed diseases not reaching the minimum number of interactions: ")
  removed_diseases = fmdf[fmdf$N_interactions <  min_n_int, ]$disease_name
  print(removed_diseases)
  fmdf = fmdf[fmdf$N_interactions >=  min_n_int, ]
}
# filename_sufix = paste0("_Split_dup_",substr(split_dup, 1, 1),"__","Min_n_int_",min_n_int,"__Remove_dis_",substr(remove_dis, 1, 1))
# if(with_icds_synonyms){
#   outdir="Results/with_icds_synonyms/"
#   filename_sufix = paste0(filename_sufix,paste0("__with_ICD10_synonyms"))
# }else{
#   outdir="Results/wo_icds_synonyms/"
# }


pdf(paste0(outdir,"heritability_histograms_new_",filename_sufix,".pdf"), width = 11, height = 8)

gg1 = ggplot(fmdf[!is.na(fmdf$h2_west), ], aes(x = reorder(disease_name, -h2_west), y = h2_west, fill = h2_west)) +
  geom_bar(stat = "identity") +
  theme_classic() +
  labs(x = "ICD10", y = "h2") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
gg1


# Heritability and spousal heritability
cfmdf = fmdf[!is.na(fmdf$h2_west), ]
gg1 = ggplot(cfmdf, aes(x = reorder(icd10, -h2_west), y = h2_west)) +
  geom_bar(stat = "identity") +
  theme_classic() +
  labs(x = "ICD10", y = "h2") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
gg1

gg2 = ggplot(cfmdf, aes(x = reorder(icd10, -h2_west), y = h2_west_Sp)) +
  geom_bar(stat = "identity") +
  theme_classic() +
  labs(x = "ICD10", y = "h2 Spousal") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
gg2

grid.arrange(gg1, gg2, ncol=1)

gg1 = ggplot(cfmdf, aes(x = reorder(icd10, -h2_west), y = h2_west)) +
  geom_bar(stat = "identity") +
  geom_bar(aes(x = reorder(icd10, -h2_west), y = h2_west_Sp), stat = "identity", fill="red") +
  theme_classic() +
  labs(x = "ICD10", y = "h2") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
gg1

gg1 <- ggplot(cfmdf, aes(x = reorder(icd10, -h2_west))) +
  geom_bar(aes(y = h2_west, fill = "h2_west"), stat = "identity", alpha = 0.4) +
  geom_bar(aes(y = h2_west_Sp, fill = "h2_west_Sp"), stat = "identity", alpha = 0.4) +
  scale_fill_manual(
    values = c("h2_west" = "grey50", "h2_west_Sp" = "black"),
    labels = c("h2_west" = "h2", "h2_west_Sp" = "h2 Sp")
    ) +
  theme_classic() +
  labs(x = "ICD10", y = "Heritability", fill = "") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Display the plot
print(gg1)

gg1 <- ggplot(cfmdf, aes(x = reorder(substr(disease_name,1,30), -h2_west)))+
  geom_bar(aes(y = h2_west, fill = "h2_west"), stat = "identity", alpha = 0.4) +
  geom_bar(aes(y = h2_west_Sp, fill = "h2_west_Sp"), stat = "identity", alpha = 0.4) +
  scale_fill_manual(
    values = c("h2_west" = "grey50", "h2_west_Sp" = "black"),
    labels = c("h2_west" = "h2", "h2_west_Sp" = "h2 Sp")
  ) +
  theme_classic() +
  labs(x = "ICD10", y = "Heritability", fill = "") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.margin = margin(t = 50, r = 50, b = 50, l = 50))

print(gg1)



dev.off()




ggplot(fmdf[!is.na(fmdf$h2_west_Sp), ], aes(x = reorder(icd10, -h2_west_Sp), y = h2_west_Sp)) +
  geom_bar(stat = "identity") +
  theme_classic() +
  labs(x = "ICD10", y = "h2 Spousal") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggplot(fmdf[!is.na(fmdf$h2_west_Sp), ], aes(x = reorder(disease_name, -h2_west_Sp), y = h2_west_Sp)) +
  geom_bar(stat = "identity") +
  theme_classic() +
  labs(x = "ICD10", y = "h2 Spousal") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Genetic
library(scales)
g = ggplot(fmdf[is.finite(fmdf$expr_vs_genetic), ], aes(Perc_genetic, h2_west, color=expr_vs_genetic, size=N_interactions)) +
  geom_point()+
  geom_text_repel(
    aes(label = substr(disease_name, 1, 28), color=expr_vs_genetic),
    hjust = 0, vjust = 1,
    size = 3,
    show.legend = FALSE
  ) +
  geom_point(data = fmdf[!is.finite(fmdf$expr_vs_genetic), ], aes(Perc_genetic, h2_west, size=N_interactions), colour = "#000CB5") +   # muted("#0111FF") -- #16B5AA green color
  geom_text_repel(data=fmdf[!is.finite(fmdf$expr_vs_genetic), ], aes(Perc_genetic, h2_west, label=substr(disease_name, 1, 27)), size=3, colour="#000CB5", box.padding=0.5) +
  scale_color_gradient2(low="#F2A90A", high="#0111FF", mid="#969696", midpoint=1, na.value = "#969696", space = "Lab", name = "Expr vs. Genetic", transform = transform_log()) +
  geom_abline(intercept = 0, slope = 0.01, color = "black", linetype = "dotted")+
  theme_minimal()+
  labs(title = "", x = "Genetic", y = "h2", size="#EIs", color="Ratio")
g
ggsave(paste0(outdir,"h2_vs_genetic",filename_sufix,".pdf"), plot = g, width = 10, height = 10)


library(ggExtra)
library(cowplot)
g = ggplot(fmdf, aes(h2_west, h2_west_Sp, color=expr_vs_genetic)) +
    geom_point()+
    geom_point(data = fmdf[!is.finite(fmdf$expr_vs_genetic), ], aes(h2_west, h2_west_Sp), colour = "#000CB5") +   # muted("#0111FF") -- #16B5AA green color
    geom_text_repel(data=fmdf[!is.finite(fmdf$expr_vs_genetic), ], aes(h2_west, h2_west_Sp, label=disease_name), size=2.5, colour="#000CB5", box.padding=0.5) +
    geom_text_repel(
      aes(label = disease_name, color=expr_vs_genetic),
      hjust = 0, vjust = 1,
      size = 3,
      show.legend = FALSE
    ) +
    scale_color_gradient2(low="#F2A90A", high="#0111FF", mid="#969696", midpoint=1, na.value = "#969696", space = "Lab", name = "Expr vs. Genetic", transform = transform_log()) +
    geom_abline(intercept = 0, slope = 1, color = "black", linetype = "dotted")+
    theme_minimal()+
    coord_fixed()+
    labs(title = "", x = "h2", y = "h2 Sp")
g
    
# Add histograms to the plot
g_with_histograms <- ggMarginal(g, type = "density", margins = "both", size = 15, fill="lightgrey", color="lightgrey")
  
# Display the plot
print(g_with_histograms)
# Save the plot
ggsave(paste0(outdir,"scatter_h2_Sp_coor_fix_plots_new_",filename_sufix,".pdf"), plot = g_with_histograms, width = 10, height = 9)

# Plotting only some disease names (the ones mentioned in the paper)
c1 = c("K51","F17", "J45", "D86", "M06", "K90", "E10","E11", "E66", "C50;D05", "C18", "F43", "K70", "K13", "K05", "C64", "A04", 
       "B95", "A41", "J02", "G20", "F31", "F20", "I67", "G93", "F33", "K70")
c2 = c("B34", "D68")
ttfmdf = fmdf[fmdf$icd10 %in% c(c1, c2), ]; ncol(ttfmdf) # 32
g = ggplot(fmdf, aes(h2_west, h2_west_Sp, color=expr_vs_genetic)) +
  geom_point()+
  geom_point(data = fmdf[!is.finite(fmdf$expr_vs_genetic), ], aes(h2_west, h2_west_Sp), colour = "#000CB5") +   # muted("#0111FF") -- #16B5AA green color
  geom_text_repel(data=ttfmdf[!is.finite(ttfmdf$expr_vs_genetic), ], aes(h2_west, h2_west_Sp, label=substr(disease_name, 5, 60)), size=3, colour="#000CB5", box.padding=0.5) +
  geom_text_repel(
    data=ttfmdf[is.finite(ttfmdf$expr_vs_genetic), ], aes(h2_west, h2_west_Sp, label = substr(disease_name, 5, 60), color=expr_vs_genetic),
    hjust = 0, vjust = 1,
    size = 3,
    show.legend = FALSE
  ) +
  scale_color_gradient2(low="#F2A90A", high="#0111FF", mid="#969696", midpoint=1, na.value = "#969696", space = "Lab", name = "Expr vs. Genetic", 
                        transform = transform_log(), labels = label_number(accuracy = 0.1)) +
  geom_abline(intercept = 0, slope = 1, color = "black", linetype = "dotted")+
  theme_minimal()+
  labs(title = "", x = "h2", y = "h2 Sp")
g

# Add histograms to the plot
g_with_histograms <- ggMarginal(g, type = "density", margins = "both", size = 15, fill="lightgrey", color="lightgrey")

# Display the plot
print(g_with_histograms)

# Save the plot
ggsave(paste0(outdir,"scatter_h2_Sp_coor_fix_plots_new_filtered_dis_names",filename_sufix,".pdf"), plot = g_with_histograms, width = 10, height = 9)



# dev.off()
pdf(file=paste0(outdir,"heritability_plots_new_",filename_sufix,".pdf"))
# COMPUTE CORRELATIONS
# Compute correlations between numeric values of this table
mdf_correlations <- compute_correlations(fmdf[, c("Perc_snps","Perc_genes","Perc_ppis","Perc_pathways",
                                                 "Perc_geneticcorr","Perc_genetic","Perc_geneexpr","expr_vs_genetic",
                                                 "h2_Agg","h2_Agg_catch","match",
                                                 "h2_west","h2_catch_west","h2_west_Sp")])



# Print the results
head(mdf_correlations)

# Save the results to a CSV file
write.csv(mdf_correlations, file=paste0(outdir,"heritability_correlations_",filename_sufix,".csv"), row.names = FALSE)

# Convert Variable1 and Variable2 to factors with the specified order
mdf_correlations$Variable1 <- factor(mdf_correlations$Variable1, levels = unique(mdf_correlations$Variable1))
mdf_correlations$Variable2 <- factor(mdf_correlations$Variable2, levels = unique(union(unique(mdf_correlations$Variable1),unique(mdf_correlations$Variable2))))

# Use dplyr to apply format to all numeric columns
library(dplyr)
mdf_correlations_not_scif_format <- mdf_correlations %>%
  mutate_all(~ if(is.numeric(.)) format(., scientific = FALSE))
mdf_correlations_not_scif_format = cbind(mdf_correlations[,c("Variable1","Variable2")], mdf_correlations_not_scif_format)

tmpdf = mdf_correlations

# Swap the content of Column1 and Column2
temp <- tmpdf$Variable1
tmpdf$Variable1 <- tmpdf$Variable2
tmpdf$Variable2 <- temp

full_mdf_correlations = rbind(mdf_correlations, tmpdf)
nrow(mdf_correlations)
nrow(full_mdf_correlations)


ggplot(full_mdf_correlations, aes(Variable1, Variable2, fill = Spearman_Correlation)) +
  geom_tile() +
  geom_point(data = subset(full_mdf_correlations, Spearman_p_value <= 0.05), aes(Variable1, Variable2), color = "black", size = 2) +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = 0, limit = c(-1, 1), space = "Lab", name = "Correlation") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
  labs(title = "title", x = "", y = "") +
  coord_fixed()

ggplot(full_mdf_correlations, aes(Variable1, Variable2, fill = Pearson_Correlation)) +
  geom_tile() +
  geom_point(data = subset(full_mdf_correlations, Pearson_p_value <= 0.05), aes(Variable1, Variable2), color = "black", size = 2) +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = 0, limit = c(-1, 1), space = "Lab", name = "Correlation") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
  labs(title = "title", x = "", y = "") +
  coord_fixed()

# Scatter plot - with match
cdf = full_mdf_correlations[full_mdf_correlations$Variable1 %in% c("Perc_snps","Perc_genes", "Perc_ppis", "Perc_pathways",
                                                                   "Perc_geneticcorr", "Perc_genetic", "Perc_geneexpr"), ]
cdf = cdf[cdf$Variable2 %in% c("h2_Agg", "h2_Agg_catch","match", "h2_west", "h2_catch_west"), ]
ggplot(cdf, aes(x = Variable1, y = Spearman_Correlation, color = Variable2, group=Variable2, shape=ifelse(Spearman_p_value <= 0.05, "Significant", "Not Significant"))) +
  geom_point() +
  # geom_line(alpha=0.5) +
  scale_shape_manual(values = c("Significant" = 8, "Not Significant" = 1)) +  # Set shapes manually
  labs(title = "Correlation Scatter Plot",
       x = "Features",
       y = "Correlation with Heritability") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Scatter plot - wo/ match
cdf = full_mdf_correlations[full_mdf_correlations$Variable1 %in% c("Perc_snps","Perc_genes", "Perc_ppis", "Perc_pathways",
                                                                   "Perc_geneticcorr", "Perc_genetic", "Perc_geneexpr",
                                                                   "expr_vs_genetic"), ]
cdf = cdf[cdf$Variable2 %in% c("h2_Agg", "h2_Agg_catch", "h2_west", "h2_catch_west"), ]
g1 = ggplot(cdf, aes(x = Variable1, y = Spearman_Correlation, color = Variable2, group=Variable2, shape=ifelse(Spearman_p_value <= 0.05, "Significant", "Not Significant"))) +
  geom_point(size=2) +
  geom_line(alpha=0.3) +
  scale_shape_manual(values = c("Significant" = 8, "Not Significant" = 1)) +  # Set shapes manually
  labs(title = "Correlation Scatter Plot",
       x = "Features",
       y = "Correlation with Heritability") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

g2 = ggplot(cdf[cdf$Spearman_p_value <= 0.15, ], aes(x = Variable1, y = Spearman_p_value, color = Variable2, group=Variable2, shape=ifelse(Spearman_p_value <= 0.05, "Significant", "Not Significant"))) +
  geom_point(size=2) +
  geom_line(alpha=0.3) +
  scale_shape_manual(values = c("Significant" = 8, "Not Significant" = 1)) +  # Set shapes manually
  labs(title = "Correlation Scatter Plot",
       x = "Features",
       y = "p-value") +
  geom_hline(yintercept = 0.05, linetype = "dashed", color = "grey") +
  geom_hline(yintercept = 0.01, linetype = "dashed", color = "grey") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

g3 = ggplot(cdf, aes(x = Variable1, y = Spearman_p_value, color = Variable2, group=Variable2, shape=ifelse(Spearman_p_value <= 0.05, "Significant", "Not Significant"))) +
  geom_point(size=2) +
  # geom_line(alpha=0.3) +
  scale_shape_manual(values = c("Significant" = 8, "Not Significant" = 1)) +  # Set shapes manually
  labs(title = "Correlation Scatter Plot",
       x = "Features",
       y = "p-value") +
  ylim(0, 0.1) +
  geom_hline(yintercept = 0.05, linetype = "dashed", color = "grey") +
  # geom_hline(yintercept = 0.01, linetype = "dashed", color = "grey") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

library(gridExtra)
grid.arrange(g1, g3, ncol=1)

west_cdf = cdf[cdf$Variable2 == "h2_west", ]

# Custom x-axis labels
custom_labels <- c("SNPs", "Genes", "PPIs", "Pathways", "Genetic Corr.",  "Genetic", "Gene Expr", "Expr vs. Genetic")
gg = ggplot(west_cdf, aes(x = Variable1, y = Spearman_Correlation, fill = Spearman_Correlation)) +
  geom_bar(stat = "identity") +
  scale_fill_gradient2(low = "#7E14FF", high = "#FF8D00", mid = "#EBEBEB", midpoint = 0, limit = c(-1, 1), space = "Lab", name = "Correlation") +
  geom_text(aes(label = round(Spearman_Correlation, 2)), vjust = 1.5, color = "black") +
  geom_text(aes(label = ifelse(Spearman_p_value <= 0.05, "*", "")), vjust = -0.5, color = "black") +
  theme_classic() +
  labs(title = "", y = "Correlation with h2", x="") +
  scale_x_discrete(labels = custom_labels) + # Apply custom labels
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
gg
ggsave(paste0(outdir,"Correlation_with_Heritability_west_",filename_sufix,".pdf"), plot = gg, width = 7, height = 6)

# Splitted:
g1 = ggplot(west_cdf[west_cdf$Variable1 %in% c("Perc_genetic", "Perc_geneexpr", "expr_vs_genetic"), ], aes(x = Variable1, y = Spearman_Correlation, fill = Spearman_Correlation)) +
  geom_bar(stat = "identity") +
  scale_fill_gradient2(low = "#7E14FF", high = "#FF8D00", mid = "#EBEBEB", midpoint = 0, limit = c(-1, 1), space = "Lab", name = "Correlation") +
  geom_text(aes(label = round(Spearman_Correlation, 2)), vjust = 1.5, color = "black") +
  geom_text(aes(label = ifelse(Spearman_p_value <= 0.05, "*", "")), vjust = -0.5, color = "black") +
  theme_classic() +
  labs(title = "", y = "Correlation with h2", x="") +
  # scale_x_discrete(labels = custom_labels) + # Apply custom labels
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
g1
g2 = ggplot(west_cdf[!west_cdf$Variable1 %in% c("Perc_genetic", "Perc_geneexpr", "expr_vs_genetic"), ], aes(x = Variable1, y = Spearman_Correlation, fill = Spearman_Correlation)) +
  geom_bar(stat = "identity") +
  scale_fill_gradient2(low = "#7E14FF", high = "#FF8D00", mid = "#EBEBEB", midpoint = 0, limit = c(-1, 1), space = "Lab", name = "Correlation") +
  geom_text(aes(label = round(Spearman_Correlation, 2)), vjust = 1.5, color = "black") +
  geom_text(aes(label = ifelse(Spearman_p_value <= 0.05, "*", "")), vjust = -0.5, color = "black") +
  theme_classic() +
  theme(legend.position = "bottom", legend.box = "horizontal") +
  
  labs(title = "", y = "Correlation with h2", x="") +
  # scale_x_discrete(labels = custom_labels) + # Apply custom labels
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
g2
gg = grid.arrange(g2, g1, ncol=2)
ggsave(paste0(outdir,"Correlation_with_Heritability_west_splitted_",filename_sufix,".pdf"), plot = gg, width = 7, height = 6)



gg = ggplot(mdf, aes(x = h2_Agg, y = h2_west)) +
  geom_point() +
  geom_text_repel(
    aes(label = disease_name),
    color = "black",
    hjust = 0, vjust = 1,
    size = 3,
    show.legend = FALSE
  ) +
  ggpubr::stat_cor(method = "spearman", label.x = 0.2, label.y = 0.9)
plot(gg)


fmdf$new_h2_Agg = fmdf$h2_Agg
fmdf$match = as.numeric(fmdf$match)

# Defining new metrics
fmdf$expr_minus_genetic = fmdf$Perc_geneexpr - fmdf$Perc_genetic
fmdf$Perc_only_expr0 = (fmdf$Perc_geneexpr-fmdf$Perc_genetic)/fmdf$Perc_geneexpr
fmdf$Perc_only_expr = fmdf$Perc_only_expr0
fmdf[fmdf$Perc_only_expr<0, ]$Perc_only_expr = 0

# Save the results to a CSV file
fmdf_sorted <- fmdf[order(-fmdf$h2_west_Sp), ]
write.csv(fmdf_sorted, file=paste0(outdir, "Top_Sp_",filename_sufix,".csv"), row.names = FALSE)

# Exploring new metric
sorted_by= "Perc_only_expr"
# sorted_by= "h2_west_Sp"
sorted_by= "expr_minus_genetic"


# Create subsets for the top 10, 15, 20, 25, 30, and 35 elements
nrow(fmdf_sorted)
# cfmdf = fmdf_sorted[!is.na(fmdf_sorted$Perc_only_expr), ]
cfmdf = fmdf_sorted[!is.na(fmdf_sorted$h2_west_Sp), ]
if(sorted_by == "Perc_only_expr"){
  cfmdf <- cfmdf[order(-cfmdf$Perc_only_expr), ]
}else if(sorted_by == "h2_west_Sp"){
  cfmdf <- cfmdf[order(-cfmdf$h2_west_Sp), ]
}else if(sorted_by == "expr_minus_genetic"){
  cfmdf <- cfmdf[order(-cfmdf$expr_minus_genetic), ]
}
head(cfmdf)
nrow(cfmdf)
top_10 <- cfmdf[1:10, ]
top_15 <- cfmdf[1:15, ]
top_20 <- cfmdf[1:20, ]
top_25 <- cfmdf[1:25, ]
top_30 <- cfmdf[1:30, ]
top_35 <- cfmdf[1:35, ]
all_elements <- cfmdf
# Create subsets for the bottom 10, 15, 20, 25, and 30 elements
bottom_10 <- cfmdf[(nrow(cfmdf)-9):nrow(cfmdf), ]
bottom_15 <- cfmdf[(nrow(cfmdf)-14):nrow(cfmdf), ]
bottom_20 <- cfmdf[(nrow(cfmdf)-19):nrow(cfmdf), ]
bottom_25 <- cfmdf[(nrow(cfmdf)-24):nrow(cfmdf), ]
bottom_30 <- cfmdf[(nrow(cfmdf)-29):nrow(cfmdf), ]
bottom_35 <- cfmdf[(nrow(cfmdf)-34):nrow(cfmdf), ]

# Add a new column to indicate the subset
top_10$Subset <- "Top 10"
top_15$Subset <- "Top 15"
top_20$Subset <- "Top 20"
top_25$Subset <- "Top 25"
top_30$Subset <- "Top 30"
top_35$Subset <- "Top 35"
bottom_10$Subset <- "Bottom 10"
bottom_15$Subset <- "Bottom 15"
bottom_20$Subset <- "Bottom 20"
bottom_25$Subset <- "Bottom 25"
bottom_30$Subset <- "Bottom 30"
bottom_35$Subset <- "Bottom 35"
all_elements$Subset <- "All Elements"

# Combine all subsets into a single dataframe
combined_df <- rbind(top_10, top_15, top_20, top_25, top_30, top_35,
                     bottom_10, bottom_15, bottom_20, bottom_25, bottom_30, bottom_35,
                     all_elements)

# Perform pairwise t-tests
t_test_results <- data.frame(Subset = character(), p_value = numeric(), stringsAsFactors = FALSE)

if(sorted_by == "Perc_only_expr"){
  for(subset_name in unique(combined_df$Subset)) {
    if(subset_name != "All Elements") {
      subset_data <- combined_df[combined_df$Subset == subset_name, "h2_west_Sp"]
      all_elements_data <- all_elements$h2_west_Sp
      p_value <- t.test(subset_data, all_elements_data)$p.value
      t_test_results <- rbind(t_test_results, data.frame(Subset = subset_name, p_value = p_value))
    }
  
  }
}else if(sorted_by == "h2_west_Sp"){
  for(subset_name in unique(combined_df$Subset)) {
    if(subset_name != "All Elements") {
      subset_data <- combined_df[combined_df$Subset == subset_name, "Perc_only_expr"]
      all_elements_data <- all_elements$Perc_only_expr
      p_value <- t.test(subset_data, all_elements_data)$p.value
      t_test_results <- rbind(t_test_results, data.frame(Subset = subset_name, p_value = p_value))
    }
  }
}else if(sorted_by == "expr_minus_genetic"){
  for(subset_name in unique(combined_df$Subset)) {
    if(subset_name != "All Elements") {
      subset_data <- combined_df[combined_df$Subset == subset_name, "h2_west_Sp"]
      all_elements_data <- all_elements$h2_west_Sp
      p_value <- t.test(subset_data, all_elements_data)$p.value
      t_test_results <- rbind(t_test_results, data.frame(Subset = subset_name, p_value = p_value))
    }
  }
}

# Add significance column based on p-value
t_test_results <- t_test_results %>%
  mutate(significance = ifelse(p_value < 0.05, "*", ""))

# Bottom 10 0.0002903159

# Merge the significance information back into the combined dataframe
combined_df <- combined_df %>%
  left_join(t_test_results, by = "Subset")

if(sorted_by == "Perc_only_expr"){
  # Plot the boxplots using ggplot2
  gg=ggplot(combined_df, aes(x = Subset, y = h2_west_Sp)) +
    geom_boxplot() +
    geom_text(data = t_test_results, aes(x = Subset, y = max(combined_df$h2_west_Sp) + 0.01, label = significance), vjust = -0.5) +
    labs(title = "Boxplots of h2_west_Sp for Different Subsets",
         x = "Subset",
         y = "h2_west_Sp") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  gg
}else if(sorted_by == "h2_west_Sp"){
  # Plot the boxplots using ggplot2
  gg = ggplot(combined_df, aes(x = Subset, y = Perc_only_expr)) +
    geom_boxplot() +
    geom_text(data = t_test_results, aes(x = Subset, y = max(combined_df$Perc_only_expr) + 0.01, label = significance), vjust = -0.5) +
    labs(title = "Boxplots of Perc_only_expr for Different Subsets",
         x = "Subset",
         y = "Perc_only_expr") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  gg
}else if(sorted_by == "expr_minus_genetic"){
  # Plot the boxplots using ggplot2
  gg=ggplot(combined_df, aes(x = Subset, y = h2_west_Sp)) +
    geom_boxplot() +
    geom_text(data = t_test_results, aes(x = Subset, y = max(combined_df$h2_west_Sp) + 0.01, label = significance), vjust = -0.5) +
    labs(title = "Boxplots of h2_west_Sp for Different Subsets",
         x = "Subset",
         y = "h2_west_Sp") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  gg
}

#### PARWISE T-TESTS
# Perform pairwise t-tests comparing top vs bottom subsets
t_test_results <- data.frame(Subset = character(), p_value = numeric(), stringsAsFactors = FALSE)

subset_pairs <- list(c("Top 10", "Bottom 10"), c("Top 15", "Bottom 15"), c("Top 20", "Bottom 20"), 
                     c("Top 25", "Bottom 25"), c("Top 30", "Bottom 30"), c("Top 35", "Bottom 35"))

# subset_pairs <- list(c("Top 15", "Bottom 15"), c("Top 20", "Bottom 20"), 
                     # c("Top 25", "Bottom 25"), c("Top 30", "Bottom 30"), c("Top 35", "Bottom 35"))

for(pair in subset_pairs) {
  # pair = subset_pairs[[6]]
  top_subset_name <- pair[1]
  bottom_subset_name <- pair[2]
  
  top_data <- combined_df[combined_df$Subset == top_subset_name, "h2_west_Sp"]
  bottom_data <- combined_df[combined_df$Subset == bottom_subset_name, "h2_west_Sp"]
  
  p_value <- t.test(top_data, bottom_data)$p.value
  t_test_results <- rbind(t_test_results, data.frame(Subset = top_subset_name, p_value = p_value))
}

# Add significance column based on p-value
t_test_results <- t_test_results %>%
  mutate(significance = ifelse(p_value < 0.05, "*", ""))

# Merge the significance information back into the combined dataframe
combined_df <- combined_df %>%
  left_join(t_test_results, by = "Subset")

# Plot the boxplots using ggplot2
ggplot(combined_df, aes(x = Subset, y = h2_west_Sp)) +
  geom_boxplot() +
  geom_text(data = t_test_results, aes(x = Subset, y = max(combined_df$h2_west_Sp) + 0.01, label = significance), vjust = -0.5) +
  labs(title = "Boxplots of h2_west_Sp for Different Subsets",
       x = "Subset",
       y = "h2_west_Sp") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


#######
# pdf(file=paste0("Results/heritability_plots_new_split_dup",split_dup,"_remove_dis_",remove_dis,".pdf"))

ggplot(combined_df, aes(x = Subset, y = h2_west_Sp)) +
  geom_boxplot() +
  labs(title = "Boxplots of h2_west_Sp for Different Subsets",
       x = "Subset",
       y = "h2_west_Sp") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))



# Plot the boxplots using ggplot2
ggplot(combined_df, aes(x = Subset, y = h2_west_Sp)) +
  geom_boxplot() +
  geom_text(data = t_test_results, aes(x = Subset, y = max(combined_df$h2_west_Sp) + 0.01, label = significance), vjust = -0.5) +
  labs(title = "Boxplots of h2_west_Sp for Different Subsets",
       x = "Subset",
       y = "h2_west_Sp") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


# CORRELATIONS
mcor = cbind(fmdf$Perc_genetic, fmdf$h2_Agg)
mcor <- mcor[complete.cases(mcor), ]

cor.test(mcor[,1], mcor[,2])
cor.test(mcor[,1], mcor[,2], method="spearman") # Does not assume a linear relationship between the variables
cor.test(mcor[,1], mcor[,2], method="kendall") # 
cosine(mcor[,1], mcor[,2]) # 0.9

fmdf$new_h2_AE[is.na(fmdf$new_h2_Agg)] <- -0.05

correlation=0.5

# h2 is the proportion of phenotypic variation that can be attributed to genetic differences
# Match
gg = ggplot(fmdf, aes(x = Perc_genetic, y = match, color=Perc_geneexpr)) +
  geom_point() +
  geom_text_repel(
    aes(label = disease_name),
    color = "black",
    hjust = 0, vjust = 1,
    size = 3,
    show.legend = FALSE
  ) +
  ggpubr::stat_cor(method = "spearman", label.x = 0.2, label.y = 0.9)
plot(gg)


gg = ggplot(fmdf, aes(x = Perc_genetic, y = match, color=Perc_geneexpr)) +
  geom_point() +
  geom_text_repel(
    aes(label = disease_name),
    color = "black",
    hjust = 0, vjust = 1,
    size = 3,
    show.legend = FALSE
  ) +
  ggpubr::stat_cor(method = "spearman", label.x = 0.2, label.y = 0.9)
plot(gg)


##### HERITABILITY vs. PERC_SNPS, PPIs.... #####
# WEST PLOTS
gg1 = ggplot(fmdf, aes(x = Perc_snps, y = h2_west)) +
  geom_point() +
  # geom_text_repel(
  #   aes(label = disease_name),
  #   color = "black",
  #   hjust = 0, vjust = 1,
  #   size = 3,
  #   show.legend = FALSE
  # ) +
  theme_classic()+
  ggpubr::stat_cor(method = "spearman", label.x = 0.2, label.y = 0.9)
  # ggpubr::stat_cor(method = "pearson", label.x = 0.2, label.y = 0.88)
# plot(gg1)

gg2 = ggplot(fmdf, aes(x = Perc_genes, y = h2_west)) +
  geom_point() +
  # geom_text_repel(
  #   aes(label = disease_name),
  #   color = "black",
  #   hjust = 0, vjust = 1,
  #   size = 3,
  #   show.legend = FALSE
  # ) +
  theme_classic()+
  ggpubr::stat_cor(method = "spearman", label.x = 0.2, label.y = 0.9)
  # ggpubr::stat_cor(method = "pearson", label.x = 0.2, label.y = 0.88)
# plot(gg2)

gg3 = ggplot(fmdf, aes(x = Perc_ppis, y = h2_west)) +
  geom_point() +
  # geom_text_repel(
  #   aes(label = disease_name),
  #   color = "black",
  #   hjust = 0, vjust = 1,
  #   size = 3,
  #   show.legend = FALSE
  # ) +
  theme_classic()+
  ggpubr::stat_cor(method = "spearman", label.x = 0.2, label.y = 0.9)
  # ggpubr::stat_cor(method = "pearson", label.x = 0.2, label.y = 0.88)
# plot(gg3)

gg4 = ggplot(fmdf, aes(x = Perc_pathways, y = h2_west)) +
  geom_point() +
  # geom_text_repel(
  #   aes(label = disease_name),
  #   color = "black",
  #   hjust = 0, vjust = 1,
  #   size = 3,
  #   show.legend = FALSE
  # ) +
  theme_classic()+
  ggpubr::stat_cor(method = "spearman", label.x = 0.2, label.y = 0.9)
  # ggpubr::stat_cor(method = "pearson", label.x = 0.2, label.y = 0.88)
# plot(gg4)

gg5 = ggplot(fmdf, aes(x = Perc_geneticcorr, y = h2_west)) +
  geom_point() +
  # geom_text_repel(
  #   aes(label = disease_name),
  #   color = "black",
  #   hjust = 0, vjust = 1,
  #   size = 3,
  #   show.legend = FALSE
  # ) +
  theme_classic()+
  ggpubr::stat_cor(method = "spearman", label.x = 0.2, label.y = 0.9)
  # ggpubr::stat_cor(method = "pearson", label.x = 0.2, label.y = 0.88)
# plot(gg5)

gg6 = ggplot(fmdf, aes(x = Perc_genetic, y = h2_west)) +
  geom_point() +
  # geom_text_repel(
  #   aes(label = disease_name),
  #   color = "black",
  #   hjust = 0, vjust = 1,
  #   size = 3,
  #   show.legend = FALSE
  # ) +
  theme_classic()+
  ggpubr::stat_cor(method = "spearman", label.x = 0.2, label.y = 0.9)
  # ggpubr::stat_cor(method = "pearson", label.x = 0.2, label.y = 0.88)
# plot(gg6)

gg7 = ggplot(fmdf, aes(x = Perc_geneexpr, y = h2_west)) +
  geom_point() +
  # geom_text_repel(
  #   aes(label = disease_name),
  #   color = "black",
  #   hjust = 0, vjust = 1,
  #   size = 3,
  #   show.legend = FALSE
  # ) +
  theme_classic()+
  ggpubr::stat_cor(method = "spearman", label.x = 0.2, label.y = 0.9)
  # ggpubr::stat_cor(method = "pearson", label.x = 0.2, label.y = 0.88)
# plot(gg7)


gg8 = ggplot(fmdf, aes(x = expr_vs_genetic, y = h2_west)) +
  geom_point() +
  # geom_text_repel(
  #   aes(label = disease_name),
  #   color = "black",
  #   hjust = 0, vjust = 1,
  #   size = 3,
  #   show.legend = FALSE
  # ) +
  theme_classic()+
  ggpubr::stat_cor(method = "spearman", label.x = 0.2, label.y = 0.9)
  # ggpubr::stat_cor(method = "pearson", label.x = 0.2, label.y = 0.88)
# plot(gg8)

grid.arrange(gg1, gg2, gg3, gg4, gg5, gg6, gg7, gg8, ncol=3)


# SNPS
gg = ggplot(fmdf, aes(x = Perc_snps, y = h2_Agg)) +
  geom_point() +
  geom_text_repel(
    aes(label = disease_name),
    color = "black",
    hjust = 0, vjust = 1,
    size = 3,
    show.legend = FALSE
  ) +
  ggpubr::stat_cor(method = "spearman", label.x = 0.2, label.y = 0.9)+
  ggpubr::stat_cor(method = "pearson", label.x = 0.2, label.y = 0.88)
plot(gg)

gg = ggplot(fmdf, aes(x = Perc_snps, y = h2_Agg_catch)) +
  geom_point() +
  geom_text_repel(
    aes(label = disease_name),
    color = "black",
    hjust = 0, vjust = 1,
    size = 3,
    show.legend = FALSE
  ) +
  ggpubr::stat_cor(method = "spearman", label.x = 0.2, label.y = 0.9)+
  ggpubr::stat_cor(method = "pearson", label.x = 0.2, label.y = 0.88)
plot(gg)

gg = ggplot(fmdf, aes(x = Perc_snps, y = h2_west)) +
  geom_point() +
  geom_text_repel(
    aes(label = disease_name),
    color = "black",
    hjust = 0, vjust = 1,
    size = 3,
    show.legend = FALSE
  ) +
  ggpubr::stat_cor(method = "spearman", label.x = 0.2, label.y = 0.9)+
  ggpubr::stat_cor(method = "pearson", label.x = 0.2, label.y = 0.88)
plot(gg)

gg = ggplot(fmdf, aes(x = Perc_snps, y = h2_catch_west)) +
  geom_point() +
  geom_text_repel(
    aes(label = disease_name),
    color = "black",
    hjust = 0, vjust = 1,
    size = 3,
    show.legend = FALSE
  ) +
  ggpubr::stat_cor(method = "spearman", label.x = 0.2, label.y = 0.9)+
  ggpubr::stat_cor(method = "pearson", label.x = 0.2, label.y = 0.88)
plot(gg)

# genes
gg = ggplot(fmdf, aes(x = Perc_genes, y = h2_Agg)) +
  geom_point() +
  geom_text_repel(
    aes(label = disease_name),
    color = "black",
    hjust = 0, vjust = 1,
    size = 3,
    show.legend = FALSE
  ) +
  ggpubr::stat_cor(method = "spearman", label.x = 0.2, label.y = 0.9)+
  ggpubr::stat_cor(method = "pearson", label.x = 0.2, label.y = 0.88)
plot(gg)

gg = ggplot(fmdf, aes(x = Perc_genes, y = h2_Agg_catch)) +
  geom_point() +
  geom_text_repel(
    aes(label = disease_name),
    color = "black",
    hjust = 0, vjust = 1,
    size = 3,
    show.legend = FALSE
  ) +
  ggpubr::stat_cor(method = "spearman", label.x = 0.2, label.y = 0.9)+
  ggpubr::stat_cor(method = "pearson", label.x = 0.2, label.y = 0.88)
plot(gg)

gg = ggplot(fmdf, aes(x = Perc_genes, y = h2_west)) +
  geom_point() +
  geom_text_repel(
    aes(label = disease_name),
    color = "black",
    hjust = 0, vjust = 1,
    size = 3,
    show.legend = FALSE
  ) +
  ggpubr::stat_cor(method = "spearman", label.x = 0.2, label.y = 0.9)+
  ggpubr::stat_cor(method = "pearson", label.x = 0.2, label.y = 0.88)
plot(gg)

gg = ggplot(fmdf, aes(x = Perc_genes, y = h2_catch_west)) +
  geom_point() +
  geom_text_repel(
    aes(label = disease_name),
    color = "black",
    hjust = 0, vjust = 1,
    size = 3,
    show.legend = FALSE
  ) +
  ggpubr::stat_cor(method = "spearman", label.x = 0.2, label.y = 0.9)+
  ggpubr::stat_cor(method = "pearson", label.x = 0.2, label.y = 0.88)
plot(gg)

# ppis
gg = ggplot(fmdf, aes(x = Perc_ppis, y = h2_Agg)) +
  geom_point() +
  geom_text_repel(
    aes(label = disease_name),
    color = "black",
    hjust = 0, vjust = 1,
    size = 3,
    show.legend = FALSE
  ) +
  ggpubr::stat_cor(method = "spearman", label.x = 0.2, label.y = 0.9)+
  ggpubr::stat_cor(method = "pearson", label.x = 0.2, label.y = 0.88)
plot(gg)

gg = ggplot(fmdf, aes(x = Perc_ppis, y = h2_Agg_catch)) +
  geom_point() +
  geom_text_repel(
    aes(label = disease_name),
    color = "black",
    hjust = 0, vjust = 1,
    size = 3,
    show.legend = FALSE
  ) +
  ggpubr::stat_cor(method = "spearman", label.x = 0.2, label.y = 0.9)+
  ggpubr::stat_cor(method = "pearson", label.x = 0.2, label.y = 0.88)
plot(gg)

gg = ggplot(fmdf, aes(x = Perc_ppis, y = h2_west)) +
  geom_point() +
  geom_text_repel(
    aes(label = disease_name),
    color = "black",
    hjust = 0, vjust = 1,
    size = 3,
    show.legend = FALSE
  ) +
  ggpubr::stat_cor(method = "spearman", label.x = 0.2, label.y = 0.9)+
  ggpubr::stat_cor(method = "pearson", label.x = 0.2, label.y = 0.88)
plot(gg)

gg = ggplot(fmdf, aes(x = Perc_ppis, y = h2_catch_west)) +
  geom_point() +
  geom_text_repel(
    aes(label = disease_name),
    color = "black",
    hjust = 0, vjust = 1,
    size = 3,
    show.legend = FALSE
  ) +
  ggpubr::stat_cor(method = "spearman", label.x = 0.2, label.y = 0.9)+
  ggpubr::stat_cor(method = "pearson", label.x = 0.2, label.y = 0.88)
plot(gg)

# geneticcorr
gg = ggplot(fmdf, aes(x = Perc_geneticcorr, y = h2_Agg)) +
  geom_point() +
  geom_text_repel(
    aes(label = disease_name),
    color = "black",
    hjust = 0, vjust = 1,
    size = 3,
    show.legend = FALSE
  ) +
  ggpubr::stat_cor(method = "spearman", label.x = 0.2, label.y = 0.9)+
  ggpubr::stat_cor(method = "pearson", label.x = 0.2, label.y = 0.88)
plot(gg)

gg = ggplot(fmdf, aes(x = Perc_geneticcorr, y = h2_Agg_catch)) +
  geom_point() +
  geom_text_repel(
    aes(label = disease_name),
    color = "black",
    hjust = 0, vjust = 1,
    size = 3,
    show.legend = FALSE
  ) +
  ggpubr::stat_cor(method = "spearman", label.x = 0.2, label.y = 0.9)+
  ggpubr::stat_cor(method = "pearson", label.x = 0.2, label.y = 0.88)
plot(gg)

gg = ggplot(fmdf, aes(x = Perc_geneticcorr, y = h2_west)) +
  geom_point() +
  geom_text_repel(
    aes(label = disease_name),
    color = "black",
    hjust = 0, vjust = 1,
    size = 3,
    show.legend = FALSE
  ) +
  ggpubr::stat_cor(method = "spearman", label.x = 0.2, label.y = 0.9)+
  ggpubr::stat_cor(method = "pearson", label.x = 0.2, label.y = 0.88)
plot(gg)

gg = ggplot(fmdf, aes(x = Perc_geneticcorr, y = h2_catch_west)) +
  geom_point() +
  geom_text_repel(
    aes(label = disease_name),
    color = "black",
    hjust = 0, vjust = 1,
    size = 3,
    show.legend = FALSE
  ) +
  ggpubr::stat_cor(method = "spearman", label.x = 0.2, label.y = 0.9)+
  ggpubr::stat_cor(method = "pearson", label.x = 0.2, label.y = 0.88)
plot(gg)

# pathways
gg = ggplot(fmdf, aes(x = Perc_genetic, y = h2_Agg)) +
  geom_point() +
  geom_text_repel(
    aes(label = disease_name),
    color = "black",
    hjust = 0, vjust = 1,
    size = 3,
    show.legend = FALSE
  ) +
  ggpubr::stat_cor(method = "spearman", label.x = 0.2, label.y = 0.9)+
  ggpubr::stat_cor(method = "pearson", label.x = 0.2, label.y = 0.88)
plot(gg)

gg = ggplot(fmdf, aes(x = Perc_genetic, y = h2_Agg_catch)) +
  geom_point() +
  geom_text_repel(
    aes(label = disease_name),
    color = "black",
    hjust = 0, vjust = 1,
    size = 3,
    show.legend = FALSE
  ) +
  ggpubr::stat_cor(method = "spearman", label.x = 0.2, label.y = 0.9)+
  ggpubr::stat_cor(method = "pearson", label.x = 0.2, label.y = 0.88)
plot(gg)

gg = ggplot(fmdf, aes(x = Perc_genetic, y = h2_west)) +
  geom_point() +
  geom_text_repel(
    aes(label = disease_name),
    color = "black",
    hjust = 0, vjust = 1,
    size = 3,
    show.legend = FALSE
  ) +
  ggpubr::stat_cor(method = "spearman", label.x = 0.2, label.y = 0.9)+
  ggpubr::stat_cor(method = "pearson", label.x = 0.2, label.y = 0.88)
plot(gg)

gg = ggplot(fmdf, aes(x = Perc_genetic, y = h2_catch_west)) +
  geom_point() +
  geom_text_repel(
    aes(label = disease_name),
    color = "black",
    hjust = 0, vjust = 1,
    size = 3,
    show.legend = FALSE
  ) +
  ggpubr::stat_cor(method = "spearman", label.x = 0.2, label.y = 0.9)+
  ggpubr::stat_cor(method = "pearson", label.x = 0.2, label.y = 0.88)
plot(gg)

# genetic
gg = ggplot(fmdf, aes(x = Perc_genetic, y = h2_Agg)) +
  geom_point() +
  geom_text_repel(
    aes(label = disease_name),
    color = "black",
    hjust = 0, vjust = 1,
    size = 3,
    show.legend = FALSE
  ) +
  ggpubr::stat_cor(method = "spearman", label.x = 0.2, label.y = 0.9)+
  ggpubr::stat_cor(method = "pearson", label.x = 0.2, label.y = 0.88)
plot(gg)

gg = ggplot(fmdf, aes(x = Perc_genetic, y = h2_Agg_catch)) +
  geom_point() +
  geom_text_repel(
    aes(label = disease_name),
    color = "black",
    hjust = 0, vjust = 1,
    size = 3,
    show.legend = FALSE
  ) +
  ggpubr::stat_cor(method = "spearman", label.x = 0.2, label.y = 0.9)+
  ggpubr::stat_cor(method = "pearson", label.x = 0.2, label.y = 0.88)
plot(gg)

gg = ggplot(fmdf, aes(x = Perc_genetic, y = h2_west)) +
  geom_point() +
  geom_text_repel(
    aes(label = disease_name),
    color = "black",
    hjust = 0, vjust = 1,
    size = 3,
    show.legend = FALSE
  ) +
  ggpubr::stat_cor(method = "spearman", label.x = 0.2, label.y = 0.9)+
  ggpubr::stat_cor(method = "pearson", label.x = 0.2, label.y = 0.88)
plot(gg)

gg = ggplot(fmdf, aes(x = Perc_genetic, y = h2_catch_west)) +
  geom_point() +
  geom_text_repel(
    aes(label = disease_name),
    color = "black",
    hjust = 0, vjust = 1,
    size = 3,
    show.legend = FALSE
  ) +
  ggpubr::stat_cor(method = "spearman", label.x = 0.2, label.y = 0.9)+
  ggpubr::stat_cor(method = "pearson", label.x = 0.2, label.y = 0.88)
plot(gg)

gg = ggplot(fmdf, aes(x = Perc_genetic, y = h2_west_Sp)) +
  geom_point() +
  geom_text_repel(
    aes(label = disease_name),
    color = "black",
    hjust = 0, vjust = 1,
    size = 3,
    show.legend = FALSE
  ) +
  ggpubr::stat_cor(method = "spearman", label.x = 0.2, label.y = 0.9)+
  ggpubr::stat_cor(method = "pearson", label.x = 0.2, label.y = 0.88)
plot(gg)


# Genexpre
gg = ggplot(fmdf, aes(x = Perc_geneexpr, y = h2_Agg)) +
  geom_point() +
  geom_text_repel(
    aes(label = disease_name),
    color = "black",
    hjust = 0, vjust = 1,
    size = 3,
    show.legend = FALSE
  ) +
  ggpubr::stat_cor(method = "spearman", label.x = 0.2, label.y = 0.9)+
  ggpubr::stat_cor(method = "pearson", label.x = 0.2, label.y = 0.88)
plot(gg)

gg = ggplot(fmdf, aes(x = Perc_geneexpr, y = h2_Agg_catch)) +
  geom_point() +
  geom_text_repel(
    aes(label = disease_name),
    color = "black",
    hjust = 0, vjust = 1,
    size = 3,
    show.legend = FALSE
  ) +
  ggpubr::stat_cor(method = "spearman", label.x = 0.2, label.y = 0.9)+
  ggpubr::stat_cor(method = "pearson", label.x = 0.2, label.y = 0.88)
plot(gg)

gg = ggplot(fmdf, aes(x = Perc_geneexpr, y = h2_west)) +
  geom_point() +
  geom_text_repel(
    aes(label = disease_name),
    color = "black",
    hjust = 0, vjust = 1,
    size = 3,
    show.legend = FALSE
  ) +
  ggpubr::stat_cor(method = "spearman", label.x = 0.2, label.y = 0.9)+
  ggpubr::stat_cor(method = "pearson", label.x = 0.2, label.y = 0.88)
plot(gg)

gg = ggplot(fmdf, aes(x = Perc_geneexpr, y = h2_catch_west)) +
  geom_point() +
  geom_text_repel(
    aes(label = disease_name),
    color = "black",
    hjust = 0, vjust = 1,
    size = 3,
    show.legend = FALSE
  ) +
  ggpubr::stat_cor(method = "spearman", label.x = 0.2, label.y = 0.9)+
  ggpubr::stat_cor(method = "pearson", label.x = 0.2, label.y = 0.88)
plot(gg)

gg = ggplot(fmdf, aes(x = Perc_geneexpr, y = h2_west_Sp)) +
  geom_point() +
  geom_text_repel(
    aes(label = disease_name),
    color = "black",
    hjust = 0, vjust = 1,
    size = 3,
    show.legend = FALSE
  ) +
  ggpubr::stat_cor(method = "spearman", label.x = 0.2, label.y = 0.9)+
  ggpubr::stat_cor(method = "pearson", label.x = 0.2, label.y = 0.88)
plot(gg)

gg = ggplot(fmdf, aes(x = expr_vs_genetic, y = h2_west_Sp)) +
  geom_point() +
  geom_text_repel(
    aes(label = disease_name),
    color = "black",
    hjust = 0, vjust = 1,
    size = 3,
    show.legend = FALSE
  ) +
  ggpubr::stat_cor(method = "spearman", label.x = 0.2, label.y = 0.9)+
  ggpubr::stat_cor(method = "pearson", label.x = 0.2, label.y = 0.88)
plot(gg)


dev.off()



####### Sp PLOTS #######

# Histograms
ggplot(fmdf, aes(x = h2_west_Sp)) +
  geom_histogram() +       # Adds points to the plot
  # geom_line(color = "blue") +        # Adds lines connecting the points
  labs(title = "Plot of Continuous Values", x = "Index", y = "Value") + # Labels
  theme_minimal()                   # Use a minimal theme for better aesthetics

ggplot(fmdf, aes(x = h2_west_Sp)) +
  geom_bar() +       # Adds points to the plot
  # geom_line(color = "blue") +        # Adds lines connecting the points
  labs(title = "Plot of Continuous Values", x = "Index", y = "Value") + # Labels
  theme_minimal()  


gg = ggplot(fmdf, aes(x = expr_vs_genetic, y = h2_west_Sp)) +
  geom_point() +
  geom_text_repel(
    aes(label = disease_name),
    color = "black",
    hjust = 0, vjust = 1,
    size = 3,
    show.legend = FALSE
  ) +
  ggpubr::stat_cor(method = "spearman", label.x = 0.2, label.y = 0.9)+
  ggpubr::stat_cor(method = "pearson", label.x = 0.2, label.y = 0.88)
plot(gg)

gg = ggplot(fmdf[fmdf$Perc_only_expr >= 0, ], aes(x = Perc_only_expr, y = h2_west_Sp)) +
  geom_point() +
  geom_text_repel(
    aes(label = disease_name),
    color = "black",
    hjust = 0, vjust = 1,
    size = 3,
    show.legend = FALSE
  ) +
  ggpubr::stat_cor(method = "spearman", label.x = 0.2, label.y = 0.9)+
  ggpubr::stat_cor(method = "pearson", label.x = 0.2, label.y = 0.88)
plot(gg)

gg = ggplot(fmdf, aes(x = Perc_only_expr, y = h2_west_Sp)) +
  geom_point() +
  geom_text_repel(
    aes(label = disease_name),
    color = "black",
    hjust = 0, vjust = 1,
    size = 3,
    show.legend = FALSE
  ) +
  ylim(0, 0.3) +
  ggpubr::stat_cor(method = "spearman", label.x = 0.2, label.y = 0.9)+
  ggpubr::stat_cor(method = "pearson", label.x = 0.2, label.y = 0.88)
plot(gg)

gg = ggplot(fmdf, aes(x = Perc_only_expr, y = h2_west_Sp)) +
  geom_point() +
  ylim(0, 0.3) +
  ggpubr::stat_cor(method = "spearman", label.x = 0.2, label.y = 0.9)+
  ggpubr::stat_cor(method = "pearson", label.x = 0.2, label.y = 0.88)
plot(gg)


stop("End of the script before fitting curves")

####### FITTING CURVES ######
pdf(file=paste0(outdir,"fitting_curves_",split_dup,"_remove_dis_",remove_dis,".pdf"))

y = fmdf[fmdf$icd10 != "C61", ]$Perc_only_expr
x = fmdf[fmdf$icd10 != "C61", ]$h2_west_Sp

# Fit a linear model
lm_model <- lm(y ~ x)

# Print summary of the model
summary(lm_model)

# Plot the data and the fitted line
plot(x, y, main = "Linear Regression")
abline(lm_model, col = "red")

# Fit a nonlinear model
nls_model <- nls(y ~ a * exp(b * x), start = list(a = 1, b = 0.5))

# Print summary of the model
summary(nls_model)

# Plot the data and the fitted curve
plot(x, y, main = "Nonlinear Regression")
lines(x, predict(nls_model, list(x = x)), col = "blue")

# Fit a polynomial model of degree 2
poly_model <- lm(y ~ poly(x, 2))

# Print summary of the model
summary(poly_model)

# Plot the data and the fitted curve
plot(x, y, main = "Polynomial Regression")
lines(x, predict(poly_model, data.frame(x = x)), col = "green")

# Fit a smoothing spline
spline_model <- smooth.spline(x, y)

# Print summary of the model
summary(spline_model)

# Plot the data and the fitted spline
plot(x, y, main = "Smoothing Spline")
lines(spline_model, col = "purple")


# Plot the data and the fitted polynomial
ggplot(fmdf, aes(x = Perc_only_expr, y = h2_west_Sp)) +
  geom_point() +
  geom_smooth(method = "lm", formula = y ~ poly(x, 2), col = "green") +
  ggtitle("Polynomial Model Fit (Degree 2)")

# Residuals for linear model
plot(lm_model, which = 1)

# Residuals for polynomial model
plot(poly_model, which = 1)

# Residuals for nonlinear model
plot(nls_model, which = 1)

# Fit a nonlinear model (e.g., exponential)
nls_model <- nls(h2_west_Sp ~ a * exp(b * Perc_only_expr), data = fmdf[fmdf$Perc_only_expr >= -1, ], start = list(a = 1, b = 0.1))

# Plot the data and the fitted nonlinear model
fmdf[fmdf$Perc_only_expr >= -1, ]$nls_pred <- predict(nls_model, newdata = fmdf[fmdf$Perc_only_expr >= -1, ])
ggplot(fmdf[fmdf$Perc_only_expr >= -1, ], aes(x = Perc_only_expr, y = h2_west_Sp)) +
  geom_point() +
  geom_line(aes(y = nls_pred), col = "blue") +
  ggtitle("Nonlinear Model Fit")

# Residuals for nonlinear model
plot(nls_model, which = 1)

# Calculate residuals
a <- 0.03474
b <- 0.51310
fmdf[fmdf$Perc_only_expr >= -1, ]$predicted <- a * exp(b * fmdf[fmdf$Perc_only_expr >= -1, ]$Perc_only_expr)
fmdf[fmdf$Perc_only_expr >= -1, ]$residuals <- fmdf[fmdf$Perc_only_expr >= -1, ]$h2_west_Sp - fmdf[fmdf$Perc_only_expr >= -1, ]$predicted

# Plot residuals
ggplot(fmdf[fmdf$Perc_only_expr >= -1, ], aes(x = Perc_only_expr, y = residuals)) +
  geom_point() +
  geom_hline(yintercept = 0, linetype = "dashed", col = "red") +
  labs(title = "Residuals of Nonlinear Regression Model",
       x = "Perc_only_expr",
       y = "Residuals")

gg = ggplot(fmdf, aes(x = expr_minus_genetic, y = h2_west_Sp)) +
  geom_point() +
  geom_text_repel(
    aes(label = disease_name),
    color = "black",
    hjust = 0, vjust = 1,
    size = 3,
    show.legend = FALSE
  ) +
  ggpubr::stat_cor(method = "spearman", label.x = 0.2, label.y = 0.9)+
  ggpubr::stat_cor(method = "pearson", label.x = 0.2, label.y = 0.88)
plot(gg)

dev.off()


#### OLD PLOTS #####
gg = ggplot(fmdf, aes(x = expr_minus_genetic, y = h2_Agg, color=Perc_geneticcorr)) +
  geom_point() +
  geom_text_repel(
    aes(label = disease_name),
    color = "black",
    hjust = 0, vjust = 1,
    size = 3,
    show.legend = FALSE
  )  +
  ggpubr::stat_cor(method = "spearman", label.x = 0.2, label.y = 0.9)
plot(gg)

# SNPS
gg = ggplot(fmdf, aes(x = Perc_genetic, y = h2_Agg)) +
  geom_point() +
  geom_text_repel(
    aes(label = disease_name),
    color = "black",
    hjust = 0, vjust = 1,
    size = 3,
    show.legend = FALSE
  ) +
  labs(x="%EIs explained by genetics", y="h2") +
  ggpubr::stat_cor(method = "spearman", label.x = 0.2, label.y = 0.9)
plot(gg)

gg = ggplot(fmdf, aes(x = Perc_genetic, y = h2_Agg, color=Perc_geneexpr)) +
  geom_point() +
  geom_text_repel(
    aes(label = disease_name),
    color = "black",
    hjust = 0, vjust = 1,
    size = 3,
    show.legend = FALSE
  ) +
  labs(x="%EIs explained by genetics", y="h2", color="%EIs expression") +
  ggpubr::stat_cor(method = "spearman", label.x = 0.2, label.y = 0.9)
plot(gg)


gg = ggplot(fmdf, aes(x = Perc_genetic, y = h2_Agg, color=Perc_geneexpr, size=N_interactions)) +
  geom_point() +
  geom_text_repel(
    aes(label = disease_name),
    color = "black",
    hjust = 0, vjust = 1,
    size = 3,
    show.legend = FALSE
  )  +
  labs(x="%EIs explained by genetics", y="h2", color="%EIs expression", size="N Int") +
  ggpubr::stat_cor(method = "spearman", label.x = 0.2, label.y = 0.9)
plot(gg)

gg = ggplot(fmdf, aes(x = Perc_genetic, y = new_h2_Agg)) +
  geom_jitter(width = 0.02, height = 0.02) +
  geom_text_repel(
    aes(label = disease_name),
    color = "black",
    hjust = 0, vjust = 1,
    size = 3,
    show.legend = FALSE
  )  +
  ggpubr::stat_cor(method = "spearman", label.x = 0.2, label.y = 0.9)
plot(gg)

gg = ggplot(fmdf, aes(x = Perc_geneexpr, y = h2_Agg, color=Perc_genetic)) +
  geom_jitter(width = 0.02, height = 0.02) +
  geom_text_repel(
    aes(label = disease_name),
    color = "black",
    hjust = 0, vjust = 1,
    size = 3,
    show.legend = FALSE
  )  +
  ggpubr::stat_cor(method = "spearman", label.x = 0.2, label.y = 0.9)
plot(gg)

