#####################################################################################
# Topological analysis of the DSN and SSN networks
#
# Author: Beatriz Urda García, 2025
######################################################################################


source("GEVtools.R")
# save_genome_medicine_final_networks()

dsn_filename = "Network_building/Defined_networks/all_diseases/DEGs/pairwise_union_cosine_distance_sDEGs_network_1FDR.txt"
ssn_filename = "Network_building/Defined_networks/all_diseases/metapatients_and_disease/metap_dis_pairwise_union_cosine_distance_sDEGs_network_1FDR.txt"
out_dsn = "all_diseases/DEGs/pairwise_union_cosine_distance_sDEGs_network_1FDR_net_analysis"
out_ssn = "all_diseases/metapatients_and_disease/metap_dis_pairwise_union_cosine_distance_sDEGs_network_1FDR_net_analysis"
dsn = read.csv2(dsn_filename, sep="\t")
head(dsn)
length(unique(union(dsn$Dis1, dsn$Dis2))) # 148

ssn = read.csv2(ssn_filename, sep="\t")
head(ssn)
length(unique(union(ssn$Dis1, ssn$Dis2))) # 511

# Are there pairs with Distance == 0? No.
nrow(dsn)
nrow(dsn[dsn$Distance == 0, ]) # 0

nrow(ssn)
nrow(ssn[ssn$Distance == 0, ]) # 0

nets = list(dsn, ssn); colnames(nets) = c("dsn", "ssn")

library(plotly)
library(ggplot2)
library(ggplotify)
library(gplots)
library(heatmaply)



net_list = list(dsn, ssn)
names(net_list) = c("final_dsn", "final_ssn")
obtain_cliques = FALSE
properties_df = data.frame()

for (k in c(1:2)){
  df = net_list[[k]]; nint_total = nrow(df)
  df_name = names(net_list[k]); print(df_name)
  df$Distance = as.numeric(df$Distance); df$pvalue = as.numeric(df$pvalue); df$adj_pvalue = as.numeric(df$adj_pvalue)
  posi <- df[df$Distance > 0,]; nint_posi = nrow(posi)   
  negat <- df[df$Distance < 0,]; nint_negat = nrow(negat) 
  n_nodes = length(unique(union(df$Dis1, df$Dis2)))
  
  print("Starting network analysis.....................")
  cdf <- data.frame(Feature = "N nodes", Value = n_nodes)
  cfeat <- c('N sign int','N pos int','N neg int','Perc pos int','Perc neg int' ) 
  cvals <- c(nint_total, nint_posi, nint_negat, (nint_posi / nint_total)*100, (nint_negat / nint_total)*100 )
  cdf <- rbind(cdf, data.frame(Feature=cfeat,Value = cvals))
  print(head(cdf))
  adf <- analyze_network(dist_filename="", out_dsn, df, posi, negat)
  cdf <- rbind(cdf, adf)
  print(head(cdf))
  print("Network analyzed.....................")
  
  
  
  graph <- graph_from_data_frame(df, directed=FALSE)
  netm <- get.adjacency(graph, attr="Distance", sparse=F)
  p = heatmap(netm, col = bluered(100), 
            scale="none", margins=c(16,16), trace="none" )
  p
  heatmap.2(netm, col = bluered(100), 
            scale="none", margins=c(16,16), trace="none" )
  # heatmap_gg <- ggplotify::as.ggplot(heatmap)
  # ggplotly(heatmap_gg)
  
  heatmaply(netm, 
            xlab = "", 
            ylab = "", 
            main = "Disease Similarity Network", 
            dendrogram = "both", col= bluered(100),
            colorbar_title = "Distance",
            layout = list(xaxis = list(title = "", tickfont = list(size = 2)),
                          yaxis = list(title = "", tickfont = list(size = 2)))
  )
  
  
  cnet_list <- list(df,posi,negat)
  
  j <- 1
  for(df in cnet_list){
    if(j == 1){
      surname = ""
    }else if(j == 2){
      surname = "_pos"
    }else if(j == 3){
      surname = "_neg"
    }
    
    try(dev.off())
    pdf(paste("Network_building/Network_analysis/Topology/",df_name,surname,"_topology_plots.pdf",sep=""), width = 13)
    nint_total = nrow(df)
    pos <- df[df$Distance > 0,]; nint_pos = nrow(pos)   
    neg <- df[df$Distance < 0,]; nint_neg = nrow(neg) 
    n_nodes = length(unique(union(df$Dis1, df$Dis2)))
    
    graph <- graph_from_data_frame(df, directed=FALSE)
    E(graph)$weight <- abs(df$Distance) 
    j = j+1
    
    # Create a distance graph with Distance as weight
    cex_axis <- 4.5
    epsilon = 0.00000001 # 1*10-8
    distance_graph = graph
    output_filename = "mynets"
    if(output_filename == "barabasi"){
      max_barabasi = max(E(graph)$weight)
      E(distance_graph)$weight = (1 - (E(graph)$weight/E(graph)$weight) + epsilon)
    }else{
      E(distance_graph)$weight = (1 - E(graph)$weight + epsilon)
    }
    
    
    connected_components <- clusters(graph)
    n_connected_components <- connected_components$no
    size_connected_components <- connected_components$csize
    if(length(size_connected_components) > 1){
      size_connected_components <- paste(size_connected_components,collapse="_")
    }
    
    deg <- degree(graph, mode="all")
    mean_deg <- mean(deg) 
    
    #Density of a graph: The proportion of present edges from all possible edges in the network
    edge_density <- edge_density(graph, loops=F)
    edge_density
    
    # Reciprocity: it only makes sense for directed graphs
    reciprocity(graph)
    
    # dyad_census(graph) # Mutual, asymmetric, and nyll node pairs
    # 2*dyad_census(graph)$mut/ecount(graph) # Calculating reciprocity
    
    # Transitivity: Transitivity measures the probability that the adjacent vertices of a vertex are connected
    transitivity <- transitivity(graph, type="global") # fot the network (global)
    local_trans <- transitivity(graph, vids=V(graph), type="local") # for each of the nodes (local)
    transitivity_table <- data.frame(nodes=names(V(graph)), transitivity=local_trans)
    print(ggplot(transitivity_table, aes(reorder(nodes,transitivity), transitivity)) + geom_bar(stat="identity") +
            theme_classic()+
            theme(axis.text.x = element_text(angle = 45,hjust=1, size=cex_axis),plot.margin = margin(0.8, 0.8, 0.8, 0.8, "cm"))+
            # ggtitle("Transitivity") + 
            # xlab("Diseases") + 
            ylab("Transitivity"))
    
    ### Diameter: longest geodesic distance (length of the shortest path between two nodes) in the network
    # With weights
    diameter <- diameter(distance_graph, directed=F)
    get_diameter(distance_graph, directed=F)
    # Without weights
    diameter_wo_weights <- diameter(graph, directed=F, weights=NA)
    get_diameter(graph, directed=F, weights=NA)
    
    ###Degree
    deg <- degree(graph, mode="all")
    plot(graph, vertex.size=deg*0.1,layout=layout.circle)
    hist(deg, breaks=1:vcount(graph)-1, main="Histogram of node degree")
    
    # Degree distribution
    deg.dist <- degree_distribution(graph, cumulative=T, mode="all")
    
    plot( x=0:max(deg), y=1-deg.dist, pch=19, cex=1.2, col="orange", 
          xlab="Degree", ylab="Cumulative Frequency")
    
    ### CENTRALITY AND CENTRALIZATION
    # Degree (number of ties) ### hacer un plot de esto con ggplot quizá!
    degree_table <- data.frame('nodes'=names(degree(graph, mode="in")), 'feature'=degree(graph, mode="in"))
    centr_degree(graph, mode="in", normalized=T)
    print(ggplot(degree_table, aes(reorder(nodes,feature), feature)) + geom_bar(stat="identity") +
            theme_classic()+
            theme(axis.text.x = element_text(angle = 45,hjust=1, size=cex_axis), plot.margin = margin(0.8, 0.8, 0.8, 0.8, "cm"))+
            # ggtitle("Transitivity") + 
            # xlab("Diseases") + 
            xlab("Diseases") + ylab("Node degree"))
    
    # Closeness (centrality based on distance to others in the graph)
    #Inverse of the node’s average geodesic distance to others in the network.
    closeness(graph, mode="all", weights=NA) 
    mean_closeness = mean(closeness(graph, mode="all", weights=NA) )
    centr_clo(graph, mode="all", normalized=T)
    degree_table <- data.frame('nodes'=names(closeness(graph, mode="all", weights=NA)), 'feature'=closeness(graph, mode="all", weights=NA))
    print(ggplot(degree_table, aes(reorder(nodes,feature), feature)) + geom_bar(stat="identity") +
            theme_classic()+
            theme(axis.text.x = element_text(angle = 45,hjust=1, size=cex_axis), plot.margin = margin(0.8, 0.8, 0.8, 0.8, "cm"))+
            # ggtitle("Transitivity") + 
            # xlab("Diseases") + 
            xlab("Diseases") + ylab("Closeness"))
    
    # Eigenvector (centrality proportional to the sum of connection centralities)
    # Values of the first eigenvector of the graph matrix.
    
    eigen_centrality(graph, directed=T, weights=NA)
    centr_eigen(graph, directed=T, normalized=T) 
    
    # Betweenness (centrality based on a broker position connecting others)
    # Number of geodesics that pass through the node or the edge.
    betweenness <- betweenness(distance_graph, directed=FALSE)
    betweenness_wo_weights <- betweenness(graph, directed=FALSE, weights=NA)
    edge_betweenness(graph, directed=T, weights=NA)
    centr_betw(distance_graph, directed=T, normalized=T)
    betweenness <- data.frame('nodes'=names(betweenness), 'feature'=betweenness)
    print(ggplot(betweenness, aes(reorder(nodes,feature), feature)) + geom_bar(stat="identity") +
            theme_classic()+
            theme(axis.text.x = element_text(angle = 45,hjust=1, size=cex_axis), plot.margin = margin(0.8, 0.8, 0.8, 0.8, "cm"))+
            ggtitle("Betweenness with weights") +
            # xlab("Diseases") + 
            xlab("Diseases") + ylab("Betweenness"))
    
    betweenness <- data.frame('nodes'=names(betweenness_wo_weights), 'feature'=betweenness_wo_weights)
    print(ggplot(betweenness, aes(reorder(nodes,feature), feature)) + geom_bar(stat="identity") +
            theme_classic()+
            theme(axis.text.x = element_text(angle = 45,hjust=1, size=cex_axis), plot.margin = margin(0.8, 0.8, 0.8, 0.8, "cm"))+
            ggtitle("Betweenness wo weights") +
            # xlab("Diseases") + 
            xlab("Diseases") + ylab("Betweenness"))
    
    mean_betweenness = mean(betweenness(distance_graph, directed=FALSE))
    mean_betweenness_wo_weights <- mean(betweenness(graph, directed=FALSE, weights=NA))
    
    ### HUBS AND AUTHORITIES
    # HUBS with weights
    hs <- hub_score(graph, weights=abs(E(graph)$weight))$vector
    hubs_table <- data.frame('nodes'=names(hs), 'feature'=hs)
    print(ggplot(hubs_table, aes(reorder(nodes,feature), feature)) + geom_bar(stat="identity") +
            theme_classic()+
            theme(axis.text.x = element_text(angle = 45,hjust=1, size=cex_axis), plot.margin = margin(0.8, 0.8, 0.8, 0.8, "cm"))+
            ggtitle("Hubs with weights") +
            # xlab("Diseases") + 
            xlab("Diseases") + ylab("Hubs"))
    # wo weights
    hs <- hub_score(graph, weights=NA)$vector
    hubs_table <- data.frame('nodes'=names(hs), 'feature'=hs)
    print(ggplot(hubs_table, aes(reorder(nodes,feature), feature)) + geom_bar(stat="identity") +
            theme_classic()+
            theme(axis.text.x = element_text(angle = 45,hjust=1, size=cex_axis), plot.margin = margin(0.8, 0.8, 0.8, 0.8, "cm"))+
            ggtitle("Hubs wo weights") +
            # xlab("Diseases") + 
            xlab("Diseases") + ylab("Hubs"))
    
    # AUTHORITIES with weights
    as <- authority_score(graph, weights=abs(E(graph)$weight))$vector
    auth_table <- data.frame('nodes'=names(as), 'feature'=as)
    print(ggplot(auth_table, aes(reorder(nodes,feature), feature)) + geom_bar(stat="identity") +
            theme_classic()+
            theme(axis.text.x = element_text(angle = 45,hjust=1, size=cex_axis), plot.margin = margin(0.8, 0.8, 0.8, 0.8, "cm"))+
            ggtitle("Authorities with weights") +
            # xlab("Diseases") + 
            xlab("Diseases") + ylab("Authority"))
    
    # wo weights
    as <- authority_score(graph, weights=NA)$vector
    auth_table <- data.frame('nodes'=names(as), 'feature'=as)
    print(ggplot(auth_table, aes(reorder(nodes,feature), feature)) + geom_bar(stat="identity") +
            theme_classic()+
            theme(axis.text.x = element_text(angle = 45,hjust=1, size=cex_axis), plot.margin = margin(0.8, 0.8, 0.8, 0.8, "cm"))+
            ggtitle("Authorities wo weights") +
            # xlab("Diseases") + 
            xlab("Diseases") + ylab("Authority"))
    
    par(mfrow=c(1,2))
    plot(graph, vertex.size=hs*20, main="Hubs")
    plot(graph, vertex.size=as*10, main="Authorities")
    
    ### DISTANCES AND PATHS
    # Average path length: the mean of the shortest distance between 
    #each pair of nodes in the network (in both directions for directed graphs).
    
    mean_distance <- mean_distance(distance_graph, directed=F)
    mean_distance_wo_weights <- mean_distance(distance_graph, directed=F, weights=NA)
    # mean_distance <- mean_distance(graph, directed=F)
    # distances(graph) # with edge weights # YIELDS AND ERROR
    distances(graph, weights=NA) # ignore weights
    
    ### FINDING COMMUNITIES
    # cliques: complete subgraphs of an undirected graph
    if(obtain_cliques == TRUE){
      graph.sym <- as.undirected(graph, mode= "collapse",
                                 edge.attr.comb=list(weight="sum", "ignore"))
      cliques(graph) # list of cliques       
      sapply(cliques(graph), length) # clique sizes
      saveRDS(largest_cliques(graph), file=paste("Network_building/Network_analysis/Topology/Cliques/",output_filename,k,'.rds', sep=""))
    }
    dev.off()
    
    properties = c(df_name, gsub("_","",surname),
                   nint_total, nint_pos, nint_neg, (nint_posi / nint_total)*100, (nint_neg / nint_total)*100,
                    n_connected_components, size_connected_components, mean_deg, edge_density, transitivity, diameter, diameter_wo_weights, 
                    mean_closeness, mean_betweenness, mean_betweenness_wo_weights, mean_distance, mean_distance_wo_weights)
    names(properties) = c("net_name","Int_type",
                          "nint_total", "nint_pos", "nint_neg", "perc_pos", "perc_neg",
                          "n_connected_components", "size_connected_components", "mean_deg", "edge_density", "transitivity", "diameter", "diameter_wo_weights", 
                          "mean_closeness", "mean_betweenness", "mean_betweenness_wo_weights", "mean_distance", "mean_distance_wo_weights")
    
    properties_df <- rbind(properties_df, properties)
    if(!"mean_closeness" %in% colnames(properties_df)){
      colnames(properties_df) = names(properties)
    }
    
    

  }
} 

fprop_df = properties_df

# Identify numeric columns
num_cols <- sapply(fprop_df, function(x) all(is.na(as.numeric(x))) == FALSE)

# Print column names of numeric columns
colnames(fprop_df)[num_cols]; length(colnames(fprop_df)[num_cols])
# ncol(fprop_df)

# Convert numeric columns to numeric
fprop_df[, num_cols] <- sapply(fprop_df[, num_cols], as.numeric)
str(fprop_df)

# Change closeness to scientific notation
# fprop_df$mean_closeness = format(fprop_df$mean_closeness, scientific = TRUE, digits = 2)
fprop_df$mean_closeness = sprintf("%.2e", fprop_df$mean_closeness)

# Remove 'mean_closeness' fromthe list of numeric columns to round
num_cols["mean_closeness"] = FALSE
# num_cols = setdiff(num_cols, which(names(fprop_df) == "mean_closeness"))

# Round numeric columns
fprop_df[, num_cols] <- round(fprop_df[, num_cols], 2)
str(fprop_df)

head(fprop_df)



fprop_df = t(fprop_df)
write.table(fprop_df,file=paste("Network_building/Network_analysis/Topology/final_networks_topology.txt",sep=""),sep="\t",row.names=T, col.names=F, quote=FALSE)


# RUN THE TOPOLOGICAL ANALYSIS FOR THE NETWORK (k=1), POS (k=2) AND NEG SUBNETWORK (k=3)
dfpos <- df[df$Distance >= 0, ]
dfneg <- df[df$Distance < 0, ]; dfneg$Distance <- abs(dfneg$Distance)
graph_list <- list(df,dfpos,dfneg)

k <- 1
for(df in graph_list){
  # head(df)
  # class(df)
  # k = 2
  pdf(paste("Network_building/Network_analysis/Topology/",output_filename,k,"_topology_plots.pdf",sep=""))
  graph <- graph_from_data_frame(df, directed=FALSE)
  E(graph)$weight <- abs(df$Distance) 
  epsilon = 0.00000001 # 1*10-8
  distance_graph = graph
  if(output_filename == "barabasi"){
    max_barabasi = max(E(graph)$weight)
    E(distance_graph)$weight = (1 - (E(graph)$weight/E(graph)$weight) + epsilon)
  }else{
    E(distance_graph)$weight = (1 - E(graph)$weight + epsilon)
  }
  
  
  connected_components <- clusters(graph)
  n_connected_components <- connected_components$no
  size_connected_components <- connected_components$csize
  if(length(size_connected_components) > 1){
    size_connected_components <- paste(size_connected_components,collapse="_")
  }
  
  deg <- degree(graph, mode="all")
  mean_deg <- mean(deg) 
  
  #Density of a graph: The proportion of present edges from all possible edges in the network
  edge_density <- edge_density(graph, loops=F)
  edge_density
  
  # Reciprocity: it only makes sense for directed graphs
  reciprocity(graph)
  
  # dyad_census(graph) # Mutual, asymmetric, and nyll node pairs
  # 2*dyad_census(graph)$mut/ecount(graph) # Calculating reciprocity
  
  # Transitivity: Transitivity measures the probability that the adjacent vertices of a vertex are connected
  transitivity <- transitivity(graph, type="global") # fot the network (global)
  local_trans <- transitivity(graph, vids=V(graph), type="local") # for each of the nodes (local)
  transitivity_table <- data.frame(nodes=names(V(graph)), transitivity=local_trans)
  print(ggplot(transitivity_table, aes(reorder(nodes,transitivity), transitivity)) + geom_bar(stat="identity") +
          theme_classic()+
          theme(axis.text.x = element_text(angle = 45,hjust=1, size=cex_axis),plot.margin = margin(0.8, 0.8, 0.8, 0.8, "cm"))+
          # ggtitle("Transitivity") + 
          # xlab("Diseases") + 
          ylab("Transitivity"))
  
  ### Diameter: longest geodesic distance (length of the shortest path between two nodes) in the network
  # With weights
  diameter <- diameter(distance_graph, directed=F)
  get_diameter(distance_graph, directed=F)
  # Without weights
  diameter_wo_weights <- diameter(graph, directed=F, weights=NA)
  get_diameter(graph, directed=F, weights=NA)
  
  ###Degree
  deg <- degree(graph, mode="all")
  plot(graph, vertex.size=deg*0.4,layout=layout.circle)
  hist(deg, breaks=1:vcount(graph)-1, main="Histogram of node degree")
  
  # Degree distribution
  deg.dist <- degree_distribution(graph, cumulative=T, mode="all")
  
  plot( x=0:max(deg), y=1-deg.dist, pch=19, cex=1.2, col="orange", 
        xlab="Degree", ylab="Cumulative Frequency")
  
  ### CENTRALITY AND CENTRALIZATION
  # Degree (number of ties) ### hacer un plot de esto con ggplot quizá!
  degree_table <- data.frame('nodes'=names(degree(graph, mode="in")), 'feature'=degree(graph, mode="in"))
  centr_degree(graph, mode="in", normalized=T)
  print(ggplot(degree_table, aes(reorder(nodes,feature), feature)) + geom_bar(stat="identity") +
          theme_classic()+
          theme(axis.text.x = element_text(angle = 45,hjust=1, size=cex_axis), plot.margin = margin(0.8, 0.8, 0.8, 0.8, "cm"))+
          # ggtitle("Transitivity") + 
          # xlab("Diseases") + 
          xlab("Diseases") + ylab("Node degree"))
  
  # Closeness (centrality based on distance to others in the graph)
  #Inverse of the node’s average geodesic distance to others in the network.
  closeness(graph, mode="all", weights=NA) 
  centr_clo(graph, mode="all", normalized=T)
  degree_table <- data.frame('nodes'=names(closeness(graph, mode="all", weights=NA)), 'feature'=closeness(graph, mode="all", weights=NA))
  print(ggplot(degree_table, aes(reorder(nodes,feature), feature)) + geom_bar(stat="identity") +
          theme_classic()+
          theme(axis.text.x = element_text(angle = 45,hjust=1, size=cex_axis), plot.margin = margin(0.8, 0.8, 0.8, 0.8, "cm"))+
          # ggtitle("Transitivity") + 
          # xlab("Diseases") + 
          xlab("Diseases") + ylab("Closeness"))
  
  # Eigenvector (centrality proportional to the sum of connection centralities)
  # Values of the first eigenvector of the graph matrix.
  
  eigen_centrality(graph, directed=T, weights=NA)
  centr_eigen(graph, directed=T, normalized=T) 
  
  # Betweenness (centrality based on a broker position connecting others)
  # Number of geodesics that pass through the node or the edge.
  betweenness <- betweenness(distance_graph, directed=FALSE)
  betweenness_wo_weights <- betweenness(graph, directed=FALSE, weights=NA)
  edge_betweenness(graph, directed=T, weights=NA)
  centr_betw(distance_graph, directed=T, normalized=T)
  betweenness <- data.frame('nodes'=names(betweenness), 'feature'=betweenness)
  print(ggplot(betweenness, aes(reorder(nodes,feature), feature)) + geom_bar(stat="identity") +
          theme_classic()+
          theme(axis.text.x = element_text(angle = 45,hjust=1, size=cex_axis), plot.margin = margin(0.8, 0.8, 0.8, 0.8, "cm"))+
          ggtitle("Betweenness with weights") +
          # xlab("Diseases") + 
          xlab("Diseases") + ylab("Betweenness"))
  
  betweenness <- data.frame('nodes'=names(betweenness_wo_weights), 'feature'=betweenness_wo_weights)
  print(ggplot(betweenness, aes(reorder(nodes,feature), feature)) + geom_bar(stat="identity") +
          theme_classic()+
          theme(axis.text.x = element_text(angle = 45,hjust=1, size=cex_axis), plot.margin = margin(0.8, 0.8, 0.8, 0.8, "cm"))+
          ggtitle("Betweenness wo weights") +
          # xlab("Diseases") + 
          xlab("Diseases") + ylab("Betweenness"))
  
  ### HUBS AND AUTHORITIES
  # HUBS with weights
  hs <- hub_score(graph, weights=abs(E(graph)$weight))$vector
  hubs_table <- data.frame('nodes'=names(hs), 'feature'=hs)
  print(ggplot(hubs_table, aes(reorder(nodes,feature), feature)) + geom_bar(stat="identity") +
          theme_classic()+
          theme(axis.text.x = element_text(angle = 45,hjust=1, size=cex_axis), plot.margin = margin(0.8, 0.8, 0.8, 0.8, "cm"))+
          ggtitle("Hubs with weights") +
          # xlab("Diseases") + 
          xlab("Diseases") + ylab("Hubs"))
  # wo weights
  hs <- hub_score(graph, weights=NA)$vector
  hubs_table <- data.frame('nodes'=names(hs), 'feature'=hs)
  print(ggplot(hubs_table, aes(reorder(nodes,feature), feature)) + geom_bar(stat="identity") +
          theme_classic()+
          theme(axis.text.x = element_text(angle = 45,hjust=1, size=cex_axis), plot.margin = margin(0.8, 0.8, 0.8, 0.8, "cm"))+
          ggtitle("Hubs wo weights") +
          # xlab("Diseases") + 
          xlab("Diseases") + ylab("Hubs"))
  
  # AUTHORITIES with weights
  as <- authority_score(graph, weights=abs(E(graph)$weight))$vector
  auth_table <- data.frame('nodes'=names(as), 'feature'=as)
  print(ggplot(auth_table, aes(reorder(nodes,feature), feature)) + geom_bar(stat="identity") +
          theme_classic()+
          theme(axis.text.x = element_text(angle = 45,hjust=1, size=cex_axis), plot.margin = margin(0.8, 0.8, 0.8, 0.8, "cm"))+
          ggtitle("Authorities with weights") +
          # xlab("Diseases") + 
          xlab("Diseases") + ylab("Authority"))
  
  # wo weights
  as <- authority_score(graph, weights=NA)$vector
  auth_table <- data.frame('nodes'=names(as), 'feature'=as)
  print(ggplot(auth_table, aes(reorder(nodes,feature), feature)) + geom_bar(stat="identity") +
          theme_classic()+
          theme(axis.text.x = element_text(angle = 45,hjust=1, size=cex_axis), plot.margin = margin(0.8, 0.8, 0.8, 0.8, "cm"))+
          ggtitle("Authorities wo weights") +
          # xlab("Diseases") + 
          xlab("Diseases") + ylab("Authority"))
  
  par(mfrow=c(1,2))
  plot(graph, vertex.size=hs*20, main="Hubs")
  plot(graph, vertex.size=as*10, main="Authorities")
  
  ### DISTANCES AND PATHS
  # Average path length: the mean of the shortest distance between 
  #each pair of nodes in the network (in both directions for directed graphs).
  
  mean_distance <- mean_distance(distance_graph, directed=F)
  # mean_distance <- mean_distance(graph, directed=F)
  # distances(graph) # with edge weights # YIELDS AND ERROR
  distances(graph, weights=NA) # ignore weights
  
  ### FINDING COMMUNITIES
  # cliques: complete subgraphs of an undirected graph
  if(obtain_cliques == TRUE){
    graph.sym <- as.undirected(graph, mode= "collapse",
                               edge.attr.comb=list(weight="sum", "ignore"))
    cliques(graph) # list of cliques       
    sapply(cliques(graph), length) # clique sizes
    saveRDS(largest_cliques(graph), file=paste("Network_building/Network_analysis/Topology/Cliques/",output_filename,k,'.rds', sep=""))
  }
  
  
  # # Community detection algorithms
  # if (input$comm_algorithm == 'greedy'){
  #   # Using greedy optimization of modularity
  #   fc <- fastgreedy.community(graph)
  #   V(graph)$community <- fc$membership
  # }else if(input$comm_algorithm == 'rand_walks'){
  #   # Using random walks
  #   fc <- cluster_walktrap(graph)
  #   V(graph)$community <- fc$membership #membership(fc)
  # }
  
  # K-core decomposition
  # The k-core is the maximal subgraph in which every node has degree of 
  # at least k. The result here gives the coreness of each vertex in the 
  # network. A node has coreness D if it belongs to a D-core but not 
  # to (D+1)-core.
  kc <- coreness(graph, mode="all")
  # plot(graph, vertex.size=kc*6, vertex.label=kc, vertex.color=colrs[kc]) # FIXX
  kc_table <- data.frame('nodes'=names(kc), 'feature'=kc)
  print(ggplot(kc_table, aes(reorder(nodes,feature), feature)) + geom_bar(stat="identity") +
          theme_classic()+
          theme(axis.text.x = element_text(angle = 45,hjust=1, size=cex_axis), plot.margin = margin(0.8, 0.8, 0.8, 0.8, "cm"))+
          # ggtitle("Authorities with weights") +
          # xlab("Diseases") + 
          xlab("Diseases") + ylab("K-core decomposition"))
  
  ### ASSORTATIVITY AND HOMOPHILY
  # assortativity_nominal(graph, V(graph)$media.type, directed=F)
  assortativity_degree <- assortativity_degree(graph, directed=F)
  
  cfeat <- c('N conn comp','Size conn comp','Mean deg',
             'edge_density','transitivity','diameter',
             'diameter_wo_weights','mean_distance','assortativity_degree') 
  cvals <- c(n_connected_components,size_connected_components,mean_deg,
             edge_density,transitivity,diameter,
             diameter_wo_weights,mean_distance,assortativity_degree)
  cdf <- data.frame(Feature=cfeat,Value = cvals)
  
  dev.off()
  write.table(cdf,file=paste("Network_building/Network_analysis/Topology/",output_filename,k,"_topology.txt",sep=""),sep="\t",row.names=F, quote=FALSE)
  k <- k+1
}
