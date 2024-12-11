# Install required packages if not already installed
suppressMessages(suppressWarnings({
  options(repos = c(CRAN = "https://cran.rstudio.com/"))
  if (!require("igraph")) install.packages("igraph")
  if (!require("WGCNA")) install.packages("WGCNA")
  if (!require("flashClust")) install.packages("flashClust")
  if (!require("ComplexHeatmap")) BiocManager::install("ComplexHeatmap")
  
  library(WGCNA)
  library(flashClust)
  library(igraph)
  library(ComplexHeatmap)
  library(jsonlite)
}))

# Function to normalize and symmetrize adjacency matrix
normalize_adj <- function(adj) {
  aux.adj <- matrix(0, nrow = nrow(adj), ncol = ncol(adj))
  for (i in 1:nrow(adj)) {
    for (j in 1:nrow(adj)) {
      aux.adj[j, i] <- (adj[j, i] + adj[i, j]) / 2
    }
  }
  adj <- aux.adj
  if (max(adj) > 1) {
    adj <- (adj - min(adj)) / (max(adj) - min(adj))
  }
  return(adj)
}

# Function to calculate entropy and refine communities
func_entr <- function(cor.data) {
  t_thr <- seq(0.01, 0.8, by = 0.01)
  ent2 <- c()
  
  for (i in seq_along(t_thr)) {
    data.thr <- ifelse(cor.data < t_thr[i], 0, 1)
    g.thr <- graph.adjacency(data.thr, mode = 'undirected', diag = FALSE)
    g.thr.bet <- betweenness(g.thr, normalized = TRUE)
    
    ent2[i] <- ifelse(max(g.thr.bet) == 0, 0, 
                      -sum(g.thr.bet[g.thr.bet > 0] * log2(g.thr.bet[g.thr.bet > 0])))
  }
  
  max_idx <- which.max(ent2)
  th <- t_thr[max_idx]
  data.thr1 <- ifelse(cor.data < th, 0, 1)
  g.thr1 <- graph.adjacency(data.thr1, mode = 'undirected', diag = FALSE)
  
  return(g.thr1)
}

# Main process function
process_communities <- function(dd, bool_value, min_limit) {
  if (min(dd[, 1]) == 0) {
    dd[, c(1, 2)] <- dd[, c(1, 2)] + 1
  }
  
  g <- graph.data.frame(dd, directed = bool_value)
  E(g)$weight <- as.numeric(dd[[3]])
  adj <- normalize_adj(as.matrix(get.adjacency(g, attr = 'weight')))
  
  # Calculate TOM and perform hierarchical clustering
  TOM <- TOMsimilarity(adj)
  dissTOM <- 1 - TOM
  dend.tree <- hclust(as.dist(dissTOM), method = "average")
  
  unmergedLabels <- cutreeDynamic(dendro = dend.tree, distM = dissTOM, deepSplit = 4, 
                                  pamRespectsDendro = FALSE, minClusterSize = min_limit)
  unmergedLabels <- as.numeric(unmergedLabels[unmergedLabels > 0])
  
  my.comms <- list()
  
  # Split communities further using fast-greedy algorithm
  for (i in seq_along(unique(unmergedLabels))) {
    adj.aux <- adj[which(unmergedLabels == i), which(unmergedLabels == i)]
    g.aux <- graph.adjacency(adj.aux, mode = 'undirected', weighted = TRUE, diag = FALSE)
    fc <- fastgreedy.community(g.aux)
    ind1 <- which(unmergedLabels == i)
    
    for (j in seq_along(unique(membership(fc)))) {
      ind2 <- which(membership(fc) == j)
      aux <- ind1[ind2]
      if (length(aux) > 2) {
        my.comms <- c(my.comms, list(aux))
      }
    }
  }
  
  # Refine communities using entropy-based function
  my.comms_final <- lapply(my.comms, function(community) {
    if (length(community) > 100) {
      g.aux1 <- func_entr(adj[community, community])
      fc1 <- fastgreedy.community(g.aux1)
      lapply(unique(membership(fc1)), function(t) {
        aux <- community[which(membership(fc1) == t)]
        if (length(aux) > 2) return(aux)
      }) %>% unlist(recursive = FALSE)
    } else {
      return(list(community))
    }
  }) %>% unlist(recursive = FALSE)
  
  my.comms_final <- unique(my.comms_final)
  
  # Further refinement and filtering
  my.comms_final <- lapply(my.comms_final, function(comm) comm - 1)
  my.comms_final <- my.comms_final[lapply(my.comms_final, length) > 2]
  
  return(my.comms_final)
}

# Main script execution
args <- commandArgs(trailingOnly = TRUE)
json_data <- args[[1]]
output_file <- args[[2]]
bool_value <- as.logical(args[[3]])
min_limit <- as.numeric(args[[4]])

dd <- fromJSON(json_data)
communities <- process_communities(dd, bool_value, min_limit)

# Output results
json_communities <- toJSON(communities, auto_unbox = TRUE)
cat(json_communities)

if (output_file != '') {
  write(communities, output_file, ncolumns = 102, append = TRUE, sep = "\t")
}
