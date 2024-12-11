library(igraph)

# Function to rescale weights to integer values
rescale <- function(x, graph, maximum, minimum) {
  numV <- vcount(graph)
  return((numV - 1) / (maximum - minimum) * (x - minimum) + 1)
}

# Main function that orchestrates the entire process
main <- function(file, b = 1.3, c = 500, i = 1, filter, threshold = 1.8, inteWeight = "yes", weighted = TRUE, dir, post, smallest, largest, b2, c2, i2) {
  input <- read.table(file, sep = '\t')
  
  # Step 1: Preprocess the input data to create a graph
  graph <- preProcessing(input, filter, threshold, integerWeight = inteWeight)
  
  # Step 2: Generate a file formatted for the MLRMCL algorithm
  generateFile(graph, weighted)
  
  # Step 3: Run the MLRMCL algorithm on the graph
  system(paste("wsl ./community_detection/mlrmcl -b ", b, " -c ", c, " -i ", i, " -o output.txt test_r.txt", sep = ""))
  
  # Step 4: Post-process the MLRMCL output to refine the clusters
  output <- postProcessing(post, smallest, largest, graph, b2, c2, i2)
  print(output)
  # Step 5: Write the final clusters to a file
  #writeFile(output, file, dir)
  
  # Clean up temporary files
  #file.remove("output.txt")
  #file.remove("test.txt")
  
  return(output)
}

# Function to preprocess the input data and create a graph
preProcessing <- function(input, method = c("quantile", "pageRank", "double"), i, integerWeight = c("yes", "no", "1")) {
  # Adjust node numbering to start from 1
  input[, 1:2] <- input[, 1:2] + 1
  
  # Filtering based on method
  if (method == "quantile") {
    filter <- input[, 3] >= quantile(input[, 3])[i]
    input <- input[filter, ]
  } else if (method == "pageRank" || method == "double") {
    graph <- make_graph(c(t(input[, 1:2])), directed = FALSE)
    pageRank <- page_rank(graph, directed = FALSE, weights = input[, 3])$vector
    index <- which(pageRank > quantile(pageRank)[i])
    input <- input[input[, 3] >= quantile(input[, 3])[i], ]
    
    if (method == "double") {
      filter <- input[, 3] >= quantile(input[, 3])[i]
      input <- input[filter, ]
    }
  }
  
  graph <- make_graph(c(t(input[, 1:2])), directed = FALSE)
  
  # Scale weights if required
  if (integerWeight == "yes") {
    maximum <- max(input[, 3])
    minimum <- min(input[, 3])
    w <- as.integer(rescale(input[, 3], graph, maximum, minimum))
    E(graph)$weight <- as.numeric(w)
  } else if (integerWeight == "no") {
    E(graph)$weight <- as.numeric(input[, 3])
  } else {
    E(graph)$weight <- 1
  }
  
  return(graph)
}

# Function to generate the input file for MLRMCL from the graph
generateFile <- function(graph, weighted = TRUE) {
  adjacency_matrix <- get.adjacency(graph, attr = "weight")
  m <- as.matrix(adjacency_matrix)
  edge_indices <- apply(m != 0, 1, which, arr.ind = TRUE)
  
  fn <- "test_r.txt"
  if (file.exists(fn)) file.remove(fn)
  
  # Write the graph details to a file
  write(c(vcount(graph), ecount(graph), ifelse(weighted, 1, 0)), file = "test_r.txt", sep = "\t")
  
  if (weighted) {
    weights <- sapply(1:nrow(m), function(i) m[i, edge_indices[[i]]])
    edge_list <- lapply(1:length(edge_indices), function(i) c(rbind(as.character(edge_indices[[i]]), weights[[i]])))
    lapply(edge_list, write, "test_r.txt", ncolumns = length(edge_list), append = TRUE, sep = "\t")
  } else {
    lapply(edge_indices, write, "test_r.txt", ncolumns = length(edge_indices), append = TRUE, sep = "\t")
  }
}

# Function to process the output from MLRMCL and refine clusters
postProcessing <- function(method = c("random", "discard", "recluster"), smallest = 3, largest = 100, graph, b2, c2, i2) {
  result <- read.table("output.txt")
  unique_clusters <- unique(result)
  output <- lapply(0:(nrow(unique_clusters) - 1), function(i) which(result == i) - 1)
  
  # Remove small clusters
  output <- output[sapply(output, length) >= smallest]
  
  if (method == "random") {
    large_clusters <- output[sapply(output, length) > largest]
    small_clusters <- output[sapply(output, length) <= largest]
    for (large in large_clusters) {
      small_clusters <- append(small_clusters, split(large, 1:ceiling(length(large) / largest)))
    }
    output <- small_clusters
  } else if (method == "discard") {
    output <- output[sapply(output, length) <= largest]
  } else if (method == "recluster") {
    output <- reclusterLargeClusters(output, graph, b2, c2, i2, smallest, largest)
  }
  
  return(output)
}

# Helper function to recluster large clusters recursively
reclusterLargeClusters <- function(clusters, graph, b2, c2, i2, smallest, largest) {
  large_clusters <- clusters[sapply(clusters, length) > largest]
  small_clusters <- clusters[sapply(clusters, length) <= largest]
  
  if (length(large_clusters) == 0) {
    return(clusters)
  }
  
  result <- list()
  for (large in large_clusters) {
    subgraph_indices <- unlist(large) + 1
    subgraph <- induced_subgraph(graph, subgraph_indices, "copy_and_delete")
    V(subgraph)$name <- as.character(subgraph_indices)
    
    generateFile(subgraph, weighted = TRUE)
    system(paste("./mlrmcl -b ", b2, " -c ", c2, " -i ", i2, " -o output2.txt test_r.txt", sep = ""))
    
    subgraph_result <- read.table("output2.txt")
    subgraph_clusters <- lapply(1:nrow(unique(subgraph_result)), function(i) as.integer(V(subgraph)$name[which(subgraph_result == i)]) - 1)
    subgraph_clusters <- subgraph_clusters[sapply(subgraph_clusters, length) > smallest]
    
    small_clusters <- append(small_clusters, subgraph_clusters[sapply(subgraph_clusters, length) <= largest])
    result <- append(result, subgraph_clusters[sapply(subgraph_clusters, length) > largest])
  }
  
  return(append(small_clusters, reclusterLargeClusters(result, graph, b2, c2, i2, smallest, largest)))
}

# Function to write the final clusters to a file
writeFile <- function(output, file, dir) {
  dir.create(dir, showWarnings = FALSE)
  fn <- paste(dir, "/", file, sep = "")
  if (file.exists(fn)) file.remove(fn)
  
  for (i in 1:length(output)) {
    cluster <- c(i, 0.5, output[[i]])
    write(cluster, file = fn, ncolumns = length(cluster), append = TRUE, sep = "\t")
  }
  
  # Clean up temporary files
  file.remove("output.txt")
  if (file.exists("output2.txt")) file.remove("output2.txt")
  #file.remove("test_r.txt")
}

# Define the parameters
file <- "C:/Users/hp/community_detection/network.dat"
b <- 1.3
c <- 500
i <- 1.8
filter <- "pageRank"  # Use "pageRank" or "double" as needed
threshold <- 2
inteWeight <- "yes"
weighted <- TRUE
dir <- "output_clusters"
post <- "random"  # Use "discard" or "recluster" as needed
smallest <- 3
largest <- 100
b2 <- 2
c2 <- 5000
i2 <- 2

# Call the main function
output <- main(file, b, c, i, filter, threshold, inteWeight, weighted, dir, post, smallest, largest, b2, c2, i2)