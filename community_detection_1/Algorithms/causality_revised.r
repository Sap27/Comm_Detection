install.packages("igraph", repos = "https://cloud.r-project.org/")
install.packages("jsonlite", repos = "https://cloud.r-project.org/")
library(igraph)
library(jsonlite)
# Function to rescale weights to integer values
rescale <- function(x, graph, maximum, minimum) {
  numV <- vcount(graph)
  return((numV - 1) / (maximum - minimum) * (x - minimum) + 1)
}

# Main function that orchestrates the entire process
main <- function(json_data, b = 1.3, c = 500, i = 1, filter, threshold = 1.8, inteWeight = "yes", weighted = TRUE, dir, post, smallest, largest, b2, c2, i2,script_dir) {
  #input <- read.table(file, sep = '\t')
  #print('hi')
  setwd(script_dir)
  input<-fromJSON(json_data)
  #print('bye')
  # Step 1: Preprocess the input data to create a graph
  graph <- preProcessing(input, filter, threshold, integerWeight = inteWeight)
  
  # Step 2: Generate a file formatted for the MLRMCL algorithm
  temp_dir <- tempdir()
  windows_test_r_file <- file.path(temp_dir, "test_r.txt")
  windows_output_file <- file.path(temp_dir, "output.txt")
  print(windows_output_file)
  # Convert Windows paths to WSL paths
  wsl_temp_dir <- system(paste("wsl wslpath -a '", temp_dir, "'", sep=""), intern = TRUE)
  wsl_test_r_file <- file.path(wsl_temp_dir, "test_r.txt")
  wsl_output_file <- file.path(wsl_temp_dir, "output.txt")
  
  generateFile(graph, weighted, windows_test_r_file)
  
  # Step 3: Run the MLRMCL algorithm on the graph
  system(paste("wsl ./mlrmcl -b", b, "-c", c, "-i", i, "-o", wsl_output_file, wsl_test_r_file, sep = " "))
  
  # Step 4: Post-process the MLRMCL output to refine the clusters
  output <- postProcessing(post, smallest, largest, graph, b2, c2, i2, windows_output_file)
  
  # Clean up temporary files
  file.remove(windows_test_r_file)
  file.remove(windows_output_file)
  
  return(output)
}

# Function to preprocess the input data and create a graph
preProcessing <- function(input, method = c("quantile", "pageRank", "double"), i, integerWeight = c("yes", "no", "1")) {
  input[, 1:2] <- input[, 1:2] + 1
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
generateFile <- function(graph, weighted = TRUE, file_path) {
  adjacency_matrix <- get.adjacency(graph, attr = "weight")
  m <- as.matrix(adjacency_matrix)
  edge_indices <- apply(m != 0, 1, which, arr.ind = TRUE)
  
  if (file.exists(file_path)) file.remove(file_path)
  write(c(vcount(graph), ecount(graph), ifelse(weighted, 1, 0)), file = file_path, sep = "\t")
  
  if (weighted) {
    weights <- sapply(1:nrow(m), function(i) m[i, edge_indices[[i]]])
    edge_list <- lapply(1:length(edge_indices), function(i) c(rbind(as.character(edge_indices[[i]]), weights[[i]])))
    lapply(edge_list, write, file_path, ncolumns = length(edge_list), append = TRUE, sep = "\t")
  } else {
    lapply(edge_indices, write, file_path, ncolumns = length(edge_indices), append = TRUE, sep = "\t")
  }
}

# Function to process the output from MLRMCL and refine clusters
postProcessing <- function(method = c("random", "discard", "recluster"), smallest = 3, largest = 100, graph, b2, c2, i2, output_file) {
  
  result <- read.table(output_file)
  unique_clusters <- unique(result)
  output <- lapply(0:(nrow(unique_clusters) - 1), function(i) which(result == i) - 1)
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

# Define the parameters
file <- "C:/Users/hp/Social_networks/football.txt"



b <- 1.3
#c <- 1000
i <-2
#filter <- "pageRank"
threshold <- 2
#inteWeight <- "no"
#weighted <- TRUE
dir <- "output_clusters"
post <- "random"
smallest <- 3
#largest <- 100
b2 <- 2
c2 <- 5000
i2 <- 2

# Call the main function
#print('hi')
args<-commandArgs(trailingOnly=TRUE)
json_data<-args[[1]]
out_path<-args[[2]]
filter<-args[[3]]
inteWeight <- args[[4]]
weighted <- as.logical(args[[5]])
#smallest<-as.logical(args[[6]])
largest<-as.numeric(args[[6]])
c<-as.numeric(args[[7]])
script_dir<-args[[8]]
communities <- main(json_data, b, c, i, filter, threshold, inteWeight, weighted, dir, post, smallest, largest, b2, c2, i2,script_dir)
#print(output)
#print(length(communities))
json_communities<-toJSON(communities,auto_unbox=TRUE)
cat(json_communities)