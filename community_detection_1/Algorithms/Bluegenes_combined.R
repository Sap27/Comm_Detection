#dir.create("raw_data")
#install.packages("igraph")
install.packages("igraph", repos = "https://cloud.r-project.org/")

options(repos = c(CRAN = "https://cran.rstudio.com/"))
#source("biocLite.R")
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(version = "3.20")
BiocManager::install(c("AnnotationDbi", "impute", "GO.db", "preprocessCore", "ComplexHeatmap"))
install.packages("WGCNA") 
install.packages("flashClust")
suppressMessages(suppressWarnings({
  library(WGCNA)
  library(flashClust)
  # Any other code where output needs to be suppressed

library(igraph)
library(WGCNA)
library(flashClust)
library(ComplexHeatmap)
args<- commandArgs(trailingOnly=TRUE)
json_data <- args[[1]]
library(jsonlite)
# Read JSON data from the file
temp <- fromJSON(json_data)
out_path <-args[[2]]
minlimit <- as.numeric(args[[3]])
bool_value<-as.logical(args[[7]])
cutsize<-as.numeric(args[[6]])
#in_file<-'C:/Users/hp/Mu_0.80/network.dat'
#out_path<-'C:/Users/hp/Mu_0.80/'
#minlimit<-21
#minClusterSize<-as.numeric(args)
#print(minClusterSize)
#in_files <- list('C:/Users/hp/Mu_0.10/network.dat', 'C:/Users/hp/Mu_0.12/network.dat', 'C:/Users/hp/Mu_0.14/network.dat', 'C:/Users/hp/Mu_0.16/network.dat', 'C:/Users/hp/Mu_0.18/network.dat', 'C:/Users/hp/Mu_0.20/network.dat', 'C:/Users/hp/Mu_0.22/network.dat', 'C:/Users/hp/Mu_0.24/network.dat', 'C:/Users/hp/Mu_0.26/network.dat', 'C:/Users/hp/Mu_0.28/network.dat', 'C:/Users/hp/Mu_0.30/network.dat', 'C:/Users/hp/Mu_0.32/network.dat', 'C:/Users/hp/Mu_0.34/network.dat', 'C:/Users/hp/Mu_0.36/network.dat', 'C:/Users/hp/Mu_0.38/network.dat', 'C:/Users/hp/Mu_0.40/network.dat', 'C:/Users/hp/Mu_0.42/network.dat', 'C:/Users/hp/Mu_0.44/network.dat', 'C:/Users/hp/Mu_0.46/network.dat', 'C:/Users/hp/Mu_0.48/network.dat', 'C:/Users/hp/Mu_0.50/network.dat', 'C:/Users/hp/Mu_0.52/network.dat', 'C:/Users/hp/Mu_0.54/network.dat', 'C:/Users/hp/Mu_0.56/network.dat', 'C:/Users/hp/Mu_0.58/network.dat', 'C:/Users/hp/Mu_0.60/network.dat', 'C:/Users/hp/Mu_0.62/network.dat', 'C:/Users/hp/Mu_0.64/network.dat', 'C:/Users/hp/Mu_0.66/network.dat', 'C:/Users/hp/Mu_0.68/network.dat', 'C:/Users/hp/Mu_0.70/network.dat', 'C:/Users/hp/Mu_0.72/network.dat', 'C:/Users/hp/Mu_0.74/network.dat', 'C:/Users/hp/Mu_0.76/network.dat', 'C:/Users/hp/Mu_0.78/network.dat', 'C:/Users/hp/Mu_0.80/network.dat','C:/Users/hp/networks/network.dat')
#in_file<-'C:/Users/hp/Documents/TripleAHC/challenge1/network.dat'
#out_paths <-list('C:/Users/hp/Mu_0.10/', 'C:/Users/hp/Mu_0.12/', 'C:/Users/hp/Mu_0.14/', 'C:/Users/hp/Mu_0.16/', 'C:/Users/hp/Mu_0.18/', 'C:/Users/hp/Mu_0.20/', 'C:/Users/hp/Mu_0.22/', 'C:/Users/hp/Mu_0.24/', 'C:/Users/hp/Mu_0.26/', 'C:/Users/hp/Mu_0.28/', 'C:/Users/hp/Mu_0.30/', 'C:/Users/hp/Mu_0.32/', 'C:/Users/hp/Mu_0.34/', 'C:/Users/hp/Mu_0.36/', 'C:/Users/hp/Mu_0.38/', 'C:/Users/hp/Mu_0.40/', 'C:/Users/hp/Mu_0.42/', 'C:/Users/hp/Mu_0.44/', 'C:/Users/hp/Mu_0.46/', 'C:/Users/hp/Mu_0.48/', 'C:/Users/hp/Mu_0.50/', 'C:/Users/hp/Mu_0.52/', 'C:/Users/hp/Mu_0.54/', 'C:/Users/hp/Mu_0.56/', 'C:/Users/hp/Mu_0.58/', 'C:/Users/hp/Mu_0.60/', 'C:/Users/hp/Mu_0.62/', 'C:/Users/hp/Mu_0.64/', 'C:/Users/hp/Mu_0.66/', 'C:/Users/hp/Mu_0.68/', 'C:/Users/hp/Mu_0.70/', 'C:/Users/hp/Mu_0.72/', 'C:/Users/hp/Mu_0.74/', 'C:/Users/hp/Mu_0.76/', 'C:/Users/hp/Mu_0.78/', 'C:/Users/hp/Mu_0.80/','C:/Users/hp/networks/')


kernel_prep <- function(path,network_name, use_weights) {
  #load(paste0(path, network_name, ".rda"))
  #temp<-read.delim(sprintf ("%s%s"  , path, 'blue_genes.rda'))
  #network_<-c("temp")
  network <- get('temp')
  #print('hi')
  G <- graph_from_data_frame(network, directed=bool_value)
  L <- graph.laplacian(as.undirected(G), normalized = TRUE, weights = use_weights)
  svd_L <- svd(L)
  #dir.create(paste0("kernels/", network_name))
  X<-sprintf ("%s%s"  , path, 'blue_genes_kernel.rda')
  #save(G, L, svd_L, file = X)
  return(list(G = G, L = L, svd_L = svd_L))
}

graph_kernel_diffusion <- function(L, svd_L, alpha){
  kernel <- svd_L$u %*% diag(exp(-alpha * svd_L$d)) %*% t(svd_L$u)
  rownames(kernel) <- rownames(L)
  return(kernel)
}

generate_clusters_1 <- function(data,out_path, minClusterSize,networks,alpha = c(1,15,2)){
  #dir.create(paste0("clusters/", data))
  #alpha = c(1,15,2)
  for (i in 1:3) {
    a=alpha[[i]]
    #dataset <- paste0(out_path, "a", a,"_","kernels",".rda")
    
    #load(dataset)
    dataset<-networks[[i]]
    colnames(dataset) <- rownames(dataset)
    #print(paste0("a",a,"_", data))
    ADJ <- get(paste0("a",a,"_", "network"))
    dissADJ <- 1-ADJ
    hierADJ <- flashClust(as.dist(dissADJ), method="average" )
    hierADJ$height <- sort(hierADJ$height)
    cluster_temp=list()
    for (k in minClusterSize){

      #Cluster detection by a constant height cut of a hierarchical clustering dendrogram
      mods_Static <- cutreeStatic(hierADJ, cutHeight=cutsize, minSize=k) 
      #the ‘Dynamic Tree’ cut, is a top-down algorithm that relies solely on the dendrogram.
      mods_Dynamic <- cutreeDynamic(hierADJ,method="tree",minClusterSize=k)
      #Dynamic Hybrid’ cut, is a bottom-up algorithm that improves the detection of outlying members of each cluster
      mods_Hybrid <- cutreeDynamic(hierADJ,distM= dissADJ, method="hybrid", minClusterSize=k,cutHeight =cutsize, deepSplit=2, pamRespectsDendro = FALSE) 
      
      clusters <- data.frame("gene" = rownames(ADJ), "static" = mods_Static, "dynamic" = mods_Dynamic, "hybrid" = mods_Hybrid)
      
      gene_mods_static <- data.frame("module" = sort(unique(clusters$static[clusters$static != 0])))
      gene_mods_dynamic <- data.frame("module" = sort(unique(clusters$dynamic[clusters$dynamic != 0])))
      gene_mods_hybrid <- data.frame("module" = sort(unique(clusters$hybrid[clusters$hybrid != 0])))
      
      if(nrow(gene_mods_hybrid) >= 1) {
        for (i in 1:nrow(gene_mods_hybrid)) {
          gene_mods_hybrid$genes[i] <- list(clusters$gene[clusters$hybrid == gene_mods_hybrid$module[i]])
          gene_mods_hybrid$size[i] <- length(unlist(gene_mods_hybrid$genes[i]))
        }
      }
      if(nrow(gene_mods_dynamic) >= 1) {
        for (i in 1:nrow(gene_mods_dynamic)) {
          gene_mods_dynamic$genes[i] <- list(clusters$gene[clusters$dynamic == gene_mods_dynamic$module[i]])
          gene_mods_dynamic$size[i] <- length(unlist(gene_mods_dynamic$genes[i]))
        }
      }
      if(nrow(gene_mods_static) >= 1) {
        for (i in 1:nrow(gene_mods_static)) {
          gene_mods_static$genes[i] <- list(clusters$gene[clusters$static == gene_mods_static$module[i]])
          gene_mods_static$size[i] <- length(unlist(gene_mods_static$genes[i]))
        }
      }
      cluster_temp[[len(cluster_temp)+1]]<-list(gene_mods_hybrid, gene_mods_static, gene_mods_dynamic)
      save(gene_mods_hybrid, gene_mods_static, gene_mods_dynamic, file = paste0(out_path, "blue_genes_clusters_a", a, "_k",k,".rda"))
    }
  }
  return (cluster_temp)
}
generate_clusters <- function(out_path, k,dataset,a,method){
    #dir.create(paste0("clusters/", data))
    alphas = c(1,15,2)
  
    #a=alpha[[i]]
    #dataset <- paste0(out_path, "a", a,"_","kernels",".rda")
    
    #load(dataset)
    #dataset<-networks[[i]]
    colnames(dataset) <- rownames(dataset)
    #print(paste0("a",a,"_", data))
    ADJ <- get(paste0("a",a,"_", "network"))
    dissADJ <- 1-ADJ
    hierADJ <- flashClust(as.dist(dissADJ), method="average" )
    hierADJ$height <- sort(hierADJ$height)
    #cluster=list()
    #Cluster detection by a constant height cut of a hierarchical clustering dendrogram
    mods_Static <- cutreeStatic(hierADJ, cutHeight=cutsize, minSize=k) 
    #the ‘Dynamic Tree’ cut, is a top-down algorithm that relies solely on the dendrogram.
    mods_Dynamic <- cutreeDynamic(hierADJ,method="tree",minClusterSize=k)
    #Dynamic Hybrid’ cut, is a bottom-up algorithm that improves the detection of outlying members of each cluster
    mods_Hybrid <- cutreeDynamic(hierADJ,distM= dissADJ, method="hybrid", minClusterSize=k,cutHeight =cutsize, deepSplit=2, pamRespectsDendro = FALSE) 
     
    clusters <- data.frame("gene" = rownames(ADJ), "static" = mods_Static, "dynamic" = mods_Dynamic, "hybrid" = mods_Hybrid)
      
    gene_mods_static <- data.frame("module" = sort(unique(clusters$static[clusters$static != 0])))
    gene_mods_dynamic <- data.frame("module" = sort(unique(clusters$dynamic[clusters$dynamic != 0])))
    gene_mods_hybrid <- data.frame("module" = sort(unique(clusters$hybrid[clusters$hybrid != 0])))
      
    if(nrow(gene_mods_hybrid) >= 1) {
        for (i in 1:nrow(gene_mods_hybrid)) {
          gene_mods_hybrid$genes[i] <- list(clusters$gene[clusters$hybrid == gene_mods_hybrid$module[i]])
          gene_mods_hybrid$size[i] <- length(unlist(gene_mods_hybrid$genes[i]))
        }
      }
    if(nrow(gene_mods_dynamic) >= 1) {
        for (i in 1:nrow(gene_mods_dynamic)) {
          gene_mods_dynamic$genes[i] <- list(clusters$gene[clusters$dynamic == gene_mods_dynamic$module[i]])
          gene_mods_dynamic$size[i] <- length(unlist(gene_mods_dynamic$genes[i]))
        }
      }
    if(nrow(gene_mods_static) >= 1) {
        for (i in 1:nrow(gene_mods_static)) {
          gene_mods_static$genes[i] <- list(clusters$gene[clusters$static == gene_mods_static$module[i]])
          gene_mods_static$size[i] <- length(unlist(gene_mods_static$genes[i]))
        }
      }
    #save(gene_mods_hybrid, gene_mods_static, gene_mods_dynamic, file = paste0(out_path, "blue_genes_clusters_a", a, "_k",k,".rda"))
    if (method=='static'){
      return((gene_mods_static))
    }
    if (method=='dynamic'){
      return(gene_mods_dynamic)
    }
    if (method=='hybrid'){
      return(gene_mods_hybrid)
    }
    #cluster<-list(hybrid=gene_mods_hybrid,static=gene_mods_static, dynamic=gene_mods_dynamic)
    
  
  
}

save_new_submissions_1 <- function(dataset,cluster,out_path, a, k, method, submission_number){
  load(paste0(out_path, "blue_genes_clusters_a", a, "_k",k,".rda"))

  assign("data_to_save", paste0("gene_mods_", method))
  network_data <- get(data_to_save)
  save(network_data, file = paste0(out_path,"submission_final_dataset", ".rda"))
}
save_new_submissions <- function(cluster,out_path ,submission_number){
  #save(cluster, file = paste0(out_path,"submission_final_dataset", ".rda"))
  return(cluster)
  }
  

DREAM_results <- function(data, n, out_path) {
  # Load the data
  #mydata <- load(paste0(out_path, "submission_final_dataset", ".rda"))
  #data <- get(mydata)
  
  # Extract genes
  genes <- lapply(data$genes, as.character)
  #print(length(genes))
  communities=list()
  
  # Check uniqueness of genes
  cat(n, length(unique(unlist(genes))) == length(unlist(genes)), "\n")
  
  # Initialize a list to store genes
  mygenes <- vector("list", length(genes))
  
  # Populate mygenes with gene indices and dummy values
  for (i in 1:length(mygenes)) {
    mygenes[[i]] <- c(i, 1, genes[[i]])
    communities[[length(communities)+1]]<-list()
    community<-list()
    for (j in 3:length(mygenes[[i]])){
        community[[length(community)+1]]<-as.numeric(mygenes[[i]][[j]])
    }
    communities[[length(communities)]]<-community
  }
  
  # Write mygenes to a file
  if (out_path!=''){
  write.table(do.call(rbind, mygenes), file = paste0(out_path, "output_blue_genes.txt"), 
              sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE, append = TRUE)}
  return(communities)
  
  # Uncomment the following section if needed
  # files <- dir(paste0("results/submission_", submission_number, "/submission"), full.names = TRUE)
  # zip(zipfile = paste0("results/submission_", submission_number, "/submission_", submission_number), files = files, zip = "zip")
}

#DREAM_results <- function(submission_number,n,out_path){
#  #results<-c("ppi1","ppi2","signal","coexpr","cancer", "homology")
#  for (l in 1:length(n)){
#    assign("mydata",load(paste0(out_path,"submission_final_dataset", ".rda")))
#    data<-get(mydata)
#    genes<-lapply(data$genes,as.character)
#    cat(n[l],length(unique(unlist(genes)))==length(unlist(genes)),"\n")
#    mygenes<-mapply(c, NA,NA,genes)
#    # made genes to my genes in below line
#    for (i in 1:length(mygenes)) { 
#      print(length(mygenes))
#      mygenes[[i]][1]<-i
#      mygenes[[i]][2]<-1
#    }
#    if (l==1){
#      lapply(mygenes, write, paste0(out_path, "output_blue_genes.txt"), append=TRUE, sep="\t",ncolumns=1000)
#    } 
  
#  }
  #files <- dir(paste0("results/submission_", submission_number, "/submission"), full.names = TRUE)
  #zip(zipfile = paste0("results/submission_", submission_number, "/submission_", submission_number), files = files, zip = "zip")
#}

#temp <- read.delim(in_file, header = FALSE)


networks <- c("temp")
for(n in networks) {
  n.tmp <- get(n)
  #print(n.tmp)
  names(n.tmp) <- c("node1", "node2", "weight") 
  assign(n, n.tmp)
}

#dir.create("r_data")
#X<-sprintf ("%s%s"  , out_path, 'blue_genes_rdata.rda')
#save(temp, file = X ) 
results=kernel_prep(out_path,temp, NA)
G=results$G
L=results$L
svd_L=results$svd_L
#load(paste0(out_path, "blue_genes_kernel", ".rda"))

a1_network <- graph_kernel_diffusion(L, svd_L, alpha=1)
a15_network <- graph_kernel_diffusion(L, svd_L, alpha=1.5)
a2_network <- graph_kernel_diffusion(L, svd_L, alpha=2)

networks<-list(a1_network,a15_network,a2_network)

#save(a1_network, file = sprintf ("%s%s"  , out_path, 'a1_kernels.rda'))
#save(a15_network, file = sprintf ("%s%s"  , out_path, 'a15_kernels.rda'))
#save(a2_network, file = sprintf ("%s%s"  , out_path, 'a2_kernels.rda'))

n<-c("temp")
num<-as.numeric(minlimit)
alpha<-as.numeric(args[[4]])
#num<-25

cluster=generate_clusters(out_path, num,networks[[2]],alpha,args[[5]])

#for (j in n){
#generate_clusters(j,out_path,c(seq(20,25)))
#clusters[[length(clusters)+1]]<-generate_clusters(j,out_path,c(seq(num-1,num+1)),networks)

#}
#save_new_submissions(cluster,out_path, "final")
communities=DREAM_results(cluster,n,out_path)
json_communities<-toJSON(communities,auto_unbox=TRUE)
cat(json_communities)
}))


