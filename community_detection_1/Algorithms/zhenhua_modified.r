#This is R code for module detection in sub-challenge one
#The code hasv been tested under R version 3.2.1. R packages required: igraph_1.0.1, plyr_1.8.3
#download the networks of sub-challenge 1 and put them in the directory data/subchallenge1/
install.packages("igraph", repos = "https://cloud.r-project.org/")
#install.packages("plyr_1.8.3", repos = "https://cloud.r-project.org/")
install.packages('dplyr', repos = "https://cloud.r-project.org/")
library(igraph)
library(dplyr)
set.seed(123)
rm(list=ls()) #remove existed objects
#data_2 <- file("stdin")
#print(data_2)

#srcfile=args[[1]]

args <- commandArgs(trailingOnly=TRUE)
json_data <- args[[1]]
library(jsonlite)
#Read JSON data from the file
data_2 <- fromJSON(json_data)
output_file=args[[2]]
max_size=as.numeric(args[[3]])
directed <- as.logical(args[[4]])
method<-as.numeric(args[[5]])

#srcfile='C:/Users/hp/Mu_0.10/network.dat'
#output_file='community_detection/output_zhenhua_1.txt'
#max_size=49
#data_2=read.table(file=srcfile,header=FALSE)
#data_2[,1]=as.character(data_2[,1])
#data_2[,2]=as.character(data_2[,2])
#data_2=as.matrix(data_2)
#g_2=graph.edgelist(data_2[,1:2],directed = FALSE)
#method<-1

g_2=graph_from_data_frame(data_2,directed = directed)

E(g_2)$weight=as.numeric(data_2[,3])

if (method==1){
one=cluster_walktrap(g_2,weights=E(g_2)$weight, steps=4, merges=TRUE, modularity = TRUE, membership=TRUE)
#export modules with size from 3-100 and assign 0.5 as weight score for all moudles
communities<-list()
module_2=data.frame(stringsAsFactors=FALSE)
for (i in 1:length(one)){
  community_members <- as.numeric(names(which(membership(one) == i)))
  if (length(names(which(membership(one)==i))) < 101 & length(names(which(membership(one)==i))) >2 ){
    module_2<-bind_rows(module_2,as.data.frame(t(c(0.5, names(which(membership(one)==i))))))
    communities[[length(communities) + 1]] <- community_members
  }
}}

#Run Infomap
if (method==2){
two=cluster_infomap(g_2,modularity = TRUE)
module_2=data.frame(stringsAsFactors=FALSE)
communities<-list()
for (i in 1:length(two)){
  community_members <- as.numeric(names(which(membership(two) == i)))
  if (length(names(which(membership(two)==i))) < max_size & length(names(which(membership(two)==i))) >2 ){
    module_2<-bind_rows(module_2,as.data.frame(t(c(0.5, as.numeric(names(which(membership(two)==i)))))))
    communities[[length(communities) + 1]] <- community_members
  }
}
#extract modules with size >100
l1_2=list()
l2_2=list()
graph2=list()
infomap_list_2=list()
data_2=as.data.frame(data_2)
#print(length(two))
for (i in 1:length(two)){
  if (length(names(which(membership(two)==i))) > max_size ) {
    l1_2[i] <- list(names(which(membership(two)==i)))
  }
} 
#print(length(l1_2))
l1_2=l1_2[!sapply(l1_2, is.null)]
#print(length(l1_2))
#save modules into l2_2 list
if (length(l1_2)>0){
for (j in 1:length(l1_2)){
  l2_2[j]<- list(data_2[data_2$V1 %in% l1_2[[j]] & data_2$V2 %in% l1_2[[j]],])     
}
#perform Infomap to modules in l2_2 list
for (k in 1:length(l2_2)){
  l2_2[[k]][,1]=as.character(l2_2[[k]][,1])
  l2_2[[k]][,2]=as.character(l2_2[[k]][,2])
  l2_2[[k]]=as.matrix(l2_2[[k]])
  graph2[k]=list(graph.edgelist(l2_2[[k]][,1:2],directed = FALSE))
  E(graph2[[k]])$weight=as.numeric(l2_2[[k]][,3])
  infomap_list_2[k]=list(cluster_infomap(graph2[[k]],modularity = TRUE))
  table(membership(infomap_list_2[[k]]))
  for (n in 1:length(infomap_list_2[[k]])){
    community_members <- as.numeric(names(which(membership(infomap_list_2) == i)))
    if (length(names(which(membership(infomap_list_2[[k]])==n))) < max_size & length(names(which(membership(infomap_list_2[[k]])==n))) >2 ){
      module_2<-bind_rows(module_2,as.data.frame(t(c(0.5,as.numeric(names(which(membership(infomap_list_2[[k]])==n)))))))
      communities[[length(communities) + 1]] <- community_members
    }
  }

}

}}

#print(communities)
json_communities<-toJSON(communities,auto_unbox=TRUE)
cat(json_communities)
if (output_file!=''){
write.table(module_2,file=output_file,sep='\t',row.names = TRUE, col.names = FALSE, na = "")}