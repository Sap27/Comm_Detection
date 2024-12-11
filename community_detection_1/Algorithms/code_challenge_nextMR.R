#Usage: Rscript code_challenge.R input_file
#This code require to istall the following package: igraph and WGCNA.
#To install igraph: install.packages("igraph")
#To install WGCNA: 
#source("http://bioconductor.org/biocLite.R") 
#biocLite(c("AnnotationDbi", "impute", "GO.db", "preprocessCore")) 
#install.packages("WGCNA")

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
library(igraph)
options(repos = c(CRAN = "https://cran.rstudio.com/"))
#source("biocLite.R")
#if (!require("BiocManager", quietly = TRUE))
#    install.packages("biocLite")
#BiocManager::install(version = "3.18")
#biocLite(c("AnnotationDbi", "impute", "GO.db", "preprocessCore", "ComplexHeatmap"))
install.packages("WGCNA") 
install.packages("flashClust")
library(WGCNA)
library(flashClust)
library(ComplexHeatmap)

args <- commandArgs(trailingOnly=TRUE)
#srcFile=args[[1]]
json_data <- args[[1]]
library(jsonlite)
# Read JSON data from the file
dd <- fromJSON(json_data)
output_file<-args[[2]]
bool_value <- as.logical(args[[4]])
min_limit<-as.numeric(args[[3]])
#args <- commandArgs(TRUE)
#srcFile ='Documents/network.csv'
#output_file='Documents/output_trial.txt'
#dd <-read.table('Documents/network.dat')
#df <- graph.data.frame(dd, directed=FALSE)
require('igraph')
require('WGCNA')

#print("installation ok")

#net=srcFile
#net = my.files[1]
#dd = read.csv(file = net ,header=F, sep=',')
#print(dim(df))
if (min(dd[,1])==0){
  dd[,c(1,2)] = dd[,c(1,2)]+1
}

#print("read file ok")

g = graph.data.frame(dd,directed=bool_value)
E(g)$weight<- as.numeric(dd[[3]])
#print(length(g))
adj = as.matrix(get.adjacency(g,attr='weight'))

#print("simmetrization start")
# simmetrizza adj se non simmetrica
aux.adj = matrix(0,nrow = dim(adj)[1],ncol = dim(adj)[2])
for (i in 1 : dim(adj)[1]){
  for (j in 1 : dim(adj)[1]){
    aux.adj[j,i] = (adj[j,i]+adj[i,j])/2
  }
}

adj = aux.adj
rm(aux.adj)

if (max(adj)>1){
  adj = (adj - min(adj)) / (max(adj) - min(adj))
}

#print("symmetrization stop")

#print("start TOM")
TOM = TOMsimilarity(adj)
dissTOM = 1-TOM

#print("end TOM")
#print("start WGCNA")
# Clustering
dend.tree = hclust(as.dist(dissTOM), method = "average")
minModuleSize = as.numeric(args[[3]])
minModuleSize=21
unmergedLabels = cutreeDynamic(dendro = dend.tree, distM = dissTOM,
                               deepSplit = 4, pamRespectsDendro = FALSE,
                               minClusterSize = minModuleSize)

unmergedLabels=as.numeric(unmergedLabels[unmergedLabels>0])

#unmergedColors = labels2colors(unmergedLabels)
my.comms = list()
my.subcomms = list()

for (i in  1 : length(sort(unique(unmergedLabels)))){
  
  #print(i)
  adj.aux  = adj[which(unmergedLabels==i),which(unmergedLabels==i)]
  g.aux = graph.adjacency(adj.aux,mode='undirected',weighted=T,diag=F)
  
  # Community detection
  fc = fastgreedy.community(g.aux)
  ind1 = which(unmergedLabels==i)
  
  for (j in 1 : length(unique(membership(fc)))) {
    #print(j)
    
    ind2 = which(membership(fc)==j)
    aux = ind1[ind2]
    #my.subcomms[,j] = aux
    if (length(aux)>2){
      
      my.subcomms[[j]] = aux
      
    
  }
ff=length(my.subcomms)  
pp=length(my.comms) 
si=pp+1
so=pp+ff
#print(pp)
#print(ff)
#print('hi')
 for (l in 1:ff){ 
  my.comms[pp+l] <-my.subcomms[l] 
}
  
}
  

  
}

#print("Stop WGCNA")

func_entr<-function(cor.data) {
  
  
  ############################
  
  t_thr<-seq(0.01,0.8,by=0.01)

  bet<-list()
  ent2<-c()
  
 
  for (i in 1:length(t_thr)) { 
    
    data.thr<-cor.data
    data.thr[data.thr<t_thr[i]]<-0
    data.thr[data.thr!=0]<-1
    
    g.thr<-graph.adjacency(data.thr,mode='undirected',diag=FALSE)
    #betweenness
    g.thr.bet<-betweenness(g.thr,normalized=T)
    
    bet[[i]]<-g.thr.bet
    
    if (max(bet[[i]])==0) { 
      ent2[i]<-0
      
    } else {
      
      ent2[i]<-(-sum(bet[[i]][bet[[i]]>0]*log2(bet[[i]][bet[[i]]>0])))
    }
    
  }
  MA=which.max(ent2)
  
  #prendo la soglia selezionata
  th=t_thr[MA]
  
  #############################
  
  data.thr1<-cor.data
  data.thr1[data.thr1<th]<-0
  data.thr1[data.thr1!=0]<-1
  
  g.thr1<-graph.adjacency(data.thr1,mode='undirected',diag=FALSE)
  ############################
  memb<-membership(fastgreedy.community(g.thr1))
  #return(th)  
  return(g.thr1)
}

ty=unique(my.comms)
my.subcomms1 = list()
my.comms_final=list()
for (j in 1 : length(ty)){
  #print(j)
  if (length(ty[[j]])>100) {
BB=adj[ty[[j]],ty[[j]]]
g.aux1 = func_entr(BB)
  
  # Community detection
  fc1 = fastgreedy.community(g.aux1)
  ind3 = ty[[j]]
  
  for (t in 1 : length(unique(membership(fc1)))) {
    
    ind2 = which(membership(fc1)==t)
    aux = ind3[ind2]
    #aux=ind2
    if (length(aux)>2){
      
      my.subcomms1[[t]] = aux
      
        
      }
      
    }
  ff=length(my.subcomms1)  
  pp=length(my.comms_final) 
  si=pp+1
  so=pp+ff

  for (l in 1:ff){ 
    my.comms_final[pp+l] <-my.subcomms1[l] 
  }
  
      
  }
else {
pp=length(my.comms_final)   
my.comms_final[pp+1]=ty[j]
}  
  }



my.comms_final1=unique(my.comms_final)

nv=vector()
num=0
for (i in 1:length(my.comms_final1)) {
  if (length(my.comms_final1[[i]])>100){
    nv[i]=i
    #print(i)
    num=num+1
  }
}

nv <- na.omit(nv)

my.subcomms11 = list()
my.comms_final2=my.comms_final1
for (j in nv){
  #print(j)
    BB=adj[my.comms_final1[[j]],my.comms_final1[[j]]]
    g.aux1 = func_entr(BB)
    
    # Community detection
    fc1 = fastgreedy.community(g.aux1)
    ind3 = my.comms_final1[[j]]
    
    for (t in 1 : length(unique(membership(fc1)))) {
      
      ind2 = which(membership(fc1)==t)
      aux = ind3[ind2]
      #aux=ind2
      if (length(aux)>2){
        
        my.subcomms11[[t]] = aux
        
        
      }
      
    }
    ff=length(my.subcomms11)  
    pp=length(my.comms_final2) 
    si=pp+1
    so=pp+ff
    
    for (l in 1:ff){ 
      my.comms_final2[pp+l] <-my.subcomms11[l] 
    }
    
    
  }

my.comms_final3=list()
for (t in 1 : length(my.comms_final2)) {
  
  
  if (length(my.comms_final2[[t]])<100){
    
    my.comms_final3[[t]] = my.comms_final2[[t]]
    
    
  }
  
}


my.comms_final4=my.comms_final3[lapply(my.comms_final3,length)>2]

my.comms_final5=unique(my.comms_final4)

my.comms_final5=unique(my.comms_final4)
for (i in  1 : length(my.comms_final5)) {
  my.comms_final5[[i]]=my.comms_final5[[i]] -1
}
#print((my.comms_final5))
communities<- list()
for (i in  1 : length(my.comms_final5)) {
qqqq=list()
qqqq=my.comms_final5[[i]]
communities[[i]] <- qqqq

fff<-c(i,0.5,qqqq)
if (output_file!=''){
write(fff,output_file, ncolumns = 102, append = TRUE, sep="\t")
}}
#print(communities)

json_communities<-toJSON(communities,auto_unbox=TRUE)

cat(json_communities)  # This should be the only output from this script when executed

#cat(json_communities)
}))

