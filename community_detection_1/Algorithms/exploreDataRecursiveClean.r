
suppressMessages(suppressWarnings({
install.packages("igraph", repos = "https://cloud.r-project.org/")

options(repos = c(CRAN = "https://cran.rstudio.com/"))
install.packages("argparser")

library("igraph")
library("argparser")

countGlob<<-0

mySplitInf <-function(gg1,X,threshold1, threshold2,community_list)
{
#print ("=================================================================================================")
# if you are processing network 3, please change 0.1 to 10 so that you can really prune data 
#print (threshold)

gg2<-delete.edges(gg1, which(E(gg1)$weight <=threshold1))
gg2<-delete.vertices(gg2,V(gg2)[igraph::degree(gg2)<1])


if (threshold2 < 1)  #only needed for network 4
{
gg2<-delete.edges(gg2, which(E(gg2)$weight >=threshold2))
gg2<-delete.vertices(gg2,V(gg2)[igraph::degree(gg2)<1])
}


#print (length(V(gg2)))

clp1 <- walktrap.community(gg2, weights= E(gg2)$weight) 
#print (modularity(clp1))
#print (transitivity(gg2))
#print (length(clp1))

if (modularity(clp1)>0.1)
{
    for (j in 1:length(clp1))
    {
		#print('hi')
		#print(j)
	
        if (length(clp1[[j]])>2 &&    length(clp1[[j]])<=75   )
        {
		community_list[[length(community_list) + 1]] <- clp1[[j]]
		#print(community_list)
	    countGlob<<-countGlob + 1
	    line<-sprintf ("%s %s" , countGlob, "0.5")
	    for (k in 1:length(clp1[[j]]))
		line<-paste ( line,  clp1[[j]][k])	    
	    #print (line)
		if (out_path!=''){
	    write(line,file=X,append=TRUE)
        }}

	else

	{

          if (length(clp1[[j]])>75)
		{
			gg3 <-induced.subgraph(gg2,clp1[[j]])
			#print(community_list)
	        community_list<-Recall(gg3,X,threshold1, threshold2,community_list)

		}

	}
	

   }
 }
membership <- membership(clp1)
#print(membership)
# Create list of lists
#community_list <- lapply(unique(membership), function(x) {
#    which(membership == x)
#})
return (community_list)

}


#### equivalent of Main##################################


#p <- argparser::arg_parser('Dream Challenge script')
#p <- add_argument(p, "--input-file", type="character", nargs=1,
        #help="input filename")
#p<- add_argument(p, "--t1", type="character", nargs=1,
        #help="lower threshold")
#p<- add_argument(p, "--t2", type="character", nargs=1,
        #help="upper threshold")




#args <- parse_args(p)
library(jsonlite)
args <- commandArgs(trailingOnly = TRUE)
t1<-as.numeric(args[[3]])
t2<-as.numeric(args[[4]])
bool_value <- as.logical(args[5])
out_path<-args[[2]]
json_data <- args[[1]]
library(jsonlite)
# Read JSON data from the file
dd <- fromJSON(json_data)


#in_file<- 'C:/Users/hp/community_detection/Social_networks/polblogs.txt'
#out_path<-'C:/Users/hp/community_detection/Social_networks/'
#t1<- 0
#t2<- 1000
#bool_value<-FALSE
#in_file <- args$input_file
#t1 <-args$t1
#t2 <-args$t2

#in_files <- list('C:/Users/hp/Mu_0.10/network.dat', 'C:/Users/hp/Mu_0.12/network.dat', 'C:/Users/hp/Mu_0.14/network.dat', 'C:/Users/hp/Mu_0.16/network.dat', 'C:/Users/hp/Mu_0.18/network.dat', 'C:/Users/hp/Mu_0.20/network.dat', 'C:/Users/hp/Mu_0.22/network.dat', 'C:/Users/hp/Mu_0.24/network.dat', 'C:/Users/hp/Mu_0.26/network.dat', 'C:/Users/hp/Mu_0.28/network.dat', 'C:/Users/hp/Mu_0.30/network.dat', 'C:/Users/hp/Mu_0.32/network.dat', 'C:/Users/hp/Mu_0.34/network.dat', 'C:/Users/hp/Mu_0.36/network.dat', 'C:/Users/hp/Mu_0.38/network.dat', 'C:/Users/hp/Mu_0.40/network.dat', 'C:/Users/hp/Mu_0.42/network.dat', 'C:/Users/hp/Mu_0.44/network.dat', 'C:/Users/hp/Mu_0.46/network.dat', 'C:/Users/hp/Mu_0.48/network.dat', 'C:/Users/hp/Mu_0.50/network.dat', 'C:/Users/hp/Mu_0.52/network.dat', 'C:/Users/hp/Mu_0.54/network.dat', 'C:/Users/hp/Mu_0.56/network.dat', 'C:/Users/hp/Mu_0.58/network.dat', 'C:/Users/hp/Mu_0.60/network.dat', 'C:/Users/hp/Mu_0.62/network.dat', 'C:/Users/hp/Mu_0.64/network.dat', 'C:/Users/hp/Mu_0.66/network.dat', 'C:/Users/hp/Mu_0.68/network.dat', 'C:/Users/hp/Mu_0.70/network.dat', 'C:/Users/hp/Mu_0.72/network.dat', 'C:/Users/hp/Mu_0.74/network.dat', 'C:/Users/hp/Mu_0.76/network.dat', 'C:/Users/hp/Mu_0.78/network.dat', 'C:/Users/hp/Mu_0.80/network.dat')
#in_file<-'C:/Users/hp/Documents/TripleAHC/challenge1/network.dat'
#out_paths <-list('C:/Users/hp/Mu_0.10/', 'C:/Users/hp/Mu_0.12/', 'C:/Users/hp/Mu_0.14/', 'C:/Users/hp/Mu_0.16/', 'C:/Users/hp/Mu_0.18/', 'C:/Users/hp/Mu_0.20/', 'C:/Users/hp/Mu_0.22/', 'C:/Users/hp/Mu_0.24/', 'C:/Users/hp/Mu_0.26/', 'C:/Users/hp/Mu_0.28/', 'C:/Users/hp/Mu_0.30/', 'C:/Users/hp/Mu_0.32/', 'C:/Users/hp/Mu_0.34/', 'C:/Users/hp/Mu_0.36/', 'C:/Users/hp/Mu_0.38/', 'C:/Users/hp/Mu_0.40/', 'C:/Users/hp/Mu_0.42/', 'C:/Users/hp/Mu_0.44/', 'C:/Users/hp/Mu_0.46/', 'C:/Users/hp/Mu_0.48/', 'C:/Users/hp/Mu_0.50/', 'C:/Users/hp/Mu_0.52/', 'C:/Users/hp/Mu_0.54/', 'C:/Users/hp/Mu_0.56/', 'C:/Users/hp/Mu_0.58/', 'C:/Users/hp/Mu_0.60/', 'C:/Users/hp/Mu_0.62/', 'C:/Users/hp/Mu_0.64/', 'C:/Users/hp/Mu_0.66/', 'C:/Users/hp/Mu_0.68/', 'C:/Users/hp/Mu_0.70/', 'C:/Users/hp/Mu_0.72/', 'C:/Users/hp/Mu_0.74/', 'C:/Users/hp/Mu_0.76/', 'C:/Users/hp/Mu_0.78/', 'C:/Users/hp/Mu_0.80/')
#t1<- 0
#t2<- 1000

#print (t1)
#print (t2)
#X<-sprintf ("%s%s"  , "C:/Users/hp/Documents/TripleAHC/challenge1/Results/", 'in_files_new.txt')

#dd <-read.table(in_file)
#gg <- graph.data.frame(dd, directed=FALSE)
#E(gg)$weight<- dd[[3]]
#remove self loops
#gg1<- simplify(gg, remove.multiple = TRUE, remove.loops = TRUE)
#	edge.attr.comb = igraph_opt("edge.attr.comb"))
#mySplitInf(gg1,X,t1, t2)
#folder_info <- file.info("C:/Users/hp/Mu_0.10")
#Sys.chmod("C:/Users/hp/Mu_0.10", "777")
#print(folder_info)

#for (i in 1:1) {
#	print(out_paths[[i]])
#	X<-sprintf ("%s%s"  , out_paths[[i]], 'new_tripleAHC.txt')
#	dd <-read.table(in_files[[i]])
#	gg <- graph.data.frame(dd, directed=FALSE)
##	#remove self loops
#	gg1<- simplify(gg, remove.multiple = TRUE, remove.loops = TRUE,
# 	edge.attr.comb = igraph_opt("edge.attr.comb"))
#	mySplitInf(gg1,X,t1, t2)
#}
#data <- fromJSON(json_data)
#print(length(data))

X<-sprintf ("%s%s"  , out_path, 'output_triple_ahc.txt')
#dd <-read.table(in_file)
#gg <- graph.data.frame(dd, directed=FALSE)
gg <- graph.data.frame(dd, directed=bool_value)

if (length(dd)==3) {
E(gg)$weight<- dd[[3]]
} else {
   E(gg)$weight<-1
}
#remove self loops
gg1<- simplify(gg, remove.multiple = TRUE, remove.loops = TRUE,
edge.attr.comb = igraph_opt("edge.attr.comb"))
community_list <- list()

community_list<-mySplitInf(gg1,X,t1, t2,community_list)
#print(community_list)
#print(clp1)
#membership_data <- split(1:vcount(gg2), membership(clp1))
#json_communities <- toJSON(clp1, auto_unbox = TRUE)
json_communities<-toJSON(community_list,auto_unbox=TRUE)
#print(json_communities)
cat(json_communities)
}))