This python module contains the codes for implementing 21 community detection algorithms.The  module community_detection_1 can be imported in an 
external script or can be accesed through the CLI. 

----Requirements---
Java version 23.01
R version >=4.2
Python Version>=3.10( python Libraries required are mentioned in the setup.py file)
WSL Environment

---How to Install--
1. Go to the directory containing the downloaded package in the terminal.
2. Run the below command

pip install -e .

---How to Run the Program--

The python module takes in an edgelist file or a networkx graph as input and returns a list of list containing the communities and optionally writes the output into .txt file when the path for the same is provided. An example usage utlising the python CLI is highlighted below. 

python -m community_detection_1.main network.dat MLRMCL --output_file output.txt --algorithm_args largest=50 
Here, community detection 1 is the name of the python package. The input file is network.dat and the algorithm chosen for this above example is MLRMCL. Further the output path is specified as output.txt and largest=50 is one of the parameter of the algorithm .

List of implemented Algorithms and corresponding  Arguements

General Algorithm Args
1. weighted (can be True or False)
2. directed (can be True or False)

1. MLRMCL
    1. largest (integer indicating the maximum allowed size of a community) 
    2. filters ("quantile", "pageRank", "double")
    3. inteWeight ("no","yes)

2. Recursive_Infomap
    1. Recursive=(True, False)
3. Recursive_Walktrap
    1. t1 (weight threshold - lower limit)
    2. t2 (weight threshold -upper limit)
4. Walktrap_infomap
    1.max_limit=integer indicating the maximum allowed size of a community,
    2.method=1 for walktrap, 2 for infomap
5. walktrap
    1.steps=(any integer)
6. spin_glass
    1. spins (corresponds to no of communities)
7. louvain
    1. resolution = (0.1-10)
8. fast_greedy
9. multiplex-louvain
    1.p=(0.1-10)(corresponds to resolution)
10. TOM_hier
    1.min_limit=(minimum size of allowed clusters)
11. Spectral_clustering
    1. n_clusters=(number of output clusters)
    2. n_components=(dimension of latent representation)
12. DSD_kernel
    1. num_com=(number of putput communities)
13. Iterative_SC
    1.min_size=(minimum size of output community),
    2.max_size=(maximum size of output community), 
    3. Max_iter=(number of iterations)
14. exp-laplacian_kernel
    1.min_limit=minimum size of output communities,
    2.alpha=(can be 1 ,15 or 2)
    3. cut_size= (0.999 when community structure not very clear, 0.99 otherwise)
15. sc_agg
    1.groupNumber=number of output communities
16. Dcut
    1.min_limit=minimum size of output community,
    2. max_limit=maximum size of output community
17. hamming_ensemble
18. SVT
    1.n_clusters=number of output communities
19. label_propagation
    1.spins=upper limit for number of communities
20. shared_neighbor
    1.limit=any integer
21. girvan_newman
    1.most_valuable_edge=(tuple (u,v))

Example for importing in an external python script
import community_detection_1.main as main
res= main.run_community_detection('network.dat','exp-laplacian_kernel')

