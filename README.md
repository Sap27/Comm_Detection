This python module contains the codes for implementing 21 community detection algorithms.The module community_detection_1 can be imported in an external script or can be accesed through the CLI.

----Requirements--- Java version 23.01 R version >=4.2 Python Version>=3.10( python Libraries required are mentioned in the setup.py file) WSL Environment

---How to Install--

Go to the directory containing the downloaded package in the terminal.
Run the below command
pip install -e .

---How to Run the Program--

The python module takes in an edgelist file or a networkx graph as input and returns a list of list containing the communities and optionally writes the output into .txt file when the path for the same is provided. An example usage utlising the python CLI is highlighted below.

python -m community_detection_1.main network.dat MLRMCL --output_file output.txt --algorithm_args largest=50 Here, community detection 1 is the name of the python package. The input file is network.dat and the algorithm chosen for this above example is MLRMCL. Further the output path is specified as output.txt and largest=50 is one of the parameter of the algorithm .

List of implemented Algorithms and corresponding Arguements

General Algorithm Args

weighted (can be True or False)

directed (can be True or False)

MLRMCL

largest (integer indicating the maximum allowed size of a community)
filters ("quantile", "pageRank", "double")
inteWeight ("no","yes)
Recursive_Infomap

Recursive=(True, False)
Recursive_Walktrap

t1 (weight threshold - lower limit)
t2 (weight threshold -upper limit)
Walktrap_infomap 1.max_limit=integer indicating the maximum allowed size of a community, 2.method=1 for walktrap, 2 for infomap

walktrap 1.steps=(any integer)

spin_glass

spins (corresponds to no of communities)
louvain

resolution = (0.1-10)
fast_greedy

multiplex-louvain 1.p=(0.1-10)(corresponds to resolution)

TOM_hier 1.min_limit=(minimum size of allowed clusters)

Spectral_clustering

n_clusters=(number of output clusters)
n_components=(dimension of latent representation)
DSD_kernel

num_com=(number of putput communities)
Iterative_SC 1.min_size=(minimum size of output community), 2.max_size=(maximum size of output community), 3. Max_iter=(number of iterations)

exp-laplacian_kernel 1.min_limit=minimum size of output communities, 2.alpha=(can be 1 ,15 or 2) 3. cut_size= (0.999 when community structure not very clear, 0.99 otherwise)

sc_agg 1.groupNumber=number of output communities

Dcut 1.min_limit=minimum size of output community, 2. max_limit=maximum size of output community

hamming_ensemble

SVT 1.n_clusters=number of output communities

label_propagation 1.spins=upper limit for number of communities

shared_neighbor 1.limit=any integer

girvan_newman 1.most_valuable_edge=(tuple (u,v))

Example for importing in an external python script import community_detection_1.main as main res= main.run_community_detection('network.dat','exp-laplacian_kernel')
