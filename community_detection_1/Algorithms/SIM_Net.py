import numpy as np
import networkx as nx
from sklearn.cluster import KMeans
from numpy.linalg import inv
from numpy import log, diag, real
from scipy.sparse import csr_matrix, identity
import os

class YakmoKMeans:
    def __init__(self, n_clusters=8, max_iter=300, tol=1e-4, init='k-means++', random_state=None):
        self.n_clusters = n_clusters
        self.max_iter = max_iter
        self.tol = tol
        self.init = init
        self.random_state = random_state
        self.centroids = None
        self.labels_ = None

    def fit(self, X):
        n_samples, n_features = X.shape
        
        # Initialize centroids
        if self.init == 'k-means++':
            self.centroids = self._kmeans_plusplus(X, self.n_clusters)
        elif self.init == 'random':
            random_indices = np.random.choice(n_samples, self.n_clusters, replace=False)
            self.centroids = X[random_indices]
        else:
            raise ValueError("Unsupported initialization method.")

        for iteration in range(self.max_iter):
            old_centroids = self.centroids.copy()
            
            # Assign labels based on closest centroid
            self.labels_ = self._assign_labels(X, self.centroids)
            
            # Recompute centroids
            self.centroids = self._compute_centroids(X, self.labels_, self.n_clusters)
            
            # Orthogonalize centroids
            self.centroids = self._orthogonalize_centroids(self.centroids)

            # Check for convergence
            if np.linalg.norm(self.centroids - old_centroids) < self.tol:
                break

    def _kmeans_plusplus(self, X, n_clusters):
        n_samples = X.shape[0]
        centroids = np.empty((n_clusters, X.shape[1]))
        centroids[0] = X[np.random.randint(n_samples)]

        for k in range(1, n_clusters):
            distances = np.min(euclidean_distances(X, centroids[:k]), axis=1)
            prob = distances / np.sum(distances)
            next_index = np.random.choice(n_samples, p=prob)
            centroids[k] = X[next_index]

        return centroids

    def _assign_labels(self, X, centroids):
        distances = euclidean_distances(X, centroids)
        return np.argmin(distances, axis=1)

    def _compute_centroids(self, X, labels, n_clusters):
        centroids = np.zeros((n_clusters, X.shape[1]))
        for k in range(n_clusters):
            cluster_points = X[labels == k]
            if len(cluster_points) > 0:
                centroids[k] = np.mean(cluster_points, axis=0)
        return centroids

    def _orthogonalize_centroids(self, centroids):
        # Gram-Schmidt process to orthogonalize the centroids
        n_clusters, n_features = centroids.shape
        ortho_centroids = centroids.copy()
        for i in range(1, n_clusters):
            for j in range(i):
                proj = np.dot(centroids[i], ortho_centroids[j]) / np.dot(ortho_centroids[j], ortho_centroids[j])
                ortho_centroids[i] -= proj * ortho_centroids[j]
        return ortho_centroids

    def predict(self, X):
        if self.centroids is None:
            raise ValueError("Model has not been fitted yet.")
        return self._assign_labels(X, self.centroids)

def read_edge_list(file_path):
    """
    Reads an edge list from a file and converts it into a NetworkX graph.
    :param file_path: Path to the edge list file.
    :return: A NetworkX graph.
    """
    G = nx.read_edgelist(file_path, nodetype=int, data=(('weight', float),))
    return G

def graph_to_adjacency_matrix(G):
    """
    Converts a NetworkX graph into an adjacency matrix.
    :param G: A NetworkX graph.
    :return: Adjacency matrix (NumPy array).
    """
    return nx.to_numpy_array(G)
def normalize_adjacency_matrix(A):
    """
    Normalize the adjacency matrix A.
    """
    A = A + csr_matrix(np.eye(A.shape[0]))  # Add self-loops
    d = np.array(A.sum(axis=1)).flatten()
    d_inv = 1. / np.sqrt(d)
    D_inv = csr_matrix(np.diag(d_inv))
    normalized_A = D_inv @ A @ D_inv
    return normalized_A
# Complexity Functions
def gacPathEntropy_bo(IminuszW, z):
    N = IminuszW.shape[0]
    sumall = N
    P = IminuszW
    temp = np.ones(N)
    
    for i in range(10):
        temp = P @ temp
        sumall += np.sum(temp)
    
    clusterComp = sumall / (N * N)
    return clusterComp

def gacPathCondEntropy_bo(P, num_i, num_j, z):
    y_ij = np.zeros((num_i + num_j, 2))  # [y_i, y_j]
    y_ij[:num_i, 0] = 1
    y_ij[num_i:, 1] = 1
    sumall = y_ij.copy()
    
    for i in range(10):
        y_ij = P @ y_ij
        sumall += y_ij
    
    L_ij = (y_ij[:num_i, 0].sum() / (num_i * num_i)) + (y_ij[num_i:, 1].sum() / (num_j * num_j))
    return L_ij

def gacZetaCondEntropy(IminuszW, cluster_i, cluster_j):
    num_i = len(cluster_i)
    num_j = len(cluster_j)

    ijGroupIndex = np.concatenate((cluster_i, cluster_j))
    subMatrix = IminuszW[np.ix_(ijGroupIndex, ijGroupIndex)]

    logZetaSelfSim = np.log(np.real(np.diag(np.linalg.inv(subMatrix))))
    L_ij = (logZetaSelfSim[:num_i].sum() / num_i) + (logZetaSelfSim[num_i:].sum() / num_j)
    
    return L_ij

def gacZetaEntropy(subIminuszW):
    logZetaSelfSim = np.log(np.real(np.diag(np.linalg.inv(subIminuszW))))
    clusterComp = logZetaSelfSim.sum() / subIminuszW.shape[0]
    return clusterComp

# Merging Function
def gacMerging_bo(graphW, initClusters, groupNumber, strDescr, z):
    numSample = graphW.shape[0]
    myInf = 1e10

    graphW = z * graphW

    if strDescr.lower() == 'zeta':
        complexity_fun = gacZetaEntropy
        conditionalComplexity_fun = gacZetaCondEntropy
    elif strDescr.lower() == 'path':
        complexity_fun = gacPathEntropy_bo
        conditionalComplexity_fun = gacPathCondEntropy_bo
    else:
        raise ValueError("GAC: Descriptor type is not supported!")

    numClusters = len(initClusters)
    if numClusters <= groupNumber:
        raise ValueError("GAC: too few initial clusters. Do not need merging!")

    clusterComp = np.zeros(numClusters)
    for i in range(numClusters):
        clusterComp[i] = complexity_fun(graphW[np.ix_(initClusters[i], initClusters[i])], z)

    affinityTab = np.full((numClusters, numClusters), np.inf)
    for j in range(numClusters):
        for i in range(j):
            ijindex = np.concatenate((initClusters[i], initClusters[j]))
            P = graphW[np.ix_(ijindex, ijindex)]
            affinityTab[i, j] = -conditionalComplexity_fun(P, len(initClusters[i]), len(initClusters[j]), z)
    
    affinityTab = affinityTab + clusterComp[:, None] + clusterComp[None, :]

    curGroupNum = numClusters

    while curGroupNum > groupNumber:
        minIndex1, minIndex2 = np.unravel_index(np.argmin(affinityTab[:curGroupNum, :curGroupNum]), (curGroupNum, curGroupNum))

        if minIndex2 < minIndex1:
            minIndex1, minIndex2 = minIndex2, minIndex1

        new_cluster = np.unique(np.concatenate((initClusters[minIndex1], initClusters[minIndex2])))

        if minIndex2 != curGroupNum - 1:
            initClusters[minIndex2] = initClusters[curGroupNum - 1]
            clusterComp[minIndex2] = clusterComp[curGroupNum - 1]
            affinityTab[:minIndex2, minIndex2] = affinityTab[:minIndex2, curGroupNum - 1]
            affinityTab[minIndex2, minIndex2 + 1:curGroupNum] = affinityTab[minIndex2 + 1:curGroupNum, curGroupNum - 1]

        initClusters[minIndex1] = new_cluster.tolist()
        initClusters.pop()

        clusterComp[minIndex1] = complexity_fun(graphW[np.ix_(new_cluster, new_cluster)], z)
        clusterComp[curGroupNum - 1] = myInf
        affinityTab[:curGroupNum, curGroupNum - 1] = myInf
        affinityTab[curGroupNum - 1, :curGroupNum] = myInf

        curGroupNum -= 1

        for groupIndex1 in range(minIndex1):
            ijindex = np.concatenate((initClusters[groupIndex1], new_cluster))
            P = graphW[np.ix_(ijindex, ijindex)]
            affinityTab[groupIndex1, minIndex1] = -conditionalComplexity_fun(P, len(initClusters[groupIndex1]), len(new_cluster), z)
        
        for groupIndex1 in range(minIndex1 + 1, curGroupNum):
            ijindex = np.concatenate((initClusters[groupIndex1], new_cluster))
            P = graphW[np.ix_(ijindex, ijindex)]
            affinityTab[minIndex1, groupIndex1] = -conditionalComplexity_fun(P, len(initClusters[groupIndex1]), len(new_cluster), z)
        
        affinityTab[:minIndex1, minIndex1] = clusterComp[:minIndex1] + clusterComp[minIndex1] + affinityTab[:minIndex1, minIndex1]
        affinityTab[minIndex1, minIndex1 + 1:curGroupNum] = clusterComp[minIndex1 + 1:curGroupNum] + clusterComp[minIndex1] + affinityTab[minIndex1, minIndex1 + 1:curGroupNum]

    clusterLabels = np.ones(numSample, dtype=int)
    for i in range(len(initClusters)):
        clusterLabels[initClusters[i]] = i + 1

    return clusterLabels

# Main Clustering Function with Input and Output Handling
def gdlCluster(graph_input, groupNumber, output_file=None):
    """
    Perform community detection using Graph Diffusion Learning.
    :param graph_input: Either a NetworkX graph or a path to an edge list .txt file.
    :param groupNumber: Desired number of groups.
    :param output_file: Optional output file path to write the communities.
    :return: List of lists representing the communities.
    """
    # Read and process the input graph
    if isinstance(graph_input, str) and os.path.isfile(graph_input):
        G = read_edge_list(graph_input)
    elif isinstance(graph_input, nx.Graph):
        G = graph_input
    else:
        raise ValueError("Invalid graph input. Provide a NetworkX graph or a valid file path to an edge list .txt file.")
    original_nodes = list(G.nodes())
    node_mapping = {node: idx for idx, node in enumerate(original_nodes)}
    reverse_mapping = {idx: node for node, idx in node_mapping.items()}
    
    # Relabel the graph nodes with new sequential integers
    relabeled_G = nx.relabel_nodes(G, node_mapping)
    S = graph_to_adjacency_matrix(relabeled_G)  # Adjacency matrix
    #S=normalize_adjacency_matrix(S)
    F = S  # Assuming F can be approximated as the adjacency matrix for initial clustering

    initial_num_clusters = min(round(len(F) / 10), min(20 * groupNumber, 1000))
    y0 = KMeans(n_clusters=initial_num_clusters, random_state=0).fit_predict(F)

    initialClusters = {i: np.where(y0 == i)[0].tolist() for i in range(np.max(y0) + 1)}
    initialClusters = list(initialClusters.values())  # Convert dictionary to list

    S = S + S.T  # Symmetrize the matrix

    clusteredLabels = gacMerging_bo(S, initialClusters, groupNumber, 'path', 0.1)

    communities = [[] for _ in range(max(clusteredLabels))]
    for idx, label in enumerate(clusteredLabels):
        communities[label - 1].append(idx)

    if output_file:
        with open(output_file, 'w') as f:
            i=1
            for community in communities:
                temp=[str(i)]+community
                f.write(" ".join(map(str, temp)) + "\n")
                i+=1

    return [[int(reverse_mapping[j]) for j in i] for i in communities]




# Example Usage
# Provide a NetworkX graph or an edge list file path
#graph_input = 'path_to_edge_list.txt'  # Replace with the actual file path or NetworkX graph
"""G = nx.Graph()  # Construct networkx graph
datafile=open('Mu_0.10/network.dat', 'r')
for line in datafile:
            g = line.strip().split("\t")
            G.add_edge(g[0], g[1], weight=float(g[2]))
graph_input=G
groupNumber = 30  # Desired number of communities
#output_file = 'output_communities.txt'  # Optional output file path

# Perform clustering
communities = gdlCluster(graph_input, groupNumber)
print("Detected Communities:", communities)"""
