import numpy as np
import networkx as nx
from sklearn.cluster import SpectralClustering, AgglomerativeClustering
from scipy.linalg import svd


def denoise_network(W, C, lambda_, beta, max_iter=1000, tol=1e-5):
    n = W.shape[0]
    
    # Initialize S as a stochastic matrix
    S = np.random.rand(n, n)
    S = S / S.sum(axis=1, keepdims=True)
    
    for iteration in range(max_iter):
        S_old = S.copy()
        
        # Update F (using SVD)
        U, _, _ = svd(S - np.eye(n))
        F = U[:, :C]
        
        # Update S
        gradient = -W + lambda_ * (np.dot(F, F.T)) + beta * S
        S = S - 0.01 * gradient
        
        # Ensure S is stochastic
        S = np.maximum(S, 0)
        S = S / S.sum(axis=1, keepdims=True)
        
        # Check for convergence
        if np.linalg.norm(S - S_old) < tol:
            break
    
    # Symmetrize S
    S = (S + S.T) / 2
    return S

def community_detection_pipeline(G, initial_clusters=28, cluster_size_threshold=50, denoise_lambda=1.0, denoise_beta=0.1):
    # Convert NetworkX graph to an adjacency matrix
    original_nodes = list(G.nodes())
    node_mapping = {node: idx for idx, node in enumerate(original_nodes)}
    reverse_mapping = {idx: node for node, idx in node_mapping.items()}
    """largest_cc = max(nx.connected_components(G), key=len)
    G = G.subgraph(largest_cc).copy()"""
    # Relabel the graph nodes with new sequential integers
    relabeled_G = nx.relabel_nodes(G, node_mapping)
    
    W = nx.to_numpy_array(relabeled_G)
    
    # Step 1: Spectral Clustering for Initial Clusters
    spectral = SpectralClustering(n_clusters=initial_clusters, affinity='precomputed', assign_labels='discretize')
    initial_labels = spectral.fit_predict(W)
    
    final_clusters = []
    
    # Step 2: Process each cluster
    for cluster_id in np.unique(initial_labels):
        cluster_indices = np.where(initial_labels == cluster_id)[0]
        if len(cluster_indices) <= cluster_size_threshold:
            # Save as a valid module
            final_clusters.append(cluster_indices.tolist())
        else:
            # Step 3: Extract subgraph and apply denoising
            subgraph = W[np.ix_(cluster_indices, cluster_indices)]
            denoised_subgraph = denoise_network(subgraph, C=initial_clusters, lambda_=denoise_lambda, beta=denoise_beta)
            
            # Step 4: Apply Agglomerative Clustering
            agglomerative = AgglomerativeClustering(n_clusters=None, distance_threshold=0.5, metric='euclidean', linkage='ward')
            agglomerative_labels = agglomerative.fit_predict(denoised_subgraph)
            
            for subcluster_id in np.unique(agglomerative_labels):
                subcluster_indices = cluster_indices[agglomerative_labels == subcluster_id]
                if 3 < len(subcluster_indices) < cluster_size_threshold:
                    final_clusters.append(subcluster_indices.tolist())
    
    return [[int(reverse_mapping[j]) for j in i] for i in final_clusters]

def read_edge_list(file_path):

    G = nx.read_edgelist(file_path, nodetype=int, data=(('weight', float),))
    return G

# Example usage:

# Option 1: Directly use a NetworkX graph
def process_network_file(network_file,weighted=True,directed=False):
    if weighted:
        if directed:
            G=nx.read_edgelist(network_file, data=(('weight', float),), nodetype=int,create_using=nx.DiGraph())
        else:
            G=nx.read_edgelist(network_file, data=(('weight', float),), nodetype=int)
    else:
        if directed:
            G=nx.read_edgelist(network_file, nodetype=int,create_using=nx.DiGraph())
        else:
            G=nx.read_edgelist(network_file, nodetype=int)
    #print(G.nodes)
    return G

"""G=process_network_file('C:/Users/hp/Mu_0.60/network.dat')
clusters = community_detection_pipeline(G)
rc=realcommunities=rc_.groundtruth('C:/Users/hp/Mu_0.60/community.dat')
print(score_benchmarks.nmi_score(clusters,rc))"""

# Option 2: Read from an edge list file
# G = read_edge_list("path_to_edge_list.txt")
# clusters = community_detection_pipeline(G)
# print(clusters)
