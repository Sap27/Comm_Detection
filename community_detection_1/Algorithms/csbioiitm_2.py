import numpy as np
import networkx as nx
from scipy.spatial.distance import pdist
from scipy.cluster.hierarchy import linkage, fcluster
import community_detection_1.community_pertub as community
from collections import defaultdict

# Construct the graph from the input file
"""def construct_graph(input_file):
    G = nx.Graph()
    with open(input_file, 'r') as datafile:
        for line in datafile:
            g = line.strip().split("\t")
            G.add_edge(int(g[0]), int(g[1]), weight=float(g[2]))
    return G"""

# Calculate modularity and get community partitions
def calculate_modularity(G, resolution):
    modpartit = community.best_partition(G, resolution=resolution)
    modularit = community.modularity(modpartit, G)
    return modpartit, modularit

# Create ensemble matrix from modularity results
def creating_ensemble_matrix(modularity_results):
    resolutions = sorted(modularity_results.keys())
    data = []
    
    for res in resolutions:
        partition = list(modularity_results[res][0].values())
        if isinstance(partition, list):  # Ensure it's a 2D array
            partition = np.array(partition).reshape(-1, 1)
        data.append(partition)
    
    data = np.hstack(data)  # Stack all the partitions horizontally
    return np.hstack((np.arange(len(data)).reshape(-1, 1), data))


# Perform hierarchical clustering based on ensemble matrix
def ensemble_hierarchical_clustering(Z, G, strength=2):
    p = Z[:, 1:]  # Use the matrix without the first column
    dist = pdist(p, 'hamming')  # Compute Hamming distance

    # Apply thresholding based on the strength parameter
    d1 = dist.copy()
    d1[d1 < (strength / 10)] = 0

    # Perform hierarchical clustering
    l = linkage(d1, method='single')
    c = fcluster(l, t=np.finfo(float).eps, criterion='distance')  # Using 'distance' as the criterion

    # Prepare the clustered nodes
    clustered_nodes = []
    for k in range(1, max(c) + 1):
        d = Z[c == k, 0].tolist()  # Get the node IDs for each cluster
        clustered_nodes.append(d)

    # Process the small remaining network
    clustered_nodes = creating_small_remaining_network(clustered_nodes, G)

    # Optionally, integrate these nodes back into clusters or handle them as needed
    """for i in small_clusters:
        clustered_nodes.append(i)  # Example integration"""

    return clustered_nodes

def creating_small_remaining_network(clustered_nodes, G):
    large_clusters = [cluster for cluster in clustered_nodes if len(cluster) >= 3]
    small_remaining_nodes = [node for cluster in clustered_nodes if len(cluster) < 3 for node in cluster]
    
    # Generate a new subnetwork from the small remaining nodes
    remaining_nodes_set = list(set(small_remaining_nodes))
    H = G.subgraph(remaining_nodes_set).copy()

    if len(H.edges)>0:
        partitions,_=calculate_modularity(H, 0.1)
        #print(partitions)1
        community_dict = defaultdict(list)
        for node, community in partitions.items():
            community_dict[community].append(node)

        # Convert dictionary values to list of lists
        small_clusters = [nodes for nodes in community_dict.values()]

        for i in small_clusters:
            large_clusters.append(i)
    return large_clusters


# Main function that ties everything together
def main(network_file='',G=None):
    #if network_file!
    #input_file = 'network.dat'
    original_nodes = list(G.nodes())
    node_mapping = {node: idx for idx, node in enumerate(original_nodes)}

    reverse_mapping = {idx: node for node, idx in node_mapping.items()}
    
    # Relabel the graph nodes with new sequential integers
    relabeled_G = nx.relabel_nodes(G, node_mapping)
    resolutions = [0.1 * i for i in range(1, 11)]

    # Step 1: Construct the graph
    #G = construct_graph(input_file)
    
    # Step 2: Calculate modularity for each resolution and store results in a dictionary
    modularity_results = {}
    for resolution in resolutions:
        modpartit, modularit = calculate_modularity(relabeled_G, resolution)
        modularity_results[resolution] = (modpartit, modularit)

    # Step 3: Create ensemble matrix based on different resolutions
    ensemble_matrix = creating_ensemble_matrix(modularity_results)
    
    # Step 4: Perform hierarchical clustering on the ensemble matrix
    final_clusters = ensemble_hierarchical_clustering(ensemble_matrix,relabeled_G)
    print(final_clusters)
    return [[int(reverse_mapping[j]) for j in i] for i in final_clusters]

"""if __name__ == "__main__":
    final_clusters = main()
    print("Final Clusters:", len(final_clusters))
"""
