import numpy as np
import networkx as nx
from scipy.spatial.distance import pdist, squareform
from sklearn_extra.cluster import KMedoids
import community_pertub as community_louvain
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
# Function to calculate modularity
# Function to calculate modularity
def calculate_modularity(G, resolution):
    modpartit = community_louvain.best_partition(G, resolution=resolution)
    modularit = community_louvain.modularity(modpartit, G)
    return modpartit, modularit

# Function to compute modularities for different resolutions
def compute_modularities(G, resolutions):
    modularity_results = {}
    for res in resolutions:
        modpartit, modularit = calculate_modularity(G, res)
        
        modularity_results[res] = (modpartit, modularit)
    return modularity_results

# Function to create ensemble matrix from modularity results
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

# Function to perform K-medoids clustering on ensemble matrix
def process_ensemble_matrix(Z, k):
    p = Z[:, 1:]
    #print(Z.shape)
    distance_matrix=np.zeros((Z.shape[0],Z.shape[0]))
    for i in range(Z.shape[0]):
        for j in range(Z.shape[0]):
            distance_matrix[i][j]=1-(list(Z[i]-Z[j]).count(0))/Z.shape[1]
    #distance_matrix = squareform(pdist(p, metric='jaccard'))
    print(distance_matrix[0])
    kmedoids = KMedoids(n_clusters=k, metric='precomputed', random_state=1).fit(distance_matrix)
    
    # Check for empty clusters
    labels = kmedoids.labels_
    print(labels)
    unique_labels = np.unique(labels)
    if len(unique_labels) < k:
        k = len(unique_labels)  # Update k to the number of non-empty clusters
    
    idx = kmedoids.labels_
    n2c = np.column_stack((Z[:, 0], idx))
    
    return n2c

# Function to recluster larger modules based on node-to-community assignments
def creating_remaining_network(G, n2c, min_cluster_size=300):
    n2c[:, 1] += 1
    community_list = [[] for _ in range(n2c[:, 1].max())]
    
    for node, com in n2c:
        community_list[com-1].append(node)
    
    large_clusters = []
    for community in community_list:
        if len(community) > min_cluster_size:
            large_clusters.extend(community)
    
    # Filter the graph to include only edges within large clusters
    filtered_network = [(u, v) for u, v in G.edges() if u in large_clusters and v in large_clusters]
    
    return np.array(filtered_network, dtype=int)

# Function to merge node-to-community assignments
def merging_clusters(n2c1, n2c2):
    if n2c1.size == 0:
        return n2c2
    
    if n2c2.size == 0:
        return n2c1
    
    a1 = np.isin(n2c2[:, 0], n2c1[:, 0], invert=True)
    n2c2_filtered = n2c2[a1]
    
    if n2c2_filtered.size == 0:
        max_label = n2c1[:, 1].max()
    else:
        max_label = n2c2_filtered[:, 1].max()
    
    n2c1[:, 1] += max_label + 1
    merged_n2c = np.vstack((n2c2_filtered, n2c1))
    
    return merged_n2c

# Main function
def main(G, resolutions=[0.1 * i for i in range(1, 11)], k=600, min_cluster_size=300):
    # Step 1: Calculate modularity for different resolutions
    modularity_results = compute_modularities(G, resolutions)
    
    # Step 2: Create ensemble matrix and run K-medoids clustering
    Z = creating_ensemble_matrix(modularity_results)
    n2c = process_ensemble_matrix(Z, k)
    #print(n2c)
    # Step 3: Reclustering larger modules
    filtered_network = creating_remaining_network(G, n2c, min_cluster_size)
    
    # Convert filtered_network to a graph for further processing
    G_filtered = nx.Graph()
    G_filtered.add_edges_from(filtered_network)
    
    # Step 4: Merge clusters (if needed)
    # Here we simulate merging with the filtered network itself
    # Adjust this part based on actual implementation needs
    final_n2c = merging_clusters(n2c, n2c)  # Placeholder for actual merging logic
    
    # Convert final_n2c to list of lists of communities
    final_communities = {}
    for node, com in final_n2c:
        if com not in final_communities:
            final_communities[com] = []
        final_communities[com].append(node)
    
    return list(final_communities.values())



# Example usage:
# G = nx.read_edgelist('path_to_edge_list.txt', delimiter='\t', nodetype=int)
G=process_network_file('network.dat')
# data = np.loadtxt('path_to_network_data.txt', dtype=int)
communities = main(G, k=28,min_cluster_size=21)
print(communities)
