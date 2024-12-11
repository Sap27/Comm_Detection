import networkx as nx
import numpy as np
from collections import defaultdict
from infomap import Infomap  # Infomap for community detection

# Load the network from Sample.txt
def load_network(file_path):
    G = nx.read_edgelist(file_path, data=(('weight', float),), nodetype=int)
    return G

# Load clusters from Clusters.txt
def load_clusters(file_path):
    clusters = []
    with open(file_path, 'r') as f:
        for line in f:
            if line.startswith("Cluster"):
                try:
                    cluster = list(map(int, line.split(":")[1].strip().strip('[]').split(',')))
                    clusters.append(cluster)
                except ValueError:
                    continue
    return clusters

# Calculate the inverse log-weighted similarity
def inverse_log_weighted_similarity(G, nodes):
    n = len(nodes)
    similarity_matrix = np.zeros((n, n))
    node_to_index = {node: i for i, node in enumerate(nodes)}

    for i, node1 in enumerate(nodes):
        neighbors1 = set(G.neighbors(node1))
        for j, node2 in enumerate(nodes):
            if i < j:
                neighbors2 = set(G.neighbors(node2))
                common_neighbors = neighbors1 & neighbors2
                if common_neighbors:
                    similarity = sum(1.0 / np.log(sum(G[nbr][x].get('weight', 1.0) for x in G.neighbors(nbr))) for nbr in common_neighbors)
                    similarity_matrix[i][j] = similarity_matrix[j][i] = similarity
    return similarity_matrix

# Sparsify the similarity matrix based on a percentile threshold
def sparsify_matrix(similarity_matrix, percentile):
    threshold = np.percentile(similarity_matrix[similarity_matrix > 0], percentile)
    sparsified_matrix = np.where(similarity_matrix < threshold, 0, similarity_matrix)
    return sparsified_matrix

# Run Infomap clustering on a graph
"""def infomap_clustering(G):
    infomap = Infomap()
    for edge in G.edges(data=True):
        infomap.addLink(edge[0], edge[1], edge[2].get('weight', 1.0))
    infomap.run()
    
    # Convert infomap.modules to a dictionary if necessary
    modules = dict(infomap.modules)
    
    clusters = defaultdict(list)
    for node_id, module_id in modules.items():
        clusters[module_id].append(node_id)
    
    return list(clusters.values())"""
def infomap_clustering(in_file='',weighted=True,directed=False,out_file='',G=None):
    if  G is None:
        G = G = nx.DiGraph() if directed else nx.Graph()
        datafile=open(in_file, 'r')
        for line in datafile:
            g = line.strip().split("\t")
            G.add_edge(g[0], g[1], weight=float(g[2]))
    if weighted and not directed:
        im = Infomap("--two-level")  # Note: No --directed flag
        for source, target,weight in G.edges(data='weight'):
            im.add_link(int(source), int(target), float(weight))
    elif weighted and directed:
        im = Infomap("--two-level --directed")  # Note: No --directed flag
        for source, target,weight in G.edges(data='weight'):
            im.add_link(int(source), int(target), float(weight))  
    elif not weighted and not directed:
        im = Infomap("--two-level")  # Note: No --directed flag
        for source, target in G.edges():
            im.add_link(int(source), int(target))     
    elif not weighted and directed:
        im = Infomap("--two-level --directed")  # Note: No --directed flag
        for source, target in G.edges():
            im.add_link(int(source), int(target)) 
    im.run()
    module_dict = {}
    for node in im.tree:
        if node.is_leaf:
            module_id = node.module_id
            if module_id not in module_dict:
                module_dict[module_id] = []
            module_dict[module_id].append(node.node_id)
    communities = list(module_dict.values())
    if out_file!='':
        with open(out_file, 'a') as file:
            for j in range(Ncla):
                tmp_nodes = np.intersect1d(block_names[j], np.unique(col1)) 
                if len(tmp_nodes) <= Max_size and len(tmp_nodes) >= Min_size:
                    file.write(f"{j} 0.5 {' '.join(map(str, tmp_nodes))}\n")

    return communities


# Function to handle recursive sparsification and clustering
def recursive_sparsify_and_cluster(G, clusters, size_threshold=100, min_size=3, percentile=40,directed=False,weighted=True):
    final_clusters = []
    for cluster in clusters:
        if len(cluster) < min_size:
            continue  # Discard clusters with less than min_size nodes
        elif len(cluster) <= size_threshold:
            final_clusters.append(cluster)
        else:
            subgraph = G.subgraph(cluster)
            similarity_matrix = inverse_log_weighted_similarity(subgraph, list(subgraph.nodes))
            sparsified_matrix = sparsify_matrix(similarity_matrix, percentile)
            sparsified_subgraph = nx.Graph()
            nodes = list(subgraph.nodes)
            for i in range(len(nodes)):
                for j in range(i + 1, len(nodes)):
                    if sparsified_matrix[i, j] > 0:
                        sparsified_subgraph.add_edge(nodes[i], nodes[j], weight=sparsified_matrix[i, j])

            # Perform Infomap clustering on the sparsified subgraph
            new_clusters = infomap_clustering(G=sparsified_subgraph,directed=directed,weighted=weighted)

            if new_clusters == [cluster]:
                final_clusters.append(cluster)
            else:
                final_clusters.extend(recursive_sparsify_and_cluster(G, new_clusters, size_threshold, min_size, percentile))
    return final_clusters

# Save the final modules to Module.txt
def save_clusters(clusters, file_path):
    with open(file_path, 'w') as f:
        for i, cluster in enumerate(clusters, 1):
            f.write(f"Cluster {i}: {cluster}\n")

def main(network_file='',output_file=None,G=None,directed=False,weighted=True,recursive=True):
    
    # Load network and initial clusters
    if network_file!='':
        G = load_network(network_file)
    cluster_assignment=infomap_clustering(G=G,directed=directed,weighted=weighted)
    #print(len(cluster_assignment))
    if recursive:
        final_clusters = recursive_sparsify_and_cluster(G, cluster_assignment,directed=directed,weighted=weighted)
    #print(len(final_clusters))
    # Save the final clusters
    else:
        final_clusters=cluster_assignment
    if output_file:
        save_clusters(final_clusters, output_file)
    return final_clusters

if __name__ == "__main__":
    # File paths
    network_file = "community_detection/network.dat"
    clusters_file = "community_detection/cluster.txt"
    output_file = "community_detection/module.txt"

    main(network_file)







