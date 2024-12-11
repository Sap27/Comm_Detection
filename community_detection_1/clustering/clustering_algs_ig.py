"""
Implementations of some graph clustering algorithms

"""
import igraph as ig
import numpy as np
import sklearn.cluster as sc

def threshold_cluster(G, threshold, node_map=[], objects=False):
    """ Calculate threshold clusters for the given similarity scores.

    Args:
        G (ig.Graph)      - the input network
        threshold (float) - the weight above which to remove edges

    Returns:
        clusters (list) - a list of lists of nodes, each sublist represents
                          a cluster
    """
    edges = []
    for edge in G.es:
        edges.append((edge.tuple[0], edge.tuple[1], edge['weight']))
    edges_to_remove = [(n1, n2) for n1, n2, w in edges if w > threshold]
    G.delete_edges(edges_to_remove)

    # hopefully the graph is disconnected now, so filter nodes into bins
    if objects: return G.clusters().subgraphs()
    else:
        clusters = [c for c in G.clusters()]
        if node_map:
            return [[node_map[n] for n in cl] for cl in clusters]
        else: return clusters


def spectral_cluster(G, n_clusters, node_map=[]):
    """ Cluster the given similarity matrix using spectral clustering.

    Assumes the given similarity network is connected.

    Args:
        G (ig.Graph)     - the input network
        n_clusters (int) - number of clusters to look for

    Returns:
        clusters (list) - a list of lists of nodes, each sublist represents
                          a cluster
    """
    # generate a numpy distance matrix from the given graph
    mat = G.get_adjacency(attribute='weight')
    dist_matrix = np.array(mat.data)

    # apply RBF kernel to generate similarity matrix from distance
    # matrix (i.e. lower DSD => higher similarity)
    sim_matrix = np.exp(-(dist_matrix) / (2 *(dist_matrix.std()) ** 2))
    print((sim_matrix.ndim))
    # now do the clustering, scikit-learn implements this
    # return a list of lists representing the clusters
    node_assignments = list(sc.spectral_clustering(sim_matrix,n_clusters=n_clusters))
    clusters = []
    for n in range(n_clusters):
        clusters.append([i for i, m in enumerate(node_assignments) if m == n])
    if node_map:
        return [[node_map[n] for n in cl] for cl in clusters]
    else: return clusters


def hierarchical_clustering(G, threshold=1.0):
    """ Hierarchical clustering using shortest path distances.

    For use as a baseline comparison against our DSD-based methods.
    """
    # TODO: not yet implemented using igraph
    pass

