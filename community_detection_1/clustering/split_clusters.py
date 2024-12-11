#!/usr/bin/env python
"""
Script for splitting large clusters in a clustering into smaller clusters,
by progressively running spectral clustering with 2 cluster centers (i.e.
finding an approximate min cut)

"""
import sys
import argparse
import community_detection_1.clustering.io_functions as io
import community_detection_1.clustering.clustering_algs_ig as cl

MAX_CL_SIZE = 100
MAX_STEP = 10

def names_to_ids(G, cluster):
    # map from vertex name to vertex ID in G
    return [v.index for v in G.vs if v["name"] in cluster]

def split_clusters_local(dsd_file,clusters,node_list,output_file):
    """parser = argparse.ArgumentParser()
    parser.add_argument("dsd_file", help="Distance (i.e. DSD) matrix for network")
    parser.add_argument("cluster_file", help="Clustering results file")
    parser.add_argument("-n", "--node_list", nargs="?",
                        help="Optionally specify a list of the nodes in\
                              the DSD file. Default is all the nodes in the\
                              graph.")
    opts = parser.parse_args()"""

    node_list = io.get_node_list(node_list)
    #clusters = io.read_clusters(cluster_file)
    G = io.build_ig_graph_from_matrix(dsd_file, False, node_list)

    clusters_to_process, final_clusters = [], []
    for cluster in clusters:
        if len(cluster) > MAX_CL_SIZE:
            clusters_to_process.append(cluster)
        else:
            final_clusters.append(cluster)

    # if all nodes have been clustered, stop looping, otherwise continue to
    # recurse on each large cluster
    step = 1
    while clusters_to_process:
        processing = clusters_to_process
        clusters_to_process = []

        for cluster in processing:
            id_cluster = names_to_ids(G, cluster)
            SG = G.subgraph(cluster)

            cluster_size = len(cluster)
            num_clusters = (int(cluster_size / float(100)) if cluster_size > 200
                                                           else 2)
            clusters = cl.spectral_cluster(SG, num_clusters)
            for cluster in clusters:
                if len(cluster) > MAX_CL_SIZE:
                    clusters_to_process.append([SG.vs[i]['name'] for i in cluster])
                else:
                    final_clusters.append([SG.vs[i]['name'] for i in cluster])
        step += 1
    if output_file!='':
        io.output_clusters(final_clusters, output_file)
    
    return [[int(j) for j in i] for i in final_clusters]

"""def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("dsd_file", help="Distance (i.e. DSD) matrix for network")
    parser.add_argument("cluster_file", help="Clustering results file")
    parser.add_argument("-n", "--node_list", nargs="?",
                        help="Optionally specify a list of the nodes in\
                              the DSD file. Default is all the nodes in the\
                              graph.")
    opts = parser.parse_args()

    node_list = io.get_node_list(opts.node_list)
    clusters = io.read_clusters(opts.cluster_file)
    G = io.build_ig_graph_from_matrix(opts.dsd_file, False, node_list)

    clusters_to_process, final_clusters = [], []
    for cluster in clusters:
        if len(cluster) > MAX_CL_SIZE:
            clusters_to_process.append(cluster)
        else:
            final_clusters.append(cluster)

    # if all nodes have been clustered, stop looping, otherwise continue to
    # recurse on each large cluster
    step = 1
    while clusters_to_process:
        processing = clusters_to_process
        clusters_to_process = []

        for cluster in processing:
            id_cluster = names_to_ids(G, cluster)
            SG = G.subgraph(cluster)

            cluster_size = len(cluster)
            num_clusters = (int(cluster_size / float(100)) if cluster_size > 200
                                                           else 2)
            clusters = cl.spectral_clustering(SG, num_clusters)
            for cluster in clusters:
                if len(cluster) > MAX_CL_SIZE:
                    clusters_to_process.append([SG.vs[i]['name'] for i in cluster])
                else:
                    final_clusters.append([SG.vs[i]['name'] for i in cluster])
        step += 1

    io.output_clusters(final_clusters, '')

if __name__ == '__main__':
    main()

"""