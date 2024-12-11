#!/usr/bin/env python
"""
Output clusters in a graph, similar to draw_clusters.py but outputs a
text file rather than an image

"""
import sys
import argparse
import community_detection_1.clustering.io_functions as io
import community_detection_1.clustering.clustering_algs_ig as cl
#g = io.build_ig_graph_from_matrix('C:/Users/hp/Downloads/tusk_dmi_code/data/DSD/6_0.dsd')
#Nodes = io.get_node_list('C:/Users/hp/Downloads/tusk_dmi_code/data/DSD/6_0.dsd')
#clusters = cl.spectral_cluster(g, n_clusters=100, node_map=Nodes)

# default to using spectral clustering
DEFAULT_ALG = 1

# options for clustering algorithm
SPECTRAL = 1
THRESHOLD = 2
HIERARCHICAL = 3
def generate_clusters_local(dsd_file,algorithm,node_list,output_file,parameter,directed=False):
    """parser = argparse.ArgumentParser()
    # parser.add_argument("network_file", help="Original network input file")
    parser.add_argument("dsd_file", help="Distance (i.e. DSD) matrix for network")
    parser.add_argument("-a", "--algorithm", nargs="?", default=DEFAULT_ALG,
                        help="The clustering algorithm to use - 1 for spectral,\
                              2 for threshold clustering, and 3 for simple\
                              shortest-path divisive hierarchical clustering.\
                              Defaults to spectral clustering.")
    parser.add_argument("-d", "--directed", action="store_true",
                        help="Flag specifying if the input represents\
                              a directed graph. Defaults to false.")
    parser.add_argument("-n", "--node_list", nargs="?",
                        help="Optionally specify a list of the nodes in\
                              the DSD file. Default is all the nodes in the\
                              graph.")
    parser.add_argument("-o", "--output_file", nargs="?", default="",
                        help="Optionally specify an output file. Output is to\
                              stdout if no file is specified.")
    parser.add_argument("-p", "--parameter", nargs="?", default='',
                        help="Specify a parameter (i.e. number of clusters,\
                              distance threshold) to be used with clustering\
                              algorithm. If none is provided, a sensible\
                              default is used.")
    opts = parser.parse_args()"""

    G = io.build_ig_graph_from_matrix(dsd_file, directed)

    nodes = io.get_node_list(node_list) if node_list else []

    algorithm = int(algorithm)
    if algorithm == SPECTRAL:
        k_val = int(parameter) if parameter else 100
        clusters = cl.spectral_cluster(G, n_clusters=k_val, node_map=nodes)
    elif algorithm == THRESHOLD:
        filter_weight = float(parameter) if parameter else 5.0
        clusters = cl.threshold_cluster(G, threshold=filter_weight,
                                              node_map=nodes)
    elif algorithm == HIERARCHICAL:
        sys.exit('Hierarchical clustering is not implemented, please choose\
                  another algorithm')
    else:
        sys.exit('Please pick a valid clustering algorithm')
    
    if output_file!='':
        io.output_clusters(clusters, output_file)
    return clusters

"""
def main():
    parser = argparse.ArgumentParser()
    # parser.add_argument("network_file", help="Original network input file")
    parser.add_argument("dsd_file", help="Distance (i.e. DSD) matrix for network")
    parser.add_argument("-a", "--algorithm", nargs="?", default=DEFAULT_ALG,
                        help="The clustering algorithm to use - 1 for spectral,\
                              2 for threshold clustering, and 3 for simple\
                              shortest-path divisive hierarchical clustering.\
                              Defaults to spectral clustering.")
    parser.add_argument("-d", "--directed", action="store_true",
                        help="Flag specifying if the input represents\
                              a directed graph. Defaults to false.")
    parser.add_argument("-n", "--node_list", nargs="?",
                        help="Optionally specify a list of the nodes in\
                              the DSD file. Default is all the nodes in the\
                              graph.")
    parser.add_argument("-o", "--output_file", nargs="?", default="",
                        help="Optionally specify an output file. Output is to\
                              stdout if no file is specified.")
    parser.add_argument("-p", "--parameter", nargs="?", default='',
                        help="Specify a parameter (i.e. number of clusters,\
                              distance threshold) to be used with clustering\
                              algorithm. If none is provided, a sensible\
                              default is used.")
    opts = parser.parse_args()

    G = io.build_ig_graph_from_matrix(opts.dsd_file, opts.directed)

    nodes = io.get_node_list(opts.node_list) if opts.node_list else []

    opts.algorithm = int(opts.algorithm)
    if opts.algorithm == SPECTRAL:
        k_val = int(opts.parameter) if opts.parameter else 100
        clusters = cl.spectral_cluster(G, n_clusters=k_val, node_map=nodes)
    elif opts.algorithm == THRESHOLD:
        filter_weight = float(opts.parameter) if opts.parameter else 5.0
        clusters = cl.threshold_cluster(G, threshold=filter_weight,
                                              node_map=nodes)
    elif opts.algorithm == HIERARCHICAL:
        sys.exit('Hierarchical clustering is not implemented, please choose\
                  another algorithm')
    else:
        sys.exit('Please pick a valid clustering algorithm')

    io.output_clusters(clusters, opts.output_file)


if __name__ == '__main__':
    main()

"""
