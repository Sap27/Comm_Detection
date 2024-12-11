"""
Generate "overlap clusters" from a standard network clustering and a
GeneCentric clustering.

"""
import sys
import argparse
import io_functions as io

def get_cluster_nodes(cluster_file):
    cluster_nodes = []
    try:
        fp = open(cluster_file, 'r')
    except IOError:
        sys.exit('Could not open file: {}'.format(cluster_file))

    for line in fp.readlines():
        l = line.rstrip().split()
        for node in l[2:]:
            cluster_nodes.append(node)

    fp.close()
    return cluster_nodes

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("gc_file", help='GeneCentric cluster file')
    parser.add_argument("cluster_file",
                        help="File containing cluster filepaths")
    opts = parser.parse_args()

    cluster_nodes = get_cluster_nodes(opts.cluster_file)
    gc_nodes = get_cluster_nodes(opts.gc_file)
    gc_clusters = io.read_clusters(opts.gc_file)
    difference_nodes = list(set(cluster_nodes) - set(gc_nodes))
    gc_clusters.append(difference_nodes)
    io.output_clusters(gc_clusters, '')


if __name__ == '__main__':
    main()
