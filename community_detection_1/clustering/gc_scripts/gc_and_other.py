"""
Add GeneCentric results to another clustering, with preference for the
GeneCentric clusters. Overlapping clusters in the other clustering are
discarded.

"""
import sys
import argparse
import io_functions as io

def resolve_clusters(gc_clusters, other_clusters):
    overlap_others = set()
    for gc_cluster in gc_clusters:
        for idx, other_cluster in enumerate(other_clusters):
            if set(gc_cluster).intersection(set(other_cluster)):
                overlap_others.add(idx)
    others_without_overlap = [oc for i, oc in enumerate(other_clusters)
                                 if i not in overlap_others]
    return gc_clusters + others_without_overlap

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("gc_clustering")
    parser.add_argument("other_clustering")
    opts = parser.parse_args()

    gc_clusters = io.read_clusters(opts.gc_clustering)
    other_clusters = io.read_clusters(opts.other_clustering)

    clusters = resolve_clusters(gc_clusters, other_clusters)
    io.output_clusters(clusters, '')

if __name__ == '__main__':
    main()

