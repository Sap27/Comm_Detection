import networkx as nx
from community_detection_1 import dsd_gen
import community_detection_1.clustering.generate_clusters as gc
import community_detection_1.clustering.split_clusters as sc
import tempfile
import os
def node_list_from_graph(G):
    node_list=list(G.nodes)
    return node_list
def adj_matrix_from_graph(G,node_list):
    adj_mat=nx.to_numpy_array(G, node_list)
    return adj_mat
    
def get_matrix_list(G,num_com,is_directed=False):
    node_list=node_list_from_graph(G)
    adj_mat=adj_matrix_from_graph(G,node_list)
    matrix_list =dsd_gen.get_matrix_components(adj_mat, node_list,is_directed=is_directed)
    #return matrix_list 
    for idx, (adj_matrix, nodelist) in enumerate(matrix_list):
        dsd_matrix = dsd_gen.cDSD(adj_matrix)
        #output_filename = "{}.dsd".format(output_prefix)
        #nodelist_filename = "{}.nodelist".format(output_prefix)
        with tempfile.NamedTemporaryFile(delete=False, mode='w', suffix='.dsd') as tmpfile:
            output_DSD = tmpfile.name
            #tmpfile.write(dsd_matrix)
            dsd_gen.write_result_to_file(dsd_matrix, output_DSD)
        with tempfile.NamedTemporaryFile(delete=False, mode='w', suffix='.nodelist') as tmpfile2:
            output_nodelist = tmpfile2.name
            # Write the nodelist to the temporary file
            with open(output_nodelist, 'w') as f:
                for node in nodelist:
                    f.write(f"{node}\n")
        """if nodelist:
                dsd_gen.write_result_to_file(np.array(nodelist), nodelist_filename, '%s')"""
    clusters=gc.generate_clusters_local(output_DSD,1,output_nodelist,'',num_com,directed=is_directed)
    communities=sc.split_clusters_local(output_DSD,clusters,output_nodelist,'')
    os.remove(output_DSD)
    os.remove(output_nodelist)
    return communities


