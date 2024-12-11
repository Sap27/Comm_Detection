import pandas as pd
import numpy as np
import os
import igraph
from igraph import *
import re
import time
from rpy2.robjects.packages import importr
import json
import subprocess
import tempfile
from infomap import Infomap
from networkx.algorithms import community
from sklearn.cluster import SpectralClustering
import argparse
import networkx as nx
import community_detection_1.realcommunities_ as rc_
import community_detection_1.score_benchmarks
from  community_detection_1.Algorithms import BigS2
from community_detection_1.Algorithms import sealangbrown_modified 
from community_detection_1.Algorithms import Tianle
from community_detection_1.Algorithms import csbioiitm_2
from  community_detection_1.Algorithms import tusk_dmi_modified
from community_detection_1.Algorithms import tsuromi_ono
from community_detection_1.Algorithms import clustering_module_csbioiitm as csbio_iitm
from community_detection_1.Algorithms import SIM_Net
from  community_detection_1.Algorithms import TeamCS

def graph_json(graph,weighted=True):
    df = nx.to_pandas_edgelist(graph)
    #print(df.head())
    edgelist=[]
    if not weighted:
        df['weight']=[1.0 for i in range(len(df))]
    #print(df.head())
    for i in range(len(df)):
        if df.loc[i]['source']!=df.loc[i]['target']:
            """if weighted:
                edgelist.append([df.loc[i]['source'],df.loc[i]['target'],df.loc[i]['weight']])
            else:
                edgelist.append([float(df.loc[i]['source']),float(df.loc[i]['target']),1.])"""
            edgelist.append([df.loc[i]['source'],df.loc[i]['target'],df.loc[i]['weight']])
    #print(edgelist[0])
    json_list=json.dumps(edgelist)
    return json_list
def convert_to_linux_path(windows_path):
    result = subprocess.run(['wsl', 'wslpath', '-a', windows_path], capture_output=True, text=True)
    if result.returncode == 0:
        linux_path = result.stdout.strip()
        return linux_path
    else:
        print(f"Error converting path: {result.stderr}")
        return None
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
# Triple AHC
#Silence of the noise & Random Walk based
def triple_ahc(t1=0.1,t2=1000,weighted=True,directed=False,network_file='',G=None):
    if network_file!='':
        G=process_network_file(network_file=network_file,directed=directed,weighted=weighted)
    json_list=graph_json(G,weighted=weighted)
    script_dir = os.path.dirname(os.path.abspath(__file__))
    r_script_path_1 = script_dir+'/Algorithms/exploreDataRecursiveClean.r'
    
    with tempfile.NamedTemporaryFile(delete=False, mode='w', suffix='.json') as tmpfile:
        json_file_path = tmpfile.name
        tmpfile.write(json_list)
    if directed:
        bool_val='TRUE'
    else:
        bool_val='False'
    cmd = [
        'Rscript', r_script_path_1,
        json_file_path,  # Passing JSON string instead of input file
        'trial.txt',
        '0.1', '1000',bool_val
        ]
    t1=time.time()
    res=subprocess.run(cmd,capture_output=True, text=True, universal_newlines=True)
    t2=time.time()
    
    if res.returncode != 0:
        print("Error running R script:", res.stderr)
    else:
        # Parse the JSON output from the R script
        #print(res.stdout)
        #print(res.stdout[660:663])

        communities = json.loads(res.stdout)
        #print(communities)
        #print('hi')
    os.remove(json_file_path)
    return [[int(i) for i in j] for j in communities],t2-t1
#triple_ahc(t1=0.1,t2=1000,weighted=True,directed=False,network_file='C:/Users/hp/Mu_0.10/network.dat')
#Zhenhua
#random walk dynamics methods (Walktrap and Infomap) to detect disease modules from biological networks
def zhenhua(max_limit=49,method=1,weighted=True,directed=False,network_file='',G=None):
    # Method 1 is walktrap and Method 2 is infomap
    if network_file!='':
        G=process_network_file(network_file=network_file,directed=directed,weighted=weighted)
    json_list=graph_json(G,weighted=weighted)
    script_dir = os.path.dirname(os.path.abspath(__file__))
    r_script_path_2=script_dir+'/Algorithms/zhenhua_modified.r'
    
    with tempfile.NamedTemporaryFile(delete=False, mode='w', suffix='.json') as tmpfile:
        json_file_path = tmpfile.name
        tmpfile.write(json_list)
    if directed:
        bool='TRUE'
    else:
        bool='False'
    cmd = ['Rscript', r_script_path_2,json_file_path, '',str(max_limit),bool,str(method)]
    t1=time.time()
    res=subprocess.run(cmd,capture_output=True, text=True, universal_newlines=True)
    t2=time.time()
    if res.returncode != 0:
        print("Error running R script:", res.stderr)
    else:
        # Parse the JSON output from the R script
        communities = json.loads(res.stdout)
        #print(len(communities))
        #print(res.stdout)
    os.remove(json_file_path)
    return communities,t2-t1
#NextMR
#Hybrid method
def nextmr(min_limit=21,weighted=True,directed=False,network_file='',G=None):
    if network_file!='':
        G=process_network_file(network_file=network_file,directed=directed,weighted=weighted)
    original_nodes = list(G.nodes())
    node_mapping = {node: idx for idx, node in enumerate(original_nodes)}
    reverse_mapping = {idx: node for node, idx in node_mapping.items()}
    # Relabel the graph nodes with new sequential integers
    relabeled_G = nx.relabel_nodes(G, node_mapping)
    json_list=graph_json(relabeled_G,weighted=weighted)
    script_dir = os.path.dirname(os.path.abspath(__file__))
    r_script_path_3=script_dir+'/Algorithms/code_challenge_nextMR.r'
    
    #print(json_list)
    with tempfile.NamedTemporaryFile(delete=False, mode='w', suffix='.json') as tmpfile:
        json_file_path = tmpfile.name
        tmpfile.write(json_list)
    if directed:
        bool='TRUE'
    else:
        bool='False'
    cmd = ['Rscript', r_script_path_3,json_file_path, '',str(min_limit),bool]
    t1=time.time()
    res=subprocess.run(cmd,capture_output=True, text=True, universal_newlines=True)
    t2=time.time()
    #print(res.stdout)
    if res.returncode != 0:
        print("Error running R script:", res.stderr)
    else:
        pattern = r"\[\[.*?\]\]"
        match = re.search(pattern, res.stdout, re.DOTALL)
        if match:
            list_of_lists_str = match.group(0)
            #communities = json.loads(res.stdout[190:])
            communities=json.loads(list_of_lists_str)
            #print(res.stdout)
            #print(communities)
        os.remove(json_file_path)
    communities=[[int(reverse_mapping[j]) for j in i] for i in communities]
    return communities,t2-t1
def blue_genes(min_limit=21,alpha=15,weighted=True,directed=False,network_file='',cut_size=0.99,G=None):
    # alpha can take values 1,15, 2 where 15 denotes 1.5
    if network_file!='':
        G=process_network_file(network_file=network_file,directed=directed,weighted=weighted)
    json_list=graph_json(G,weighted=weighted)
    script_dir = os.path.dirname(os.path.abspath(__file__))
    r_script_path_4=script_dir+'/Algorithms/Bluegenes_combined.r'
    #json_list=graph_json(G)
    #print(json_list)
    with tempfile.NamedTemporaryFile(delete=False, mode='w', suffix='.json') as tmpfile:
        json_file_path = tmpfile.name
        tmpfile.write(json_list)
    if directed:
        bool='TRUE'
    else:
        bool='False'
    cmd = ['Rscript', r_script_path_4,json_file_path, '',str(min_limit),str(alpha),'static',str(cut_size),bool]
    t1=time.time()
    res=subprocess.run(cmd,capture_output=True, text=True, universal_newlines=True)
    t2=time.time()
    if res.returncode != 0:
        print("Error running R script:", res.stderr)
    else:
        #print(res.stdout[19:])
        #communities = json.loads(res.stdout[19:])
        #print(communities)
        pattern = r"\[\[.*?\]\]"
        match = re.search(pattern, res.stdout, re.DOTALL)
        if match:
            list_of_lists_str = match.group(0)
            #communities = json.loads(res.stdout[190:])
            communities=json.loads(list_of_lists_str)
            print(res.stdout)
            #print(communities)
        os.remove(json_file_path)
    return communities,t2-t1

def causality(network_file='',G=None,weighted=True,directed=False,filters='pageRank',inteWeight='no',largest=50):
    if network_file!='':
        G=process_network_file(network_file=network_file,directed=directed,weighted=weighted)
    json_list=graph_json(G,weighted=weighted)
    
    #json_list=graph_json(G)
    #print(json_list)
    script_dir = os.path.dirname(os.path.abspath(__file__))
    r_script_path_5=script_dir+'/Algorithms/causality_revised.r'
    with tempfile.NamedTemporaryFile(delete=False, mode='w', suffix='.json') as tmpfile:
        json_file_path = tmpfile.name
        tmpfile.write(json_list)
    if weighted:
        bool='TRUE'
    else:
        bool='False'
    cmd = ['Rscript', r_script_path_5,json_file_path, '',filters,inteWeight,bool,str(largest),str(len(G.nodes)),script_dir]
    t1=time.time()
    res=subprocess.run(cmd,capture_output=True, text=True, universal_newlines=True)
    t2=time.time()
    print(res.stdout)
    if res.returncode != 0:
        print("Error running R script:", res.stderr)
    else:
        #print(res.stdout[19:])
        #communities = json.loads(res.stdout[19:])
        #print(communities)
        pattern = r"\[\[.*?\]\]"
        match = re.search(pattern, res.stdout, re.DOTALL)
        if match:
            list_of_lists_str = match.group(0)
            #communities = json.loads(res.stdout[190:])
            communities=json.loads(list_of_lists_str)
            #print(res.stdout)
            #print(communities)
        os.remove(json_file_path)
    return communities,t2-t1

    
# Team Tusk
# Double Spectral Approach to DREAM 11 Subchallenge 3
def tuskdmi(num_com=28,G=None,network_file='',weighted=True,directed=False):
    if network_file!='':
        G=process_network_file(network_file=network_file,directed=directed,weighted=weighted)
    if not weighted:
        for u, v in G.edges():
            G[u][v]['weight'] = 1.
    t1=time.time()
    clusters=tusk_dmi_modified.get_matrix_list(G,num_com,is_directed=directed)
    t2=time.time()
    return clusters,t2-t1

#Tsuromi-ono
# Hybrid
def tsuromi_ono_local(network_file='',threads=4,max_limit=49,min_limit=21,G=None,weighted=True,directed=False):
    if network_file!='':
        G=process_network_file(network_file=network_file,directed=directed,weighted=weighted)
    t1=time.time()
    output=tsuromi_ono.run_dcut(threads=4, max_module_size=49, min_module_size=21,graph=G)
    t2=time.time()
    return output,t2-t1

#csbio-iitm
#Louvain
def csbio_iitm_louvain(network_file='',resolution = 0.1,G=None,weighted=True,directed=False):
    if network_file!='':
        G=process_network_file(network_file=network_file,directed=directed,weighted=weighted)
    if directed:
        G=G.to_undirected()
    t1=time.time()
    communities=csbio_iitm.main(resolution=resolution,G=G)
    t2=time.time()
    return communities,t2-t1
#csbio-iitm
#hierarchical louvain
def csbio_iitm_louvain2(network_file='',G=None,weighted=True,directed=False):
    if network_file!='':
        G=process_network_file(network_file=network_file,directed=directed,weighted=weighted)
    if directed:
        G=G.to_undirected()
    t1=time.time()
    communities=csbioiitm_2.main(G=G)
    t2=time.time()
    return communities,t2-t1
#SealangBrown
#Shared Neighbor Clustering for Disease Module Identification
def sealangbrown(G=None,network_file='',limit=50,weighted=True,directed=False):
    if network_file!='':
        G=process_network_file(network_file=network_file,directed=directed,weighted=weighted)
    t1=time.time()
    com=sealangbrown_modified.SNcluster(graph=G,limit=50,directed=directed)
    t2=time.time()
    return com,t2-t1
# Tianle
def tianle(G=None,weighted=True,network_file='',directed=False,n_clusters=30):
    if network_file!='':
        G=process_network_file(network_file=network_file,directed=directed,weighted=weighted)
    """if not weighted:
        for u, v in G.edges():
            G[u][v]['weight'] = 1."""
    t1=time.time()
    communities=Tianle.sol_dream11('', M=None,G=G,directed=directed, n_clusters=n_clusters)
    t2=time.time()
    return communities,t2-t1
# Big-S2
def big_s2(G=None,weighted=True,network_file='',directed=False,min_size=21,max_size=49, Max_iter=15, M=5):
    if network_file!='':
        G=process_network_file(network_file=network_file,directed=directed,weighted=weighted)
    t1=time.time()
    communities=BigS2.multi_stage_score_local(N=2, in_file='',out_file='', Max_size=max_size, Min_size=min_size, Max_iter=Max_iter, M=M,G=G)
    t2=time.time()
    return [[k for k in j] for j in communities],t2-t1
# Team CS and Aleph(Infomap)
def teamcs_aleph(G=None,weighted=True,directed=False,network_file='',recursive=True):
    if network_file!='':
        G=process_network_file(network_file=network_file,directed=directed,weighted=weighted)
    t1=time.time()
    clust=TeamCS.main(network_file='',output_file=None,G=G,directed=directed,weighted=weighted,recursive=recursive)
    t2=time.time()
    return clust,t2-t1
# SIM_NET
def sim_net(G=None,weighted=True,directed=False,network_file='',groupNumber=28):
    if network_file!='':
        G=process_network_file(network_file=network_file,directed=directed,weighted=weighted)
    t1=time.time()
    communities = SIM_Net.gdlCluster(graph_input=G, groupNumber=groupNumber)
    t2=time.time()
    return communities,t2-t1
# Luminex
def luminex(weighted=True,directed=False,network_file='',G=None,p=1):
    if network_file=='':
        
        if not weighted:
            with tempfile.NamedTemporaryFile(delete=False,suffix='.txt') as temp_file:
                nx.write_edgelist(G, temp_file.name)
                network_file=temp_file.name
                network_file = network_file.replace("\\", "/")
        else:
            with tempfile.NamedTemporaryFile(delete=False,suffix='.txt') as temp_file:
                nx.write_edgelist(G, temp_file.name,data=['weight'])
                
                network_file=temp_file.name
                network_file = network_file.replace("\\", "/")
    network_file=convert_to_linux_path(network_file)
    
    script_dir = os.path.dirname(os.path.abspath(__file__))
    os.chdir(script_dir)
    t1=time.time()
    res=(subprocess.run(f'wsl ./molti-console -p {p} {network_file}',capture_output=True, text=True, universal_newlines=True))
    #res=(subprocess.run(f'wsl {executable_path} -p {p} {network_file}',capture_output=True, text=True, universal_newlines=True))
    t2=time.time()
    print(res.stdout)
    if res.returncode != 0:
        print("Error running:", res.stderr)
    else:
        communities =[[int(j) for j in i] for i in json.loads(res.stdout)]
        print(communities)
    return communities,t2-t1
# Spectral Clustering
def spectral_clustering(weighted=True,directed=False,network_file='',G=None,n_clusters=30,n_components=10):
    if network_file!='':
        G=process_network_file(network_file=network_file,directed=directed,weighted=weighted)
    original_nodes = list(G.nodes())
    node_mapping = {node: idx for idx, node in enumerate(original_nodes)}
    reverse_mapping = {idx: node for node, idx in node_mapping.items()}
    # Relabel the graph nodes with new sequential integers
    relabeled_G = nx.relabel_nodes(G, node_mapping)
    adj_mat=nx.adjacency_matrix(relabeled_G)
    t1=time.time()
    sc=SpectralClustering(n_clusters=n_clusters,affinity='precomputed',n_components=n_components)
    t2=time.time()
    clusters=sc.fit_predict(adj_mat)
    cluster_list=[[] for _ in range(sc.n_clusters)]
    for idx,cluster_label in enumerate(clusters):
        cluster_list[cluster_label].append(idx)
    cluster_list=[[int(reverse_mapping[j]) for j in i] for i in cluster_list]
    return cluster_list,t2-t1
    


#Basic Walktrap
def walktrap(G=None,weighted=True,directed=False,network_file='',steps=10):
    if network_file!='':
        G=process_network_file(network_file=network_file,directed=directed,weighted=weighted)
    """if directed:
        G=G.to_undirected()"""
    original_nodes = list(G.nodes())
    node_mapping = {node: idx for idx, node in enumerate(original_nodes)}
    reverse_mapping = {idx: node for node, idx in node_mapping.items()}
    
    # Relabel the graph nodes with new sequential integers
    G_relabelled = nx.relabel_nodes(G, node_mapping)
    if directed:
        g = igraph.Graph(directed=True)
        g.add_vertices(list(G_relabelled.nodes))
        g.add_edges(list(G_relabelled.edges))
    else:
        g=igraph.Graph.from_networkx(G_relabelled)
    t1=time.time()
    if weighted:
        wtrap = g.community_walktrap(weights=g.es["weight"],steps = steps)
    else:
        wtrap = g.community_walktrap( steps = steps)
    t2=time.time()
    clust = wtrap.as_clustering()
    walktrap_communities=[list(clust[i]) for i in range(len(clust))]
    walktrap_communities=[[int(reverse_mapping[j]) for j in community] for community in walktrap_communities]
    return walktrap_communities,t2-t1

# Spring_glass
# Here, spin refers to the upper limit for number of communities
def label_propogation(G=None,weighted=True,directed=False,network_file='',spins=None):
    if network_file!='':
        G=process_network_file(network_file=network_file,directed=directed,weighted=weighted)
    if directed:
        G=G.to_undirected()
    original_nodes = list(G.nodes())
    node_mapping = {node: idx for idx, node in enumerate(original_nodes)}
    reverse_mapping = {idx: node for node, idx in node_mapping.items()}
    
    # Relabel the graph nodes with new sequential integers
    G_relabelled = nx.relabel_nodes(G, node_mapping)
    
    g=igraph.Graph.from_networkx(G_relabelled)
    t1=time.time()
    if weighted:
        label = g.community_label_propagation(weights=g.es["weight"])
    else:
        label = g.community_label_propagation()
    t2=time.time()
    #clust = label.as_clustering()
    label_communities=[list(label[i]) for i in range(len(label))]
    label_communities=[[int(reverse_mapping[j]) for j in community] for community in label_communities]
    return label_communities,t2-t1

def fast_greedy(G=None,weighted=True,directed=False,network_file=''):
    if network_file!='':
        G=process_network_file(network_file=network_file,directed=directed,weighted=weighted)
    if directed:
        G=G.to_undirected()
    original_nodes = list(G.nodes())
    node_mapping = {node: idx for idx, node in enumerate(original_nodes)}
    reverse_mapping = {idx: node for node, idx in node_mapping.items()}
    
    # Relabel the graph nodes with new sequential integers
    G_relabelled = nx.relabel_nodes(G, node_mapping)
    
    g=igraph.Graph.from_networkx(G_relabelled)
    t1=time.time()
    if weighted:
        label = g.community_fastgreedy(weights=g.es["weight"])
    else:
        label = g.community_fastgreedy()
    t2=time.time()
    label = label.as_clustering()
    label_communities=[list(label[i]) for i in range(len(label))]
    label_communities=[[int(reverse_mapping[j]) for j in community] for community in label_communities]
    return label_communities,t2-t1

def spin_glass(G=None,weighted=True,directed=False,network_file='',spins=25):
    if network_file!='':
        G=process_network_file(network_file=network_file,directed=directed,weighted=weighted)
    if directed:
        G=G.to_undirected()
    original_nodes = list(G.nodes())
    print(G.edges)
    node_mapping = {node: idx for idx, node in enumerate(original_nodes)}
    reverse_mapping = {idx: node for node, idx in node_mapping.items()}
    
    # Relabel the graph nodes with new sequential integers
    G_relabelled = nx.relabel_nodes(G, node_mapping)
    
    g=igraph.Graph.from_networkx(G_relabelled)
    t1=time.time()
    if weighted:
        label = g.community_spinglass(weights=g.es["weight"],spins=spins)
    else:
        label = g.community_spinglass(spins=spins)
    t2=time.time()
    #label = label.as_clustering()
    label_communities=[list(label[i]) for i in range(len(label))]
    label_communities=[[int(reverse_mapping[j]) for j in community] for community in label_communities]
    return label_communities,t2-t1

#def save_clusters(communities, file_path):
    


def run_community_detection(input_data, method='',  output_file='', **kwargs):
    """
    Central function to run community detection algorithms.

    Args:
        input_data: The input data, either a path to a DGEList file or a NetworkX graph.
        method (str): The community detection method to use ('louvain', 'leiden', etc.).
        save_clusters (bool): Whether to save the results to a file.
        output_file (str): Path to the output file where clusters will be saved.
        kwargs: Additional keyword arguments for the specific algorithm.

    Returns:
        clusters: The detected community clusters.
    """

    # Check if the input is a path to a DGEList file or a NetworkX graph
    if isinstance(input_data, str):
        # Assuming DGEList is provided as a file path and needs to be converted to a NetworkX graph
        if os.path.exists(input_data):
            # Code to convert DGEList file to NetworkX graph
            G = process_network_file(input_data)
        else:
            raise FileNotFoundError(f"The file {input_data} does not exist.")
    elif isinstance(input_data, nx.Graph):
        G = input_data
    else:
        raise ValueError("input_data must be a file path to a DGEList or a NetworkX graph.")

    if method == 'triple_ahc':
        clusters,t = triple_ahc(G=G, **kwargs)
    elif method == 'zhenhua':
        clusters,t = zhenhua(G=G, **kwargs)
    elif method=='nextmr':
        clusters,t=nextmr(G=G,**kwargs)
    elif method=='blue_genes':
        clusters,t=blue_genes(G=G,**kwargs)
    elif method=='causality':
        clusters,t=causality(G=G,**kwargs)
    elif method=='tuskdmi':
        clusters,t=tuskdmi(G=G,**kwargs)
    elif method=='tsuromi_ono':
        clusters,t=tsuromi_ono_local(G=G,**kwargs)
    elif method=='sealang_brown':
        clusters,t=sealangbrown(G=G,**kwargs)
    elif method=='tianle':
        clusters,t=tianle(G=G,**kwargs)
    elif method=='csbio_iitm':
        clusters,t=csbio_iitm_louvain(G=G,**kwargs)
    elif method=='csbio_iitm_hier':
        clusters,t=csbio_iitm_louvain2(G=G,**kwargs)
    elif method=='big_s2':
        clusters,t=big_s2(G=G,**kwargs)
    elif method=='sim_net':
        clusters,t=sim_net(G=G,**kwargs)
    elif method=='luminex':
        clusters,t=luminex(G=G,**kwargs)
    elif method=='teamcs_aleph':
        clusters,t=teamcs_aleph(G=G,**kwargs)
    elif method=='spectral_clustering':
        clusters,t=spectral_clustering(G=G,**kwargs)
    elif method=='walktrap':
        clusters,t=walktrap(G=G,**kwargs)
    elif method=='label_propagation':
        clusters,t=label_propogation(G=G,**kwargs)
    elif method=='fast_greedy':
        clusters,t=fast_greedy(G=G,**kwargs)
    elif method=='spin_glass':
        clusters,t=spin_glass(G=G,**kwargs)
    else:
        raise ValueError(f"Unsupported algorithm: {algorithm}")

    # Optionally save the clusters to a file
    if output_file!='':
        #if output_file
        with open(output_file, 'w') as f:
            for i, cluster in enumerate(communities, 1):
                f.write(f"Cluster {i}: {cluster}\n")
    print(clusters)
    return clusters

def main():
    parser = argparse.ArgumentParser(description='Community Detection CLI')
    parser.add_argument('input_data', type=str, help='Path to the input data file')
    parser.add_argument('method', type=str, help='Algorithm to use for community detection')
    parser.add_argument('--output_file', type=str, default='', help='Path to the output file')
    
    # Add a placeholder for algorithm-specific arguments
    parser.add_argument('--algorithm_args', nargs='*', help='Algorithm-specific arguments')
    args = parser.parse_args()
    algorithm_args = {}
    if args.algorithm_args:
        for arg in args.algorithm_args:
            key, value = arg.split('=')
            algorithm_args[key] = value
    run_community_detection(
        input_data=args.input_data,
        output_file=args.output_file,
        method=args.method,
        **algorithm_args
    )

if __name__ == '__main__':
    main()
"""G=process_network_file('C:/Users/hp/community_detection/network.dat')
network_file_1='E:/networksProject/Mu-0.10/R_1/network.dat'"""

"""networks_1=['adjnoun','dolphins','football','internet','karate_new','polblogs','polbooks','power','Email bw 1000 people dig','musae-git']
root_directory='C:/Users/hp/community_detection/Social_networks/'
network_file=root_directory+'hep-weighted'+'.txt'
print(network_file)
G=process_network_file(network_file)"""
"""communities,t=run_community_detection(network_file='',threads=4,max_limit=49,min_limit=21,G=G,weighted=True,directed=False)
print(t)
rc1=rc_.groundtruth('Mu_0.10/community.dat')
nmi=score_benchmarks.nmi_score(rc1,communities)
print(nmi)"""
"""communities,t=luminex(weighted=True,directed=False,network_file='',G=G1,p=1)
realcommunities=rc_.groundtruth('C:/Users/hp/Mu_0.10/community.dat')
nmi=score_benchmarks.nmi_score(realcommunities,communities)
print(nmi)
print(len(communities))
print(t)"""
"""G_dir=G.to_directed()
communities=label_propogation(G=G,weighted=True,directed=False,network_file='')"""
#print(communities)
"""communities=run_community_detection(input_data='C:/Users/hp/Mu_0.10/network.dat',method='blue_genes')
print(communities)"""
"""# label_propogation
clusters_label_propagation = g.community_label_propagation()
#k-clique
#def k_clique_
communities = list(k_clique_communities(G, 3))
    
# Convert sets to lists for easier processing and return
communities= [list(community) for community in communities]
communities=[[int(j) for j in community] for community in communities]"""


"""networks_1=['adjnoun','dolphins','football','internet','karate','polblogs','polbooks','power','Email bw 1000 people dig','musae-git']
networks_2=['cond-mat-weighted','lesmis-weighted','hep-weighted','netscience-weighted']
networks_4=['celegansneural']
root_directory='C:/Users/hp/community_detection/Social_networks/'


#triple_ahc(t1=0.1,t2=1000,weighted=True,directed=False,network_file='',G=None)
for i in range(len(networks_1)):
    print(i)
    network_file_curr=f'{root_directory}{networks_1[i]}.txt'
    G=process_network_file(network_file_curr,weighted=False,directed=False)
    #t1=time.time()
    pc=triple_ahc(t1=0.1,t2=1000,weighted=False,directed=False,network_file='',G=G)
    #t2=time.time()
    #
    #times[0][i]=t2-t1
    #ccbar_scores[0][i]=len(pc)/len(ground_truths_1[i])


#Spectral Clustering
"""