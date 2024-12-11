import numpy as np
import scipy.sparse as sp
import pandas as pd
from scipy.sparse.linalg import eigsh
from sklearn.cluster import KMeans
import networkx as nx
from sklearn.impute import SimpleImputer


def score(adj, K):
    # Compute the largest K eigenvalues and corresponding eigenvectors of the adjacency matrix
    values, vectors = eigsh(adj, k=K, which='LM')

    m = adj.shape[0]
    ratio = np.zeros((m, K - 1))

    for j in range(1, K):
        ratio[:, j - 1] = vectors[:, j] / vectors[:, 0]

    # Thresholding the ratio to prevent extreme values
    thresh = np.log(m)
    ratio = np.clip(ratio, -thresh, thresh)
    #ratio_cleaned = ratio[~np.isnan(ratio).any(axis=1)]
    imputer = SimpleImputer(strategy='mean')  # You can use 'median', 'most_frequent', etc.
    ratio_imputed = imputer.fit_transform(ratio)
    kmeans = KMeans(n_clusters=K, n_init=100, random_state=0)
    labels = kmeans.fit_predict(ratio_imputed)

    return labels

def multi_stage_score_local(N, in_file, out_file, Max_size, Min_size, Max_iter, M,G=None,directed=False,weighted=True):
    if G is None:
        """edge_input = np.loadtxt(in_file)
        if edge_input.shape[1] < 3:
            values = np.ones(len(edge_input) * 2)
        else:
            values = np.hstack([edge_input[:, 2], edge_input[:, 2]])

        col1 = np.hstack([edge_input[:, 0] , edge_input[:, 1] ]).astype(int)
        col2 = np.hstack([edge_input[:, 1] , edge_input[:, 0] ]).astype(int)

        net_size = np.max(col1) + 1
        weighted_adj = sp.coo_matrix((values, (col1, col2)), shape=(net_size, net_size))"""
        G = G = nx.DiGraph() if directed else nx.Graph()
        datafile=open(in_file, 'r')
        for line in datafile:
            g = line.strip().split("\t")
            G.add_edge(g[0], g[1], weight=float(g[2]))
    original_nodes = list(G.nodes())
    node_mapping = {node: idx for idx, node in enumerate(original_nodes)}
    reverse_mapping = {idx: node for node, idx in node_mapping.items()}
    
    # Relabel the graph nodes with new sequential integers
    relabeled_G = nx.relabel_nodes(G, node_mapping)
    
    weighted_adj=nx.adjacency_matrix(relabeled_G)
    weighted_adj = weighted_adj.astype(np.float64)
    weighted_adj=sp.coo_matrix(weighted_adj)
    col1=np.array(relabeled_G.nodes)
    
    

    labels = score(weighted_adj, N)
    size_community = np.array([np.sum(labels == i) for i in range(N)])
    block_matrix = [weighted_adj.tocsr()[labels == i, :][:, labels == i] for i in range(N)]
    block_names = [np.where(labels == i)[0] for i in range(N)]

    iter_count = 1
    while np.max(size_community) > Max_size and iter_count <= Max_iter:
        new_block_matrix = []
        new_block_names = []
        new_size_community = []
        
        for i in range(len(block_names)):
            if size_community[i] > Max_size:
                S = min(M, np.ceil(size_community[i] / Max_size).astype(int))
                sub_labels = score(block_matrix[i], S)
                tmp_size_community = np.array([np.sum(sub_labels == j) for j in range(S)])
                tmp_block_matrix = [block_matrix[i].tocsr()[sub_labels == j, :][:, sub_labels == j] for j in range(S)]
                tmp_block_names = [block_names[i][sub_labels == j] for j in range(S)]
                
                new_block_matrix.extend(tmp_block_matrix)
                new_block_names.extend(tmp_block_names)
                new_size_community.extend(tmp_size_community)
            else:
                new_block_matrix.append(block_matrix[i])
                new_block_names.append(block_names[i])
                new_size_community.append(size_community[i])
        
        iter_count += 1
        size_community = np.array(new_size_community)
        block_matrix = new_block_matrix
        block_names = new_block_names
    Ncla = len(size_community)
    if out_file!='':
        with open(out_file, 'a') as file:
            for j in range(Ncla):
                tmp_nodes = np.intersect1d(block_names[j], np.unique(col1)) 
                if len(tmp_nodes) <= Max_size and len(tmp_nodes) >= Min_size:
                    file.write(f"{j} 0.5 {' '.join(map(str, tmp_nodes))}\n")
    communities=[]
    for j in range(Ncla):
        tmp_nodes = np.intersect1d(block_names[j], np.unique(col1))
        if len(tmp_nodes) <= Max_size and len(tmp_nodes) >= Min_size:
            original_labels = [reverse_mapping[idx] for idx in tmp_nodes]
            communities.append(original_labels)
    
    return [[int(j)  for j in i]for i in communities]



# Example usage
#communities=multi_stage_score_local(N=2, in_file='community_detection/network.dat', out_file='community_detection/bigs2_.txt', Max_size=49, Min_size=31, Max_iter=15, M=5)
#print(communities)
#def  multi_stage_score_local(N,path_to_score,in_file,out_file,Max_size,Min_size,Max_iter=15,M=5):
"'"