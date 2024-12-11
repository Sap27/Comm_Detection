import numpy as np
from scipy.stats import spearmanr, scoreatpercentile
from scipy.sparse.linalg import svds
from sklearn.utils.extmath import randomized_svd
from sklearn.cluster import AgglomerativeClustering
from sklearn.neighbors import kneighbors_graph
from sklearn.decomposition import PCA

import time
import _pickle as pickle
import networkx as nx

### SVT feature extraction
def svt_feature(mat, M = None, include_diagonal = False, svd_k = 3, svd_maxiter = 100, svt_delta = 1.5, e=0.0001, svt_maxiter=100, save_feature=False, feature_filename=''):    
    if M is None or include_diagonal:
        idx1, idx2 = mat.nonzero()
        M = mat[idx1, idx2]

    
    else:
        idx1 = np.array(M[:,0], dtype='int')
        idx2 = np.array(M[:,1], dtype='int')
        M = M[:,2]
    
    tic = time.time()
    # Efficient svd implementation from scipy
    U,s,V = svds(mat, k = svd_k, maxiter = svd_maxiter)
    #U, s, V = randomized_svd(mat, n_components=svd_k, n_iter=5, random_state=None)
    i = 0
    err = 100    # Initial error set a large number
    # First svd approximation
    Y = U.dot(np.diag(s)).dot(V)
    while err > e and i <= svt_maxiter:
        Y[idx1, idx2] += svt_delta * (M - Y[idx1, idx2])
        U,s,V = svds(Y, k = svd_k, maxiter=svd_maxiter)
        Y = U.dot(np.diag(s)).dot(V)
        i += 1
        err = np.sum((M - Y[idx1, idx2])**2) / np.sum(M**2)
        #print ('iteration:', i, ' relative error:', err)
        if i % 10 == 0:
            c = np.corrcoef(U.dot(np.diag(np.sqrt(s))))
            #print ('correlation:', np.corrcoef(c[idx1, idx2], M)[0,1])
    #print ('SVT time:', time.time()-tic)
    #print ('Top ', svd_k, ' singlar values:', np.sqrt(s))
    X = U.dot(np.diag(np.sqrt(s)))
    #print(X.shape)
    if save_feature:
        pickle.dump(X, open(feature_filename, 'wb'), -1) 
    
    """pca = PCA(n_components=30)
    X = pca.fit_transform(X)"""
     
    return X

### Second-level clustering
def sub_cluster(X, nodes, module_ass, module_size = 40, linkage='ward', constraint = False, n_neighbors=100):
    if constraint:
        connectivity = kneighbors_graph(X[nodes], n_neighbors=n_neighbors, include_self=False)
        ward_res = AgglomerativeClustering(n_clusters = nodes.shape[0]//module_size, connectivity=connectivity, linkage=linkage).fit(X[nodes])
    else:
        ward_res = AgglomerativeClustering(n_clusters = nodes.shape[0]//module_size, linkage=linkage).fit(X[nodes])
    label = ward_res.labels_
    c = np.bincount(label)
    #print c,np.sum(c), np.where(np.logical_and(c>2, c<101))[0].shape
    num = max(module_ass) + 1
    for i in np.where(np.logical_and(c>2, c<101))[0]:
        module_ass[nodes[label == i]] = num
        num += 1
    for i in np.where(c>100)[0]:
        module_ass = sub_cluster(X, nodes[label == i], module_ass)
    return module_ass

### Module discovery based on learned features
def module_disc(X, output_filename, n_neighbors = 100, n_clusters=1000, linkage='ward', module_size = 40, constraint2 = False, n_neighbors2=100):
    # Initial clustering with connectivity constraint
    #print("Compute structured hierarchical clustering...")
    connectivity = kneighbors_graph(X, n_neighbors=n_neighbors, include_self=False)    
    ward = AgglomerativeClustering(n_clusters=n_clusters, connectivity=connectivity, linkage=linkage).fit(X)
    label = ward.labels_
    #print("Number of points: %i" % label.size)
    
    c = np.bincount(label) 
    nodes = np.array(range(label.shape[0]))
    module_ass = np.zeros(label.shape)
    num = 1
    # Identify modules of size 3 to 100 from initial clustering result
    for i in np.where(np.logical_and(c>2, c<101))[0]:
        module_ass[label == i] = num
        num += 1
    #print (output_filename)
    #print ('Initial clustering: ', np.where(np.logical_and(c>2, c<101))[0].shape[0], ' clusters')
    #print  (np.where(module_ass!=0)[0].shape[0], ' out of ', module_ass.shape[0], ' nodes included')
    #print ('Percentage: ', 1.0 * np.where(module_ass!=0)[0].shape[0] / module_ass.shape[0])
    #print ('Average module size:', 1.0 * np.where(module_ass!=0)[0].shape[0] / np.where(np.logical_and(c>2, c<101))[0].shape[0])
    for i in np.where(c>100)[0]:
        module_ass = sub_cluster(X, nodes[label == i], module_ass, module_size = module_size, linkage=linkage, constraint = constraint2, n_neighbors=n_neighbors2)
    
    c = np.bincount(np.array(module_ass, dtype='int'))
    #print (output_filename)
    #print ('Final result:')
    num_clus = np.where(np.logical_and(c>2, c<101))[0].shape[0] - np.logical_and(c[0]>2, c[0]<101)
    print ('Total clusters:', num_clus)
    #print  (np.where(module_ass!=0)[0].shape[0], ' out of ', module_ass.shape[0], ' nodes included')
    # print 'Excluded:', np.where(module_ass == 0)[0].shape[0]
    #print ('Percentage:', 1.0 * np.where(module_ass != 0)[0].shape[0] / module_ass.shape[0])
    #print ('Average module size', 1.0 * np.where(module_ass!=0)[0].shape[0] / num_clus)
    if output_filename!='':
        with open(output_filename, 'w') as f:
            for i in range(1, int(max(module_ass))+1):
                f.write(str(i) + '\t0.5\t')
                for j in np.where(module_ass==i)[0]:
                    f.write(str(j)+'\t')
                f.write('\n')
    
    communities=[]
    for i in range(1, int(max(module_ass))+1):
        communities.append([])
        for j in np.where(module_ass==i)[0]:
            communities[-1].append(j)
    return communities



### The main function
def sol_dream11(output_filename, input_filename='', M=None, include_diagonal = False, G=None, directed = False, svd_k = 50, svd_maxiter = 100, svt_delta = 1.5, e=0.0001, svt_maxiter=100, save_feature=False, feature_filename='', n_neighbors = 100, n_clusters=1000, linkage='ward', module_size = 40, constraint2 = False, n_neighbors2=100):
    if M is None and G is None:
        M = np.loadtxt(open(input_filename, 'r'), delimiter='\t')
    if G is None:
        """idx1 = np.array(M[:,0], dtype='int')
        idx2 = np.array(M[:,1], dtype='int')
        num = np.max(np.hstack((idx1, idx2)))+(1)
        mat = np.zeros((num,num))
        mat[idx1, idx2] = M[:,2]
        if not directed:                          
            mat[idx2, idx1] = M[:,2]
        mat[range(num), range(num)] = 1       """ 
        G = G = nx.DiGraph() if directed else nx.Graph()
        #datafile=open(in_file, 'r')
        for g in M:
            G.add_edge(g[0], g[1], weight=float(g[2]))
        original_nodes = list(G.nodes())
        node_mapping = {node: idx for idx, node in enumerate(original_nodes)}
        for edge in M:
            edge[0]=node_mapping[edge[0]]
            edge[1]=node_mapping[edge[1]]

        reverse_mapping = {idx: node for node, idx in node_mapping.items()}
    
        # Relabel the graph nodes with new sequential integers
        relabeled_G = nx.relabel_nodes(G, node_mapping)
        A=nx.adjacency_matrix(relabeled_G)
        mat=A.toarray()
    else:
        #A=nx.adjacency_matrix(G)
        original_nodes = list(G.nodes())
        node_mapping = {node: idx for idx, node in enumerate(original_nodes)}
        reverse_mapping = {idx: node for node, idx in node_mapping.items()}
        #print(node_mapping)
        # Relabel the graph nodes with new sequential integers
        relabeled_G = nx.relabel_nodes(G, node_mapping)
        mat=nx.adjacency_matrix(relabeled_G)
        """mat=B.toarray()         
        mat = mat.astype(float)  """
        #mat=B      
    X = svt_feature(mat, M = M, include_diagonal = include_diagonal, svd_k = svd_k, svd_maxiter = svd_maxiter, svt_delta = svt_delta, e=e, svt_maxiter=svt_maxiter, save_feature=save_feature, feature_filename=feature_filename)
    communities=module_disc(X, output_filename, n_neighbors = n_neighbors, n_clusters=n_clusters, linkage=linkage, module_size = module_size, constraint2 = constraint2, n_neighbors2=n_neighbors2)                          
    return [[int(reverse_mapping[j]) for j in i] for i in communities]

### Default parameter settings used for final submissions. 
### subchallenge 1
"""print ('PPI1:')
sol_dream11('submitted/1_ppi_480.txt', input_filename='subchallenge1/1_ppi_anonym_v2.txt', n_clusters=480, save_feature=True, feature_filename='submitted/PPI1_feature.pkl')

print ('PPI2:')
sol_dream11('submitted/2_ppi_800.txt', input_filename='subchallenge1/2_ppi_anonym_v2.txt', n_clusters=800, save_feature=True, feature_filename='submitted/PPI2_feature.pkl')
"""
"""print ('SIGNAL:')
signal = np.loadtxt(open("C:/Users/hp/Downloads/sol_dream11 (1)/sol_dream11/subchallenge1/3_signal_anonym_directed_v3.txt", "r"), delimiter='\t')
a = scoreatpercentile(signal[:,2], 75)
b = scoreatpercentile(signal[:,2], 25)
idx = signal[:,2] > a + 3*(a-b)
signal[idx,2] = 3 * (signal[idx, 2] - np.min(signal[idx,2])) / (np.max(signal[idx,2]) - np.min(signal[idx,2])) + np.min(signal[idx,2])
signal[:,2] = 4*(signal[:,2]-np.min(signal[:,2])) / (np.max(signal[:,2]) - np.min(signal[:,2])) + np.min(signal[:,2])
sol_dream11('C:/Users/hp/Downloads/sol_dream11 (1)/sol_dream11/submitted/3_signal_260.txt', M=signal, directed = True, svd_k=30, n_clusters=260, save_feature=True, feature_filename='C:/Users/hp/Downloads/sol_dream11 (1)/sol_dream11/submitted/SIGNAL_feature.pkl')

print ('COEXPR:')
sol_dream11('submitted/4_coexpr_340.txt', input_filename='C:/Users/hp/Downloads/sol_dream11 (1)/sol_dream11/subchallenge1/4_coexpr_anonym_v2.txt', n_clusters=340, save_feature=True, feature_filename='submitted/COEXPR_feature.pkl')

print ('CANCER:')
sol_dream11('submitted/5_cancer_250.txt', input_filename='subchallenge1/5_cancer_anonym_v2.txt', n_clusters=250, save_feature=True, feature_filename='submitted/CANCER_feature.pkl')
"""
#print ('HOMOLOGY:')
#homology = np.loadtxt(open("C:/Users/hp/Downloads/Tianle(2nd)/sol_dream11/subchallenge1/cora_edgelist_test.txt", "r"), delimiter='\t')
"""a = scoreatpercentile(homology[:,2], 75)
b = scoreatpercentile(homology[:,2], 25)
idx = homology[:,2] > a + 1.5*(a-b)
homology[idx,2] = 3 * (homology[idx, 2] - np.min(homology[idx,2])) / (np.max(homology[idx,2]) - np.min(homology[idx,2])) + np.min(homology[idx,2])
homology[:,2] = 3*(homology[:,2]-np.min(homology[:,2])) / (np.max(homology[:,2]) - np.min(homology[:,2])) + np.min(homology[:,2])
"""
#sol_dream11('C:/Users/hp/Downloads/Tianle(2nd)/sol_dream11/submitted/cora_edgelist.txt', M=homology, n_clusters=7, save_feature=True, feature_filename='C:/Users/hp/Downloads/Tianle(2nd)/sol_dream11/submitted//HOMOLOGY_feature.pkl')


### subchallenge 2

## Load the data
"""ppi1 = np.loadtxt(open("subchallenge2/1_ppi_anonym_aligned_v2.txt", "r"), delimiter='\t')
ppi2 = np.loadtxt(open("subchallenge2/2_ppi_anonym_aligned_v2.txt", "r"), delimiter='\t')
signal = np.loadtxt(open("subchallenge2/3_signal_anonym_aligned_directed_v3.txt", "r"), delimiter='\t')
coexpr = np.loadtxt(open("subchallenge2/4_coexpr_anonym_aligned_v2.txt", "r"), delimiter='\t')
cancer = np.loadtxt(open("subchallenge2/5_cancer_anonym_aligned_v2.txt", "r"), delimiter='\t')
homology = np.loadtxt(open("subchallenge2/6_homology_anonym_aligned_v2.txt", "r"), delimiter='\t')

idx = np.hstack((ppi1[:,0],ppi2[:,0],signal[:,0],coexpr[:,0], cancer[:,0], homology[:,0], ppi1[:,1],ppi2[:,1],signal[:,1],coexpr[:,1], cancer[:,1], homology[:,1]))
num = int(max(idx)) + 1
mat = np.zeros((num,num))

M = ppi1
mat[np.array(M[:,0], dtype='uint'), np.array(M[:,1],dtype='uint')] = M[:,2]
mat[np.array(M[:,1], dtype='uint'), np.array(M[:,0],dtype='uint')] = M[:,2]
M = ppi2
mat[np.array(M[:,0], dtype='uint'), np.array(M[:,1],dtype='uint')] += M[:,2]
mat[np.array(M[:,1], dtype='uint'), np.array(M[:,0],dtype='uint')] += M[:,2]
#M = signal
#mat[np.array(M[:,0], dtype='uint'), np.array(M[:,1],dtype='uint')] += M[:,2]
M = coexpr
mat[np.array(M[:,0], dtype='uint'), np.array(M[:,1],dtype='uint')] += M[:,2]
mat[np.array(M[:,1], dtype='uint'), np.array(M[:,0],dtype='uint')] += M[:,2]
M = cancer
mat[np.array(M[:,0], dtype='uint'), np.array(M[:,1],dtype='uint')] += M[:,2]
mat[np.array(M[:,1], dtype='uint'), np.array(M[:,0],dtype='uint')] += M[:,2]
#M = homology
#mat[np.array(M[:,0], dtype='uint'), np.array(M[:,1],dtype='uint')] += M[:,2]
#mat[np.array(M[:,1], dtype='uint'), np.array(M[:,0],dtype='uint')] += M[:,2]
mat[xrange(num), xrange(num)] = 1
sol_dream11('submitted/sub2_1000.txt', mat=mat, svd_k = 100, n_clusters=1000, save_feature=True, feature_filename='submitted/sub2_k100_i100_feature.pkl')"""