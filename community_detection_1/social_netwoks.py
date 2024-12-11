import numpy as np
import networkx as nx
import pandas as pd
import re
from collections import defaultdict
import time
import score_benchmarks
import community_detection_main
import json
from collections import defaultdict
def parse_graph_file(file_path):
    with open(file_path, 'r') as file:
        content = file.read()

    # Find all nodes in the file
    nodes = re.findall(r'node\s+\[\s+id\s+(\d+)\s+label\s+"([^"]+)"\s+value\s+(\d+)\s+\]', content)

    # Create a dictionary to store the node information
    node_info = {}
    for node in nodes:
        node_id, label, value = node
        node_info[int(node_id)] = {
            'label': label,
            'value': int(value)
        }

    return node_info


def parse_graph_file1(file_path):
    with open(file_path, 'r') as file:
        content = file.read()

    # Use regex to find all nodes in the file
    nodes = re.findall(r'node\s*\[\s*id\s+(\d+)\s+label\s+"([^"]+)"\s+value\s+(\d+)\s+source\s+"([^"]+)"\s*\]', content)

    # Create a dictionary to store the node information
    node_info = {}
    for node in nodes:
        node_id, label, value, source = node
        node_info[int(node_id)] = {
            'label': label,
            'value': int(value),
            'source': source
        }

    return node_info




def group_nodes_by_community(file_path,source=False):
    if source:
        node_info=parse_graph_file1(file_path)
    else:
        node_info=parse_graph_file(file_path)
    communities = defaultdict(list)
    
    for node_id, info in node_info.items():
        communities[info['value']].append(node_id)
    community_list = list(communities.values())
    
    return community_list
#print(group_nodes_by_community(file_path))



# Define the path to your text file
file_path = 'C:/Users/hp/community_detection/Social_networks/Email_bw_1000_people_dig_ground_truth.txt'

# Dictionary to store communities
communities_dict = {}

# Read the file and populate the dictionary
with open(file_path, 'r') as file:
    for line in file:
        node, community_id = map(int, line.split())
        if community_id not in communities_dict:
            communities_dict[community_id] = []
        communities_dict[community_id].append(node)

# Convert the dictionary to a list of lists
communities_e = list(communities_dict.values())

file_path = 'C:/Users/hp/community_detection/Social_networks/cora_ground_truth.txt'

# Dictionary to store communities
communities_dict = {}

# Read the file and populate the dictionary
with open(file_path, 'r') as file:
    for line in file:
        node, community_id = map(int, line.split())
        if community_id not in communities_dict:
            communities_dict[community_id] = []
        communities_dict[community_id].append(node)

# Convert the dictionary to a list of lists
communities_e1 = list(communities_dict.values())
print(communities_e1)


networks_1=['adjnoun','dolphins','football','polbooks','polblogs','karate','Email_bw_1000_people_dig','cora']
networks_1_=['internet','power']
networks_1_2=['musae_git']
networks_2=['cond-mat-weighted','lesmis-weighted','hep-weighted','netscience-weighted','astro-dig-weighted']
networks_4=['celegansneural']
networks_3=['Email-Enron-dig','Wikivote-dig']
root_directory='C:/Users/hp/community_detection/Social_networks/'

ground_truths_1=[group_nodes_by_community(f'{root_directory}{networks_1[i]}.gml') if i!=4 else group_nodes_by_community(f'{root_directory}{networks_1[i]}.gml',source=True)  for i in range(len(networks_1)-2)]
ground_truths_1.append(communities_e)
ground_truths_1.append(communities_e1)

#ground_truths_1.append(communities_g)

#triple_ahc(t1=0.1,t2=1000,weighted=True,directed=False,network_file='',G=None)
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
nmi_scores=np.zeros((8,16))
ccbar_scores=np.zeros((8,16))
times=np.zeros((8,16))
"""#triple_ahc(t1=0.1,t2=1000,weighted=True,directed=False,network_file='',G=None)

# triple ahc
for i in range(len(networks_1)):
    network_file_curr=f'{root_directory}{networks_1[i]}.txt'
    #print(i)
    G=process_network_file(network_file_curr,weighted=False,directed=False)
    #print(ground_truths_1[i])
    pc,t=community_detection_main.triple_ahc(t1=0.1,t2=1000,weighted=False,directed=False,network_file='',G=G)
    nmi_scores[i][0]=score_benchmarks.nmi_score(ground_truths_1[i],pc)
    times[i][0]=t
    ccbar_scores[i][0]=len(pc)/len(ground_truths_1[i])
print(nmi_scores[:,0])"""
# zhenhua(walktrap)
"""for i in range(len(networks_1)):
    network_file_curr=f'{root_directory}{networks_1[i]}.txt'
    #print(i)
    l_curr=[len(j) for j in ground_truths_1[i]]
    G=process_network_file(network_file_curr,weighted=False,directed=False)
    #print(ground_truths_1[i])
    pc,t=community_detection_main.zhenhua(max_limit=max(l_curr),method=1,weighted=False,directed=False,network_file='',G=G)
    print(len(pc))
    nmi_scores[i][1]=score_benchmarks.nmi_score(ground_truths_1[i],pc)
    times[i][1]=t
    ccbar_scores[i][1]=len(pc)/len(ground_truths_1[i])
print(nmi_scores[:,1])"""

"""# zhenhua(infomap)
for i in range(len(networks_1)):
    network_file_curr=f'{root_directory}{networks_1[i]}.txt'
    #print(i)
    l_curr=[len(j) for j in ground_truths_1[i]]
    G=process_network_file(network_file_curr,weighted=False,directed=False)
    t1=time.time()
    #print(ground_truths_1[i])
    pc=community_detection_main.zhenhua(max_limit=max(l_curr),method=2,weighted=False,directed=False,network_file='',G=G)
    t2=time.time()
    nmi_scores[i][2]=score_benchmarks.nmi_score(ground_truths_1[i],pc)
    times[i][2]=t2-t1
    ccbar_scores[i][2]=len(pc)/len(ground_truths_1[i])
print(nmi_scores[:,2])"""

"""# nextmr
for i in range(len(networks_1)):
    network_file_curr=f'{root_directory}{networks_1[i]}.txt'
    if i!=5:
        l_curr=[len(j) for j in ground_truths_1[i]]
        G=process_network_file(network_file_curr,weighted=False,directed=False)
        #print(ground_truths_1[i])
        pc,t=community_detection_main.nextmr(min_limit=min(l_curr),weighted=False,directed=False,network_file='',G=G)
        nmi_scores[i][2]=score_benchmarks.nmi_score(ground_truths_1[i],pc)
        times[i][2]=t
        ccbar_scores[i][2]=len(pc)/len(ground_truths_1[i])
print(nmi_scores[:,2])"""

# bluegenes
"""for i in range(0,len(networks_1)):
    print(i)
    if i!=4 and i!=6:
        network_file_curr=f'{root_directory}{networks_1[i]}.txt'
        l_curr=[len(j) for j in ground_truths_1[i]]
        G=process_network_file(network_file_curr,weighted=False,directed=False)
        #print(ground_truths_1[i])
        if i==6:
            print(len(G.nodes))
        pc,t=community_detection_main.blue_genes(min_limit=min(l_curr),alpha=1,weighted=False,directed=False,network_file='',G=G)
        nmi_scores[i][3]=score_benchmarks.nmi_score(ground_truths_1[i],pc)
        times[i][3]=t
        ccbar_scores[i][3]=len(pc)/len(ground_truths_1[i])
print(nmi_scores[:,3])"""
# tuskdmi
"""for i in range(3,len(networks_1)):
    if i!=4 and i!=6:
        network_file_curr=f'{root_directory}{networks_1[i]}.txt'
        print(i,'hi')
        l_curr=[len(j) for j in ground_truths_1[i]]
        G=process_network_file(network_file_curr,weighted=False,directed=False)
        #print(G.nodes)
        #print(ground_truths_1[i])
        if i==6:
            print(G.nodes)
        pc,t=community_detection_main.tuskdmi(num_com=len(l_curr),G=G,network_file='',weighted=False,directed=False)
        nmi_scores[i][3]=score_benchmarks.nmi_score(ground_truths_1[i],pc)
        times[i][3]=t
        ccbar_scores[i][3]=len(pc)/len(ground_truths_1[i])
print(nmi_scores[:,3])
for i in range(len(networks_1)):
        #if i!=4 and i!=6 :
        network_file_curr=f'{root_directory}{networks_1[i]}.txt'
        #print(i,'hi')
        l_curr=[len(j) for j in ground_truths_1[i]]
        G=process_network_file(network_file_curr,weighted=False,directed=False)
        #print(G.nodes)
        #print(ground_truths_1[i])
        pc,t=community_detection_main.tsuromi_ono_local(network_file='',threads=4,max_limit=len(l_curr),min_limit=len(l_curr),G=G,weighted=False,directed=False)
        nmi_scores[i][4]=score_benchmarks.nmi_score(ground_truths_1[i],pc)
        times[i][4]=t
        ccbar_scores[i][4]=len(pc)/len(ground_truths_1[i])
print(nmi_scores[:,4])
for i in range(len(networks_1)):
        #if i!=4 and i!=6 :
        network_file_curr=f'{root_directory}{networks_1[i]}.txt'
        #print(i,'hi')
        l_curr=[len(j) for j in ground_truths_1[i]]
        G=process_network_file(network_file_curr,weighted=False,directed=False)
        #print(G.nodes)
        #print(ground_truths_1[i])
        pc,t=community_detection_main.csbio_iitm_louvain(network_file='',resolution = 0.1,G=G,weighted=False,directed=False)
        nmi_scores[i][5]=score_benchmarks.nmi_score(ground_truths_1[i],pc)
        times[i][5]=t
        ccbar_scores[i][5]=len(pc)/len(ground_truths_1[i])
print(nmi_scores[:,5])
for i in range(len(networks_1)):
        #if i!=4 and i!=6 :
        network_file_curr=f'{root_directory}{networks_1[i]}.txt'
        #print(i,'hi')
        l_curr=[len(j) for j in ground_truths_1[i]]
        G=process_network_file(network_file_curr,weighted=False,directed=False)
        #print(G.nodes)
        #print(ground_truths_1[i])
        pc,t=community_detection_main.sealangbrown(G=G,network_file='',limit=50,weighted=False,directed=False)
        nmi_scores[i][6]=score_benchmarks.nmi_score(ground_truths_1[i],pc)
        times[i][6]=t
        ccbar_scores[i][6]=len(pc)/len(ground_truths_1[i])
print(nmi_scores[:,6])
for i in range(0,len(networks_1)):
    if i!=5 and i!=1:
        network_file_curr=f'{root_directory}{networks_1[i]}.txt'
        print(i,'hi')
        l_curr=[len(j) for j in ground_truths_1[i]]
        G=process_network_file(network_file_curr,weighted=False,directed=False)
        #print(G.nodes)
        #print(ground_truths_1[i])
        pc,t=community_detection_main.tianle(G=G,weighted=True,network_file='',directed=False,n_clusters=len(l_curr))
        nmi_scores[i][7]=score_benchmarks.nmi_score(ground_truths_1[i],pc)
        times[i][7]=t
        ccbar_scores[i][7]=len(pc)/len(ground_truths_1[i])
print(nmi_scores[:,7])

for i in range(1,len(networks_1)):
    if i!=0 and i!=4:
        network_file_curr=f'{root_directory}{networks_1[i]}.txt'
        print(i,'hi')
        l_curr=[len(j) for j in ground_truths_1[i]]
        G=process_network_file(network_file_curr,weighted=False,directed=False)
        #print(G.nodes)
        #print(ground_truths_1[i])
        pc,t=community_detection_main.big_s2(G=G,weighted=False,network_file='',directed=False,min_size=min(l_curr),max_size=max(l_curr), Max_iter=15, M=5)
        print(pc)
        nmi_scores[i][8]=score_benchmarks.nmi_score(ground_truths_1[i],pc)
        times[i][8]=t
        ccbar_scores[i][8]=len(pc)/len(ground_truths_1[i])
print(nmi_scores[:,8])
for i in range(0,len(networks_1)):
        #if i!=0 and i!=4:
        network_file_curr=f'{root_directory}{networks_1[i]}.txt'
        print(i,'hi')
        l_curr=[len(j) for j in ground_truths_1[i]]
        G=process_network_file(network_file_curr,weighted=False,directed=False)
        #print(G.nodes)
        #print(ground_truths_1[i])
        pc,t=community_detection_main. teamcs_aleph(G=G,weighted=False,directed=False,network_file='',recursive=True)
        #print(pc)
        nmi_scores[i][9]=score_benchmarks.nmi_score(ground_truths_1[i],pc)
        times[i][9]=t
        ccbar_scores[i][9]=len(pc)/len(ground_truths_1[i])
print(nmi_scores[:,9])
for i in range(0,len(networks_1)):
    if i!=2 :
        network_file_curr=f'{root_directory}{networks_1[i]}.txt'
        print(i,'hi')
        l_curr=[len(j) for j in ground_truths_1[i]]
        G=process_network_file(network_file_curr,weighted=False,directed=False)
        #print(G.nodes)
        #print(ground_truths_1[i])
        pc,t=community_detection_main. sim_net(G=G,weighted=False,directed=False,network_file='',groupNumber=len(l_curr))
        print(len(pc))
        print(len(ground_truths_1[i]))
        nmi_scores[i][10]=score_benchmarks.nmi_score(ground_truths_1[i],pc)
        times[i][10]=t
        ccbar_scores[i][10]=len(pc)/len(ground_truths_1[i])
print(nmi_scores[:,10])"""

for i in range(0,len(networks_1)):
        #if i!=0:
        network_file_curr=f'community_detection/Social_networks/{networks_1[i]}.txt'
        print(i,'hi')
        l_curr=[len(j) for j in ground_truths_1[i]]
        if i not in [6,7]:
                G=process_network_file(network_file_curr,weighted=False,directed=False)
                pc,t=community_detection_main.luminex(weighted=False,directed=False,network_file=network_file_curr,G=G,p=1)
        #print(len(pc))
        #print(len(ground_truths_1[i]))
        else:
                G=process_network_file(network_file_curr,weighted=False,directed=True)
                pc,t=community_detection_main.luminex(weighted=False,directed=True,network_file=network_file_curr,G=G,p=1)
        nmi_scores[i][11]=score_benchmarks.nmi_score(ground_truths_1[i],pc)
        times[i][11]=t
        ccbar_scores[i][11]=len(pc)/len(ground_truths_1[i])
#print(nmi_scores[:,11])
print(nmi_scores[:,11])
print(ccbar_scores[:,11])
print(times[:,11])
"""
for i in range(0,len(networks_1)):
        #if i!=0:
        network_file_curr=f'{root_directory}{networks_1[i]}.txt'
        print(i,'hi')
        l_curr=[len(j) for j in ground_truths_1[i]]
        G=process_network_file(network_file_curr,weighted=False,directed=False)
        #print(G.nodes)
        #print(ground_truths_1[i])
        pc,t=community_detection_main.spectral_clustering(weighted=False,directed=False,network_file='',G=G,n_clusters=len(l_curr))
        
        nmi_scores[i][12]=score_benchmarks.nmi_score(ground_truths_1[i],pc)
        times[i][12]=t
        ccbar_scores[i][12]=len(pc)/len(ground_truths_1[i])
print(nmi_scores[:,12])

for i in range(0,len(networks_1)):
        #if i!=0:
        network_file_curr=f'{root_directory}{networks_1[i]}.txt'
        print(i,'hi')
        l_curr=[len(j) for j in ground_truths_1[i]]
        G=process_network_file(network_file_curr,weighted=False,directed=False)
        #print(G.nodes)
        #print(ground_truths_1[i])
        pc,t=community_detection_main.walktrap(G=G,weighted=False,directed=False,network_file='',steps=10)
        
        nmi_scores[i][13]=score_benchmarks.nmi_score(ground_truths_1[i],pc)
        times[i][13]=t
        ccbar_scores[i][13]=len(pc)/len(ground_truths_1[i])
print(nmi_scores[:,13])"""
"""for i in range(0,len(networks_1)):
        #if i!=0:
        network_file_curr=f'{root_directory}{networks_1[i]}.txt'
        print(i,'hi')
        l_curr=[len(j) for j in ground_truths_1[i]]
        G=process_network_file(network_file_curr,weighted=False,directed=False)
        #print(G.nodes)
        #print(ground_truths_1[i])
        pc,t=community_detection_main.label_propogation(G=G,weighted=False,directed=False,network_file='')
        if i==0:
            print(pc)
        nmi_scores[i][14]=score_benchmarks.nmi_score(ground_truths_1[i],pc)
        times[i][14]=t
        ccbar_scores[i][14]=len(pc)/len(ground_truths_1[i])
print(nmi_scores[:,14])"""

"""for i in range(0,len(networks_1)):
    print(i)
    if i!=4 and i!=6:
        network_file_curr=f'{root_directory}{networks_1[i]}.txt'
        l_curr=[len(j) for j in ground_truths_1[i]]
        G=process_network_file(network_file_curr,weighted=False,directed=False)
        #print(ground_truths_1[i])
        if i==6:
            print(len(G.nodes))
        pc,t=community_detection_main.blue_genes(min_limit=min(l_curr),alpha=1,weighted=False,directed=False,network_file='',G=G)
        nmi_scores[i][15]=score_benchmarks.nmi_score(ground_truths_1[i],pc)
        times[i][15]=t
        ccbar_scores[i][15]=len(pc)/len(ground_truths_1[i])
    else:
        nmi_scores[i][15]=None
        times[i][15]=None
        ccbar_scores[i][15]=None
print(nmi_scores[:,15])"""