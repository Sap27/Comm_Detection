import csbioiitm.community_csbio as community
import networkx as nx
def compute_modularity(datafile, resolution_parameter,G=None):
    if G is None:
        G = nx.Graph()  # Construct networkx graph
        for line in datafile:
            g = line.strip().split("\t")
            G.add_edge(g[0], g[1], weight=float(g[2]))

    mod_partition = community.best_partition(G, resolution=float(resolution_parameter))
    modularity_value = community.modularity(mod_partition, G)
    
    return G, mod_partition, modularity_value

from collections import defaultdict

def identify_core_modules(G, partition, output_path=None):
    dict_of_comm = defaultdict(list)
    for node, comm_id in partition.items():
        dict_of_comm[comm_id].append(node)

    results = []
    ii = 0
    overall_node_count = 0

    for i, nodes in dict_of_comm.items():
        if len(nodes) < 3:
            overall_node_count += len(nodes)
        elif len(nodes) > 100:
            core_comm = {}
            for kk in nodes:
                edgelist = G.edges(kk)
                overall_degree = len(edgelist)
                indegree = sum(1 for edges in edgelist if edges[0] in nodes and edges[1] in nodes)
                outdegree = overall_degree - indegree
                core_comm[kk] = outdegree
            core_community = sorted(core_comm.items(), key=lambda x: x[1])
            
            new_core = [k for k, v in core_community if v < 50]
            results.append((ii+1, new_core))
            
        else:
            results.append((ii+1, nodes))
        
        ii += 1

    if output_path:
        with open(output_path, 'w') as f:
            for res in results:
                f.write(f"{res[0]}\t{','.join(res[1])}\n")
    return results

# Example usage:
# G, partition, _ = compute_modularity(open('network.dat', 'r'), 0.1)
# result = identify_core_modules(G, partition)
# print(result)

# Example usage:
# with open('network.dat', 'r') as file:
#     G, partition, modularity_score = compute_modularity(file, 0.1)
#     print("Modularity score:", modularity_score)




