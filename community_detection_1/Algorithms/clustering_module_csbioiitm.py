import networkx as nx
import community_detection_1.community_pertub as community

# Load Graph from file


# Calculate modularity and get community partitions
def calculate_modularity(G, resolution):
    modpartit = community.best_partition(G, resolution=resolution)
    modularit = community.modularity(modpartit, G)
    return modpartit, modularit

# Identify core modules based on the described criteria
def identify_core_modules(G, modpartit):
    core_modules = []
    overallnodecount = 0

    for community_id in set(modpartit.values()):
        nodes = [node for node, comm_id in modpartit.items() if comm_id == community_id]
        if len(nodes) < 3:
            overallnodecount += len(nodes)
        elif len(nodes) > 100:
            corecomm = {}
            for node in nodes:
                edgelist = G.edges(node)
                overall_degree = len(G.edges(node))
                indegree = sum(1 for edge in edgelist if edge[0] in nodes and edge[1] in nodes)
                outdegree = overall_degree - indegree
                corecomm[node] = outdegree

            corecommunity = sorted(corecomm.items(), key=lambda x: x[1])
            newcore = [node for node, outdegree in corecommunity if outdegree < 50]

            if len(newcore) < 3:
                newcore = [node for node, outdegree in corecommunity if outdegree < 100]

            if len(newcore) < 3:
                newcore = [node for node, _ in corecommunity[:10]]

            if len(newcore) > 100:
                newcore = [node for node, outdegree in corecommunity if outdegree < 25]

            if len(newcore) > 100:
                newcore = [node for node, _ in corecommunity[:50]]

            core_modules.append(newcore)
            overallnodecount += len(corecommunity) - len(newcore)
        else:
            core_modules.append(nodes)

    return core_modules

# Main function that ties everything together
def main(resolution=0.1, input_file='',output_modularity=None, output_core_modules=None,G=None):
    
    """if input_file!='':
        G = G = nx.read_edgelist(input_file, data=(('weight', float),), nodetype=int)"""

    modpartit, modularit = calculate_modularity(G, resolution)
    core_modules = identify_core_modules(G, modpartit)

    if output_modularity:
        with open(output_modularity, 'w') as f:
            for node, community in modpartit.items():
                f.write(f"{node}\t{community}\n")
    else:
        modularity_output = [[node, community] for node, community in modpartit.items()]
        #print(f"Modularity: {modularit}")
        #print("Modularity output:", modularity_output)

    if output_core_modules:
        with open(output_core_modules, 'w') as f:
            for idx, module in enumerate(core_modules, start=1):
                f.write(f"{idx}\t1\t" + "\t".join(module) + "\n")
    return [[int(j) for j in i] for i in core_modules]

if __name__ == "__main__":
    import sys
    if len(sys.argv) < 4:
        #print("Usage: python clustering_module.py <graph_file> <resolution> [output_modularity_file] [output_core_modules_file]")
        sys.exit(1)

    graph_file = sys.argv[1]
    resolution = float(sys.argv[2])
    output_modularity_file = sys.argv[3] if len(sys.argv) > 3 else None
    output_core_modules_file = sys.argv[4] if len(sys.argv) > 4 else None

    main(graph_file, resolution, output_modularity_file, output_core_modules_file)
