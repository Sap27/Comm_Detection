from infomap import Infomap
import networkx as nx
def do_infomap(in_file='',weighted=True,directed=False,out_file='',G=None):
    if  G is None:
        G = G = nx.DiGraph() if directed else nx.Graph()
        datafile=open(in_file, 'r')
        for line in datafile:
            g = line.strip().split("\t")
            G.add_edge(g[0], g[1], weight=float(g[2]))
    if weighted and not directed:
        im = Infomap("--two-level")  # Note: No --directed flag
        for source, target,weight in G.edges(data='weight'):
            im.add_link(int(source), int(target), float(weight))
    elif weighted and directed:
        im = Infomap("--two-level --directed")  # Note: No --directed flag
        for source, target,weight in G.edges(data='weight'):
            im.add_link(int(source), int(target), float(weight))  
    elif not weighted and not directed:
        im = Infomap("--two-level")  # Note: No --directed flag
        for source, target in G.edges():
            im.add_link(int(source), int(target))     
    elif not weighted and directed:
        im = Infomap("--two-level --directed")  # Note: No --directed flag
        for source, target in G.edges():
            im.add_link(int(source), int(target)) 
    im.run()
    module_dict = {}
    for node in im.tree:
        if node.is_leaf:
            module_id = node.module_id
            if module_id not in module_dict:
                module_dict[module_id] = []
            module_dict[module_id].append(node.node_id)
    communities = list(module_dict.values())
    if out_file!='':
        with open(out_file, 'a') as file:
            for j in range(Ncla):
                tmp_nodes = np.intersect1d(block_names[j], np.unique(col1)) 
                if len(tmp_nodes) <= Max_size and len(tmp_nodes) >= Min_size:
                    file.write(f"{j} 0.5 {' '.join(map(str, tmp_nodes))}\n")

    return communities
"""G1=nx.karate_club_graph()
communities=do_infomap(G=G1,weighted=False,directed=False)
print(communities)"""