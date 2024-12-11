import networkx as nx
import community
import sys
###construct graph
G=nx.Graph() #Construct networkx graph
datafile=open(sys.argv[1],'r')
for line in datafile:
    g=line.split("\t")
    G.add_edges_from([(g[0],g[1])],weight=float(g[2][:-1]))

parameter=sys.argv[2]
###Modularity
modpartit = community.best_partition(G,resolution=float(parameter))
modularit, mod = community.modularity(modpartit,G)
#print "Modularity score:",round(modularit,2)
for i in modpartit.keys():
    print (str(i)+"\t"+str(modpartit[i]))
