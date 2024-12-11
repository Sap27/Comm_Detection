import sys
import networkx as nx
from collections import OrderedDict
import operator
ff=open(sys.argv[2])
G=nx.Graph()
for i in ff:
    g=i.split("\t")
    G.add_edges_from([(g[0],g[1])])
f=open(sys.argv[1])
from collections import defaultdict
dictofcomm=defaultdict(list)
for j in f:
    jj=j.split("\t")
    dictofcomm[jj[1][:-1]].append(jj[0])
ii=0
overallnodecount=0
for i in dictofcomm.keys():
    if len(dictofcomm[i])<3:
        overallnodecount += len(dictofcomm[i])
    elif len(dictofcomm[i])>100:
        corecomm={}
        for kk in dictofcomm[i]:
           edgelist=G.edges(kk)
           overall_degree=len(G.edges(kk))
           indegree=0
           for edges in edgelist:
               if edges[0] in dictofcomm[i] and edges[1] in dictofcomm[i]:
                   indegree += 1
           outdegree=overall_degree-indegree
           corecomm[kk]=outdegree
        corecommunity = sorted(corecomm.items(), key=lambda x: x[1])
        print (ii+1,"\t",1,"\t")
        newcore=[]
        for j in corecommunity:#community size
            if j[1]< 50:#out-degree cutoff
                newcore.append(j[0])
        diffofbothcore=len(corecommunity)-len(newcore)
        overallnodecount += diffofbothcore
        newcore2=[]
        if len(newcore) < 3:
            for j in corecommunity:#community size
                if j[1]< 100:#out-degree cutoff
                    newcore2.append(j[0])
            newcore=newcore2
        if len(newcore) < 3:
            newcore=[]
            for j in corecommunity[:10]:
                newcore.append(j[0])
        newcore1=[] 
        if len(newcore) > 100:
            for j in corecommunity:#community size
                if j[1]< 25:#out-degree cutoff
                    newcore1.append(j[0])
            newcore=newcore1
        if len(newcore) > 100:
            newcore=[]
            for j in corecommunity[:50]:
                newcore.append(j[0])
        for k in newcore:               
            print (k,"\t")
        print ('\n')
        ii=ii+1
    else:
        print (ii+1,"\t",1,"\t")
        for j in dictofcomm[str(i)]:
            print (j,"\t")
        print ('\n')
        ii=ii+1

