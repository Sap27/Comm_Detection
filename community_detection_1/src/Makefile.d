ClustNSeeToDream.o: ClustNSeeToDream.c Utils.h Partition.h
Compute.o: Compute.c Utils.h Graph.h Partition.h Louvain.h \
 EdgesComposition.h
EdgesComposition.o: EdgesComposition.c EdgesComposition.h Partition.h \
 Graph.h Utils.h
ExtractPartition.o: ExtractPartition.c Utils.h Partition.h
ExtractSubGraph.o: ExtractSubGraph.c Utils.h Graph.h Partition.h
FilterGraph.o: FilterGraph.c Utils.h Graph.h Partition.h
FilterPartition.o: FilterPartition.c Utils.h Partition.h
Graph.o: Graph.c Graph.h Utils.h Partition.h
Louvain.o: Louvain.c Louvain.h Utils.h Graph.h Partition.h
Partition.o: Partition.c Partition.h Utils.h
StandardToClustNSee.o: StandardToClustNSee.c Utils.h Partition.h
Utils.o: Utils.c Utils.h
