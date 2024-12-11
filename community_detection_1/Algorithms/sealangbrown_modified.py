import numpy as np
import networkx as nx

def SNcluster(Filename=None, graph=None, limit=10, groupAllNodes=True, mergeFileList=False, weighted=True, output_filename=None,directed=False):
    """Clusters nodes in a network edges list or NetworkX graph and returns output as a list of lists.
    
    Arguments:
    Filename -- Path to the input file containing the edges (optional if graph is provided).
    graph -- Direct input of the graph as a NetworkX graph.
    limit -- Maximum size of clusters.
    groupAllNodes -- Whether to attempt grouping all nodes into clusters.
    mergeFileList -- Whether to merge multiple files into one network (if Filename is a list).
    weighted -- Whether the graph is weighted.
    output_filename -- Optional path to save the output clusters.
    
    Returns:
    clusters -- List of lists, where each sublist contains nodes in a cluster.
    """

    def mergeLim(ng, nn, n1, n2, limit):
        """Determine if groups should merge based on a size penalty."""
        if n1 + n2 > limit:
            return False
        a = 140
        g = 0.2 * float(limit)
        return ng > nn * 0.5 * (1 + np.exp(-1 * np.power(-1 * (n1 + n2 - limit - 0.001) / a, g)))

    def mergeCol(r1, r2, A, G, W, weighted):
        """Merge column r2 into r1 and update matrices A, G, and W."""
        A[:, r1] += A[:, r2]
        A[np.nonzero(A[:, r1]), r1] = 1
        A[np.nonzero(A[:, r2]), r2] = 0
        G[:, r1] += G[:, r2]
        G[np.nonzero(G[:, r2]), r2] = 0
        if weighted:
            W[:, r1] += W[:, r2]
            W[np.nonzero(W[:, r1]), r1] = 1
            W[np.nonzero(W[:, r2]), r2] = 0
        return A, G, W

    def getOrder(A, G):
        """Get order of nodes to be processed based on their degree in A."""
        degrees = np.sum(A, axis=0)
        active_nodes = np.nonzero(np.diag(G))[0]
        return active_nodes[np.argsort(-degrees[active_nodes])]

    edges = []

    # If graph parameter is provided, use it to extract edges
    if graph:
        # Extract edges and weights from the NetworkX graph
        for u, v, data in graph.edges(data=True):
            weight = data['weight'] if 'weight' in data else 1.0
            edges.append([u, v, weight])
    else:
        # Import edges from file
        if mergeFileList:  # Merge multiple files into one network
            for file in Filename:
                InFile = file + ".txt"
                with open(InFile, 'rt') as f:
                    for line in f:
                        edges.append(line.strip().split('\t'))
                print(f"Processing file: {file}")
        else:
            InFile = Filename + ".txt"
            with open(InFile, 'rt') as f:
                for line in f:
                    edges.append(line.strip().split('\t'))
            print(f"Processing file: {Filename}")

    # Read file or graph, count nodes, create number/name mapping
    nodes = 0
    nameToNum = dict()
    numToName = dict()

    for row in edges:
        for nName in row[0:2]:  # Get 1st and 2nd entries
            if nName not in nameToNum:
                nameToNum[nName] = nodes
                numToName[nodes] = nName
                nodes += 1

    A = np.eye(nodes)  # Create adjacency matrix A - initialize as Identity
    W = np.eye(nodes)  # Weight matrix
    G = np.eye(nodes)  # Membership matrix

    # Read edges into adjacency matrix
    for row in edges:
        try:
            i = nameToNum[row[0]]
            j = nameToNum[row[1]]
            weight = float(row[2])
        except IndexError:
            print(f"Warning: Input file format exception. Possibly blank or incomplete row. Row: {row}")
        else:
            if not directed:
                A[min(i, j), max(i, j)] = 1  # Add edge
                W[min(i, j), max(i, j)] = weight 
            else:
                A[i,j]=1
                W[i,j]=weight

    groupCount, nodeCount = np.shape(A)
    checkNodes = np.nonzero(np.diag(G))  # Get active nodes

    # Step 1: First join degree 1 nodes
    degreeOneNodes = True
    while degreeOneNodes:
        degreeOneNodes = False
        for r in np.nditer(checkNodes):
            r0 = int(r)
            if np.sum(A[r0, :]) + np.sum(A[:, r0]) == 3:  # If degree 1
                R = A[:, r0] - G[:, r0]  # Get external neighbors
                try:
                    nbr = int(np.nonzero(R)[0])  # Find index of neighbor
                except TypeError:
                    nbr = int(max(np.nonzero(A[r0, :])[0]))  # Find index of neighbor
                if 1 + np.sum(A[:, nbr]) <= limit:
                    A, G, W = mergeCol(nbr, r0, A, G, W, weighted)
                    checkNodes = np.nonzero(np.diag(G))
                    degreeOneNodes = True
                    break

    # Step 2: Perform other merges
    lastGroupCount = -55        
    while lastGroupCount != groupCount:
        lastGroupCount = groupCount
        checkNodes = np.nonzero(np.diag(G))  # Only check nodes that have not been merged
        for r in np.nditer(checkNodes):
            r1 = int(r)
            merge = True  
            while merge:
                merge = False
                mb1 = G[:, r1]  # Get membership of r1
                r1n = A[:, r1] - mb1  # Get external neighbors r1
                n1 = np.sum(mb1)  # Sum nodes in r1
                
                if np.sum(r1n) > 0:  # Prevent iterating over empty
                    for rb in np.nditer(np.nonzero(r1n)):
                        r2 = int(rb)
                        r2n = A[:, r2]  # Get all neighbors of r2
                        mb2 = G[:, r2]  # Get all members of r2
                        n2 = np.sum(mb2)  # Sum nodes in r2
                        if weighted:
                            w2n = W[:, r2]
                            ng = (np.sum(w2n * r2n * (r1n + mb1)) + np.sum(W[:, r1] * r2n * (r1n + mb1))) / 2
                            nn = np.sum(w2n * (r2n - mb2))
                        else:
                            ng = np.sum(r2n * (r1n + mb1))  # Count neighbors that are in group
                            nn = np.sum(r2n - mb2)  # Total neighbors minus corresponding row in G
                        
                        if mergeLim(ng, nn, n1, n2, limit):
                            A, G, W = mergeCol(r1, r2, A, G, W, weighted)
                            groupCount -= 1
                            merge = True
                            break

    if groupAllNodes:  # Merge smaller groups into larger nodes
        singleNodes = True
        while singleNodes:
            checkNodes = getOrder(A, G)
            for r in np.nditer(checkNodes):
                singleNodes = False
                r1 = int(r)
                mb1 = G[:, r1]
                n1 = np.sum(mb1)  # Sum nodes in r1
                if n1 <= 3 and n1 > 0:
                    lastNg = 0
                    lnbr = 0
                    R = A[:, r1] - mb1  # Get external neighbors
                    if np.sum(R) > 0:
                        for ngbr in np.nditer(np.nonzero(R)):
                            r2 = int(ngbr)
                            r2n = A[:, r2]
                            mb2 = G[:, r2]
                            n2 = np.sum(mb2)
                            ng = np.sum(r2n * mb1)  # Neighbors that are in group
                            if ng > lastNg and n1 + n2 <= limit:
                                nbr = r2  # Most common neighbor
                                singleNodes = True
                            elif n2 > lnbr and n1 + n2 <= limit and ng == lastNg:
                                nbr = r2
                                singleNodes = True
                            lastNg = ng
                            lnbr = n2

                    if singleNodes:
                        A, G, W = mergeCol(nbr, r1, A, G, W, weighted)
                        groupCount -= 1
                        checkNodes = getOrder(A, G)

    clusters = []
    checkNodes = np.nonzero(np.diag(G))
    for col in np.nditer(checkNodes):
        colG = G[:, int(col)]
        if sum(colG) > 0:
            cluster_nodes = [numToName[int(nNum)] for nNum in np.nditer(np.nonzero(colG))]
            clusters.append(cluster_nodes)
    
    # Write to file if filename is provided
    if output_filename:
        with open(output_filename, 'w', newline='') as f:
            for gNum, cluster in enumerate(clusters, 1):
                f.write(f"{gNum}\t1.0\t" + "\t".join(cluster) + "\n")

    print(f"{len(clusters)} Groups created from {nodeCount} nodes.")
    return [[int(j) for j in i] for i in clusters]
