import networkx as nx
import numpy as np
import subprocess
import tempfile
import os

def rescale(weights, graph, maximum, minimum):
    numV = graph.number_of_nodes()
    return ((numV - 1) / (maximum - minimum)) * (weights - minimum) + 1

def preProcessing(input_data, method='pageRank', quantile_level=0.75, integerWeight='no'):
    input_data[:, :2] += 1  # Adjust node numbering to start from 1
    
    if method == 'quantile':
        threshold = np.quantile(input_data[:, 2], quantile_level)
        input_data = input_data[input_data[:, 2] >= threshold]
    elif method in ['pageRank', 'double']:
        graph = nx.Graph()
        graph.add_weighted_edges_from(input_data[:, :3])
        pagerank_scores = nx.pagerank(graph, weight='weight')
        threshold = np.quantile(list(pagerank_scores.values()), quantile_level)
        input_data = input_data[input_data[:, 2] >= threshold]
        
        if method == 'double':
            threshold = np.quantile(input_data[:, 2], quantile_level)
            input_data = input_data[input_data[:, 2] >= threshold]

    graph = nx.Graph()
    graph.add_weighted_edges_from(input_data[:, :3])
    
    if integerWeight == 'yes':
        max_weight = input_data[:, 2].max()
        min_weight = input_data[:, 2].min()
        scaled_weights = rescale(input_data[:, 2], graph, max_weight, min_weight)
        for idx, (u, v, data) in enumerate(graph.edges(data=True)):
            data['weight'] = scaled_weights[idx]
    elif integerWeight == 'no':
        for u, v, data in graph.edges(data=True):
            data['weight'] = float(data['weight'])
    else:
        for u, v, data in graph.edges(data=True):
            data['weight'] = 1
    
    return graph

def run_mlrmcl(b, c, i, input_file, output_file):
    

    cmd = f"wsl ./community_detection/mlrmcl -b {b} -c {c} -i {i} -o {output_file} {input_file}"
    res = subprocess.run(cmd, capture_output=True, text=True, shell=True)
    print(res.stdout)
    if res.returncode != 0:
        print(f"Error running mlrmcl: {res.stderr}")
        return None
    return output_file

def postProcessing(clusters, method='random', smallest=3, largest=100, graph=None, b2=None, c2=None, i2=None):
    filtered_clusters = [c for c in clusters if len(c) >= smallest]
    
    if method == 'random':
        large_clusters = [c for c in filtered_clusters if len(c) > largest]
        small_clusters = [c for c in filtered_clusters if len(c) <= largest]
        
        for cluster in large_clusters:
            chunks = [cluster[i:i + largest] for i in range(0, len(cluster), largest)]
            small_clusters.extend(chunks)
        
        return small_clusters
    elif method == 'discard':
        return [c for c in filtered_clusters if len(c) <= largest]
    elif method == 'recluster':
        return reclusterLargeClusters(filtered_clusters, graph, b2, c2, i2, smallest, largest)
    else:
        return filtered_clusters

def reclusterLargeClusters(clusters, graph, b2, c2, i2, smallest, largest):
    large_clusters = [c for c in clusters if len(c) > largest]
    small_clusters = [c for c in clusters if len(c) <= largest]
    
    if not large_clusters:
        return clusters
    
    result = []
    for cluster in large_clusters:
        subgraph = graph.subgraph(cluster).copy()
        
        # Use a temporary file for subgraph input data
        with tempfile.NamedTemporaryFile(delete=False, mode='w', suffix='.txt') as tmp_input:
            tmp_input_path = tmp_input.name
            tmp_input.write(f"{subgraph.number_of_nodes()}\t{subgraph.number_of_edges()}\t1\n")
            for u, v, data in subgraph.edges(data=True):
                tmp_input.write(f"{u}\t{v}\t{data['weight']}\n")
        
        # Use a temporary file for subgraph output data
        with tempfile.NamedTemporaryFile(delete=False, mode='w', suffix='.txt') as tmp_output:
            tmp_output_path = tmp_output.name
        
        run_mlrmcl(b2, c2, i2, input_file=tmp_input_path, output_file=tmp_output_path)
        
        subgraph_result = np.loadtxt(tmp_output_path, dtype=int)
        
        subgraph_clusters = [np.where(subgraph_result == i)[0].tolist() for i in np.unique(subgraph_result)]
        subgraph_clusters = [c for c in subgraph_clusters if len(c) > smallest]
        
        small_clusters.extend([c for c in subgraph_clusters if len(c) <= largest])
        result.extend([c for c in subgraph_clusters if len(c) > largest])
        
        os.remove(tmp_input_path)  # Clean up temporary input file
        os.remove(tmp_output_path)  # Clean up temporary output file
    
    return small_clusters + reclusterLargeClusters(result, graph, b2, c2, i2, smallest, largest)

def main(file, b=1.3, c=500, i=1, filter='pageRank', threshold=0.75, inteWeight='yes', weighted=True, dir='output_clusters', post='random', smallest=3, largest=100, b2=None, c2=None, i2=None):
    input_data = np.loadtxt(file, delimiter='\t')
    
    # Step 1: Preprocess the input data to create a graph
    graph = preProcessing(input_data, method=filter, quantile_level=threshold, integerWeight=inteWeight)
    
    # Step 2: Generate a file formatted for the community detection algorithm
    """with tempfile.NamedTemporaryFile(delete=False, mode='w', suffix='.txt') as temp_file:
        filename = temp_file.name
        
        # Write number of nodes and edges
        temp_file.write(f"{graph.number_of_nodes()}\t{graph.number_of_edges()}\t{1 if weighted else 0}\n")
        
        # Iterate over the edges and write them to the file
        for u, v, data in graph.edges(data=True):
            if weighted:
                # If the graph is weighted, write the weight of the edge
                weight = data.get('weight', 1)  # Default weight to 1 if not present
                temp_file.write(f"{u}\t{v}\t{weight}\t")
            else:
                # If the graph is not weighted, just write the edge
                temp_file.write(f"{u}\t{v}\t")"""
    #print(len(graph.nodes))
    with open("community_detection/test_file.txt","w") as f:
        f.write(f"{graph.number_of_nodes()}\t{graph.number_of_edges()}\t{1 if weighted else 0}\n")
        for u, v, data in graph.edges(data=True):
            if weighted:
                # If the graph is weighted, write the weight of the edge
                weight = data.get('weight', 1)  # Default weight to 1 if not present
                f.write(f"{u}\t{v}\t{weight}\n")
            else:
                # If the graph is not weighted, just write the edge
                f.write(f"{u}\t{v}\n")

    # Step 3: Run the mlrmcl community detection algorithm
    """with tempfile.NamedTemporaryFile(delete=False, mode='w', suffix='.txt') as tmp_output:
        output_file_path = tmp_output.name"""
    output_file_path='community_detection/output3.txt'
    
    run_mlrmcl(b, c, i, input_file="community_detection/test_file.txt", output_file=output_file_path)
    
    # Read the result from the temporary output file
    result = np.loadtxt(output_file_path, dtype=int)
    
    # Convert result to list of lists (clusters)
    clusters = [np.where(result == i)[0].tolist() for i in np.unique(result)]
    
    # Step 4: Post-process the mlrmcl output to refine the clusters
    output = postProcessing(clusters, method=post, smallest=smallest, largest=largest, graph=graph, b2=b2, c2=c2, i2=i2)
    
    # Print the final clusters instead of writing to a file
    print("Final Clusters:")
    for cluster in output:
        print(cluster)
    
    # Clean up temporary files
    #os.remove(filename)
    os.remove(output_file_path)
    
    return output

# Define the parameters
file = "C:/Users/hp/community_detection/network.dat"
b = 1.3
c = 1000
i = 1.8
filter = "pageRank"  # Use "pageRank" or "double" as needed
threshold = 0.75
inteWeight = "yes"
weighted = True
dir = "output_clusters"
post = "random"  # Use "discard" or "recluster" as needed
smallest = 3
largest = 100
b2 = 2
c2 = 5000
i2 = 2

# Call the main function and print the output
output = main(file, b, c, i, filter, threshold, inteWeight, weighted, dir, post, smallest, largest, b2, c2, i2)
print(output)
