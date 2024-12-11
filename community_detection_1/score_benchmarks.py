import numpy as np
import community_detection_1.realcommunities_ as rc
"""def count(l1,l2):
    s=0
    for i in l1:
        if i in l2:
            s+=1
    return s
def nmi_score(realcommunities,predictedcommunities):
    confusionmatrix=[]
    for i in range(len(realcommunities)):
        t=[]
        for j in range(len(predictedcommunities)):
            t.append(count(realcommunities[i],predictedcommunities[j]))
        confusionmatrix.append(t)
    
    confusionmatrix=np.array(confusionmatrix)
    #print(confusionmatrix[29])
    nmi=0
    num=0
    n = sum(len(community) for community in realcommunities)
    print(n)
    confusionmatrixt=np.transpose(confusionmatrix)
    for i in range(len(realcommunities)):
        for j in range(len(predictedcommunities)):
            if confusionmatrix[i][j]!=0:
                num+=confusionmatrix[i][j]*(np.log(confusionmatrix[i][j]*n/(confusionmatrix[i][:].sum()*confusionmatrixt[j][:].sum())))
    d1=0
    for i in range(len(realcommunities)):
        if confusionmatrix[i][:].sum() !=0:
        #print(len(realcommunities[i]))
            d1+=(confusionmatrix[i][:].sum())*(np.log(confusionmatrix[i][:].sum()/n))

    d2=0
    for j in range(len(predictedcommunities)):
        if confusionmatrixt[j][:].sum()!=0:
            d2+=(confusionmatrixt[j][:].sum())*(np.log(confusionmatrixt[j][:].sum()/n))
    nmi=-2*(num/(d1+d2))
    return nmi"""
import numpy as np

def count(l1, l2):
    """Count the number of common elements in two lists."""
    return sum(1 for i in l1 if i in l2)

def nmi_score(real_communities, predicted_communities):
    # Number of elements (nodes)
    
    
    # Build the confusion matrix
    confusion_matrix = []
    for real_community in real_communities:
        t = []
        for predicted_community in predicted_communities:
            t.append(count(real_community, predicted_community))
        confusion_matrix.append(np.array(t))
    
    confusion_matrix = np.array(confusion_matrix)
    confusion_matrix_t = np.transpose(confusion_matrix)
    #n = sum(len(community) for community in predicted_communities)
    n = sum(sum(ele) for ele in confusion_matrix)
    #print(n1,n)
    # Calculate the numerator
    num = 0
    for i in range(len(real_communities)):
        for j in range(len(predicted_communities)):
            if confusion_matrix[i][j] != 0:
                num += (confusion_matrix[i][j] * 
                        np.log(confusion_matrix[i][j] * n / 
                               (confusion_matrix[i][:].sum() * confusion_matrix_t[j][:].sum())))
    
    # Calculate the denominators
    d1 = 0
    for i in range(len(real_communities)):
        if confusion_matrix[i][:].sum() != 0:
            d1 += (confusion_matrix[i][:].sum() * 
                   np.log(confusion_matrix[i][:].sum() / n))

    d2 = 0
    for j in range(len(predicted_communities)):
        if confusion_matrix_t[j][:].sum() != 0:
            d2 += (confusion_matrix_t[j][:].sum() * 
                   np.log(confusion_matrix_t[j][:].sum() / n))
    
    # Compute NMI
    nmi = -2 * (num / (d1 + d2))
    return nmi

"""def parse_cluster_file(file_path):
    clusters = []  # This will store the list of lists
    current_cluster = []

    with open(file_path, 'r') as file:
        for line in file:
            line = line.strip()
            
            if line.startswith("ClusterID"):
                if current_cluster is not []:
                    clusters.append(current_cluster)  # Save the current cluster
                current_cluster = []  # Start a new cluster
            
            elif line and line != '||':  # Avoid empty lines or the '||' separator
                current_cluster.append(int(line))
        
        # Append the last cluster after the loop ends
        if current_cluster is not []:
            clusters.append(current_cluster)
    
    return clusters

# Example usage
file_path = 'community_detection/src/outpu1.txt'
communities = parse_cluster_file(file_path)
print(communities)
nmi=score_benchmarks(rc.groundtruth('Mu_0.10/community.dat'),communities)


"""