import numpy as np

#########################################################################
###-----------------2.CLUSTERING ANALYSIS - K-MEANS-------------------###
#########################################################################

def k_means(data, k, max_iterations=100):
    # Randomly initialize centroids
    centroids = data[np.random.choice(data.shape[0], k, replace=False)]

    for _ in range(max_iterations):
        # Assign each data point to the closest centroid
        distances = np.linalg.norm(data - centroids[:, np.newaxis], axis=2)
        labels = np.argmin(distances, axis=0)

        # Update centroids based on the mean of assigned points
        new_centroids = np.array([data[labels == i].mean(axis=0) for i in range(k)])

        # Check for convergence
        if np.all(centroids == new_centroids):
            break

        centroids = new_centroids

    return labels, centroids



def cluster_gene_extraction(cleaned_data, labels, output_file):
    # Find the cluster with 102 elements
    desired_cluster_size = 102

    # Count the number of members in each cluster
    unique_labels, counts = np.unique(labels, return_counts=True)
    cluster_sizes = dict(zip(unique_labels, counts))

    # Find the cluster label with the desired size
    chosen_cluster_label = [cluster for cluster, size in cluster_sizes.items() if size == desired_cluster_size]

    # Check if a cluster with the desired size is found
    if not chosen_cluster_label:
        print(f"No cluster found with {desired_cluster_size} members.")
    else:
        # Use the first cluster found with the desired size
        chosen_cluster_label = chosen_cluster_label[0]

        # Extract gene names from the chosen cluster
        genes_in_cluster = cleaned_data[cleaned_data.iloc[:, -1] == chosen_cluster_label].iloc[:, 1].tolist()

        # Write gene names to the specified output file
        with open(output_file, 'w') as file:
            for gene_name in genes_in_cluster:
                file.write(f"{gene_name}\n")

        print(f"Gene names in the cluster with {desired_cluster_size} members have been written to {output_file}.")



def count_rows(file_path):
    with open(file_path, 'r') as file:
        num_rows = sum(1 for line in file)
    return num_rows
