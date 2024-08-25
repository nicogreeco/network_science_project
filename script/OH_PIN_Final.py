import networkx as nx
import gzip
import json
import numpy as np
import datetime
from itertools import combinations
import sys
import os
import heapq

def count_edges(adj_matrix):
    edge_count = sum(sum(row) for row in adj_matrix)
    return edge_count // 2  # Each edge is counted twice in an undirected grap

def load_network_to_matrix(file_path):
    # Step 1: Identify all unique nodes
    nodes = set()
    with gzip.open(file_path, 'rt') as f:
        for i, line in enumerate(f):
            if i == 0:
                continue  # Skip the first line
            node1, node2, w = line.strip().split()
            nodes.add(node1)
            nodes.add(node2)
    
    # Step 2: Create a mapping from node to index and index to node
    node_list = list(nodes)
    node_index = {node: idx for idx, node in enumerate(node_list)}
    index_node = {idx: node for node, idx in node_index.items()}
    n = len(node_list)
    
    # Step 3: Initialize the adjacency matrix
    adj_matrix = [[0] * n for _ in range(n)]
    
    # Step 4: Fill the adjacency matrix with edges
    with gzip.open(file_path, 'rt') as f:
        for i, line in enumerate(f):
            if i == 0:
                continue  # Skip the first line
            node1, node2, w = line.strip().split()
            i, j = node_index[node1], node_index[node2]
            adj_matrix[i][j] = 1
            adj_matrix[j][i] = 1
    
    return adj_matrix, index_node


# Step 1


## Convert a list of sets containing indexes to a list of sets containing nodes.
def convert_index_sets_to_node_sets(index_sets, index_node):
    node_sets = [{index_node[idx] for idx in index_set} for index_set in index_sets]
    return node_sets


## Generate a B_cluster for a given edge, including the edge's vertices
## and their common neighbors. Optimized to minimize set operations.
def generate_B_cluster_3(adj_matrix, edge):

    u_idx, v_idx = edge
    
    # Get neighbors of u and v from the adjacency matrix
    neighbors_u = {idx for idx, is_neighbor in enumerate(adj_matrix[u_idx]) if is_neighbor}
    neighbors_v = {idx for idx, is_neighbor in enumerate(adj_matrix[v_idx]) if is_neighbor}
    
    # Utilize intersection directly for performance improvement
    return {u_idx, v_idx}.union(neighbors_u & neighbors_v)

## Initializes M_clusters based on the description of the algorithm.
def initialize_M_clusters_3(adj_matrix):

    C_set = []  # List of sets to keep track of clusters
    n = len(adj_matrix)
    
    # Iterate over the upper triangle of the adjacency matrix to find edges
    for i in range(n):
        for j in range(i + 1, n):
            if adj_matrix[i][j] == 1:
                new_B_cluster = generate_B_cluster_3(adj_matrix, (i, j))
                new_C_set = []
                for cluster in C_set:
                    if not cluster <= new_B_cluster:
                        new_C_set.append(cluster)
                
                # Add the new B_cluster if it is not a subset of any cluster in the current C_set
                if not any(new_B_cluster <= cluster for cluster in new_C_set):
                    new_C_set.append(new_B_cluster)
                
                C_set = new_C_set

    return C_set


# Step 2

## Calculate the overlapping score OS between two clusters, optimized to use precomputed sizes.
def overlapping_score(cluster1, cluster2):

    intersection_size = len(cluster1 & cluster2)
    return (intersection_size ** 2) / (len(cluster1) * len(cluster2))

## Optimizes the merging of clusters by reducing the number of set operations and improving the merge loop logic.
def merge_clusters_4(C_set, os_threshold=0.5):
    # Initial score calculation for all pairs
    scores_heap = []
    for cluster1, cluster2 in combinations(C_set, 2):
        score = overlapping_score(cluster1, cluster2)
        if score > os_threshold:
            # Push scores as negative to create a max-heap using heapq which is by default a min-heap
            heapq.heappush(scores_heap, (-score, cluster1, cluster2))

    while scores_heap:
        # Pop the pair with the highest score
        max_score, cluster1, cluster2 = heapq.heappop(scores_heap)
        max_score = -max_score  # Convert back to positive score
        print(f'current max score: {max_score}' )

        # Merge the two clusters with the highest score
        merged_cluster = cluster1 | cluster2
        new_C_set = [cluster for cluster in C_set if cluster != cluster1 and cluster != cluster2]
        new_C_set.append(merged_cluster)

        # Update the scores heap: remove scores involving the merged clusters
        new_scores_heap = []
        for score, c1, c2 in scores_heap:
            if c1 != cluster1 and c1 != cluster2 and c2 != cluster1 and c2 != cluster2:
                new_scores_heap.append((score, c1, c2))
        scores_heap = new_scores_heap  # Reassign the filtered heap

        # Recreate the heap since the structure might be destroyed
        heapq.heapify(scores_heap)

        # Calculate scores for the new merged cluster with all other clusters
        for cluster in new_C_set[:-1]:  # Exclude the newly added merged cluster itself
            new_score = overlapping_score(merged_cluster, cluster)
            if new_score > os_threshold:
                heapq.heappush(scores_heap, (-new_score, merged_cluster, cluster))

        C_set = new_C_set

        if not scores_heap:
            break

    return C_set

# Step 3

def set_default(obj):
    if isinstance(obj, set):
        return list(obj)
    raise TypeError


def weighted_clustering_coefficient_3(adj_matrix, cluster1, cluster2, neighbors):

    actual_weight_sum = 0
    for node1 in cluster1:
        for node2 in cluster2:
            if adj_matrix[node1, node2] == 1:
                common_neighbors = len(neighbors[node1] & neighbors[node2])
                actual_weight_sum += common_neighbors
    possible_edges = len(cluster1) * len(cluster2)
    return 0.0 if possible_edges == 0 else actual_weight_sum / possible_edges


## Calculate the lambda_C value for a cluster in an adjacency matrix.
def calculate_lambda_C_3(adj_matrix, cluster):
    internal_edges = np.sum(adj_matrix[np.ix_(list(cluster), list(cluster))]) / 2

    # Calculate external edges: sum the adjacency matrix over the cluster and the rest of the graph
    # Subtract internal edges (already counted in the cluster's sum) to get only external connections.
    total_edges_to_cluster = np.sum(adj_matrix[list(cluster), :])  # Sum all edges connected to nodes in the cluster
    external_edges = total_edges_to_cluster - 2 * internal_edges  # Subtract double the internal edges
    return 0.0 if external_edges == 0 else internal_edges / external_edges

## Final hierarchial Step
def assemble_into_lambda_modules_4(id, adj_matrix, C_list, lambda_thresholds, index_node):
    adj_matrix = np.array(adj_matrix)
    clusters = {i: set(clust) for i, clust in enumerate(C_list)}
    neighbors = {i: set(np.nonzero(adj_matrix[i])[0]) for i in range(len(adj_matrix))}
    
    # Initialize lambda values for each cluster
    lambda_values = {i: calculate_lambda_C_3(adj_matrix, clust) for i, clust in clusters.items()}

    print(f'Lambda values are {lambda_thresholds}')
    for lambda_th in lambda_thresholds:
        CCV = {(i, j): weighted_clustering_coefficient_3(adj_matrix, clusters[i], clusters[j], neighbors)
               for i in clusters for j in clusters if i < j}

        # Initialize max-heap
        heap = [(-val, key) for key, val in CCV.items()]
        heapq.heapify(heap)

        while heap and -heap[0][0] > 0:
            max_val, max_pair = heapq.heappop(heap)
            max_val = -max_val  # Convert back to positive value
            i, j = max_pair
            if CCV.get((i, j), -1) == -1:
                continue  # Skip if already processed
            if lambda_values[i] >= lambda_th and lambda_values[j] >= lambda_th:
                CCV[max_pair] = -1  # Mark as processed
                print(f'Clusters: {i} and {j} lambdas >= {lambda_th}', flush=True)
            else:
                print(f'Max clustering coeff: {max_val}, clusters: {i} and {j} merged', flush=True)
                clusters[i].update(clusters[j])  # Merge clusters
                lambda_values[i] = calculate_lambda_C_3(adj_matrix, clusters[i])  # Update lambda value for the new merged cluster
                del clusters[j]  # Remove merged cluster
                del lambda_values[j]  # Remove lambda value for the removed cluster
                CCV = {key: val for key, val in CCV.items() if j not in key}
                for k in list(clusters):
                    if k != i:
                        if clusters[k] <= clusters[i]:
                            del clusters[k]
                            del lambda_values[k]  # Also remove lambda value for deleted cluster
                            CCV = {key: val for key, val in CCV.items() if k not in key}
                        else:
                            new_key = tuple(sorted((i, k)))
                            new_val = weighted_clustering_coefficient_3(adj_matrix, clusters[i], clusters[k], neighbors)
                            CCV[new_key] = new_val
                            heapq.heappush(heap, (-new_val, new_key))  # Push new value to max-heap

        print(f'Lambda threshold: {lambda_th}, number of clusters: {len(clusters)}')

        # Save to JSON file
        clusters_list = [list(v) for v in clusters.values()]
        clusters_list_names = convert_index_sets_to_node_sets(clusters_list, index_node)
        json_data = json.dumps(clusters_list_names, default=set_default)
        with open(f'./{id}_new_OHPIN_clusters_{lambda_th}.json', 'w') as json_file:
            json_file.write(json_data)

    return [list(v) for v in clusters.values()]



# total function

def OH_PIN(id, graph, lambda_thresholds, index_node, os_threshold=0.5):
    M_clusters = initialize_M_clusters_3(graph)
    print(f'Step 1 completed, {len(M_clusters)} M clusters', flush=True)
    C_set=merge_clusters_4(M_clusters, os_threshold=os_threshold)
    print(f'Step 2 completed, {len(C_set)} non overlapping (%50) clusters', flush=True)
    C_set=assemble_into_lambda_modules_4(id, graph, C_set, lambda_thresholds=lambda_thresholds, index_node=index_node)
    return(C_set)

# Main

## Usage: python OH_PIN_Final.py <input graph.txt.gz>
def main(argv):
    if len(argv) < 2:
        sys.stderr.write("Usage: %s <input graph>" % (argv[0],))
        return 1
    graph_path = argv[1]   
    print(datetime.datetime.now(), flush=True)
    adj_matrix, index_node = load_network_to_matrix(graph_path)
    print(f'Network "{graph_path}" has been loaded', flush=True)
    print(f'The graph has {len(adj_matrix)} nodes and {count_edges(adj_matrix)} edges', flush=True)
    print(f'Starting the clustering', flush=True)
    id=os.path.basename(graph_path).split('.')[0]
    C_set=OH_PIN(id, adj_matrix, lambda_thresholds=[0.05, 0.1, 0.2, 0.4, 1], index_node=index_node)
    print(f'Total number of clusters: {len(C_set)}')
    print(datetime.datetime.now())
    

    
if __name__ == "__main__":
    sys.exit(main(sys.argv))
    