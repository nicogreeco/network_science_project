import networkx as nx
import csv
import sys
import gzip
import json
from multiprocessing import Pool

DEBUG = False

# This funcion reads the graph structure from the input file
def load_network_weighted(file_path):

    G = nx.Graph()
    with gzip.open(file_path, 'rt') as file:  # 'rt' mode for reading text
        for line in file:
            parts = line.strip().split()
            if len(parts) == 3:
                node1, node2, w = parts
                if node1 != 'protein1' and node2 != 'protein2':
                    G.add_edge(node1, node2, weight=1)
            else:
                continue  # skip lines that do not match the expected format

    return G

# This method keeps removing edges from Graph until one of the connected components of Graph splits into two
# compute the edge betweenness
def cmtyGirvanNewmanStep(G, num_cores):
    if DEBUG:
        print("Running cmtyGirvanNewmanStep method ...")
    init_ncomp = nx.number_connected_components(G)  # no of components
    ncomp = init_ncomp
    while ncomp <= init_ncomp:
        # Parallel edge betweenness centrality calculation
        bw = parallel_edge_betweenness_centrality(G, num_cores, weight='weight')
        # find the edge with max centrality
        max_ = max(bw.values())
        # find the edge with the highest centrality and remove all of them if there is more than one!
        for k, v in bw.items():
            if v == max_:
                G.remove_edge(*k)  # remove the central edge
        ncomp = nx.number_connected_components(G)  # recalculate the no of components

# Parallel calculation of edge betweenness centrality
def parallel_edge_betweenness_centrality(G, num_cores, weight=None):
    nodes = list(G.nodes())
    chunk_size = len(nodes) // num_cores
    
    node_chunks = [nodes[i:i + chunk_size] for i in range(0, len(nodes), chunk_size)]
    
    with Pool(processes=num_cores) as pool:
        results = pool.starmap(compute_betweenness, [(G, chunk, weight) for chunk in node_chunks])
    
    # Combine results from all chunks
    betweenness = results[0]
    for result in results[1:]:
        for edge in result:
            betweenness[edge] += result[edge]
    
    return betweenness

def compute_betweenness(G, nodes, weight):
    return nx.edge_betweenness_centrality_subset(G, nodes, list(G.nodes()), weight=weight)

# This method compute the modularity of current split
def girvanNewmanGetModularity(G, deg_, m_):
    New_A = nx.adjacency_matrix(G)
    New_deg = updateDeg(New_A, G.nodes())
    # Let's compute Q
    comps = nx.connected_components(G)  # list of components
    print('No of communities in decomposed G: {}'.format(nx.number_connected_components(G)))
    print('No of communities remaining edges in G: {}'.format(G.number_of_edges()))
    Mod = 0  # Modularity of a given partitioning
    for c in comps:
        EWC = 0  # no of edges within a community
        RE = 0  # no of random edges
        for u in c:
            EWC += New_deg[u]
            RE += deg_[u]  # count the probability of a random edge
        Mod += (EWC - RE * RE / (2 * m_))
    Mod = Mod / (2 * m_)
    if DEBUG:
        print("Modularity: {}".format(Mod))
    return Mod

def updateDeg(A, nodes):
    deg_dict = {}
    n = len(nodes)
    B = A.sum(axis=1)
    for i, node_id in enumerate(nodes):
        deg_dict[node_id] = B[i]
    return deg_dict

def set_default(obj):
    if isinstance(obj, set):
        return list(obj)
    raise TypeError

# This method runs Girvan-Newman algorithm and finds the best community split by maximizing the modularity measure
def runGirvanNewman(G, Orig_deg, m_, num_cores):
    # Let's find the best split of the graph
    BestQ = 0.0
    Q = 0.0
    while True:
        cmtyGirvanNewmanStep(G, num_cores)
        Q = girvanNewmanGetModularity(G, Orig_deg, m_)
        print("Modularity of decomposed G: {}".format(Q))
        if Q > BestQ:
            BestQ = Q
            Bestcomps = list(nx.connected_components(G))  # Best Split
            print("Identified components: {}".format(Bestcomps))
        if G.number_of_edges() == 0:
            break
    if BestQ > 0.0:
        print("Max modularity found (Q): {} and number of communities: {}".format(BestQ, len(Bestcomps)))
        print("Graph communities: {}".format(Bestcomps))
    else:
        print("Max modularity (Q):", BestQ)


## Function modified to output and save not only the highest mosularity partition but also 
## several partition with n of clusters similiar to target_clusters. As 
def runGirvanNewmanWithSizes(G, Orig_deg, m_, target_clusters, save_folder, num_cores):
    BestQ = 0.0
    Q = 0.0
    # Dictionary to store the community structures for desired cluster counts
    cluster_snapshots = {}
    target_clusters = sorted(target_clusters, reverse=True)  # Sort descending to handle thresholds properly

    while True:
        cmtyGirvanNewmanStep(G, num_cores)
        Q = girvanNewmanGetModularity(G, Orig_deg, m_)
        # print("Modularity of decomposed G: {}".format(Q))

        current_components = list(nx.connected_components(G))
        num_components = len(current_components)

        if Q > BestQ:
            BestQ = Q
            Bestcomps = current_components
            # print("Identified components: {}".format(Bestcomps))

        # Check if current number of components is close to or has just passed any target thresholds
        while target_clusters and num_components >= target_clusters[-1]:
            close_target = target_clusters.pop()
            cluster_snapshots[close_target] = current_components

            json_data = json.dumps(current_components, default=set_default)

            with open(f'{save_folder}/clusters_Girvan-Newman_{close_target}.json', 'w') as json_file:
                json_file.write(json_data)
            print(f"Snapshot at {close_target} clusters taken.")

        if G.number_of_edges() == 0 or not target_clusters:
            break

    if BestQ > 0.0:
        print("Max modularity found (Q): {} and number of communities: {}".format(BestQ, len(Bestcomps)))
        with open(f'{save_folder}/clusters_Girvan-Newman_Max_Modularity_{BestQ}.json', 'w') as json_file:
            json_file.write(json_data)
        print(f"Snapshot at {close_target} clusters taken.")
    else:
        print("Max modularity (Q):", BestQ)

    return cluster_snapshots


## Usage: python GirvanNewmanParallel.py <input graph.txt.gz> <num_cores>
def main(argv):
    if len(argv) < 3:
        sys.stderr.write("Usage: %s <input graph> <num_cores>\n" % (argv[0],))
        return 1
    graph_fn = argv[1]
    num_cores = int(argv[2])

    G = load_network_weighted(graph_fn)

    print(G, flush=True)

    n = G.number_of_nodes()  # |V|
    A = nx.adjacency_matrix(G)  # adjacency matrix

    m_ = 0.0  # the weighted version for the number of edges
    for i in range(0, n):
        for j in range(0, n):
            m_ += A[i, j]
    m_ = m_ / 2.0

    print("m: {}".format(m_))

    # calculate the weighted degree for each node
    Orig_deg = updateDeg(A, G.nodes())

    # run Girvan-Newman algorithm
    runGirvanNewmanWithSizes(G, Orig_deg, m_, [1000, 500, 350, 280, 275, 220], save_folder='/home/nicola/internship/PPI/511145_clustering', num_cores=num_cores)

if __name__ == "__main__":
    sys.exit(main(sys.argv))

