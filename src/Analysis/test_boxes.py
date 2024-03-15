import networkx as nx
import numpy as np

# Sample adjacency matrix (replace this with your actual adjacency matrix)
adjacency_matrix = np.array([[0, 1, 1, 0],
                             [1, 0, 1, 0],
                             [1, 1, 0, 1],
                             [0, 0, 1, 0]])

# Step 1: Calculate the degree of each node
degrees = np.sum(adjacency_matrix, axis=1)

# Step 2: Construct arrays representing the degrees of connected nodes
deg1, deg2 = [], []
for i in range(adjacency_matrix.shape[0]):
    for j in range(i + 1, adjacency_matrix.shape[0]):
        if adjacency_matrix[i, j] == 1:
            deg1.append(degrees[i])
            deg2.append(degrees[j])

# Step 3: Calculate the degree assortativity using the Pearson correlation coefficient
degree_assortativity_np = np.corrcoef(deg1, deg2)[0, 1]

print("Degree assortativity coefficient (using numpy):", degree_assortativity_np)

G = nx.from_numpy_array(adjacency_matrix)

# Calculate the degree assortativity
degree_assortativity = nx.degree_assortativity_coefficient(G)

print("Degree assortativity coefficient (using nx):", degree_assortativity)


def degree_assortativity_coefficient(adj_matrix):
    degrees = np.sum(adj_matrix, axis=1)
    total_degree = np.sum(degrees)
    degrees_squared = np.sum(degrees ** 2)
    sum_edges = np.trace(adj_matrix)
    sum_degree_product = np.sum(np.outer(degrees, degrees))

    assortativity = (sum_edges / total_degree) - (sum_degree_product / (total_degree ** 2))

    return assortativity


degree_assortativity_own = degree_assortativity_coefficient(adjacency_matrix)
print("Degree Assortativity coefficient (own): ", degree_assortativity_own)
