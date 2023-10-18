import copy

import numpy as np
import networkx as nx

nodes = 10
state_dict = {i: np.random.randint(0, 2, (3)) for i in range(nodes)}

G = nx.erdos_renyi_graph(nodes, p=0.5)
nx.set_node_attributes(G, state_dict, name="state")
current_states = np.array([attribute.get("state") for node, attribute in G.nodes(data=True)])
a = [set() for v in range(3)]
b = copy.deepcopy(a)
b[1] = {2, 3, 4}

print(a)
print(b)
