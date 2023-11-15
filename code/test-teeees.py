import copy

import numpy as np
import networkx as nx

G = nx.erdos_renyi_graph(10, 0.1)
state_dict = {i: np.zeros(3) for i in range(10)}
nx.set_node_attributes(G, state_dict, name="state")
nononodes = [0,1]
nodes = [1,2,3]
b = G.nodes(data=True)[3]["state"][nononodes]
print(b)

