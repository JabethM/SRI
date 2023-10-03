import numpy as np
import random
import networkx as nx
import matplotlib.pyplot as plt


class SIRS():
    """
    This is the class where an arbitrary network is used as the environment for a disease which follows SIRS dynamics
    is simulated on. Upon initialising the network we need to define the number of nodes, the type of network we want
    to use (initialisation) and the probability of transition between states (p1, p2 and p3).
    """

    def __init__(self, nodes, initialisation=(0, 0.1), p1=0.5, p2=0.1, p3=0.015):
        """

        :param nodes: total number of nodes in our system :param initialisation: A tuple where the first parameter
        specifies the type of network and any additional parameters are values that may be needed for the
        construction of specified network type
            - initialisation = (0, a) - erdos reyni network where:
                                            'a' is a fixed probability that any given node is connected to another

            - initialisation = (1, a, b) - watts_strogatz network where:
                                            'a' is the number of neighbouring nodes each node is connected to in a ring
                                              topology,
                                            'b' is the probability of rewiring each edge

            - initialisation = (2, a) - barabasi_albert network where:
                                            'a' is the number edges to attach from a new node to an existing node

            - initialisation = (3, a, b) - stochastic_block network where:
                                            'a' is populations for each 'cluster' in the graph. The sum of the
                                              populations of all the clusters must sum to the total number of networks.
                                            'b' is a nxn symmetric matrix (where n is the number of clusters) that
                                              represents the probability of edge connections between nodes from one
                                              cluster to another

        :param p1: fixed probability that any given susceptible node that shares an edge with an infected node becomes
        infected as well

        :param p2: fixed probability that any given infected node recovers
        :param p3: fixed probability that any given recovered node becomes susceptible
        """
        self.dt = 0.1
        self.time = 0
        self.p1 = p1
        self.p2 = p2
        self.p3 = p3
        self.P = [self.p1, self.p2, self.p3]

        self.num_nodes = nodes
        self.G = None

        self.initialise_graph(initialisation)

        self.set_disease()
        self.pos = None

    def initialise_graph(self, init_tuple):
        # Random
        if init_tuple[0] == 0:
            self.G = nx.erdos_renyi_graph(self.num_nodes, init_tuple[1])
        elif init_tuple[0] == 1:
            self.G = nx.watts_strogatz_graph(self. num_nodes, init_tuple[1], init_tuple[2])
        # Scale Free
        elif init_tuple[0] == 2:
            self.G = nx.barabasi_albert_graph(self.num_nodes, init_tuple[1])
        # High Modularity
        elif init_tuple[0] == 3:
            assert (sum(init_tuple[1]) == self.num_nodes)
            assert (len(init_tuple[2]) == len(init_tuple[1]))
            self.G = nx.stochastic_block_model(init_tuple[1], init_tuple[2])

        return self.G

    def set_disease(self):
        """
        Sets the initial state of the system by assigning all nodes to a state of susceptibility aside from one which
        we assign as patient zero.

        :return: the instance network self.G
        """
        infected_node = random.randint(0, self.num_nodes - 1)
        state_dict = {i: 0 for i in range(self.num_nodes)}
        state_dict[infected_node] = 1

        nx.set_node_attributes(self.G, state_dict, name="state")
        return self.G

    def iterate(self):
        H = self.G.copy()
        probs = np.random.random(self.num_nodes)
        states = np.array([attribute.get("state") for node, attribute in H.nodes(data=True)])

        for i in range(len(states)):
            transform_bool = False

            if states[i] == 0:
                nbs = [n for n in H.neighbors(i)]
                for j in nbs:
                    if states[j] == 1:
                        transform_bool = True
                        break
            else:
                transform_bool = True

            if transform_bool and probs[i] < self.P[states[i]]:
                self.G.nodes[i]["state"] = (states[i] + 1) % 3

    def colors(self):
        colors = ['green', 'red', 'blue']
        node_colors = [colors[data.get("state")] for node, data in self.G.nodes(data=True)]
        return node_colors

    def low_dim_draw(self):
        node_colors = self.colors()
        nx.draw(self.G, with_labels=True, node_color=node_colors, edge_color='gray')
        plt.show()

    def mid_dim_draw(self):
        if self.pos is None:
            self.pos = nx.spring_layout(self.G, dim=3)
        test = sorted(self.G)
        node_xyz = np.array([self.pos[v] for v in sorted(self.G)])
        edge_xyz = np.array([(self.pos[u], self.pos[v]) for u, v in self.G.edges()])

        # Create the 3D figure
        fig = plt.figure()
        ax = fig.add_subplot(111, projection="3d")

        # Plot the nodes - alpha is scaled by "depth" automatically
        node_colors = self.colors()
        ax.scatter(*node_xyz.T, s=100, c=node_colors, ec="w")

        for vizedge in edge_xyz:
            ax.plot(*vizedge.T, color="tab:gray")

        # Turn gridlines off
        # ax.grid(False)
        # Suppress tick labels
        # for dim in (ax.xaxis, ax.yaxis, ax.zaxis):
        #     dim.set_ticks([])
        # Set axes labels
        ax.set_xlabel("x")
        ax.set_ylabel("y")
        ax.set_zlabel("z")

        fig.tight_layout()
        plt.show()
        return


def main():
    e = SIRS(50, 0)

    for i in range(10):
        e.iterate()
        e.mid_dim_draw()
