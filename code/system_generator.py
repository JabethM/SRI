import numpy as np
import random
import networkx as nx
import matplotlib.pyplot as plt


class network():
    def __init__(self, nodes, initialisation=0, p1=0.5, p2=0.1, p3=0.015):
        self.dt = 0.1
        self.time = 0
        self.p1 = p1
        self.p2 = p2
        self.p3 = p3
        self.P = [self.p1, self.p2, self.p3]

        self.num_nodes = nodes
        self.G = None

        if initialisation == 0:
            self.G = nx.erdos_renyi_graph(self.num_nodes, 0.1)
        elif initialisation == 1:
            self.G = nx.erdos_renyi_graph(self.num_nodes, 0.5)

        self.set_disease()
        self.pos = None

    def set_disease(self):
        random_state = [random.randint(0, 2) for i in range(self.num_nodes)]
        state_dict = {i: random_state[i] for i in range(self.num_nodes)}
        nx.set_node_attributes(self.G, state_dict, name="state")

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
    e = network(50, 0)

    for i in range(10):
        e.iterate()
        e.mid_dim_draw()

