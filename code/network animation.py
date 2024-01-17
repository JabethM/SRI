import matplotlib;
import networkx as nx

matplotlib.use("TkAgg")

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from SIRS_single import SIRS
import mpl_toolkits.mplot3d.art3d as art3d

dt = None
sim = None
tMax = None
node_plt = None
node_xyz = None


def init():
    global dt
    global tMax
    global sim
    global node_plt
    node_colors = sim.colors()
    n_c = triple_array(node_colors)
    node_plt = v_color(node_plt, n_c)
    # for i in range(len(node_plt)):
    #     node_plt[i].set_color(node_colors[i])

    return node_plt


def animate(i):
    global sim
    global node_plt

    while i * dt < tMax:
        sim.iterate()
        node_colors = sim.colors()
        n_c = triple_array(node_colors)
        node_plt = v_color(node_plt, n_c)
        # for i in range(len(node_plt)):
        #    node_plt[i].set_color(node_colors[i])

        return node_plt


def run(nodes, typ, p1, p2, p3):
    global fig
    global ax
    global dt
    global tMax
    global sim
    global node_plt
    global node_xyz

    dt = 0.1
    tMax = 1000

    sim = SIRS(nodes, typ, variants=2, ps=[[p1, p2, p3], [0, 0, 0]])

    fig = plt.figure(figsize=(10, 6))

    ax = fig.add_subplot(111, projection="3d")

    sim.pos = nx.spring_layout(sim.G, dim=3)

    node_xyz = np.array([sim.pos[v] for v in sorted(sim.G)])
    edge_xyz = np.array([(sim.pos[u], sim.pos[v]) for u, v in sim.G.edges()])

    node_plt = np.array([[plt.Circle((node_xyz[i][0], node_xyz[i][1]), 0.02, color='tab:gray'),
                          plt.Circle((node_xyz[i][1], node_xyz[i][2]), 0.02, color='tab:gray'),
                          plt.Circle((node_xyz[i][0], node_xyz[i][2]), 0.02, color='tab:gray')]
                         for i in range(nodes)])

    for i in range(len(node_plt)):
        axes = ['x', 'y', 'z']
        for j in range(3):
            ax.add_patch(node_plt[i][j])
            art3d.pathpatch_2d_to_3d(node_plt[i][j], z=node_xyz[i][(j + 2) % 3], zdir=axes[(j + 2) % 3])

    for edge in edge_xyz:
        ax.plot(*edge.T, color='tab:gray')
    X_data = node_xyz[:, 0]
    Y_data = node_xyz[:, 1]
    Z_data = node_xyz[:, 2]

    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_zlabel("z")
    ratio = 5
    ax.set_box_aspect([np.ptp(X_data) * ratio, np.ptp(Y_data) * ratio, np.ptp(Z_data) * ratio])

    fig.canvas.mpl_connect('draw_event', updates)

    ani = animation.FuncAnimation(fig, animate,
                                  interval=0.00001, repeat=False, init_func=init)
    plt.show()
    return


def color(n_plt, col):
    n_plt.set_color(col)
    return n_plt


def v_color(n, c):
    n = np.array(n)
    return np.vectorize(color)(n, c)


def distance(p1, p2):
    return (p1 - p2) ** 2


def v_distance(p1, p2):
    return np.vectorize(distance)(p1, p2)


def set_ad(node, alpha):
    node.set_alpha(abs(alpha))
    return


def v_set_ad(n, a):
    return np.vectorize(set_ad)(n, a)


def triple_array(single):
    triple = 3
    z = np.shape(node_plt)

    repeat = np.tile(single, (1, triple))

    a = np.shape(single)
    if len(a) > 1:
        b = np.shape(single)[0] * 3
        c = np.shape(single)[1]

        repeat = np.reshape(repeat, (b, c))
        tripled = list(map(tuple, repeat))
    else:
        repeat = np.repeat(single, triple)
        tripled = repeat.reshape((np.shape(node_plt)))
    return tripled


def modified_sigmoid(x):
    modifier = 1.1
    sig = (1 / (modifier + np.exp(-103.0 * (x - 0.5)))) + ((1 - 1 / modifier) / 2)
    return sig


def update_opacities(event):
    global node_plt
    # Get the camera position in 3D space
    # cam_pos = ax.get_proj().proj_transform.transform([ax.azim, ax.elev, 1])

    alpha = ax.azim * np.pi / 180.
    beta = ax.elev * np.pi / 180.
    n = np.array([np.cos(alpha) * np.sin(beta),
                  np.sin(alpha) * np.cos(beta),
                  np.sin(beta)])

    p_d = v_distance(node_xyz, n)
    p_d = np.linalg.norm(p_d, axis=1)
    normal = lambda x: (x - np.min(p_d)) / (np.max(p_d) - np.min(p_d))
    v_norm = lambda l: np.vectorize(normal)(l)

    p_d = v_norm(p_d)
    p_d = 1 - p_d
    p_d = modified_sigmoid(p_d)
    p_d = triple_array(p_d)
    v_set_ad(node_plt, p_d)


def updates(event):
    update_opacities(event)


def run_configs(config):
    num_of_nodes = config.get("nodes", 50)
    initialisation = config.get("init", (0, 0.1))
    p1 = config.get("p1", 0.5)
    p2 = config.get("p2", 0.1)
    p3 = config.get("p3", 0.015)

    run(num_of_nodes, initialisation, p1, p2, p3)


if __name__ == '__main__':
    # probs = [[0.2, 0.015, 0.015, 0.015],
    #         [0.015, 0.2, 0.015, 0.015],
    #         [0.015, 0.015, 0.2, 0.015],
    #         [0.015, 0.015, 0.015, 0.2]]

    # configuration = {"nodes": 100, "init": (3, [20, 40, 27, 13], probs)}
    # print(probs)
    configuration = {"nodes": 100, "init": (2, 2)}
    # configuration = {"node": 50, "init": (1, 5, 0.1)}
    # configuration = {}
    run_configs(config=configuration)
