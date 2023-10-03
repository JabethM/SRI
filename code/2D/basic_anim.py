import matplotlib;

matplotlib.use("TkAgg")

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from twoD_SIRS import sim

fig = None
ax = None
board_plt = None
dt = None

dim = None
sim_instance = None
tMax = None
sim_type = None


def init():
    global dt
    global tMax
    global sim_instance
    global board_plt

    dt = 0.1
    tMax = 1000

    sim_instance = sim(dt, (dim, dim), sim_type)
    board = sim_instance.board
    board_plt = v_colouring(board, board_plt)

    return board_plt


def animate(i):
    global board_plt
    while i * dt < tMax:
        board = sim_instance.iterate()
        board_plt = v_colouring(board, board_plt)
        return board_plt


def run(d=50, typ=0):
    global fig
    global ax
    global dim
    global dt
    global board_plt
    global sim_type

    fig, ax = plt.subplots()
    sim_type = typ
    dim = d
    board_plt = [
        [plt.Rectangle(xy=(i, j), width=1, height=1, fill=True, color='tab:gray') for i in range(dim)] for j in
        range(dim)]

    board_plt = np.array(board_plt)

    f = lambda x: ax.add_patch(x)
    vf = lambda x: np.vectorize(f)(x)

    vf(board_plt)

    ax.set_aspect('equal')
    ax.set_xlim(0, dim)
    ax.set_ylim(0, dim)
    ani = animation.FuncAnimation(fig, animate,
                                  interval=0.00001, repeat=False, init_func=init)
    plt.show()
    return


def colouring(cell, patch):
    colour_options = ['g', 'r', 'w']
    a = cell
    patch.set_color(colour_options[int(cell)])
    return patch


def v_colouring(c, p):
    return np.vectorize(colouring)(c, p)

run()
