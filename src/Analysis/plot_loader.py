import argparse
import pickle
import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import interp1d


def adjust_lightness(color, amount=0.5):
    import matplotlib.colors as mc
    import colorsys
    try:
        c = mc.cnames[color]
    except:
        c = color
    c = colorsys.rgb_to_hls(*mc.to_rgb(c))
    return colorsys.hls_to_rgb(c[0], max(0, min(1, amount * c[1])), c[2])


def read_pickle_file(file_path, regime):
    """
    regime 0: dataFiles
    regime 1: connectionStatsData
    :param file_path:
    :param regime:
    :return:
    """
    with open(file_path, 'rb') as f:
        data = pickle.load(f)
    if regime == 0:
        return data['x'], data['y']
    if regime == 1:
        return data['time'], data['stats']


def plot_PlotPickle_files(file_paths):
    plt.figure(figsize=(8, 6))
    mean_x = np.array([])
    mean_y = np.array([])
    func = None
    base_colors = ['tab:red', 'tab:orange', 'tab:blue', 'tab:green']
    for f, file_path in enumerate(file_paths):
        plot_x, plot_y = read_pickle_file(file_path, 0)
        number_of_variants = plot_y.shape[0] // 2

        if f == 0:
            mean_x = plot_x
            mean_y = plot_y
        else:
            func = interp1d(plot_x, plot_y, axis=-1, fill_value='extrapolate')
            interpolated_y = func(mean_x)
            mean_y += interpolated_y

        for v in range(2 * number_of_variants):
            chosen_colour = base_colors[v]
            plt.plot(plot_x, plot_y[v], color=adjust_lightness(chosen_colour, 1.7))
    number_of_variants = mean_y.shape[0] // 2
    mean_y = mean_y / (f + 1)
    for v in range(number_of_variants):
        letter = (chr(ord('A') + v))
        plt.plot(mean_x, mean_y[v], color=adjust_lightness(base_colors[v], 0.8),
                 label=f"Infected mean for variant {letter}")
        plt.plot(mean_x, mean_y[v + number_of_variants],
                 color=adjust_lightness(base_colors[v + number_of_variants], 0.8),
                 label=f"Susceptible mean for variant {letter}")

    plt.xlabel('Time')
    plt.ylabel('Number of Infected')
    plt.title('Infectious Curve plot')
    plt.legend()
    plt.grid(True)
    plt.show()


def plot_ConnectionPickle_files(connection_paths):
    plt.figure(figsize=(8, 6))
    avg_time = np.array([])

    avg_most = np.array([])
    avg_uq = np.array([])
    avg_mean = np.array([])
    avg_lq = np.array([])
    avg_least = np.array([])
    func = None

    for c, c_path in enumerate(connection_paths):
        connection_x, connection_y = read_pickle_file(c_path, 1)

        def interpolate_connections(old_x, new_x, new_y, axis):
            func = interp1d(new_x, new_y, axis=axis, fill_value='extrapolate')
            interpolated_y = func(old_x)
            return interpolated_y

        if c == 0:
            avg_time = np.array(connection_x)
            avg_most, avg_uq, avg_mean, avg_lq, avg_least = np.array(connection_y)


        else:
            for a, avg in enumerate([avg_most, avg_uq, avg_mean, avg_lq, avg_least]):
                avg += interpolate_connections(avg_time, connection_x, connection_y[a], 0)


    c_colours = ['tab:blue', 'tab:orange']
    styles = ['-.', '--', '-']

    length = c + 1
    avg_most = avg_most / length
    avg_uq = avg_uq / length
    avg_mean = avg_mean / length
    avg_lq = avg_lq / length
    avg_least = avg_least / length

    fig, ax1 = plt.subplots()
    ax1.set_xlabel('Time')
    ax1.set_ylabel('Number of Connections')

    for var in range(2):
        ax1.plot(avg_time, avg_most[:, var], color=adjust_lightness(c_colours[var], 1),
                 ls=styles[0])
        ax1.plot(avg_time, avg_mean[:, var], color=adjust_lightness(c_colours[var], 1),
                 ls=styles[1])
        ax1.plot(avg_time, avg_least[:, var], color=adjust_lightness(c_colours[var], 1),
                 ls=styles[2])

    colour_names = ['Variant A', 'Variant B']
    for cc, col in enumerate(c_colours):
        ax1.plot(np.NaN, np.NaN, c=c_colours[cc], label=colour_names[cc])

    ax2 = ax1.twinx()
    style_names = ['Averaged Most Connected Nodes', 'Averaged Mean Node connections',
                   'Averaged Least Connected Nodes']
    for ss, sty in enumerate(styles):
        ax2.plot(np.NaN, np.NaN, ls=styles[ss],
                 label=style_names[ss], c='black')
    ax2.get_yaxis().set_visible(False)

    plt.title('Average # of connections vs time')
    ax1.legend(loc=1)
    ax2.legend(loc=5)
    ax1.grid(True)
    plt.show()
    return


def plot_CombinedPickle_files(plot_paths, connection_paths):

    plt.figure(figsize=(8, 6))
    avg_time = np.array([])

    avg_most = np.array([])
    avg_uq = np.array([])
    avg_mean = np.array([])
    avg_lq = np.array([])
    avg_least = np.array([])

    avg_x = np.array([])
    avg_y = np.array([])
    func = None
    first = True
    for p_path, c_path in zip(plot_paths, connection_paths):
        plot_x, plot_y = read_pickle_file(p_path, 0)
        num_of_variants = plot_y.shape[0]//2
        plot_y = plot_y[:num_of_variants]
        connection_x, connection_y = read_pickle_file(c_path, 1)

        def interpolate_connections(old_x, new_x, new_y, axis):
            func = interp1d(new_x, new_y, axis=axis, fill_value='extrapolate')
            interpolated_y = func(old_x)
            return interpolated_y

        if first is True:
            avg_time = np.array(connection_x)
            avg_most, avg_uq, avg_mean, avg_lq, avg_least = np.array(connection_y)

            avg_x = plot_x
            avg_y = plot_y
            first = False
        else:
            for a, avg in enumerate([avg_most, avg_uq, avg_mean, avg_lq, avg_least]):
                avg += interpolate_connections(avg_time, connection_x, connection_y[a], 0)
            avg_y += interpolate_connections(avg_x, plot_x, plot_y, -1)
    p_colours = ['red', 'purple']
    c_colours = ['tab:blue', 'tab:orange']

    length = len(plot_paths)
    avg_most = avg_most / length
    avg_uq = avg_uq / length
    avg_mean = avg_mean / length
    avg_lq = avg_lq / length
    avg_least = avg_least / length
    avg_y = avg_y / length

    fig, ax1 = plt.subplots()
    ax2 = ax1.twinx()
    ax1.set_xlabel('Time')
    ax1.set_ylabel('Number of Connections', color='peru')
    legend_plots = []

    for var in range(2):
        a = ax1.plot(avg_time, avg_most[:, var], color=adjust_lightness(c_colours[var], 1),
                     label=f'Average number of Connections for variant {chr(ord("A") + var)}')
        b = ax1.plot(avg_time, avg_mean[:, var], color=adjust_lightness(c_colours[var], 1))
        c = ax1.plot(avg_time, avg_least[:, var], color=adjust_lightness(c_colours[var], 1))
        legend_plots.extend([a, b, c])


    ax2.set_ylabel('Number of Infected', color='cornflowerblue')

    a = ax2.plot(avg_x, avg_y[0], color=adjust_lightness(c_colours[0], 1.7),
                 label=f"Average number of nodes infected by A")
    b = ax2.plot(avg_x, avg_y[1], color=adjust_lightness(c_colours[1], 1.7),
                 label=f"Average number of nodes infected by B")
    legend_plots.extend([a, b])

    lines2, labels2 = ax2.get_legend_handles_labels()
    ax2.legend(lines2, labels2, loc='center right')

    ax2.annotate('Averaged Maximum Degree for A', xy=(0.35, 0.91), xycoords='axes fraction', color='tab:orange')
    ax2.annotate('Averaged Maximum Degree for B', xy=(0.35, 0.71), xycoords='axes fraction', color=adjust_lightness(c_colours[0], 0.8))
    ax2.annotate('Averaged Mean Degree for A', xy=(0.35, 0.195), xycoords='axes fraction', color='tab:orange')
    ax2.annotate('Averaged Mean Degree for B', xy=(0.35, 0.08), xycoords='axes fraction', color=adjust_lightness(c_colours[0], 0.8))

    ax1.grid(True)
    plt.show()

    return







if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Plot data from pickle files')
    parser.add_argument('file_paths', metavar='FILE', type=str, nargs='+',
                        help='path to pickle files')
    args = parser.parse_args()

    #plot_PlotPickle_files(args.file_paths)
    plot_ConnectionPickle_files(args.file_paths)