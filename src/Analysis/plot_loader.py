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


def read_pickle_file(file_path):
    with open(file_path, 'rb') as f:
        data = pickle.load(f)
    return data['x'], data['y']


def plot_pickle_files(file_paths):
    plt.figure(figsize=(8, 6))
    mean_x = np.array([])
    mean_y = np.array([])
    func = None
    base_colors = ['tab:red', 'tab:orange', 'tab:blue', 'tab:green']
    for f, file_path in enumerate(file_paths):
        plot_x, plot_y = read_pickle_file(file_path)
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


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Plot data from pickle files')
    parser.add_argument('file_paths', metavar='FILE', type=str, nargs='+',
                        help='path to pickle files')
    args = parser.parse_args()

    plot_pickle_files(args.file_paths)
