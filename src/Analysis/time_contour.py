import pickle
import sys
import call_plot_loader
import numpy as np
import scipy.sparse as sp
import matplotlib.pyplot as plt

def read_file(file_path, ext_type):
    # ext_type [0] - pickle file
    # ext_type [1] - npz file
    if ext_type == 0:
        with open(file_path, 'rb') as f:
            data = pickle.load(f)
        return data['time'], data['nodes']
    elif ext_type == 1:
        data = sp.load_npz(file_path)
        return data
    else:
        raise ValueError


def find_closest_index(arr, x):
    # Find the index where x would be inserted to maintain sorted order
    idx = np.searchsorted(arr, x)

    # Handle cases where x is outside the range of arr
    if idx == 0:
        return 0
    elif idx == len(arr):
        return len(arr) - 1

    # Determine the index of the closest value
    left_diff = abs(arr[idx - 1] - x)
    right_diff = abs(arr[idx] - x)

    # Return the index with the minimum difference
    return idx - 1 if left_diff < right_diff else idx


def delete_row_col(matrix, x):
    #    if not sp.issparse(matrix):
    #        raise ValueError("Input matrix must be sparse.")

    #    if x < 0 or x >= matrix.shape[0]:
    #        raise ValueError("Invalid row/column index.")

    matrix[x, :] = 0  # Set all elements in row x to zero
    matrix[:, x] = 0  # Set all elements in column x to zero

    return matrix


def degree_distribution(imm_files, adj_files, t0, t1, var):
    avg_density_before_A = np.array([])
    avg_density_after_A = np.array([])

    avg_density_before_B = np.array([])
    avg_density_after_B = np.array([])
    average_list = [[avg_density_before_A, avg_density_after_A], [avg_density_before_B, avg_density_after_B]]

    for imm, adj in zip(imm_files, adj_files):
        time_arr, nodes = read_file(imm, 0)

        before_idx = find_closest_index(time_arr, t0)
        after_idx = find_closest_index(time_arr, t1)

        for variant in range(var):
            before_adjacency = read_file(adj, 1).toarray()
            after_adjacency = before_adjacency.copy()

            before_nodes = nodes[before_idx][variant]
            after_nodes = nodes[after_idx][variant]

            before_adjacency = delete_row_col(before_adjacency, before_nodes)

            after_adjacency = delete_row_col(after_adjacency, after_nodes)

            density_before = np.count_nonzero(before_adjacency, axis=1)
            density_after = np.count_nonzero(after_adjacency, axis=1)

            densities = [np.bincount(density_before[density_before > 0]), np.bincount(density_after[density_after > 0])]
            averages = average_list[variant]
            for i in range(len(densities)):
                max_length = max(len(averages[i]), len(densities[i]))

                padded_array1 = np.pad(averages[i], (0, max_length - len(averages[i])), 'constant')
                padded_array2 = np.pad(densities[i], (0, max_length - len(densities[i])), 'constant')

                averages[i] = padded_array1 + padded_array2
            average_list[variant][0] = averages[0]
            average_list[variant][1] = averages[1]

    num_of_files = len(imm_files)
    normalise = lambda avg: avg / num_of_files

    average_list[0] = list(map(normalise, average_list[0]))
    average_list[1] = list(map(normalise, average_list[1]))

    max_length = max(len(array) for sublist in average_list for array in sublist)
    average_list = [[np.pad(array, (0, max_length - len(array)), 'constant') for array in sublist]
                    for sublist in average_list]

    fig, ax = plt.subplots(1, 2, figsize=(10, 6), gridspec_kw={'wspace': 0.5})

    base_colours = ['blue', 'orange']
    styles = ['-.', '-']

    for i in range(2): #range(len(average_list)):
        #ax1.plot(np.arange(len(average_list[i][0])), average_list[i][0], color=base_colours[i], ls=styles[0])
        #ax1.plot(np.arange(len(average_list[i][1])), average_list[i][1], color=base_colours[i], ls=styles[1])
        bins = range(len(average_list[i][0]) + 1)
        ax[i].hist(range(len(average_list[i][0])), bins=bins, weights=average_list[i][0], color=base_colours[0],
                 histtype='barstacked', rwidth=0.8)
        bins = range(len(average_list[i][1]) + 1)
        ax[i].hist(range(len(average_list[i][1])), bins=bins, weights=average_list[i][1], color=base_colours[1],
                 histtype='barstacked', rwidth=0.4)

    color_names = [f't={t0}', f't={t1}']

    for cc, col in enumerate(base_colours):
        ax[0].plot(np.NaN, np.NaN, c=base_colours[cc], label=color_names[cc])
        ax[1].plot(np.NaN, np.NaN, c=base_colours[cc], label=color_names[cc])


    for i in range(2):
        ax[i].set_xlabel('Degree number')
        ax[i].set_ylabel('Number of nodes')
        ax[i].legend(loc=1)

        ax[i].grid(True)
        ax[i].set_title(f'Variant {chr(ord("A") + i)}')
    fig.suptitle("Degree Distribution of nodes at two key times")
    plt.show()
    return 2


def lcc(imm_files, adj_files):
    return


def degree_assortativity(imm_files, adj_files):
    return


if __name__ == '__main__':
    if len(sys.argv) != 4:
        print(f"Usage: python {sys.argv[0]} <input_file> <anomalies_file> <regime>")
        sys.exit(1)

    input_file = sys.argv[1]
    anom_file = sys.argv[2]
    regime = int(sys.argv[3])

    immune_files = call_plot_loader.generate_file_paths(input_file, anom_file)
    adjacency_files = [i.split('immuneData')[0] + 'networkAdjacencySparse' +
                       i.split('immuneData')[1].split('pickle')[0] + 'npz'
                       for i in immune_files]
    if regime == 0:
        before = 15
        after = 20
        number_of_variants = 2
        degree_distribution(immune_files, adjacency_files, before, after, number_of_variants)
    elif regime == 1:
        lcc(immune_files, adjacency_files)
    elif regime == 2:
        degree_assortativity(immune_files, adjacency_files)
    else:
        print("Regime OutOfBounds")
        sys.exit(1)
