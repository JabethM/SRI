import pickle
import sys
import call_plot_loader
import numpy as np
import scipy.sparse as sp
import matplotlib.pyplot as plt
import networkx as nx
import os
from scipy.interpolate import interp1d


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


def find_consecutive_zeros(arr):
    # Find indices where the array is equal to zero
    zero_indices = np.where(arr == 0)[0]

    # Initialize variables to track consecutive zeros
    start_index = None
    consecutive_zeros = 0
    zero_ranges = []

    # Iterate through the zero indices
    for index in zero_indices:
        if start_index is None:
            # If this is the first zero encountered, mark its index
            start_index = index
            consecutive_zeros = 1
        elif index == start_index + consecutive_zeros:
            # If this zero is consecutive to the previous one, increment the count
            consecutive_zeros += 1
        else:
            # If the zero is not consecutive, record the range of consecutive zeros
            zero_ranges.append((start_index, consecutive_zeros))
            # Reset start_index and consecutive_zeros for the new range
            start_index = index
            consecutive_zeros = 1

    # Append the last range of consecutive zeros if any
    if start_index is not None:
        zero_ranges.append((start_index, consecutive_zeros))

    return zero_ranges


def remove_ranges(arr, ranges):
    mask = np.zeros_like(arr, dtype=bool)
    for start, length in ranges:
        mask[start:start + length] = True
    return arr[~mask]


def add_consecutive_numbers(arr):
    # Create an array with the original numbers and numbers incremented by 1
    incremented = np.concatenate([arr, arr + 1])
    # Sort the array
    incremented.sort()
    # Remove duplicates
    result = np.unique(incremented)
    return result


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

    for i in range(2):  # range(len(average_list)):
        # ax1.plot(np.arange(len(average_list[i][0])), average_list[i][0], color=base_colours[i], ls=styles[0])
        # ax1.plot(np.arange(len(average_list[i][1])), average_list[i][1], color=base_colours[i], ls=styles[1])
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


def lcc(imm_files, adj_files, t0, t1, var):
    avg_density_before_A = np.array([])
    avg_density_after_A = np.array([])

    avg_density_before_B = np.array([])
    avg_density_after_B = np.array([])
    average_list = [[avg_density_before_A, avg_density_after_A], [avg_density_before_B, avg_density_after_B]]
    count = 0
    for imm, adj in zip(imm_files[:3], adj_files[:3]):
        #if count == 1:
        #    break
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

            graphs = [nx.from_numpy_array(before_adjacency), nx.from_numpy_array(after_adjacency)]

            density_before = np.sort([len(component) for component in nx.connected_components(graphs[0])])
            density_after = np.sort([len(component) for component in nx.connected_components(graphs[1])])

            # density_before = np.count_nonzero(before_adjacency, axis=1)
            # density_after = np.count_nonzero(after_adjacency, axis=1)

            densities = [np.bincount(density_before[density_before > 0]), np.bincount(density_after[density_after > 0])]

            averages = average_list[variant]
            for i in range(len(densities)):
                max_length = max(len(averages[i]), len(densities[i]))

                padded_array1 = np.pad(averages[i], (0, max_length - len(averages[i])), 'constant')
                padded_array2 = np.pad(densities[i], (0, max_length - len(densities[i])), 'constant')

                averages[i] = padded_array1 + padded_array2
            average_list[variant][0] = averages[0]
            average_list[variant][1] = averages[1]
        count += 1

    num_of_files = len(imm_files)
    normalise = lambda avg: avg / num_of_files

    average_list[0] = list(map(normalise, list(map(np.array, average_list[0]))))
    average_list[1] = list(map(normalise, list(map(np.array, average_list[1]))))

    max_length = max(len(array) for sublist in average_list for array in sublist)
    average_list = [[np.pad(array, (0, max_length - len(array)), 'constant') for array in sublist]
                    for sublist in average_list]

    fig, ax = plt.subplots(1, 2, figsize=(10, 6), gridspec_kw={'wspace': 0.5})

    base_colours = ['blue', 'orange']
    styles = ['-.', '-']

    for i in range(2):  # range(len(average_list)):

        bins = range(len(average_list[i][0]) + 1)
        ax[i].stairs(average_list[i][0], bins, linewidth=5)

        bins = range(len(average_list[i][1]) + 1)
        ax[i].stairs(average_list[i][1], bins, linewidth=3)


    color_names = [f't={t0}', f't={t1}']

    for cc, col in enumerate(base_colours):
        ax[0].plot(np.NaN, np.NaN, c=base_colours[cc], label=color_names[cc])
        ax[1].plot(np.NaN, np.NaN, c=base_colours[cc], label=color_names[cc])

    for i in range(2):
        ax[i].set_xlabel('Size of components')
        ax[i].set_ylabel('Number of components (Averaged over multiple trials)')
        ax[i].legend(loc=1)
        ax[i].set_yscale('log')
        #ax[i].set_xscale('log')
        ax[i].grid(True)
        ax[i].set_title(f'Variant {chr(ord("A") + i)}')
    fig.suptitle("Averaged size distribution of Network Components at key times")
    plt.show()
    return


def assortativity_helper(adj, nodes, loop):
    assort = [np.array([]), np.array([])]
    for t in range(loop):
        full_t_adjacency = read_file(adj, 1).toarray()
        removed_nodes = nodes[t][0]  # 0 being variant A
        reduced_t_adjacency_A = delete_row_col(full_t_adjacency, removed_nodes)

        full_t_adjacency = read_file(adj, 1).toarray()
        removed_nodes = nodes[t][1]  # 1 being variant B
        reduced_t_adjacency_B = delete_row_col(full_t_adjacency, removed_nodes)

        G_A = nx.from_numpy_array(reduced_t_adjacency_A)
        G_B = nx.from_numpy_array(reduced_t_adjacency_B)

        d_A = nx.degree_assortativity_coefficient(G_A)
        d_B = nx.degree_assortativity_coefficient(G_B)

        assort[0] = np.append(assort[0], d_A)
        assort[1] = np.append(assort[1], d_B)
    return assort


def degree_assortativity(imm_files, adj_files):
    avg_time = np.array([])
    all_assort = [[np.array([0]), np.array([0])]]
    count = 0
    for imm, adj in zip(imm_files, adj_files):
        time_arr, nodes = read_file(imm, 0)

        if count == 0:
            avg_time = time_arr
            all_assort[0][0], all_assort[0][1] = assortativity_helper(adj, nodes, len(time_arr))
        else:
            temp_assortativity = assortativity_helper(adj, nodes, len(time_arr))

            func = interp1d(time_arr, temp_assortativity, axis=-1, fill_value='extrapolate')
            interpolated_y = func(avg_time)

            all_assort.extend([interpolated_y])

        count += 1
    all_assort = np.array(list(map(np.array, all_assort)))
    num_of_files = count
    normalise = lambda avg: avg / num_of_files

    avg_assort = np.sum(np.nan_to_num(all_assort), axis=0)
    avg_assort = normalise(avg_assort)

    max_assort = np.max(np.nan_to_num(all_assort), axis=0)
    min_assort = np.min(np.nan_to_num(all_assort), axis=0)


    base_colours = ['blue', 'orange']
    color_names = [f'Variant A', f'Variant B']

    styles = ['-.', '-']

    for var in range(2):  # range(len(average_list)):
        for i in range(len(avg_time)):
            a_mean = avg_assort[var][i]
            a_plus = max_assort[var][i] - a_mean
            a_minus = a_mean - min_assort[var][i]
            asymmetric_error = np.array([[a_minus], [a_plus]])

            plt.errorbar(avg_time[i], a_mean, yerr=asymmetric_error, fmt='o', capsize=5,
                         mfc=base_colours[var], ecolor=base_colours[var], mec=base_colours[var])

    for cc, col in enumerate(base_colours):
        plt.plot(np.NaN, np.NaN, c=base_colours[cc], label=color_names[cc])

    plt.xlabel('Time')
    plt.ylabel('Degree Assortativity of the network (as computed by nx)')
    plt.legend(loc=2)
    plt.grid(True)
    plt.title("Network Assortativity.vs.Time")

    plt.show()


if __name__ == '__main__':
    if len(sys.argv) != 4:
        print(f"Usage: python {sys.argv[0]} <input_file> <anomalies_file> <regime>")
        sys.exit(1)

    input_file = sys.argv[1]
    anom_file = sys.argv[2]
    regime = int(sys.argv[3])

    immune_files = call_plot_loader.generate_file_paths(input_file, anom_file, resolution=25)
    adjacency_files = [i.split('immuneData')[0] + 'networkAdjacencySparse' +
                       i.split('immuneData')[1].split('pickle')[0] + 'npz'
                       for i in immune_files]
    number_of_variants = 2
    before = 15
    after = 20

    if regime == 0:
        degree_distribution(immune_files, adjacency_files, before, after, number_of_variants)
    elif regime == 1:
        lcc(immune_files, adjacency_files, before, after, number_of_variants)
    elif regime == 2:

        degree_assortativity(immune_files, adjacency_files)
    else:
        print("Regime OutOfBounds")
        sys.exit(1)
