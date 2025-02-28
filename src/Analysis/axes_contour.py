import scipy.sparse as sp
import numpy as np
import matplotlib.pyplot as plt
import os
from argparse import ArgumentParser
import sys
from call_plot_loader import filter_numbers


def extract_data_from_csv(csv_file, z_axis):
    """
    # z axis can be several things:
    # Peak time - Start time for A[0]
    # Peak time - Start time for B [1]

    # Disease A peak infection number [2]
    # Disease B peak infection number [3]

    # Disease A hit 25% [4]
    # Disease B hit 25% [5]

    # Disease A hit  50% [6]
    # Disease B hit  50% [7]

    # Unique Infections of A [8]
    # Unique Infections of B [9]

    # Unique Recoveries of A [10]
    # Unique Recoveries of B [11]

    # Cumulative Infections of A [12]
    # Cumulative Infections of B [13]

    :param csv_file:
    :param z_axis:
    :return:
    """
    with open(csv_file, 'r') as f:
        lines = f.readlines()
        start_time = (lines[1].strip().split(',')[1:])
        lowerQ_infection_time = (lines[2].strip().split(',')[1:])
        median_infection_time = (lines[3].strip().split(',')[1:])

        peak_times = (lines[4].strip().split(',')[1:])
        peak_infections = lines[5].strip().split(',')[1:]
        cumulative_infections = lines[6].strip().split(',')[1:]
        total_unique_infections = lines[7].strip().split(',')[1:]
        total_unique_recoveries = (lines[8].strip().split(',')[1:])

        strlist_to_float = lambda a, d: [float(a.replace('"', '')), float(d.replace('"', ''))]

        if z_axis < 2:
            st = np.array(strlist_to_float(*start_time))
            pt = np.array(strlist_to_float(*peak_times))
            key_stat = pt - st
        elif z_axis < 4:
            key_stat = np.array(strlist_to_float(*peak_infections))
        elif z_axis < 6:
            st = np.array(strlist_to_float(*start_time))
            key_stat = np.array(strlist_to_float(*lowerQ_infection_time)) - st
        elif z_axis < 8:
            st = np.array(strlist_to_float(*start_time))
            key_stat = np.array(strlist_to_float(*median_infection_time)) - st
        elif z_axis < 10:
            key_stat = np.array(strlist_to_float(*total_unique_infections))
        elif z_axis < 12:
            key_stat = np.array(strlist_to_float(*total_unique_recoveries))
        elif z_axis < 14:
            key_stat = np.array(strlist_to_float(*cumulative_infections))

        z_axis = z_axis % 2

        return key_stat[z_axis]


def format_data(directory, anomalies, resolution, x_axis, y_axis, z_axis_value, doubles=False, directory_addition=None, anomalies_addition=None):

    x = []
    y = []
    z = []
    data_files = [file for file in os.listdir(directory) if file.startswith("dataFile") and file.endswith(".csv")]
    safe_files = [filter_numbers(j, k, anomalies) for j in range(resolution) for k in range(resolution)]
    if doubles:
        df_2 = [file for file in os.listdir(directory_addition) if file.startswith("dataFile")
                and file.endswith(".csv")]
        sf_2 = [[-repeat - 1 for repeat in filter_numbers(j, k, anomalies_addition)] for j in range(resolution) for k in
                range(resolution)]
        safe_files = [a + b for a, b in zip(safe_files, sf_2)]
    for f, file_sets in enumerate(safe_files):
        j = f // resolution
        k = f % resolution

        x_val, y_val = j, k
        x_val = x_axis[x_val]
        y_val = y_axis[y_val]

        sample_size = len(file_sets)
        z_val = 0

        for r, repeats in enumerate(file_sets):
            if repeats < 0:
                repeats = abs(repeats + 1)
                loop_dir = directory_addition
            else:
                loop_dir = directory

            filename = f'dataFile_{repeats}.{j}.{k}.csv'

            filepath = os.path.join(loop_dir, filename)
            temp_z = extract_data_from_csv(filepath, z_axis_value)

            if ((4 <= z_axis_value <= 7) or (1 <= z_axis_value <= 2)) and temp_z < 0:
                temp_z = 100
            z_val += temp_z

        z_val = z_val / sample_size
        x.append(x_val)
        y.append(y_val)
        z.append(z_val)
    return np.array(x), np.array(y), np.array(z)


if __name__ == '__main__':

    doubles = False
    directory_addition = None
    anomalies_addition = None
    if len(sys.argv) == 4:
        directory = sys.argv[1]
        data_mode = int(sys.argv[2])
        anomalies = sys.argv[3]

    elif len(sys.argv) == 6:
        doubles = True
        directory = sys.argv[1]
        directory_addition = sys.argv[2]

        data_mode = int(sys.argv[3])

        anomalies = sys.argv[4]
        anomalies_addition = sys.argv[5]

    else:
        print(f"Usage: {sys.argv[0]} <directory_of_data> <data_mode> <anomalies_file>")
        sys.exit(1)

    if data_mode > 3 or data_mode < 1:
        raise ValueError

    resolution = 10

    # z axis can be several things:
    # Peak time - Start time for A[0]
    # Peak time - Start time for B [1]

    # Disease A peak infection number [2]
    # Disease B peak infection number [3]

    # Disease A hit 25% [4]
    # Disease B hit 25% [5]

    # Disease A hit  50% [6]
    # Disease B hit  50% [7]

    # Unique Infections of A [8]
    # Unique Infections of B [9]

    # Unique Recoveries of A [10]
    # Unique Recoveries of B [11]

    # Cumulative Infections of A [12]
    # Cumulative Infections of B [13]

    z_axis_value = 1

    infective_range = np.linspace(5, 15, resolution)
    time_delay_range = np.linspace(0, 15, resolution)
    relationship_range = np.linspace(0, 1, resolution)
    axes = (time_delay_range, relationship_range, infective_range)
    names = ("Time Delay", "Relationship Coefficient", "R0 Value")
    titles = ("Contour Plot of Time Delay vs Relationship Coefficient",
              "Contour Plot of Relationship Coefficient vs R0 Value",
              "Contour Plot of R0 Value vs Time Delay")
    z_axis_names = ("Peak time (since start of infection) for A",
                    "Peak time (since start of infection) for B",
                    "Disease A peak infection number",
                    "Disease B peak infection number",
                    "Time (since start of infection) when A hit 25% of population",
                    "Time (since start of infection) when B hit 25% of population",
                    "Time (since start of infection) when A hit 50% of population",
                    "Time (since start of infection) when B hit 50% of population",
                    "Number of nodes only infected by A",
                    "Number of nodes only infected by B",
                    "Number of nodes recovered from A",
                    "Number of nodes recovered from B",
                    "Total cumulative infections of A",
                    "Total cumulative infections of B")
    x_axis = axes[data_mode % 3]
    y_axis = axes[(data_mode + 1) % 3]

    if not os.path.exists(directory):
        print(f"{directory} does not exist")

    x, y, z = format_data(directory, anomalies, resolution, x_axis, y_axis, z_axis_value, doubles=doubles,
                          directory_addition=directory_addition, anomalies_addition=anomalies_addition)

    # Changing transmission time to R0
    if data_mode % 3 == 2:
        x = 15 / x
    elif (data_mode + 1) % 3 == 2:
        y = 15 / y

    plt.tricontourf(x, y, z, levels=20, cmap='viridis')
    plt.colorbar(label=z_axis_names[z_axis_value])
    plt.xlabel(names[data_mode % 3])
    plt.ylabel(names[(data_mode + 1) % 3])
    plt.title(titles[data_mode % 3])
    plt.grid()
    plt.show()
