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


if __name__ == '__main__':

    if len(sys.argv) != 4:
        print(f"Usage: {sys.argv[0]} <directory_of_data> <mode> <anomalies_file>")
        sys.exit(1)

    directory = sys.argv[1]
    mode = int(sys.argv[2])
    anomalies = sys.argv[3]

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

    z_axis_value = 3


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
    x_axis = axes[mode % 3]
    y_axis = axes[(mode + 1) % 3]

    if not os.path.exists(directory):
        print(f"{directory} does not exist")

    x = []
    y = []
    z = []
    data_files = [file for file in os.listdir(directory) if file.startswith("dataFile") and file.endswith(".csv")]
    safe_files = [filter_numbers(j, k, anomalies) for j in range(resolution) for k in range(resolution)]

    for f, file_sets in enumerate(safe_files):
        j = f // resolution
        k = f % resolution

        x_val, y_val = j, k
        x_val = x_axis[x_val]
        y_val = y_axis[y_val]

        sample_size = len(file_sets)
        z_val = 0

        for r, repeats in enumerate(file_sets):
            filename = f'dataFile_{repeats}.{j}.{k}.csv'
            filepath = os.path.join(directory, filename)

            z_val += extract_data_from_csv(filepath, z_axis_value)

        z_val = z_val / sample_size
        x.append(x_val)
        y.append(y_val)
        z.append(z_val)


    x = np.array(x)
    y = np.array(y)
    z = np.array(z)


    # Changing transmission time to R0
    if mode % 3 == 2:
        x = 15 / x
    elif (mode + 1) % 3 == 2:
        y = 15 / y

    plt.tricontourf(x, y, z, levels=20, cmap='viridis')
    plt.colorbar(label=z_axis_names[z_axis_value])
    plt.xlabel(names[mode % 3])
    plt.ylabel(names[(mode + 1) % 3])
    plt.title(titles[mode % 3])
    plt.grid()
    plt.show()
