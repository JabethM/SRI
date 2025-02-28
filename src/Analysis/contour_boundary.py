import os
import matplotlib.pyplot as plt
from sklearn.linear_model import LinearRegression
import numpy as np
from axes_contour import format_data
import sys
from scipy.interpolate import griddata, interp1d

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

    z_axis_value = 5

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

    if data_mode % 3 == 2:
        x = 15 / x
    elif (data_mode + 1) % 3 == 2:
        y = 15 / y

    plt.tricontourf(x, y, z, levels=20, cmap='viridis')
    plt.colorbar(label=z_axis_names[z_axis_value])


    z_midpoint = np.max(z) / 2

    X = x.reshape(10, 10)
    Y = x.reshape(10, 10)

    interp_y = griddata((x, z), y, (z_midpoint, z_midpoint), method='linear')
    interp_x = griddata((y, z), x, (interp_y, z_midpoint), method='linear')


    plt.xlabel(names[data_mode % 3])
    plt.ylabel(names[(data_mode + 1) % 3])
    plt.title(titles[data_mode % 3])

    plt.plot(interp_x, interp_y, color='red', linestyle='--', label='Fitted Line')
    plt.grid()
    plt.show()

