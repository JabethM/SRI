import os
import matplotlib.pyplot as plt
import numpy as np
from axes_contour import format_data
import sys
from scipy.interpolate import RegularGridInterpolator
from sklearn.preprocessing import PolynomialFeatures
from sklearn.linear_model import LinearRegression
from regression import normalise

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
    inversion = False
    z_axis_value = 7
    level = 5000

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

    xs, ys, zs = format_data(directory, anomalies, resolution, x_axis, y_axis, z_axis_value, number_of_repeats=25, doubles=doubles,
                             directory_addition=directory_addition, anomalies_addition=anomalies_addition)

    x_interp_axis = x_axis
    y_interp_axis = y_axis
    if data_mode % 3 == 2:
        xs = 15 / xs
        x_interp_axis = 15 / x_axis
    elif (data_mode + 1) % 3 == 2:
        ys = 15 / ys
        y_interp_axis = 15 / y_axis

    if inversion:
        plt.tricontourf(ys, xs, zs, levels=20, cmap='viridis')
    else:
        plt.tricontourf(xs, ys, zs, levels=20, cmap='viridis', extend='both')

    plt.colorbar(label=z_axis_names[z_axis_value])
    z_grid = zs.reshape(10, 10)
    # Interpolate the function z=f(x, y)

    interp_func = RegularGridInterpolator((x_interp_axis, y_interp_axis), z_grid, method="cubic")

    # Define a function to find x given y for z=50
    def find_x_for_z_50(y):
        # Use scipy's optimize to find the x value where z is closest to 50 for a given y
        from scipy.optimize import root_scalar, minimize_scalar
        objective = lambda x: interp_func([[x, y]])[0] - level
        try:
            root = root_scalar(objective, bracket=[min(x_interp_axis), max(x_interp_axis)])
        except ValueError:
            return np.nan

        # res = minimize_scalar(objective, method='bounded',
        #                      bounds=(0.5, 1))  # Adjust bounds as per your x range

        if root.converged:
            return root.root
        else:
            return np.nan


    def chebyshev_nodes(a, b, n):
        k = np.arange(1, n + 1)
        x = np.cos((2 * k - 1) * np.pi / (2 * n))
        return 0.5 * (a + b) + 0.5 * (b - a) * x


    y_range = chebyshev_nodes(min(ys), max(ys), len(ys))
    # Find x values corresponding to z=50 for a range of y values
    #y_range = np.linspace(ys.min(), ys.max(), 100)
    removals = []
    x_range = []
    for idx, y in enumerate(y_range):
        a = find_x_for_z_50(y)
        if not np.isnan(a):
            x_range.append(find_x_for_z_50(y))
        else:
            removals.append(idx)

    x_range = np.array(x_range)

    # Polynomial regression
    #degree = 4  # Degree of the polynomial
    #y_range = np.delete(y_range, removals)
    #coefficients = np.polyfit(x_range, y_range, degree)
    #print(coefficients)
    # Create the polynomial function
    #polynomial = np.poly1d(coefficients)

    # Generate x values for plotting the curve
    x_values = np.linspace(min(xs), max(xs), 100)

    # Calculate corresponding y values using the polynomial function
    #y_values = polynomial(x_values)

    #idys = np.where(((min(ys) < y_values) & (y_values < max(ys))))
    #x_out = x_values[idys]
    #y_out = y_values[idys]

    M2Coeff = [-2.25736785,   30.28240281, -113.8008649,   179.09496867,  -94.4154052 ]
    #M3Coeff = [3.82221478e-05, -1.49879364e-03, 2.34403774e-02, -1.90909918e-01, 9.90811191e-01]
    M2_xs, M2_ys, f2 = normalise(x_values, ys, M2Coeff)
    #M3_xs, M3_ys, f3 = normalise(x_values, ys, M3Coeff)

    plt.plot(M2_xs, M2_ys, color='red', linestyle='--', label='Peak-Infection Generated Curve')
    #plt.plot(M3_xs, M3_ys, color='red', linestyle='--', label='Peak-Infection Generated Curve')
    plt.xlabel(names[data_mode % 3])
    plt.ylabel(names[(data_mode + 1) % 3])
    """if inversion:
        plt.plot(y_out, x_out, color='red', linestyle='--', label=str(polynomial).replace('x','y'))
        plt.ylabel(names[data_mode % 3])
        plt.xlabel(names[(data_mode + 1) % 3])
    else:
        plt.plot(x_out, y_out, color='red', linestyle='--', label=str(polynomial))
        plt.xlabel(names[data_mode % 3])
        plt.ylabel(names[(data_mode + 1) % 3])"""
    #print(polynomial)
    plt.title(titles[data_mode % 3])
    plt.legend()
    plt.grid()
    plt.show()
