import os
import argparse
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats
from call_plot_loader import filter_numbers
from axes_contour import extract_data_from_csv


def parse_arguments():
    parser = argparse.ArgumentParser(description='Perform regression analysis on data files in a specified directory.')
    parser.add_argument('directory', type=str, help='Path to the directory containing data files.')

    parser.add_argument("mode", type=int, help="1 - Rel.vs.R0, 2 - Time.vs.Rel, 3 - R0.vs.Time")
    parser.add_argument("fixedAxis", type=int, help="Which axis in the mode are we keeping fixed; 0 or 1")
    parser.add_argument("fixedPoint", type=int, help="Which point (within fixed axis) are we examining")
    parser.add_argument("anomaliesFile", type=str, help="Path directory for the anomalies text file")
    return parser.parse_args()


# plt.plot(second_digit, intercept2 + slope2 * np.array(second_digit), color='blue', label='Second Peak Regression')
def normalise(xs, y_standard, coeffs):
    f = np.poly1d(coeffs)
    ys = f(xs)
    idys = np.where(((min(y_standard) < ys) & (ys < max(y_standard))))
    xs = xs[idys]
    ys = ys[idys]
    return xs, ys, f


def calculate_rmse(func, points):
    # Extract x and y values from the points
    x_values = np.array([point[0] for point in points])
    y_values = np.array([point[1] for point in points])

    # Calculate predicted y values using the function
    predicted_y_values = np.array([func(xs) for xs in x_values])
    # Calculate residuals
    residuals = y_values - predicted_y_values
    # Calculate mean squared error
    mse = np.mean(residuals ** 2)
    rmse = np.sqrt(mse)
    return rmse


if __name__ == '__main__':

    args = parse_arguments()

    if not os.path.exists(args.directory):
        print("Directory does not exist")
        exit()

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
    axis_resolution = 5
    file_repeats = 25

    infective_range = np.linspace(5, 15, axis_resolution)
    time_delay_range = np.linspace(0, 15, axis_resolution)
    relationship_range = np.linspace(0, 1, axis_resolution)
    axes = (time_delay_range, relationship_range, infective_range)
    names = ("Time Delay", "Relationship Coefficient", "R0 Value")
    z_axis_names = ("Peak time (since start of infection) for A (unit)",
                    "Peak time (since start of infection) for B (unit)",
                    "Disease A peak infection number",
                    "Disease B peak infection number",
                    "Time (since start of infection) when A hit 25% of population (unit)",
                    "Time (since start of infection) when B hit 25% of population (unit)",
                    "Time (since start of infection) when A hit 50% of population (unit)",
                    "Time (since start of infection) when B hit 50% of population (unit)",
                    "Number of nodes only infected by A",
                    "Number of nodes only infected by B",
                    "Number of nodes recovered from A",
                    "Number of nodes recovered from B",
                    "Total cumulative infections of A",
                    "Total cumulative infections of B")
    # Parse command-line arguments

    data_files = [file for file in os.listdir(args.directory) if file.startswith("dataFile") and file.endswith(".csv")]

    if args.fixedAxis == 0:
        safe_files = [filter_numbers(args.fixedPoint, k, args.anomaliesFile, resolution=file_repeats) for k in
                      range(axis_resolution)]
    elif args.fixedAxis == 1:
        safe_files = [filter_numbers(j, args.fixedPoint, args.anomaliesFile, resolution=file_repeats) for j in
                      range(axis_resolution)]
    else:
        raise ValueError

    x_axis = axes[(args.mode + int(not bool(int(args.fixedAxis)))) % 3]
    x_axis = 15 / x_axis
    y_axis = []
    for jk, file_sets in enumerate(safe_files):
        temp_y = []

        if args.fixedAxis == 0:
            j = args.fixedPoint
            k = jk
        else:
            j = jk
            k = args.fixedPoint

        for r, repeats in enumerate(file_sets):
            filename = f'dataFile_{repeats}.{j}.{k}.csv'
            filepath = os.path.join(args.directory, filename)

            if filename not in data_files:
                continue
            temp_y.append(extract_data_from_csv(filepath, z_axis_value))
        y_axis.append(temp_y)

    x_reg = []
    y_reg = []
    y_mean_points = []
    first = True
    repeats_lower_bound = min(list(map(len, y_axis)))
    for x, y_values in zip(x_axis, y_axis):
        x_reg.extend([x] * len(y_values))
        y_reg.extend(y_values)

        y_mean = sum(y_values) / len(y_values)
        y_mean_points.append(y_mean)
        y_plus = max(y_values) - y_mean
        y_minus = y_mean - min(y_values)
        asymmetric_error = np.array([[y_minus], [y_plus]])
        if first:
            plt.errorbar(x, y_mean, yerr=asymmetric_error, fmt='o', capsize=5, mfc='red', ecolor='black', mec='red',
                         label=f'Averaged y-value over {repeats_lower_bound} repeats')
            first = False
        else:
            plt.errorbar(x, y_mean, yerr=asymmetric_error, fmt='o', capsize=5, mfc='red', ecolor='black', mec='red')

        plt.xlabel(names[(args.mode + int(not bool(int(args.fixedAxis)))) % 3])
        plt.ylabel(z_axis_names[z_axis_value])
        plt.grid(True)

    degree = 2

    coefficients = np.polyfit(x_reg, y_reg, degree)
    print(coefficients)
    # Create the polynomial function
    polynomial = np.poly1d(coefficients)

    x_values = np.linspace(min(x_axis), max(x_axis), 100)
    y_values = polynomial(x_values)
    idys = np.where(((min(y_reg) < y_values) & (y_values < max(y_reg))))
    x_values_plot = x_values[idys]
    y_values_plot = y_values[idys]
    # Perform linear regression
    # slope1, intercept1, r_value1, p_value1, std_err1 = stats.linregress(x_reg, y_reg)
    # slope2, intercept2, r_value2, p_value2, std_err2 = stats.linregress(second_digit, second_peak_times)

    # Plot regression lines
    # plt.plot(x_values_plot, y_values_plot, color='green', label='r_{A,B} = 0.5')
    plt.plot(x_values_plot, y_values_plot, color='green', label=str(polynomial))

    print("BA RMSE:")
    print(calculate_rmse(polynomial, list(zip(x_axis, y_mean_points))))

    # coef1 = [  -57.33598945,   694.77176826, -3309.7876702 , 7740.83690261, 646.298354  ]
    # coef5 = [  -553.91978891,   5090.30065674, -17965.54255285,  30192.31975994, -13383.7954425 ]
    # coeff9 = [  462.33289747, -2773.11639881  ,3127.84919905  ,8261.01880247,-7743.9243068 ]
    # x_values5, er_ys5, f5 = normalise(x_values, y_reg, coef5)
    # x_values1, er_ys1, f1 = normalise(x_values, y_reg, coef1)
    # x_values9, er_ys9, f9 = normalise(x_values, y_reg, coeff9)

    # plt.plot(x_values1, er_ys1, label='r_{A,B} = 0.1')
    # plt.plot(x_values5, er_ys5, label='r_{A,B} = 0.5 for ER')
    # plt.plot(x_values9, er_ys9, label='r_{A,B} = 0.9')

    # print("ER RMSE:")
    # print(calculate_rmse(f5, list(zip(x_axis, y_mean_points))))

    plt.legend()
    plt.title("Regression Plot of Key Output Value vs changing input parameter ")
    plt.show()

    # print("First Peak Regression: slope =", slope1, ", intercept =", intercept1)
    # print("Second Peak Regression: slope =", slope2, ", intercept =", intercept2)
