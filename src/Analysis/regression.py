import os
import argparse
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats


# Function to parse command-line arguments
def parse_arguments():
    parser = argparse.ArgumentParser(description='Perform regression analysis on data files in a specified directory.')
    parser.add_argument('directory', type=str, help='Path to the directory containing data files.')
    return parser.parse_args()


# Parse command-line arguments
args = parse_arguments()

# Find relevant files in the specified directory
file_paths = [os.path.join(args.directory, file) for file in os.listdir(args.directory) if
              file.startswith('dataFile_') and file.endswith('.csv')]

if not file_paths:
    print("No matching CSV files found in the specified directory.")
    exit()

# Lists to store extracted data
second_digit = []
first_peak_times = []
second_peak_times = []

# Process each file
for file_path in file_paths:
    # Extract second digit from filename
    second_digit.append(int(file_path.split('_')[1][1:]))

    # Read CSV file
    df = pd.read_csv(file_path)

    # Extract peak times
    peak_times = df.iloc[0]['Peak_Times'].split(',')
    first_peak_times.append(float(peak_times[0]))
    second_peak_times.append(float(peak_times[1]))

# Plot the data
plt.scatter(second_digit, first_peak_times, label='First Peak Time')
plt.scatter(second_digit, second_peak_times, label='Second Peak Time')
plt.xlabel('Second Digit')
plt.ylabel('Peak Time')
plt.title('Peak Times vs. Second Digit')
plt.legend()
plt.grid(True)

# Perform linear regression
slope1, intercept1, r_value1, p_value1, std_err1 = stats.linregress(second_digit, first_peak_times)
slope2, intercept2, r_value2, p_value2, std_err2 = stats.linregress(second_digit, second_peak_times)

# Plot regression lines
plt.plot(second_digit, intercept1 + slope1 * np.array(second_digit), color='red', label='First Peak Regression')
plt.plot(second_digit, intercept2 + slope2 * np.array(second_digit), color='blue', label='Second Peak Regression')

plt.legend()

plt.show()

print("First Peak Regression: slope =", slope1, ", intercept =", intercept1)
print("Second Peak Regression: slope =", slope2, ", intercept =", intercept2)