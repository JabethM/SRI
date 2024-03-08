import scipy.sparse as sp
import numpy as np
import matplotlib.pyplot as plt
import os
from argparse import ArgumentParser
import sys
from call_plot_loader import filter_numbers


def extract_data_from_csv(csv_file):
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

        pt = strlist_to_float(*peak_times)
        pi = strlist_to_float(*peak_infections)
        ci = strlist_to_float(*cumulative_infections)
        tui = strlist_to_float(*total_unique_infections)

        return pt, pi, ci, tui

if len(sys.argv) != 3:
    print(f"Usage: {sys.argv[0]} <directory of data> <mode>")
    sys.exit(1)


directory = sys.argv[1]
mode = int(sys.argv[2])

infective_range = np.linspace(5, 15, 10)
time_delay_range = np.linspace(0, 15, 10)
relationship_range = np.linspace(0, 1, 10)
axes = (time_delay_range, relationship_range, infective_range)
names = ("Time Delay", "Relationship Coefficient", "R0 Value")
titles = ("Contour Plot of Time Delay vs Relationship Coefficient",
          "Contour Plot of Relationship Coefficient vs R0 Value",
          "Contour Plot of R0 Value vs Time Delay")
x_axis = axes[mode % 3]
y_axis = axes[(mode + 1) % 3]

if not os.path.exists(directory):
    print(f"{directory} does not exist")

x = []
y = []
z_ancestor = []
z_descendent = []

for filename in os.listdir(directory):
    if filename.startswith("datafile") and filename.endswith(".csv"):
        filepath = os.path.join(directory, filename)

        x_val, y_val = map(int, filename.split('.')[1:3])
        x_val = x_axis[x_val]
        y_val = y_axis[y_val]

        peak_inf_vals = extract_data_from_csv(filepath)[3]
        x.append(x_val)
        y.append(y_val)
        z_ancestor.append(peak_inf_vals[0])
        z_descendent.append(peak_inf_vals[1])

x = np.array(x)
y = np.array(y)
z_ancestor = np.array(z_ancestor)
z_descendent = np.array(z_descendent)

if mode % 3 == 2:
    x = 15 / x
elif (mode + 1) % 3 == 2:
    y = 15/ y

plt.tricontourf(x, y, z_descendent, levels=20, cmap='viridis')
plt.colorbar(label='Unique Infections')
plt.xlabel(names[mode % 3])
plt.ylabel(names[(mode + 1) % 3])
plt.title(titles[mode % 3])
plt.grid()
plt.show()