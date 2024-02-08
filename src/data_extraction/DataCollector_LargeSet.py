from ..simulation.SIR_system import SIR
from ..simulation.SIR_system import Variant
from argparse import ArgumentParser
import numpy as np
import networkx as nx
import os
import json
import uuid
import shutil
import csv
import matplotlib.pyplot as plt
import scipy.sparse as sp

par = ArgumentParser()

current_directory = os.getcwd()
config_file = "config_sim.json"
simulation_folder = os.path.join('scripts', 'SIR_simulation_files')
sim_relative_path = os.path.relpath(simulation_folder, current_directory)
default_config_file_path = os.path.join(sim_relative_path, config_file)

par.add_argument("-f", "--config-file", type=str, default=default_config_file_path, help="Path for JSON config file")
args = par.parse_args()

config_file_path = args.config_file

print("Initalising...")
execute = SIR(config_file=config_file_path)
num_of_nodes = execute.num_nodes
number_of_variants = execute.num_of_variants
end_steps = 50000

# Set Up
# ------------------------------------------------------------------------
# Peak infected nodes; number and time
peak_infection = np.zeros((number_of_variants, 2))

# Total time of epidemic: Beginning time, End time, Boolean: True when beginning has started
total_infection_time = np.zeros((number_of_variants, 3))

# Rates of infection & recovery
variant_rate_properties = np.zeros((number_of_variants, 2))

# Network Adjacency Matrix
adj_matrix_sparse = nx.adjacency_matrix(execute.G)

# Relationship Matrix
relationship_matrix = np.zeros((number_of_variants, number_of_variants))

# Infected and Recovered populations at key points: Start of infection, Peak Infection
data_sets = [[[], []] for _ in range(number_of_variants)]

plot_y = np.zeros((number_of_variants, 1))
plot_x = [0]
time = np.array([0])
end = False

print("Running Sim...")
while not end:
    current_time = execute.step_run()
    previous_time = execute.old_time

    infected_sets = execute.infected_set
    recovered_sets = execute.recovered_set
    current_infection_num = np.array(list(map(len, infected_sets)))
    current_recovered_num = np.array(list(map(len, recovered_sets)))

    plot_y = np.hstack((plot_y, current_infection_num.reshape(-1, 1)))
    plot_x.append(current_time)

    # Calculate Peak Time
    # ------------------------------------------------------------------------
    peak_inf_mask = current_infection_num > peak_infection[:, 1]
    peak_infection[peak_inf_mask, 1] = current_infection_num[peak_inf_mask]
    peak_infection[:, 0] = np.where(peak_inf_mask, current_time, peak_infection[:, 0])

    # Note Beginning and End of an Epidemic
    # ------------------------------------------------------------------------
    # > list of diseases that haven't "started" (AND) list of diseases who have infected more than 0 people

    infection_start_condition = (total_infection_time[:, -1] == 0) & (current_infection_num > 0)

    # > Note down the current time, mark disease as having started (aka 1)
    total_infection_time[infection_start_condition, 0] = previous_time
    total_infection_time[infection_start_condition, -1] = 1

    # > list of diseases that have started (AND) list of diseases that have less than 1 infected nodes left
    infection_end_condition = (total_infection_time[:, -1] == 1) & (current_infection_num <= 0)
    total_infection_time[infection_end_condition, 1] = current_time
    total_infection_time[infection_end_condition, -1] = 2

    # Store infected and recovered population multi_variant_files
    # ------------------------------------------------------------------------
    for i in range(len(data_sets)):
        # > if a disease has just started, mark down the infected and recovered
        data_sets[i][0] = list(recovered_sets[i]) if infection_start_condition[i] \
            else data_sets[i][0]
        data_sets[i][1] = list(recovered_sets[i]) if peak_inf_mask[i] else data_sets[i][1]

    # Collect Infection and Recovery Rates of each disease
    # ------------------------------------------------------------------------
    rates_collected = False
    if (execute.variant_count >= number_of_variants) and not rates_collected:
        variant_rate_properties = execute.rates
        relationship_matrix = Variant.relation_matrix

        rates_collected = True

    # End conditions
    # ------------------------------------------------------------------------
    end = execute.end

    check = 0
    if rates_collected and (execute.steps >= end_steps):
        is_undying = np.where(variant_rate_properties[:, 0] > variant_rate_properties[:, 1])[0]
        carrier_majority = current_infection_num > 0.9 * (num_of_nodes - current_recovered_num)
        leader_of_the_pack = np.argmax(current_infection_num)
        end = end or (carrier_majority[leader_of_the_pack] and is_undying[leader_of_the_pack])

        check += 1
        if check >= 100:
            check = 0
            end_steps += 1000
print()
print("Data Dump...")
# Data Dump
# ------------------------------------------------------------------------
data_folder_name = str(uuid.uuid4())
data_output_path = os.path.join("run_data", data_folder_name)

if not os.path.exists(data_output_path):
    os.makedirs(data_output_path)

# > family tree dump (+ rates)
family_tree = execute.variants[0].to_json(execute.variants[0])
with open(os.path.join(data_output_path, "family_tree.json"), "w") as file:
    json.dump(family_tree, file, indent=4)


def np_data_dump(file_name, array, headers):
    with open(file_name, mode='w', newline='') as file:
        writer = csv.writer(file)
        writer.writerow(headers)
        writer.writerows(array)


# > peak infection time dump (+ rates of infection and recovery)
np_data_dump(os.path.join(data_output_path, 'peak_infections.csv'), peak_infection, ['Time', 'Number of Infections'])

# > infection Start and End time
np_data_dump(os.path.join(data_output_path, 'infection_duration.csv'), total_infection_time, ['Start Time', 'End Time', 'Status'])

# > variant rate properties
np_data_dump(os.path.join(data_output_path, 'rate_properties.csv'), variant_rate_properties, ['Rate of Infection per day',
                                                                                  'Rate of Recovery per dat'])

# > relationship matrix dump
np.savetxt(os.path.join(data_output_path, 'relationship_matrix.csv'), relationship_matrix, delimiter=',', fmt='%.5f')

# > network adjacency matrix dump
sp.save_npz(os.path.join(data_output_path, 'network_adjacency_sparse.npz'), adj_matrix_sparse)

# > plot against time dump
for v in range(number_of_variants):
    plt.plot(plot_x, plot_y[v], label=(chr(ord('A') + v)))
plt.xlabel('Time (days)')
plt.ylabel('Number of Infected')
plt.legend()
plt.savefig(os.path.join(data_output_path, 'Infections_vs_Time_plot.png'))


# > recovered multi_variant_files at key points (for network reconstruction)
with open(os.path.join(data_output_path, 'recovered_data.csv'), 'w', newline='') as csvfile:
    writer = csv.writer(csvfile)

    for header, recovery_data in enumerate(data_sets):
        writer.writerow([f"-{chr(ord('A') + header)}"])

        for j, rec in enumerate(recovery_data):
            if j == 0:
                row_header = 'Start'
            else:
                row_header = 'Peak'
            writer.writerow([f"{row_header}:", *rec])

shutil.copy2(config_file_path, data_output_path)
