from ..simulation.SIR_system import SIR
import numpy as np
import networkx as nx
import scipy.sparse as sp
import os
import csv
import matplotlib.pyplot as plt
import sys

if len(sys.argv) != 3:
    print("Usage: DataCollector_LargeSet.py <configuration_file> <output_folder>")
    sys.exit(1)

config_file_path = sys.argv[1]
destination_folder = sys.argv[2]

if not os.path.exists(destination_folder):
    os.makedirs(destination_folder)


print("Initalising...")
execute = SIR(config_file=config_file_path)

if not execute.deterministic_spawning:
    print("Configuration file not in deterministic regime. Set [deterministic_spawning: true] or use DataCollector_Random.py")
    sys.exit(1)

num_of_nodes = execute.num_nodes
number_of_variants = execute.num_of_variants
end_steps = 50000

peak_infection = np.zeros((number_of_variants, 2))
total_infections = np.zeros(number_of_variants)
total_unique_infections = np.zeros(number_of_variants)

adj_matrix_sparse = nx.adjacency_matrix(execute.G)
data_sets = [[] for _ in range(number_of_variants)]


plot_y = np.zeros((number_of_variants, 1))
plot_x = [0]
time = np.array([0])
end = False
print("Running Sim...")

while not end:
    current_time = execute.step_run()

    infected_sets = execute.infected_set
    recovered_sets = execute.recovered_set
    current_infection_num = np.array(list(map(len, infected_sets)))
    current_recovered_num = np.array(list(map(len, recovered_sets)))

    plot_y = np.hstack((plot_y, current_infection_num.reshape(-1, 1)))
    plot_x.append(current_time)

    peak_inf_mask = current_infection_num > peak_infection[:, 1]
    peak_infection[peak_inf_mask, 1] = current_infection_num[peak_inf_mask]
    peak_infection[:, 0] = np.where(peak_inf_mask, current_time, peak_infection[:, 0])

    total_infections += current_infection_num

    infection_start_condition = total_infections > 0

    for i in range(len(data_sets)):
        # > if a disease has just started, mark down the infected and recovered
        data_sets[i] = list(recovered_sets[i]) if infection_start_condition[i] \
            else data_sets[i]

    end = execute.end
    if (execute.variant_count >= 2) and (execute.steps >= end_steps):
        sim_rates = execute.rates
        is_undying = np.where(sim_rates[:, 0] > sim_rates[:, 1])[0]
        carrier_majority = current_infection_num > 0.9 * (num_of_nodes - current_recovered_num)
        leader_of_the_pack = np.argmax(current_infection_num)
        end = end or (carrier_majority[leader_of_the_pack] and is_undying[leader_of_the_pack])

print()
print("Data Dump...")

total_unique_infections += current_recovered_num

data_output_path = destination_folder
os.makedirs(data_output_path, exist_ok=True)
point_number = config_file_path.split("_")[-1]

sp.save_npz(os.path.join(data_output_path, f'networkAdjacencySparse_{point_number}.npz'), adj_matrix_sparse)

output_filename = f'datafile_{point_number}.csv'
with open(os.path.join(data_output_path, output_filename), mode='w', newline='') as file:
    writer = csv.writer(file)
    writer.writerow(['Title', 'Values'])
    writer.writerow(['Peak_Times', ', '.join(map(str, peak_infection[:, 0]))])
    writer.writerow(['Peak_Infections', ', '.join(map(str, peak_infection[:, 1]))])
    writer.writerow(['Total Infections', ', '.join(map(str, total_infections))])
    writer.writerow(['Total_Unique_Infections', ', '.join(map(str, total_unique_infections))])

for v in range(execute.num_of_variants):
    plt.plot(plot_x, plot_y[v], label=(chr(ord('A') + v)))

with open(os.path.join(data_output_path, f'recoveredData_{point_number}.csv'), 'w', newline='') as csvfile:
    writer = csv.writer(csvfile)
    for header, recovery_data in enumerate(data_sets):
        writer.writerow([f"-{chr(ord('A') + header)}"])
        writer.writerow(recovery_data)
            

plt.xlabel('Time (days)')
plt.ylabel('Number of Infected')
plt.legend()
plt.savefig(os.path.join(data_output_path, f'InfectionsVsTimePlot_{point_number}.png'))
