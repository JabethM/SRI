from src.simulation.SIR_system import SIR
# from ..simulation.SIR_system import SIR
import numpy as np
import networkx as nx
import scipy.sparse as sp
import os
import csv
import matplotlib.pyplot as plt
import sys
import pickle

if len(sys.argv) != 4:
    print(f"Usage: {sys.argv[0]} <configuration_file> <output_folder> <name_suffix>")
    sys.exit(1)

config_file_path = sys.argv[1]
destination_folder = sys.argv[2]
name_suffix = sys.argv[3]

if not os.path.exists(destination_folder):
    os.makedirs(destination_folder)

print("Initalising...")
execute = SIR(config_file=config_file_path)

if not execute.deterministic_spawning:
    print(
        "Configuration file not in deterministic regime. Set [deterministic_spawning: true] or use "
        "DataCollector_Random.py")
    sys.exit(1)

plot = False

num_of_nodes = execute.num_nodes
number_of_variants = execute.num_of_variants
end_steps = 50000
lagged_time_count = 0

# Infection Statistics
peak_infection = np.zeros((number_of_variants, 2))
total_infections = np.zeros(number_of_variants)
unique_recoveries = np.zeros(number_of_variants)
unique_infections = [set() for _ in range(number_of_variants)]
population_percentile_times = np.zeros((4, number_of_variants)) - 1   # Time when disease starts, infects 25% and then 50%
population_stages = np.zeros(number_of_variants)

# Network Statistics
adj_matrix_sparse = nx.adjacency_matrix(execute.G)
immune_data = []
rec_data_time = []

# Connection Statistics
upper_Qs = []
lower_Qs = []
most_connections_list = []
least_connections_list = []
mean_connections_list = []
connection_time = []

plot_y = np.zeros((2 * number_of_variants, 1))
plot_y[number_of_variants:, 0] = num_of_nodes
plot_x = [0]
time = np.array([0])
end = False
print("Running Sim...")

while not end:
    current_time = execute.step_run()
    current_variant_count = execute.variant_count

    if current_variant_count > lagged_time_count:
        population_percentile_times[0][lagged_time_count] = current_time
        lagged_time_count += 1

    infected_sets = execute.infected_set
    recovered_sets = execute.recovered_set
    current_infection_num = np.array(list(map(len, infected_sets)))
    current_recovered_num = np.array(list(map(len, recovered_sets)))
    current_susceptible_num = num_of_nodes - (current_infection_num + current_recovered_num)

    # Diseases that have infected over 25% of the total population
    quarter_pop_mask = (current_infection_num > 0.25 * num_of_nodes) & (population_stages == 0)
    median_pop_mask = (current_infection_num > 0.5 * num_of_nodes) & (population_stages == 1)
    end_of_disease_mask = (current_infection_num == 0) & (population_percentile_times[0] >= 0)

    population_percentile_times[1][quarter_pop_mask] = current_time
    population_stages[quarter_pop_mask] = 1

    population_percentile_times[2][median_pop_mask] = current_time
    population_stages[median_pop_mask] = 2

    population_percentile_times[3][end_of_disease_mask] = current_time

    current_y_plots = np.vstack((current_infection_num.reshape(-1, 1), current_susceptible_num.reshape(-1, 1)))
    plot_y = np.hstack((plot_y, current_y_plots))
    plot_x.append(current_time)

    peak_inf_mask = current_infection_num > peak_infection[:, 1]
    peak_infection[peak_inf_mask, 1] = current_infection_num[peak_inf_mask]
    peak_infection[:, 0] = np.where(peak_inf_mask, current_time, peak_infection[:, 0])

    total_infections += current_infection_num
    unique_infections = [set1.union(set2) for set1, set2 in zip(unique_infections, infected_sets)]

    if (execute.steps % 500) == 0:
        non_negative_counts = np.sum(execute.all_current_nbs >= 0, axis=2)
        upper_quartile, lower_quartile, most_connections, least_connections = np.percentile(non_negative_counts,
                                                                                            [75, 25, 100, 0], axis=1)
        average_connections = np.mean(non_negative_counts, axis=1)

        upper_Qs.append(upper_quartile)
        lower_Qs.append(lower_quartile)
        most_connections_list.append(most_connections)
        least_connections_list.append(least_connections)
        mean_connections_list.append(average_connections)
        connection_time.append(current_time)

    if (execute.steps % 1000) == 0:
        immune_data.append(list(map(list, recovered_sets)))
        rec_data_time.append(current_time)
    if (execute.steps % 10000) == 0:
        print(current_susceptible_num)
    end = execute.end
    if (current_variant_count >= 2) and (execute.steps >= end_steps):
        sim_rates = execute.rates
        carrier_majority = current_infection_num < 0.01 * num_of_nodes
        leader_of_the_pack = np.argmax(current_infection_num)
        end = end or (carrier_majority[leader_of_the_pack])

print()
print("Data Dump...")

end_time_blanks = (current_infection_num != 0) | (population_stages < 2)
population_percentile_times[3][end_time_blanks] = current_time
unique_recoveries += current_recovered_num
unique_infections = map(len, unique_infections)
connection_statistics = [most_connections_list, upper_Qs, mean_connections_list, lower_Qs, least_connections_list]

plot_dict = {'x': plot_x, 'y': plot_y}
connection_dict = {'stats': connection_statistics, 'time': connection_time}
immune_dict = {'nodes': immune_data, 'time': rec_data_time}

data_output_path = destination_folder
os.makedirs(data_output_path, exist_ok=True)

sp.save_npz(os.path.join(data_output_path, f'networkAdjacencySparse_{name_suffix}.npz'), adj_matrix_sparse)

output_filename = f'dataFile_{name_suffix}.csv'
with open(os.path.join(data_output_path, output_filename), mode='w', newline='') as file:
    writer = csv.writer(file)
    writer.writerow(['Title', 'Values'])
    writer.writerow(['Start Time', ', '.join(map(str, population_percentile_times[0]))])
    writer.writerow(['25% prevalence rate time', ', '.join(map(str, population_percentile_times[1]))])
    writer.writerow(['50% prevalence rate time', ', '.join(map(str, population_percentile_times[2]))])
    writer.writerow(['Peak_Times', ', '.join(map(str, peak_infection[:, 0]))])
    writer.writerow(['Peak_Infections', ', '.join(map(str, peak_infection[:, 1]))])
    writer.writerow(['Total Infections', ', '.join(map(str, total_infections))])
    writer.writerow(['Total Unique Infections', ', '.join(map(str, unique_infections))])
    writer.writerow(['Total Unique Recoveries', ', '.join(map(str, unique_recoveries))])

with open(os.path.join(data_output_path, f'immuneData_{name_suffix}.pickle'), 'wb') as f:
    pickle.dump(immune_dict, f)

with open(os.path.join(data_output_path, f'connectionStatsData_{name_suffix}.pickle'), 'wb') as f:
    pickle.dump(connection_dict, f)

with open(os.path.join(data_output_path, f'plotData_{name_suffix}.pickle'), 'wb') as f:
    pickle.dump(plot_dict, f)

if plot:
    for v in range(2 * number_of_variants):
        if v < execute.num_of_variants:
            label = "Infected by " + (chr(ord('A') + v))
        else:
            label = "Susceptible to " + (chr(ord('A') + v - number_of_variants))
        plt.plot(plot_x, plot_y[v], label=label)
    plt.xlabel('Time (days)')
    plt.ylabel('Number of nodes')
    plt.legend()
    plt.savefig(os.path.join(data_output_path, f'InfectionsVsTimePlot_{name_suffix}.png'))
