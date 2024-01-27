from scripts.SIR_simulation_files.SIR_system import SIR
from scripts.SIR_simulation_files.SIR_system import Variant
import os
import numpy as np
import networkx as nx
import json
import matplotlib

current_directory = os.getcwd()
config_file = "config_sim.json"
parent_folder = os.path.join('..', 'SIR_simulation_files')

relative_path = os.path.relpath(parent_folder, current_directory)
file_path = os.path.join(relative_path, config_file)

execute = SIR(config_file=file_path)

number_of_variants = execute.num_of_variants
end_time = execute.end_time

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

time = np.array([0])
end = False

while not end:
    current_time = execute.step_run()
    previous_time = execute.old_time

    infected_sets = execute.infected_set
    recovered_sets = execute.recovered_set
    current_infection_num = np.array(list(map(len, infected_sets)))
    current_recovered_num = np.array(list(map(len, recovered_sets)))

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

    # Store infected and recovered population data
    # ------------------------------------------------------------------------
    for i in range(len(data_sets)):
        # > if a disease has just started, mark down the infected and recovered
        data_sets[i][0] = [infected_sets[i].copy(), recovered_sets[i].copy()] if infection_start_condition[i] \
            else data_sets[i][0]
        data_sets[i][1] = [infected_sets[i].copy(), recovered_sets[i].copy()] if peak_inf_mask[i] else data_sets[i][1]

    # Collect Infection and Recovery Rates of each disease
    # ------------------------------------------------------------------------
    rates_collected = False
    if (execute.variant_count >= number_of_variants) and not rates_collected:
        variant_rate_properties = execute.rates
        rates_collected = True

        relationship_matrix = Variant.relation_matrix
    # ------------------------------------------------------------------------

    if current_time >= end_time:
        end = True

# Any continuing disease will be marked as not having finished
undying_infections_condition = (total_infection_time[:, -1] == 1) & (current_infection_num > 0)
total_infection_time[undying_infections_condition, 1] = np.nan
total_infection_time[undying_infections_condition, -1] = 2

assert(np.all(total_infection_time[:, -1] == 2))

# Saving the variant family tree
family_tree = execute.variants[0]

with open("family_tree.json", "w") as file:
    json.dump(family_tree, file)

# Saving the variant relationship matrix
relationship_matrix = Variant.relation_matrix

#### Outputs ####
## For Each Run: ##
# - Variant Relationship Matrix
# - Variant Family Tree
# - Network adjacency matrix
## For Each Disease: ##
# - Peak Infection Time
# - Total infection time/ Length of infection
# - Infection and Recovery Rates
# - Shortened / Contracted Adjacency matrix for each peak and t_0
