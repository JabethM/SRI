from scripts.SIR_simulation_files.SIR_system import SIR
from scripts.SIR_simulation_files.SIR_system import Variant
import os
import numpy as np
import networkx as nx
import json

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

time = np.array([0])
end = False

while not end:
    current_time = execute.step_run()
    infected_sets = execute.infected_set
    current_infection_num = np.array(list(map(len, infected_sets)))

    # Calculate Peak Time
    # ------------------------------------------------------------------------
    mask = current_infection_num > peak_infection[:, 1]
    peak_infection[mask, 1] = current_infection_num[mask]
    peak_infection[:, 0] = np.where(mask, current_time, peak_infection[:, 0])

    # Note Beginning and End of an Epidemic
    # ------------------------------------------------------------------------
    # > Are there diseases that haven't started? And are there infected individuals of that disease?
    infection_start_condition = (total_infection_time[:, -1] == 0) & current_infection_num > 0

    # > Note down the current time, mark disease as having started (aka 1)
    total_infection_time[infection_start_condition, 0] = current_time
    total_infection_time[infection_start_condition, -1] = 1

    # > Are there diseases that have started? And are there no more infected individuals of that disease
    infection_end_condition = (total_infection_time[:, -1] == 1) & current_infection_num <= 0
    total_infection_time[infection_end_condition, 1] = current_time
    total_infection_time[infection_end_condition, -1] = 2

    # Collect Infection and Recovery Rates of each disease
    # ------------------------------------------------------------------------
    rates_collected = False
    if execute.variant_count > number_of_variants and not rates_collected:
        variant_rate_properties = execute.P
        rates_collected = True
    # ------------------------------------------------------------------------


    np.apply_along_axis(len, 1, infected_sets)
    if current_time >= end_time:
        end = True

# Any continuing disease will be marked as not having finished
undying_infections_condition = (total_infection_time[:, -1] == 1) & current_infection_num > 0
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
