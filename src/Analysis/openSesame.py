import os
import sys
import numpy as np
import pickle

infective_range = np.linspace(5, 15, 10)
time_delay_range = np.linspace(0, 15, 10)
relationship_range = np.linspace(0, 1, 10)
independent_variables = [time_delay_range, relationship_range, infective_range]
variable_names = ["delay", "relation", "R0"]

# Mode defines which config files will be made
# Mode1 : R0.vs.Relation
# Mode2: R0.vs.TimeDelay
# Mode3: Relation.vs.TimeDelay

if len(sys.argv) != 2:
    print(f"Usage: python {sys.argv[0]} <data_path>")
    sys.exit(1)

destination_folder = sys.argv[1]

if not os.path.exists(destination_folder):
    print(f"{destination_folder} does not exist")
    sys.exit(1)

with open(f'{sys.argv[1]}', 'rb') as f:
    loaded_data = pickle.load(f)

print(loaded_data)

