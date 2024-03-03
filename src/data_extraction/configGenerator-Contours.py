import json
import os
import sys
import numpy as np
import uuid

infective_range = np.linspace(5, 15, 10)
time_delay_range = np.linspace(0, 15, 10)
relationship_range = np.linspace(0, 1, 10)
independent_variables = [time_delay_range, relationship_range, infective_range]
variable_names = ["delay", "relation", "R0"]

# Mode defines which config files will be made
# Mode1 : R0.vs.Relation
# Mode2: R0.vs.TimeDelay
# Mode3: Relation.vs.TimeDelay

if len(sys.argv) != 4:
    print(f"Usage: python {sys.argv[0]} <template_file> <destination_folder> <mode>")
    sys.exit(1)

source_file = sys.argv[1]
destination_folder = sys.argv[2]
mode = int(sys.argv[3])

if not os.path.exists(destination_folder):
    os.makedirs(destination_folder)

with open(source_file, 'r') as file:
    data = json.load(file)


for var1_idx, var1 in enumerate(independent_variables[mode % 3]):
    for var2_idx, var2 in enumerate(independent_variables[(mode + 1) % 3]):
            new_file_name = f'{variable_names[mode % 3]}_{variable_names[(mode+1) % 3]}-{var1_idx}{var2_idx}.json'
            new_file_name = os.path.join(destination_folder, new_file_name)

            new_data = data.copy()
            new_data["internal_params"]["deterministic_spawning"]["deterministic_mode"] = True

            if mode == 1:
                new_data["internal_params"]["deterministic_spawning"]["relation_matrix"] = \
                    [[1, var1], [var1, 1]]
                new_data["internal_params"]["deterministic_spawning"]["child_transmission_T"] = [var2]
            elif mode == 2:
                new_data["internal_params"]["deterministic_spawning"]["child_transmission_T"] = [var1]
                new_data["internal_params"]["deterministic_spawning"]["time_separation"] = [var2]
            elif mode == 3:
                new_data["internal_params"]["deterministic_spawning"]["time_separation"] = [var1]
                new_data["internal_params"]["deterministic_spawning"]["relation_matrix"] = \
                    [[1, var2], [var2, 1]]


            with open(new_file_name, 'w') as new_file:
                json.dump(new_data, new_file, indent=4)
print(destination_folder)