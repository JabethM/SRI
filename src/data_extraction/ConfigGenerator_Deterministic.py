import json
import os
import sys
import numpy as np
import uuid

infective_range = np.linspace(5, 25, 20)
time_delay_range = np.linspace(0, 15, 20)
relationship_range = np.linspace(0, 1, 20)

if len(sys.argv) != 3:
    print("Usage: python script.py <template_file> <destination_folder>")
    sys.exit(1)

source_file = sys.argv[1]
destination_folder = sys.argv[2]

run_id = str(uuid.uuid4())
destination_folder = os.path.join(destination_folder, run_id)

if not os.path.exists(destination_folder):
    os.makedirs(destination_folder)

with open(source_file, 'r') as file:
    data = json.load(file)

initial_transmission_rate = data["input_params"]["transmission_time"]

for rel_idx, rel in enumerate(relationship_range):
    rel_directory = f'Relationship_Weight_{rel_idx}'
    rel_directory = os.path.join(destination_folder, rel_directory)
    for time_idx, time in enumerate(time_delay_range):
        time_directory = f'Time_delay_{time_idx}'
        time_directory = os.path.join(rel_directory, time_directory)
        os.makedirs(time_directory, exist_ok=True)
        for inf_idx, inf in enumerate(infective_range):
            new_file_name = f'transmission_rate_{inf_idx}.json'
            new_file_name = os.path.join(time_directory, new_file_name)

            new_data = data.copy()
            new_data["internal_params"]["deterministic_spawning"]["deterministic_mode"] = True

            new_data["internal_params"]["deterministic_spawning"]["child_transmission_T"] = [inf]
            new_data["internal_params"]["deterministic_spawning"]["time_separation"] = time
            new_data["internal_params"]["deterministic_spawning"]["relation_matrix"] = \
                [[1, rel], [rel, 1]]

            with open(new_file_name, 'w') as new_file:
                json.dump(new_data, new_file, indent=4)
print()
print(destination_folder)