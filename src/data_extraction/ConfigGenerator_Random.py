import json
import os
import sys
import uuid


if len(sys.argv) != 3:
    print("Usage: python ConfigGenerator_Random.py <source_file> <destination_folder>")
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

for i in range(0, 10):
    new_data = data.copy()
    new_data["input_params"]["transmission_time"] = initial_transmission_rate + i

    new_file_name = f'config_{i}.json'
    new_file_path = os.path.join(destination_folder, new_file_name)

    with open(new_file_path, 'w') as new_file:
        json.dump(new_data, new_file, indent=4)
    
    print(f'File "{new_file_name}" created with .transmissions value: {new_data["input_params"]["transmission_time"]}')
