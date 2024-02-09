from ..simulation.SIR_system import SIR
import numpy as np
import os
import csv
import matplotlib.pyplot as plt
import sys
import uuid

if len(sys.argv) != 3:
    print("Usage: DataCollector_LargeSet.py <configuration_file> <output_folder>")
    sys.exit(1)

config_file_path = sys.argv[1]
destination_folder = sys.argv[2]


print("Initalising...")
execute = SIR(config_file=config_file_path)
num_of_nodes = execute.num_nodes
number_of_variants = execute.num_of_variants
end_steps = 50000

peak_infection = np.zeros((number_of_variants, 2))
total_infections = np.zeros(number_of_variants)
total_unique_infections = np.zeros(number_of_variants)

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

    end = execute.end
    if (execute.variant_count >= 2) and (execute.steps >= end_steps):
        sim_rates = execute.rates
        is_undying = np.where(sim_rates[:, 0] > sim_rates[:, 1])[0]
        carrier_majority = current_infection_num > 0.9 * (num_of_nodes - current_recovered_num)
        leader_of_the_pack = np.argmax(current_infection_num)
        end = end or (carrier_majority[leader_of_the_pack] and is_undying[leader_of_the_pack])

data_folder_name = str(uuid.uuid4())
data_output_path = os.path.join(destination_folder, data_folder_name)
os.makedirs(data_output_path, exist_ok=True)

total_unique_infections += current_recovered_num
point_number = config_file_path.split("_")[-1]
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
plt.xlabel('Time (days)')
plt.ylabel('Number of Infected')
plt.legend()
plt.savefig(os.path.join(data_output_path, 'Infections_vs_Time_plot.png'))
