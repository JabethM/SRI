from scripts.SIR_simulation_files.SIR_system import SIR
from argparse import ArgumentParser
import numpy as np
import os
import csv
import matplotlib.pyplot as plt

par = ArgumentParser()

par.add_argument("-f", "--config-file", type=str, help="Path for JSON config file", required=True)
args = par.parse_args()
config_file_path = args.config_file

data_output_folder = os.path.join(os.path.dirname(config_file_path), f"output_files")
os.makedirs(data_output_folder, exist_ok=True)

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

total_unique_infections += current_recovered_num
point_number = config_file_path.split("_")[-1]
output_filename = f'datafile_{point_number}.csv'
with open(os.path.join(data_output_folder, output_filename), mode='w', newline='') as file:
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
plt.savefig(os.path.join(data_output_folder, 'Infections_vs_Time_plot.png'))
