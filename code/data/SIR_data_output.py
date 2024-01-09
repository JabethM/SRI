import numpy as np
from ..SIR_simulation_files.SIR_system import SIR
import matplotlib.pyplot as plt


end = False
seed = 875485
num_of_variants = 16
np.random.seed(seed)
probs = tuple([(np.random.randint(0, 1000), np.random.randint(0, 1000)) for _ in range(1)])

execute = SIR(100, initialisation=(0, 0.1), variants=num_of_variants, probabilities=probs, seed=seed)

num_of_infected = np.zeros((num_of_variants, 1))
num_of_recovered = np.zeros((num_of_variants, 1))

while not end:
    end = execute.step_run()
    infected_sets = execute.infected_set
    recovered_sets = execute.recovered_set

    assert(len(infected_sets) == len(num_of_infected))
    assert(len(recovered_sets) == len(num_of_recovered))

    # lists containing the number of infected and recovered nodes for each variant
    infected_lengths = np.array([set_length for set_length in map(len, infected_sets)])
    recovered_lengths = np.array([set_length for set_length in map(len, recovered_sets)])

    num_of_infected = np.hstack((num_of_infected, infected_lengths.reshape(-1, 1)))
    num_of_recovered = np.hstack((num_of_recovered, recovered_lengths.reshape(-1, 1)))

    if execute.time > 300:
        end = True

count = len(num_of_infected[0])
x = np.arange(count)

for v in range(num_of_variants):
    plt.plot(x, num_of_infected[v], label=(chr(ord('A') + v)))

execute.root.print_tree(execute.root)

plt.xlabel('Time')
plt.ylabel('Number of infected')

plt.legend()
plt.show()
