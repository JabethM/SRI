import numpy as np
import networkx as nx
from SIR_system import SIR
import matplotlib.pyplot as plt


end = False
seed = 39
num_of_variants = 7
np.random.seed(seed)
probs = tuple([(np.random.randint(0, 1000), np.random.randint(0, 600)) for _ in range(7)])

execute = SIR(50, initialisation=(0, 0.1), variants=num_of_variants, probabilities=probs, seed=seed)

num_of_infected = np.zeros((7, 1))
num_of_recovered = np.zeros((7, 1))

while not end:
    end = execute.step_run()
    infected_sets = execute.infected_set
    recovered_sets = execute.recovered_set

    assert(len(infected_sets) == len(num_of_infected))
    assert(len(recovered_sets) == len(num_of_infected))

    infected_lengths = np.array([set_length for set_length in map(len, infected_sets)])
    recovered_lengths = np.array([set_length for set_length in map(len, recovered_sets)])

    num_of_infected = np.hstack((num_of_infected, infected_lengths.reshape(-1, 1)))
    num_of_recovered = np.hstack((num_of_recovered, recovered_lengths.reshape(-1, 1)))

count = len(num_of_infected[0])
x = np.arange(count)

for v in range(num_of_variants):
    plt.plot(x, num_of_infected[v], label=('variant' + str(v)))

plt.xlabel('Time')
plt.ylabel('Number of infected')

plt.legend()
plt.show()
