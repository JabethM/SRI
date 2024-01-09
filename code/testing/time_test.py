import sys

print(sys.executable)

"""from line_profiler import LineProfiler
from ..SIR_simulation_files.SIR_system import SIR
import numpy as np

if __name__ == '__main__':
    seed = 38204518
    np.random.seed(seed)
    total_variants = 7
    probs = tuple([(np.random.randint(0, 1000), np.random.randint(0, 1000)) for _ in range(1)])
    execute = SIR(50, initialisation=(0, 0.1), variants=total_variants, probabilities=probs, seed=seed)

    profiler = LineProfiler()
"""
