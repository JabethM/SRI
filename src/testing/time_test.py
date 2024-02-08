from line_profiler import LineProfiler
import datetime

from scripts.SIR_simulation_files.SIR_system import SIR
import numpy as np

if __name__ == '__main__':
    seed = 87548
    np.random.seed(seed)
    total_variants = 7
    probs = tuple([(np.random.randint(0, 1000), np.random.randint(0, 1000)) for _ in range(1)])
    execute = SIR(50, initialisation=(0, 0.1), variants=total_variants, probabilities=probs, seed=seed)

    print("time")
    profiler = LineProfiler()
    profiler.add_function(SIR.iterate)
    profiler.add_function(SIR.run)

    profiler.add_function(SIR.choose_outcome)
    profiler.add_function(SIR.choice_probability_setup)
    profiler.add_function(SIR.node_neighbours)
    profiler.add_function(SIR.set_node_infected)
    profiler.add_function(SIR.set_node_recovery)

    profiler.run("execute.run()")
    current_datetime = datetime.datetime.now()
    formatted_datetime = current_datetime.strftime("%m-%d_%H-%M-%S")
    filename = f"ProfileData_{formatted_datetime}.lprof"
    profiler.dump_stats(filename)

    profiler.print_stats()


