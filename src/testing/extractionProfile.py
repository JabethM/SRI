from line_profiler import LineProfiler
import datetime

from src.simulation.SIR_system import SIR
from src.data_extraction.DataCollector_Deterministic import DCmain

import os

if __name__ == '__main__':
    print("hello")
    input_folder = "default_configuration_test.json"
    output_folder = os.path.join("src", "testing", "output_folder")
    suffix = "tester"
    print("time")
    profiler = LineProfiler()
    profiler.add_function(SIR.iterate)
    profiler.add_function(SIR.step_run)

    profiler.add_function(SIR.choose_outcome)
    profiler.add_function(SIR.choice_probability_setup)
    profiler.add_function(SIR.node_neighbours)
    profiler.add_function(SIR.set_node_infected)
    profiler.add_function(SIR.set_node_recovery)
    profiler.add_function(DCmain)
    profiler.run("DCmain(input_folder, output_folder, suffix)")
    current_datetime = datetime.datetime.now()
    formatted_datetime = current_datetime.strftime("%m-%d_%H-%M-%S")
    filename = f"ProfileData_{formatted_datetime}.lprof"
    profiler.dump_stats(filename)

    profiler.print_stats()


