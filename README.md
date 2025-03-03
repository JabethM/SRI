# MPhys Project
This project explores the spread of infectious diseases on networks using computational modeling techniques. The focus looks at how diseases on different networks successfully propagate when there are other diseases with nonzero cross immunity coefficients actively spreading through the same network.

## Features

- Network-based disease modeling: Implements various network topologies (e.g., random graphs, scale-free networks, small-world networks).

- Epidemic simulations: Uses the Gillespie algorithm alongside the stochastic SIR model to simulate outbreaks.

- Data-driven analysis: Examines correlations between network properties infection dynamics and how these change both in time and with varying infectious rates.

- Parameter tuning & scenario testing: Allows for controlled variation of parameters such as transmission rates, recovery rates, and vaccination strategies (preset recovered/resistant/immune nodes).

## Technologies

- Programming Language: Python

- Libraries: NetworkX, NumPy, Matplotlib


## Usage

Using `src/simulation/SIR_system.py` one could import and use the main bulk of the simulation alongside `SIR.step_run()` to perform their own data analysis. Simply 
```python
from SIR_system import SIR

sim = SIR(config_file=json_configuration_file)
while True:
    sim.step_run()
```

For data collection purposes once could simply run
```commandline
python3 src/data_extraction/DataCollector_Deterministic.py config.json output_folder file_name_prefix
```
or `src/data_extraction/DataCollector_Random.py` instead


## Results

For results please consult MPhys.pdf

