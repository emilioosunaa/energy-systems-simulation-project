# Coupled Simulation of a Power and Gas Network

This project implements a coupled simulation of a 4-bus power grid and a 4-node gas network using the Newton-Raphson method. It was developed as part of the Energy Systems Simulation course at RWTH Aachen University.

![System Diagram](power_and_gas_system.png)

## Overview

The simulation models the interdependency between a power grid and a gas network:

- A **gas-fired turbine** at Bus 2 consumes natural gas from the gas network and generates electrical power for the grid.
- A **compressor station** between nodes 2 and 3 in the gas network is powered electrically from Bus 3 of the power grid.

This coupling requires solving both systems together, as gas demand affects power generation and compressor power consumption feeds back into the electrical load.

## Project Structure

| File | Description |
|---|---|
| `coupled_simulation.ipynb` | Main notebook running gas, power, and co-simulation scenarios |
| `power_flow.py` | Newton-Raphson power flow solver for a 4-bus system |
| `gas_network.py` | Newton-Raphson gas network solver for a 4-node pipeline system |
| `power_and_gas_system.png` | Diagram of the coupled system |

## Gas Network

The gas network consists of 4 nodes connected by 3 pipelines, with a compressor station between nodes 2 and 3. Steady-state volumetric flow rates are computed using the 1D Isothermal Euler Equations:

- **Node 1**: Supply node (50 bar)
- **Nodes 2 & 3**: Demand nodes (grouped, 60 sm³/s combined)
- **Node 4**: Turbine supply node
- **Compressor**: Between nodes 2 and 3, compression ratio r = 1.1

## Power Grid

The power grid is a 4-bus system with 5 transmission lines. The power flow is solved using the Newton-Raphson method with the full Jacobian matrix:

- **Bus 1**: Slack bus (1.05 p.u., 115.5 kV)
- **Bus 2**: PV bus (generator, 1.02 p.u.)
- **Bus 3**: PQ bus (200 MW, 50 MVAR load)
- **Bus 4**: PQ bus (100 MW, 20 MVAR load)

Line parameters: R = 109 mΩ/km, X = 1.2 mH/km, C = 9.5 nF/km at 50 Hz.

## Simulation Scenarios

1. **Independent simulations**: Gas network and power flow solved separately.
2. **Coupled simulation**: Compressor power from the gas network is added as load at Bus 3, and turbine gas consumption determines power generation at Bus 2.

## Requirements

- Python 3
- NumPy
- SciPy
- Jupyter Notebook

## Usage

```bash
jupyter notebook coupled_simulation.ipynb
```

Or run the solvers independently:

```bash
python power_flow.py
python gas_network.py
```

## Author

Emilio Osuna Aguilar
