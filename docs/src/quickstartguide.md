# Quick Start Guide

## Pre-requisites

- **Julia**: HolisticElectricityModel is written in Julia. Instructions to install Julia can be found in the [julia website](https://julialang.org/downloads/).  If this is your first time using Julia visit the official [Getting started with Julia](https://julialang.org/learning/).

- **Solvers**: Running HEM requires access to LP and MIP solvers. Because different choices can be made in this regard, we recommend setting up your computational environment based on one of the environments under [runner](https://github.com/nrel-hem/HolisticElectricityModel.jl/tree/main/runner). Currently, it is recommended to use Gurobi for HPC and Xpress for local runs.

- **Data**: The data for which the simulation will be run needs to be pre-populated and placed in the local directory. Detailed instructions on how to generate the run inputs can be found in the [HolisticElectricityModelData.jl](https://github.com/nrel-hem/HolisticElectricityModelData.jl) repository.


## Running Your First Model

### Running on NREL HPC

- Set up the model configurations using the hem_config.yaml file in HolisticElectricityModel.jl/script/configs/. Make sure that the solver is 'Gurobi'.

- Load the Gurobi module

```bash
> module load gurobi/11.0.0
```
    and then use the `hem.jl` script for your model runs.
  
    For example, to run in interactive mode on a debug node:
```bash
> salloc -N 1 -t 60 --account=hem --partition=debug
# ... Wait for a compute node
> cd HolisticElectricityModel.jl
> module load gurobi/11.0.0
> julia --project=runner/Gurobi script/hem.jl
```

### Running on NREL Laptop

- Ensure that Xpress is installed and activated.

- Set up the model configurations using the hem_config.yaml file in HolisticElectricityModel.jl/script/configs/. Make sure that the solver is 'Xpress'.

- Activate the runner/Xpress environment and run the hem.jl script.

```bash
$ cd HolisticElectricityModel.jl
$ julia --project=runner/Xpress script/hem.jl
```