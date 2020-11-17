# HolisticElectricityModel.jl

The Holistic Electricity Model (HEM) is a computational framework for analyzing electricity systems in their entirety, from the points of view of all key stakeholders.

## Start-Up Instructions

- Clone this repository
- Clone the [HolisticElectricityModel-Data repository](https://github.nrel.gov/HEM/HolisticElectricityModel-Data)

After this, you have some options:
- [Computational Environment](#computational-environment): Local or NREL High-Performance Computing (HPC)
- [Run Style](#run-style): REPL or command line

### Computational Environment

#### Local Set-up

- [Install the Xpress solver](https://github.nrel.gov/dcutler/fico-xpress)
- Each time vefore a model run, make sure the driver.jl file is set to use the XpressSolver

#### NREL HPC - Eagle Set-Up

- [Download](https://julialang.org/downloads/) julia linux binaries
- Extract them
- SCP to install on Eagle:
    ```
    scp -r /Users/$USER/Downloads/julia-1.5.2-linux-x86_64.tar.gz $USER@ed1.hpc.nrel.gov:/home/$USER/
    ssh ed1.hpc.nrel.gov
    cd ~
    tar -xvzf julia-1.5.2-linux-x86_64.tar.gz
    ```
- Add lines analogous to these to your ~/.bashrc file on Eagle:
    ```
    # User specific aliases and functions
    alias julia="/home/ehale/julia-1.5.2/bin/julia"
    ```
- Once on a login node, make it so the required packages will load by running the first part of the REPL code:
    ```bash
    > module load gurobi
    > julia
    ```

    ```julia
    julia> ]
    pkg> activate .
    pkg> instantiate
    pkg> build Gurobi
    ```
- Each time before a model run, be sure to load the Gurobi module:
    ```
    module load gurobi
    ```
    and make sure the driver.jl file is set to use the GurobiSolver
  
    For example, to run in interactive mode on a debug node:
    ```bash
    > salloc -N 1 -t 60 --account=mpafess --partition=debug
    # ... Wait for a compute node
    > cd HolisticElectricityModel.jl
    > module load gurobi
    > julia --project=. script/driver.jl
    ```


### Run Style

#### REPL

Open the REPL from the HolisticElectrictyModel.jl directory. Then:

```julia
julia> ]
pkg> activate .
# Hit Backspace
julia> include("script/driver.jl")
```

#### Command Line

```bash
> cd ~/HolisticElectricityModel.jl
> julia --project=. script/driver.jl
```
