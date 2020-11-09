# HolisticElectricityModel.jl

The Holistic Electricity Model (HEM) is a computational framework for analyzing electricity systems in their entirety, from the points of view of all key stakeholders.

## Start-Up Instructions

- Clone this repository
- Clone the [HolisticElectricityModel-Data repository](https://github.nrel.gov/HEM/HolisticElectricityModel-Data)

After this, you have some options:
- [Computational Environment](#computational-environment): Local or NREL High-Performance Computing (HPC)
- [Run Style](#run-style): REPL or command line

### Computational Environment

#### Local

- [Install the Xpress solver](https://github.nrel.gov/dcutler/fico-xpress)

#### NREL HPC - Eagle

- [Download](https://julialang.org/downloads/) julia linux binaries
- Extract them
- SCP to install on Eagle:
    ```
    scp -r /Users/dsigler/Downloads/julia-1.5.2 dsigler@ed1.hpc.nrel.gov:/home/dsigler/julia-1.5.2
    ```
- Add lines analogous to these to your ~/.bashrc file on Eagle:
    ```
    # User specific aliases and functions
    alias julia="/home/dsigler/julia-1.5.2/bin/julia"
    ```
- Once on a login node, make it so the required packages will load by running the first part of the REPL code:
    ```julia
    julia> ]
    pkg> activate .
    ```
- Each time before a model run, be sure to load the Gurobi module:
    ```
    module load gurobi
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
> julia --project=. script\driver.jl
```
