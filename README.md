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
- Each time before a model run, make sure the driver.jl file is set to use the XpressSolver

#### NREL HPC - Eagle Set-Up

- Go to https://julialang.org/downloads and copy the link to the latest Julia Linux image.
- ssh to Eagle and download the image.
    ```bash
    $ ssh ed1.hpc.nrel.gov
    $ wget https://julialang-s3.julialang.org/bin/linux/x64/1.8/julia-1.8.1-linux-x86_64.tar.gz
    $ tar -xzf julia-1.8.1-linux-x86_64.tar.gz
    ```
- Add lines analogous to these to your ~/.bashrc file on Eagle:
    ```
    # User specific aliases and functions
    alias julia="/home/ehale/julia-1.5.2/bin/julia"
    ```
- Reload the file so that you can run `julia`.
    ```bash
    $ source ~/.bashrc
    ```
- Once on a login node, make it so the required packages will load by running the first part of the REPL code:
    ```bash
    > module load gurobi/9.1.2
    > julia --project=test
    ```

    ```julia
    julia> ]
    pkg> instantiate
    pkg> build Gurobi
    ```
- Each time before a model run, be sure to load the Gurobi module:
    ```
    module load gurobi/9.1.2
    ```
    and make sure the driver.jl file is set to use the GurobiSolver
  
    For example, to run in interactive mode on a debug node:
    ```bash
    > salloc -N 1 -t 60 --account=mpafess --partition=debug
    # ... Wait for a compute node
    > cd HolisticElectricityModel.jl
    > module load gurobi
    > julia --project=test script/driver.jl
    ```


### Run Style

Note that a few solver packages are installed in the `test` project and not in the main HolisticElectricityModel
package. The HolisticElectricityModel team used those solvers for test and development. The intention is to allow
users to install any compatible solver in their own environment.

#### REPL

Open the REPL from the HolisticElectrictyModel.jl directory. Then:

```julia
julia> ]
pkg> activate test
# Hit Backspace
julia> include("script/driver.jl")
```

#### Command Line

```bash
> cd ~/HolisticElectricityModel.jl
> julia --project=test script/driver.jl
```
