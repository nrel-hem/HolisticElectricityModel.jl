# HolisticElectricityModel.jl

The Holistic Electricity Model (HEM) is a computational framework for analyzing electricity systems in their entirety, from the points of view of all key stakeholders.

## Start-Up Instructions

- Clone this repository
- Clone the [HolisticElectricityModelData.jl repository](https://github.nrel.gov/HEM/HolisticElectricityModelData.jl)

After this, you have some options:
- [Computational Environment](#computational-environment): Local or NREL High-Performance Computing (HPC)
- [Run Style](#run-style): REPL or command line

### Computational Environment

#### Local Set-up

##### Open-Source Solver
TODO. Ipopt.

##### Xpress
If you have access to the Xpress solver, then you may choose to use it instead of Ipopt.
- [Install the Xpress solver](https://github.nrel.gov/dcutler/fico-xpress)
- Use the driver_xpress.jl file for your model runs

#### NREL HPC - Kestrel Set-Up
First, [set up an HPC account with a username and password](https://www.nrel.gov/hpc/user-accounts.html).

- ssh to Kestrel:
    ```bash
    $ ssh <username>@kestrel.hpc.nrel.gov
    # Enter your password, and accept the terms and conditions when prompted to do so.
    ```
- Navigate to the HEM project directory:
    ```bash
    $ cd /kfs2/projects/hem/HolisticElectricityModel.jl
    ```
Note: Logan's directory is set-up for a Quick Start run. We will assume one has navigated to the below subdirectory for the remainder of our Quick Start instructions.
    ```bash
    $ cd /kfs2/projects/hem/Github/Graham/HolisticElectricityModel.jl
    ```

- Load Julia v1.7.2, so one can call `julia` from the command line:
    ```bash
    $ module load julia/1.7.2
    ```
As a result of the above, when calling `julia` from the command line, version 1.7.2 will automatically be used.
- Load Gurobi, so the JuMP code in the HEM package can use the Gurobi solver:
    ```bash
    $ module load gurobi
    ```

- Set-up the Julia project (start Julia by calling `julia` from the command line):
    ```julia
    julia> ]
    pkg> activate .
    pkg> instantiate
    ```

- If there are issues with `instantiate`, then run:
    ```julia
    julia> ]
    pkg> update
    pkg> precompile
    ```
before retrying the above `instantiate` call.

Return to the bash command line (hit `<Backspace>` to exit the package manager; enter `exit()` or hit `<CTRL + d>` from the Julia REPL).

- Grab a debug node for an hour:
    ```bash
    $ srun --time=1:00:00 --account=hem --nodes=1 --partition=debug --pty $SHELL
    ```

- Run HEM using the Gurobi solver (do so from /kfs2/projects/hem/Github/Graham/HolisticElectricityModel.jl):
    ```bash
    $ julia --project=. scripts/driver_gurobi.jl
    ```
If you prefer running HEM from the Julia REPL, then call `julia --project=.` from the bash command line to start Julia with the correct project, then replace the above with the following:
    ```julia
    julia> include("scripts/driver_gurobi.jl")
    ```
