# HolisticElectricityModel.jl

The Holistic Electricity Model (HEM) is a computational framework for analyzing electricity systems in their entirety, from the points of view of all key stakeholders.

## Workflow

1. Prepare input data using the [HolisticElectricityModelData.jl repository](https://github.nrel.gov/HEM/HolisticElectricityModelData.jl).
2. Manually edit the `scripts/hem_config.yaml` configuration file, specifying options for one's desired HEM model run.
3. Run the `scripts/hem.jl` script.

To support this workflow, one must first set-up the necessary Julia environment, which we outline below.

## Start-Up Instructions

1. Clone this repository
2. Clone the [HolisticElectricityModelData.jl repository](https://github.nrel.gov/HEM/HolisticElectricityModelData.jl)

### Computational Environment
HEM may be run locally or in NREL's High-Performance Computing environment, Kestrel (if one has access to do so).

#### Local Set-up
##### Xpress
If you have access to the Xpress solver, then you may choose to use it instead of Ipopt.
1. [Install the Xpress solver](https://github.nrel.gov/dcutler/fico-xpress)
2. Specify "Xpress" as the solver in the `hem_config.yaml` configuration file.

#### NREL HPC - Kestrel Set-Up
First, [set up an HPC account with a username and password](https://www.nrel.gov/hpc/user-accounts.html).

1. ssh to Kestrel:
    ```bash
    $ ssh <username>@kestrel.hpc.nrel.gov
    # Enter your password, and accept the terms and conditions when prompted to do so.
    ```

2. Load Julia v1.7.2, so one can call `julia` from the command line:
    ```bash
    $ module load julia/1.7.2
    ```
As a result of the above, when calling `julia` from the command line, version 1.7.2 will automatically be used.

3. Load Gurobi, so the JuMP code in the HEM package can use the Gurobi solver:
    ```bash
    $ module load gurobi
    ```

4. Set-up the Julia project (start Julia by calling `julia` from the command line):
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

5. Return to the bash command line (hit `<Backspace>` to exit the package manager; enter `exit()` or hit `<CTRL + d>` from the Julia REPL).

6. Submit a job to run HEM. Below, we grab a debug node for an hour, then run HEM in the following step:
    ```bash
    $ srun --time=1:00:00 --account=hem --nodes=1 --partition=debug --pty $SHELL
    ```

7. Run HEM using the Gurobi solver (do so from the top level of the HolisticElectricityModel.jl directory):
    ```bash
    $ julia --project=. scripts/hem.jl
    ```
- If you prefer running HEM from the Julia REPL, then call `julia --project=.` from the bash command line to start Julia with the correct project, then replace the above with the following:
    ```julia
    julia> include("scripts/hem.jl")
    ```

#### Adjusting HEM's Parameters
The above 'Quick Start' guide uses the default HEM configuration. HEM is designed to facilitate comparisons arising from different market structures, time scales, balancing areas, net metering policies, rate designs, etc. As such, users are encouraged to try out a variety of scenarios. To do so, edit the `scripts/hem_config.yaml` configuration file before running `scripts/hem.jl`. See the comments in the configuration file for a list of valid options for each variable.

#### Tests
Run tests
    ```bash
    $ julia --project=. tests/runtests.jl
    ```