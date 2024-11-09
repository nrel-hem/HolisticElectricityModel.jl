# HolisticElectricityModel.jl

The Holistic Electricity Model (HEM) is a computational framework for analyzing electricity systems in their entirety, from the points of view of all key stakeholders.

## Pre-requisites

### Julia

- Go to https://julialang.org/downloads and copy the link to the latest Julia Linux image.
- Download the image.
    ```bash
    $ ssh ed1.hpc.nrel.gov
    $ wget https://julialang-s3.julialang.org/bin/linux/x64/1.8/julia-1.8.1-linux-x86_64.tar.gz
    $ tar -xzf julia-1.8.1-linux-x86_64.tar.gz
    ```
- Add lines analogous to these to your ~/.bashrc file:
    ```
    # User specific aliases and functions
    alias julia="/home/<user_name>/julia-1.5.2/bin/julia"
    ```
- Reload the file so that you can run `julia`.
    ```bash
    $ source ~/.bashrc
    ```


## Start-Up Instructions

- Clone this repository
- Clone the [HolisticElectricityModelData.jl repository](https://github.nrel.gov/HEM/HolisticElectricityModelData.jl) and follw the instructions in the README to get the latest input data files.
- Within the HolisticElectrictyModelData.jl directory, the `parse_inputs()` function needs to be executed to populate the `runs/` directory with the parsed input data for running a test case. This can be done by running `src/input_data_parsing_enhanced.jl` (after uncommenting the line that calls the `parse_inputs()` function). Another option is to run the `scripts/check_input_parse.jl` script (after commenting all the code after the call to `parse_inputs()`).

After this, you have some options:
- [Computational Environment](#computational-environment): Local or NREL High-Performance Computing (HPC)
- [Run Style](#run-style): REPL or command line

### Computational Environment

#### Local Set-up

- [Install the Xpress solver](https://github.nrel.gov/MSOC/fico-xpress)
- Install the required project dependencies.
    ```cd HolisticElectricityModel.jl
    julia --project
    julia> ]
    pkg> dev ../HolisticElectricityModelData.jl
    pkg> instantiate
    pkg> build Xpress
    ```
- Use the driver_xpress.jl file for your model runs
    ```bash
    > julia --project script/driver_xpress.jl
    ``` 

#### NREL HPC - Kestrel Set-Up

- Run the prereuisite steps to install Julia if it is not installed.
- Once on a login node, make it so the required packages will load by running the first part of the REPL code:
    ```bash
    > module load gurobi/9.1.2
    > julia --project
    ```

    ```julia
    julia> ]
    pkg> dev ../HolisticElectricityModelData.jl
    pkg> instantiate
    pkg> build Gurobi
    ```
- Each time before a model run, be sure to load the Gurobi module:
    ```
    module load gurobi/9.1.2
    ```
    and use the `driver_gurobi.jl` script for your model runs.
  
    For example, to run in interactive mode on a debug node:
    ```bash
    > salloc -N 1 -t 60 --account=hem --partition=debug
    # ... Wait for a compute node
    > cd HolisticElectricityModel.jl
    > module load gurobi/9.1.2
    > julia --project script/driver.jl
    ```


### Run Style

Note that a few solver packages are installed in the `runner` project and not in the main HolisticElectricityModel
package. The HolisticElectricityModel team used those solvers for test and development. The intention is to allow
users to install any compatible solver in their own environment.

To set up the `runner` environment before running in either of the manners described below, complete the following
set-up steps.

First, if you don't have Gurobi installed on your system, set the GUROBI_JL_SKIP_LIB_CHECK environment variable:

*Linux/Mac*
```bash
$ export GUROBI_JL_SKIP_LIB_CHECK=1
```

*Windows*
```
> set GUROBI_JL_SKIP_LIB_CHECK=1
```

Or, if you don't have Xpress installed on your system, set the XPRESS_JL_SKIP_LIB_CHECK environment variable:

*Linux/Mac*
```bash
$ export XPRESS_JL_SKIP_LIB_CHECK=1
```

*Windows*
```
> set XPRESS_JL_SKIP_LIB_CHECK=1
```

Second, install HolisticElectricityModel.jl and HolisticElectricityModelData.jl in the runner project:

```bash
> julia
```
```julia
julia> ]
pkg> activate runner
pkg> dev . ../HolisticElectricityModelData.jl
```

#### REPL

Open the REPL from the HolisticElectrictyModel.jl directory. Then:

```julia
julia> ]
pkg> activate runner
# Hit Backspace
julia> include("script/driver_xpress.jl") # or include("script/driver_gurobi.jl")
```

or, if you specify the runner project when you open the REPL:

```bash
> cd ~/HolisticElectricityModel.jl
> julia --project=runner
```

and if the runner environment is already set up then you can simply:

```julia
julia> include("script/driver_xpress.jl") # or include("script/driver_gurobi.jl")
```

#### Command Line

```bash
> cd ~/HolisticElectricityModel.jl
> julia --project=runner script/driver_xpress.jl # or script/driver_gurobi.jl
```

### Tests

Specify a scenario you want to test by setting `driver_name` as desired at the top of `test/test_driver.jl` and then editing the driver file accordingly. Currently the model outputs available to test against are in the `test/driver_outputs/ba_1_base_2018_future_2_ipps_1` directory created with input parameters:

```
ba = ["p13"]
ba_len = length(ba)
base_year = 2018
future_years = [2019, 2020]
future_years_len = length(future_years)
ipp_number = 1
```

The option names printed in the `test/driver_outputs/ba_1_base_2018_future_2_ipps_1/Results_...` directory names indicate what sets of HEMOptions and RegulatoryOptions can be tested against.

Then, to run all tests:
```bash
> julia --project=test
julia> ]
pkg> dev . ../HolisticElectricityModelData.jl
# Hit Backspace
julia> include("test/runtests.jl")
```

NREL Software Record SWR-21-62 "HEM (HolisticElectricityModel.jl)"
