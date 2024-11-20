# HolisticElectricityModel.jl

The Holistic Electricity Model (HEM) is a computational framework for analyzing electricity systems in their entirety, from the points of view of all key stakeholders.

## Start-Up Instructions

- Clone this repository
- Clone the [HolisticElectricityModelData.jl repository](https://github.nrel.gov/HEM/HolisticElectricityModelData.jl)


### Installation

#### NREL HPC Instructions
(* Skip this step if Julia is already installed in the HPC *)
- Go to https://julialang.org/downloads and copy the link to the latest supported Julia Linux image.
- ssh to NREL HPC (currently Kestrel) and download the image.
    ```bash
    $ ssh kestrel.hpc.nrel.gov
    $ wget https://julialang-s3.julialang.org/bin/linux/x64/1.8/julia-<x.y.z>-linux-x86_64.tar.gz
    $ tar -xzf julia-<x.y.z>-linux-x86_64.tar.gz
    ```
- Add lines analogous to these to your ~/.bashrc file on the HPC:
    ```
    # User specific aliases and functions
    alias julia="/home/<user_name>/julia-<x.y.z>/bin/julia"
    ```
- Reload the file so that you can run `julia`.
    ```bash
    $ source ~/.bashrc
    ```
- Once on a login node, make it so the required packages will load:
    ```bash
    > module load gurobi
    > julia --project=runner/Gurobi
    ```

    ```julia
    julia> ]
    pkg> dev --local ../HolisticElectricityModelData.jl
    pkg> instantiate
    pkg> build Gurobi
    ```

#### NREL Laptop Instructions

- Go to https://julialang.org/downloads/ and install Julia based on the OS.
- [Install the Xpress solver](https://github.nrel.gov/MSOC/fico-xpress)
- Update the .bashrc file similar to the HPC (depending on the OS).
- Load julia and set up the packages.
    ```bash
    > julia --project=runner/Xpress
    ```

    ```julia
    julia> ]
    pkg> dev --local ../HolisticElectricityModelData.jl
    pkg> instantiate
    ```

### Performing a Model Run

Running HEM requires access to LP and MIP solvers. Because different choices can 
be made in this regard, we recommend setting up your computational environment 
based on one of the environments under [runner](https://github.com/nrel-hem/HolisticElectricityModel.jl/tree/main/runner). Currently, it is recommended to use Gurobi for HPC and Xpress for local runs.

#### Running on NREL HPC

- Set up the model configurations using the hem_config.yaml file in HolisticElectricityModel.jl/script/configs/. Make sure that the solver is 'Gurobi'.

- Load the Gurobi module

    ```
    module load gurobi
    ```
    and then use the `hem.jl` script for your model runs.
  
    For example, to run in interactive mode on a debug node:
    ```bash
    > salloc -N 1 -t 60 --account=hem --partition=debug
    # ... Wait for a compute node
    > cd HolisticElectricityModel.jl
    > module load gurobi
    > julia --project=runner/Gurobi script/hem.jl
    ```

#### Running on NREL Laptop

- Ensure that Xpress is installed and activated.

- Set up the model configurations using the hem_config.yaml file in HolisticElectricityModel.jl/script/configs/. Make sure that the solver is 'Xpress'.

- Activate the runner/Xpress environment and run the hem.jl script.

    ```bash
    > cd HolisticElectricityModel.jl
    > module load gurobi
    > julia --project=runner/Xpress script/hem.jl
    ```

#### Other Environments

The key consideration for running in non-NREL computational environments is which
solver will be used for linear programs (LPs) and mixed integer programs (MIPs). 
The HEM team recommends using a commercial solver; for example, we use either 
Xpress or Gurobi. If you would prefer to use a different solver, please let us 
know and/or try implementing on your own. The key steps to supporting another 
solver are:

- Add a new environment specification under [runner](https://github.com/nrel-hem/HolisticElectricityModel.jl/tree/main/runner)
- Create a new `HEMSolver` in [src/solvers.jl](https://github.com/nrel-hem/HolisticElectricityModel.jl/blob/main/src/solvers.jl)
- Add an entry to the get_optimizer_for_solver in [src/solvers.jl](https://github.com/nrel-hem/HolisticElectricityModel.jl/blob/main/src/solvers.jl)
- Edit `script/hem.jl` accordingly

A nonlinear programming solver is also required, but we've generally found the 
open source solver Ipopt to be sufficient for our purposes. Thus, Ipopt is 
installed by default in our runner projects. That said, if you'd like to try 
other options, feel free to explore and please do reach out if you learn anything 
interesting or would like to make a contribution.


### Running Tests

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