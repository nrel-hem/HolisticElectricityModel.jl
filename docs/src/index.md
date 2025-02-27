# HolisticElectricityModel.jl

```@meta
CurrentModule = HolisticElectricityModel
```

## Overview

`HolisticElectricityModel.jl` (HEM) is a [`Julia`](http://www.julialang.org) package
that provides a computational framework for analyzing electricity systems in their entirety, from the points of view of all key stakeholders. Within HEM, various stakeholders operate as independent agents with different constraints, objectives, and decision-making processes. HEM iteratively captures the interactions between the different agents and the influence of their decisions on the evolution of the system under study year by year. 

*Guo, Nongchao, Elaine Hale, Jessica Lau, Jenny Heeter, Mark O’Malley, Bryan Palmintier, and Shanti Pless. 2023. “Integrated, Multi-Stakeholder Analysis of Electricity System Structures: Methodology and a Case Study.” IEEE Transactions on Energy Markets, Policy and Regulation 1 (4): 237–47. https://doi.org/10.1109/TEMPR.2023.3277345.*

## Installation

### NREL HPC Instructions
(*Skip this step if Julia is already installed on the HPC*)
- Go to https://julialang.org/downloads and copy the link to the latest supported Julia Linux image (currently Julia 1.10.0).
- ssh to NREL HPC (currently Kestrel) and download the image.
```bash
$ ssh kestrel.hpc.nrel.gov
$ wget https://julialang-s3.julialang.org/bin/linux/x64/1.8/julia-1.10.0-linux-x86_64.tar.gz
$ tar -xzf julia-1.10.0-linux-x86_64.tar.gz
```
- Add lines analogous to these to your ~/.bashrc file on the HPC:
```bash
# User specific aliases and functions
alias julia="/home/<user_name>/julia-1.10.0/bin/julia"
```
- Reload the file so that you can run `julia`.
```bash
$ source ~/.bashrc
```
- Once on a login node, make it so the required packages will load. Also, load the latest suported Gurobi solver (currently 11.0.0):
```bash
> module load gurobi/11.0.2
> julia --project=runner/Gurobi
```

```julia
julia> ]
pkg> dev .
pkg> instantiate
pkg> build Gurobi
```

### NREL Laptop Instructions

- Go to https://julialang.org/downloads/ and install Julia based on the OS.
- [Install the Xpress solver](https://github.nrel.gov/MSOC/fico-xpress)
- Update the .bashrc file similar to the HPC (depending on the OS).
- Load julia and set up the packages.
```bash
$ julia --project=runner/Xpress
```

```julia
julia> ]
pkg> dev .
pkg> instantiate
pkg> build Xpress
```