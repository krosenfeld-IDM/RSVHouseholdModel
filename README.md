# RSVHouseholdModel

Julia 0.6 code for running the RSV Household transmission model (see https://www.biorxiv.org/content/10.1101/569335v1.abstract for preprint).

The main files for inferring parameters, running simulations and running forecasts are:

* Fit_HHModel_EMalg.jl --- This contains the code used to define the simulation parameters, the log-likelihood with respect to hospitalisation data, and the EM algorithm method.
* Fit_HHModel_EMalg.jl --- This contains additional definitions for running the simulation with vaccination and uses multiple processes for averaging over future simulations.
* RSV_household_demo.jl --- This demonstrates the model running with parameters inferred from the EM algorithm (maximum likelihood estimates are saved in DATA folder), and compares to data.

The MATLAB® code used to generate the publication plots, and the outputed data behind these, are located in the Plots folder.

## Setup

Open `julia` from the terminal and install the dependencies:
```
julia> Pkg.add("Distributions"); Pkg.add("DifferentialEquations"); Pkg.add("Plots"); Pkg.add("JLD"); Pkg.add("Hiccup"); Pkg.add("BlackBoxOptim")
```
