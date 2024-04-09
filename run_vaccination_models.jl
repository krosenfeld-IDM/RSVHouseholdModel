#This runs the forecasting of the model using inferred parameters
#Default is to involve 5 processors to accelerate forecasting by simulating forwards in parallel before pooling result
using Revise
using DifferentialEquations
using Distributed

# #Dependencies for running simulations in parallel
# addprocs(7)
# @everywhere cd(pwd())
# @everywhere include("Basic_params.jl");
# @everywhere include("Basic_vac_params.jl");
# @everywhere include("SetupStates_HHModel.jl");
# @everywhere include("ParamDefn_HHModel.jl");
# @everywhere include("Params_vacc_model.jl");
# @everywhere include("ReadInData_HHModel.jl");
# @everywhere include("LoadOptimalParams.jl");
# @everywhere include("HH_distrib.jl");
# @everywhere include("AgeDistributionModel_HMModel.jl");
# @everywhere include("ConstructStateFlow_HHModel.jl");
# @everywhere include("VacHouseholds.jl");
# @everywhere include("VectorFieldForHHModel.jl");
# @everywhere include("FutureSimulations.jl");
# include("PlotsForModel.jl");

#Dependencies for running simulations in serial
cd(pwd())
include("Basic_params.jl");
include("Basic_vac_params.jl");
include("SetupStates_HHModel.jl");
include("ParamDefn_HHModel.jl");
include("Params_vacc_model.jl");
include("ReadInData_HHModel.jl");
include("LoadOptimalParams.jl");
include("HH_distrib.jl");
include("AgeDistributionModel_HMModel.jl");
include("ConstructStateFlow_HHModel.jl");
include("VacHouseholds.jl");
include("VectorFieldForHHModel.jl");
include("FutureSimulations.jl");
include("PlotsForModel.jl");


# RepeatedForecasts = MonteCarloProblem(prob_vac,
#                                       output_func=output_func_for_rsv_inf_endpoint,
#                                       prob_func=prob_func_for_rsv_sims)

RepeatedForecasts = EnsembleProblem(prob_vac,
                                output_func=output_func_for_rsv_inf_endpoint,
                                prob_func=prob_func_for_rsv_sims)

#Comment/uncomment to choose which set of scenarios to loop over
include("Loopovervaccinationscenarios.jl");
# include("numberofvaccines.jl");
# include("loopfortrueinfs.jl")
