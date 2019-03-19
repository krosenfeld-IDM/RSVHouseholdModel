#Demonstration code for model

cd(pwd())
# Basic_params.jl contains the basic params of RSV, change SchoolsOnlyR_0 variable to 0,0.5,1,1.5 as desired
#Also defines a lookback variable which determines the time window over which RSV hospitalisations were binned (default 1 week)
include("Basic_params.jl");
# SetupStates_HHModel.jl contains code to automatically generate all household configurations and produce a
# lookup table for household configurations
include("SetupStates_HHModel.jl");
# ParamDefn_HHModel.jl contains code to that defines a struct containing all the parameters and other
# data types that are needed to run the household model. The target inference parameters get updated here.
include("ParamDefn_HHModel.jl");
# ReadInData_HHModel.jl reads in all the data required to initialise model and determine the yearly varying rates of
# household growth, population turnover and the empirical person type distribution. Also reads in all hospitalisation data and converts
include("ReadInData_HHModel.jl");
#HH_distrib.jl uses root-finding to convert the yearly household size distributions into turnover rates
include("HH_distrib.jl");
#AgeDistributionModel_HMModel.jl uses the empirical person type distribution to construct yearly conditional distributions
#e.g. probability of age category given household size and sharing the household with an U1
include("AgeDistributionModel_HMModel.jl");
#ConstructStateFlow_HHModel.jl sets up the matrix representation of the epidemic model by specifying the rate of transistion between
#any two possible household configurations
include("ConstructStateFlow_HHModel.jl");
#VectorFieldForHHModel.jl defines a vector field for the household model, the callback times where the yearly parameters change and
#incidence rates
include("VectorFieldForHHModel.jl");
#PlotsForModel.jl defines a number of plotting options for solutions to the household RSV transmission model
include("PlotsForModel.jl");
#ObjFunctions_HHModel.jl includes all the functions for turning solutions into a log-likelihood
include("ObjFunctions_HHModel.jl");
#EMalgs_HHModel.jl contains the E-step and M-step methods for EM algorithm inference
include("EMalgs_HHModel.jl");



# Load optimal parameters for chosen schools R_0 value chosen as 0,0.5,1,1.5
include("LoadOptimalParams.jl");

# Solve the model
neg_LL_opt,sol_rsv = neg_objective_for_params(Θ_Opt)

#Plot the model
x,y,z,p = create_comparison_plot(sol_rsv)
display(p)
