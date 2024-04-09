#This runs the inference algorithm using the KCH hospitalisation data to fit for the unobserved
#random seasonality and parameters using an EM algorithm

cd(pwd())
include("Basic_params.jl");
include("SetupStates_HHModel.jl");
include("ParamDefn_HHModel.jl");
include("ReadInData_HHModel.jl");
include("HH_distrib.jl");
include("AgeDistributionModel_HMModel.jl");
include("ConstructStateFlow_HHModel.jl");
include("VectorFieldForHHModel.jl");
include("PlotsForModel.jl");
include("ObjFunctions_HHModel.jl");
include("EMalgs_HHModel.jl");




#Code for EM algorithm --- starting point chosen from optimising on a limited data set
(Θ_init_opt,res_init) = perform_initial_M_step_det_bb(:adaptive_de_rand_1_bin_radiuslimited,200,200)

CurrBestNegLL = Inf
Θ_Opt = Θ_init_opt
Z_ξ_opt = P_ModelParams.Z_ξ
Z_ϕ_opt = P_ModelParams.Z_ϕ
ξ̄_opt = P_ModelParams.ξ̄
ϕ̄_opt = P_ModelParams.ϕ̄
σ_ξ_opt = P_ModelParams.σ_ξ
σ_ϕ_opt = P_ModelParams.σ_ϕ
ρ_ξϕ_opt = P_ModelParams.ρ_ξϕ
P_currbest = P_ModelParams

#loop over restart algorithm
for i = 1:10
    println("Starting iteration $(i) of the EM algorithm")
    Perform_E_step(100)
    (Θ_Opt_poss,neg_LL_opt_poss) = Perform_M_step(600, 100)

    if neg_LL_opt_poss < CurrBestNegLL
        Θ_Opt = Θ_Opt_poss
        CurrBestNegLL = neg_LL_opt_poss
        Z_ξ_opt = P_ModelParams.Z_ξ
        Z_ϕ_opt = P_ModelParams.Z_ϕ
        ξ̄_opt = P_ModelParams.ξ̄
        ϕ̄_opt = P_ModelParams.ϕ̄
        σ_ξ_opt = P_ModelParams.σ_ξ
        σ_ϕ_opt = P_ModelParams.σ_ϕ
        ρ_ξϕ_opt = P_ModelParams.ρ_ξϕ
        P_currbest = P_ModelParams
        println("Found improvement on iteration $(i) of the EM algorithm; params = ",Θ_Opt,", LL = ",CurrBestNegLL)
    else
        println("No improvement on iteration $(i) of the EM algorithm, ending with params = ",Θ_Opt,", LL = ",CurrBestNegLL)
        P_ModelParams = P_currbest
        break
    end

    P_ModelParams = P_currbest
end


P_ModelParams = P_currbest
Perform_E_step(100)
neg_LL_opt,sol_rsv = neg_objective_for_params(Θ_Opt)
x_0_end = sol_rsv[:,end]

# #Julia formats for saving
using DataFrames,IterableTables
@save "./DATA/Opt_params_schools_R_$(SchoolsOnlyR_0).jld" Θ_Opt CurrBestNegLL
@save "./DATA/Latent_stat_schools_R_$(SchoolsOnlyR_0).jld" P_ModelParams.ξ̄ P_ModelParams.ϕ̄ P_ModelParams.σ_ξ P_ModelParams.σ_ϕ P_ModelParams.ρ_ξϕ
Z_array = [P_ModelParams.Z_ξ[2:16]';P_ModelParams.Z_ϕ[2:16]']
@save "./DATA/Latent_Vars_schools_R_$(SchoolsOnlyR_0).jld" Z_array
@save "./DATA/x_0_for_subsequent_sims_R_$(SchoolsOnlyR_0).jld" x_0_end

# # MATLAB formats
using MAT
t_opt,true_inc,model_pred,plt = create_comparison_plot(sol_rsv)
file = matopen("./DATA/predictions_for_plotting_model_schools_R_$(SchoolsOnlyR_0).mat","w")
write(file,"pred_times",t_opt)
write(file,"true_incidence",true_inc)
write(file,"model_pred",model_pred)
close(file)
