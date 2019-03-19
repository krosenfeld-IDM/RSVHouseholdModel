# Load best parameters from the inference stage
P_ModelParams.b_S = SchoolsContactRate;
P_ModelParams.b_A = 0.
Θ_Opt = load("./DATA/Opt_params_schools_R_$(SchoolsOnlyR_0).jld")["Θ_Opt"]
P_ModelParams.σ_ξ = load("./DATA/Latent_stat_schools_R_$(SchoolsOnlyR_0).jld")["P_ModelParams.σ_ξ"]
P_ModelParams.σ_ϕ = load("./DATA/Latent_stat_schools_R_$(SchoolsOnlyR_0).jld")["P_ModelParams.σ_ϕ"]
P_ModelParams.ρ_ξϕ = load("./DATA/Latent_stat_schools_R_$(SchoolsOnlyR_0).jld")["P_ModelParams.ρ_ξϕ"]
P_ModelParams.ξ̄ = load("./DATA/Latent_stat_schools_R_$(SchoolsOnlyR_0).jld")["P_ModelParams.ξ̄"]
P_ModelParams.Z_ξ[:] = P_ModelParams.ξ̄
P_ModelParams.ϕ̄ = load("./DATA/Latent_stat_schools_R_$(SchoolsOnlyR_0).jld")["P_ModelParams.ϕ̄"]
P_ModelParams.Z_ϕ[:] = P_ModelParams.ϕ̄
z = load("./DATA/Latent_Vars_schools_R_$(SchoolsOnlyR_0).jld")["Z_array"]
P_ModelParams.b_C = Θ_Opt[1]
P_ModelParams.β = Θ_Opt[2]
P_ModelParams.α = 1/Θ_Opt[3]
P_ModelParams.τ = Θ_Opt[4]
P_ModelParams.Z_ξ[2:16] = z[1,:]
P_ModelParams.Z_ϕ[2:16] = z[2,:]
SetMixingMatrix!(P_ModelParams)


#Load initial condition from the optimised process
x_0 = load("DATA/x_0_for_subsequent_sims_R_$(SchoolsOnlyR_0).jld")["x_0_end"]
x_0[1:M_a] = 0.
x_0[find(x_0 .< 0)] = 0.


#Define the multivariate normal distribution that defines the seasonality
SeasonalDistrib = MvNormal([P_ModelParams.ξ̄,P_ModelParams.ϕ̄],[P_ModelParams.σ_ξ P_ModelParams.ρ_ξϕ;P_ModelParams.ρ_ξϕ P_ModelParams.σ_ϕ])
