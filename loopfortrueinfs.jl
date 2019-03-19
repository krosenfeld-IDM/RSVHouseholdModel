# Loop for finding true infecteds
#Loop over vaccination scenarios

@everywhere HHCoverage = 0.
@everywhere generate_replacement_rate!(P_VacModel,100)
sol_repvac = solve(RepeatedForecasts,num_monte = 500,CVODE_BDF(linear_solver=:GMRES),abstol = 1e-3,reltol = 1e-3)
@save "./DATA/HH_inf_baseline_R_$(SchoolsOnlyR_0).jld" sol_repvac.u

# @save "Baseline_incidence.jld" sol_repvac.u
# tot = [sum(sol_repvac.u[i][end]) for i =1:500]

@everywhere CurrentSolidProtectionDuration = 15.
sol_repvac = solve(RepeatedForecasts,num_monte = 500,CVODE_BDF(linear_solver=:GMRES),abstol = 1e-3,reltol = 1e-3)
@save "./DATA/solid_prot_inf_$(CurrentSolidProtectionDuration)_HHCov_$(HHCoverage)_MABCov_$(MABCoverage)__R_$(SchoolsOnlyR_0).jld" sol_repvac.u

@everywhere CurrentSolidProtectionDuration = 30.
sol_repvac = solve(RepeatedForecasts,num_monte = 500,CVODE_BDF(linear_solver=:GMRES),abstol = 1e-3,reltol = 1e-3)
@save "./DATA/solid_prot_inf_$(CurrentSolidProtectionDuration)_HHCov_$(HHCoverage)_MABCov_$(MABCoverage)__R_$(SchoolsOnlyR_0).jld" sol_repvac.u

@everywhere CurrentSolidProtectionDuration = 45.
sol_repvac = solve(RepeatedForecasts,num_monte = 500,CVODE_BDF(linear_solver=:GMRES),abstol = 1e-3,reltol = 1e-3)
@save "./DATA/solid_prot_inf_$(CurrentSolidProtectionDuration)_HHCov_$(HHCoverage)_MABCov_$(MABCoverage)__R_$(SchoolsOnlyR_0).jld" sol_repvac.u

@everywhere CurrentSolidProtectionDuration = 60.
sol_repvac = solve(RepeatedForecasts,num_monte = 500,CVODE_BDF(linear_solver=:GMRES),abstol = 1e-3,reltol = 1e-3)
@save "./DATA/solid_prot_inf_$(CurrentSolidProtectionDuration)_HHCov_$(HHCoverage)_MABCov_$(MABCoverage)__R_$(SchoolsOnlyR_0).jld" sol_repvac.u

@everywhere CurrentSolidProtectionDuration = 75.
sol_repvac = solve(RepeatedForecasts,num_monte = 500,CVODE_BDF(linear_solver=:GMRES),abstol = 1e-3,reltol = 1e-3)
@save "./DATA/solid_prot_inf_$(CurrentSolidProtectionDuration)_HHCov_$(HHCoverage)_MABCov_$(MABCoverage)__R_$(SchoolsOnlyR_0).jld" sol_repvac.u

@everywhere CurrentSolidProtectionDuration = 90.
sol_repvac = solve(RepeatedForecasts,num_monte = 500,CVODE_BDF(linear_solver=:GMRES),abstol = 1e-3,reltol = 1e-3)
@save "./DATA/solid_prot_inf_$(CurrentSolidProtectionDuration)_HHCov_$(HHCoverage)_MABCov_$(MABCoverage)__R_$(SchoolsOnlyR_0).jld" sol_repvac.u

# @everywhere CurrentSolidProtectionDuration = 120.
# sol_repvac = solve(RepeatedForecasts,num_monte = 500,CVODE_BDF(linear_solver=:GMRES),abstol = 1e-3,reltol = 1e-3)
# @save "solid_prot_inf_$(CurrentSolidProtectionDuration)_HHCov_$(HHCoverage)_MABCov_$(MABCoverage)__R_$(SchoolsOnlyR_0).jld" sol_repvac.u
# gr()
# BaseLine = load("Baseline_incidence.jld")["sol_repvac.u"]
# solid_30 = load("solid_prot_inf_30.0_vs2.jld")["sol_repvac.u"]
#
# BaseMean = mean([sum(BaseLine[i][end]) for i = 1:500])
# HHVacMean = mean([sum(sol_repvac.u[i][end]) for i = 1:10])
# title!("Forecast RSV Hospitalisation: Mean = $(BaseMean)")
#
# solid_30 = load("solid_prot_inf_30.0_HHCov_$(HHCoverage)_MABCov_$(MABCoverage)__R_$(SchoolsOnlyR_0).jld")["sol_repvac.u"]
# t,x,p = plot_forecast_solution(solid_30[1:20],30.);display(p)
# solid_30_mean = mean([sum(solid_30[i][end]) for i = 1:500])
# title!("30. days solid prot.: Mean imp. = $(BaseMean - solid_30_mean)")
#
# solid_60 = load("solid_prot_inf_60.0_vs2.jld")["sol_repvac.u"]
# t,x,p = plot_forecast_solution(solid_60[1:20],30.);display(p)
# solid_60_mean = mean([sum(solid_60[i][end]) for i = 1:500])
# title!("60. days solid prot.: Mean imp. = $(100*(BaseMean - solid_60_mean)/BaseMean)%")
#
# solid_90 = load("./DATA/solid_prot_inf_90.0_HHCov_$(HHCoverage)_MABCov_$(MABCoverage)__R_$(SchoolsOnlyR_0).jld")["sol_repvac.u"]
# t,x,p = plot_forecast_solution(solid_90[1:20],30.);display(p)
# solid_90_mean = mean([sum(solid_90[i][end]) for i = 1:500])
# title!("90. days solid prot.: Mean imp. = $(100*(BaseMean - solid_90_mean)/BaseMean)%")
#
# solid_120 = load("./DATA/solid_prot_inf_120.0_vs2.jld")["sol_repvac.u"]
# t,x,p = plot_forecast_solution(solid_120[1:20],30.);display(p)
# solid_120_mean = mean([sum(solid_120[i][end]) for i = 1:500])
# title!("120. days solid prot.: Mean imp. = $(100*(BaseMean - solid_120_mean)/BaseMean)%")
#
# # solid_120 = load("solid_prot_inf_120.0_vs2.jld")["sol_repvac.u"]
# t,x,p = plot_forecast_solution(sol_repvac.u[1:10],30.);display(p)
# hh_inf_mean = mean([sum(sol_repvac.u[i][end]) for i = 1:10])
# title!("HH vaccination: Mean imp. = $(100*(BaseMean - hh_inf_mean)/BaseMean)%")
#

@everywhere FittedMeanDuration = 1/P_ModelParams.α
@everywhere CurrentSolidProtectionDuration = 0.
@everywhere HHCoverage = 0.
@everywhere generate_replacement_rate!(P_VacModel,100)

@everywhere P_ModelParams.α = 1/(FittedMeanDuration + 15.)
sol_repvac = solve(RepeatedForecasts,num_monte = 500,CVODE_BDF(linear_solver=:GMRES),abstol = 1e-3,reltol = 1e-3)
@save "./DATA/mean_prot_inf_$((1/P_ModelParams.α) - FittedMeanDuration)_HHCov_$(HHCoverage)_MABCov_$(MABCoverage)__R_$(SchoolsOnlyR_0).jld" sol_repvac.u

@everywhere P_ModelParams.α = 1/(FittedMeanDuration + 30.)
sol_repvac = solve(RepeatedForecasts,num_monte = 500,CVODE_BDF(linear_solver=:GMRES),abstol = 1e-3,reltol = 1e-3)
@save "./DATA/mean_prot_inf_$((1/P_ModelParams.α) - FittedMeanDuration)_HHCov_$(HHCoverage)_MABCov_$(MABCoverage)__R_$(SchoolsOnlyR_0).jld" sol_repvac.u

@everywhere P_ModelParams.α = 1/(FittedMeanDuration + 45.)
sol_repvac = solve(RepeatedForecasts,num_monte = 500,CVODE_BDF(linear_solver=:GMRES),abstol = 1e-3,reltol = 1e-3)
@save "./DATA/mean_prot_inf_$((1/P_ModelParams.α) - FittedMeanDuration)_HHCov_$(HHCoverage)_MABCov_$(MABCoverage)__R_$(SchoolsOnlyR_0).jld" sol_repvac.u

@everywhere P_ModelParams.α = 1/(FittedMeanDuration + 60.)
sol_repvac = solve(RepeatedForecasts,num_monte = 500,CVODE_BDF(linear_solver=:GMRES),abstol = 1e-3,reltol = 1e-3)
@save "./DATA/mean_prot_inf_$((1/P_ModelParams.α) - FittedMeanDuration)_HHCov_$(HHCoverage)_MABCov_$(MABCoverage)__R_$(SchoolsOnlyR_0).jld" sol_repvac.u

@everywhere P_ModelParams.α = 1/(FittedMeanDuration + 75.)
sol_repvac = solve(RepeatedForecasts,num_monte = 500,CVODE_BDF(linear_solver=:GMRES),abstol = 1e-3,reltol = 1e-3)
@save "./DATA/mean_prot_inf_$((1/P_ModelParams.α) - FittedMeanDuration)_HHCov_$(HHCoverage)_MABCov_$(MABCoverage)__R_$(SchoolsOnlyR_0).jld" sol_repvac.u

@everywhere P_ModelParams.α = 1/(FittedMeanDuration + 90.)
sol_repvac = solve(RepeatedForecasts,num_monte = 500,CVODE_BDF(linear_solver=:GMRES),abstol = 1e-3,reltol = 1e-3)
@save "./DATA/mean_prot_inf_$((1/P_ModelParams.α) - FittedMeanDuration)_HHCov_$(HHCoverage)_MABCov_$(MABCoverage)__R_$(SchoolsOnlyR_0).jld" sol_repvac.u
# t,x,p = plot_forecast_solution(sol_repvac.u[1:20],30.);display(p)

@everywhere HHCoverage = 0.25
@everywhere generate_replacement_rate!(P_VacModel,100)
@everywhere P_ModelParams.α = 1/(FittedMeanDuration)
@everywhere CurrentSolidProtectionDuration = 0.
sol_repvac = solve(RepeatedForecasts,num_monte = 500,CVODE_BDF(linear_solver=:GMRES),abstol = 1e-3,reltol = 1e-3)
@save "./DATA/solid_prot_inf_$(CurrentSolidProtectionDuration)_HHCov_$(HHCoverage)_MABCov_$(MABCoverage)__R_$(SchoolsOnlyR_0).jld" sol_repvac.u

@everywhere HHCoverage = 0.5
@everywhere generate_replacement_rate!(P_VacModel,100)
@everywhere P_ModelParams.α = 1/(FittedMeanDuration)
@everywhere CurrentSolidProtectionDuration = 0.
sol_repvac = solve(RepeatedForecasts,num_monte = 500,CVODE_BDF(linear_solver=:GMRES),abstol = 1e-3,reltol = 1e-3)
@save "./DATA/solid_prot_inf_$(CurrentSolidProtectionDuration)_HHCov_$(HHCoverage)_MABCov_$(MABCoverage)__R_$(SchoolsOnlyR_0).jld" sol_repvac.u

@everywhere HHCoverage = 0.75
@everywhere generate_replacement_rate!(P_VacModel,100)
@everywhere P_ModelParams.α = 1/(FittedMeanDuration)
@everywhere CurrentSolidProtectionDuration = 0.
sol_repvac = solve(RepeatedForecasts,num_monte = 500,CVODE_BDF(linear_solver=:GMRES),abstol = 1e-3,reltol = 1e-3)
@save "./DATA/solid_prot_inf_$(CurrentSolidProtectionDuration)_HHCov_$(HHCoverage)_MABCov_$(MABCoverage)__R_$(SchoolsOnlyR_0).jld" sol_repvac.u

@everywhere HHCoverage = 1.
@everywhere generate_replacement_rate!(P_VacModel,100)
@everywhere P_ModelParams.α = 1/(FittedMeanDuration)
@everywhere CurrentSolidProtectionDuration = 0.
sol_repvac = solve(RepeatedForecasts,num_monte = 500,CVODE_BDF(linear_solver=:GMRES),abstol = 1e-3,reltol = 1e-3)
@save "./DATA/solid_prot_inf_$(CurrentSolidProtectionDuration)_HHCov_$(HHCoverage)_MABCov_$(MABCoverage)__R_$(SchoolsOnlyR_0).jld" sol_repvac.u

#------------------

@everywhere HHCoverage = 0.25
@everywhere generate_replacement_rate!(P_VacModel,100)
@everywhere P_ModelParams.α = 1/(FittedMeanDuration)
@everywhere CurrentSolidProtectionDuration = 15.
sol_repvac = solve(RepeatedForecasts,num_monte = 500,CVODE_BDF(linear_solver=:GMRES),abstol = 1e-3,reltol = 1e-3)
@save "./DATA/solid_prot_inf_$(CurrentSolidProtectionDuration)_HHCov_$(HHCoverage)_MABCov_$(MABCoverage)__R_$(SchoolsOnlyR_0).jld" sol_repvac.u

@everywhere HHCoverage = 0.5
@everywhere generate_replacement_rate!(P_VacModel,100)
@everywhere P_ModelParams.α = 1/(FittedMeanDuration)
@everywhere CurrentSolidProtectionDuration = 15.
sol_repvac = solve(RepeatedForecasts,num_monte = 500,CVODE_BDF(linear_solver=:GMRES),abstol = 1e-3,reltol = 1e-3)
@save "./DATA/solid_prot_inf_$(CurrentSolidProtectionDuration)_HHCov_$(HHCoverage)_MABCov_$(MABCoverage)__R_$(SchoolsOnlyR_0).jld" sol_repvac.u

@everywhere HHCoverage = 0.75
@everywhere generate_replacement_rate!(P_VacModel,100)
@everywhere P_ModelParams.α = 1/(FittedMeanDuration)
@everywhere CurrentSolidProtectionDuration = 15.
sol_repvac = solve(RepeatedForecasts,num_monte = 500,CVODE_BDF(linear_solver=:GMRES),abstol = 1e-3,reltol = 1e-3)
@save "./DATA/solid_prot_inf_$(CurrentSolidProtectionDuration)_HHCov_$(HHCoverage)_MABCov_$(MABCoverage)__R_$(SchoolsOnlyR_0).jld" sol_repvac.u

@everywhere HHCoverage = 1.
@everywhere generate_replacement_rate!(P_VacModel,100)
@everywhere P_ModelParams.α = 1/(FittedMeanDuration)
@everywhere CurrentSolidProtectionDuration = 15.
sol_repvac = solve(RepeatedForecasts,num_monte = 500,CVODE_BDF(linear_solver=:GMRES),abstol = 1e-3,reltol = 1e-3)
@save "./DATA/solid_prot_inf_$(CurrentSolidProtectionDuration)_HHCov_$(HHCoverage)_MABCov_$(MABCoverage)__R_$(SchoolsOnlyR_0).jld" sol_repvac.u


#------------------
@everywhere CurrentSolidProtectionDuration = 30.


@everywhere HHCoverage = 0.25
@everywhere generate_replacement_rate!(P_VacModel,100)
@everywhere P_ModelParams.α = 1/(FittedMeanDuration)
sol_repvac = solve(RepeatedForecasts,num_monte = 500,CVODE_BDF(linear_solver=:GMRES),abstol = 1e-3,reltol = 1e-3)
@save "./DATA/solid_prot_inf_$(CurrentSolidProtectionDuration)_HHCov_$(HHCoverage)_MABCov_$(MABCoverage)__R_$(SchoolsOnlyR_0).jld" sol_repvac.u

@everywhere HHCoverage = 0.5
@everywhere generate_replacement_rate!(P_VacModel,100)
@everywhere P_ModelParams.α = 1/(FittedMeanDuration)
sol_repvac = solve(RepeatedForecasts,num_monte = 500,CVODE_BDF(linear_solver=:GMRES),abstol = 1e-3,reltol = 1e-3)
@save "./DATA/solid_prot_inf_$(CurrentSolidProtectionDuration)_HHCov_$(HHCoverage)_MABCov_$(MABCoverage)__R_$(SchoolsOnlyR_0).jld" sol_repvac.u

@everywhere HHCoverage = 0.75
@everywhere generate_replacement_rate!(P_VacModel,100)
@everywhere P_ModelParams.α = 1/(FittedMeanDuration)
sol_repvac = solve(RepeatedForecasts,num_monte = 500,CVODE_BDF(linear_solver=:GMRES),abstol = 1e-3,reltol = 1e-3)
@save "./DATA/solid_prot_inf_$(CurrentSolidProtectionDuration)_HHCov_$(HHCoverage)_MABCov_$(MABCoverage)__R_$(SchoolsOnlyR_0).jld" sol_repvac.u

@everywhere HHCoverage = 1.
@everywhere generate_replacement_rate!(P_VacModel,100)
@everywhere P_ModelParams.α = 1/(FittedMeanDuration)
sol_repvac = solve(RepeatedForecasts,num_monte = 500,CVODE_BDF(linear_solver=:GMRES),abstol = 1e-3,reltol = 1e-3)
@save "./DATA/solid_prot_inf_$(CurrentSolidProtectionDuration)_HHCov_$(HHCoverage)_MABCov_$(MABCoverage)__R_$(SchoolsOnlyR_0).jld" sol_repvac.u


#------------------
@everywhere CurrentSolidProtectionDuration = 45.


@everywhere HHCoverage = 0.25
@everywhere generate_replacement_rate!(P_VacModel,100)
@everywhere P_ModelParams.α = 1/(FittedMeanDuration)
sol_repvac = solve(RepeatedForecasts,num_monte = 500,CVODE_BDF(linear_solver=:GMRES),abstol = 1e-3,reltol = 1e-3)
@save "./DATA/solid_prot_inf_$(CurrentSolidProtectionDuration)_HHCov_$(HHCoverage)_MABCov_$(MABCoverage)__R_$(SchoolsOnlyR_0).jld" sol_repvac.u

@everywhere HHCoverage = 0.5
@everywhere generate_replacement_rate!(P_VacModel,100)
@everywhere P_ModelParams.α = 1/(FittedMeanDuration)
sol_repvac = solve(RepeatedForecasts,num_monte = 500,CVODE_BDF(linear_solver=:GMRES),abstol = 1e-3,reltol = 1e-3)
@save "./DATA/solid_prot_inf_$(CurrentSolidProtectionDuration)_HHCov_$(HHCoverage)_MABCov_$(MABCoverage)__R_$(SchoolsOnlyR_0).jld" sol_repvac.u

@everywhere HHCoverage = 0.75
@everywhere generate_replacement_rate!(P_VacModel,100)
@everywhere P_ModelParams.α = 1/(FittedMeanDuration)
sol_repvac = solve(RepeatedForecasts,num_monte = 500,CVODE_BDF(linear_solver=:GMRES),abstol = 1e-3,reltol = 1e-3)
@save "./DATA/solid_prot_inf_$(CurrentSolidProtectionDuration)_HHCov_$(HHCoverage)_MABCov_$(MABCoverage)__R_$(SchoolsOnlyR_0).jld" sol_repvac.u

@everywhere HHCoverage = 1.
@everywhere generate_replacement_rate!(P_VacModel,100)
@everywhere P_ModelParams.α = 1/(FittedMeanDuration)
sol_repvac = solve(RepeatedForecasts,num_monte = 500,CVODE_BDF(linear_solver=:GMRES),abstol = 1e-3,reltol = 1e-3)
@save "./DATA/solid_prot_inf_$(CurrentSolidProtectionDuration)_HHCov_$(HHCoverage)_MABCov_$(MABCoverage)__R_$(SchoolsOnlyR_0).jld" sol_repvac.u

#------------------
@everywhere CurrentSolidProtectionDuration = 60.


@everywhere HHCoverage = 0.25
@everywhere generate_replacement_rate!(P_VacModel,100)
@everywhere P_ModelParams.α = 1/(FittedMeanDuration)
sol_repvac = solve(RepeatedForecasts,num_monte = 500,CVODE_BDF(linear_solver=:GMRES),abstol = 1e-3,reltol = 1e-3)
@save "./DATA/solid_prot_inf_$(CurrentSolidProtectionDuration)_HHCov_$(HHCoverage)_MABCov_$(MABCoverage)__R_$(SchoolsOnlyR_0).jld" sol_repvac.u

@everywhere HHCoverage = 0.5
@everywhere generate_replacement_rate!(P_VacModel,100)
@everywhere P_ModelParams.α = 1/(FittedMeanDuration)
sol_repvac = solve(RepeatedForecasts,num_monte = 500,CVODE_BDF(linear_solver=:GMRES),abstol = 1e-3,reltol = 1e-3)
@save "./DATA/solid_prot_inf_$(CurrentSolidProtectionDuration)_HHCov_$(HHCoverage)_MABCov_$(MABCoverage)__R_$(SchoolsOnlyR_0).jld" sol_repvac.u

@everywhere HHCoverage = 0.75
@everywhere generate_replacement_rate!(P_VacModel,100)
@everywhere P_ModelParams.α = 1/(FittedMeanDuration)
sol_repvac = solve(RepeatedForecasts,num_monte = 500,CVODE_BDF(linear_solver=:GMRES),abstol = 1e-3,reltol = 1e-3)
@save "./DATA/solid_prot_inf_$(CurrentSolidProtectionDuration)_HHCov_$(HHCoverage)_MABCov_$(MABCoverage)__R_$(SchoolsOnlyR_0).jld" sol_repvac.u

@everywhere HHCoverage = 1.
@everywhere generate_replacement_rate!(P_VacModel,100)
@everywhere P_ModelParams.α = 1/(FittedMeanDuration)
sol_repvac = solve(RepeatedForecasts,num_monte = 500,CVODE_BDF(linear_solver=:GMRES),abstol = 1e-3,reltol = 1e-3)
@save "./DATA/solid_prot_inf_$(CurrentSolidProtectionDuration)_HHCov_$(HHCoverage)_MABCov_$(MABCoverage)__R_$(SchoolsOnlyR_0).jld" sol_repvac.u


#------------------
@everywhere CurrentSolidProtectionDuration = 75.


@everywhere HHCoverage = 0.25
@everywhere generate_replacement_rate!(P_VacModel,100)
@everywhere P_ModelParams.α = 1/(FittedMeanDuration)
sol_repvac = solve(RepeatedForecasts,num_monte = 500,CVODE_BDF(linear_solver=:GMRES),abstol = 1e-3,reltol = 1e-3)
@save "./DATA/solid_prot_inf_$(CurrentSolidProtectionDuration)_HHCov_$(HHCoverage)_MABCov_$(MABCoverage)__R_$(SchoolsOnlyR_0).jld" sol_repvac.u

@everywhere HHCoverage = 0.5
@everywhere generate_replacement_rate!(P_VacModel,100)
@everywhere P_ModelParams.α = 1/(FittedMeanDuration)
sol_repvac = solve(RepeatedForecasts,num_monte = 500,CVODE_BDF(linear_solver=:GMRES),abstol = 1e-3,reltol = 1e-3)
@save "./DATA/solid_prot_inf_$(CurrentSolidProtectionDuration)_HHCov_$(HHCoverage)_MABCov_$(MABCoverage)__R_$(SchoolsOnlyR_0).jld" sol_repvac.u

@everywhere HHCoverage = 0.75
@everywhere generate_replacement_rate!(P_VacModel,100)
@everywhere P_ModelParams.α = 1/(FittedMeanDuration)
sol_repvac = solve(RepeatedForecasts,num_monte = 500,CVODE_BDF(linear_solver=:GMRES),abstol = 1e-3,reltol = 1e-3)
@save "./DATA/solid_prot_inf_$(CurrentSolidProtectionDuration)_HHCov_$(HHCoverage)_MABCov_$(MABCoverage)__R_$(SchoolsOnlyR_0).jld" sol_repvac.u

@everywhere HHCoverage = 1.
@everywhere generate_replacement_rate!(P_VacModel,100)
@everywhere P_ModelParams.α = 1/(FittedMeanDuration)
sol_repvac = solve(RepeatedForecasts,num_monte = 500,CVODE_BDF(linear_solver=:GMRES),abstol = 1e-3,reltol = 1e-3)
@save "./DATA/solid_prot_inf_$(CurrentSolidProtectionDuration)_HHCov_$(HHCoverage)_MABCov_$(MABCoverage)__R_$(SchoolsOnlyR_0).jld" sol_repvac.u


#------------------
@everywhere CurrentSolidProtectionDuration = 90.


@everywhere HHCoverage = 0.25
@everywhere generate_replacement_rate!(P_VacModel,100)
@everywhere P_ModelParams.α = 1/(FittedMeanDuration)
sol_repvac = solve(RepeatedForecasts,num_monte = 500,CVODE_BDF(linear_solver=:GMRES),abstol = 1e-3,reltol = 1e-3)
@save "./DATA/solid_prot_inf_$(CurrentSolidProtectionDuration)_HHCov_$(HHCoverage)_MABCov_$(MABCoverage)__R_$(SchoolsOnlyR_0).jld" sol_repvac.u

@everywhere HHCoverage = 0.5
@everywhere generate_replacement_rate!(P_VacModel,100)
@everywhere P_ModelParams.α = 1/(FittedMeanDuration)
sol_repvac = solve(RepeatedForecasts,num_monte = 500,CVODE_BDF(linear_solver=:GMRES),abstol = 1e-3,reltol = 1e-3)
@save "./DATA/solid_prot_inf_$(CurrentSolidProtectionDuration)_HHCov_$(HHCoverage)_MABCov_$(MABCoverage)__R_$(SchoolsOnlyR_0).jld" sol_repvac.u

@everywhere HHCoverage = 0.75
@everywhere generate_replacement_rate!(P_VacModel,100)
@everywhere P_ModelParams.α = 1/(FittedMeanDuration)
sol_repvac = solve(RepeatedForecasts,num_monte = 500,CVODE_BDF(linear_solver=:GMRES),abstol = 1e-3,reltol = 1e-3)
@save "./DATA/solid_prot_inf_$(CurrentSolidProtectionDuration)_HHCov_$(HHCoverage)_MABCov_$(MABCoverage)__R_$(SchoolsOnlyR_0).jld" sol_repvac.u

@everywhere HHCoverage = 1.
@everywhere generate_replacement_rate!(P_VacModel,100)
@everywhere P_ModelParams.α = 1/(FittedMeanDuration)
sol_repvac = solve(RepeatedForecasts,num_monte = 500,CVODE_BDF(linear_solver=:GMRES),abstol = 1e-3,reltol = 1e-3)
@save "./DATA/solid_prot_inf_$(CurrentSolidProtectionDuration)_HHCov_$(HHCoverage)_MABCov_$(MABCoverage)__R_$(SchoolsOnlyR_0).jld" sol_repvac.u
