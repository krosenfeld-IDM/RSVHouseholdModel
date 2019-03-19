# Objective functions for minimisation --- these are the likelihoods for both the E and M steps

# Set up objective functions
function notinbounds(P::HH_RSVModelParameters) #Checks the positivity of the appropriate parameters
    return any((P.b_C <0,P.b_S<0,P.b_A<0,P.τ<0,P.α<=0.,P.ϵ <0,P.EffHHSizePower<0 || P.EffHHSizePower>1))
end

function isinpredtimes(times,prediction_times)
    return [in(times[i],prediction_times) for i = 1:length(times)]
end

function relevant_pred_times_from_sol(sol,prediction_times)
    return prediction_times[(prediction_times .>= sol.t[1] ).*(prediction_times .<= sol.t[end])]
end

function relevant_pred_times_from_time_interval(time_span,prediction_times)
    return prediction_times[(prediction_times .>= time_span[1] ).*(prediction_times .<= time_span[end])]
end

function find_model_incidence(sol,R,prediction_times)
    return [diff(sol[i,isinpredtimes(sol.t,prediction_times)]).*R[2:end] for i = 1:16]
end

function find_KDHSS_ratio(prediction_times)
    return [(Ratio(t) + 1.)/Ratio(t) for t in prediction_times]#Include people from outside KDHSS
end

function find_model_incidence_using_interpolation(sol,prediction_times)
    pred_times_in_sol = relevant_pred_times_from_sol(sol,prediction_times)
    R = find_KDHSS_ratio(pred_times_in_sol)
    return [diff([sol(t)[i] for t in pred_times_in_sol]).*R[2:end] for i = 1:16]
end

function find_model_incidence_using_interpolation_age_cat(sol,prediction_times,age_cat)
    pred_times_in_sol = relevant_pred_times_from_sol(sol,prediction_times)
    R = find_KDHSS_ratio(pred_times_in_sol)
    return diff([sol(t)[age_cat] for t in pred_times_in_sol]).*R[2:end]
end

function ismodelprednegative(ModelPred)
    return any((any(v -> v .< -0.001,ModelPred[i]) for i = 1:length(ModelPred)))
end

function make_non_neg!(X)
    for i = 1:length(X)
        for j = 1:length(X[i])
            if X[i][j] <0.
                X[i][j] = 0.
            end
        end
    end
    return nothing
end

function get_trueincidencedata(sol,prediction_times)
    time_bool = (prediction_times .>= sol.t[1]).&(prediction_times .<= sol.t[end])
    return [Agg_RSVData[i][time_bool][2:end] for i = 1:16] #First element left off because of diff operation
end

function get_trueincidencedata_age_cat(sol,prediction_times,age_cat::Int64)
    time_bool = (prediction_times .>= sol.t[1]).&(prediction_times .<= sol.t[end])
    return Agg_RSVData[age_cat][time_bool][2:end] #First element left off because of diff operation
end

function neg_LL_for_RSV_data_HHModel(sol)
   tot_loss = 0.0
   if any((s.retcode != :Success for s in sol)) ||  notinbounds(sol.prob.p)
     tot_loss = Inf
   else
      ModelIncidence = find_model_incidence_using_interpolation(sol,pred_times)
      if ismodelprednegative(ModelIncidence)
          tot_loss = Inf
      else
          make_non_neg!(ModelIncidence)
          D = get_trueincidencedata(sol,pred_times)
          l = length(D[1])
          for n = 1:16
              for i = 1:l
                  tot_loss -= logpdf(Poisson(ModelIncidence[n][i]+1e-5),D[n][i])
              end
          end
      end
   end
   return tot_loss
end


#E step functions -- for each year using a specified intial condition

function Neg_objective_for_latent_vars_one_season(Z,season_num,new_x_0)
    if Z[1] < 0 || ~(Z[2] > -180. && Z[2] < 180)
        return Inf
    else
        # Set the season's latent variables
        P_ModelParams.Z_ξ[season_num] = Z[1]
        P_ModelParams.Z_ϕ[season_num] = Z[2]
        t_0 = SeasonStartTimes[season_num]
        t_1 = SeasonEndTimes[season_num]
        tspan = (t_0,t_1)
        times_to_saveat = relevant_pred_times_from_time_interval(tspan,pred_times)
        prob_season = ODEProblem(F,new_x_0,tspan,P_ModelParams)
        sol_one_season = solve(prob_season,CVODE_BDF(linear_solver=:GMRES),abstol = 1e-3, reltol = 1e-3,saveat = times_to_saveat,callback = cb_set_annual_changes)        #Return log joint probability
        m = [P_ModelParams.ξ̄,P_ModelParams.ϕ̄]
        Σ = [P_ModelParams.σ_ξ P_ModelParams.ρ_ξϕ;P_ModelParams.ρ_ξϕ P_ModelParams.σ_ϕ]
        return neg_LL_for_RSV_data_HHModel(sol_one_season) - logpdf(MvNormal(m, Σ),Z), sol_one_season
    end
end

#M step functions




function neg_objective_for_params(Θ::Vector{Float64})
    if any(Θ .< 0)
        return Inf
    else
        #Preliminary solution to find an initial condition for 01-09-2001 (beginning of season 2)
        #Adjust the variable ParameterizedFunction
        P_ModelParams.b_C = Θ[1];
        P_ModelParams.b_S = SchoolsContactRate;
        P_ModelParams.b_A = 0.#Θ[2]
        P_ModelParams.β = Θ[2]
        P_ModelParams.β_S = 0.
        P_ModelParams.β_C = 0.
        SetMixingMatrix!(P_ModelParams) #This updates with within population parameter changes
        P_ModelParams.α = 1/Θ[3] #Convert from rate to mean
        P_ModelParams.τ = Θ[4]
        P_ModelParams.ϵ = 10.
        P_ModelParams.Z_ξ[1] = P_ModelParams.ξ̄
        P_ModelParams.Z_ϕ[1] = P_ModelParams.ϕ̄

        x_0_eq = get_equilibrium_state(SeasonStartDates[2])[1]
        t_0 = SeasonStartTimes[2]
        t_1 = SeasonEndTimes[16]
        tspan_rsv = (t_0,t_1)
        prob_rsv = ODEProblem(F,x_0_eq,tspan_rsv,P_ModelParams);
        sol_rsv = solve(prob_rsv,CVODE_BDF(linear_solver=:GMRES),abstol = 1e-3, reltol = 1e-3,callback = cb_set_annual_changes)
        return neg_LL_for_RSV_data_HHModel(sol_rsv),sol_rsv
    end
end



function neg_objective_for_params_initial(Θ_initial::Vector{Float64})
    if any(Θ_initial[1:5] .< 0) #Only the phase can be negative
        return Inf
    else
        #Preliminary solution to find an initial condition for 01-09-2001 (beginning of season 2)
        #Adjust the variable ParameterizedFunction
        P_ModelParams.b_C = Θ_initial[1];
        P_ModelParams.b_S = SchoolsContactRate;
        P_ModelParams.b_A = 0.#Θ_initial[2]
        P_ModelParams.β = Θ_initial[2]
        P_ModelParams.β_S = 0.
        P_ModelParams.β_C = 0.
        SetMixingMatrix!(P_ModelParams) #This updates with within population parameter changes
        P_ModelParams.α = 1/Θ_initial[3] #Convert from rate to mean
        P_ModelParams.τ = Θ_initial[4]
        P_ModelParams.ϵ = 10.
        P_ModelParams.ξ̄ = Θ_initial[5]
        P_ModelParams.ϕ̄ = Θ_initial[6]
        P_ModelParams.Z_ξ = P_ModelParams.ξ̄*ones(18)
        P_ModelParams.Z_ϕ = P_ModelParams.ϕ̄*ones(18)
        # set_matrices_for_year!(P_ModelParams,-1)
        t_0 = Float64((SeasonStartDates[5] - Date(2000,1,1)).value)
        t_1 = Float64((SeasonEndDates[9] - Date(2000,1,1)).value)
        tspan_init = (t_0,t_1)
        times_to_saveat = relevant_pred_times_from_time_interval(tspan_init,pred_times)

        #Get equilibrium and draw initial state from that
        x_0_new,sol_eq = get_equilibrium_state(SeasonStartDates[5]) #Gets an equilibrium solution for beginning of 04-05 RSV season
        prob_init = ODEProblem(F,x_0_new,tspan_init,P_ModelParams);
        sol_init = solve(prob_init,CVODE_BDF(linear_solver=:GMRES),callback = cb_set_annual_changes,saveat = times_to_saveat,abstol = 1e-3, reltol = 1e-3)
        #Now fit to the regular bit of the data
        return neg_LL_for_RSV_data_HHModel(sol_init)
    end
end

function sol_for_params_initial(Θ_initial::Vector{Float64})
    if any(Θ_initial[1:5] .< 0) #Only the phase can be negative
        return Inf
    else
        x_0 = [zeros(M_a);construct_susceptible_HH_distribution_from_year(N_HEachYear[5],d1)];
        neg_LL = 0.;
        #Preliminary solution to find an initial condition for 01-09-2001 (beginning of season 2)
        #Adjust the variable ParameterizedFunction
        P_ModelParams.b_C = Θ_initial[1];
        P_ModelParams.b_S = SchoolsContactRate;
        P_ModelParams.b_A = 0.#Θ_initial[2]
        P_ModelParams.β = Θ_initial[2]
        P_ModelParams.β_S = 0.
        P_ModelParams.β_C = 0.
        SetMixingMatrix!(P_ModelParams) #This updates with within population parameter changes
        P_ModelParams.α = 1/Θ_initial[3] #Convert from rate to mean
        P_ModelParams.τ = Θ_initial[4]
        P_ModelParams.ϵ = 10.
        P_ModelParams.ξ̄ = Θ_initial[5]
        P_ModelParams.ϕ̄ = Θ_initial[6]
        P_ModelParams.Z_ξ = P_ModelParams.ξ̄*ones(18)
        P_ModelParams.Z_ϕ = P_ModelParams.ϕ̄*ones(18)
        # set_matrices_for_year!(P_ModelParams,-1)
        t_0 = Float64((SeasonStartDates[5] - Date(2000,1,1)).value)
        t_1 = Float64((SeasonEndDates[9] - Date(2000,1,1)).value)
        tspan_init = (t_0,t_1)
        times_to_saveat = relevant_pred_times_from_time_interval(tspan_init,pred_times)

        #Get equilibrium and draw initial state from that
        x_0_new,sol_eq = get_equilibrium_state(SeasonStartDates[5]) #Gets an equilibrium solution for beginning of 04-05 RSV season
        prob_init = ODEProblem(F,x_0_new,tspan_init,P_ModelParams);
        sol_init = solve(prob_init,CVODE_BDF(linear_solver=:GMRES),callback = cb_set_annual_changes,saveat = times_to_saveat,abstol = 1e-3, reltol = 1e-3)
        #Now fit to the regular bit of the data
        return sol_init,sol_eq
    end
end





function sol_for_latent_vars_one_season(Z,season_num,new_x_0)
    if Z[1] < 0 || ~(Z[2] > -180. && Z[2] < 180)
        return Inf
    else
        # Set the season's latent variables
        P_ModelParams.Z_ξ[season_num] = Z[1]
        P_ModelParams.Z_ϕ[season_num] = Z[2]
        t_0 = SeasonStartTimes[season_num]
        t_1 = SeasonEndTimes[season_num]
        tspan = (t_0,t_1)
        times_to_saveat = relevant_pred_times_from_time_interval(tspan,pred_times)
        prob_season = ODEProblem(F,new_x_0,tspan,P_ModelParams)
        sol_one_season = solve(prob_season,CVODE_BDF(linear_solver=:GMRES),abstol = 1e-3, reltol = 1e-3,saveat = times_to_saveat,callback = cb_set_annual_changes)        #Return log joint probability
        m = [P_ModelParams.ξ̄,P_ModelParams.ϕ̄]
        Σ = [P_ModelParams.σ_ξ P_ModelParams.ρ_ξϕ;P_ModelParams.ρ_ξϕ P_ModelParams.σ_ϕ]
        return  sol_one_season
    end
end
