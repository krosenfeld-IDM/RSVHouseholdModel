# EM algorithms
using Optim,BlackBoxOptim

# Define M step --- Optimise parameters
SearchRangeForInitialParams = [(0.,1.),(0.,1.),(10.,90.),(0.,1.),(0.,1.),(-30.,30.)]
SearchRangeForLatentVars = [(0.,20.),(-90.,90.)]
SearchRangeForModelParams = [(0.,.5),(0.,0.5),(10.,90.),(0.,0.3)]

function perform_initial_M_step_det_optim(Θ_initial::Vector{Float64},n_steps,OptimSolver)
    #Optimise the negative log likelihood without variation in the seasonality
    result = optimize(neg_objective_for_params_initial, Θ_initial, OptimSolver,Optim.Options(iterations = n_steps,show_trace = true))
    Θ_intial_Opt = Optim.minimizer(result)
    Neg_LL_opt = Optim.minimum(result)
    return Θ_intial_Opt,Neg_LL_opt
end

function perform_initial_M_step_det_bb(BBOptimSolver,maxfuncs,iters_nm)
    println("Starting initial optimisation for parameters on limited set")
    opt1 = bbsetup(neg_objective_for_params_initial;
        Method=BBOptimSolver,
        SearchRange = SearchRangeForInitialParams,
        MaxFuncEvals = maxfuncs,
        PopulationSize = 100,
        TraceMode = :verbose)
    res = bboptimize(opt1)
    Θ_init_opt = best_candidate(res)
    println("Constructing simplex for NelderMead algorithm")
    result = optimize(neg_objective_for_params_initial, Θ_init_opt, NelderMead(),Optim.Options(iterations = iters_nm,show_trace = true))
    Θ_init_opt = Optim.minimizer(result)
    # Neg_LL_opt_init = Optim.minimum(result)
    P_ModelParams.σ_ξ = 100.
    P_ModelParams.σ_ϕ = 10000.
    P_ModelParams.ρ_ξϕ = 0.
    return Θ_init_opt
end


Perform_E_step = function(iter_steps_nm)
    #1) Equilibriate past data
    P_ModelParams.Z_ξ[1] = P_ModelParams.ξ̄
    P_ModelParams.Z_ϕ[1] = P_ModelParams.ϕ̄
    new_x_0 = get_equilibrium_state(SeasonStartDates[2])[1]
    # Map to an objective function for subsequent seasons
    for num = 2:16
        Obj(Z) = Neg_objective_for_latent_vars_one_season(Z,num,new_x_0)[1]
        println("Optimising for seasonality on season $(num)")
        result = optimize(Obj, [P_ModelParams.ξ̄,P_ModelParams.ϕ̄], NelderMead(),Optim.Options(iterations = iter_steps_nm))
        Z_opt = Optim.minimizer(result)
        Neg_LL_opt = Optim.minimum(result)
        #Optimise the MAP for the latent variables season by season
       sol_latent_var = Neg_objective_for_latent_vars_one_season(Z_opt,num,new_x_0)[2]
       new_x_0 = sol_latent_var[:,end]
    end
    # Plot_Solution()
    return nothing
end


function Perform_M_step(funcevals_bb::Int64,iters_nm::Int64)
    Obj(Θ) = neg_objective_for_params(Θ)[1]
    #Optimise the negative log likelihood using current latent variable estimates
    println("Parameter variable optimisation for whole data set")
    opt1 = bbsetup(Obj;
        Method=:adaptive_de_rand_1_bin_radiuslimited,
        SearchRange = SearchRangeForModelParams,
        MaxFuncEvals = funcevals_bb,
        PopulationSize = 40,
        TraceMode = :verbose)
    res_params = bboptimize(opt1)
    Θ_opt = best_candidate(res_params)
    Neg_LL_opt = best_fitness(res_params)
    #Use the global optimisation as starting point for local optimisation
    println("Creating simplex for NM algorithm")
    result = optimize(Obj, Θ_opt, NelderMead(),Optim.Options(iterations = iters_nm,show_trace = true))
    Θ_opt = Optim.minimizer(result)
    Neg_LL_opt = Optim.minimum(result)

    # solve the MLE estimators for the latent variable hyper parameters directly
    Z_array = [P_ModelParams.Z_ξ[2:16]';P_ModelParams.Z_ϕ[2:16]'];
    fit_Z = fit_mle(MvNormal,Z_array)
    P_ModelParams.ξ̄ = fit_Z.μ[1]
    P_ModelParams.ϕ̄ = fit_Z.μ[2]
    P_ModelParams.σ_ξ = fit_Z.Σ.mat[1,1]
    P_ModelParams.ρ_ξϕ = fit_Z.Σ.mat[1,2]
    P_ModelParams.σ_ϕ = fit_Z.Σ.mat[2,2]
    #Return log joint probability
    m = [P_ModelParams.ξ̄,P_ModelParams.ϕ̄]
    Σ = [P_ModelParams.σ_ξ P_ModelParams.ρ_ξϕ;P_ModelParams.ρ_ξϕ P_ModelParams.σ_ϕ]
    neg_ll_for_seasons = 0.
    for i = 2:16
        neg_ll_for_seasons -= logpdf(MvNormal(m, Σ),[P_ModelParams.Z_ξ[i],P_ModelParams.Z_ϕ[i]])
    end
    # Plot_Solution()
    return Θ_opt,(Neg_LL_opt + neg_ll_for_seasons)
end

function Perform_final_optimisation(iters::Int64,opt_method,Θ_0)
    function Obj(Θ)
        Θ_test = Θ_0
        Θ_test[2:4] = Θ
        return neg_objective_for_params(Θ_test)[1]
    end
    # Obj(Θ) = neg_objective_for_params(Θ)[1]
    #Optimise the negative log likelihood using current latent variable estimates
    println("Final optimisation scheme")

    result = optimize(Obj, Θ_0[2:4], opt_method,Optim.Options(iterations = iters,show_trace = true))
    Θ_1 = Optim.minimizer(result)
    Neg_LL_opt = Optim.minimum(result)

    # solve the MLE estimators for the latent variable hyper parameters directly
    Z_array = [P_ModelParams.Z_ξ[2:16]';P_ModelParams.Z_ϕ[2:16]'];
    fit_Z = fit_mle(MvNormal,Z_array)
    P_ModelParams.ξ̄ = fit_Z.μ[1]
    P_ModelParams.ϕ̄ = fit_Z.μ[2]
    P_ModelParams.σ_ξ = fit_Z.Σ.mat[1,1]
    P_ModelParams.ρ_ξϕ = fit_Z.Σ.mat[1,2]
    P_ModelParams.σ_ϕ = fit_Z.Σ.mat[2,2]
    #Return log joint probability
    m = [P_ModelParams.ξ̄,P_ModelParams.ϕ̄]
    Σ = [P_ModelParams.σ_ξ P_ModelParams.ρ_ξϕ;P_ModelParams.ρ_ξϕ P_ModelParams.σ_ϕ]
    neg_ll_for_seasons = 0.
    for i = 2:16
        neg_ll_for_seasons -= logpdf(MvNormal(m, Σ),[P_ModelParams.Z_ξ[i],P_ModelParams.Z_ϕ[i]])
    end
    # Plot_Solution()
    return result,(Neg_LL_opt + neg_ll_for_seasons)
end

function Perform_final_optimisation(iters::Int64,Θ_0,dim,lower,upper)
    function Obj(Θ)
        Θ_test = Θ_0
        Θ_test[dim] = Θ
        return neg_objective_for_params(Θ_test)[1]
    end

    #Optimise the negative log likelihood using current latent variable estimates
    println("Final univariate optimisation scheme")

    result = optimize(Obj, lower, upper)
    Θ_1 = Optim.minimizer(result)
    Neg_LL_opt = Optim.minimum(result)

    # solve the MLE estimators for the latent variable hyper parameters directly
    Z_array = [P_ModelParams.Z_ξ[2:16]';P_ModelParams.Z_ϕ[2:16]'];
    fit_Z = fit_mle(MvNormal,Z_array)
    P_ModelParams.ξ̄ = fit_Z.μ[1]
    P_ModelParams.ϕ̄ = fit_Z.μ[2]
    P_ModelParams.σ_ξ = fit_Z.Σ.mat[1,1]
    P_ModelParams.ρ_ξϕ = fit_Z.Σ.mat[1,2]
    P_ModelParams.σ_ϕ = fit_Z.Σ.mat[2,2]
    #Return log joint probability
    m = [P_ModelParams.ξ̄,P_ModelParams.ϕ̄]
    Σ = [P_ModelParams.σ_ξ P_ModelParams.ρ_ξϕ;P_ModelParams.ρ_ξϕ P_ModelParams.σ_ϕ]
    neg_ll_for_seasons = 0.
    for i = 2:16
        neg_ll_for_seasons -= logpdf(MvNormal(m, Σ),[P_ModelParams.Z_ξ[i],P_ModelParams.Z_ϕ[i]])
    end
    # Plot_Solution()
    return result,(Neg_LL_opt + neg_ll_for_seasons)
end

function Perform_M_step_nm_only(iters_nm::Int64)
    Obj(Θ) = neg_objective_for_params(Θ)[1]
    #Use the global optimisation as starting point for local optimisation
    println("Creating simplex for NM algorithm")
    result = optimize(Obj, Θ_Opt, NelderMead(),Optim.Options(iterations = iters_nm,show_trace = true))
    Θ_opt = Optim.minimizer(result)
    Neg_LL_opt = Optim.minimum(result)

    # solve the MLE estimators for the latent variable hyper parameters directly
    Z_array = [P_ModelParams.Z_ξ[2:16]';P_ModelParams.Z_ϕ[2:16]'];
    fit_Z = fit_mle(MvNormal,Z_array)
    P_ModelParams.ξ̄ = fit_Z.μ[1]
    P_ModelParams.ϕ̄ = fit_Z.μ[2]
    P_ModelParams.σ_ξ = fit_Z.Σ.mat[1,1]
    P_ModelParams.ρ_ξϕ = fit_Z.Σ.mat[1,2]
    P_ModelParams.σ_ϕ = fit_Z.Σ.mat[2,2]
    #Return log joint probability
    m = [P_ModelParams.ξ̄,P_ModelParams.ϕ̄]
    Σ = [P_ModelParams.σ_ξ P_ModelParams.ρ_ξϕ;P_ModelParams.ρ_ξϕ P_ModelParams.σ_ϕ]
    neg_ll_for_seasons = 0.
    for i = 2:16
        neg_ll_for_seasons -= logpdf(MvNormal(m, Σ),[P_ModelParams.Z_ξ[i],P_ModelParams.Z_ϕ[i]])
    end
    # Plot_Solution()
    return Θ_opt,(Neg_LL_opt + neg_ll_for_seasons)
end
