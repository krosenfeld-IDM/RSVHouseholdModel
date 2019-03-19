#Plotting functions for household epidemic

#Declare the plotting backend
# plotlyjs()
gr()

function CreateBulkPlot(sol)
    S1_vec = States[:,1];
    I1_vec = States[:,2];
    R1_vec = States[:,3];
    S2_vec = States[:,4];
    I2_vec = States[:,5];
    R2_vec = States[:,6];

    S1 = [vecdot(sol[(M_a+1):end,n],S1_vec) for n = 1:length(sol.t)]
    I1 = [vecdot(sol[(M_a+1):end,n],I1_vec) for n = 1:length(sol.t)]
    R1 = [vecdot(sol[(M_a+1):end,n],R1_vec) for n = 1:length(sol.t)]
    S2 = [vecdot(sol[(M_a+1):end,n],S2_vec) for n = 1:length(sol.t)]
    I2 = [vecdot(sol[(M_a+1):end,n],I2_vec) for n = 1:length(sol.t)]
    R2 = [vecdot(sol[(M_a+1):end,n],R2_vec) for n = 1:length(sol.t)]
    plt = plot(sol.t,S1,lab="S1",reuse = false);
    plot!(plt,sol.t,I1,lab="I1");
    plot!(plt,sol.t,R1,lab="R1");
    plot!(plt,sol.t,S2,lab="S2");
    plot!(plt,sol.t,I2,lab="I2");
    plot!(plt,sol.t,R2,lab="R2");
    return plt
end

function CreatePlotHHSizeRestriction(sol,HHSize)
    S1_vec = States[:,1].*(N_vect .== HHSize);
    I1_vec = States[:,2].*(N_vect .== HHSize);
    R1_vec = States[:,3].*(N_vect .== HHSize);
    S2_vec = States[:,4].*(N_vect .== HHSize);
    I2_vec = States[:,5].*(N_vect .== HHSize);
    R2_vec = States[:,6].*(N_vect .== HHSize);

    S1 = [vecdot(sol[(M_a+1):end,n],S1_vec) for n = 1:length(sol.t)]
    I1 = [vecdot(sol[(M_a+1):end,n],I1_vec) for n = 1:length(sol.t)]
    R1 = [vecdot(sol[(M_a+1):end,n],R1_vec) for n = 1:length(sol.t)]
    S2 = [vecdot(sol[(M_a+1):end,n],S2_vec) for n = 1:length(sol.t)]
    I2 = [vecdot(sol[(M_a+1):end,n],I2_vec) for n = 1:length(sol.t)]
    R2 = [vecdot(sol[(M_a+1):end,n],R2_vec) for n = 1:length(sol.t)]
    plt = plot(sol.t,S1,lab="S1",reuse = false);
    plot!(plt,sol.t,I1,lab="I1");
    plot!(plt,sol.t,R1,lab="R1");
    plot!(plt,sol.t,S2,lab="S2");
    plot!(plt,sol.t,I2,lab="I2");
    plot!(plt,sol.t,R2,lab="R2");
    return plt
end

function CreateBulkPlot_just_u1s(sol)
    S1_vec = States[:,1];
    I1_vec = States[:,2];
    R1_vec = States[:,3];

    S1 = [vecdot(sol[(M_a+1):end,n],S1_vec) for n = 1:length(sol.t)]
    I1 = [vecdot(sol[(M_a+1):end,n],I1_vec) for n = 1:length(sol.t)]
    R1 = [vecdot(sol[(M_a+1):end,n],R1_vec) for n = 1:length(sol.t)]

    plt = plot(sol.t,S1,lab="S1",reuse = false);
    plot!(plt,sol.t,I1,lab="I1");
    plot!(plt,sol.t,R1,lab="R1");
    return plt
end

function CreatePlotHHSizeRestriction_just_u1s(sol,HHSize)
    S1_vec = States[:,1].*(N_vect .== HHSize);
    I1_vec = States[:,2].*(N_vect .== HHSize);
    R1_vec = States[:,3].*(N_vect .== HHSize);

    S1 = [vecdot(sol[(M_a+1):end,n],S1_vec) for n = 1:length(sol.t)]
    I1 = [vecdot(sol[(M_a+1):end,n],I1_vec) for n = 1:length(sol.t)]
    R1 = [vecdot(sol[(M_a+1):end,n],R1_vec) for n = 1:length(sol.t)]
    plt = plot(sol.t,S1,lab="S1",reuse = false);
    plot!(plt,sol.t,I1,lab="I1");
    plot!(plt,sol.t,R1,lab="R1");
    return plt
end


function CreateBulkPlot_just_o1s(sol)
    S2_vec = States[:,4];
    I2_vec = States[:,5];
    R2_vec = States[:,6];

    S2 = [vecdot(sol[(M_a+1):end,n],S2_vec) for n = 1:length(sol.t)]
    I2 = [vecdot(sol[(M_a+1):end,n],I2_vec) for n = 1:length(sol.t)]
    R2 = [vecdot(sol[(M_a+1):end,n],R2_vec) for n = 1:length(sol.t)]

    plt = plot(sol.t,S2,lab="S2",reuse = false);
    plot!(plt,sol.t,I2,lab="I2");
    plot!(plt,sol.t,R2,lab="R2");
    return plt
end

function CreatePlotHHSizeRestriction_just_01s(sol,HHSize)
    S2_vec = States[:,4].*(N_vect .== HHSize);
    I2_vec = States[:,5].*(N_vect .== HHSize);
    R2_vec = States[:,6].*(N_vect .== HHSize);

    S2 = [vecdot(sol[(M_a+1):end,n],S2_vec) for n = 1:length(sol.t)]
    I2 = [vecdot(sol[(M_a+1):end,n],I2_vec) for n = 1:length(sol.t)]
    R2 = [vecdot(sol[(M_a+1):end,n],R2_vec) for n = 1:length(sol.t)]
    plt = plot(sol.t,S2,lab="S2",reuse = false);
    plot!(plt,sol.t,I2,lab="I2");
    plot!(plt,sol.t,R2,lab="R2");
    return plt
end

function create_pop_size_plot(sol)
    d1,d2 = size(sol)
    N1 = [vecdot(sol[(M_a+1):end,n],N1_vec) for n = 1:length(sol.t)]
    N2 = [vecdot(sol[(M_a+1):end,n],N2_vec) for n = 1:length(sol.t)]
    plt = plot(sol.t,N1,lab="N1",reuse=false);
    plot!(plt,sol.t,N2,lab="N2");
    return plt
end

function create_pop_size_plot_for_o1s_or_u1s(sol,n)
    d1,d2 = size(sol)
    N1 = [vecdot(sol[(M_a+1):end,n],N1_vec) for n = 1:length(sol.t)]
    N2 = [vecdot(sol[(M_a+1):end,n],N2_vec) for n = 1:length(sol.t)]
    if n == 1
        plt = plot(sol.t,N1,lab="N1",reuse=false);
    elseif n == 2
        plt = plot(sol.t,N2,lab="N2",reuse=false);
    end
    return plt
end


function run_model_with_initial_params(Θ_initial)
    x_0 = [zeros(M_a);construct_susceptible_HH_distribution_from_year(N_HEachYear[5],d1)];
    P_ModelParams.b_C = 0.#Θ_initial[1]
    P_ModelParams.b_S = Θ_initial[1]
    P_ModelParams.b_A = Θ_initial[2]
    SetMixingMatrix!(P_ModelParams) #This updates with within population parameter changes
    P_ModelParams.α = 1/Θ_initial[3] #Convert from rate to mean
    P_ModelParams.τ = Θ_initial[4]
    P_ModelParams.ϵ = 10.
    P_ModelParams.ξ̄ = Θ_initial[5]
    P_ModelParams.ϕ̄ = Θ_initial[6]
    P_ModelParams.Z_ξ = P_ModelParams.ξ̄*ones(17)
    P_ModelParams.Z_ϕ = P_ModelParams.ϕ̄*ones(17)
    set_matrices_for_season!(4,P_ModelParams)
    #This updates with each year (= season + 1) because of changes in the household recycling rate
    prob_preliminary = ODEProblem(F,x_0,(Float64((SeasonStartDates[4] - Date(2000,1,1)).value),Float64((SeasonStartDates[5] - Date(2000,1,1)).value)),P_ModelParams);
    sol_beginning = solve(prob_preliminary,CVODE_BDF(linear_solver=:GMRES),abstol = 1e-3, reltol = 1e-2)
    a = sol_beginning[1:M_a,end]
    z_0_new = rescale_endstate(sol_beginning[(M_a+1):end,end],6,N_HEachYear)
    z_0_new[z_0_new .< 0.] = 0.
    t_sol = []
    cum_incidence_sol = zeros(M_a)
    #Go through each subsequent season
    for n = 5:9
        set_matrices_for_season!(n,P_ModelParams)
        prob_season = ODEProblem(F,[a;z_0_new],(Float64((SeasonStartDates[n] - Date(2000,1,1)).value),Float64((SeasonEndDates[n] - Date(2000,1,1)).value)),P_ModelParams)
        valid_prediction_times = pred_times[[(pred_times[i] >= prob_season.tspan[1])&&(pred_times[i] <= prob_season.tspan[2]) for i = 1:length(pred_times) ]]
        sol_to_end = solve(prob_season,CVODE_BDF(linear_solver=:GMRES),abstol = 1e-3,reltol = 1e-2,saveat = valid_prediction_times)
        a = sol_to_end[1:M_a,end]
        b = sol_to_end[(M_a+1):end,end]
        z_0_new = rescale_endstate(b,n+2,N_HEachYear)
        z_0_new[z_0_new .< 0.] = 0.
        t_sol = [t_sol;sol_to_end.t[2:end]]
        inc = sol_to_end[1:M_a,2:end]
        cum_incidence_sol = [cum_incidence_sol inc]
    end
    return t_sol,cum_incidence_sol[:,2:end]
end

function create_incidence_plot_u1(sol,age_cat::Int64)
    times = relevant_pred_times_from_sol(sol,pred_times)
    incidence = find_model_incidence_using_interpolation(sol,pred_times)
    plot(sol.t[2:end],incidence[age_cat])
end

function create_true_incidence_plot_age_cat(sol,age_cat::Int64)
    times = relevant_pred_times_from_sol(sol,pred_times)
    incidence = get_trueincidencedata_age_cat(sol,pred_times,age_cat)
    scatter(times,incidence)
end

function create_comparison_plot_age_cat(sol,age_cat::Int64)
    times = relevant_pred_times_from_sol(sol,pred_times)
    true_incidence = get_trueincidencedata_age_cat(sol,pred_times,age_cat)
    model_incidence = find_model_incidence_using_interpolation_age_cat(sol,pred_times,age_cat)
    scatter(times,true_incidence)
    plot!(times[2:end],model_incidence)
end

function create_comparison_plot(sol)
    times = relevant_pred_times_from_sol(sol,pred_times)
    R = find_KDHSS_ratio(times)
    date_times = Dates.Day.(round.(Int64,times)) + Date(2000,1,1)
    true_incidence = get_trueincidencedata_age_cat(sol,pred_times,1)
    model_incidence = find_model_incidence_using_interpolation_age_cat(sol,pred_times,1)
    TI = zeros(length(true_incidence),16)
    MI = zeros(length(true_incidence),16)
    TI[:,1] = true_incidence
    MI[:,1] = model_incidence
    for n =2:16
        x = get_trueincidencedata_age_cat(sol,pred_times,n)
        y = find_model_incidence_using_interpolation_age_cat(sol,pred_times,n)
        TI[:,n] = x
        MI[:,n] = y
        true_incidence += x
        model_incidence += y
    end
    model_incidence.*R[2:end]
    # scatter(times,true_incidence)
    # plot!(times[2:end],model_incidence)
    plt = scatter(date_times[2:end],true_incidence,lab = "Incidence")
    plot!(plt,date_times[2:end],model_incidence, lab = "Model pred.",lw=3)
    title!(plt,"Max HH size = $(MaxHouseholdSize), max num U1s per HH = $(MaxNumberOfU1s)")
    return times,TI,MI,plt
end

function create_comparison_plot(sol,Agg_scale::Int64)
    times = relevant_pred_times_from_sol(sol,pred_times)
    R = find_KDHSS_ratio(times)
    date_times = Dates.Day.(round.(Int64,times)) + Date(2000,1,1)
    true_incidence = get_trueincidencedata_age_cat(sol,pred_times,1)
    model_incidence = find_model_incidence_using_interpolation_age_cat(sol,pred_times,1)
    TI = zeros(length(true_incidence),16)
    MI = zeros(length(true_incidence),16)
    TI[:,1] = true_incidence
    MI[:,1] = model_incidence
    for n =2:16
        x = get_trueincidencedata_age_cat(sol,pred_times,n)
        y = find_model_incidence_using_interpolation_age_cat(sol,pred_times,n)
        TI[:,n] = x
        MI[:,n] = y
        true_incidence += x
        model_incidence += y
    end
    model_incidence.*R[2:end]
    date_times = date_times[Agg_scale:Agg_scale:length(date_times)]
    true_incidence = movsum(true_incidence,Agg_scale)[2:end]
    model_incidence = movsum(model_incidence,Agg_scale)[2:end]
    # scatter(times,true_incidence)
    # plot!(times[2:end],model_incidence)
    plt = scatter(date_times,true_incidence,lab = "Incidence")
    plot!(plt,date_times,model_incidence, lab = "Model pred.",lw=3)
    title!(plt,"Max HH size = $(MaxHouseholdSize), max num U1s per HH = $(MaxNumberOfU1s)")
    return times,TI,MI,plt
end

function plot_forecast_solution(sol,ts)
    times = sol.t[1]:ts:sol.t[end]
    date_times = [Date(2000,1,1) + Dates.Day(round(Int64,t)) for t in times]
    R = Ratio(sol.t[1])
    model_incidence_array = [diff([sol(t)[i] for t in times]).*R for i = 1:16]
    MI = zeros(length(model_incidence_array[1]),16)
    for n =1:16
        MI[:,n] = model_incidence_array[n]
    end
    # scatter(times,true_incidence)
    # plot!(times[2:end],model_incidence)
    plt = plot(date_times[2:end],sum(MI,2),lab = "",lw=1.,color=:red)
    title!(plt,"Forecast RSV Hospitalisation")
    return collect(times),MI,plt
end

function plot_forecast_solution(sol,ts,age_cat::Int64)
    times = sol.t[1]:ts:sol.t[end]
    date_times = [Date(2000,1,1) + Dates.Day(round(Int64,t)) for t in times]
    R = Ratio(sol.t[1])
    model_incidence_array = [diff([sol(t)[i] for t in times]).*R for i = 1:16]
    MI = zeros(length(model_incidence_array[1]),16)
    for n =1:16
        MI[:,n] = model_incidence_array[n]
    end
    # scatter(times,true_incidence)
    # plot!(times[2:end],model_incidence)
    plt = plot(date_times[2:end],MI[:,age_cat],lab = "Forecast Incidence for age cat $(age_cat)")
    title!(plt,"Forecast RSV Hospitalisation")
    return collect(times),MI,plt
end
function plot_forecast_solution(inc::Vector{Array{Float64,1}},ts)
    times = FutureSeasonStartTimes[1]:ts:FutureSeasonEndTimes[end]
    date_times = [Date(2000,1,1) + Dates.Day(round(Int64,t)) for t in times]
    R = Ratio(FutureSeasonStartTimes[1])
    model_incidence_array = [diff([v[i] for v in inc]).*R for i = 1:16]
    MI = zeros(length(model_incidence_array[1]),16)
    for n =1:16
        MI[:,n] = model_incidence_array[n]
    end
    # scatter(times,true_incidence)
    # plot!(times[2:end],model_incidence)
    plt = plot(date_times[2:end],sum(MI,2),lab = "")
    title!(plt,"Forecast RSV Hospitalisation")
    return collect(times),MI,plt
end
function plot_forecast_solution(inc_array::Vector{Array{Array{Float64,1},1}},ts)
    times = FutureSeasonStartTimes[1]:ts:FutureSeasonEndTimes[end]
    date_times = [Date(2000,1,1) + Dates.Day(round(Int64,t)) for t in times]
    R = Ratio(FutureSeasonStartTimes[1])
    inc = inc_array[1]
    model_incidence_array = [diff([v[i] for v in inc]).*R for i = 1:16]
    MI = zeros(length(model_incidence_array[1]),16)
    for n =1:16
        MI[:,n] = model_incidence_array[n]
    end
    # scatter(times,true_incidence)
    # plot!(times[2:end],model_incidence)
    plt = plot(date_times[2:end],sum(MI,2),lab = "",color = :grey,lw=0.5)
    title!(plt,"Forecast RSV Hospitalisation")
    # mean = sum(MI,2)
    for i = 2:length(inc_array)
        n = Float64(i)
        inc = inc_array[i]
        model_incidence_array = [diff([v[i] for v in inc]).*R for i = 1:16]
        MI = zeros(length(model_incidence_array[1]),16)
        for n =1:16
            MI[:,n] = model_incidence_array[n]
        end
        # scatter(times,true_incidence)
        # plot!(times[2:end],model_incidence)
        plot!(plt,date_times[2:end],sum(MI,2),lab = "",color = :grey,lw=0.5)
        # mean += ((n-1.)/n)*mean + (sum(MI,2)/n)
    end
    # plot!(plt,date_times[2:end],mean,lab = "",color = :red,lw=2)

    return collect(times),MI,plt
end
