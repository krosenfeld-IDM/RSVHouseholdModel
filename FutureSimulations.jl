#Alterations to the vector field in vaccination
#NB: Default is to not rescale incidence by hospitalisation probability
#There is an option to redraw random seasonality realisations or load from seed, default is load from seed
#NB: population is treated as stationary with the same population structure as 2016

function F(du,u,P::HH_RSV_VaccinationModelParameters,t)
    #Get the numbers of households in each state and the cumulative incidence
    q = @view du[1:M_a];
    p = @view du[(M_a+1):end];
    z = @view u[(M_a+1):end];
    # Seasonality --- check what season it is, and if season has changed then update
    current_season = find_season_from_time(t)
    # set_HH_growthrate!(z,P)
    if ~(P.Pre.season_curr == current_season)
        P.Pre.season_curr = current_season
    end
    # P.Pre.S = P.Χ(t)
    Seasonality = Χ(t,P)
    # # Calculate incidence rates at different ages and different household configurations
    #check if the state and time have been updated, if they have then recalculate the transmission rates and update the current state
    #and time, otherwise don't
    if ~((P.Pre.x_curr == u) && (P.Pre.t_curr == t))
        TransmissionIncidence!(z,P)
        P.Pre.x_curr = u
        P.Pre.t_curr = t
    end

    # # Vector field for hospitalisation rate amongst U1s and O1s and Age specific hospitalisation rates due to within household infection
    CondAgeDistrib = CondAgeOfInfectedU1(1/P.α,CurrentSolidProtectionDuration,MABCoverage)
    RateOfU1Incidence = mean(P.Pre.S1_A)*mean(P.Pre.λ_A[P.U1_cats])
    for i = 1:M_a
        if P.U1_cats[i]
            # q[i] = Seasonality*(RateOfU1Incidence*CondAgeDistrib[i] + P.Pre.r_h[i])*P.HR[i]
            q[i] = Seasonality*(RateOfU1Incidence*CondAgeDistrib[i] + P.Pre.r_h[i])
        else
            # q[i] = Seasonality*(P.Pre.S_A[i]*P.Pre.λ_A[i] + P.Pre.r_h[i])*P.HR[i]
            q[i] = Seasonality*(P.Pre.S_A[i]*P.Pre.λ_A[i] + P.Pre.r_h[i])
        end
    end
    # # # Vector field for change in household configurations  (1 - AvMaternalProtU1s(1/P.α))
    # p .= P.CRE_mat*z .+ P.Pre.S*P.τ*P.IntInf_mat*z .+ P.Pre.S*P.ExtInfU1_mat*P.Pre.λ_state_U1s .+ P.Pre.S*P.ExtInfO1_mat*P.Pre.λ_state_O1s;
    p .= P.CRE_mat*z .+
        P.ReplacementRate_mat*z .+
        (1 - AvMaternalProtU1s(1/P.α,CurrentSolidProtectionDuration,MABCoverage))*Seasonality*P.τ*P.IntInfU1_mat*z .+
        Seasonality*P.τ*P.IntInfO1_mat*z .+
        Seasonality*P.ExtInfU1_mat*P.Pre.λ_state_U1s .+
        Seasonality*P.ExtInfO1_mat*P.Pre.λ_state_O1s
    return nothing
end


#Set the random seasonality parameters
function create_random_seasonality_series(d)
    x = Vector{Tuple{Array{Float64,1},Array{Float64,1}}}(NumOfMCReps)
    for n = 1:NumOfMCReps
        ξ_series = zeros(length(FutureSeasonStartDates))
        ϕ_series = zeros(length(FutureSeasonStartDates))
        for i = 1:length(FutureSeasonStartDates)
            ThisSeason = rand(d)
            ξ_series[i] = ThisSeason[1]
            ϕ_series[i] = ThisSeason[2]
        end
        x[n] = (ξ_series,ϕ_series)
    end
    return x
end

#Comment/uncomment here to either draw and save a new set of realisations for random seasonality or load from seed
# RandSeasonalitySeries = create_random_seasonality_series(SeasonalDistrib)
# @save "Seed_random_seasonality_$(SchoolsOnlyR_0).jld" RandSeasonalitySeries

RandSeasonalitySeries = load("Seed_random_seasonality_$(SchoolsOnlyR_0).jld")["RandSeasonalitySeries"]

function set_seasonality!(P::HH_RSV_VaccinationModelParameters,i,RandSeasonalitySeries)
    P.Z_ξ = RandSeasonalitySeries[i][1]
    P.Z_ϕ = RandSeasonalitySeries[i][2]
    return nothing
end

function set_seasonality!(prob,i,RandSeasonalitySeries)
    prob.p.Z_ξ = RandSeasonalitySeries[i][1]
    prob.p.Z_ϕ = RandSeasonalitySeries[i][2]
    return nothing
end

function set_seasonality(prob,i,RandSeasonalitySeries)
    prob.p.Z_ξ = RandSeasonalitySeries[i][1]
    prob.p.Z_ϕ = RandSeasonalitySeries[i][2]
    return prob
end

#Overload the between household calculations

function InplaceStateToHouseholdConversion!(z,P::HH_RSV_VaccinationModelParameters)
    #Calculate the household size to household size transmission rate
    A_mul_B!(P.Pre.S1_H,P.ConvMatrix_S1_H,z)
    A_mul_B!(P.Pre.I1_H,P.ConvMatrix_I1_H,z)
    A_mul_B!(P.Pre.N1_H,P.ConvMatrix_N1,z)
    A_mul_B!(P.Pre.N2_0_H,P.ConvMatrix_N2_0,z)
    A_mul_B!(P.Pre.N2_1_H,P.ConvMatrix_N2_1,z)
    return nothing
end

function InPlaceHouseholdToAgeConversion!(z,P::HH_RSV_VaccinationModelParameters)
    #Number of infecteds and susceptibles in each conditional category for each household size
    P.Pre.N1_H[1] = 1.
    P.Pre.N2_1_H[1] = 1.
    #Convert into a prediction about the numbers of infecteds and susceptibles in each age category
    P.Pre.I1_A .= mean(P.Pre.I1_H)*ones(sum(P.U1_cats))
    A_mul_B!(P.Pre.I2_A,P.Inf_AgivState,z)
    P.Pre.I_A .= [P.Pre.I1_A;inf_2*P.Pre.I2_A] #Over 1s are less infectious
    P.Pre.S1_A .=  (1-AvMaternalProtU1s(1/P.α,CurrentSolidProtectionDuration,MABCoverage))*mean(P.Pre.S1_H)*ones(sum(P.U1_cats)) #Accounts for possibly being maternally protected
    A_mul_B!(P.Pre.S2_A,P.Sus_AgivState,z)
    P.Pre.S_A .= [P.Pre.S1_A;sus_2*P.Pre.S2_A]#Effective susceptibility is lower for O1s
    P.Pre.N1_A .= mean(P.Pre.S1_H)*ones(sum(P.U1_cats))
    A_mul_B!(P.Pre.N2_A,P.N_AgivState,z)
    P.Pre.N_A .= [P.Pre.N1_A;P.Pre.N2_A]
    return nothing
end

function InPlaceAgedepForceOfInfection!(z,P::HH_RSV_VaccinationModelParameters)
    #Force of infection at the level of the population outside of homestead
    (P.Pre.N,P.Pre.N2) = get_population_size(z)
    A_mul_B!(P.Pre.λ_A,P.MixingMatrix,(P.Pre.I_A + P.ϵ)/P.Pre.N) #note the scaling by 1/total population size
    return nothing
end

function InPlaceExternalForceOfInfOnHouseholdCalculation!(z,P::HH_RSV_VaccinationModelParameters)
    #Expected number of infectious contact for each housesize
    λ_A_U1s = @view P.Pre.λ_A[U1_cats]
    λ_A_O1s = @view P.Pre.λ_A[.~U1_cats]
    P.Pre.λ_H_U1s = P.U1s_HU_giv_A*λ_A_U1s
    P.Pre.λ_H_U1s[1] = 0.;
    P.Pre.λ_H_0_O1s = P.O1s_HU_giv_A_0*λ_A_O1s
    P.Pre.λ_H_1_O1s = P.O1s_HU_giv_A_1*λ_A_O1s
    P.Pre.λ_H_1_O1s[1] = 0.
    return nothing
end


function InPlaceInternalForceOfInfOnHouseholdCalculation!(z,P::HH_RSV_VaccinationModelParameters)
    # #Within households
    # #Household type specific infection rates
    P.Pre.r1 = P.τ*CondAgeOfInfectedU1(1/P.α,CurrentSolidProtectionDuration)*(1-AvMaternalProtU1s(1/P.α,CurrentSolidProtectionDuration,MABCoverage))*sum(P.ConvMatrix_InfRate1_H*z);
    P.Pre.r2 = sus_2*P.τ*P.Inf_Rate_AgivState*z
    P.Pre.r_h = [P.Pre.r1;P.Pre.r2]
    return nothing
end

function RescaleHouseholdToStates!(z,P::HH_RSV_VaccinationModelParameters)
    # #Give rescaled household states vectors
    for i = 1:d1
        P.Pre.λ_state_U1s[i] = P.Pre.λ_H_U1s[N_vect[i]]*z[i]
        if P.U1_state[i]
            P.Pre.λ_state_O1s[i] = P.Pre.λ_H_1_O1s[N_vect[i]]*z[i]
        else
            P.Pre.λ_state_O1s[i] = P.Pre.λ_H_0_O1s[N_vect[i]]*z[i]
        end
    end
    return nothing
end

function TransmissionIncidence!(z,P::HH_RSV_VaccinationModelParameters)
    InplaceStateToHouseholdConversion!(z,P)
    InPlaceHouseholdToAgeConversion!(z,P)
    InPlaceAgedepForceOfInfection!(z,P)
    InPlaceExternalForceOfInfOnHouseholdCalculation!(z,P)
    InPlaceInternalForceOfInfOnHouseholdCalculation!(z,P)
    RescaleHouseholdToStates!(z,P)
    return nothing
end



#Declare the main objects to for simulation
set_matrices_for_year!(P_ModelParams,100)
P_VacModel = create_RSV_vaccination_model(P_ModelParams)
generate_replacement_rate!(P_VacModel,100)
prob_vac = ODEProblem(F,x_0,(FutureSeasonStartTimes[1],FutureSeasonEndTimes[10]),P_VacModel)


function output_func_for_rsv_vac_timeseries(sol,i)
    times = sol.t[1]:7.:sol.t[end]
    u = [sol(t)[1:16] for t in times]
    (u,false)
end

function output_func_for_rsv_vac_endpoint(sol,i)
    u = sol(sol.t[end])[1:16]
    (u,false)
end

function output_func_for_rsv_inf_endpoint(sol,i)
    u = sol(sol.t[end])[1:M_a]
    (u,false)
end

function calculate_rate_of_vaccinations(z,P::HH_RSV_VaccinationModelParameters,year)
    if year < 1
        index = 1
    end
    if year >=1 && year <= length(N_HEachYear)
        index = year
    end
    if year > length(N_HEachYear)
        index = length(N_HEachYear)
    end
#Find the household size dependent per-adult replacement rate
    μ_H = FindReplacementRatesForYear(N_HEachYear,index,η_1,MaxNumberOfU1s)
    birth_rate_hh = μ_H[N_vect].*N2_vec.*z
    return sum(birth_rate_hh),HHCoverage*vecdot(birth_rate_hh,N2_vec-1)
end

function output_func_for_rsv_number_vaccines(sol,i)
    times = sol.t[1]:1.:sol.t[end]
    num_vacs = [calculate_rate_of_vaccinations(sol(t)[(M_a+1):end],sol.prob.p,100) for t in times]
    u = (sum([x[1] for x in num_vacs]) , sum([x[2] for x in num_vacs])  )
    (u,false)
end

#This fuction defines the random seasonality for each simulation
prob_func_for_rsv_sims = (prob,i,repeat)-> set_seasonality(prob,i,RandSeasonalitySeries)
