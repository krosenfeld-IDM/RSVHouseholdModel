#Constructs the vector field (the RHS of the fundamental dynamical equation of the model) for the household model

function set_conversion_matrices!(P::HH_RSVModelParameters)
    InvEffHHSize = (max.(1,N_vect-1)).^(-P.EffHHSizePower);
    ConvMatrix_I1_H = zeros(MaxHouseholdSize,d1)
    ConvMatrix_S1_H = zeros(MaxHouseholdSize,d1)
    ConvMatrix_S2_0_H = zeros(MaxHouseholdSize,d1)
    ConvMatrix_S2_1_H = zeros(MaxHouseholdSize,d1)
    ConvMatrix_I2_0_H = zeros(MaxHouseholdSize,d1)
    ConvMatrix_I2_1_H = zeros(MaxHouseholdSize,d1)
    ConvMatrix_InfRate1_H = zeros(MaxHouseholdSize,d1)
    ConvMatrix_InfRate2_0_H = zeros(MaxHouseholdSize,d1)
    ConvMatrix_InfRate2_1_H = zeros(MaxHouseholdSize,d1)
    ConvMatrix_N1 = zeros(MaxHouseholdSize,d1)
    ConvMatrix_N2_0 = zeros(MaxHouseholdSize,d1)
    ConvMatrix_N2_1 = zeros(MaxHouseholdSize,d1)
    for i = 1:MaxHouseholdSize
        for j = 1:d1
            if N_vect[j] == i
                ConvMatrix_I1_H[i,j] = States[j,2];
                ConvMatrix_S1_H[i,j] = States[j,1];
                ConvMatrix_S2_0_H[i,j] = States[j,4]*(~U1_state[j]);
                ConvMatrix_S2_1_H[i,j] = States[j,4]*(U1_state[j]);
                ConvMatrix_I2_0_H[i,j] = States[j,5]*(~U1_state[j]);
                ConvMatrix_I2_1_H[i,j] = States[j,5]*(U1_state[j]);
                ConvMatrix_InfRate1_H[i,j] = States[j,1]*(States[j,2] + inf_2*States[j,5])*InvEffHHSize[j];
                ConvMatrix_InfRate2_0_H[i,j] = (States[j,4]*(States[j,2] + inf_2*States[j,5])*InvEffHHSize[j])*(~U1_state[j]);
                ConvMatrix_InfRate2_1_H[i,j] = (States[j,4]*(States[j,2] + inf_2*States[j,5])*InvEffHHSize[j])*(U1_state[j]);
                ConvMatrix_N1[i,j] = N1_vec[j];
                ConvMatrix_N2_0[i,j] = N2_vec[j]*(~U1_state[j])
                ConvMatrix_N2_1[i,j] = N2_vec[j]*(U1_state[j])
            end
        end
    end
    ConvMatrix_InfRate1_H = sparse(ConvMatrix_InfRate1_H)
    ConvMatrix_InfRate2_0_H = sparse(ConvMatrix_InfRate2_0_H)
    ConvMatrix_InfRate2_1_H = sparse(ConvMatrix_InfRate2_1_H)
    ConvMatrix_I1_H = sparse(ConvMatrix_I1_H)
    ConvMatrix_S1_H = sparse(ConvMatrix_S1_H)
    ConvMatrix_S2_0_H = sparse(ConvMatrix_S2_0_H)
    ConvMatrix_S2_1_H = sparse(ConvMatrix_S2_1_H)
    ConvMatrix_I2_0_H = sparse(ConvMatrix_I2_0_H)
    ConvMatrix_I2_1_H = sparse(ConvMatrix_I2_1_H)
    P.ConvMatrix_S1_H = ConvMatrix_S1_H
    P.ConvMatrix_I1_H = ConvMatrix_I1_H
    P.ConvMatrix_N1 = ConvMatrix_N1
    P.ConvMatrix_N2_0 = ConvMatrix_N2_0
    P.ConvMatrix_N2_1 = ConvMatrix_N2_1
    P.ConvMatrix_InfRate1_H = ConvMatrix_InfRate1_H
    P.ConvMatrix_I2_0_H = ConvMatrix_I2_0_H
    P.ConvMatrix_I2_1_H = ConvMatrix_I2_1_H
    P.ConvMatrix_S2_0_H = ConvMatrix_S2_0_H
    P.ConvMatrix_S2_1_H = ConvMatrix_S2_1_H
    P.ConvMatrix_InfRate2_0_H = ConvMatrix_InfRate2_0_H
    P.ConvMatrix_InfRate2_1_H = ConvMatrix_InfRate2_1_H
    return nothing
end
set_conversion_matrices!(P_ModelParams)

get_population_size(z) = (vecdot(N_vect,z),vecdot(N2_vec,z))

function set_conditionalmatrices_thatchangeeachyear!(P::HH_RSVModelParameters,year::Int64)
    if year < 1
        year = 1
    end
    if year > length(N_HEachYear)
        year = length(N_HEachYear)
    end
    #Construct the distributions for age given household type and being O1 and if in home with an U1
    Y1 = P_AgivHU[year][.~U1_cats,:,1]
    Y2 = P_AgivHU[year][.~U1_cats,:,2]
    s1,s2 = size(Y2)
    for j = 2:s2
        Y2[:,j] = normalize!(Y2[:,j],1)
    end
    Inf_AgivState = sparse(Y1*P.ConvMatrix_I2_0_H + Y2*P.ConvMatrix_I2_1_H)
    Sus_AgivState = sparse(Y1*P.ConvMatrix_S2_0_H + Y2*P.ConvMatrix_S2_1_H)
    N_AgivState = sparse(Y1*P.ConvMatrix_N2_0 + Y2*P.ConvMatrix_N2_1)
    Inf_Rate_AgivState = sparse(Y1*P.ConvMatrix_InfRate2_0_H + Y2*P.ConvMatrix_InfRate2_1_H)
    #Construct the distributions for household type given age
    U1s_HU_giv_A = P_HUgivA[year][U1_cats,:,2]';
    O1s_HU_giv_A_0 = P_HUgivA[year][.~U1_cats,:,1]'
    O1s_HU_giv_A_1 = P_HUgivA[year][.~U1_cats,:,2]'
    P.Inf_AgivState = Inf_AgivState
    P.Sus_AgivState = Sus_AgivState
    P.N_AgivState = N_AgivState
    P.Inf_Rate_AgivState = Inf_Rate_AgivState
    P.U1s_HU_giv_A = U1s_HU_giv_A
    P.O1s_HU_giv_A_0 = O1s_HU_giv_A_0
    P.O1s_HU_giv_A_1 = O1s_HU_giv_A_1
    return nothing
end

function InplaceStateToHouseholdConversion!(z,P::HH_RSVModelParameters)
    #Calculate the household size to household size transmission rate
    A_mul_B!(P.Pre.S1_H,P.ConvMatrix_S1_H,z)
    A_mul_B!(P.Pre.I1_H,P.ConvMatrix_I1_H,z)
    A_mul_B!(P.Pre.N1_H,P.ConvMatrix_N1,z)
    A_mul_B!(P.Pre.N2_0_H,P.ConvMatrix_N2_0,z)
    A_mul_B!(P.Pre.N2_1_H,P.ConvMatrix_N2_1,z)
    return nothing
end

function InPlaceHouseholdToAgeConversion!(z,P::HH_RSVModelParameters)
    #Number of infecteds and susceptibles in each conditional category for each household size
    P.Pre.N1_H[1] = 1.
    P.Pre.N2_1_H[1] = 1.
    #Convert into a prediction about the numbers of infecteds and susceptibles in each age category
    P.Pre.I1_A .= mean(P.Pre.I1_H)*ones(sum(P.U1_cats))
    A_mul_B!(P.Pre.I2_A,P.Inf_AgivState,z)
    P.Pre.I_A .= [P.Pre.I1_A;inf_2*P.Pre.I2_A] #Over 1s are less infectious
    P.Pre.S1_A .=  (1-AvMaternalProtU1s(1/P.α,CurrentSolidProtectionDuration))*mean(P.Pre.S1_H)*ones(sum(P.U1_cats)) #Accounts for possibly being maternally protected
    A_mul_B!(P.Pre.S2_A,P.Sus_AgivState,z)
    P.Pre.S_A .= [P.Pre.S1_A;sus_2*P.Pre.S2_A]#Effective susceptibility is lower for O1s
    P.Pre.N1_A .= mean(P.Pre.S1_H)*ones(sum(P.U1_cats))
    A_mul_B!(P.Pre.N2_A,P.N_AgivState,z)
    P.Pre.N_A .= [P.Pre.N1_A;P.Pre.N2_A]
    return nothing
end

function InPlaceAgedepForceOfInfection!(z,P::HH_RSVModelParameters)
    #Force of infection at the level of the population outside of homestead
    (P.Pre.N,P.Pre.N2) = get_population_size(z)
    A_mul_B!(P.Pre.λ_A,P.MixingMatrix,(P.Pre.I_A + P.ϵ)/P.Pre.N) #note the scaling by 1/total population size
    return nothing
end

function InPlaceExternalForceOfInfOnHouseholdCalculation!(z,P::HH_RSVModelParameters)
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


function InPlaceInternalForceOfInfOnHouseholdCalculation!(z,P::HH_RSVModelParameters)
    # #Within households
    # #Household type specific infection rates
    P.Pre.r1 = P.τ*CondAgeOfInfectedU1(1/P.α,CurrentSolidProtectionDuration)*(1-AvMaternalProtU1s(1/P.α,CurrentSolidProtectionDuration))*sum(P.ConvMatrix_InfRate1_H*z);
    P.Pre.r2 = sus_2*P.τ*P.Inf_Rate_AgivState*z
    P.Pre.r_h = [P.Pre.r1;P.Pre.r2]
    return nothing
end

function RescaleHouseholdToStates!(z,P::HH_RSVModelParameters)
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

function TransmissionIncidence!(z,P::HH_RSVModelParameters)
    InplaceStateToHouseholdConversion!(z,P)
    InPlaceHouseholdToAgeConversion!(z,P)
    InPlaceAgedepForceOfInfection!(z,P)
    InPlaceExternalForceOfInfOnHouseholdCalculation!(z,P)
    InPlaceInternalForceOfInfOnHouseholdCalculation!(z,P)
    RescaleHouseholdToStates!(z,P)
    return nothing
end

function set_HH_growthrate!(z,P::HH_RSVModelParameters)
    A_mul_B!(P.Pre.Num_H,P.ConvMatrix_N_H,z)
    for i = 1:d1
        P.Pre.HH_growth[i] = P.RateOfChangeOfHHs[N_vect[i]]*z[i]/P.Pre.Num_H[N_vect[i]]
    end
    return nothing
end

function set_matrices_for_year!(P,year)
    set_conditionalmatrices_thatchangeeachyear!(P,year) #This updates with each year (= season + 1) because of the changes in the conditional age distribution
    generate_replacement_rate!(P,year+1)
    set_annual_rate_of_change_of_household_size!(P,year)
    return nothing
end

function F(du,u,P,t)
    #Get the numbers of households in each state and the cumulative incidence
    q = @view du[1:M_a];
    p = @view du[(M_a+1):end];
    z = @view u[(M_a+1):end];
    # Seasonality --- check what season it is, and if season has changed then update
    current_season = find_season_from_time(t)
    set_HH_growthrate!(z,P)
    if ~(P.Pre.season_curr == current_season)
        P.Pre.season_curr = current_season
    end
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
    CondAgeDistrib = CondAgeOfInfectedU1(1/P.α,CurrentSolidProtectionDuration)
    RateOfU1Incidence = mean(P.Pre.S1_A)*mean(P.Pre.λ_A[P.U1_cats])
    for i = 1:M_a
        if P.U1_cats[i]
            q[i] = Seasonality*(RateOfU1Incidence*CondAgeDistrib[i] + P.Pre.r_h[i])*P.HR[i]
        else
            q[i] = Seasonality*(P.Pre.S_A[i]*P.Pre.λ_A[i] + P.Pre.r_h[i])*P.HR[i]
        end
    end
    # # # Vector field for change in household configurations
    p .=  P.Pre.HH_growth .+
        P.CRE_mat*z .+
        P.ReplacementRate_mat*z .+
        (1 - AvMaternalProtU1s(1/P.α,CurrentSolidProtectionDuration))*Seasonality*P.τ*P.IntInfU1_mat*z .+
        Seasonality*P.τ*P.IntInfO1_mat*z .+
        Seasonality*P.ExtInfU1_mat*P.Pre.λ_state_U1s .+
        Seasonality*P.ExtInfO1_mat*P.Pre.λ_state_O1s
    return nothing
end


function create_callback_set_for_annual_changes()
    cb_set = CallbackSet(ContinuousCallback((u,t,integrator) -> t - YearStartTimes[1],int -> set_matrices_for_year!(int.sol.prob.p,1)))
    for i = 2:(length(YearStartDates))
        cb_set = CallbackSet(cb_set,ContinuousCallback((u,t,integrator) -> t - YearStartTimes[i],int -> set_matrices_for_year!(int.sol.prob.p,i)))
    end
    return cb_set
end
cb_set_annual_changes = create_callback_set_for_annual_changes()


function construct_susceptible_HH_distribution_from_year(X::Matrix{Int64},d1::Int64)
    z_0 = zeros(d1);
    s1,s2 = size(X)
    for i = 1:s1
        for j = 1:s2
            N_C = i-1;
            N_A = j-1;
            row = [N_C,0,0,N_A,0,0]';
            G = ismember_row(States,row);
            if any(G)
                z_0[G] += Float64(X[i,j]);
            end
        end
    end
    return z_0
end

#Function for running the model from the past to get an equilibrium intial condition for fitting to data
function get_equilibrium_state(date_for_solving_until)
    #Create date 10 years in the past from 2000
    date_in_past = Date(Dates.year(date_for_solving_until)-10,1,1)
    t_0 = Float64((date_in_past - Date(2000,1,1)).value)
    t_1 = Float64((date_for_solving_until - Date(2000,1,1)).value)
    tspan_eq = (t_0,t_1)
    C_0 = zeros(M_a)
    z_0 = construct_susceptible_HH_distribution_from_year(N_HEachYear[min(max(1,find_year(date_in_past)),18)],d1)
    x_0 = [C_0;z_0]
    set_matrices_for_year!(P_ModelParams,find_year(date_in_past))
    prob_eq = ODEProblem(F,x_0,tspan_eq,P_ModelParams);
    sol_eq = solve(prob_eq,CVODE_BDF(linear_solver=:GMRES),callback = cb_set_annual_changes,abstol = 1e-2, reltol = 1e-3)
    x_0_new = sol_eq[:,end]
    x_0_new[x_0_new .< 0] = 0.
    return x_0_new,sol_eq
end

function rescale_endstate(x_end,season_num,N_H)
    # x_end = sol[(M_a+1):end,end]
    if season_num < length(N_H)
        for i = 1:MaxHouseholdSize
            for j = 0:min(MaxNumberOfU1s,MaxHouseholdSize-i)
                n_o1 = i
                n_u1 = j
                RelevantStates = find((N1_vec .== n_u1).*(N2_vec .== n_o1))
                NumOfHHs = sum(x_end[RelevantStates])
                RescaleFactor = N_H[season_num][n_u1+1,n_o1+1]/NumOfHHs
                x_end[RelevantStates] .*= RescaleFactor
            end
        end
    end
    return x_end
end


function sol_for_params(Θ::Vector{Float64})
    if any(Θ .< 0)
        return Inf
    else
        #Preliminary solution to find an initial condition for 01-09-2001 (beginning of season 2)
        #Adjust the variable ParameterizedFunction
        P_ModelParams.b_C = 0.
        P_ModelParams.b_S = Θ[1]
        P_ModelParams.b_A = Θ[2]
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
        return sol_rsv
    end
end
