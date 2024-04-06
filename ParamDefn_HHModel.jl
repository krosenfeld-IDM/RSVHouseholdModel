
# Define the Household model --- all parameters relevant to simuation are collected in a struct: HH_RSVModelParameters
# Declare Age mesh for background demography model
# Declare the season start and finish dates and the year start and finish dates
# Declare the conditional age distribution functions

using Base.MathConstants
using SparseArrays
using LinearAlgebra

#Age mesh and categories
MeshAgesEndPts = [[Float64(i) for i = range(30.4,365.25,12)]
;[Float64(i) for i = range(2*365.25,365.25*17,16)]
;18*365.25;150*365.25]
MeshAgesStartPts = [0;MeshAgesEndPts[1:(end-1)]]

U1_cats = MeshAgesEndPts .<= 365.25*1
C_cats = (MeshAgesEndPts .< 365.25*5)
S_cats = (MeshAgesEndPts .>= 365.25*5).&(MeshAgesEndPts .< 365.25*18)
A_cats = (MeshAgesEndPts .>= 365.25*18)
O1_cats = MeshAgesEndPts .> 365.25
M_a = length(MeshAgesEndPts) # Number of age nodes

# Define the start and end points of each season
SeasonStartDates = [Date(y-1,9,1) for y = 2001:2017]
SeasonEndDates = [Date(y,9,1) for y = 2001:2017]
YearStartDates = [Date(y,1,1) for y = 2000:2017]
YearEndDates = [Date(y,12,31) for y = 2000:2017]
YearStartTimes = [Float64((d - Date(2000,1,1)).value) for d in YearStartDates]
SeasonStartTimes = [Float64((d - Date(2000,1,1)).value) for d in SeasonStartDates]
SeasonEndTimes = [Float64((d - Date(2000,1,1)).value) for d in SeasonEndDates]




# @time searchsortedfirst(SeasonEndDates,Date(2020,1,1))

FindSeason(date_t) = searchsortedfirst(SeasonEndDates,date_t)
find_year(date_t) = Dates.Year(date_t).value - 1999
find_season_from_time(t) = FindSeason(Date(2000,1,1)+Dates.Day(round(Int64,t)))
find_year_from_time(t) = find_year(Date(2000,1,1)+Dates.Day(round(Int64,t)))
time_to_date(t) = Date(2000,1,1)+Dates.Day(round(Int64,t))
find_actual_year_from_time(t) = Dates.Year(time_to_date(t)).value

function Χ(t,Z_ξ::Vector{Float64},Z_ϕ::Vector{Float64})#Time varying seasonality
    date_t = Date(2000,1,1)+Dates.Day(round(Int64,t))
    season = FindSeason(date_t)
    ξ = Z_ξ[season]
    ϕ = Z_ϕ[season]
    # return 1/(1 + exp(-ξ*cospi(2*(t - ϕ)/365.25)))
    return exp(ξ*cospi(2*(t - ϕ)/365.25))
end


#Maternal protection functions

function NoMatProtection(a,M)
    L = 365.25
    if a <= L
        return (1 - exp(-a/M))/(1-exp(-L/M))
    else
        return 1.
    end
end

function NoMatProtection(a,M,S)
    L = 365.25
    if a <= S
        return 0.
    end
    if a > S && a<=L
        return (1 - exp(-(a-S)/M))/(1-exp(-(L-S)/M))
    end
    if a > L
        return 1.
    end
end

function NoMatProtection(a,M,S,coverage)
    return coverage*NoMatProtection(a,M,S) + (1-coverage)*NoMatProtection(a,M,0.)
end

function MatProtection(a,M)
    return 1 - NoMatProtection(a,M)
end

function MatProtection(a,M,S)
    return 1 - NoMatProtection(a,M,S)
end

function MatProtection(a,M,S,coverage)
    return 1 - NoMatProtection(a,M,S,coverage)
end

function AvMaternalProtU1s(M)
    return (M/365.25)
end

exp_ratio_func(x) = e^x/(1-e^x)

function AvMaternalProtU1s(M,S) #Form if mat. prot. is  conditioned on ending at age 1
    T = 365.25
    return (S + M - (T-S)*exp_ratio_func(-(T-S)/M) )/T
end


function AvMaternalProtU1s(M,S,coverage) #Form if mat. prot. is  conditioned on ending at age 1
    return (1-coverage)*AvMaternalProtU1s(M,0) + coverage*AvMaternalProtU1s(M,S)
end

#Conditional age distributions
function prob_of_being_susceptible_conditional_on_age_interval(M,S,a0,a1)
    return 1 - ((M/(a1-a0))*(exp(-a0/M) - exp(-a1/M)))
end

function CondAgeOfInfectedU1(M)
    return [prob_of_being_in_age_int_conditional_on_being_infected(M,0.,MeshAgesStartPts[i],MeshAgesEndPts[i]) for i = 1:sum(U1_cats) ]
end

exp_ratio_func2(x,y) = e^x/(1-e^y)

function prob_integral_func(M,S,a0,a1)
    T = 365.25
    x = max(a0-S,0.)
    y = max(a1-S,0.)
    return ((y-x)/(1-exp(-(T-S)/M))) + M*exp_ratio_func2(-y/M,-T/M) -M*exp_ratio_func2(-x/M,-T/M)
end


function prob_of_being_in_age_int_conditional_on_being_infected(M,S,a0,a1,coverage)
    T = 365.25
    return ((1-coverage)*prob_integral_func(M,0.,a0,a1) + coverage*prob_integral_func(M,S,a0,a1))/(T*(1-AvMaternalProtU1s(M,S,coverage)))
end

function prob_of_being_in_age_int_conditional_on_being_infected(M,S,a0,a1)
    T = 365.25
    return prob_of_being_in_age_int_conditional_on_being_infected(M,S,a0,a1,0.)
end

function CondAgeOfInfectedU1(M)
    return [prob_of_being_in_age_int_conditional_on_being_infected(M,0.,MeshAgesStartPts[i],MeshAgesEndPts[i]) for i = 1:sum(U1_cats) ]
end

function CondAgeOfInfectedU1(M,S)
    return [prob_of_being_in_age_int_conditional_on_being_infected(M,S,MeshAgesStartPts[i],MeshAgesEndPts[i]) for i = 1:sum(U1_cats) ]
end

function CondAgeOfInfectedU1(M,S,coverage)
    return [prob_of_being_in_age_int_conditional_on_being_infected(M,S,MeshAgesStartPts[i],MeshAgesEndPts[i],coverage) for i = 1:sum(U1_cats) ]
end

CondAgeDistrib = CondAgeOfInfectedU1(19.,30,1) #Just a placeholder
RateOfU1Incidence = 0.

function MeanHR_U1(M,S)
    return dot(CondAgeOfInfectedU1(M,S),HR1)
end

#Inplace allocation
mutable struct PreAllocated
    S1_H::Vector{Float64}
    I1_H::Vector{Float64}
    Num_H::Vector{Float64}
    N1_H::Vector{Float64}
    N2_0_H::Vector{Float64}
    N2_1_H::Vector{Float64}
    I1_A::Vector{Float64}
    I2_A::Vector{Float64}
    I_A::Vector{Float64} #Over 1s are less infectious
    S1_A::Vector{Float64}
    S2_A::Vector{Float64}
    S_A::Vector{Float64}
    N1_A::Vector{Float64}
    N2_A::Vector{Float64}
    N_A::Vector{Float64}
    λ_A::Vector{Float64}
    λ_H_U1s::Vector{Float64}
    λ_H_0_O1s::Vector{Float64}
    λ_H_1_O1s::Vector{Float64}
    λ_state_U1s::Vector{Float64}
    λ_state_O1s::Vector{Float64}
    N::Float64
    N1::Float64
    N2::Float64
    NC::Float64
    NS::Float64
    NA::Float64
    I::Float64
    I1::Float64
    IC::Float64
    IS::Float64
    IA::Float64
    I2::Float64
    S::Float64
    S1::Float64
    SC::Float64
    SS::Float64
    SA::Float64
    S2::Float64
    r1::Vector{Float64}
    r2::Vector{Float64}
    r_h::Vector{Float64}
    x_curr::Vector{Float64}
    t_curr::Float64
    season_curr::Int64
    year_curr::Int64
    HH_growth::Vector{Float64}
    weighted_inf::Float64
end



mutable struct HH_RSVModelParameters
    b_C::Float64
    b_S::Float64
    b_A::Float64
    β::Float64
    β_S::Float64
    β_C::Float64
    τ::Float64
    ξ̄::Float64
    ϕ̄::Float64
    σ_ξ::Float64
    σ_ϕ::Float64
    ρ_ξϕ::Float64
    α::Float64
    ϵ::Float64
    EffHHSizePower::Float64;
    Χ::Function
    Z_ξ::Vector{Float64}
    Z_ϕ::Vector{Float64}
    HR::Vector{Float64}
    ConvMatrix_S1_H::SparseMatrixCSC{Float64,Int64}
    ConvMatrix_I1_H::SparseMatrixCSC{Float64,Int64}
    ConvMatrix_N1::SparseMatrixCSC{Float64,Int64}
    ConvMatrix_N2_0::SparseMatrixCSC{Float64,Int64}
    ConvMatrix_N2_1::SparseMatrixCSC{Float64,Int64}
    ConvMatrix_InfRate1_H::SparseMatrixCSC{Float64,Int64}
    ConvMatrix_I2_0_H::SparseMatrixCSC{Float64,Int64}
    ConvMatrix_I2_1_H::SparseMatrixCSC{Float64,Int64}
    ConvMatrix_S2_0_H::SparseMatrixCSC{Float64,Int64}
    ConvMatrix_S2_1_H::SparseMatrixCSC{Float64,Int64}
    ConvMatrix_InfRate2_0_H::SparseMatrixCSC{Float64,Int64}
    ConvMatrix_InfRate2_1_H::SparseMatrixCSC{Float64,Int64}
    ConvMatrix_N_H::SparseMatrixCSC{Float64,Int64}
    Inf_AgivState::Array{Float64,2}
    Sus_AgivState::Array{Float64,2}
    N_AgivState::Array{Float64,2}
    Inf_Rate_AgivState::Array{Float64,2}
    U1s_HU_giv_A::Array{Float64,2}
    O1s_HU_giv_A_0::Array{Float64,2}
    O1s_HU_giv_A_1::Array{Float64,2}
    MixingMatrix::Array{Float64,2}
    CRE_mat::SparseMatrixCSC{Float64,Int64}
    ReplacementRate_mat::SparseMatrixCSC{Float64,Int64}
    IntInfU1_mat::SparseMatrixCSC{Float64,Int64}
    IntInfO1_mat::SparseMatrixCSC{Float64,Int64}
    ExtInfU1_mat::SparseMatrixCSC{Float64,Int64}
    ExtInfO1_mat::SparseMatrixCSC{Float64,Int64}
    U1_cats::BitArray{1}
    C_cats::BitArray{1}
    S_cats::BitArray{1}
    A_cats::BitArray{1}
    U1_state::BitArray{2}
    Pre::PreAllocated
    # RateOfChangeOfHHs::Vector{Array{Float64,1}}
    RateOfChangeOfHHs::Vector{Float64}
end


function create_RSV_model()
    PreAllocatedVects = PreAllocated(Vector{Float64}(undef, MaxHouseholdSize),
    Vector{Float64}(undef, MaxHouseholdSize),
    Vector{Float64}(undef, MaxHouseholdSize),
    Vector{Float64}(undef, MaxHouseholdSize),
    Vector{Float64}(undef, MaxHouseholdSize),
    Vector{Float64}(undef, MaxHouseholdSize),
    Vector{Float64}(undef, sum(U1_cats)),
    Vector{Float64}(undef, sum(.~U1_cats)),
    Vector{Float64}(undef, M_a), #Over 1s are less infectious
    Vector{Float64}(undef, sum(U1_cats)),
    Vector{Float64}(undef, sum(.~U1_cats)),
    Vector{Float64}(undef, M_a),
    Vector{Float64}(undef, sum(U1_cats)),
    Vector{Float64}(undef, sum(.~U1_cats)),
    Vector{Float64}(undef, M_a),
    Vector{Float64}(undef, M_a),
    Vector{Float64}(undef, sum(U1_cats)),
    Vector{Float64}(undef, sum(.~U1_cats)),
    Vector{Float64}(undef, sum(.~U1_cats)),
    Vector{Float64}(undef, d1),
    Vector{Float64}(undef, d1),
    Float64(1),
    Float64(1),
    Float64(1),
    Float64(1),
    Float64(1),
    Float64(1),
    Float64(1),
    Float64(1),
    Float64(1),
    Float64(1),
    Float64(1),
    Float64(1),
    Float64(1),
    Float64(1),
    Float64(1),
    Float64(1),
    Float64(1),
    Float64(1),
    Vector{Float64}(undef, M_a),
    Vector{Float64}(undef, M_a),
    Vector{Float64}(undef, M_a),
    Vector{Float64}(undef, d1+M_a),
    Float64(1),
    -99,
    -99,
    Vector{Float64}(undef, d1),
    Float64(1))

    I_z = I(10)
    I_z_sp = sparse(I_z)
    I_z_vect = [zeros(MaxHouseholdSize) for i=1:18]

    P_ModelParams = HH_RSVModelParameters(1.,1.,1.,1.,1.,1.,
    0.1,
    0.1,60.,1.,1.,0.,
    1/30.,10.,0.,
    t -> 1/(1 + exp(-0. * cospi(2*(t - 0.)/365.25))),
    zeros(18),zeros(18),zeros(30),
    I_z_sp,I_z_sp,I_z_sp,I_z_sp,I_z_sp,I_z_sp,
    I_z_sp,I_z_sp,I_z_sp,I_z_sp,I_z_sp,I_z_sp,I_z_sp,
    I_z,I_z,I_z,I_z,I_z,I_z,I_z,I_z,
    I_z_sp,I_z_sp,I_z_sp,I_z_sp,I_z_sp,I_z_sp,
    U1_cats,C_cats,S_cats,A_cats,
    U1_state,
    PreAllocatedVects,
    zeros(MaxHouseholdSize))

    return P_ModelParams
end

P_ModelParams = create_RSV_model()


function SetMixingMatrix!(P::HH_RSVModelParameters)
    U1Mixing = zeros(M_a,M_a);
    U1Mixing[U1_cats,:] .= 1.;
    U1Mixing[:,U1_cats] .= 1.;
    ChildMixing = zeros(M_a,M_a);
    ChildMixing[C_cats,C_cats] .= 1.;
    SchoolMixing = zeros(M_a,M_a);
    SchoolMixing[S_cats,S_cats] .= 1.;
    AdultMixing = zeros(M_a,M_a);
    AdultMixing[A_cats,A_cats] .= 1.;
    HomogeneousMixing = zeros(M_a,M_a)
    HomogeneousMixing[.~U1_cats,.~U1_cats] .= 1.;
    P.MixingMatrix = P.b_S*SchoolMixing .+ P.b_A*AdultMixing + P.β*HomogeneousMixing + P.b_C*U1Mixing;
    return nothing
end

function SetSeasonalityFunction!(P::HH_RSVModelParameters,season::Int64)
    ξ = P.Z_ξ[season]
    ϕ = P.Z_ϕ[season]
    # P.Χ = t -> 1/(1 + exp(-ξ*cospi(2*(t - ϕ)/365.25)))
    P.Χ = t -> exp(ξ*cospi(2*(t - ϕ)/365.25))
end

function Χ(t,P::HH_RSVModelParameters)#Time varying seasonality
    season = find_season_from_time(t)
    if season < 1
        season = 1
    end
    if season > length(P.Z_ξ)
        season = length(P.Z_ξ)
    end
    ξ = P.Z_ξ[season]
    ϕ = P.Z_ϕ[season]
    # return 1/(1 + exp(-ξ*cospi(2*(t - ϕ)/365.25)))
    return exp(ξ*cospi(2*(t - ϕ)/365.25))
end
