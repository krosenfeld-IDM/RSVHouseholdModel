#Extra parameters and definitions for vaccination forecasting models

# Define the endpoints of each future season for vaccination
FutureSeasonStartDates = [Date(y-1,9,1) for y = 2017:2036]
FutureSeasonEndDates = [Date(y,9,1) for y = 2017:2036]
FutureYearStartDates = [Date(y,1,1) for y = 2017:2036]
FutureYearEndDates = [Date(y,12,31) for y = 2017:2036]
FutureYearStartTimes = [Float64((d - Date(2000,1,1)).value) for d in FutureYearStartDates]
FutureSeasonStartTimes = [Float64((d - Date(2000,1,1)).value) for d in FutureSeasonStartDates]
FutureSeasonEndTimes = [Float64((d - Date(2000,1,1)).value) for d in FutureYearEndDates]

find_future_season(date_t) = min(searchsortedfirst(FutureSeasonEndDates,date_t),length(FutureSeasonEndDates))
find_future_season_from_time(t) = find_future_season(Date(2000,1,1)+Dates.Day(round(Int64,t)))

function NoMatProtection(a,M,SolidProtDuration)
    if a < SolidProtDuration
        return 0.
    else
        return 1 - exp(-(a-SolidProtDuration)/M);
    end
end

function CondAgeOfInfectedU1(M,SolidProtDuration)
    return normalize!(map(NoMatProtection,0.5*(MeshAgesStartPts[U1_cats]+MeshAgesEndPts[U1_cats]),M.*ones(M_a)[U1_cats],SolidProtDuration.*ones(M_a)[U1_cats]),1)
end

function MeanHR_U1(M,SolidProtDuration)
    return dot(CondAgeOfInfectedU1(M,SolidProtDuration),HR1)
end

mutable struct HH_RSV_VaccinationModelParameters
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

function create_RSV_vaccination_model(P::HH_RSVModelParameters)
    PreAllocatedVects = PreAllocated(Vector{Float64}(MaxHouseholdSize),
    Vector{Float64}(MaxHouseholdSize),
    Vector{Float64}(MaxHouseholdSize),
    Vector{Float64}(MaxHouseholdSize),
    Vector{Float64}(MaxHouseholdSize),
    Vector{Float64}(MaxHouseholdSize),
    Vector{Float64}(sum(U1_cats)),
    Vector{Float64}(sum(.~U1_cats)),
    Vector{Float64}(M_a), #Over 1s are less infectious
    Vector{Float64}(sum(U1_cats)),
    Vector{Float64}(sum(.~U1_cats)),
    Vector{Float64}(M_a),
    Vector{Float64}(sum(U1_cats)),
    Vector{Float64}(sum(.~U1_cats)),
    Vector{Float64}(M_a),
    Vector{Float64}(M_a),
    Vector{Float64}(sum(U1_cats)),
    Vector{Float64}(sum(.~U1_cats)),
    Vector{Float64}(sum(.~U1_cats)),
    Vector{Float64}(d1),
    Vector{Float64}(d1),
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
    Vector{Float64}(M_a),
    Vector{Float64}(M_a),
    Vector{Float64}(M_a),
    Vector{Float64}(d1+M_a),
    Float64(1),
    -99,
    -99,
    Vector{Float64}(d1),
    Float64(1))

    P_VacModel = HH_RSV_VaccinationModelParameters(
        P.b_C,
        P.b_S,
        P.b_A,
        P.β,
        P.β_S,
        P.β_C,
        P.τ,
        P.ξ̄,
        P.ϕ̄,
        P.σ_ξ,
        P.σ_ϕ,
        P.ρ_ξϕ,
        P.α,
        P.ϵ,
        P.EffHHSizePower,
        P.Χ,
        P.ξ̄*ones(length(FutureSeasonStartTimes)),
        P.ϕ̄*ones(length(FutureSeasonStartTimes)),
        P.HR,
        P.ConvMatrix_S1_H,
        P.ConvMatrix_I1_H,
        P.ConvMatrix_N1,
        P.ConvMatrix_N2_0,
        P.ConvMatrix_N2_1,
        P.ConvMatrix_InfRate1_H,
        P.ConvMatrix_I2_0_H,
        P.ConvMatrix_I2_1_H,
        P.ConvMatrix_S2_0_H,
        P.ConvMatrix_S2_1_H,
        P.ConvMatrix_InfRate2_0_H,
        P.ConvMatrix_InfRate2_1_H,
        P.ConvMatrix_N_H,
        P.Inf_AgivState,
        P.Sus_AgivState,
        P.N_AgivState,
        P.Inf_Rate_AgivState,
        P.U1s_HU_giv_A,
        P.O1s_HU_giv_A_0,
        P.O1s_HU_giv_A_1,
        P.MixingMatrix,
        P.CRE_mat,
        P.ReplacementRate_mat,
        P.IntInfU1_mat,
        P.IntInfO1_mat,
        P.ExtInfU1_mat,
        P.ExtInfO1_mat,
        P.U1_cats,
        P.C_cats,
        P.S_cats,
        P.A_cats,
        P.U1_state,
        PreAllocatedVects,
        zeros(MaxHouseholdSize) )

    return P_VacModel
end

P_VacModel = create_RSV_vaccination_model(P_ModelParams)

function Χ(t,P::HH_RSV_VaccinationModelParameters)#Time varying seasonality
    season = find_future_season_from_time(t)
    if season < 1
        season = 1
    end
    if season > length(P.Z_ξ)
        season = length(P.Z_ξ)
    end
    ξ = P.Z_ξ[season]
    ϕ = P.Z_ϕ[season]
    return exp(ξ*cospi(2*(t - ϕ)/365.25))
end
