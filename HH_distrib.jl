# Calculate the population turnover rates that give the desired stationary household distributions

using Roots,LinearAlgebra
function HH_distrib_of_U1s(μ::Real,η::Real,H::Int64,MaxNumU1s::Int64)
    π_H = zeros(H+1)
    for n_u1 = 0:min(MaxNumU1s,H-1)
        π_H[n_u1+1] = ((μ/η)^n_u1)*binomial(H,n_u1)
    end
    return normalize!(π_H,1)
end

mean_num_U1s_in_HH(μ::Real,η::Real,H::Int64,MaxNumU1s::Int64) = dot(collect(0:H),HH_distrib_of_U1s(μ,η,H,MaxNumU1s))

function find_replacement_rate_matching_mean(η::Float64,H::Int64,MN_U1s::Real,MaxNumU1s::Int64)
    F(μ) = mean_num_U1s_in_HH(μ,η,H,MaxNumU1s) - MN_U1s;
    return fzero(F,0.0000000000001,10.)
end

function get_empirical_distrib_of_u1s_for_each_HHSize(N_H::Vector{Array{Int64,2}},year::Int64,η::Real,MaxNumU1s::Int64)
    Household_size_numU1s = zeros(MaxHouseholdSize,MaxNumU1s+1)
    for i = 1:MaxHouseholdSize
        for j = max(i-MaxNumU1s,1):i
            n_O1s = j
            n_U1s = i - j
            Household_size_numU1s[i,n_U1s+1] = N_HEachYear[year][n_U1s+1,n_O1s+1]
        end
        Household_size_numU1s[i,:] = normalize!(Household_size_numU1s[i,:],1)
    end
    return Household_size_numU1s
end

function FindReplacementRatesForYear(N_H::Vector{Array{Int64,2}},year::Int64,η::Real,MaxNumU1s::Int64)
    Household_size_numO1s = get_empirical_distrib_of_u1s_for_each_HHSize(N_H,year,η,MaxNumU1s)
    MeanNumO1sPerHousehold = Household_size_numO1s*collect(0:MaxNumU1s)
    μ_H = zeros(MaxHouseholdSize)
    for i = 2:MaxHouseholdSize
        μ_H[i] = find_replacement_rate_matching_mean(η,i,MeanNumO1sPerHousehold[i],MaxNumU1s::Int64)
    end
    return μ_H
end
FindReplacementRatesForYear(N_HEachYear,1,η_1,MaxNumberOfU1s)/
    maximum(FindReplacementRatesForYear(N_HEachYear,1,η_1,MaxNumberOfU1s))

function find_numbers_of_HH_by_size(N_H::Vector{Array{Int64,2}},MaxHHSize::Int64,MaxU1s::Int64)
    N_H_size = [zeros(Int64,MaxHHSize) for i = 1:length(N_H)]
    for year = 1:length(N_H)
        X = N_H[year]
        d1,d2 = size(X)
        for i = 1:d1
            for j = 2:d2
                if i+j-2 <= MaxHHSize && i-1 <= MaxU1s
                    N_H_size[year][i+j-2] += X[i,j]
                end
            end
        end
    end
    return N_H_size
end

function set_annual_rate_of_change_of_household_size!(P::HH_RSVModelParameters,year::Int64)
    if year>= 1 && year<= (length(N_HEachYear)-1)
        N_size = find_numbers_of_HH_by_size(N_HEachYear,MaxHouseholdSize,MaxNumberOfO1s)
        X = N_size[year + 1] - N_size[year]
        P.RateOfChangeOfHHs = X/365.25
        return nothing
    else
        P.RateOfChangeOfHHs = zeros(MaxHouseholdSize)
        return nothing
    end
end
set_annual_rate_of_change_of_household_size!(P_ModelParams,-1)
