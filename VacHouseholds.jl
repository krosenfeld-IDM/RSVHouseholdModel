#Alterations to the flow between household configurations due to vaccination

V_HH_S2_out = Dict()
for i = 1:d1
    StatesOut = States[i,:] + [1,0,0,-States[i,4],0,States[i,4]-1]
    if any(ismember_row(States,StatesOut'))
        k = findfirst(!iszero, ismember_row(States,StatesOut'))
        push!(V_HH_S2_out,i => k)
    end
end

V_HH_I2_out = Dict()
for i = 1:d1
    StatesOut = States[i,:] + [1,0,0,-States[i,4],-1,States[i,4]]
    if any(ismember_row(States,StatesOut'))
        k = findfirst(!iszero, ismember_row(States,StatesOut'))
        push!(V_HH_I2_out,i => k)
    end
end

V_HH_R2_out = copy(V_HH_S2_out)

D_HH_S2_out = Dict()
for i = 1:d1
    StatesOut = States[i,:] + [1,0,0,-1,0,0]
    if any(ismember_row(States,StatesOut'))
        k = findfirst(!iszero,ismember_row(States,StatesOut'))
        push!(D_HH_S2_out,i => k)
    end
end

D_HH_I2_out = Dict()
for i = 1:d1
    StatesOut = States[i,:] + [1,0,0,0,-1,0]
    if any(ismember_row(States,StatesOut'))
        k = findfirst(!iszero,ismember_row(States,StatesOut'))
        push!(D_HH_I2_out,i => k)
    end
end

D_HH_R2_out = Dict()
for i = 1:d1
    StatesOut = States[i,:] + [1,0,0,0,0,-1]
    if any(ismember_row(States,StatesOut'))
        k = findfirst(!iszero,ismember_row(States,StatesOut'))
        push!(D_HH_R2_out,i => k)
    end
end

function generate_replacement_rate!(P::HH_RSV_VaccinationModelParameters,year)
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
    ReplacementRate_mat = zeros(Float64,d1,d1)
    for k = 1:d1
        ReplacementRate_mat[k,k] += -μ_H[N_vect[k]]*(States[k,4]+States[k,5]+States[k,6])*(N2_vec[k] > max(1,N_vect[k] - MaxNumberOfU1s) );
        if get(D_HH_S2_out,k,false) != false
            ReplacementRate_mat[k,D_HH_S2_out[k]] += (1-HHCoverage)*μ_H[N_vect[k]]*States[k,4]*(N2_vec[k] > 1);
        end
        if get(D_HH_I2_out,k,false) != false
            ReplacementRate_mat[k,D_HH_I2_out[k]] += (1-HHCoverage)*μ_H[N_vect[k]]*States[k,5]*(N2_vec[k] > 1);
        end
        if get(D_HH_R2_out,k,false) != false
            ReplacementRate_mat[k,D_HH_R2_out[k]] += (1-HHCoverage)*μ_H[N_vect[k]]*States[k,6]*(N2_vec[k] > 1);
        end
        if get(V_HH_S2_out,k,false) != false
            ReplacementRate_mat[k,V_HH_S2_out[k]] += (HHCoverage)*μ_H[N_vect[k]]*States[k,4]*(N2_vec[k] > 1);
        end
        if get(V_HH_I2_out,k,false) != false
            ReplacementRate_mat[k,V_HH_I2_out[k]] += (HHCoverage)*μ_H[N_vect[k]]*States[k,5]*(N2_vec[k] > 1);
        end
        if get(V_HH_R2_out,k,false) != false
            ReplacementRate_mat[k,V_HH_R2_out[k]] += (HHCoverage)*μ_H[N_vect[k]]*States[k,6]*(N2_vec[k] > 1);
        end
    end
    P.ReplacementRate_mat = sparse(ReplacementRate_mat')
    return nothing
end
generate_replacement_rate!(P_VacModel,100)
