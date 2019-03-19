#Defines how the different states transition into one another

#Define the fast sparse conversion matrices for going from household states to household sizes
function AddConvMatricesToParams!(P::HH_RSVModelParameters)
    ConvMatrix_I1_H = zeros(MaxHouseholdSize,d1)
    ConvMatrix_S1_H = zeros(MaxHouseholdSize,d1)
    ConvMatrix_S2_0_H = zeros(MaxHouseholdSize,d1)
    ConvMatrix_S2_1_H = zeros(MaxHouseholdSize,d1)
    ConvMatrix_I2_0_H = zeros(MaxHouseholdSize,d1)
    ConvMatrix_I2_1_H = zeros(MaxHouseholdSize,d1)
    ConvMatrix_N1 = zeros(MaxHouseholdSize,d1)
    ConvMatrix_N2_0 = zeros(MaxHouseholdSize,d1)
    ConvMatrix_N2_1 = zeros(MaxHouseholdSize,d1)
    ConvMatrix_N_H = zeros(MaxHouseholdSize,d1)

    for i = 1:MaxHouseholdSize
        for j = 1:d1
            if N_vect[j] == i
                ConvMatrix_I1_H[i,j] = States[j,2];
                ConvMatrix_S1_H[i,j] = States[j,1];
                ConvMatrix_S2_0_H[i,j] = States[j,4]*(~U1_state[j]);
                ConvMatrix_S2_1_H[i,j] = States[j,4]*(U1_state[j]);
                ConvMatrix_I2_0_H[i,j] = States[j,5]*(~U1_state[j]);
                ConvMatrix_I2_1_H[i,j] = States[j,5]*(U1_state[j]);
                ConvMatrix_N1[i,j] = N1_vec[j];
                ConvMatrix_N2_0[i,j] = N2_vec[j]*(~U1_state[j])
                ConvMatrix_N2_1[i,j] = N2_vec[j]*(U1_state[j])
                ConvMatrix_N_H[i,j] = 1.
            end
        end
    end
    ConvMatrix_I1_H = sparse(ConvMatrix_I1_H)
    ConvMatrix_S1_H = sparse(ConvMatrix_S1_H)
    ConvMatrix_S2_0_H = sparse(ConvMatrix_S2_0_H)
    ConvMatrix_S2_1_H = sparse(ConvMatrix_S2_1_H)
    ConvMatrix_I2_0_H = sparse(ConvMatrix_I2_0_H)
    ConvMatrix_I2_1_H = sparse(ConvMatrix_I2_1_H)
    ConvMatrix_N_H = sparse(ConvMatrix_N_H)
    ConvMatrix_N1 = sparse(ConvMatrix_N1)
    ConvMatrix_N2_0 = sparse(ConvMatrix_N2_0)
    ConvMatrix_N2_1 = sparse(ConvMatrix_N2_1)
    P.ConvMatrix_S1_H = ConvMatrix_S1_H
    P.ConvMatrix_I1_H = ConvMatrix_I1_H
    P.ConvMatrix_N1 = ConvMatrix_N1
    P.ConvMatrix_N2_0 = ConvMatrix_N2_0
    P.ConvMatrix_N2_1 = ConvMatrix_N2_1
    P.ConvMatrix_N_H = ConvMatrix_N_H
    return nothing
end
AddConvMatricesToParams!(P_ModelParams)
#Find the states that flow into each state
# Declare useful ismember function
ismember_row = function(mat,row)
  s1,s2 = size(mat)
  B = Bool[mat[i,:]' == row for i = 1:s1]
end

#Generate the household changes matrices

function GenerateConstRateEventMatrix!(P::HH_RSVModelParameters)

  # Aging
  is_AgingS1_in = Bool[any(ismember_row(States,States[i,:]' + [1,0,0,-1,0,0]')) for i = 1:d1]
  AgingS1_in = zeros(Int64,d1)
  for i = 1:d1
    if is_AgingS1_in[i]
    AgingS1_in[i]  = (find(ismember_row(States,States[i,:]' + [1,0,0,-1,0,0]'))[1])
    end
  end

  is_AgingI1_in = Bool[any(ismember_row(States,States[i,:]' + [0,1,0,0,-1,0]')) for i = 1:d1]
  AgingI1_in = zeros(Int64,d1)
  for i = 1:d1
    if is_AgingI1_in[i]
    AgingI1_in[i]  = (find(ismember_row(States,States[i,:]' + [0,1,0,0,-1,0]'))[1])
    end
  end

  is_AgingR1_in = Bool[any(ismember_row(States,States[i,:]' + [0,0,1,0,0,-1]')) for i = 1:d1]
  AgingR1_in = zeros(Int64,d1)
  for i = 1:d1
    if is_AgingR1_in[i]
    AgingR1_in[i]  = (find(ismember_row(States,States[i,:]' + [0,0,1,0,0,-1]'))[1])
    end
  end


    # Recovery
    is_RecR1_in = Bool[any(ismember_row(States,States[i,:]' + [0,1,-1,0,0,0]')) for i = 1:d1]
    RecR1_in = zeros(Int64,d1)
    for i = 1:d1
      if is_RecR1_in[i]
      RecR1_in[i]  = (find(ismember_row(States,States[i,:]' + [0,1,-1,0,0,0]'))[1])
      end
    end

    is_RecR2_in = Bool[any(ismember_row(States,States[i,:]' + [0,0,0,0,1,-1]')) for i = 1:d1]
    RecR2_in = zeros(Int64,d1)
    for i = 1:d1
      if is_RecR2_in[i]
      RecR2_in[i]  = (find(ismember_row(States,States[i,:]' + [0,0,0,0,1,-1]'))[1])
      end
    end

    #Reversion
    is_RevR1_in = Bool[any(ismember_row(States,States[i,:]' + [-1,0,1,0,0,0]')) for i = 1:d1]
    RevR1_in = zeros(Int64,d1)
    for i = 1:d1
      if is_RevR1_in[i]
      RevR1_in[i]  = (find(ismember_row(States,States[i,:]' + [-1,0,1,0,0,0]'))[1])
      end
    end

    is_RevR2_in = Bool[any(ismember_row(States,States[i,:]' + [0,0,0,-1,0,1]')) for i = 1:d1]
    RevR2_in = zeros(Int64,d1)
    for i = 1:d1
      if is_RevR2_in[i]
      RevR2_in[i]  = (find(ismember_row(States,States[i,:]' + [0,0,0,-1,0,1]'))[1])
      end
    end
    #Constant rate events
    CRE_mat = zeros(Float64,d1,d1)
    for k = 1:d1
      CRE_mat[k,k] += - γ_1*States[k,2] - γ_2*States[k,5] - ν*(States[k,3] + States[k,6]) - η_1*(States[k,1] +States[k,2] + States[k,3]);
      if is_RecR1_in[k]
        CRE_mat[k,RecR1_in[k]] += γ_1*States[RecR1_in[k],2];
      end
      if is_RecR2_in[k]
        CRE_mat[k,RecR2_in[k]] += γ_2*States[RecR2_in[k],5];
      end
      if is_RevR1_in[k]
        CRE_mat[k,RevR1_in[k]] += ν*States[RevR1_in[k],3];
      end
      if is_RevR2_in[k]
        CRE_mat[k,RevR2_in[k]] += ν*States[RevR2_in[k],6];
      end
      if is_AgingS1_in[k]
        CRE_mat[k,AgingS1_in[k]] += η_1*States[AgingS1_in[k],1];
      end
      if is_AgingI1_in[k]
        CRE_mat[k,AgingI1_in[k]] += η_1*States[AgingI1_in[k],2];
      end
      if is_AgingR1_in[k]
        CRE_mat[k,AgingR1_in[k]] += η_1*States[AgingR1_in[k],3];
      end
    end
    P.CRE_mat = sparse(CRE_mat)
    return nothing
end
GenerateConstRateEventMatrix!(P_ModelParams)

# Infection
function GenerateInfectionEventMatrices!(P::HH_RSVModelParameters)
    #Houshold size effects
    InvEffHHSize = (max.(1,N_vect-1)).^(-P.EffHHSizePower);
    HHInfRate1 = (States[:,2] + inf_2*States[:,5]).*InvEffHHSize;
    HHInfRate2 = (States[:,2] + inf_2*States[:,5]).*InvEffHHSize;

    is_InfI1_in = Bool[any(ismember_row(States,States[i,:]' + [1,-1,0,0,0,0]')) for i = 1:d1]
    InfI1_in = zeros(Int64,d1)
    for i = 1:d1
      if is_InfI1_in[i]
      InfI1_in[i]  = (find(ismember_row(States,States[i,:]' + [1,-1,0,0,0,0]'))[1])
      end
    end

    is_InfI2_in = Bool[any(ismember_row(States,States[i,:]' + [0,0,0,1,-1,0]')) for i = 1:d1]
    InfI2_in = zeros(Int64,d1)
    for i = 1:d1
      if is_InfI2_in[i]
      InfI2_in[i]  = (find(ismember_row(States,States[i,:]' + [0,0,0,1,-1,0]'))[1])
      end
    end

    IntInfU1_mat = zeros(Float64,d1,d1)
    for k = 1:d1
      IntInfU1_mat[k,k] += -InvEffHHSize[k]*States[k,1]*(States[k,2] + inf_2*States[k,5]);
      if is_InfI1_in[k]
        IntInfU1_mat[k,InfI1_in[k]] += InvEffHHSize[InfI1_in[k]]*States[InfI1_in[k],1]*(States[InfI1_in[k],2] + inf_2*States[InfI1_in[k],5]);
      end
    end

    IntInfO1_mat = zeros(Float64,d1,d1)
    for k = 1:d1
      IntInfO1_mat[k,k] += -InvEffHHSize[k]*sus_2*States[k,4]*(States[k,2] + inf_2*States[k,5]);
      if is_InfI2_in[k]
        IntInfO1_mat[k,InfI2_in[k]] += InvEffHHSize[InfI2_in[k]]*sus_2*States[InfI2_in[k],4]*(States[InfI2_in[k],2] + inf_2*States[InfI2_in[k],5]);
      end
    end

    ExtInfU1_mat = zeros(Float64,d1,d1)
    for k = 1:d1
      ExtInfU1_mat[k,k] += -States[k,1];
      if is_InfI1_in[k]
        ExtInfU1_mat[k,InfI1_in[k]] += States[InfI1_in[k],1];
      end
    end

    ExtInfO1_mat = zeros(Float64,d1,d1)
    for k = 1:d1
      ExtInfO1_mat[k,k] += -sus_2*States[k,4];
      if is_InfI2_in[k]
        ExtInfO1_mat[k,InfI2_in[k]] += sus_2*States[InfI2_in[k],4];
      end
    end
    # P.IntInf_mat = sparse((1 - AvMaternalProtU1s(1/P.α))*IntInfU1_mat + IntInfO1_mat)
    P.IntInfU1_mat = sparse(IntInfU1_mat)
    P.IntInfO1_mat = sparse(IntInfO1_mat)
    P.ExtInfU1_mat = sparse(ExtInfU1_mat)
    P.ExtInfO1_mat = sparse(ExtInfO1_mat)
    return nothing
end

GenerateInfectionEventMatrices!(P_ModelParams)


#Generate replacement rate
#Replacements --- replace an adult with a new susceptible baby
is_D_S2_in = Bool[any(ismember_row(States,States[i,:]' + [-1,0,0,1,0,0]')) for i = 1:d1]
D_S2_in = zeros(Int64,d1)
for i = 1:d1
  if is_D_S2_in[i]
  D_S2_in[i]  = (find(ismember_row(States,States[i,:]' + [-1,0,0,1,0,0]'))[1])
  end
end

is_D_I2_in = Bool[any(ismember_row(States,States[i,:]' + [-1,0,0,0,1,0]')) for i = 1:d1]
D_I2_in = zeros(Int64,d1)
for i = 1:d1
  if is_D_I2_in[i]
  D_I2_in[i]  = (find(ismember_row(States,States[i,:]' + [-1,0,0,0,1,0]'))[1])
  end
end

is_D_R2_in = Bool[any(ismember_row(States,States[i,:]' + [-1,0,0,0,0,1]')) for i = 1:d1]
D_R2_in = zeros(Int64,d1)
for i = 1:d1
  if is_D_R2_in[i]
  D_R2_in[i]  = (find(ismember_row(States,States[i,:]' + [-1,0,0,0,0,1]'))[1])
  end
end

function generate_replacement_rate!(P::HH_RSVModelParameters,year)
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
        if is_D_S2_in[k]
        ReplacementRate_mat[k,D_S2_in[k]] += μ_H[N_vect[D_S2_in[k]]]*States[D_S2_in[k],4]*(N2_vec[D_S2_in[k]] > 1);
        end
        if is_D_I2_in[k]
        ReplacementRate_mat[k,D_I2_in[k]] += μ_H[N_vect[D_I2_in[k]]]*States[D_I2_in[k],5]*(N2_vec[D_I2_in[k]] > 1);
        end
        if is_D_R2_in[k]
        ReplacementRate_mat[k,D_R2_in[k]] += μ_H[N_vect[D_R2_in[k]]]*States[D_R2_in[k],6]*(N2_vec[D_R2_in[k]] > 1);
        end
    end
    P.ReplacementRate_mat = sparse(ReplacementRate_mat)
    return nothing
end
generate_replacement_rate!(P_ModelParams,initial_year_index)
