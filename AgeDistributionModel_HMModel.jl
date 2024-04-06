#Use the KHDSS data to construct conditional age distributions

function SetupConditionalDistributions(YearlyJointDistributions::Vector{Array{Float64,3}})
    P_AHU = [YearlyJointDistributions[i][:,1:MaxHouseholdSize,:]/sum(YearlyJointDistributions[i][:,1:MaxHouseholdSize,:]) for i = 1:18]
    P_AgivHU = [zeros(size(P_AHU[1])) for i = 1:18]
    P_HUgivA = [zeros(size(P_AHU[1])) for i = 1:18]
    P_UgivH = [zeros(size(P_AHU[1][1,:,:])) for i = 1:18]


    P_HU = [sum(P_AHU[n],dims=1) for n = 1:18]
    P_A = [sum(P_AHU[n],dims=2:3) for n = 1:18]
    P_H = [sum(P_AHU[n],dims=[1,3]) for n = 1:18]
    for n = 1:18, j = 1:MaxHouseholdSize, k = 1:2
        P_UgivH[n][j,k] = P_HU[n][1,j,k]/P_H[n][1,j,1]
    end
    #Calculate conditional distributions
    for n = 1:18
        for i = 1:30
            for j = 1:MaxHouseholdSize
                for k = 1:2
                    P_AgivHU[n][i,j,k] = P_AHU[n][i,j,k]/P_HU[n][1,j,k]
                    P_HUgivA[n][i,j,k] = P_AHU[n][i,j,k]/P_A[n][i,1,1]
                    if isnan(P_AgivHU[n][i,j,k])
                        P_AgivHU[n][i,j,k] = 0.
                    end
                    if isnan(P_HUgivA[n][i,j,k])
                        P_HUgivA[n][i,j,k] = 0.
                    end
                end
            end
        end
    end
    return P_AgivHU,P_HUgivA ,P_HU,P_A,P_H,P_UgivH
end

P_AgivHU,P_HUgivA,P_HU,P_A,P_H,P_UgivH =SetupConditionalDistributions(YearlyJointDistributions)
