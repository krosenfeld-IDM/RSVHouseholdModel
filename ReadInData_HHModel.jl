# Read in data and organise incidence into age category and weeks
using JLD,Dates
YearlyJointDistributions = load("input/EachYearJointDistribs_new.jld","people_distribution")

DailyRSV = load("input/DailyIFATbyAgeCat.jld","DailyRSVByAgeCat")
days = load("input/DailyIFATbyAgeCat.jld","x")
N_HEachYear = load("input/N_HEachYearvs3.jld","N_HEachYear")

for year = 1:length(N_HEachYear)
    N_HEachYear[year][:,2] += N_HEachYear[year][:,1]
    N_HEachYear[year][:,1] .= 0
end
# Aggregate the RSV case data
movsum = function(x::Vector{Float64},DT::Int64)
  y = zeros(size(x))
  for i = 1:length(x)
    if i <= DT
      y[i] = sum(x[1:i])
    else
      y[i] = sum(x[(i-DT+1):i])
    end
  end
  return y
end

AggregateRSVData = function(DailyRSV::Vector{Array{Float64,1}},days::StepRange{Date,Base.Dates.Day},L::Base.Dates.Day)
  #Find first day of data
    S = sum(DailyRSV)
    FirstDayData = days[findfirst(S .> 0)]
    #aggregate data
    SummedDailyRSV = [ movsum(DailyRSV[i],L.value) for i = 1:16 ]
    SummedDays = (days[1]+L):lookback:days[end]
    F_logical = [days[i] ∈ SummedDays for i = 1:length(days)]
    G_logical = days .> FirstDayData
    H_logical = F_logical.*G_logical
    SummedDailyRSV = [SummedDailyRSV[i][H_logical] for i = 1:16]
    pred_times = [Float64((days[H_logical][i] - Date(2000,1,1)).value) for i = 1:sum(H_logical)]
    return SummedDailyRSV,pred_times
end

Agg_RSVData,pred_times = AggregateRSVData(DailyRSV,days,lookback)
pred_dates = Date(2000,1,1) + Dates.Day.(pred_times)

#Hospitalisation estimates
function SetHospitalisationRate!(P::HH_RSVModelParameters)
  HR1 = [32.76*31.2,33.07*31.2,28.6*21.9,28.6*20.74,28.6*18.86,20.0*12.27,20.0*9.4,20.0*10.76,13.0*9.1,13.0*12.11,13.0*9.87,7.6*6.7];
  HR1 = HR1/10000
  HR2 = mean([([7.11,7.78,7.34,4.13,4.1])*(7.6/10000); ([10.52,20,14.84,13.95,10.18,2.36,8.41])*(2./10000)])
  HR2 = [HR2;[3.76,1.08,0.19]*(2.0/10000)]
  HR3 = zeros(14)
  P.HR = [HR1;HR2;HR3]
  return nothing
end
SetHospitalisationRate!(P_ModelParams)
# Fitted function for ratio of residents to non-residents
function Ratio(t::Float64)
    t -= Float64((Date(2002,4,22) - Date(2000,1,1)).value) # shift to the fitting time of the function
    if t > 0
        return 1.2409370942399072 +
        0.002241544326043194 * t +
        -2.4542080692642405e-6 * t^2 +
        9.446287443394756e-10 * t^3 +
        -1.552135315860274e-13 * t^4 +
        9.097449985032835e-18 * t^5
    else
        return 1.
    end
end
