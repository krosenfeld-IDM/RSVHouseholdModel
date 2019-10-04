#Basic model definitions
#Fundamental model set-up parameters

using Distributions,DifferentialEquations,Plots
lookback = Dates.Day(7);
MaxHouseholdSize = 10;
MaxNumberOfO1s = 10;
MaxNumberOfU1s = 2;
initial_year_index = 1;

#basic dynamic parameters
const γ_1 = (1.0/9.)
const γ_2 = 0.25
const ν = (1.0/180.0)
const sus_2 = 0.75
const inf_2 = 0.5
const η_1 = 1/365.25;
CurrentSolidProtectionDuration = 0.
SchoolsOnlyR_0 = 1.5
SchoolsContactRate = SchoolsOnlyR_0*γ_2/(sus_2*inf_2)
