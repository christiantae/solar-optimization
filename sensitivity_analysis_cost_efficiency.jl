using JuMP
using Cbc
using CSV
using DataFrames
using AmplNLWriter

PATH_TO_SOLVERS = ENV["ERE291_SOLVERS"]

# Get data from csv files
HomePowerDataFrame = CSV.read("home_load_and_gen_hourly.csv"; delim=',', types=[Float64, Float64]);
# Energy consumed by the home [kWh]
EnergyUsed = HomePowerDataFrame[:1]; #averaged across all days of the year
#EnergyUsed = [0.104166667, 0.270833333, 0.4375, 0.604166667, 0.770833333, 0.9375, 1.104166667, 1.270833333, 1.4375, 1.604166667, 1.770833333, 1.9375, 2.104166667, 2.270833333, 2.4375, 2.604166667, 2.770833333, 2.9375, 3.104166667, 3.270833333, 3.4375, 3.604166667, 3.770833333, 2.9375]#kWh
#SunAngleDataFrame = CSV.read("cos_zenith_angle_vermont.csv"; delim=',', types=[Float64, Float64]);
#CosZenithAngle = SunAngleDataFrame[:2]; #currently getting the January 1st data
IrradianceDataFrame = CSV.read("solar_irradiance_vermont.csv"; delim=',', types=[Float64, Float64, Float64]);
#DiffuseHorizontalIrradiance = IrradianceDataFrame[:1]; #averaged across all days of the year
#DirectNormalIrradiance = IrradianceDataFrame[:2]; #averaged across all days of the year
#GlobalHorizontalIrradiance = DirectNormalIrradiance.*CosZenithAngle + DiffuseHorizontalIrradiance;
# The solar power hitting location at every timestep in the average day [W/m^2]
GlobalHorizontalIrradiance = IrradianceDataFrame[:3];
#GlobalHorizontalIrradiance = [0 0 0 0 0.705479452 11.59452055 34.25068493 64.72054795 106.6849315
#        147.3479452 177.0123288 191.1369863 196.3671233 184.1808219 158.0054795 117.8136986 79.77123288
#        43.21232877 16.58630137 2.078082192 0 0 0 0]
# The energy hitting location during each timestep of the average day [Wh/m^2]

###### SETS #########

nTimeSteps = length(EnergyUsed);
nTypes = 4;
TIME = 1:nTimeSteps;
SolarCellType = 1:nTypes;

###### PARAMETERS and DATA ###############

# The size of one time step [hours]
TimeStepSize = 24/nTimeSteps;
# The first time step when peak energy price is charged
PeakStartTime = 13/TimeStepSize + 1;
# The last time step when peak energy price is charged
PeakEndTime = 21/TimeStepSize;
DaysInYear = 365;


#Roof Parameters
RoofArea = 140 #m^2

# Energy market values [$/kWh]
EnergyPricePeak = 0.24769;
EnergyPriceOffPeak = 0.10558;
EnergyGenCredit = 0.21668;

#Technology-Specific Data
#[MonoCrystallline, PolyCrystalline, ThinFilm, MultiJunction]
#SolarCellEfficiency = [0.2, 0.15, 0.11, 0.41]; #Fraction

SolarCellArea = [2, 2, 2, 2]; #m^2
SolarCellOutput = [300, 300, 300, 300]; #Rated Watt Output of one panel under STC

numberpanels = zeros(50,50)
objective = zeros(50,50)
selected = zeros(50,50)

delta_ii = (20-1)/49
delta_jj = (0.45-0.05)/49

YearsInLifetime = 25;

#for ii = 1: 50
    #SelectCost = 1 + (ii-1)*delta_ii;
    #for jj = 1:50
#    SelectEfficiency = 0.05 + (jj-1)*delta_jj
    SelectCost = 5
    SelectEfficiency = 0.2
        m = Model(solver=AmplNLSolver(joinpath(PATH_TO_SOLVERS,"knitro"), ["outlev=1"]))

        @variable(m, NumberPanels >= 0)
        @variable(m, EnergyGenerated[t=TIME]) #[kWh]
        @variable(m, EnergyPurchasedIfImporting[t=TIME]) #[kWh]
        @variable(m, EnergySoldIfExporting[t=TIME]) #[kWh]
        @variable(m, Exporting[t=TIME], Bin)
        @variable(m, CostOfImported[t=TIME]) #$
        @variable(m, CreditOfExported[t=TIME]) #$
        @variable(m, CostOfImportedIfNoSolar[t=TIME]) #$
        @variable(m, SolarEnergyOut[t=TIME]) #kWh/m^2
        @variable(m, SolarEnergyIn[t=TIME])
        @variable(m, SelectType[s=SolarCellType], Bin)
        #@variable(m, SelectEfficiency >=0)
        @variable(m, SelectArea >= 0)
        #@variable(m, SelectCost >= 0)
        @variable(m, SelectOutput >= 0)

        #@constraint(m, SelectEfficiency == sum(SelectType[s]*SolarCellEfficiency[s] for s=SolarCellType))
        @constraint(m, SelectArea == sum(SelectType[s]*SolarCellArea[s] for s=SolarCellType))
        #@constraint(m, SelectCost == sum(SelectType[s]*SolarCellCost[s] for s=SolarCellType))
        @constraint(m, SelectOutput == sum(SelectType[s]*SolarCellOutput[s] for s=SolarCellType))
        # The energy hitting location during each timestep of the average day [kWh/m^2]
        @constraint(m, [t=TIME], SolarEnergyIn[t] == GlobalHorizontalIrradiance[t]*TimeStepSize/1000)
        #SolarEnergyIn = GlobalHorizontalIrradiance*TimeStepSize*1000;
        @constraint(m, [t=TIME], SolarEnergyOut[t] == SolarEnergyIn[t]*SelectEfficiency);

        @NLconstraint(m, [t=TIME], EnergyGenerated[t] == SolarEnergyOut[t]*NumberPanels*SelectArea) #kWh
        @constraint(m, [t=TIME], EnergyGenerated[t] >= EnergyUsed[t] - 100000 * (1 - Exporting[t]))
        @constraint(m, [t=TIME], EnergyUsed[t] >= EnergyGenerated[t] - 100000 * Exporting[t])
        @constraint(m, [t=TIME], EnergyPurchasedIfImporting[t] == EnergyUsed[t] - EnergyGenerated[t])
        @constraint(m, [t=TIME], EnergySoldIfExporting[t] == EnergyGenerated[t] - EnergyUsed[t])

        @constraint(m, NumberPanels*SelectArea <= RoofArea)
        @NLconstraint(m, [t=TIME; t<=PeakEndTime && t>=PeakStartTime], CostOfImportedIfNoSolar[t] == EnergyUsed[t] * EnergyPricePeak)
        @NLconstraint(m, [t=TIME; t>PeakEndTime || t<PeakStartTime], CostOfImportedIfNoSolar[t] == EnergyUsed[t] * EnergyPriceOffPeak)
        @NLconstraint(m, [t=TIME; t<=PeakEndTime && t>=PeakStartTime], CostOfImported[t] == EnergyPurchasedIfImporting[t] * (1-Exporting[t]) * EnergyPricePeak)
        @NLconstraint(m, [t=TIME; t>PeakEndTime || t<PeakStartTime], CostOfImported[t] == EnergyPurchasedIfImporting[t] * (1-Exporting[t]) * EnergyPriceOffPeak)
        @NLconstraint(m, [t=TIME], CreditOfExported[t] == EnergySoldIfExporting[t] * Exporting[t] * EnergyGenCredit)
        @constraint(m, sum(SelectType[s] for s=SolarCellType) == 1)

        @NLobjective(m, Max, (sum(CostOfImportedIfNoSolar[t] for t=TIME) - sum(CostOfImported[t] - CreditOfExported[t] for t=TIME))*DaysInYear*YearsInLifetime - SelectCost*SelectOutput*NumberPanels)

        solve(m)

        type = getvalue(SelectType)
        panels = getvalue(NumberPanels)
        generated = getvalue(EnergyGenerated)
        purchased = getvalue(EnergyPurchasedIfImporting)
        sold = getvalue(EnergySoldIfExporting)
        exporting = getvalue(Exporting)

        #selected[ii,jj] = type
#        numberpanels[ii,jj] = panels
#        objective[ii,jj] = getobjectivevalue(m)
getobjectivevalue(m)

#    end
#end
