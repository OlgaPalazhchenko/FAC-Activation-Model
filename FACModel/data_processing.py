'''
Created on Oct 22, 2017
ksdjksljgd
@author: opalazhc
'''

dgdsgdsg
import pht_model
import lepreau_data as ld
import sg_heattransfer as SGHX
import numpy as np
import csv
import matplotlib.pyplot as plt
import composition as c
from matplotlib import rc

rc('mathtext', default='regular')

RealTimeHeatTransfer = "no"
ConstantRate = "yes"
Activation = "no"
PlotOutput = "yes"
OutputLogging = "yes"
FullLoop = "no"

Solubility = []
IronConcentration = []
Loading_time = []
TotalInnerLoading = []
TotalOuterLoading = []
RIHT = [] # monitored with time 
OutletTemperatures1 = []
OutletTemperatures2 = []
# StreamOutletTemperatures = [] # monitored with time 
TemperatureProfile = []

SimulationYears = 1  # years
SimulationHours = SimulationYears * 8760

# load initial chemistry for full/half loop
pht_model.initial_chemistry(FullLoop)
# 1.52 m is the u-bend arc length of an average SG tube
Default_Tube = SGHX.closest(1.52 * 100)

import time
start_time = time.time()

for j in range(SimulationHours):
    In = pht_model.PHT_FAC(ld.InletFeeder, ld.FuelChannel, RealTimeHeatTransfer, Activation, ConstantRate, j)
    Co = pht_model.PHT_FAC(ld.FuelChannel, ld.OutletFeeder, RealTimeHeatTransfer, Activation, ConstantRate, j)
    Ou = pht_model.PHT_FAC(
        ld.OutletFeeder, ld.SteamGenerator[Default_Tube], RealTimeHeatTransfer, Activation, ConstantRate, j
        )
    
    if FullLoop == "yes":
        InletInput = ld.InletFeeder_2
    else:
        InletInput = ld.InletFeeder
    
    if RealTimeHeatTransfer == "no":
        Sg = pht_model.PHT_FAC(
            ld.SteamGenerator[Default_Tube], InletInput, RealTimeHeatTransfer, Activation, ConstantRate, j
            )
            
    else:
        # Set input concentrations for all SG zones to be same as input of first (Outlet output)
        BulkOutletOutput = [
                Ou.Section1.Bulk.FeTotal, Ou.Section1.Bulk.NiTotal, Ou.Section1.Bulk.CoTotal, Ou.Section1.Bulk.CrTotal
                ]
        SteamGeneratorTubes = []

        for Zone in SGHX.selected_tubes:
            BulkSGInput = [Zone.Bulk.FeTotal, Zone.Bulk.NiTotal, Zone.Bulk.CoTotal, Zone.Bulk.CrTotal]
            for x, y in zip(BulkSGInput, BulkOutletOutput):
                x[0] = y[Ou.Section1.NodeNumber - 1]
            
            Sg_tube = pht_model.PHT_FAC(Zone, ld.InletFeeder, RealTimeHeatTransfer, Activation, ConstantRate, j)   
            SteamGeneratorTubes.append(Sg_tube)
            
    if FullLoop == "yes":
<<<<<<< HEAD
        In_2 = pht_model.PHT_FAC(ld.InletFeeder_2, ld.FuelChannel_2, RealTimeHeatTransfer, Activation, ConstantRate, j)
        Co_2 = pht_model.PHT_FAC(ld.FuelChannel_2, ld.OutletFeeder_2, RealTimeHeatTransfer, Activation, ConstantRate, j)
=======
        In_2 = pht_model.PHT_FAC(ld.InletFeeder_2, ld.FuelChannel_2, RealTimeHeatTransfer, Activation, j)
        Co_2 = pht_model.PHT_FAC(ld.FuelChannel_2, ld.OutletFeeder_2, RealTimeHeatTransfer, Activation, j)
>>>>>>> branch 'SGHeatTransfer' of https://github.com/OlgaPalazhchenko/FAC-Activation-Model.git
        Ou_2 = pht_model.PHT_FAC(
<<<<<<< HEAD
            ld.OutletFeeder_2, ld.SteamGenerator_2[SGHX.tube_number[0]], RealTimeHeatTransfer, Activation, ConstantRate,
            j
=======
            ld.OutletFeeder_2, ld.SteamGenerator_2[Default_Tube], RealTimeHeatTransfer, Activation, j
>>>>>>> branch 'SGHeatTransfer' of https://github.com/OlgaPalazhchenko/FAC-Activation-Model.git
            )
        Sg_2 = pht_model.PHT_FAC(
<<<<<<< HEAD
            ld.SteamGenerator_2[SGHX.tube_number[0]], ld.InletFeeder, RealTimeHeatTransfer, Activation, ConstantRate, j
=======
            ld.SteamGenerator_2[Default_Tube], ld.InletFeeder, RealTimeHeatTransfer, Activation, j
>>>>>>> branch 'SGHeatTransfer' of https://github.com/OlgaPalazhchenko/FAC-Activation-Model.git
            )
    
    
    if j % 8759 == 0:  # yearly  
        # parameters tracked with time 
        T_RIH = (SGHX.energy_balance(
            ld.SteamGenerator[Default_Tube].NodeNumber - 1, ld.SteamGenerator[Default_Tube].InnerOxThickness,
            ld.SteamGenerator[Default_Tube].OuterOxThickness, j
            ) - 273.15)
        
        for Zone in ld.SteamGenerator:
                Zone.Bulk.FeSatFe3O4 = c.iron_solubility(Zone)
            
        if OutputLogging == "yes":
            RIHT.append(T_RIH)
            
            Temperature1 = (
                ld.SteamGenerator[SGHX.tube_number[0]].PrimaryBulkTemperature[21] - 273.15
                           )
            OutletTemperatures1.append(Temperature1)
            
            if len(SGHX.selected_tubes) > 1:
                Temperature2 = (
                    ld.SteamGenerator[SGHX.tube_number[0]].PrimaryBulkTemperature[21] - 273.15
                               )
                OutletTemperatures2.append(Temperature2)
            
            
            Loading_time.append(ld.SteamGenerator[SGHX.tube_number[0]].SolutionOxide.FeSatFe3O4)
    else:
        None

          
# FACRate = sum(
#     ld.UnitConverter(Ou.Section1, "Corrosion Rate Grams", "Corrosion Rate Micrometers", None, Ou.Section1.CorrRate,
#     None, None, None, None)
#     ) / (Ou.Section1.NodeNumber)


# parameters at the end of run 
for Zone in SGHX.selected_tubes:
    x = ld.UnitConverter(
    Zone, "Grams per Cm Squared", "Grams per M Squared", None, None, Zone.InnerIronOxThickness, None, None, None
    )
    y = ld.UnitConverter(
    Zone, "Grams per Cm Squared", "Grams per M Squared", None, None, Zone.OuterFe3O4Thickness, None, None, None
    )
    
#     totalloading = [i + j for i, j in zip(x, y)]
#     z = sum(totalloading[11:len(totalloading)]) / (len(totalloading) - 11)
    
    TotalInnerLoading.append(x)
    TotalOuterLoading.append(y)
    Solubility.append(Zone.SolutionOxide.FeSatFe3O4)
    IronConcentration.append(Zone.SolutionOxide.FeTotal)
    Temperature_C = [i - 273.15 for i in Zone.PrimaryBulkTemperature]
    TemperatureProfile.append(Temperature_C)
    
Data = [SGHX.ubends, TotalInnerLoading, TotalOuterLoading, Solubility, IronConcentration, TemperatureProfile]
Labels = [
    "U-bend length (m)", "Inner Loading (g/m^2)", "Outer Loading (g/m^2)", "Solubility (mol/kg)", "S/O [Fe] (mol/kg)",
    "Temperature Profile (oC)"]
    
csvfile = "RIHTOutput.csv"
with open(csvfile, "w") as output:
    writer = csv.writer(output, lineterminator='\n')
    writer.writerow(['RIHT (oC)'])
    writer.writerow(RIHT)
    writer.writerow([''])
     
    for i, j in zip(Labels, Data):
        writer.writerow([i])
        if i == "U-bend length (m)":
            writer.writerow(j)
        else:
            writer.writerows(j)
        writer.writerow([''])
        
    writer.writerow(['Outlet Streams (oC)'])
    writer.writerow(OutletTemperatures1)
    writer.writerow(OutletTemperatures2)
    
    writer.writerow([''])
    writer.writerow(['FeSat over time (mol/kg)'])
    writer.writerows(Loading_time)
    
end_time = time.time()
delta_time = end_time - start_time
 
hours = delta_time // 3600
temp = delta_time - 3600 * hours
minutes = delta_time // 60
seconds = delta_time - 60 * minutes
print('%d:%d:%d' % (hours, minutes, seconds))


def property_log10(Element, Interface):
    Sat = []  # x
    Bulk = []  # y
    SolutionOxide = []  # z
    
    if RealTimeHeatTransfer == "yes":
        SteamGenerator = SteamGeneratorTubes[0].Section1
    else:
        SteamGenerator = Sg.Section1
    
    # only in main 4 PHTS sections, not counting SG Zones, can be changed to include all, if needed 
    for Section in [In.Section1, Co.Section1, Ou.Section1, SteamGenerator]:
        if Element == "Fe":
            Concentrations = [Section.SolutionOxide.FeSatFe3O4, Section.Bulk.FeTotal, Section.SolutionOxide.FeTotal]
        elif Element == "Ni":
            Concentrations = [Section.SolutionOxide.NiSatFerrite, Section.Bulk.NiTotal, Section.SolutionOxide.NiTotal]
        elif Element == "Cr":
            Concentrations = [Section.SolutionOxide.CrSat, Section.Bulk.CrTotal, Section.SolutionOxide.CrTotal]
        elif Element == "Co":
            Concentrations = [
                Section.SolutionOxide.CoSatFerrite, Section.Bulk.CoTotal, Section.SolutionOxide.CoTotal
                ]
        else:
            None
            
        x = np.log10(Concentrations[0])
        y = np.log10(Concentrations[1])
        z = np.log10(Concentrations[2])
    
        Sat.append(x)
        Bulk.append(y)
        SolutionOxide.append(z)
        
    Sat = [j for i in Sat for j in i]
    Bulk = [j for i in Bulk for j in i]
    SolutionOxide = [j for i in SolutionOxide for j in i]
     
    if Interface == "Sat":
        return Sat
    elif Interface == "Bulk":
        return Bulk
    elif Interface == "SO":
        return SolutionOxide
    else:
        return None


def oxide_loading(Layer, RealTimeHeatTransfer):
    Oxide = []
    if RealTimeHeatTransfer == "yes":
        SteamGenerator = SteamGeneratorTubes[0].Section1
    else:
        SteamGenerator = Sg.Section1
        
    for Section in [In.Section1, Co.Section1, Ou.Section1, SteamGenerator]:
        if Layer == "Inner":
            x = Section.InnerOxThickness
        elif Layer == "Outer":
            x = Section.OuterOxThickness
        elif Layer == "Cobalt":
            x = Section.CoThickness
        elif Layer == "Nickel":
            x = Section.NiThickness
        else:
            None 
        Oxide.append(x)
    Oxide = [j for i in Oxide for j in i]
    # unpacks 4 separate lists into one list of oxide thickness down entire PHT
    ConvertedOxide = ld.UnitConverter(
        Section, "Grams per Cm Squared", "Grams per M Squared", None, None, Oxide, None, None, None
        )
    return ConvertedOxide        


def plot_output():
    LoopDistance = []
    for Section in [ld.InletFeeder, ld.FuelChannel, ld.OutletFeeder, ld.SteamGenerator[Default_Tube]]:
        x = Section.Length.magnitude
        LoopDistance.append(x)
    LoopDistance = [j for i in LoopDistance for j in i]
    TotalLoopDistance = [i / 100 for i in np.cumsum(LoopDistance)]  # Distance down length of PHTS [m]
    
    fig1 = plt.figure()
    
    ax1 = fig1.add_subplot(221)
    ax1.plot(TotalLoopDistance, property_log10("Fe", "SO"), marker='o', color='r', label='S/O Interface')
    ax1.plot(TotalLoopDistance, property_log10("Fe", "Sat"), marker='o', color='0', label='Magnetite Solubility')
    ax1.plot(TotalLoopDistance, property_log10("Fe", "Bulk"), marker='o', color='0.50', label='Bulk Coolant')
    ax1.set_xlabel('Distance (m)')
    ax1.set_ylabel('$Log_{10}$[Fe Concentration (mol/kg)]')
    # ax1.legend()
  
    ax2 = fig1.add_subplot(222)
    ax2.plot(TotalLoopDistance, property_log10("Ni", "SO"), marker='o', color='b', label='S/O Interface')
    ax2.plot(TotalLoopDistance, property_log10("Ni", "Sat"), marker='o', color='0', label='Ni-Ferrite Solubility')
    ax2.plot(TotalLoopDistance, property_log10("Ni", "Bulk"), marker='o', color='0.50', label='Bulk Coolant') 
    ax2.set_xlabel('Distance (m)')
    ax2.set_ylabel('$Log_{10}$[Ni Concentration (mol/kg)]')
    # ax2.legend()
  
    ax3 = fig1.add_subplot(223)
    ax3.plot(TotalLoopDistance, property_log10("Co", "SO"), marker='^', color='m', label='S/O Interface')
    ax3.plot(TotalLoopDistance, property_log10("Co", "Sat"), marker='^', color='0', label='Co-Ferrite Solubility')
    ax3.plot(TotalLoopDistance, property_log10("Co", "Bulk"), marker='^', color='0.50', label='Bulk Coolant') 
    ax3.set_xlabel('Distance (m)')
    ax3.set_ylabel('$Log_{10}$[Co Concentration (mol/kg)]')
    # ax3.legend()
  
  
    ax4 = fig1.add_subplot(224)
    ax4.plot(TotalLoopDistance, property_log10("Cr", "Sat"), marker='*', color='0', label='Chromite Solubility')
    ax4.plot(TotalLoopDistance, property_log10("Cr", "Bulk"), marker='*', color='0.50', label='Bulk Coolant') 
    ax4.set_xlabel('Distance (m)')
    ax4.set_ylabel('$Log_{10}$[Cr Concentration (mol/kg)]')
      
    plt.tight_layout()
    plt.show()
    
    fig2, ax1 = plt.subplots()
    # ax1 = fig2.add_subplot(221)
    ax1.plot(
        TotalLoopDistance, oxide_loading("Inner", RealTimeHeatTransfer), linestyle=None, marker='o', color='0.50',
        label='Inner Oxide'
        )
    ax1.plot(
        TotalLoopDistance, oxide_loading("Outer", RealTimeHeatTransfer), linestyle=None, marker='o', color='k',
        label='Outer Oxide'
        )
    # ax1.axis([51,69, 0, 30])
    ax1.set_xlabel('Distance (m)')
    ax1.set_ylabel('Oxide Layer Loadings (${g/m^2}$)')
                  
    ax2 = ax1.twinx()
    ax2.plot(
        TotalLoopDistance, oxide_loading("Nickel", RealTimeHeatTransfer), linestyle=None, marker='o', color='c',
        label='Nickel'
        )
    ax2.plot(
        TotalLoopDistance, oxide_loading("Cobalt", RealTimeHeatTransfer), linestyle=None, marker='o', color='m',
        label='Cobalt'
        )
    ax2.set_ylabel('Ni, Co, and Cr Loadings (${g/m^2}$)', rotation=270, labelpad=20)
    lines, labels = ax1.get_legend_handles_labels()
    lines2, labels2 = ax2.get_legend_handles_labels()
    ax2.legend(lines + lines2, labels + labels2, loc=0)        
    plt.axis([4, 69, 0, 0.01])
    plt.tight_layout()
    plt.show()
    
if PlotOutput == "yes":
    plot_output()

