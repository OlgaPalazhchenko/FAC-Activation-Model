'''
Created on Oct 22, 2017

@author: opalazhc
'''
import pht_model
import lepreau_data as ld
import sg_heattransfer as SGHX
import numpy as np
import csv
import matplotlib.pyplot as plt
import composition as c
from matplotlib import rc

rc('mathtext', default='regular')

RealTimeHeatTransfer = "yes"
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

SimulationYears = 9  # years
SimulationHours = SimulationYears * 8760

import time
start_time = time.time()

for j in range(SimulationHours):
    In = pht_model.PHT_FAC(ld.Inlet, ld.Core, RealTimeHeatTransfer, Activation, j)
    Co = pht_model.PHT_FAC(ld.Core, ld.Outlet, RealTimeHeatTransfer, Activation, j)
    Ou = pht_model.PHT_FAC(ld.Outlet, ld.SGZones[SGHX.tube_number[0]], RealTimeHeatTransfer, Activation, j)
    
    if FullLoop == "yes":
        InletInput = ld.Inlet_2
    else:
        InletInput = ld.Inlet
    
    if RealTimeHeatTransfer == "no":
        Sg = pht_model.PHT_FAC(ld.SGZones[SGHX.tube_number[0]], InletInput, RealTimeHeatTransfer, Activation, j)
            
    else:
        # Set input concentrations for all SG zones to be same as input of first (Outlet output)
        BulkOutletOutput = [
                Ou.Section1.Bulk.FeTotal, Ou.Section1.Bulk.NiTotal, Ou.Section1.Bulk.CoTotal, Ou.Section1.Bulk.CrTotal
                ]
        SteamGeneratorTubes = []

        for Zone in SGHX.desired_tubes:
            BulkSGInput = [Zone.Bulk.FeTotal, Zone.Bulk.NiTotal, Zone.Bulk.CoTotal, Zone.Bulk.CrTotal]
            for x, y in zip(BulkSGInput, BulkOutletOutput):
                x[0] = y[Ou.Section1.NodeNumber - 1]
            
            Sg_tube = pht_model.PHT_FAC(Zone, ld.Inlet, RealTimeHeatTransfer, Activation, j)   
            SteamGeneratorTubes.append(Sg_tube)
            
    if FullLoop == "yes":
        In_2 = pht_model.PHT_FAC(ld.Inlet_2, ld.Core, RealTimeHeatTransfer, Activation, j)
    
    if OutputLogging == "yes":
        if j % 3000 == 0:  # yearly
            
            # parameters tracked with time 
            T_RIH = (SGHX.energy_balance(
                ld.SGZones[SGHX.tube_number[0]].NodeNumber - 1, ld.SGZones[SGHX.tube_number[0]].InnerOxThickness,
                ld.SGZones[SGHX.tube_number[0]].OuterOxThickness, j
                ) - 273.15)
            RIHT.append(T_RIH)
            
            for Zone in ld.SGZones:
                Zone.Bulk.FeSatFe3O4 = c.iron_solubility(Zone)
            
            Temperature1 = (
                ld.SGZones[SGHX.tube_number[0]].PrimaryBulkTemperature[ld.SGZones[SGHX.tube_number[0]].NodeNumber - 1]
                           - 273.15
                           )
            Temperature2 = (
                ld.SGZones[SGHX.tube_number[0]].PrimaryBulkTemperature[ld.SGZones[SGHX.tube_number[0]].NodeNumber - 1]
                           - 273.15
                           )
            OutletTemperatures1.append(Temperature1)
            OutletTemperatures2.append(Temperature2)
            
            Loading_time.append(ld.SGZones[SGHX.tube_number[0]].SolutionOxide.FeSatFe3O4)
    else:
        None
            
# FACRate = sum(
#     ld.UnitConverter(Ou.Section1, "Corrosion Rate Grams", "Corrosion Rate Micrometers", None, Ou.Section1.CorrRate,
#     None, None, None, None)
#     ) / (Ou.Section1.NodeNumber)


# parameters at the end of run 
for Zone in SGHX.desired_tubes:
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
    
Data = [SGHX.desired_ubends, TotalInnerLoading, TotalOuterLoading, Solubility, IronConcentration, TemperatureProfile]
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
    
    # only in main 4 PHTS sections, not counting SG Zones, can be changed to include all, if needed 
    for Section in [In.Section1, Co.Section1, Ou.Section1, ld.SGZones[SGHX.tube_number[0]]]:
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
    for Section in [ld.Inlet, ld.Core, ld.Outlet, ld.SGZones[12]]:
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

