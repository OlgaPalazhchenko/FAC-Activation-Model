'''
Created on Oct 22, 2017
@author: opalazhc
'''
import pht_model
import lepreau_data as ld
import sg_heattransfer as SGHX
import composition as c
import numpy as np
import thermochemistry_and_constants as nc
import csv

import matplotlib.pyplot as plt
from matplotlib import rc
rc('mathtext', default='regular')

# preset corrosion rate instead of calculated from FAC model
ConstantRate = "yes"
Activation = "no"
PlotOutput = "yes"
OutputLogging = "yes"
Loop = "half"
ElementTracking = "no"

# 1.52 m is the u-bend arc length of an average SG tube
Default_Tube = SGHX.closest_ubend(1.52 * 100)

          
SimulationYears = 15  # years
SimulationHours = SimulationYears * 876


if OutputLogging == "yes":
    Solubility = []
    IronConcentration = []
    Loading_time = []
    TotalInnerLoading = []
    TotalOuterLoading = []
    TotalOxide = []
    RIHT = [] # monitored with time 
    OutletTemperatures1 = [] 
    OutletTemperatures2 = []
    pht_SteamFraction = []
    # StreamOutletTemperatures = [] # monitored with time 
    TemperatureProfile = []
    TotalDistance = []
    Years = []
    for i in range((SimulationYears + 1) * 4):
        Years.append((i / 4) + SGHX.YearStartup)


# load initial chemistry for full/half loop
pht_model.initial_chemistry(Loop)


import time
start_time = time.time()


def sg_heat_transfer(Outlet, InletInput, j):
    Tubes = []
    # Set input concentrations for all SG zones to be same as output of outlet feeder
    BulkOutletActivityOutput = [
        Outlet.Bulk.Co60, Outlet.Bulk.Co58, Outlet.Bulk.Fe59, Outlet.Bulk.Fe55, Outlet.Bulk.Cr51, Outlet.Bulk.Mn54,
        Outlet.Bulk.Ni63
        ]
    
    BulkOutletOutput = [Outlet.Bulk.FeTotal, Outlet.Bulk.NiTotal, Outlet.Bulk.CoTotal, Outlet.Bulk.CrTotal]
    
    for Tube in SGHX.selected_tubes:
        BulkSGInput = [Tube.Bulk.FeTotal, Tube.Bulk.NiTotal, Tube.Bulk.CoTotal, Tube.Bulk.CrTotal]
        BulkSGActivityInput = [
        Tube.Bulk.Co60, Tube.Bulk.Co58, Tube.Bulk.Fe59, Tube.Bulk.Fe55, Tube.Bulk.Cr51, Tube.Bulk.Mn54, Tube.Bulk.Ni63
        ]
        for x, y, z, q in zip(BulkSGInput, BulkOutletOutput, BulkSGActivityInput, BulkOutletActivityOutput):
            x[0] = y[Outlet.NodeNumber - 1]
            z[0] = q[Outlet.NodeNumber - 1]
        
        w = pht_model.PHT_FAC(Tube, InletInput, ElementTracking, Activation, ConstantRate, j)   
        Tubes.append(w)
    return Tubes


for j in range(SimulationHours):
    In = pht_model.PHT_FAC(ld.InletFeeder, ld.FuelChannel, ElementTracking, Activation, ConstantRate, j)
    Co = pht_model.PHT_FAC(ld.FuelChannel, ld.OutletFeeder, ElementTracking, Activation, ConstantRate, j)
    Ou = pht_model.PHT_FAC(
        ld.OutletFeeder, ld.SteamGenerator[Default_Tube], ElementTracking, Activation, ConstantRate, j
        ) 
    if Loop == "full":
        InletInput = ld.InletFeeder_2
    else:
        InletInput = ld.InletFeeder
    
    Sg = pht_model.PHT_FAC(ld.SteamGenerator[Default_Tube], InletInput, ElementTracking, Activation, ConstantRate, j)
    SteamGeneratorTubes = sg_heat_transfer(Ou.Section1, InletInput, j)

#     print (
#         ld.UnitConverter(Ou.Section1, "Corrosion Rate Grams", "Corrosion Rate Micrometers", None, Ou.Section1.CorrRate,
#     None, None, None, None)
#            )
    if Loop == "full":
        In_2 = pht_model.PHT_FAC(ld.InletFeeder_2, ld.FuelChannel_2, ElementTracking, Activation, ConstantRate, j)
        Co_2 = pht_model.PHT_FAC(ld.FuelChannel_2, ld.OutletFeeder_2, ElementTracking, Activation, ConstantRate, j)
        Ou_2 = pht_model.PHT_FAC(
            ld.OutletFeeder_2, ld.SteamGenerator_2[SGHX.tube_number[0]], ElementTracking, Activation, ConstantRate,
            j
            )
        Sg_2 = pht_model.PHT_FAC(
            ld.SteamGenerator_2[SGHX.tube_number[0]], ld.InletFeeder, ElementTracking, Activation, ConstantRate, j
            )   
    
    if j % (2190 / nc.TIME_STEP) == 0:  # twice a year  
        
        if j ==0:
            x_pht = 0.002 # PHT steam fraction for "clean" boiler

        
        #pass cleaned and uncleaned tubes into heat transfer function
#         CleanedInnerOxide = SteamGeneratorTubes[0].Section1.InnerOxThickness
#         CleanedOuterOxide = SteamGeneratorTubes[0].Section1.OuterOxThickness
#         
#         UncleanedInnerOxide = Sg.Section1.InnerOxThickness
#         UncleanedOuterOxide = Sg.Section1.OuterOxThickness
        
        # parameters tracked with time 
        T_RIH = (SGHX.energy_balance(ld.SteamGenerator[Default_Tube].NodeNumber - 1, x_pht, j) - 273.15)
        
#         T_RIH = (
#             SGHX.energy_balance(ld.SteamGenerator[Default_Tube].NodeNumber - 1, UncleanedInnerOxide,
#                                 UncleanedOuterOxide, CleanedInnerOxide, CleanedOuterOxide, x_pht, j) - 273.15
#                  )

        print (1983 + j / (8760 / nc.TIME_STEP), x_pht, T_RIH)
        x_pht = SGHX.pht_steam_quality(T_RIH + 273.15, j)
        
        InletInput.PrimaryBulkTemperature = [T_RIH + 273.15] * InletInput.NodeNumber
        for Section in ld.HalfLoop:
            Section.Bulk.FeSatFe3O4 = c.iron_solubility(Section, None)
            
        if OutputLogging == "yes":
            RIHT.append(T_RIH)
            pht_SteamFraction.append(x_pht)
            Temperature1 = (
                ld.SteamGenerator[SGHX.tube_number[0]].PrimaryBulkTemperature[21] - 273.15
                           )
            OutletTemperatures1.append(Temperature1)
            
            # if multiple tubes arc/total lengths are run
            if len(SGHX.selected_tubes) > 1:
                Temperature2 = (
                    ld.SteamGenerator[SGHX.tube_number[1]].PrimaryBulkTemperature[21] - 273.15
                               )
                OutletTemperatures2.append(Temperature2)
            
    else:
        None
          
for Zone in [SGHX.selected_tubes[0], ld.SteamGenerator[Default_Tube]]:
    x = ld.UnitConverter(
    Zone, "Grams per Cm Squared", "Grams per M Squared", None, None, Zone.InnerIronOxThickness, None, None, None
    )
    y = ld.UnitConverter(
    Zone, "Grams per Cm Squared", "Grams per M Squared", None, None, Zone.OuterFe3O4Thickness, None, None, None
    )
    z = [i / 100 for i in Zone.Distance]
    q = [i + j for i, j in zip(x, y)]
#     z = sum(totalloading[11:len(totalloading)]) / (len(totalloading) - 11)
    d = ld.UnitConverter(Zone, "Mol per Kg", "Grams per Cm Cubed", Zone.SolutionOxide.FeSatFe3O4, None, None, None,
                         nc.FeMolarMass, None)
    e = ld.UnitConverter(Zone, "Mol per Kg", "Grams per Cm Cubed", Zone.SolutionOxide.FeTotal, None, None, None,
                         nc.FeMolarMass, None)
    
    TotalDistance.append(z)
    TotalInnerLoading.append(x)
    TotalOuterLoading.append(y)
    TotalOxide.append(q)
    Solubility.append(d) # Zone.SolutionOxide.FeSatFe3O4
    IronConcentration.append(e) # Zone.SolutionOxide.FeTotal
    Temperature_C = [i - 273.15 for i in Zone.PrimaryBulkTemperature]
    TemperatureProfile.append(Temperature_C)
    
Data = [SGHX.TubeLengths, TotalDistance, TotalInnerLoading, TotalOuterLoading, TotalOxide, Solubility, IronConcentration,
        TemperatureProfile]
Labels = [
    "U-bend length (cm)", "Distance (m)", "Inner Loading (g/m^2)", "Outer Loading (g/m^2)", "Total Oxide (g/m^2)",
    "Solubility (mol/kg)", "S/O [Fe] (mol/kg)", "Temperature Profile (oC)"]

RIHT_delta = [x-y for x, y in zip (RIHT[1:], RIHT)]   

csvfile = "RIHTOutput.csv"
with open(csvfile, "w") as output:
    writer = csv.writer(output, lineterminator='\n')
    writer.writerow(['RIHT (oC) and year'])
    writer.writerow(RIHT)
    writer.writerow(Years)
    writer.writerow(['Delta RIHT (oC'])
    writer.writerow(RIHT_delta)
    writer.writerow(['Steam fraction'])
    writer.writerow(pht_SteamFraction)
    writer.writerow([''])
    
    for i, j in zip(Labels, Data):
        writer.writerow([i])
        if i == "U-bend length (cm)":
            writer.writerow(j)
        else:
            writer.writerows(j)
        writer.writerow([''])
        
    writer.writerow(['Outlet Streams (oC)'])
    writer.writerow(OutletTemperatures1)
    writer.writerow(OutletTemperatures2)
    
    
end_time = time.time()
delta_time = end_time - start_time
 
hours = delta_time // 3600
temp = delta_time - 3600 * hours
minutes = delta_time // 60
seconds = delta_time - 60 * minutes
print('%d:%d:%d' % (hours, minutes, seconds))



if Loop == "full" or Loop == "half":
    Sections = [In.Section1, Co.Section1, Ou.Section1, SteamGeneratorTubes[0].Section1]
else:
    Sections = [In.Section1, Co.Section1, Ou.Section1, Sg.Section1]
    
    
def property_log10(Element, Interface):
    Sat = []  # x
    Bulk = []  # y
    SolutionOxide = []  # z
    
    # only in main 4 PHTS sections, not counting SG Zones, can be changed to include all, if needed 
    for Section in Sections:
        if Element == "Fe":
            Concentrations = [Section.SolutionOxide.FeSatFe3O4, Section.Bulk.FeTotal, Section.SolutionOxide.FeTotal]
        elif Element == "Ni":
            Concentrations = [Section.Bulk.NiSatFerrite, Section.Bulk.NiTotal, Section.SolutionOxide.NiTotal]
        elif Element == "Cr":
            Concentrations = [Section.Bulk.CrSat, Section.Bulk.CrTotal, Section.SolutionOxide.CrTotal]
        elif Element == "Co":
            Concentrations = [
                Section.Bulk.CoSatFerrite, Section.Bulk.CoTotal, Section.SolutionOxide.CoTotal
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


def oxide_loading(Layer):
    Oxide = []
    
    for Section in Sections:
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


def activity_volumetric(Isotope):
    VolumetricActivity = []
    
    for Section in Sections:
        if Isotope == "Co60":
            x = Section.Bulk.Co60
        elif Isotope == "Co58":
            x = Section.Bulk.Co58
        elif Isotope == "Fe55":
            x = Section.Bulk.Fe55
        elif Isotope == "Fe59":
            x = Section.Bulk.Fe59
        elif Isotope == "Cr51":
            x = Section.Bulk.Cr51
        elif Isotope == "Mn54":
            x = Section.Bulk.Mn54
        elif Isotope == "Ni63":
            x = Section.Bulk.Ni63
        else:
            None
        VolumetricActivity.append(x)
    
    VolumetricActivity = [j for i in VolumetricActivity for j in i]
    
    return VolumetricActivity

   
def plot_output():
    LoopDistance = []
    if Loop == "full" or Loop == "half":
        Sections = [ld.InletFeeder, ld.FuelChannel, ld.OutletFeeder, SteamGeneratorTubes[0].Section1]
    else:
        Sections = [ld.InletFeeder, ld.FuelChannel, ld.OutletFeeder, ld.SteamGenerator[Default_Tube]]
    
    for Section in Sections:
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
    ax1.plot(
        TotalLoopDistance, oxide_loading("Inner"), linestyle=None, marker='o', color='0.50',
        label='Inner Oxide'
        )
    ax1.plot(
        TotalLoopDistance, oxide_loading("Outer"), linestyle=None, marker='o', color='k',
        label='Outer Oxide'
        )
    # ax1.axis([51,69, 0, 30])
    ax1.set_xlabel('Distance (m)')
    ax1.set_ylabel('Oxide Layer Loadings (${g/m^2}$)')
                  
    ax2 = ax1.twinx()
    ax2.plot(
        TotalLoopDistance, oxide_loading("Nickel"), linestyle=None, marker='o', color='c',
        label='Nickel'
        )
    ax2.plot(
        TotalLoopDistance, oxide_loading("Cobalt"), linestyle=None, marker='o', color='m',
        label='Cobalt'
        )
    ax2.set_ylabel('Ni, Co, and Cr Loadings (${g/m^2}$)', rotation=270, labelpad=20)
    lines, labels = ax1.get_legend_handles_labels()
    lines2, labels2 = ax2.get_legend_handles_labels()
    ax2.legend(lines + lines2, labels + labels2, loc=0)        
    plt.axis([4, 71, 0, 10])
    plt.tight_layout()
    plt.show()
    
    
    fig3, ax1 = plt.subplots()
    ax1.plot(
        TotalLoopDistance, activity_volumetric("Co60"), linestyle=None, marker='o', color='0.50',
        label='Co60'
        )
    ax1.plot(
        TotalLoopDistance, activity_volumetric("Fe55"), linestyle=None, marker='o', color='m',
        label='Co58'
        )
    ax1.plot(
        TotalLoopDistance, activity_volumetric("Cr51"), linestyle=None, marker='o', color='c',
        label='Cr51'
        )
    ax1.plot(
        TotalLoopDistance, activity_volumetric("Fe59"), linestyle=None, marker='o', color='g',
        label='Fe59'
        )
    
    ax1.set_xlabel('Distance (m)')
    ax1.set_ylabel('Coolant Activity (${uCi/m^3}$)')
    plt.tight_layout()
    plt.show()
    
if PlotOutput == "yes":
    plot_output()