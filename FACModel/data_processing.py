'''
Created on Oct 22, 2017
@author: opalazhc
'''

import lepreau_data as ld
import sg_heattransfer as SGHX
import composition as c
import numpy as np
import thermochemistry_and_constants as nc
import csv
import iteration as it
import pht_model

import matplotlib.pyplot as plt
from matplotlib import rc
rc('mathtext', default='regular')

PlotOutput = "yes"

# for j in range(SimulationHours):
#     In = pht_model.PHT_FAC(ld.InletFeeder, ld.FuelChannel, ElementTracking, Activation, ConstantRate, j)
#     Co = pht_model.PHT_FAC(ld.FuelChannel, ld.OutletFeeder, ElementTracking, Activation, ConstantRate, j)
#     Ou = pht_model.PHT_FAC(
#         ld.OutletFeeder, ld.SteamGenerator[Default_Tube], ElementTracking, Activation, ConstantRate, j
#         ) 
#     if Loop == "full":
#         InletInput = ld.InletFeeder_2
#     else:
#         InletInput = ld.InletFeeder
#     
#     Sg = pht_model.PHT_FAC(ld.SteamGenerator[Default_Tube], InletInput, ElementTracking, Activation, ConstantRate, j)
#     SteamGeneratorTubes = sg_heat_transfer(Ou.Section1, InletInput, j)
#     
#     if Loop == "full":
#         In_2 = pht_model.PHT_FAC(ld.InletFeeder_2, ld.FuelChannel_2, ElementTracking, Activation, ConstantRate, j)
#         Co_2 = pht_model.PHT_FAC(ld.FuelChannel_2, ld.OutletFeeder_2, ElementTracking, Activation, ConstantRate, j)
#         Ou_2 = pht_model.PHT_FAC(
#             ld.OutletFeeder_2, ld.SteamGenerator_2[SGHX.tube_number[0]], ElementTracking, Activation, ConstantRate,
#             j
#             )
#         Sg_2 = pht_model.PHT_FAC(
#             ld.SteamGenerator_2[SGHX.tube_number[0]], ld.InletFeeder, ElementTracking, Activation, ConstantRate, j
#             )
#         SteamGeneratorTubes_2 = sg_heat_transfer(Ou_2.Section1, ld.InletFeeder, j)   
#     
#     if j % (219) == 0:  # 2190 h * 10 = 4x a year  
#         
#         if j ==0:
#             x_pht = 0.0245 # PHT steam fraction for "clean" boiler
# 
#         if SGHX.SGFastMode == "yes":
#             #pass cleaned and uncleaned tubes into heat transfer function
#             CleanedInnerOxide = SteamGeneratorTubes[0].Section1.InnerOxThickness
#             CleanedOuterOxide = SteamGeneratorTubes[0].Section1.OuterOxThickness
#              
#             UncleanedInnerOxide = Sg.Section1.InnerOxThickness
#             UncleanedOuterOxide = Sg.Section1.OuterOxThickness
#         else:
#             CleanedInnerOxide = None
#             CleanedOuterOxide = None
# 
#             UncleanedInnerOxide = None
#             UncleanedOuterOxide = None
#             
#         # parameters tracked with time 
#         T_RIH = (SGHX.energy_balance(
#             ld.SteamGenerator[Default_Tube].NodeNumber - 1, UncleanedInnerOxide, UncleanedOuterOxide,
#             CleanedInnerOxide, CleanedOuterOxide, x_pht, j, SGHX.SGFastMode) - 273.15)
# 
# 
#         print (SGHX.YearStartup + j / (8760 / nc.TIME_STEP), x_pht, T_RIH)
#         x_pht = SGHX.pht_steam_quality(T_RIH + 273.15, j)
#         
#         InletInput.PrimaryBulkTemperature = [T_RIH + 273.15] * InletInput.NodeNumber
#         for Section in ld.HalfLoop:
#             Section.Bulk.FeSatFe3O4 = c.iron_solubility_SB(Section)
#             
#         if OutputLogging == "yes":
# #             Time.append(j)
#             # S/O concentrations are updated individually in code, e.g., self.Section1.FeTotal
#             # Bulk concentrations are not, so same location is pointed to, so all places self.Section1.Bulk.FeTotal
#             # was used change to updated values, even previously appended values
#             # new list needs to be created for bulk concentrations only
# #             x = list(In.Section1.Bulk.FeTotal)
#             
# #             SGOxide.append(Sg.Section1.OuterFe3O4Thickness)
# #             OutletOxide.append(Ou.Section1.InnerIronOxThickness)
# #             InletBulkConcentration.append(x)
# #             OutletSOConcentration.append(Ou.Section1.SolutionOxide.FeTotal)
#             OutletCorrosionRate.append(Ou.Section1.CorrRate)
#             
#             RIHT.append(T_RIH)
#             pht_SteamFraction.append(x_pht)
#             Temperature1 = (
#                 ld.SteamGenerator[SGHX.tube_number[0]].PrimaryBulkTemperature[21] - 273.15
#                            )
#             OutletTemperatures1.append(Temperature1)
#             
#             # if multiple tubes arc/total lengths are run
#             if len(SGHX.selected_tubes) > 1:
#                 Temperature2 = (
#                     ld.SteamGenerator[SGHX.tube_number[1]].PrimaryBulkTemperature[21] - 273.15
#                                )
#                 OutletTemperatures2.append(Temperature2)
#             
#     else:
#         None


def purification_csv():
    
    OutletCorrosionRate_uma = []
    InletBulkConcentration_gcm3 = []
    OutletSOConcentration_gcm3 = []
    AvgCorrRate = []
     
    for Rate, Conc1, Conc2 in zip(OutletCorrosionRate, InletBulkConcentration, OutletSOConcentration):
        x = ld.UnitConverter(
        Ou.Section1, "Corrosion Rate Grams", "Corrosion Rate Micrometers", None, Rate, None, None, None, None
        )
        y = ld.UnitConverter(
            In.Section1, "Mol per Kg", "Grams per Cm Cubed", Conc1, None, None, None, nc.FeMolarMass, None)
         
        z = ld.UnitConverter(
            Ou.Section1, "Mol per Kg", "Grams per Cm Cubed", Conc2, None, None, None, nc.FeMolarMass, None
            )
        q = sum(x) / Ou.Section1.NodeNumber
         
        OutletCorrosionRate_uma.append(x)
        InletBulkConcentration_gcm3.append(y)
        OutletSOConcentration_gcm3.append(z)
        AvgCorrRate.append(q)
     
    OutletSolubility_gcm3 = ld.UnitConverter(
            Ou.Section1, "Mol per Kg", "Grams per Cm Cubed", Ou.Section1.SolutionOxide.FeSatFe3O4, None, None, None,
            nc.FeMolarMass, None
            )
     
    DissolutionRate = []
    InnerOxideGrowth = []
    DeltaOx = []
    for Concentration, Rate in zip(OutletSOConcentration_gcm3, OutletCorrosionRate):
        diss = [nc.KdFe3O4 * (x - y) for x, y in zip (OutletSolubility_gcm3, Concentration)]
        growth = [i * it.Diffusion(Section, "Fe") / Section.FractionFeInnerOxide for i in Rate]
        delta = [x - y for x, y in zip(growth, diss)]
        DissolutionRate.append(diss)
        InnerOxideGrowth.append(growth)
        DeltaOx.append(delta)
     
    csvfile = "PurificationOutput.csv"
    with open(csvfile, "w") as output:
        writer = csv.writer(output, lineterminator='\n')
        writer.writerow(['1.67 %'])
        writer.writerow([''])
        writer.writerow(['Outlet Corrosion Rate (g/cm^2 s and um/a)'])
        writer.writerows(OutletCorrosionRate_uma)
    #     writer.writerows(OutletCorrosionRate)
         
        writer.writerow([''])
        writer.writerow(['Time (s)'])
        writer.writerow(Time)
        writer.writerow(['Average Corrosion Rate (um/a)'])
        writer.writerow(AvgCorrRate)
     
        writer.writerow([''])
        writer.writerow(['Inlet Bulk Concentration (mol/kg and g/cm^3)'])
    #     writer.writerows(InletBulkConcentration)
    #     writer.writerow([''])
        writer.writerows(InletBulkConcentration_gcm3)
         
    #     writer.writerow([''])
    #     writer.writerow(['Outlet Oxide (g/cm2)'])
    #     writer.writerows(OutletOxide)
          
        writer.writerow([''])
        writer.writerow(['Outlet S/O Concentration (g/cm^3)'])
        writer.writerows(OutletSOConcentration_gcm3)
        writer.writerow([''])
         
    #     writer.writerow(['Outlet Solubility (g/cm^3)'])
    #     writer.writerow(OutletSolubility_gcm3)
         
    #     writer.writerow([''])
        writer.writerow(['Dissolution Rate (g/cm^2 s)'])
        writer.writerows(DissolutionRate)
    #     writer.writerow([''])
         
    #     writer.writerow([''])
    #     writer.writerow(['Inner Oxide Growth from Corrosion (g/cm^2 s)'])
    #     writer.writerows(InnerOxideGrowth)
         
        writer.writerow([''])
        writer.writerow(['Growth - Dissolution (g/cm^2 s)'])
        writer.writerows(DeltaOx)
        writer.writerow([''])
    #     writer.writerow(['SG Oxide (g/cm2)'])
    #     writer.writerows(SGOxide)
    #     writer.writerow([''])
    

def RIHT_csv(InletFeeder, FuelChannel, OutletFeeder, SteamGenerator, FileName):
    #would have to select tubes from SG2 as well (in SG module)
    InnerOxide_SteamGeneratorTubes = []
    OuterOxide_SteamGeneratorTubes = []
    TotalOxide_SteamGeneratorTubes = []
    TotalDistance = []
    TemperatureProfile = []
    kp_Tdependent = []
    OutletCorrosionRate_um = [] 
    
#     FeSolubility_SteamGeneratorTubes = []
#     FeConcentration_SteamGeneratorTubes = []

    if OutletFeeder == ld.OutletFeeder_2:
        Output = pht_model.output_1
    elif OutletFeeder == ld.OutletFeeder:
        Output = pht_model.output_2 
    
    SelectedTubes = SGHX.tube_picker(SGHX.Method, SteamGenerator)[0]
    
    for Tube in SelectedTubes:
        x = ld.UnitConverter(
        Tube, "Grams per Cm Squared", "Grams per M Squared", None, None, Tube.InnerIronOxThickness, None, None, None
        )
        y = ld.UnitConverter(
        Tube, "Grams per Cm Squared", "Grams per M Squared", None, None, Tube.OuterFe3O4Thickness, None, None, None
        )
        z = [i / 100 for i in Tube.Distance]
        q = [i + j for i, j in zip(x, y)]
#         d = ld.UnitConverter(Tube, "Mol per Kg", "Grams per Cm Cubed", Tube.SolutionOxide.FeSatFe3O4, None, None, None,
#                              nc.FeMolarMass, None)
#         e = ld.UnitConverter(Tube, "Mol per Kg", "Grams per Cm Cubed", Tube.SolutionOxide.FeTotal, None, None, None,
#                              nc.FeMolarMass, None)
        f = Tube.KpFe3O4electrochem
            
        TotalDistance.append(z)
        InnerOxide_SteamGeneratorTubes.append(x)
        OuterOxide_SteamGeneratorTubes.append(y)
        TotalOxide_SteamGeneratorTubes.append(q)
        
    #     FeSolubility_SteamGeneratorTubes.append(d) # Tube.SolutionOxide.FeSatFe3O4
    #     FeConcentration_SteamGeneratorTubes.append(e) # Tube.SolutionOxide.FeTotal
        Temperature_C = [i - 273.15 for i in Tube.PrimaryBulkTemperature]
        TemperatureProfile.append(Temperature_C)
        kp_Tdependent.append(f)

    
    OutletCorrosionRate = Output[0]
    RIHT = Output[2]
    pht_SteamFraction = Output[3]
    OutletTemperatures1 = Output[4]
    OutletTemperatures2 = Output[5]
    
    # for outlet feeder that connects to current steam generator
    for Rate in OutletCorrosionRate:
        x = ld.UnitConverter(
        OutletFeeder, "Corrosion Rate Grams", "Corrosion Rate Micrometers", None, Rate, None, None, None, None
        )
        OutletCorrosionRate_um.append(x)
    
    Years = []
    for i in range((pht_model.SimulationYears + 1) * 4):
        Years.append((i / 4) + SGHX.YearStartup)
    

    Data = [
        SGHX.TubeLengths, TotalDistance, InnerOxide_SteamGeneratorTubes, OuterOxide_SteamGeneratorTubes,
            TotalOxide_SteamGeneratorTubes, OutletCorrosionRate_um, TemperatureProfile, kp_Tdependent
            ]
    Labels = [
        "U-bend length (cm)", "Distance (m)", "Inner Loading (g/m^2)", "Outer Loading (g/m^2)", "Total Oxide (g/m^2)",
        "Outlet Corrosion Rate (um/a)", "Temperature Profile (oC)", "kp (cm/s"]
     
    RIHT_delta = [x-y for x, y in zip (RIHT[1:], RIHT)]   
     
    csvfile = FileName
    with open(csvfile, "w") as output:
        writer = csv.writer(output, lineterminator='\n')
        writer.writerow(['RIHT (oC) and year'])
        writer.writerow(RIHT)
        writer.writerow(Years)
        writer.writerow(['Delta RIHT (oC)'])
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

RIHT_csv(ld.InletFeeder, ld.FuelChannel, ld.OutletFeeder_2, ld.SteamGenerator_2, "RIHTOutputSG2.csv")

if pht_model.Loop == "full":
    RIHT_csv(ld.InletFeeder_2, ld.FuelChannel_2, ld.OutletFeeder, ld.SteamGenerator, "RIHTOutputSG1.csv")  
else:
    None


LoopDistance = []
SectionsHalfLoop = [
    pht_model.InletFeeder_1_Loop1.Section1, pht_model.FuelChannel_1_Loop1.Section1,
    pht_model.OutletFeeder_2_Loop1.Section1, pht_model.SteamGeneratorTube_2_Loop1.Section1
    ] 

if pht_model.Loop == "half":
    Sections = SectionsHalfLoop
else:
    Sections = (
        SectionsHalfLoop + [
            pht_model.InletFeeder_2_Loop1.Section1, pht_model.FuelChannel_2_Loop1.Section1,
            pht_model.OutletFeeder_1_Loop1.Section1, pht_model.SteamGeneratorTube_1_Loop1.Section1
            ]
                )
           
for Section in Sections:
    x = Section.Length.magnitude
    LoopDistance.append(x)

LoopDistance = [j for i in LoopDistance for j in i]
TotalLoopDistance = [i / 100 for i in np.cumsum(LoopDistance)]  # Distance down length of PHTS [m]

 
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
            
        x = Concentrations[0]
        y = Concentrations[1]
        z = Concentrations[2]
    
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