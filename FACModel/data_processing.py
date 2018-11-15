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
import pandas as pd
from datetime import date, timedelta, datetime

import matplotlib.pyplot as plt
from matplotlib import rc
rc('mathtext', default='regular')

PlotOutput = "yes"


def purification_csv(InletFeeder, FuelChannel, OutletFeeder, SteamGenerator):
    Output = pht_model.output
    
    OutletCorrosionRate = Output[0]
    Time = Output[4]
#     InletBulkConcentration = Output[5]
#     OutletBulkConcentration = Output[6]
#     FuelChannelBulkConcentration = Output[7]
#     SteamGeneratorBulkConcentration = Output[8]
#     InletSolubility = Output[9]
    
    OutletCorrosionRate_uma = []
    InletBulkConcentration_gcm3 = []
    OutletBulkConcentration_gcm3 = []
    FuelChannelBulkConcentration_gcm3 = []
#     SteamGeneratorBulkConcentration_gm3 = []
#     InletSolubility_gm3 = []
    AvgCorrRate = []
    
    for Rate in OutletCorrosionRate:
        x = ld.UnitConverter(
        OutletFeeder, "Corrosion Rate Grams", "Corrosion Rate Micrometers", None, Rate, None, None, None, None
        )
        q = sum(x) / OutletFeeder.NodeNumber
        
        OutletCorrosionRate_uma.append(x)
        AvgCorrRate.append(q)
   
    for Conc1, Conc2, Conc3 in zip(
        InletBulkConcentration, OutletBulkConcentration, FuelChannelBulkConcentration):
        
        y = ld.UnitConverter(
            OutletFeeder, "Mol per Kg", "Grams per Cm Cubed", Conc1, None, None, None, nc.FeMolarMass, None
            )
         
        z = ld.UnitConverter(
            OutletFeeder, "Mol per Kg", "Grams per Cm Cubed", Conc2, None, None, None, nc.FeMolarMass, None
            )
    
        
        w = ld.UnitConverter(
            OutletFeeder, "Mol per Kg", "Grams per Cm Cubed", Conc3, None, None, None, nc.FeMolarMass, None
            )
        
#         p = ld.UnitConverter(
#             OutletFeeder, "Mol per Kg", "Grams per Cm Cubed", Conc4, None, None, None, nc.FeMolarMass, None
#             )
#         
#         u = ld.UnitConverter(
#             OutletFeeder, "Mol per Kg", "Grams per Cm Cubed", Conc5, None, None, None, nc.FeMolarMass, None
#             )
        
         
        
        InletBulkConcentration_gcm3.append(y)
        OutletBulkConcentration_gcm3.append(z)
        FuelChannelBulkConcentration_gcm3.append(w)
#         SteamGeneratorBulkConcentration_gm3.append(p)
#         InletSolubility_gm3.append(u)
        
    
#     DissolutionRate = []
#     InnerOxideGrowth = []
#     DeltaOx = []
#     for Concentration, Rate in zip(OutletSOConcentration_gcm3, OutletFeeder.CorrRate):
#         diss = [nc.KdFe3O4 * (x - y) for x, y in zip (OutletSolubility_gcm3, Concentration)]
#         growth = [i * it.Diffusion(Section, "Fe") / Section.FractionFeInnerOxide for i in Rate]
#         delta = [x - y for x, y in zip(growth, diss)]
#         DissolutionRate.append(diss)
#         InnerOxideGrowth.append(growth)
#         DeltaOx.append(delta)
     
    csvfile = "PurificationOutput.csv"
    with open(csvfile, "w") as output:
        writer = csv.writer(output, lineterminator='\n')
        writer.writerow(['Outlet Corrosion Rate (g/cm^2 s and um/a)'])
        writer.writerows(OutletCorrosionRate_uma)
    #     writer.writerows(OutletCorrosionRate)
         
        writer.writerow([''])
        writer.writerow(['Time (s)'])
        writer.writerow(Time)
        writer.writerow(['Average Corrosion Rate (um/a)'])
        writer.writerow(AvgCorrRate)
     
        writer.writerow([''])
        writer.writerow(['Inlet Bulk Concentration (g/cm^3)'])
        writer.writerows(InletBulkConcentration_gcm3)
        
#         writer.writerow([''])
#         writer.writerow(['Inlet Bulk Concentration (mol/kg)'])
#         writer.writerows(InletBulkConcentration)
    #     writer.writerow([''])
    #     writer.writerow(['Outlet Oxide (g/cm2)'])
    #     writer.writerows(OutletOxide)
        
        writer.writerow([''])
        writer.writerow(['Fuel Channel Bulk Concentration (g/cm^3)'])
        writer.writerows(FuelChannelBulkConcentration_gcm3)
        writer.writerow([''])
        
        writer.writerow([''])
        writer.writerow(['Outlet Bulk Concentration (g/cm^3)'])
        writer.writerows(OutletBulkConcentration_gcm3)
        writer.writerow([''])
        
#         writer.writerow(['SG Bulk Concentration (g/cm^3)'])
#         writer.writerows(SteamGeneratorBulkConcentration_gm3)
#         writer.writerow([''])
        
#         writer.writerow(['Outlet Bulk Concentration (mol/kg)'])
#         writer.writerows(OutletBulkConcentration)
#         writer.writerow([''])
         
#         writer.writerow(['Inlet Solubility (g/cm^3)'])
#         writer.writerows(InletSolubility_gm3)
#         writer.writerow(['Outlet Solubility (mol/kg)'])
#         writer.writerow(OutletFeeder.SolutionOxide.FeSatFe3O4)
         

def for_input_file_csv(InletFeeder, FuelChannel, OutletFeeder, SteamGenerator, FileName1):
    
    #would have to select tubes from SG2 as well (in SG module)
    InnerOxide = []
    OuterOxide = []
#     TotalOxide = []
#     TotalArea = []
    TotalGrams = []
    CorrosionRate = []
    kp_Tdependent = []
    CorrosionRate_um = []
    SG_grams_all_tubes = []
    TubeSize = []
    
    FeSolubility_SteamGeneratorTubes = []
    FeConcentration = []
    Output = pht_model.output
            
    AllSections = [InletFeeder, FuelChannel, OutletFeeder] + SteamGenerator
            
    for Section in AllSections:
        q = Section.InnerIronOxLoading
        w = Section.OuterFe3O4Loading
        r = Section.SolutionOxide.FeTotal
        
        if Section in SteamGenerator:
            x = [x + y for x, y in zip(Section.InnerIronOxLoading, Section.OuterFe3O4Loading)] # total oxide
            y = [np.pi * i * j for i, j in zip(Section.Length.magnitude, Section.Diameter)] # pi * D* L
            z = [i * j for i, j in zip(x, y)] # grams
            bundle_grams = [i * Section.TubeNumber for i in z]
            f = Section.TubeNumber
            TubeSize.append(f)
            k = sum(bundle_grams)
            SG_grams_all_tubes.append(k)

        InnerOxide.append(q)
        OuterOxide.append(w)
        FeConcentration.append(r)
#         CorrosionRate.append(r)
#         TotalOxide.append(x)
#         TotalGrams.append(z) 
        
    DividerPlateLeakage = Output[1]
    x_pht = Output[2]
    Date = Output[3]
    
    
    AverageTube_Fe_SO_Concentration = ld.UnitConverter(
        SteamGenerator[SGHX.Default_Tube], "Mol per Kg", "Grams per Cm Cubed",
        SteamGenerator[SGHX.Default_Tube].SolutionOxide.FeTotal, None, None, None, nc.FeMolarMass, None
        )
    
    AverageTube_BulkTemp_C = ld.UnitConverter(SteamGenerator[SGHX.Default_Tube], "Kelvin", "Celsius", None, None, None,
                                              None, None, SteamGenerator[SGHX.Default_Tube].PrimaryBulkTemperature
        )
    
    AverageTube_Solubility = ld.UnitConverter(
        SteamGenerator[SGHX.Default_Tube], "Mol per Kg", "Grams per Cm Cubed",
        SteamGenerator[SGHX.Default_Tube].SolutionOxide.FeSatFe3O4, None, None, None, nc.FeMolarMass, None
        )
    
    Data = [InnerOxide, OuterOxide]    
    Labels = ["Inner Oxide", "Outer Oxide"]
      
    csvfile = FileName1
    with open(csvfile, "w") as output:
        writer = csv.writer(output, lineterminator='\n')
        
        writer.writerow(['Start Date'])
        writer.writerow([Date[0]])
        writer.writerow([''])
        
        writer.writerow(['End Date'])
        writer.writerow([Date[len(Date) - 1]])
        writer.writerow([''])
        
        writer.writerow(['Bundle Size'])
        writer.writerow(TubeSize)
        writer.writerow([''])
         
        for i, j in zip(Labels, Data):
            writer.writerow([i])
            writer.writerows(j)
            writer.writerow([''])
        
        writer.writerow(['Total SG Oxide per Bundle (g)'])
        writer.writerow(SG_grams_all_tubes)
        writer.writerow([''])
         
        writer.writerow(['Sludge Loading (g/cm^2)'])
        writer.writerow(SteamGenerator[pht_model.Default_Tube].SludgeLoading)
        writer.writerow([''])
        
        writer.writerow(['Divider Plate Leakage'])
        writer.writerow([DividerPlateLeakage])
        writer.writerow([''])
        
        writer.writerow(['SteamQuality'])
        writer.writerow([x_pht])
        writer.writerow([''])
        
        writer.writerow(['Fe Tot (S/O) (g/cm^3)'])
        writer.writerows(FeConcentration)
        writer.writerow([''])
        
        writer.writerow(['Average Tube Temp Profile (oC)'])
        writer.writerow(AverageTube_BulkTemp_C)
        writer.writerows([''])
        
        writer.writerow(['Average Tube Concentration Profile'])
        writer.writerow(AverageTube_Fe_SO_Concentration)
        writer.writerows([''])
        
        writer.writerow(['Average Tube Distance'])
        writer.writerow(SteamGenerator[SGHX.Default_Tube].Distance)
        writer.writerows([''])
        
        writer.writerow(['SG Solublity'])
        writer.writerow(AverageTube_Solubility)
      
    return None   


# if loop run in half configuration, these will be identical
for_input_file_csv(ld.InletFeeder, ld.FuelChannel, ld.OutletFeeder_2, ld.SteamGenerator_2, "OutputSG2.csv")
# purification_csv(ld.InletFeeder, ld.FuelChannel, ld.OutletFeeder_2, ld.SteamGenerator_2)

if pht_model.Loop == "full":
    for_input_file_csv(
        ld.InletFeeder_2, ld.FuelChannel_2, ld.OutletFeeder, ld.SteamGenerator, "OutputSG1.csv")
    purification_csv(ld.InletFeeder_2, ld.FuelChannel_2, ld.OutletFeeder, ld.SteamGenerator) 


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
            x = Section.InnerOxLoading
        elif Layer == "Outer":
            x = Section.OuterOxLoading
        elif Layer == "Cobalt":
            x = Section.CoLoading
        elif Layer == "Nickel":
            x = Section.NiLoading
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

     
def concentration_and_loading_plots():   
    
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
    concentration_and_loading_plots()