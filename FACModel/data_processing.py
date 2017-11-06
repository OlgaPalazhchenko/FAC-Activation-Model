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
from matplotlib import rc
rc('mathtext', default='regular')

RealTimeHeatTransfer = "no"
Activation ="no"
PlotOutput = "no"

AverageColdLegLoading = []
RIHT = []
StreamOutletTemperatures = []

SimulationYears = 1  # years
SimulationHours = SimulationYears * 100

import time
start_time = time.time()

for j in range(SimulationHours):
    if RealTimeHeatTransfer == "no":
    
        In = pht_model.PHT_FAC(ld.Inlet, ld.Core, RealTimeHeatTransfer, j)
        Co = pht_model.PHT_FAC(ld.Core, ld.Outlet, RealTimeHeatTransfer, j)
        Ou = pht_model.PHT_FAC(ld.Outlet, ld.SteamGenerator, RealTimeHeatTransfer, j)
        Sg = pht_model.PHT_FAC(ld.SteamGenerator, ld.Inlet, RealTimeHeatTransfer, j)
    
    elif RealTimeHeatTransfer =="yes":
        In = pht_model.PHT_FAC(ld.Inlet, ld.Core, RealTimeHeatTransfer, j)
        Co = pht_model.PHT_FAC(ld.Core, ld.Outlet, RealTimeHeatTransfer, j)
        Ou = pht_model.PHT_FAC(ld.Outlet, ld.SteamGenerator, RealTimeHeatTransfer, j)
        Sg = pht_model.PHT_FAC(ld.SteamGenerator, ld.Inlet, RealTimeHeatTransfer, j)
        
        #Set input concentrations for all SG zones to be same as input of first (Outlet output)
        for Zone in ld.SGZones[1:len(ld.SGZones)]:
            BulkOutletOutput= [
                Ou.Section1.Bulk.FeTotal, Ou.Section1.Bulk.NiTotal, Ou.Section1.Bulk.CoTotal, Ou.Section1.Bulk.CrTotal
                ]
            BulkSGInput = [Zone.Bulk.FeTotal, Zone.Bulk.NiTotal, Zone.Bulk.CoTotal, Zone.Bulk.CrTotal]
            for x,y in zip(BulkSGInput,BulkOutletOutput):
                x[0] = y[Ou.Section1.NodeNumber-1]
                
        Sg = pht_model.PHT_FAC(ld.SteamGenerator, ld.Inlet, j)
        Sg1 = pht_model.PHT_FAC(ld.SG_Zone1, ld.Inlet, j)
        Sg2 = pht_model.PHT_FAC(ld.SG_Zone2, ld.Inlet, j)
        Sg3 = pht_model.PHT_FAC(ld.SG_Zone3, ld.Inlet, j)

    if j % 20 == 0:  # yearly
        
        TotalLoadingSG = [x + y for x, y in zip(Sg.Section1.OuterOxThickness, Sg.Section1.InnerOxThickness)]
        ConvertedLoading = ld.UnitConverter(
            Sg.Section1, "Grams per Cm Squared", "Grams per M Squared", None, None, TotalLoadingSG, None, None, None
            )
        AverageColdLeg = sum(ConvertedLoading[11:len(ConvertedLoading)]) \
            / (len(ConvertedLoading) - 11)
  
        T_x = (SGHX.energy_balance(21, j)) - 273.15
        
        OutletTemperatures = []
        for Zone in ld.SGZones:
            Temperature = Zone.PrimaryBulkTemperature[Zone.NodeNumber - 1] - 273.15
            OutletTemperatures.append(Temperature)
        
        StreamOutletTemperatures.append(OutletTemperatures)
        AverageColdLegLoading.append(AverageColdLeg)
        RIHT.append(T_x)
            

FACRate = sum(
    ld.UnitConverter(Ou.Section1, "Corrosion Rate Grams", "Corrosion Rate Micrometers", None, Ou.Section1.CorrRate,
    None, None, None, None)
    ) / (Ou.Section1.NodeNumber)

csvfile = "RIHTOutput.csv"
with open(csvfile, "w") as output:
    writer = csv.writer(output, lineterminator='\n')
    writer.writerow(["average cold leg loading"])
    writer.writerows([AverageColdLegLoading])

    writer.writerow(['Outlet Streams'])
    writer.writerows(StreamOutletTemperatures)
    
    writer.writerow(['RIHT'])
    writer.writerow(RIHT) 
    
    writer.writerow(['Inner Oxide Layer [g/cm^2]'])
    writer.writerow(Sg.Section1.InnerOxThickness)
    writer.writerow(['Outer Oxide Layer [g/cm^2]'])
    writer.writerow(Sg.Section1.OuterOxThickness)   

    # oxide profiel for all sections
    # if RealTimeHeatTransfer == "yes":
        

end_time = time.time()
delta_time = end_time - start_time

hours = delta_time // 3600
temp = delta_time - 3600 * hours
minutes = delta_time // 60
seconds = delta_time - 60 * minutes
print('%d:%d:%d' % (hours, minutes, seconds))




LoopDistance = (In.Section1.Length.magnitude
                + Co.Section1.Length.magnitude
                + Ou.Section1.Length.magnitude
                + Sg.Section1.Length.magnitude)
TotalLoopDistance = [i / 100 for i in np.cumsum(LoopDistance)]  # Distance down length of PHTS [m]
 
BulkFeSat = [np.log10(i) for i in
             In.Section1.Bulk.FeSatFe3O4
             + Co.Section1.Bulk.FeSatFe3O4
             + Ou.Section1.Bulk.FeSatFe3O4
             + Sg.Section1.Bulk.FeSatFe3O4]
BulkFeTotal = [np.log10(i) for i in
               In.Section1.Bulk.FeTotal
               + Co.Section1.Bulk.FeTotal
               + Ou.Section1.Bulk.FeTotal
               + Sg.Section1.Bulk.FeTotal]
SolutionOxideFeTotal = [np.log10(i) for i in In.Section1.SolutionOxide.FeTotal + Co.Section1.SolutionOxide.FeTotal + \
                        Ou.Section1.SolutionOxide.FeTotal + Sg.Section1.SolutionOxide.FeTotal]
SolutionOxideFeSat = [np.log10(i) for i in
                      In.Section1.SolutionOxide.FeSatFe3O4
                      + Co.Section1.SolutionOxide.FeSatFe3O4
                      + Ou.Section1.SolutionOxide.FeSatFe3O4
                      + Sg.Section1.SolutionOxide.FeSatFe3O4]
 
BulkNiTotal = [np.log10(i) for i in
               In.Section1.Bulk.NiTotal
               + Co.Section1.Bulk.NiTotal
               + Ou.Section1.Bulk.NiTotal + Sg.Section1.Bulk.NiTotal]
SolutionOxideNiTotal = [np.log10(i) for i in
                        In.Section1.SolutionOxide.NiTotal + Co.Section1.SolutionOxide.NiTotal + \
                        Ou.Section1.SolutionOxide.NiTotal + Sg.Section1.SolutionOxide.NiTotal]
SolutionOxideNiSat = [np.log10(i) for i in
                      In.Section1.SolutionOxide.NiSatFerrite + Co.Section1.SolutionOxide.NiSatFerrite + \
                        Ou.Section1.SolutionOxide.NiSatFerrite + Sg.Section1.SolutionOxide.NiSatFerrite]
 
BulkCoTotal = [np.log10(i) for i in
               In.Section1.Bulk.CoTotal
               + Co.Section1.Bulk.CoTotal
               + Ou.Section1.Bulk.CoTotal
               + Sg.Section1.Bulk.CoTotal]
SolutionOxideCoTotal = [np.log10(i) for i in
                        In.Section1.SolutionOxide.CoTotal + Co.Section1.SolutionOxide.CoTotal + \
                        Ou.Section1.SolutionOxide.CoTotal + Sg.Section1.SolutionOxide.CoTotal]
SolutionOxideCoSat = [np.log10(i) for i in
                      In.Section1.SolutionOxide.CoSatFerrite
                      + Co.Section1.SolutionOxide.CoSatFerrite
                      + Ou.Section1.SolutionOxide.CoSatFerrite
                      + Sg.Section1.SolutionOxide.CoSatFerrite]
 
 
BulkCrTotal = [np.log10(i) for i in
               In.Section1.Bulk.CrTotal
               + Co.Section1.Bulk.CrTotal
               + Ou.Section1.Bulk.CrTotal
               + Sg.Section1.Bulk.CrTotal]
 
SolutionOxideCrSat = [np.log10(i) for i in
                      In.Section1.SolutionOxide.CrSat
                      + Co.Section1.SolutionOxide.CrSat
                      + Ou.Section1.SolutionOxide.CrSat
                      + Sg.Section1.SolutionOxide.CrSat]
 
InnerOxide = []
InnerOxideThicknesses = [
    In.Section1.InnerOxThickness, Co.Section1.InnerOxThickness, Ou.Section1.InnerOxThickness, 
    Sg.Section1.InnerOxThickness
    ]
for Thickness, Sect in zip (InnerOxideThicknesses, ld.Sections):
    z = ld.UnitConverter(Sect, "Grams per Cm Squared", "Grams per M Squared", None, None, Thickness, None, None, None)
    InnerOxide.append(z)
InnerOxideThickness = InnerOxide[0] + InnerOxide[1] + InnerOxide[2] + InnerOxide[3]
 
OuterOxide = []
OuterOxideThicknesses = [
    In.Section1.OuterOxThickness, Co.Section1.OuterOxThickness, Ou.Section1.OuterOxThickness, 
    Sg.Section1.OuterOxThickness
    ]
for Thickness, Sect in zip (OuterOxideThicknesses, ld.Sections):
    q = ld.UnitConverter(Sect, "Grams per Cm Squared", "Grams per M Squared", None, None, Thickness, None, None, None)
    OuterOxide.append(q)
OuterOxideThickness = OuterOxide[0] + OuterOxide[1] + OuterOxide[2] + OuterOxide[3]
 
Cobalt = []
CoThicknesses = [In.Section1.CoThickness, Co.Section1.CoThickness, Ou.Section1.CoThickness, Sg.Section1.CoThickness]
for Thickness, Sect in zip (CoThicknesses, ld.Sections):
    c = ld.UnitConverter(Sect, "Grams per Cm Squared", "Grams per M Squared", None, None, Thickness, None, None, None)
    Cobalt.append(c)
CobaltThickness = Cobalt[0] + Cobalt[1] + Cobalt[2] + Cobalt[3]
 
Nickel = []
NiThicknesses = [In.Section1.NiThickness, Co.Section1.NiThickness, Ou.Section1.NiThickness, Sg.Section1.NiThickness]
for Thickness, Sect in zip (NiThicknesses, ld.Sections):
    n = ld.UnitConverter(Sect, "Grams per Cm Squared", "Grams per M Squared", None, None, Thickness, None, None, None)
    Nickel.append(n)
NickelThickness = Nickel[0] + Nickel[1] + Nickel[2] + Nickel[3]

if PlotOutput == "yes":
    fig1 = plt.figure()
    ax1 = fig1.add_subplot(221)
    ax1.plot(TotalLoopDistance, SolutionOxideFeTotal, marker='o', color='r', label='S/O Interface')
    ax1.plot(TotalLoopDistance, SolutionOxideFeSat, marker='o', color='0', label='Magnetite Solubility')
    ax1.plot(TotalLoopDistance, BulkFeTotal, marker='o', color='0.50', label='Bulk Coolant')
    ax1.set_xlabel('Distance (m)')
    ax1.set_ylabel('$Log_{10}$[Fe Concentration (mol/kg)]')
    # ax1.legend()
 
    ax2 = fig1.add_subplot(222)
    ax2.plot(TotalLoopDistance, SolutionOxideNiTotal, marker='o', color='b', label='S/O Interface')
    ax2.plot(TotalLoopDistance, SolutionOxideNiSat, marker='o', color='0', label='Ni-Ferrite Solubility')
    ax2.plot(TotalLoopDistance, BulkNiTotal, marker='o', color='0.50', label='Bulk Coolant') 
    ax2.set_xlabel('Distance (m)')
    ax2.set_ylabel('$Log_{10}$[Ni Concentration (mol/kg)]')
    # ax2.legend()
 
    ax3 = fig1.add_subplot(223)
    ax3.plot(TotalLoopDistance, SolutionOxideCoTotal, marker='^', color='m', label='S/O Interface')
    ax3.plot(TotalLoopDistance, SolutionOxideCoSat, marker='^', color='0', label='Co-Ferrite Solubility')
    ax3.plot(TotalLoopDistance, BulkCoTotal, marker='^', color='0.50', label='Bulk Coolant') 
    ax3.set_xlabel('Distance (m)')
    ax3.set_ylabel('$Log_{10}$[Co Concentration (mol/kg)]')
    # ax3.legend()
 
 
    ax4 = fig1.add_subplot(224)
    ax4.plot(TotalLoopDistance, SolutionOxideCrSat, marker='*', color='0', label='Chromite Solubility')
    ax4.plot(TotalLoopDistance, BulkCrTotal, marker='*', color='0.50', label='Bulk Coolant') 
    ax4.set_xlabel('Distance (m)')
    ax4.set_ylabel('$Log_{10}$[Cr Concentration (mol/kg)]')
     
    plt.tight_layout()
    plt.show()
     
    fig2, ax1 = plt.subplots()
    # ax1 = fig2.add_subplot(221)
    ax1.plot(TotalLoopDistance, InnerOxideThickness, linestyle=None, marker='o', color='0.50', label='Inner Oxide')
    ax1.plot(TotalLoopDistance, OuterOxideThickness, linestyle=None, marker='o', color='k', label='Outer Oxide')
    # ax1.axis([51,69, 0, 30])
    ax1.set_xlabel('Distance (m)')
    ax1.set_ylabel('Oxide Layer Loadings (${g/m^2}$)')
                
    ax2 = ax1.twinx()
    ax2.plot(TotalLoopDistance, NickelThickness, linestyle=None, marker='o', color='c', label='Nickel')
    ax2.plot(TotalLoopDistance, CobaltThickness, linestyle=None, marker='o', color='m', label='Cobalt')
    ax2.set_ylabel('Ni, Co, and Cr Loadings (${g/m^2}$)', rotation=270, labelpad=20)
    lines, labels = ax1.get_legend_handles_labels()
    lines2, labels2 = ax2.get_legend_handles_labels()
    ax2.legend(lines + lines2, labels + labels2, loc=0)        
    # plt.axis([51, 69, 0, 30])
    plt.tight_layout()
    plt.show()