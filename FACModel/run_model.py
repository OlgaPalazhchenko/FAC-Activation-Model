import numpy as np
import constants as nc
import lepreau_data as ld
import composition as c
import rk_4
import activities as a
import electrochemistry as e
import iteration as it
import sg_heattransfer as SGHX
import matplotlib.pyplot as plt
from matplotlib import rc
rc('mathtext', default='regular')

# #Initial Temperatures and T-dependent parametrs in SG zones
for Zone in ld.SGZones:
    Zone.PrimaryBulkTemperature = SGHX.temperature_profile(
        Zone, Zone.InnerOxThickness, Zone.OuterOxThickness, 
        SGHX.MassFlow_h.magnitude - 0.03 * SGHX.MassFlow_h.magnitude, 0
        )
    Zone.DensityH2O = [ld.Density("water", "PHT", i) for i in Zone.PrimaryBulkTemperature]
    Zone.ViscosityH2O = [ld.Viscosity("water", "PHT", i) for i in Zone.PrimaryBulkTemperature]
    Zone.NernstConstant = [x * (2.303 * nc.R / (2 * nc.F)) for x in Zone.PrimaryBulkTemperature]

# #Initial Concentrations
for Section in ld.Sections:    
    # #Temperature-dependent parameters            
    Section.NernstConstant = [x * (2.303 * nc.R / (2 * nc.F)) for x in Section.PrimaryBulkTemperature]
    
    Section.DensityH2O = [ld.Density("water", "PHT", x) for x in Section.PrimaryBulkTemperature]
    Section.ViscosityH2O = [ld.Viscosity("water", "PHT", x) for x in Section.PrimaryBulkTemperature]
    
    Section.FractionFeInnerOxide = c.FractionMetalInnerOxide(Section, "Fe")
    Section.FractionNiInnerOxide = c.FractionMetalInnerOxide(Section, "Ni")
    
    [
        Section.Bulk.Co60, Section.Bulk.Co58, Section.Bulk.Fe55, Section.Bulk.Fe59, Section.Bulk.Mn54, 
        Section.Bulk.Cr51, Section.Bulk.Ni63
        ] = [[1e-10] * Section.NodeNumber] * 7
    
    Interfaces = [Section.MetalOxide, Section.SolutionOxide, Section.Bulk]
    for Interface in Interfaces:
        # makes into array with appropriate # of nodes for that PHTS section
        Interface.ConcentrationH2 = [nc.H2 * nc.H2Density / nc.H2MolarMass] * Section.NodeNumber
        Interface.ConcentrationH = c.bulkpH_calculator(Section)  # from system pH calculation
        
        # Concentration/Saturation Input [mol/kg]
        Interface.FeTotal = c.IronSolubility(Section)
        Interface.FeSatFe3O4 = [1 * i for i in c.IronSolubility(Section)]
        Interface.NiTotal = [i * 1 for i in Section.SolubilityNi]
        Interface.NiSatFerrite = [i * 1 for i in Section.SolubilityNi]
        Interface.NiSatMetallicNi = [i * 1 for i in Section.SolubilityNi]
    
        Interface.CoTotal = [i * 1 for i in Section.SolubilityCo]
        Interface.CoSatFerrite = [i * 1 for i in Section.SolubilityCo]
        Interface.CrTotal = [i * 1 for i in Section.SolubilityCr]
        Interface.CrSat = [i * 1 for i in Section.SolubilityCr]
        
        # Initial RIHT temperature 
        if Section == ld.Inlet:
            Section.PrimaryBulkTemperature = [SGHX.energy_balance(ld.SteamGenerator.NodeNumber - 1, 0)]\
            * Section.NodeNumber
        
        if Section != ld.Outlet:
            if Interface == Section.SolutionOxide:
                Interface.FeTotal = [i * 0.8 for i in c.IronSolubility(Section)]
        
        if Section == ld.Core and Interface == Section.MetalOxide:
            Interface.FeTotal = [0] * Section.NodeNumber
        
        if Section == ld.Outlet and Interface == Section.MetalOxide:
            # From Cook's thesis - experimental corrosion rate measurements and calcs
            Interface.FeTotal = [0.00000026] * Section.NodeNumber 
        
        if Section == ld.Inlet or Section == ld.Outlet or Section == ld.Core:
            Interface.NiTotal = [0] * Section.NodeNumber
        
    Section.SolutionOxide.MixedPotential, Section.SolutionOxide.EqmPotentialFe3O4 = e.ECP(Section)
    
    if Section == ld.Core:
        Section.CorrRate = [0] * Section.NodeNumber
    else:
        Section.CorrRate, Section.MetalOxide.MixedPotential = it.CorrosionRate(Section)
        
    Section.KpFe3O4electrochem, Section.KdFe3O4electrochem, Section.SolutionOxide.FeSatFe3O4,\
    Section.MetalOxide.ConcentrationH = e.ElectrochemicalAdjustment(
        Section, Section.SolutionOxide.EqmPotentialFe3O4, Section.SolutionOxide.MixedPotential,
        Section.MetalOxide.MixedPotential, Section.SolutionOxide.FeTotal, Section.SolutionOxide.FeSatFe3O4,
        Section.Bulk.FeSatFe3O4, Section.SolutionOxide.ConcentrationH
        )


class PHT_FAC():
    def __init__(self, ActiveSection, OutgoingSection, j):  # j = overall time step
        self.Section1 = ActiveSection
        self.Section2 = OutgoingSection
            
        # BulkActivities = [self.Section1.Bulk.Co60, self.Section1.Bulk.Co58, self.Section1.Bulk.Fe59, self.Section1.Bulk.Fe55, \
        #                  self.Section1.Bulk.Mn54, self.Section1.Bulk.Cr51, self.Section1.Bulk.Ni63]
        
        # BulkActivities2 = [self.Section2.Bulk.Co60, self.Section2.Bulk.Co58, self.Section2.Bulk.Fe59, self.Section2.Bulk.Fe55, \
        #                  self.Section2.Bulk.Mn54, self.Section2.Bulk.Cr51, self.Section2.Bulk.Ni63]
        
        # SurfaceActivities = [self.Section1.Bulk.Co60, self.Section1.Bulk.Co58, self.Section1.Bulk.Fe59, self.Section1.Bulk.Fe55, \
        #                  self.Section1.Bulk.Mn54, self.Section1.Bulk.Cr51, self.Section1.Bulk.Ni63]
        
        # Tags = ["Co60", "Co58", "Fe59", "Fe55", "Mn54", "Cr51", "Ni63"]
        
        BulkConcentrations = [
            self.Section1.Bulk.FeTotal, self.Section1.Bulk.NiTotal, self.Section1.Bulk.CoTotal, 
            self.Section1.Bulk.CrTotal
            ]
        BulkConcentrations2 = [
            self.Section2.Bulk.FeTotal, self.Section2.Bulk.NiTotal, self.Section2.Bulk.CoTotal, 
            self.Section2.Bulk.CrTotal
            ]
        
        SolutionOxideConcentrations = [
            self.Section1.SolutionOxide.FeTotal, self.Section1.SolutionOxide.NiTotal,
            self.Section1.SolutionOxide.CoTotal, self.Section1.SolutionOxide.CrTotal
            ]
        Saturations = [
            self.Section1.SolutionOxide.FeSatFe3O4, self.Section1.SolutionOxide.NiSatFerrite, 
            self.Section1.SolutionOxide.CoSatFerrite
            ]
        
        for i in range(self.Section1.NodeNumber):
            if i > 0: 
                for r, q in zip (BulkConcentrations, SolutionOxideConcentrations):
                    r[i] = rk_4.Spatial(
                        q[i - 1], r[i - 1], ld.MassTransfer(self.Section1)[i - 1], self.Section1.Diameter[i - 1],
                        self.Section1.Velocity[i - 1], self.Section1.Length.magnitude[i - 1]
                        )

                # Exponential decay of bulk particulate at start of section as function of distance + removal due to 
                # deposition and erosion source
                self.Section1.BigParticulate[i] = rk_4.particulate(
                    self.Section1, self.Section1.BigParticulate[0], self.Section1.Diameter[i], 
                    self.Section1.DensityH2O[i], self.Section1.Velocity[i], self.Section1.Distance[i]
                    )
                self.Section1.SmallParticulate[i] = rk_4.particulate(
                    Section, self.Section1.SmallParticulate[0], self.Section1.Diameter[i], self.Section1.DensityH2O[i],
                    self.Section1.Velocity[i], self.Section1.Distance[i]
                    )
                                
#                 for x,y in zip (BulkActivities, Tags):
#                     x[i] = a.BulkActivity(self.Section1, self.Section2, x[0], self.Section1.CorrRate, self.Section1.Bulk.FeSatFe3O4, \
#                                 self.Section1.InnerOxThickness, self.Section1.OuterOxThickness, self.Section1.OuterFe3O4Thickness, \
#                                 self.Section1.NiThickness, self.Section1.CoThickness, self.Section1.InnerIronOxThickness, \
#                                 y, self.Section1.BigParticulate, self.Section1.SmallParticulate, j, i)
                    
                # Inlet header purification system
                if self.Section1 == ld.Inlet:      
                    if i == 3:
                        # for x,y in zip(BulkConcentrations, BulkActivities):
                        for x in BulkConcentrations: 
                            x[i] = 0.59 * x[i - 1]
                            # y[i] = 0.59*y[i-1]
                    
                        self.Section1.SmallParticulate[i] = 0.59 * self.Section1.SmallParticulate[i - 1]
                        self.Section1.BigParticulate[i] = 0  # 0.45 um filter removes everything over this size   
                    
                    if i > 3:  # decay for the rest of the inlet section depends on purification system 
                        self.Section1.BigParticulate[i] = rk_4.Particulate(
                            Section, self.Section1.BigParticulate[3], self.Section1.Diameter[i],
                            self.Section1.DensityH2O[i], self.Section1.Velocity[i], self.Section1.Distance[i]
                            )
                        self.Section1.SmallParticulate[i] = rk_4.Particulate(
                            Section, self.Section1.SmallParticulate[3], self.Section1.Diameter[i],
                            self.Section1.DensityH2O[i], self.Section1.Velocity[i], self.Section1.Distance[i]
                            ) 
                        
#                         for x,y in zip (BulkActivities, Tags):
#                             #print (y, x[i], i, self.Section1.NodeNumber)
#                             x[i] = a.BulkActivity(self.Section1, self.Section2, x[3], self.Section1.CorrRate, self.Section1.Bulk.FeSatFe3O4, \
#                                     self.Section1.InnerOxThickness, self.Section1.OuterOxThickness, self.Section1.OuterFe3O4Thickness, \
#                                     self.Section1.NiThickness, self.Section1.CoThickness, self.Section1.InnerIronOxThickness, y, \
#                                     self.Section1.BigParticulate, self.Section1.SmallParticulate, j, i)
                                        
        end = self.Section1.NodeNumber - 1
        for x, y in zip(BulkConcentrations, BulkConcentrations2):
        # for x,y,z,q in zip(BulkConcentrations,BulkConcentrations2, BulkActivities, BulkActivities2):
            y[0] = x[end]
            # q[0] = z[end]
        
        self.Section2.BigParticulate[0] = self.Section1.BigParticulate[end]
        self.Section2.SmallParticulate[0] = self.Section1.SmallParticulate[end]       
        
        # #Stellite wear bulk input term for cobalt and chromium    
        if self.Section1 == ld.Core:
            self.Section2.Bulk.CoTotal[0] = self.Section1.Bulk.CoTotal[end] + nc.CobaltWear
            self.Section2.Bulk.CrTotal[0] = self.Section1.Bulk.CrTotal[end] + nc.CobaltWear * (nc.FractionCr_Stellite / nc.FractionCo_Stellite)
   
#         else:
#             ##Surface activities 
#             self.Section1.MetalOxide.Co60 = a.SurfaceActivity(self.Section1, self.Section1.CorrRate, self.Section1.SolutionOxide.FeSatFe3O4, \
#                                                           self.Section1.InnerOxThickness, self.Section1.Bulk.Co60, j, "Co60")
#             
#             self.Section1.MetalOxide.Co58 = a.SurfaceActivity(self.Section1, self.Section1.CorrRate, self.Section1.SolutionOxide.FeSatFe3O4, \
#                                                           self.Section1.InnerOxThickness, self.Section1.Bulk.Co58, j, "Co58")
#             
#             self.Section1.MetalOxide.Fe59 = a.SurfaceActivity(self.Section1, self.Section1.CorrRate, self.Section1.SolutionOxide.FeSatFe3O4, \
#                                                           self.Section1.InnerOxThickness, self.Section1.Bulk.Fe59, j, "Fe59")
#             
#             self.Section1.MetalOxide.Fe55 = a.SurfaceActivity(self.Section1, self.Section1.CorrRate, self.Section1.SolutionOxide.FeSatFe3O4, \
#                                                           self.Section1.InnerOxThickness, self.Section1.Bulk.Fe55, j, "Fe55")
#             
#             self.Section1.MetalOxide.Mn54 = a.SurfaceActivity(self.Section1, self.Section1.CorrRate, self.Section1.SolutionOxide.FeSatFe3O4, \
#                                                           self.Section1.InnerOxThickness, self.Section1.Bulk.Mn54, j, "Mn54")
#             
#             self.Section1.MetalOxide.Cr51 = a.SurfaceActivity(self.Section1, self.Section1.CorrRate, self.Section1.SolutionOxide.FeSatFe3O4, \
#                                                           self.Section1.InnerOxThickness, self.Section1.Bulk.Cr51, j, "Cr51")
#             
#             self.Section1.MetalOxide.Ni63 = a.SurfaceActivity(self.Section1, self.Section1.CorrRate, self.Section1.SolutionOxide.FeSatFe3O4, \
#                                                           self.Section1.InnerOxThickness, self.Section1.Bulk.Ni63, j, "Ni63")
            # #    
        
        # #SG heat transfer 
#         if self.Section1 in ld.SGZones:   
#             self.Section1.PrimaryBulkTemperature = SGHX.temperature_profile(self.Section1, self.Section1.InnerOxThickness, self.Section1.OuterOxThickness)
#             self.Section1.Bulk.FeSatFe3O4 = c.IronSolubility(self.Section1) 
            
#         if self.Section1 == ld.SG_Zone1:    
#             if j %699==0: #(17520)
#                 print ([i-273.15 for i in self.Section1.PrimaryBulkTemperature],j)
#                 print (self.Section1.SolutionOxide.FeSatFe3O4, "FeSat S/O")
#                 print(self.Section1.SolutionOxide.FeTotal, "FeTot S/O")
#                 print()
        
        # RIHT  
#         if self.Section1 == ld.Inlet:
#             self.Section1.PrimaryBulkTemperature = SGHX.energy_balance(21, self.Section1.NodeNumber)
#             self.Section1.Bulk.FeSatFe3O4 = c.IronSolubility(self.Section1)

        
        # Deposit thickness around PHTS
        # self.Section1.DepositThickness = a.Deposition(self.Section1, self.Section1.BigParticulate, self.Section1.SmallParticulate, j)
        # #
                        
        # #rk_4 oxide thickness calculation (no spalling)
        rk_4.oxidegrowth(self.Section1, Saturations, BulkConcentrations)

                
        # Spalling    
        self.Section1.ElapsedTime, self.Section1.SpallTime = rk_4.Spall(
            self.Section1, j, self.Section1.ElapsedTime, self.Section1.SpallTime
            )

# TotalTime = []
# SACo60 = [] 
# SACo58 = [] 
# SAFe59 = []
# SAFe55 = []
# SAMn54 = []
# SACr51 = []
# SANi63 = []       

import time
start_time = time.time()
for j in range(8760 * 9):  # nc.SimulationDuration
    I = PHT_FAC(ld.Inlet, ld.Core, j)
    C = PHT_FAC(ld.Core, ld.Outlet, j)
    O = PHT_FAC(ld.Outlet, ld.SteamGenerator, j)
    SG = PHT_FAC(ld.SteamGenerator, ld.Inlet, j)
    if j % 8759 == 0:  # 8759
        totalthicknessSG = [x + y for x, y in zip(SG.Section1.OuterOxThickness, SG.Section1.InnerOxThickness)]
        convertedtotalthicknessSG = ld.UnitConverter(
            SG.Section1, "Grams per Cm Squared", "Grams per M Squared", None, None, totalthicknessSG, None, None, None
            )
  
        RIHT = (SGHX.energy_balance(21, j))
        print(RIHT - 273.15, "RIHT")
        print(ld.SteamGenerator.PrimaryBulkTemperature[21] - 273.15, ld.SG_Zone1.PrimaryBulkTemperature[21] - 273.15, \
                  ld.SG_Zone2.PrimaryBulkTemperature[21] - 273.15, ld.SG_Zone3.PrimaryBulkTemperature[21] - 273.15)
        print(SG.Section1.SolutionOxide.FeTotal, "FeTotal")
        print(SG.Section1.SolutionOxide.FeSatFe3O4, "FeSat")
        print (sum(convertedtotalthicknessSG[11:22]) / (22 - 11), "Average SG cold leg [g/m^2] oxide thickness", j, "time [h]")
#         print (sum(ld.UnitConverter(O.Section1, "Corrosion Rate Grams", "Corrosion Rate Micrometers", None, O.Section1.CorrRate, None, None, None, None))\
#                    /(O.Section1.NodeNumber), "Average outlet corrosion rate [um/a]")
        print()
#     TotalTime.append(j/7000)# converts to EFPY
#     SACo60.append(sum(S.Section1.MetalOxide.Co60)/S.Section1.NodeNumber)
#     SACo58.append(sum(S.Section1.MetalOxide.Co58)/S.Section1.NodeNumber)
#     SAFe59.append(sum(S.Section1.MetalOxide.Fe59)/S.Section1.NodeNumber)
#     SAFe55.append(sum(S.Section1.MetalOxide.Fe55)/S.Section1.NodeNumber)
#     SAMn54.append(sum(S.Section1.MetalOxide.Mn54)/S.Section1.NodeNumber)
#     SACr51.append(sum(S.Section1.MetalOxide.Cr51)/S.Section1.NodeNumber)
#     SANi63.append(sum(S.Section1.MetalOxide.Ni63)/S.Section1.NodeNumber)
    
end_time = time.time()
delta_time = end_time - start_time

hours = delta_time // 3600
temp = delta_time - 3600 * hours
minutes = delta_time // 60
seconds = delta_time - 60 * minutes
print('%d:%d:%d' % (hours, minutes, seconds))

# Corrosion= []
# Rates = [I.Section1.CorrRate, C.Section1.CorrRate, O.Section1.CorrRate, SG.Section1.CorrRate]
# for Rate, Section in zip (Rates, ld.Sections):
#     x = ld.UnitConverter(Section, "Corrosion Rate Grams", "Corrosion Rate Micrometers", None, Rate, None, None, None, None)
#     Corrosion.append(x)
# CorrosionRate = Corrosion[0] + Corrosion[1] + Corrosion[2] + Corrosion[3]
# print (CorrosionRate)
# print(-np.log10(O.Section1.SolutionOxide.ConcentrationH[0]), "Outlet pH")
# print([i-273.15 for i in SG.Section1.PrimaryBulkTemperature], "Temperature profile")
# print()
print (SG.Section1.InnerOxThickness, "InnerOx")
print(SG.Section1.OuterOxThickness, "OuterOx")

# ActivityTime = TotalTime[0::3]
# Co60SGActivity = SACo60[0::3]
# Co58SGActivity = SACo58[0::3]
# Fe59SGActivity = SAFe59[0::3]
# Fe55SGActivity = SAFe55[0::3]
# Mn54SGActivity = SAMn54[0::3]
# Cr51SGActivity = SACr51[0::3]
# Ni63SGActivity = SANi63[0::3]

LoopDistance = I.Section1.Length.magnitude + C.Section1.Length.magnitude + O.Section1.Length.magnitude + SG.Section1.Length.magnitude
TotalLoopDistance = [i / 100 for i in np.cumsum(LoopDistance)]  # Distance down length of PHTS [m]

BulkFeSat = [np.log10(i) for i in I.Section1.Bulk.FeSatFe3O4 + C.Section1.Bulk.FeSatFe3O4 + O.Section1.Bulk.FeSatFe3O4 + SG.Section1.Bulk.FeSatFe3O4]
BulkFeTotal = [np.log10(i) for i in I.Section1.Bulk.FeTotal + C.Section1.Bulk.FeTotal + O.Section1.Bulk.FeTotal + SG.Section1.Bulk.FeTotal]
SolutionOxideFeTotal = [np.log10(i) for i in I.Section1.SolutionOxide.FeTotal + C.Section1.SolutionOxide.FeTotal + \
                        O.Section1.SolutionOxide.FeTotal + SG.Section1.SolutionOxide.FeTotal]
SolutionOxideFeSat = [np.log10(i) for i in I.Section1.SolutionOxide.FeSatFe3O4 + C.Section1.SolutionOxide.FeSatFe3O4 + \
                        O.Section1.SolutionOxide.FeSatFe3O4 + SG.Section1.SolutionOxide.FeSatFe3O4]

BulkNiTotal = [np.log10(i) for i in I.Section1.Bulk.NiTotal + C.Section1.Bulk.NiTotal + O.Section1.Bulk.NiTotal + SG.Section1.Bulk.NiTotal]
SolutionOxideNiTotal = [np.log10(i) for i in I.Section1.SolutionOxide.NiTotal + C.Section1.SolutionOxide.NiTotal + \
                        O.Section1.SolutionOxide.NiTotal + SG.Section1.SolutionOxide.NiTotal]
SolutionOxideNiSat = [np.log10(i) for i in I.Section1.SolutionOxide.NiSatFerrite + C.Section1.SolutionOxide.NiSatFerrite + \
                        O.Section1.SolutionOxide.NiSatFerrite + SG.Section1.SolutionOxide.NiSatFerrite]

BulkCoTotal = [np.log10(i) for i in I.Section1.Bulk.CoTotal + C.Section1.Bulk.CoTotal + O.Section1.Bulk.CoTotal + SG.Section1.Bulk.CoTotal]
SolutionOxideCoTotal = [np.log10(i) for i in I.Section1.SolutionOxide.CoTotal + C.Section1.SolutionOxide.CoTotal + \
                        O.Section1.SolutionOxide.CoTotal + SG.Section1.SolutionOxide.CoTotal]
SolutionOxideCoSat = [np.log10(i) for i in I.Section1.SolutionOxide.CoSatFerrite + C.Section1.SolutionOxide.CoSatFerrite + \
                        O.Section1.SolutionOxide.CoSatFerrite + SG.Section1.SolutionOxide.CoSatFerrite]


BulkCrTotal = [np.log10(i) for i in I.Section1.Bulk.CrTotal + C.Section1.Bulk.CrTotal + O.Section1.Bulk.CrTotal + SG.Section1.Bulk.CrTotal]

SolutionOxideCrSat = [np.log10(i) for i in I.Section1.SolutionOxide.CrSat + C.Section1.SolutionOxide.CrSat + O.Section1.SolutionOxide.CrSat + \
                      SG.Section1.SolutionOxide.CrSat]


InnerOxide = []
InnerOxideThicknesses = [I.Section1.InnerOxThickness, C.Section1.InnerOxThickness, O.Section1.InnerOxThickness, SG.Section1.InnerOxThickness]
for Thickness, Sect in zip (InnerOxideThicknesses, ld.Sections):
    z = ld.UnitConverter(Sect, "Grams per Cm Squared", "Grams per M Squared", None, None, Thickness, None, None, None)
    InnerOxide.append(z)
InnerOxideThickness = InnerOxide[0] + InnerOxide[1] + InnerOxide[2] + InnerOxide[3]

OuterOxide = []
OuterOxideThicknesses = [I.Section1.OuterOxThickness, C.Section1.OuterOxThickness, O.Section1.OuterOxThickness, SG.Section1.OuterOxThickness]
for Thickness, Sect in zip (OuterOxideThicknesses, ld.Sections):
    q = ld.UnitConverter(Sect, "Grams per Cm Squared", "Grams per M Squared", None, None, Thickness, None, None, None)
    OuterOxide.append(q)
OuterOxideThickness = OuterOxide[0] + OuterOxide[1] + OuterOxide[2] + OuterOxide[3]

Cobalt = []
CoThicknesses = [I.Section1.CoThickness, C.Section1.CoThickness, O.Section1.CoThickness, SG.Section1.CoThickness]
for Thickness, Sect in zip (CoThicknesses, ld.Sections):
    c = ld.UnitConverter(Sect, "Grams per Cm Squared", "Grams per M Squared", None, None, Thickness, None, None, None)
    Cobalt.append(c)
CobaltThickness = Cobalt[0] + Cobalt[1] + Cobalt[2] + Cobalt[3]

Nickel = []
NiThicknesses = [I.Section1.NiThickness, C.Section1.NiThickness, O.Section1.NiThickness, SG.Section1.NiThickness]
for Thickness, Sect in zip (NiThicknesses, ld.Sections):
    n = ld.UnitConverter(Sect, "Grams per Cm Squared", "Grams per M Squared", None, None, Thickness, None, None, None)
    Nickel.append(n)
NickelThickness = Nickel[0] + Nickel[1] + Nickel[2] + Nickel[3]


# LogCorrosion= [np.log10(i) for i in Corrosion]


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

# fig3, ax1 = plt.subplots()
# #ax1 = fig2.add_subplot(221)
# ax1.plot(ActivityTime, Co60SGActivity, linestyle ="-", marker='o', color='0.50', label='Co-60')
# ax1.plot(ActivityTime, Co58SGActivity, linestyle ="-", marker='o', color='b', label='Co-58')
# ax1.plot(ActivityTime, Fe59SGActivity, linestyle ="-", marker='o', color='r', label='Fe-59')
# ax1.plot(ActivityTime, Fe55SGActivity, linestyle ="-", marker='o', color='g', label='Fe-55')
# ax1.plot(ActivityTime, Mn54SGActivity, linestyle ="-", marker='o', color='m', label='Mn-54')
# ax1.plot(ActivityTime, Cr51SGActivity, linestyle ="-", marker='o', color='y', label='Cr-51')
# ax1.plot(ActivityTime, Ni63SGActivity, linestyle ="-", marker='o', color='k', label='Ni-63')
#  
# #ax1.axis([51,69, 0, 30])
# ax1.set_xlabel('EFPY')
# ax1.set_ylabel('Surface Activity (${mCi/m^2}$)')
#             
# plt.tight_layout()
# plt.show()


