import constants as nc
import lepreau_data as ld
import composition as c
import rk_4
#import activities as a
import electrochemistry as e
import iteration as it
import sg_heattransfer as SGHX

# Initial temperatures in SG zones
[SecondarySidePressure, RemainingPHTMassFlow, DividerPlateMassFlow] = SGHX.station_events(1983)

for Zone in ld.SGZones:
    Zone.PrimaryBulkTemperature = SGHX.temperature_profile(
        Zone, Zone.InnerOxThickness, Zone.OuterOxThickness, RemainingPHTMassFlow, SecondarySidePressure, 1983
        )

# Initial concentrations
for Section in ld.Sections:    
    # #Temperature-dependent parameters            
    Section.NernstConstant = [x * (2.303 * nc.R / (2 * nc.F)) for x in Section.PrimaryBulkTemperature]
    
    Section.DensityH2O = [ld.Density("water", "PHT", x, SecondarySidePressure) for x in Section.PrimaryBulkTemperature]
    Section.ViscosityH2O = [ld.Viscosity("water", "PHT", x) for x in Section.PrimaryBulkTemperature]
    
    Section.FractionFeInnerOxide = c.fraction_metal_inner_oxide(Section, "Fe")
    Section.FractionNiInnerOxide = c.fraction_metal_inner_oxide(Section, "Ni")
    
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
        Interface.FeTotal = c.iron_solubility(Section)
        Interface.FeSatFe3O4 = [1 * i for i in Interface.FeTotal]
        Interface.NiTotal = [i * 1 for i in Section.SolubilityNi]
        Interface.NiSatFerrite = [i * 1 for i in Section.SolubilityNi]
        Interface.NiSatMetallicNi = [i * 1 for i in Section.SolubilityNi]
    
        Interface.CoTotal = [i * 1 for i in Section.SolubilityCo]
        Interface.CoSatFerrite = [i * 1 for i in Section.SolubilityCo]
        Interface.CrTotal = [i * 1 for i in Section.SolubilityCr]
        Interface.CrSat = [i * 1 for i in Section.SolubilityCr]
        
        # Initial RIHT temperature 
        if Section == ld.Inlet:
            Section.PrimaryBulkTemperature = [SGHX.energy_balance(ld.SGZones[0].NodeNumber - 1, 0)]\
            * Section.NodeNumber
        
        if Section != ld.Outlet:
            if Interface == Section.SolutionOxide:
                Interface.FeTotal = [i * 0.8 for i in c.iron_solubility(Section)]
        
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
        
    [
        Section.KpFe3O4electrochem, Section.KdFe3O4electrochem, Section.SolutionOxide.FeSatFe3O4,
        Section.MetalOxide.ConcentrationH
        ] \
        = e.ElectrochemicalAdjustment(
        Section, Section.SolutionOxide.EqmPotentialFe3O4, Section.SolutionOxide.MixedPotential,
        Section.MetalOxide.MixedPotential, Section.SolutionOxide.FeTotal, Section.SolutionOxide.FeSatFe3O4,
        Section.Bulk.FeSatFe3O4, Section.SolutionOxide.ConcentrationH
        )


class PHT_FAC():
    def __init__(self, ActiveSection, OutgoingSection, RealTimeHeatTransfer, j):  # j = overall time step
        self.Section1 = ActiveSection
        self.Section2 = OutgoingSection
        
#         if Activation == "yes":    
#             BulkActivities = [
#                 self.Section1.Bulk.Co60, self.Section1.Bulk.Co58, self.Section1.Bulk.Fe59, self.Section1.Bulk.Fe55,
#                 self.Section1.Bulk.Mn54, self.Section1.Bulk.Cr51, self.Section1.Bulk.Ni63
#                 ]
#                
#             BulkActivities2 = [
#                 self.Section2.Bulk.Co60, self.Section2.Bulk.Co58, self.Section2.Bulk.Fe59, self.Section2.Bulk.Fe55,
#                 self.Section2.Bulk.Mn54, self.Section2.Bulk.Cr51, self.Section2.Bulk.Ni63
#                 ]
#                
# #             SurfaceActivities = [
# #                 self.Section1.Bulk.Co60, self.Section1.Bulk.Co58, self.Section1.Bulk.Fe59, self.Section1.Bulk.Fe55,
# #                 self.Section1.Bulk.Mn54, self.Section1.Bulk.Cr51, self.Section1.Bulk.Ni63
# #                 ]
#                
#             Tags = ["Co60", "Co58", "Fe59", "Fe55", "Mn54", "Cr51", "Ni63"]
#             
#             for x, y in zip (BulkActivities, Tags):
#                 x = a.BulkActivity(self.Section1, x[0], y, j)
         
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
                for x, y in zip (BulkConcentrations, SolutionOxideConcentrations):
                    x[i] = rk_4.Spatial(
                        y[i - 1], x[i - 1], ld.MassTransfer(self.Section1)[i - 1], self.Section1.Diameter[i - 1],
                        self.Section1.Velocity[i - 1], self.Section1.Length.magnitude[i - 1]
                        )

                # Exponential decay of bulk particulate at start of section as function of distance + removal due to 
                # deposition and erosion source
                self.Section1.BigParticulate[i] = rk_4.particulate(
                    self.Section1, self.Section1.BigParticulate[0], self.Section1.Diameter[i],
                    self.Section1.DensityH2O[i], self.Section1.Velocity[i], self.Section1.Distance[i]
                    )
#                 self.Section1.SmallParticulate[i] = rk_4.particulate(
#                     Section, self.Section1.SmallParticulate[0], self.Section1.Diameter[i], self.Section1.DensityH2O[i],
#                     self.Section1.Velocity[i], self.Section1.Distance[i]
#                     )

                    
                # Inlet header purification system
                if self.Section1 == ld.Inlet:      
                    if i == 3:
                        # for x,y in zip(BulkConcentrations, BulkActivities):
                        for x in BulkConcentrations: 
                            x[i] = 0.59 * x[i - 1]
                            # y[i] = 0.59*y[i-1]
                    
#                         self.Section1.SmallParticulate[i] = 0.59 * self.Section1.SmallParticulate[i - 1]
#                         self.Section1.BigParticulate[i] = 0  # 0.45 um filter removes everything over this size   
                    
#                     if i > 3:  # decay for the rest of the inlet section depends on purification system 
#                         self.Section1.BigParticulate[i] = rk_4.particulate(
#                             Section, self.Section1.BigParticulate[3], self.Section1.Diameter[i],
#                             self.Section1.DensityH2O[i], self.Section1.Velocity[i], self.Section1.Distance[i]
#                             )
#                         self.Section1.SmallParticulate[i] = rk_4.particulate(
#                             Section, self.Section1.SmallParticulate[3], self.Section1.Diameter[i],
#                             self.Section1.DensityH2O[i], self.Section1.Velocity[i], self.Section1.Distance[i]
#                             )
#                         for x,y in zip (BulkActivities, Tags):
#                             x = a.BulkActivity(self.Section1, x[3], y, j)
                                        
        end = self.Section1.NodeNumber - 1
        for x, y in zip(BulkConcentrations, BulkConcentrations2):
        # for x,y,z,q in zip(BulkConcentrations,BulkConcentrations2, BulkActivities, BulkActivities2):
            y[0] = x[end]
            # q[0] = z[end]
        
#         self.Section2.BigParticulate[0] = self.Section1.BigParticulate[end]
#         self.Section2.SmallParticulate[0] = self.Section1.SmallParticulate[end]       
        
        # #Stellite wear bulk input term for cobalt and chromium    
        if self.Section1 == ld.Core:
            self.Section2.Bulk.CoTotal[0] = self.Section1.Bulk.CoTotal[end] + nc.CobaltWear
            self.Section2.Bulk.CrTotal[0] = (self.Section1.Bulk.CrTotal[end]
                                             + nc.CobaltWear * (nc.FractionCr_Stellite / nc.FractionCo_Stellite))
   
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

        # SG heat transfer 
        if RealTimeHeatTransfer == "yes":
            if self.Section1 in ld.SGZones:   
                self.Section1.PrimaryBulkTemperature = SGHX.temperature_profile(
                    self.Section1, self.Section1.InnerOxThickness, self.Section1.OuterOxThickness
                    )
                self.Section1.Bulk.FeSatFe3O4 = c.iron_solubility(self.Section1) 
          
        
            # RIHT  
            if self.Section1 == ld.Inlet:
                self.Section1.PrimaryBulkTemperature = SGHX.energy_balance(21, j)
                self.Section1.Bulk.FeSatFe3O4 = c.iron_solubility(self.Section1)

        # Deposit thickness around PHTS
        # self.Section1.DepositThickness = a.Deposition(self.Section1, self.Section1.BigParticulate, self.Section1.SmallParticulate, j)
                        
        # #rk_4 oxide thickness calculation (no spalling)
        rk_4.oxidegrowth(self.Section1, Saturations, BulkConcentrations, ElementTracking = "no")
        
        # Spalling    
        self.Section1.ElapsedTime, self.Section1.SpallTime = rk_4.Spall(
            self.Section1, j, self.Section1.ElapsedTime, self.Section1.SpallTime, ElementTracking= "no"
            )


