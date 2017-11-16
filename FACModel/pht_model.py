import constants as nc
import lepreau_data as ld
import composition as c
import rk_4
#import activities as a
import electrochemistry as e
import iteration as it
import sg_heattransfer as SGHX

# Initial Temperatures and T-dependent parametrs in SG zones
SecondarySidePressure = 4.593

# for Zone in ld.SGZones:
#     Zone.PrimaryBulkTemperature = SGHX.temperature_profile(
#         Zone, Zone.InnerOxThickness, Zone.OuterOxThickness,
#         SGHX.MassFlow_h.magnitude - 0.03 * SGHX.MassFlow_h.magnitude, SecondarySidePressure, 1983
#         )
#     Zone.DensityH2O = [
#         ld.Density("water", "PHT", i, SecondarySidePressure) for i in Zone.PrimaryBulkTemperature
#         ]
#     Zone.ViscosityH2O = [ld.Viscosity("water", "PHT", i) for i in Zone.PrimaryBulkTemperature]
#     Zone.NernstConstant = [x * (2.303 * nc.R / (2 * nc.F)) for x in Zone.PrimaryBulkTemperature]

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
        ] = e.ElectrochemicalAdjustment(
        Section, Section.SolutionOxide.EqmPotentialFe3O4, Section.SolutionOxide.MixedPotential,
        Section.MetalOxide.MixedPotential, Section.SolutionOxide.FeTotal, Section.SolutionOxide.FeSatFe3O4,
        Section.Bulk.FeSatFe3O4, Section.SolutionOxide.ConcentrationH
        )


class PHT_FAC():
    def __init__(self, ActiveSection, OutgoingSection, RealTimeHeatTransfer, j):  # j = overall time step
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
#                 self.Section1.BigParticulate[i] = rk_4.particulate(
#                     self.Section1, self.Section1.BigParticulate[0], self.Section1.Diameter[i],
#                     self.Section1.DensityH2O[i], self.Section1.Velocity[i], self.Section1.Distance[i]
#                     )
#                 self.Section1.SmallParticulate[i] = rk_4.particulate(
#                     Section, self.Section1.SmallParticulate[0], self.Section1.Diameter[i], self.Section1.DensityH2O[i],
#                     self.Section1.Velocity[i], self.Section1.Distance[i]
#                     )
                                
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

     
# InnerOxide = []
# InnerOxideThicknesses = [I.Section1.InnerOxThickness, C.Section1.InnerOxThickness, O.Section1.InnerOxThickness, SG.Section1.InnerOxThickness]
# for Thickness, Sect in zip (InnerOxideThicknesses, ld.Sections):
#     z = ld.UnitConverter(Sect, "Grams per Cm Squared", "Grams per M Squared", None, None, Thickness, None, None, None)
#     InnerOxide.append(z)
# InnerOxideThickness = InnerOxide[0] + InnerOxide[1] + InnerOxide[2] + InnerOxide[3]
# 
# OuterOxide = []
# OuterOxideThicknesses = [I.Section1.OuterOxThickness, C.Section1.OuterOxThickness, O.Section1.OuterOxThickness, SG.Section1.OuterOxThickness]
# for Thickness, Sect in zip (OuterOxideThicknesses, ld.Sections):
#     q = ld.UnitConverter(Sect, "Grams per Cm Squared", "Grams per M Squared", None, None, Thickness, None, None, None)
#     OuterOxide.append(q)
# OuterOxideThickness = OuterOxide[0] + OuterOxide[1] + OuterOxide[2] + OuterOxide[3]
# 
# Cobalt = []
# CoThicknesses = [I.Section1.CoThickness, C.Section1.CoThickness, O.Section1.CoThickness, SG.Section1.CoThickness]
# for Thickness, Sect in zip (CoThicknesses, ld.Sections):
#     c = ld.UnitConverter(Sect, "Grams per Cm Squared", "Grams per M Squared", None, None, Thickness, None, None, None)
#     Cobalt.append(c)
# CobaltThickness = Cobalt[0] + Cobalt[1] + Cobalt[2] + Cobalt[3]
# 
# Nickel = []
# NiThicknesses = [I.Section1.NiThickness, C.Section1.NiThickness, O.Section1.NiThickness, SG.Section1.NiThickness]
# for Thickness, Sect in zip (NiThicknesses, ld.Sections):
#     n = ld.UnitConverter(Sect, "Grams per Cm Squared", "Grams per M Squared", None, None, Thickness, None, None, None)
#     Nickel.append(n)
# NickelThickness = Nickel[0] + Nickel[1] + Nickel[2] + Nickel[3]
# 

# 
# fig2, ax1 = plt.subplots()
# # ax1 = fig2.add_subplot(221)
# ax1.plot(TotalLoopDistance, InnerOxideThickness, linestyle=None, marker='o', color='0.50', label='Inner Oxide')
# ax1.plot(TotalLoopDistance, OuterOxideThickness, linestyle=None, marker='o', color='k', label='Outer Oxide')
# # ax1.axis([51,69, 0, 30])
# ax1.set_xlabel('Distance (m)')
# ax1.set_ylabel('Oxide Layer Loadings (${g/m^2}$)')
#            
# ax2 = ax1.twinx()
# ax2.plot(TotalLoopDistance, NickelThickness, linestyle=None, marker='o', color='c', label='Nickel')
# ax2.plot(TotalLoopDistance, CobaltThickness, linestyle=None, marker='o', color='m', label='Cobalt')
# ax2.set_ylabel('Ni, Co, and Cr Loadings (${g/m^2}$)', rotation=270, labelpad=20)
# lines, labels = ax1.get_legend_handles_labels()
# lines2, labels2 = ax2.get_legend_handles_labels()
# ax2.legend(lines + lines2, labels + labels2, loc=0)        
# # plt.axis([51, 69, 0, 30])
# plt.tight_layout()
# plt.show()

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


