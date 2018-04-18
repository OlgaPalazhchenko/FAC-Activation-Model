import thermochemistry_and_constants as nc
import lepreau_data as ld
import composition as c
import rk_4
import activities as a
import electrochemistry as e
import iteration as it
import sg_heattransfer as SGHX

from operator import itemgetter

def initial_chemistry(Loop):
    
    [SecondarySidePressure, RemainingPHTMassFlow, DividerPlateMassFlow] = SGHX.station_events(
        SGHX.YearStartup, x_pht= 0.002
        )
    # initial temperatures in steam generator(s)
    
    if Loop == "full":
        # if full loop selected, both steam generators iterated through for initial temperatures (87 bundles each)
        Sections = ld.FullLoop
#         for Zone in ld.SteamGenerator_2:
#             Zone.PrimaryBulkTemperature = SGHX.temperature_profile(
#             Zone, Zone.InnerOxThickness, Zone.OuterOxThickness, RemainingPHTMassFlow, SecondarySidePressure, 1983
#             )       
    elif Loop == "half":
        Sections = ld.HalfLoop
    
    elif Loop == "single tube":
        Sections = ld.SingleTubeLoop
    else:
        None
   
    if Loop == "full" or Loop == "half":
        for Zone in ld.SteamGenerator:
            Zone.PrimaryBulkTemperature = SGHX.temperature_profile(
            Zone, Zone.InnerOxThickness, Zone.OuterOxThickness, RemainingPHTMassFlow, SecondarySidePressure,
            x_pht=0.002, calendar_year=SGHX.YearStartup
            )
    else:
        ld.SteamGenerator[57].PrimaryBulkTemperature = SGHX.temperature_profile(
            ld.SteamGenerator[57], ld.SteamGenerator[57].InnerOxThickness, ld.SteamGenerator[57].OuterOxThickness,
            RemainingPHTMassFlow, SecondarySidePressure, x_pht=0.2, calendar_year=SGHX.YearStartup
            )
        
    # initial concentrations
    for Section in Sections:
        # temperature-dependent parameters            
        Section.NernstConstant = [x * (2.303 * nc.R / (2 * nc.F)) for x in Section.PrimaryBulkTemperature]
        
        Section.DensityH2O = [
            nc.density("PHT", x, SecondarySidePressure) for x in Section.PrimaryBulkTemperature
            ]
        Section.ViscosityH2O = [
            nc.viscosity("PHT", x, SecondarySidePressure) for x in Section.PrimaryBulkTemperature
            ]
        
        Section.FractionFeInnerOxide = c.fraction_metal_inner_oxide(Section, "Fe")
        Section.FractionNiInnerOxide = c.fraction_metal_inner_oxide(Section, "Ni")
        
        Interfaces = [Section.MetalOxide, Section.SolutionOxide, Section.Bulk]
        for Interface in Interfaces:
            # makes into array with appropriate # of nodes for that PHTS section
            Interface.ConcentrationH2 = [nc.H2 * nc.H2Density / nc.H2MolarMass] * Section.NodeNumber
            Interface.ConcentrationH = c.bulkpH_calculator(Section)  # from system pH calculation
            
            # concentration/Saturation Input [mol/kg]
            Interface.FeTotal = c.iron_solubility_SB(Section)
            Interface.FeSatFe3O4 = [1 * i for i in Interface.FeTotal]
            Interface.NiTotal = [i * 1 for i in Section.SolubilityNi]
            Interface.NiSatFerrite = [i * 1 for i in Section.SolubilityNi]
            Interface.NiSatMetallicNi = [i * 1 for i in Section.SolubilityNi]
        
            Interface.CoTotal = [i * 1 for i in Section.SolubilityCo]
            Interface.CoSatFerrite = [i * 1 for i in Section.SolubilityCo]
            Interface.CrTotal = [i * 1 for i in Section.SolubilityCr]
            Interface.CrSat = [i * 1 for i in Section.SolubilityCr]
            
            [
                Interface.Co60, Interface.Co58, Interface.Fe55, Interface.Fe59, Interface.Mn54, Interface.Cr51,
                Interface.Ni63
                ] = [[0] * Section.NodeNumber] * 7
             
            # initial RIHT temperature 
            if Section in ld.InletSections:
                Section.PrimaryBulkTemperature = (
                    [SGHX.energy_balance(21, ld.SteamGenerator[0].InnerOxThickness,
                                         ld.SteamGenerator[0].OuterOxThickness, ld.SteamGenerator[0].InnerOxThickness,
                                         ld.SteamGenerator[0].OuterOxThickness, x_pht=0.002, j=0,
                                         SGFastMode=SGHX.SGFastMode)] * Section.NodeNumber
                                                )
            
            if Section in ld.FuelSections and Interface == Section.MetalOxide:
                Interface.FeTotal = [0] * Section.NodeNumber
            
            if Section in ld.OutletSections and Interface == Section.MetalOxide:
                # from Cook's thesis - experimental corrosion rate measurements and calcs
                Interface.FeTotal = [0.00000026] * Section.NodeNumber 
            
            if Section not in ld.SteamGenerator and Section not in ld.SteamGenerator_2:
                Interface.NiTotal = [0] * Section.NodeNumber
            
        Section.SolutionOxide.MixedPotential, Section.SolutionOxide.EqmPotentialFe3O4 = e.ECP(Section)
        
        if Section in ld.FuelSections:
            Section.CorrRate == [0] * Section.NodeNumber
        else:
            Section.CorrRate, Section.MetalOxide.MixedPotential = it.FAC_solver(Section, ConstantRate="yes")
            
        
        Section.KpFe3O4electrochem = [nc.KpFe3O4] * Section.NodeNumber
        Section.KdFe3O4electrochem = [nc.KdFe3O4] * Section.NodeNumber   
#         [
#             Section.KpFe3O4electrochem, Section.KdFe3O4electrochem, Section.SolutionOxide.FeSatFe3O4,
#             Section.MetalOxide.ConcentrationH
#             ] \
#             = e.ElectrochemicalAdjustment(
#             Section, Section.SolutionOxide.EqmPotentialFe3O4, Section.SolutionOxide.MixedPotential,
#             Section.MetalOxide.MixedPotential, Section.SolutionOxide.FeTotal, Section.SolutionOxide.FeSatFe3O4,
#             Section.Bulk.FeSatFe3O4, Section.SolutionOxide.ConcentrationH
#             )

Tags = ["Co60", "Co58", "Fe59", "Fe55", "Mn54", "Cr51", "Ni63"]


class PHT_FAC():
    def __init__(self, ActiveSection, OutgoingSection, ElementTracking, Activation, ConstantRate, j):  
        # j = time step
        self.Section1 = ActiveSection
        self.Section2 = OutgoingSection
        
        if Activation == "yes":    
            BulkActivities = [
                self.Section1.Bulk.Co60, self.Section1.Bulk.Co58, self.Section1.Bulk.Fe59, self.Section1.Bulk.Fe55,
                self.Section1.Bulk.Mn54, self.Section1.Bulk.Cr51, self.Section1.Bulk.Ni63
                ]
    
            BulkActivities2 = [
                self.Section2.Bulk.Co60, self.Section2.Bulk.Co58, self.Section2.Bulk.Fe59, self.Section2.Bulk.Fe55,
                self.Section2.Bulk.Mn54, self.Section2.Bulk.Cr51, self.Section2.Bulk.Ni63
                ]
             
            SurfaceActivities = [
                self.Section1.MetalOxide.Co60, self.Section1.MetalOxide.Co58, self.Section1.MetalOxide.Fe59,
                self.Section1.MetalOxide.Fe55, self.Section1.MetalOxide.Mn54, self.Section1.MetalOxide.Cr51,
                self.Section1.MetalOxide.Ni63
                ]
         
            # Deposit thickness around PHTS only calculated if activity transport is being tracked
            self.Section1.DepositThickness = a.deposition(self.Section1, j)
            
            # Exponential decay of bulk particulate at start of section as function of distance + removal due to 
            # deposition and addition from erosion
            self.Section1.BigParticulate = a.particulate(self.Section1, self.Section1.BigParticulate[0])
            self.Section1.SmallParticulate = a.particulate(self.Section1, self.Section1.SmallParticulate[0])
            
            # self.Section1.NodeNumber - 1 = last node of section
            self.Section2.BigParticulate[0] = self.Section1.BigParticulate[self.Section1.NodeNumber - 1]
            self.Section2.SmallParticulate[0] = self.Section1.SmallParticulate[self.Section1.NodeNumber - 1]
            
        
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
            # solves for bulk volumetric activities and connects activity concentrations b/w PHT sections
            if Activation == "yes":
                self.Section1.Bulk.Co60[i] = a.bulk_activity(self.Section1, self.Section1.Bulk.Co60[0], "Co60", j, i)
                self.Section1.Bulk.Fe59[i] = a.bulk_activity(self.Section1, self.Section1.Bulk.Fe59[0], "Fe59", j, i)
#                 for x, y in zip (BulkActivities, Tags):
#                     x[i] = a.bulk_activity(self.Section1, x[0], y, j, i)
                    
                for x, y, z in zip (BulkActivities, SurfaceActivities, Tags):
                    y[i] = a.surface_activity(self.Section1, x[i], j, i, z)
                    
            for x, y in zip (BulkConcentrations, SolutionOxideConcentrations):
                if i > 0:
                    x[i] = rk_4.spatial(
                        self.Section1, y[i - 1], x[i - 1], ld.MassTransfer(self.Section1)[i], self.Section1.Diameter[i],
                        self.Section1.Velocity[i], self.Section1.Length.magnitude[i], i
                        )
                
            # Inlet header purification systemc comes off of only one inlet feeder header in a full figure-of-8 loop
            if self.Section1 == ld.InletFeeder and i == 2:
                purificationfactor = 1900 / (1900 + 12)     
#                 if j > 168: 
#                     for x in BulkConcentrations:
#                         x[i] =  purificationfactor * x[i - 1]
                if Activation == "yes":
                    for y in BulkActivities:
                        y[i] = purificationfactor * y[i - 1]
            
            elif self.Section1 in ld.InletSections and i > 2:
                if Activation == "yes":
                    self.Section1.BigParticulate = a.particulate(self.Section1, self.Section1.BigParticulate[3])
                    self.Section1.SmallParticulate = a.particulate(self.Section1, self.Section1.SmallParticulate[3])
                    
                    for x,y in zip (BulkActivities, Tags):
                        x[i] = a.bulk_activity(self.Section1, x[3], y, j, i)
            else:
                None
                
       
        # Connects output of one PHT section to input of subsequent section 
        for x, y in zip(BulkConcentrations, BulkConcentrations2):
            y[0] = x[self.Section1.NodeNumber - 1]
        
        if Activation == "yes":
            for x, y in zip(BulkActivities, BulkActivities2):
                y[0] = x[self.Section1.NodeNumber - 1]
        
  
        # Stellite wear bulk input term for cobalt and chromium    
        if self.Section1 in ld.FuelSections:
            self.Section2.Bulk.CoTotal[0] = self.Section1.Bulk.CoTotal[self.Section1.NodeNumber - 1] + nc.CobaltWear
            self.Section2.Bulk.CrTotal[0] = (self.Section1.Bulk.CrTotal[self.Section1.NodeNumber - 1]
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


        # RK4 oxide thickness calculation (no spalling)
        rk_4.oxide_layers(
            self.Section1, ConstantRate, Saturations, BulkConcentrations, ElementTracking, j, SGHX.SGFastMode)
                
        # Spalling    
        self.Section1.ElapsedTime, self.Section1.SpallTime = rk_4.spall(
            self.Section1, j, self.Section1.ElapsedTime, self.Section1.SpallTime, ElementTracking
            )
