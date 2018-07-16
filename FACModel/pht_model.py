import thermochemistry_and_constants as nc
import lepreau_data as ld
import composition as c
import rk_4
import activities as a
import electrochemistry as e
import iteration as it
import sg_heattransfer as SGHX
import numpy as np
from datetime import date, timedelta

# from operator import itemgetter

Loop = "half"
ConstantRate = "yes" # preset corrosion rate instead of calculated from FAC model
Purification = "yes"
Activation = "no"
OutputLogging = "yes"
ElementTracking = "no"

if Loop == "half":
    Sections = ld.HalfLoop
elif Loop == "full":
    Sections = ld.FullLoop
else:
    None

def initial_conditions():
    
    # initial temperatures in steam generator(s)
    RIHT1 = SGHX.energy_balance(ld.SteamGenerator, x_pht = 0.01, j = 0, SGFastMode = SGHX.SGFastMode)
    ld.InletFeeder.PrimaryBulkTemperature = [RIHT1] * ld.InletFeeder.NodeNumber
        
    RIHT2 = SGHX.energy_balance(ld.SteamGenerator_2, x_pht = 0.01, j = 0, SGFastMode = SGHX.SGFastMode)
    ld.InletFeeder_2.PrimaryBulkTemperature = [RIHT2] * ld.InletFeeder_2.NodeNumber
  
        
    for Section in Sections:
        Section.KpFe3O4electrochem = [nc.KpFe3O4] * Section.NodeNumber
        Section.KdFe3O4electrochem = [nc.KdFe3O4] * Section.NodeNumber   
               
        # temperature-dependent parameters            
        Section.NernstConstant = [x * (2.303 * nc.R / (2 * nc.F)) for x in Section.PrimaryBulkTemperature]
        
        Section.DensityH2O = [
            nc.densityH2O_liquid(x, nc.PrimarySidePressure) for x in Section.PrimaryBulkTemperature
            ]
        Section.ViscosityH2O = [
            nc.viscosityH2O_liquid(x, nc.PrimarySidePressure) for x in Section.PrimaryBulkTemperature
            ]
        
        Section.FractionFeInnerOxide = c.fraction_metal_inner_oxide(Section, "Fe")
        Section.FractionNiInnerOxide = c.fraction_metal_inner_oxide(Section, "Ni")
        
        # initial concentrations
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
            Section.CorrRate, Section.MetalOxide.MixedPotential = it.FAC_solver(Section, "yes", j=0)
        
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


class PHTS():
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
                if j >= 168: 
                    for x in BulkConcentrations:
                        x[i] =  purificationfactor * x[i - 1]
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


FACRate_OutletFeeder = []
RIHT_average = []
RIHT_InletFeeder1 = []
RIHT_InletFeeder2 = []
YM = []
OutletTemperature_Bundle_1 = [] 
OutletTemperature_Bundle_2 = []
pht_SteamFraction = []

def output_time_logging(FACRate, RIHT_avg, RIHT1, RIHT2, x, Temperature1, Temperature2):
        
    FACRate_OutletFeeder.append(FACRate)
    RIHT_InletFeeder1.append(RIHT1)
    RIHT_InletFeeder2.append(RIHT2)
    RIHT_average.append(RIHT_avg)
    pht_SteamFraction.append(x)
#     OutletTemperature_Bundle_1.append(Temperature1)
#     OutletTemperature_Bundle_2.append(Temperature2)
    
    return (
        FACRate_OutletFeeder, RIHT_InletFeeder1, RIHT_InletFeeder2, RIHT_average, pht_SteamFraction,
        OutletTemperature_Bundle_1, OutletTemperature_Bundle_2
        ) 


def sg_heat_transfer(Outlet, Inlet, SelectedTubes, j):
    Tubes = []
    # Set input concentrations for all SG zones to be same as output of outlet feeder
    BulkOutletActivityOutput = [
        Outlet.Bulk.Co60, Outlet.Bulk.Co58, Outlet.Bulk.Fe59, Outlet.Bulk.Fe55, Outlet.Bulk.Cr51, Outlet.Bulk.Mn54,
        Outlet.Bulk.Ni63
        ]
    
    BulkOutletOutput = [Outlet.Bulk.FeTotal, Outlet.Bulk.NiTotal, Outlet.Bulk.CoTotal, Outlet.Bulk.CrTotal]
    
    for Tube in SelectedTubes:
        BulkSGInput = [Tube.Bulk.FeTotal, Tube.Bulk.NiTotal, Tube.Bulk.CoTotal, Tube.Bulk.CrTotal]
        BulkSGActivityInput = [
        Tube.Bulk.Co60, Tube.Bulk.Co58, Tube.Bulk.Fe59, Tube.Bulk.Fe55, Tube.Bulk.Cr51, Tube.Bulk.Mn54, Tube.Bulk.Ni63
        ]
        for x, y, z, q in zip(BulkSGInput, BulkOutletOutput, BulkSGActivityInput, BulkOutletActivityOutput):
            x[0] = y[Outlet.NodeNumber - 1]
            z[0] = q[Outlet.NodeNumber - 1]
        
        w = PHTS(Tube, Inlet, ElementTracking, Activation, ConstantRate, j)   
        Tubes.append(w)
    return Tubes


# just a tube number generator (number of the tube that has closest u-bend arc length to the avg. 1.52 m length)
Default_Tube = SGHX.closest_ubend(1.52 * 100)
SimulationYears = 25 # years
SimulationHours = SimulationYears * 876 # 851


import time
start_time = time.time()

# load initial chemistry for full/half loop
initial_conditions()

for j in range(SimulationHours):
   
    if Loop == "full":
        InletInput = ld.InletFeeder_2
    # for 1/2 of single figure-of-eight loop, SG flow returned to same inlet header as that at start of loop
    # for full loop, 
    else:
        InletInput = ld.InletFeeder
    # half loop iterates through 4 sections (inlet feeder 1, fuel channel 1, outlet feeder 2, sg 2, back to inlet 1)        
    
    InletFeeder_1_Loop1 = PHTS(ld.InletFeeder, ld.FuelChannel, ElementTracking, Activation, ConstantRate, j)
    
    FuelChannel_1_Loop1 = PHTS(ld.FuelChannel, ld.OutletFeeder_2, ElementTracking, Activation, ConstantRate, j)
    
    OutletFeeder_2_Loop1 = PHTS(
    ld.OutletFeeder_2, ld.SteamGenerator_2[Default_Tube], ElementTracking, Activation, ConstantRate, j
    )
    
    SteamGeneratorTube_2_Loop1 = PHTS(
        ld.SteamGenerator_2[Default_Tube], InletInput, ElementTracking, Activation, ConstantRate, j
        )
    
    SelectedTubes, SelectedTubeNumbers = SGHX.tube_picker(SGHX.Method, ld.SteamGenerator_2)
    SteamGeneratorTubes_2 = sg_heat_transfer(ld.OutletFeeder_2, InletInput, SelectedTubes, j)
    
    if Loop == "full":
    
        InletFeeder_2_Loop1 = PHTS(ld.InletFeeder_2, ld.FuelChannel_2, ElementTracking, Activation, ConstantRate, j)
        
        FuelChannel_2_Loop1 = PHTS(ld.FuelChannel_2, ld.OutletFeeder, ElementTracking, Activation, ConstantRate, j)
        
        OutletFeeder_1_Loop1 = PHTS(
        ld.OutletFeeder, ld.SteamGenerator[Default_Tube], ElementTracking, Activation, ConstantRate, j
        )
        
        SteamGeneratorTube_1_Loop1 = PHTS(
            ld.SteamGenerator[Default_Tube], ld.InletFeeder, ElementTracking, Activation, ConstantRate, j
            )
        
        SelectedTubes, SelectedTubeNumbers = SGHX.tube_picker(SGHX.Method, ld.SteamGenerator)
        SteamGeneratorTubes_1 = sg_heat_transfer(ld.OutletFeeder, ld.InletFeeder, SelectedTubes, j)
    
    else:
        None

    # parameters tracked/updated with time
    if j % (73) == 0:  # 73 h * 10 = 12 x a year  
        
        if j ==0:
            x_pht = 0.01 # PHT steam fraction assumed for "clean" boiler

        if Loop == "full":
            RIHT_1 = SGHX.energy_balance(ld.SteamGenerator, x_pht, j, SGHX.SGFastMode) - 273.15
            RIHT_2 = SGHX.energy_balance(ld.SteamGenerator_2, x_pht, j, SGHX.SGFastMode) - 273.15
            ld.InletFeeder_2.PrimaryBulkTemperature = [RIHT_2 + 273.15] * ld.InletFeeder_2.NodeNumber
        
        else:
            #SG 2 is in the half and full loop configurations (default steam generator)
            RIHT_2 = SGHX.energy_balance(ld.SteamGenerator_2, x_pht, j, SGHX.SGFastMode) - 273.15
            RIHT_1 = RIHT_2 * 1  # output logging function needs value for both RIHT's and inlet header 1 needs input 

        # in full loop mode, this is calculated based on SG 1 output
        # in half loop mode, this is based on SG 2 output (RIHT_1 set equal to RIHT_2)
        ld.InletFeeder.PrimaryBulkTemperature  = [RIHT_1 + 273.15] * ld.InletFeeder.NodeNumber
        
        # in half loop mode, these are equal, so avg = RIHT_1 = RIHT_2
        T_RIH_average = (RIHT_1 + RIHT_2) / 2
        x_pht = SGHX.pht_steam_quality(T_RIH_average + 273.15, j)
        
        # core and outlet temperatures currently not being updated, but all sections called for continuity
        for Section in Sections:
            # update solubility based on new temperatures in steam generators and inlet heades/feeders
            Section.Bulk.FeSatFe3O4 = c.iron_solubility_SB(Section)        

        #final node temperature of first 2 of selected bundles (first two bundles from all 87 if not run in "fast mode")
        # currently tracking only one steam generator...second set of Temp1 and Temp2 can be added for SG1 if needed
        Temperature1 = (
            ld.SteamGenerator_2[SelectedTubeNumbers[0]].PrimaryBulkTemperature[21] - 273.15
                       )

        if len(SelectedTubes) > 1: # provided that the selected bundle array is not just one element
            Temperature2 = (
                ld.SteamGenerator_2[SelectedTubeNumbers[1]].PrimaryBulkTemperature[21] - 273.15
                           )
        else:
            Temperature2 = None
            
        # optional preview of RIHT and primary-side steam quality
        start = SGHX.YearStartup
        delta = timedelta(hours = j * nc.TIME_STEP)
        CalendarYear = start + delta
        
        Year_Month = (CalendarYear.year, CalendarYear.month)
        
        print (Year_Month, x_pht, RIHT_1)#, RIHT_2)
#         SelectedTubes = SGHX.tube_picker(SGHX.Method, ld.SteamGenerator_2)[1]
#          
#         print ('cleaned', ld.SteamGenerator_2[SelectedTubes[0]].OuterFe3O4Thickness)
#         print ('normal', ld.SteamGenerator_2[Default_Tube].OuterFe3O4Thickness)
        output = output_time_logging(
            OutletFeeder_2_Loop1.Section1.CorrRate, T_RIH_average, RIHT_1, RIHT_2, x_pht, Temperature1,
            Temperature2
            )
            
    else:
        None
    
    
end_time = time.time()
delta_time = end_time - start_time

 
hours = delta_time // 3600
temp = delta_time - 3600 * hours
minutes = delta_time // 60
seconds = delta_time - 60 * minutes
print('%d:%d:%d' % (hours, minutes, seconds))
