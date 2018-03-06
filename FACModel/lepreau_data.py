import numpy as np 
import thermochemistry_and_constants as nc
import csv


def UnitConverter(Section, UnitInput, UnitOutput, Concentration, Rate, Oxide, OxideDensity, MolarMass, Temperature):
    # Converts temperature from Celsius to Kelvin
    if UnitInput == "Celsius" and UnitOutput == "Kelvin":
        return [i + 273.15 for i in Temperature]

    # Converts surface area from /cm^2 to /m^2
    elif UnitInput == "Grams per Cm Squared" and UnitOutput == "Grams per M Squared":
        return [i * (100 ** 2) for i in Oxide]

    # Converts coolant activity from Bq/cm^3 to microCi/m^3
    elif UnitInput == "Bq" and UnitOutput == "Ci":
        # return Concentration/(3.7*10**10)
        return [i / (3.7 * 10 ** 10) for i in Concentration]
    
    # Converts concentrations between mol/kg to g/cm^3 and vice versa
    # [mol_solute/kg_coolant] *[1 kg/1000 g coolant]* [g/mol solute] = [g_solute]/[g_coolant]* [g/cm^3 coolant]
    # = [g_solute/cm^3 coolant]
    elif UnitInput == "Mol per Kg" and UnitOutput == "Grams per Cm Cubed":
        return [(x / 1000) * (MolarMass * y) for x, y in zip(Concentration, Section.DensityH2O)]

    elif UnitInput == "Grams per Cm Cubed" and UnitOutput == "Mol per Kg":
    # [g_solute/cm^3 coolant]/([g/cm^3 coolant]*[g/mol solute]) = [mol_solute/g coolant]*[1000 g/1 kg coolant]
    # = [mol solute/kg coolant]
        return [x * 1000 / (MolarMass * y) for x, y in zip(Concentration, Section.DensityH2O)]

    # Converts corrosion rates from [g/cm^2 s] to micrometers per annum [um/yr]
    elif UnitInput == "Corrosion Rate Grams" and UnitOutput == "Corrosion Rate Micrometers":
    # ([g/cm^2 s] / [g/cm^3] *[10000 um/cm]) = ([cm/s] * [3600 s/h] * [24 h/d] * [365 d/yr]) = [um/yr]
        if Section in InletSections or Section in OutletSections:
            return [i * ((1 / nc.FeDensity) * 3600 * 24 * 365 * 10000) for i in Rate]
        
        elif Section in SteamGenerator or Section in SteamGenerator_2:
            return [i * ((1 / nc.Alloy800Density) * 3600 * 24 * 365 * 10000) for i in Rate]
        
        else:
            return None
    
    # Converts from oxide mass per surface area to oxide thickness [g/cm^2] to [um]
    # assuming uniform distribution on surface   
    elif UnitInput == "Oxide Thickness Grams" and UnitOutput == "Oxide Thickness Micrometers":
        return [x * 10000 / y for x, y in zip(Oxide, OxideDensity)]

    else:
        None

SizingParameters = open('SizingParameters.txt', 'r')      
SizingParametersReader = list(csv.reader(SizingParameters, delimiter=','))  


class Interface():
    def __init__(self):
        self.ConcentrationH2 = None
        self.ConcentrationH = None
        self.FeSatFe3O4 = None
        self.FeTotal = None
        self.ConcentrationFe2 = None
        # self.ConcentrationFeOH = None
        self.ConcentrationFeOH2 = None
        # self.ConcentrationFeOH3 = None

        self.NiSatFerrite = None
        self.NiSatMetallicNi = None
        self.NiTotal = None
        self.ConcentrationNi2 = None
        self.ConcentrationNiOH = None
        self.CrSat = None
        self.CoSatFerrite = None
        self.CrTotal = None
        self.CoTotal = None

        self.Co60 = None
        self.Co58 = None
        self.Fe55 = None
        self.Fe59 = None
        self.Mn54 = None
        self.Cr51 = None
        self.Ni63 = None

        self.MixedPotential = None
        self.EqmPotentialFe3O4 = None


class Section():  # Defining each primary heat transport section as a class 
    def __init__(self, j, RowStart, RowEnd): 
        self.RowStart = RowStart
        self.RowEnd = RowEnd
        self.j = j
        self.NodeNumber = RowEnd

        self.MetalOxide = Interface()
        self.SolutionOxide = Interface()
        self.Bulk = Interface()

        # General section parameter input (interface specification not necessary)
        self.Diameter = None
        self.OuterDiameter = None
        self.TubeThickness = None
        self.Velocity = None
        self.Length = nc.SGParameters() 
        self.Distance = None 
        self.PrimaryBulkTemperature = None
        # self.PrimaryWallTemperature = None
        # self.SecondaryWallTemperature = None
        self.SecondaryBulkTemperature = None
        self.NernstConstant = None

        self.DensityH2O = None
        self.ViscosityH2O = None
        self.SolubilityFe = None
        self.SolubilityCr = None
        self.SolubilityCo = None
        self.SolubilityNi = None

        self.InnerOxThickness = None
        self.InnerIronOxThickness = None
        self.OuterFe3O4Thickness = None
        self.NiThickness = None
        self.OuterOxThickness = None
        self.CoThickness = None
        self.OxThickness = None

        self.StandardEqmPotentialFe = [float(SizingParametersReader[j + 13][i])\
                                       for i in range(self.RowStart, self.RowEnd)]
        self.StandardEqmPotentialFe3O4red = [float(SizingParametersReader[j + 14][i]) \
                                            for i in range(self.RowStart, self.RowEnd)]
        self.StandardEqmPotentialFe3O4oxid = [float(SizingParametersReader[j + 15][i]) \
                                            for i in range(self.RowStart, self.RowEnd)]
        self.StandardEqmPotentialH2 = [float(SizingParametersReader[j + 16][i]) \
                                    for i in range(self.RowStart, self.RowEnd)]
        self.StandardEqmPotentialNi = [float(SizingParametersReader[j + 17][i]) \
                                    for i in range(self.RowStart, self.RowEnd)]

        self.km = None
        self.KpFe3O4electrochem = None
        self.KdFe3O4electrochem = None
        self.CorrRate = None

        self.DepositThickness = None
        self.Particle = []
        self.FractionFeInnerOxide = None
        self.FractionNiInnerOxide = None

        self.ElapsedTime = None
        self.SpallTime = []
        self.TubeNumber = None

 
# creates the PHTS sections and their methods based on the Sections class template/blueprint
InletFeeder = Section(21, 0, 7)
InletFeeder_2 = Section(21, 0, 7)

FuelChannel = Section(42, 0, 12)
FuelChannel_2 = Section(42, 0, 12)

OutletFeeder = Section(63, 0, 9)
OutletFeeder_2 = Section(63, 0, 9)

# each steam generator split into 87 distinct bundles based on the tube bend arc lengths
SteamGenerator = [Section(84, 0, 22) for each in range(87)]
SteamGenerator_2 = [Section(84, 0, 22) for each in range(87)]

# only one feeder and fuel channel in each section for now 
InletSections = [InletFeeder, InletFeeder_2]
OutletSections = [OutletFeeder, OutletFeeder_2]
FuelSections = [FuelChannel, FuelChannel_2]
SteamGeneratorSections = [SteamGenerator, SteamGenerator_2]

for InletPiping in InletSections:
    InletPiping.Diameter = [44.3, 50, 106, 5.68, 5.68, 5.68, 5.68]
    InletPiping.Velocity = [1530, 1200, 270, 985, 985, 985, 985]
    InletPiping.Length.magnitude = [477.6, 281.8, 78.6, 350, 350, 350, 350]
    # Solubility (mol/kg)
    InletPiping.SolubilityNi = [2.71098E-09] * InletPiping.NodeNumber
    InletPiping.SolubilityCo = [3.77742E-09] * InletPiping.NodeNumber
    InletPiping.SolubilityCr = [8.81E-11] * InletPiping.NodeNumber
    InletPiping.PrimaryBulkTemperature = UnitConverter(
    InletPiping, "Celsius", "Kelvin", None, None, None, None, None, [266] * InletPiping.NodeNumber
    )

for Channel in FuelSections:
    Channel.Diameter = [1.3] * Channel.NodeNumber 
    Channel.Velocity = [883.08, 890.66, 900.3, 910.64, 920.97, 932.68, 945.08, 958.17, 973.32, 989.16, 1073.89, 1250.92]
    Channel.Length.magnitude = [49.5] * Channel.NodeNumber
    
    Channel.SolubilityNi = [
        2.71098E-09, 2.66196E-09, 2.60372E-09, 2.54535E-09, 2.29445E-09, 2.2137E-09, 1.91294E-09, 1.82788E-09,
        1.54081E-09, 1.54384E-09, 1.54584E-09, 1.54584E-09
        ]
    Channel.SolubilityCo = [
        3.77742E-09, 3.51525E-09, 3.20371E-09, 2.8915E-09, 2.60041E-09, 2.3011E-09, 2.08369E-09, 1.86963E-09,
        1.67578E-09, 1.53187E-09, 1.43654E-09, 1.43654E-09
        ]
    Channel.SolubilityCr = [
        8.81E-11, 9.61E-11, 1.01E-10, 9.40E-11, 8.69E-11, 7.98E-11, 7.28E-11, 6.57E-11, 5.86E-11, 5.16E-11, 4.69E-11,
        4.69E-11
        ]
    Channel.PrimaryBulkTemperature = UnitConverter(
        Channel, "Celsius", "Kelvin", None, None, None, None, None, [266.55, 270.48, 275.15, 279.83, 284.51, 289.19,
                                                                   293.87, 298.54, 303.22, 307.9, 310, 310]
                                                   )

for OutletPiping in OutletSections: 
    OutletPiping.Diameter = [6.4, 6.4, 6.4, 6.4, 8.9, 8.9, 8.9, 116, 40.8] 
    OutletPiping.Velocity = [1619, 1619, 1619, 1619, 857, 857, 857, 306, 1250]
    OutletPiping.Length.magnitude = [17, 3.5, 139.5, 432, 225.5, 460.3, 460.3, 400, 100]
    OutletPiping.SolubilityNi = [1.54584E-09] * OutletPiping.NodeNumber
    OutletPiping.SolubilityCo = [1.44E-09] * OutletPiping.NodeNumber
    OutletPiping.SolubilityCr = [4.84E-11] * OutletPiping.NodeNumber
    OutletPiping.PrimaryBulkTemperature = UnitConverter(
    OutletPiping, "Celsius", "Kelvin", None, None, None, None, None, [310] * OutletPiping.NodeNumber
    )

# assumed that these u-bends/straight leg lengths are = for all steam generators 
u_bend = []
straight_u_bend_section = [9.5] * 82 + [7.95] + [6.14] + [3.94] + [1.4] + [0]  # in.
straight_u_bend_section = [i * 2.54 for i in straight_u_bend_section]  # in. to cm
radius_decrease = [0.475] * 81 + [0.17] + [0.06] + [-0.16] + [-0.32] + [0.25]
radii = list(range(87))

for i, j in zip (radii, radius_decrease):
    if i == 0:
        # convert from in. to cm 
        # 2piR = circumferance of circle, piR = half of circle
        x = 43.225
        u_bend.append(x * 2.54 * np.pi)
    # in. # each u-bend radius is 0.475 in. shorter than the previous (slightly deviates at the end)
    x = x - j
    y = (x * 2.54) * np.pi
    u_bend.append(y)
    
u_bend_total = [x + y for x, y in zip(u_bend, straight_u_bend_section)]

hot_leg_length = [41, 157, 93, 189, 96, 99, 96, 99]
cold_leg_length = [136, 95, 99, 97, 178, 29, 92, 66, 42, 32]

number_tubes = [8, 11, 14, 15, 18, 19, 20, 21, 22, 23, 24, 25, 26, 26, 28, 28, 28, 29, 28, 31, 32, 33, 34, 35, 36, 35,
                36, 37, 36, 39, 38, 39, 40, 41, 40, 36, 40, 43, 42, 43, 44, 43, 44, 45, 42, 45, 46, 45, 46, 47, 46, 47,
                48, 47, 48, 43, 48, 48, 49, 50, 49, 50, 51, 50, 51, 50, 51, 50, 51, 52, 50, 46, 51, 50, 53, 52, 51, 52,
                53, 52, 53, 52, 51, 50, 50, 48, 47]


def steam_generator_properties(SteamGenerator):
    
    for Zone, length, i in zip(SteamGenerator, u_bend_total, number_tubes):
        # u-bend split into 4 nodes, but length is just a float 
        Zone.Length.magnitude = hot_leg_length + [length / 4] * 4 + cold_leg_length
        Zone.TubeNumber = i
        
        Zone.Diameter = [1.368] * Zone.NodeNumber
        Zone.Velocity = [
            533.002, 533.001, 533, 531, 524, 517, 511, 506, 506, 502, 498, 494, 491, 489, 487, 484, 483, 481, 480, 479,
            476, 474
            ]
        
        Zone.SolubilityNi = [
            1.5452E-09, 1.5452E-09, 1.5452E-09, 1.54453E-09, 1.54189E-09, 1.78271E-09, 1.84719E-09, 1.9062E-09,
            1.96011E-09, 2.22698E-09, 2.27478E-09, 2.31567E-09, 2.35035E-09, 2.3821E-09, 2.41091E-09, 2.59037E-09,
            2.60733E-09, 2.62118E-09, 2.63802E-09, 2.66147E-09, 2.68978E-09, 2.71747E-09
            ]
        Zone.SolubilityCo = [
            1.46729E-09, 1.46729E-09, 1.46729E-09, 1.49896E-09, 1.62443E-09, 1.75595E-09, 1.91821E-09, 2.06673E-09,
            2.2024E-09, 2.35035E-09, 2.5275E-09, 2.67907E-09, 2.80762E-09, 2.92529E-09, 3.0321E-09, 3.13232E-09,
            3.22305E-09, 3.2971E-09, 3.38716E-09, 3.51258E-09, 3.66401E-09, 3.81211E-09
            ]
        Zone.SolubilityCr = [
            4.84E-11, 4.84E-11, 4.84E-11, 4.99E-11, 5.61E-11, 6.20E-11, 6.73E-11, 7.22E-11, 7.67E-11, 8.10E-11,
            8.52E-11, 8.88E-11, 9.18E-11, 9.46E-11, 9.71E-11, 9.94E-11, 1.01E-10, 1.03E-10, 1.00E-10, 9.62E-11,
            9.15E-11, 8.70E-11
                             ]
        
        Zone.PrimaryBulkTemperature = UnitConverter(
            Zone, "Celsius", "Kelvin", None, None, None, None, None,
            [310.002, 310.001, 310, 308.97, 304.89, 301.02, 297.48, 294.24, 291.28, 288.42, 285.65, 283.28, 281.27,
             279.43, 277.76, 276.22, 274.86, 273.75, 272.4, 270.52, 268.25, 266.03]
            )
        
        Zone.Length.label = ["PHT boiling"] * 2 \
        + [None] * 6 \
        + ["u-bend"] * 4 \
        + [None] * 5 \
        + ["preheater start"] \
        + ["preheater"] * 3\
        + ["thermal plate"]

steam_generator_properties(SteamGenerator)
steam_generator_properties(SteamGenerator_2)

# Combines PHT sections and SG Zones (in the event each zone will be tracked for oxide growth/heat transfer)
FullLoop = InletSections + FuelSections + OutletSections + SteamGenerator + SteamGenerator_2
HalfLoop = [InletSections[0], FuelSections[0], OutletSections[0]] + SteamGeneratorSections[0]
SingleTubeLoop = [InletSections[0], FuelSections[0], OutletSections[0], SteamGenerator[57]]

for Section in FullLoop:
    # Particulate #[mg/kg] (ppm)
    Section.SmallParticulate = [0] * Section.NodeNumber
    Section.BigParticulate = [0] * Section.NodeNumber

    # Oxide thicknesses [g/cm^2]
    if Section in SteamGenerator or Section in SteamGenerator_2:
        Section.OuterFe3O4Thickness = [1.3E-4] * Section.NodeNumber
        Section.NiThickness = [1.3E-4] * Section.NodeNumber
        Section.OuterOxThickness = [1 * x + 1 * y for x, y in zip(Section.OuterFe3O4Thickness, Section.NiThickness)]

        Section.TubeThickness = 0.113
        Section.OuterDiameter = [x + 2 * Section.TubeThickness for x in Section.Diameter]

    if Section in OutletSections or Section in InletSections:
        Section.OuterFe3O4Thickness = [2.5E-4] * Section.NodeNumber
        Section.NiThickness = [0] * Section.NodeNumber
        Section.OuterOxThickness = [i * 1 for i in Section.OuterFe3O4Thickness]

    Section.InnerOxThickness = [2.5E-4] * Section.NodeNumber
    Section.InnerIronOxThickness = [i * 1 for i in Section.InnerOxThickness]
    Section.CoThickness = [0] * Section.NodeNumber
    
    if Section in FuelSections:
        Section.InnerOxThickness = [0] * Section.NodeNumber
        Section.OuterFe3O4Thickness = [0] * Section.NodeNumber
        Section.NiThickness = [0] * Section.NodeNumber
        Section.OuterOxThickness = [i * 1 for i in Section.OuterFe3O4Thickness]
    
    Section.OxThickness = [x + y for x, y in zip(Section.InnerOxThickness, Section.OuterOxThickness)]

    Section.Distance = np.cumsum(Section.Length.magnitude)

    
def ReynoldsNumber(Section, Diameter):
    # Diameter is an input due to difference in desired dimension (e.g., inner, outer, hydraulic, etc.)
    # [cm/s][cm][g/cm^3]/([g/cm s]
    Reynolds = [x * y * q / z  
                for x, y, z, q in zip(Section.Velocity, Diameter, Section.ViscosityH2O, Section.DensityH2O)]
    return Reynolds


def MassTransfer(Section):
    Schmidt = [x / (y * nc.FeDiffusivity) for x, y in zip(Section.ViscosityH2O, Section.DensityH2O)]
    Reynolds = ReynoldsNumber(Section, Section.Diameter)
    Sherwood = [0.0165 * (x ** 0.86) * (y ** 0.33) for x, y in zip(Reynolds, Schmidt)]  
    # Berger & Hau for straight pipe, single phase, fully developed (turbulent) flow
# 
#     SurfaceRoughness = 0.000075  # [m]
#     # [unitless]
#     HydraulicResistance = [(1.8 * np.log10((6.9 / x) + (SurfaceRoughness / (3.75 * y)) ** 1.11)) ** (-2) for x, y in
#                            zip(Reynolds, Section.Diameter)]
#     
#     Sherwood = [(x/8) * y * z / (1.07 + np.sqrt(x / 8) * ((z**0.667) - 1)) for x, y, z in zip(
#         HydraulicResistance, Reynolds, Schmidt)]
#     
#     r_elbow = 0.382 #0.09225 # [m]
#     
#     GeometryFactor= [0.68 + (1.2 - 0.044 * np.log(x)) * np.exp(-0.065 * r_elbow / y) + 0.58/(np.log(z + 2.5)) for
#                       x, y, z in zip(Reynolds, Section.Diameter, Schmidt)]
    GeometryFactor = [1.3] * Section.NodeNumber
    h_BH =  [z* nc.FeDiffusivity * x / y for x, y, z in zip(Sherwood, Section.Diameter, GeometryFactor)] # [cm/s]
    
#     if Section in OutletSections:
#         h_BH[1] = h_BH[1] * GeometryFactor[1]
    
    return h_BH
