import numpy as np
import csv 
import constants as nc


class SGParameters():
    def __init__(self):
        self.magnitude = None
        self.label = None
        self.unit = None
        self.steam_quality = None

PrimarySidePressure = 10  # MPa
SecondarySidePressure = 4.593  # MPa

SizingParameters = open('SizingParameters.txt', 'r')      
SizingParametersReader = list(csv.reader(SizingParameters, delimiter=','))  
# Assign file data to list, Reader[row][column], delimiter = comma 

n_IAPWS = [
    0.14632971213167, -0.84548187169114, -0.37563603672040 * 10 ** 1, 0.33855169168385 * 10 ** 1, -0.95791963387872, 
    0.15772038513228, -0.16616417199501 * 10 ** (-1), 0.81214629983568 * 10 ** (-3),
    0.28319080123804 * 10 ** (-3), -0.60706301565874 * 10 ** (-3), -0.18990068218419 * 10 ** (-1), 
    -0.32529748770505 * 10 ** (-1), -0.21841717175414 * 10 ** (-1), -0.52838357969930 * 10 ** (-4),
    - 0.47184321073267 * 10 ** (-3), -0.30001780793026 * 10 ** (-3), 0.47661393906987 * 10 ** (-4), 
    -0.44141845330846 * 10 ** (-5), -0.72694996297594 * 10 ** (-15), -0.31679644845054 * 10 ** (-4),
    - 0.28270797985312 * 10 ** (-5), -0.85205128120103 * 10 ** (-9), -0.22425281908000 * 10 ** (-5), 
    -0.65171222895601 * 10 ** (-6), -0.14341729937924 * 10 ** (-12), -0.40516996860117 * 10 ** (-6),
    - 0.12734301741641 * 10 ** (-8), -0.17424871230634 * 10 ** (-9), -0.68762131295531 * 10 ** (-18), 
    0.14478307828521 * 10 ** (-19), 0.26335781662795 * 10 ** (-22), -0.11947622640071 * 10 ** (-22),
    0.18228094581404 * 10 ** (-23), -0.93537087292458 * 10 ** (-25)
    ]

I_IAPWS = [0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 3, 3, 3, 4, 4, 4, 5, 8, 8, 21, 23, 29, 30, 31, 32]
J_IAPWS = [
    -2, -1, 0, 1, 2, 3, 4, 5, -9, -7, -1, 0, 1, 3 , -3, 0, 1, 3, 17, -4, 0, 6, -5, -2, 10, -8, -11, -6, -29, -31, -38, 
    -39, -40, -41
    ]

n_backwards_IAPWS = [
    - 0.23872489924521 * 10 ** 3, 0.40421188637945 * 10 ** 3, 0.11349746881718 * 10 ** 3, -0.58457616048039 * 10 ** 1,
    -0.15285482413140 * 10 ** (-3), -0.10866707695377 * 10 ** (-5), -0.13391744872602 * 10 ** 2,
    0.43211039183559 * 10 ** 2, -0.54010067170506 * 10 ** 2, 0.30535892203916 * 10 ** 2, -0.65964749423638 * 10 ** 1, 
    0.93965400878363 * 10 ** (-2), 0.11573647505340 * 10 ** (-6), -0.25858641282073 * 10 ** (-4), 
    -0.40644363084799 * 10 ** (-8), 0.66456186191635 * 10 ** (-7), 0.80670734103027 * 10 ** (-10), 
    -0.93477771213947 * 10 ** (-12), 0.58265442020601 * 10 ** (-14), -0.15020185953503 * 10 ** (-16)
    ]

I_backwards_IAPWS = [0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 2, 2, 3, 3, 4, 5, 6]

J_backwards_IAPWS = [0, 1, 2, 6, 22, 32, 0, 1, 2, 3, 4, 10, 32, 10, 32, 10, 32, 32, 32, 32] 

R_IAPWS = 0.461526  # [kJ/kg K] 

i_1 = [0, 1, 2, 3, 0, 1, 2, 3, 5, 0, 1, 2, 3, 4, 0, 1, 0, 3, 4, 3, 5]
j_1 = [0, 0, 0, 0, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 3, 3, 4, 4, 5, 6, 6]
H_1 = [
    0.520094, 0.0850895, -1.08374, -0.289555, 0.222531, 0.999115, 1.88797, 1.26613, 0.120573, -0.281378, -0.906851, 
    -0.772479, -0.489837, -0.25704, 0.161913, 0.257399, - 0.0325372, 0.0698452, 0.00872102, -0.00435673, -5.93E-04
    ]


class ThermoDimensionlessInput():
    def __init__(self):
        self.p_ref = 16.53  # [MPa] reducing pressure [MPa]
        self.T_ref = 1386  # [K] reducing temperature [K]
        self.rho_ref = 0.3220  # [g/cm^3]
        self.mu_ref = 1  # [micro Pa s]
        self.lambda_ref = 1  # [mW/mK]
        self.enthalpy = 2500  # [J/g]


ReferenceValues = ThermoDimensionlessInput()


def Enthalpy(side, Temperature, SecondarySidePressure):
    if side == "PHT" or side == "PHTS" or side == "phts" or side == "pht":
        p = PrimarySidePressure  # MPa
    else:
        p = SecondarySidePressure  # MPa secondary side

    ratio_pressures = p / ReferenceValues.p_ref  # p/p*, reduced pressure [unitless]
    ratio_temperatures = ReferenceValues.T_ref / Temperature  # T/T*, reduced temperature [unitless]

    Gibbs_T = sum(
        [x * ((7.1 - ratio_pressures) ** y) * z * ((ratio_temperatures - 1.222) ** (z - 1))
        for x, y, z in zip(n_IAPWS, I_IAPWS, J_IAPWS)]
        )
    return ratio_temperatures * Gibbs_T * R_IAPWS * Temperature  # [J/g mol]*[K] = [J/g]


def TemperaturefromEnthalpy(side, Enthalpy, SecondarySidePressure):
    if side == "PHT" or side == "PHTS" or side == "phts" or side == "pht":
        p = PrimarySidePressure
    else:
        p = SecondarySidePressure  # MPa secondary side
    p_ref = 1  # MPa
    T_ref = 1  # K
    
    ratio_pressures = p / p_ref  # p/p*, reduced pressure [unitless]
    ratio_enthalpies = Enthalpy / ReferenceValues.enthalpy
    
    ratio_temmperatures = sum(
        [x * ((ratio_pressures) ** y) * (ratio_enthalpies + 1) ** z
        for x, y, z in zip(n_backwards_IAPWS, I_backwards_IAPWS, J_backwards_IAPWS)]
        )
    
    return ratio_temmperatures * T_ref

    
def Density(species, side, Temperature, SecondarySidePressure):
    if side == "phts" or side == "PHTS" or side == "PHT":
        p = PrimarySidePressure
    elif side == "shts" or side == "SHTS" or side == "SHT":
        p = SecondarySidePressure 
    else: print ("Error: HTS side not specified")
        
    ratio_pressures = p / ReferenceValues.p_ref
    ratio_temperatures = ReferenceValues.T_ref / Temperature  # ratio is reverse here from all other functions 
    
    if species == "water":    
        Gibbs_p = sum(
            [-x * y * ((7.1 - ratio_pressures) ** (y - 1)) * (ratio_temperatures - 1.222) ** z
            for x, y, z in zip(n_IAPWS, I_IAPWS, J_IAPWS)]
            )
      
        return (1 / (R_IAPWS * Temperature * ratio_pressures * (Gibbs_p / (p * 1000)))) * 1000 / (100 ** 3)  # [g/cm^3] 
        
        
def Viscosity(species, side, Temperature):
    T_ref = 647.096  # K
    rho_ref = ReferenceValues.rho_ref
    mu_ref = ReferenceValues.mu_ref
    
    T_rel = Temperature / T_ref
    rho_rel = Density(species, side, Temperature, SecondarySidePressure) / rho_ref
    
    H_0 = [1.67752, 2.20462, 0.6366564, -0.241605]
    terms = [0, 1, 2, 3]
    summation = sum([x / (T_rel) ** y for x, y in zip(H_0, terms)])
    
    mu_0 = 100 * np.sqrt(T_rel) / summation
    
    summation1 = sum([x * (((1 / T_rel) - 1) ** y) * ((rho_rel - 1) ** z) for x, y, z in zip(H_1, i_1, j_1)])
    mu_1 = np.exp(rho_rel * summation1)
    mu_2 = 1
    # [microPa s] --> [g/cm s]
    return (mu_ref * (mu_0 * mu_1 * mu_2) / 1000000) * 10  # [g/cm s]


def HeatCapacity(side, Temperature, SecondarySidePressure):
    if side == "PHT" or side == "PHTS" or side == "phts" or side == "pht":
        p = PrimarySidePressure
    else:
        p = SecondarySidePressure
    
    ratio_pressures = p / ReferenceValues.p_ref  # p/p*, reduced pressure [unitless]
    ratio_temperatures = ReferenceValues.T_ref / Temperature  # T/T*, reduced temperature [unitless]
    
    Gibbs_TT = sum(
        [x * ((7.1 - ratio_pressures) ** y) * z * (z - 1) * ((ratio_temperatures - 1.222) ** (z - 2))
        for x, y, z in zip(n_IAPWS, I_IAPWS, J_IAPWS)]
        )
    
    return (-ratio_temperatures ** 2) * Gibbs_TT * R_IAPWS  # [kJ/kg K] ([J/g K] or *nc.H2OMolarMass for [J/mol K]
 
 
def thermal_conductivityH2O(side, Temperature, SecondarySidePressure):
#     if side == "PHT" or side == "PHTS" or side == "phts" or side=="pht":
#         p = PrimarySidePressure
#     else:
#         p = SecondarySidePressure
    T_ref = 647.096  # [K]
#     p_ref = 22.064 #[MPa]

    ratio_densities = Density("water", side, Temperature, SecondarySidePressure) / ReferenceValues.rho_ref
    ratio_temperatures = Temperature / T_ref
    # ratio_pressures = p/p_ref

    L_IAPWS = [2.443221E-03, 1.323095E-02, 6.770357E-03, -3.454586E-03, 4.096266E-04]
    k_IAPWS = [0, 1, 2, 3, 4]
    summation = sum([x / (ratio_temperatures ** y) for x, y in zip(L_IAPWS, k_IAPWS)])

    lambda_0 = np.sqrt(ratio_temperatures) / summation

    i_IAPWS = [0, 1, 2, 3, 4]
    j_IAPWS = [0, 1, 2, 3, 4, 5]
    L1_IAPWS = [
        [1.60397357, -0.646013523, 0.111443906, 0.102997357, -0.050412363, 0.006098593], 
        [2.33771842, -2.78843778, 1.53616167, -0.463045512, 0.083282702, -0.007192012],
        [2.19650529, -4.54580785, 3.55777244, -1.40944978, 0.275418278, -0.020593882], 
        [-1.21051378, 1.60812989, -0.621178141, 0.071637322, 0, 0],
        [-2.720337, 4.57586331, -3.18369245, 1.1168348, -0.19268305, 0.012913842]
        ]

    summation1 = []
    for row in L1_IAPWS:
        y = sum([x * (ratio_densities - 1) ** y for x, y in zip(row, j_IAPWS)])
        summation1.append(y)

    lambda_1 = np.exp(
        ratio_densities * sum([(((1 / ratio_temperatures) - 1) ** x) * y for x, y in zip(i_IAPWS, summation1)])
        )

    return ((lambda_0 * lambda_1 * ReferenceValues.lambda_ref) / 1000) / 100  # [W/m K] 


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
        if Section == Inlet or Section == Outlet:
            return [i * ((1 / nc.FeDensity) * 3600 * 24 * 365 * 10000) for i in Rate]
        elif Section in SGZones:
            return [x * y for x, y in zip(Rate, [(1 / nc.Alloy800Density) * 3600 * 24 * 365 * 10000] 
                                          * Section.NodeNumber)]
        else:
            return None
    
    # Converts from oxide mass per surface area to oxide thickness [g/cm^2] to [um]
    # assuming uniform distribution on surface   
    elif UnitInput == "Oxide Thickness Grams" and UnitOutput == "Oxide Thickness Micrometers":
        return [x * 10000 / y for x, y in zip(Oxide, OxideDensity)]

    else:
        None


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
        self.Length = SGParameters() 
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

 
# Creates the 4 PHTS sections and their methods based on the Sections class template/blueprint
Inlet = Section(21, 0, 7)
Core = Section(42, 0, 12)
Outlet = Section(63, 0, 9)
# Steam generator split into 87 zones based on the distinct tube bend arc lengths
SGZones = [Section(84, 0, 22) for each in range(87)]

# Diameter [cm]
Inlet.Diameter = [44.3, 50, 106, 5.68, 5.68, 5.68, 5.68] 
Core.Diameter = [1.3] * Core.NodeNumber 
Outlet.Diameter = [6.4, 6.4, 6.4, 6.4, 8.9, 8.9, 8.9, 116, 40.8] 

# Velocity [cm/s]
Inlet.Velocity = [1530, 1200, 270, 985, 985, 985, 985]
Core.Velocity = [883.08, 890.66, 900.3, 910.64, 920.97, 932.68, 945.08, 958.17, 973.32, 989.16, 1073.89, 1250.92]
Outlet.Velocity = [1619, 1619, 1619, 1619, 857, 857, 857, 306, 1250]

# Length [cm]
Inlet.Length.magnitude = [477.6, 281.8, 78.6, 350, 350, 350, 350] 
Core.Length.magnitude = [49.5] * Core.NodeNumber 
Outlet.Length.magnitude = [17, 3.5, 139.5, 432, 225.5, 460.3, 460.3, 400, 100]

# Solubility (mol/kg)
Inlet.SolubilityNi = [2.71098E-09] * Inlet.NodeNumber
Inlet.SolubilityCo = [3.77742E-09] * Inlet.NodeNumber
Inlet.SolubilityCr = [8.81E-11] * Inlet.NodeNumber

Core.SolubilityNi = [
    2.71098E-09, 2.66196E-09, 2.60372E-09, 2.54535E-09, 2.29445E-09, 2.2137E-09, 1.91294E-09, 1.82788E-09, 1.54081E-09, 
    1.54384E-09, 1.54584E-09, 1.54584E-09
    ]
Core.SolubilityCo = [
    3.77742E-09, 3.51525E-09, 3.20371E-09, 2.8915E-09, 2.60041E-09, 2.3011E-09, 2.08369E-09, 1.86963E-09, 1.67578E-09, 
    1.53187E-09, 1.43654E-09, 1.43654E-09
    ]
Core.SolubilityCr = [
    8.81E-11, 9.61E-11, 1.01E-10, 9.40E-11, 8.69E-11, 7.98E-11, 7.28E-11, 6.57E-11, 5.86E-11, 5.16E-11, 4.69E-11, 
    4.69E-11
    ]

Outlet.SolubilityNi = [1.54584E-09] * Outlet.NodeNumber
Outlet.SolubilityCo = [1.44E-09] * Outlet.NodeNumber
Outlet.SolubilityCr = [4.84E-11] * Outlet.NodeNumber
 
u_bend = []
straight_u_bend_section = [9.5] * 82 + [7.95] + [6.14] + [3.94] + [1.4] + [0] # in.
straight_u_bend_section = [i * 2.54 for i in straight_u_bend_section] # in. to cm
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

hot_leg_length = [74.05, 74.05, 148.1, 113.96, 113.96, 113.96, 113.96, 113.96]
cold_leg_length = [113.96, 113.96, 113.96, 113.96, 113.96, 74.05, 75, 74.05, 42.05, 32]

number_tubes = [8, 11, 14, 15, 18, 19, 20, 21, 22, 23, 24, 25, 26, 26, 28, 28, 28, 29, 28, 31, 32, 33, 34, 35, 36, 35,
                36, 37, 36, 39, 38, 39, 40, 41, 40, 36, 40, 43, 42, 43, 44, 43, 44, 45, 42, 45, 46, 45, 46, 47, 46, 47,
                48, 47, 48, 43, 48, 48, 49, 50, 49, 50, 51, 50, 51, 50, 51, 50, 51, 52, 50, 46, 51, 50, 53, 52, 51, 52,
                53, 52, 53, 52, 51, 50, 50, 48, 47]

# not including last appended "SteamGenerator" zone
for Zone, length, i in zip(SGZones, u_bend_total, number_tubes):
    # u-bend split into 4 nodes
    Zone.Length.magnitude = hot_leg_length + [length / 4] * 4 + cold_leg_length
    Zone.TubeNumber = i
    
    Zone.Diameter = [1.368] * Zone.NodeNumber
    Zone.Velocity = [
        533.002, 533.001, 533, 531, 524, 517, 511, 506, 506, 502, 498, 494, 491, 489, 487, 484, 483, 481, 480, 479, 476, 
        474
        ]
    
    Zone.SolubilityNi = [
        1.5452E-09, 1.5452E-09, 1.5452E-09, 1.54453E-09, 1.54189E-09, 1.78271E-09, 1.84719E-09, 1.9062E-09, 1.96011E-09, 
        2.22698E-09, 2.27478E-09, 2.31567E-09, 2.35035E-09, 2.3821E-09, 2.41091E-09, 2.59037E-09, 2.60733E-09, 
        2.62118E-09, 2.63802E-09, 2.66147E-09, 2.68978E-09, 2.71747E-09
        ]
    Zone.SolubilityCo = [
        1.46729E-09, 1.46729E-09, 1.46729E-09, 1.49896E-09, 1.62443E-09, 1.75595E-09, 1.91821E-09, 2.06673E-09,
        2.2024E-09, 2.35035E-09, 2.5275E-09, 2.67907E-09, 2.80762E-09, 2.92529E-09, 3.0321E-09, 3.13232E-09,
        3.22305E-09, 3.2971E-09, 3.38716E-09, 3.51258E-09, 3.66401E-09, 3.81211E-09
        ]
    Zone.SolubilityCr = [
        4.84E-11, 4.84E-11, 4.84E-11, 4.99E-11, 5.61E-11, 6.20E-11, 6.73E-11, 7.22E-11, 7.67E-11, 8.10E-11, 8.52E-11, 
        8.88E-11, 9.18E-11, 9.46E-11, 9.71E-11, 9.94E-11, 1.01E-10, 1.03E-10, 1.00E-10, 9.62E-11, 9.15E-11, 8.70E-11
                         ]
    
    Zone.PrimaryBulkTemperature = UnitConverter(
        Zone, "Celsius", "Kelvin", None, None, None, None, None, 
        [310.002, 310.001, 310, 308.97, 304.89, 301.02, 297.48, 294.24, 291.28, 288.42, 285.65, 283.28, 281.27, 279.43, 
         277.76, 276.22, 274.86, 273.75, 272.4, 270.52, 268.25, 266.03]
        )
    
    Zone.Length.label = ["PHT boiling"] * 2 \
    + [None] * 6 \
    + ["u-bend"] * 4 \
    + [None] * 5 \
    + ["preheater start"] \
    + ["preheater"] * 3

# Temperature [Celsius]
Inlet.PrimaryBulkTemperature = UnitConverter(
    Inlet, "Celsius", "Kelvin", None, None, None, None, None, [266] * Inlet.NodeNumber
    )
Core.PrimaryBulkTemperature = UnitConverter(
    Inlet, "Celsius", "Kelvin", None, None, None, None, None, [266.55, 270.48, 275.15, 279.83, 284.51, 289.19, 293.87, 
    298.54, 303.22, 307.9, 311, 311]
    )
Outlet.PrimaryBulkTemperature = UnitConverter(
    Inlet, "Celsius", "Kelvin", None, None, None, None, None, [310] * Outlet.NodeNumber
    )

# Combines PHT sections and SG Zones (in the event each zone will be tracked for oxide growth/heat transfer)
Sections = [Inlet, Core, Outlet] + SGZones

for Section in Sections:
    # Particulate #[mg/kg] (ppm)
    Section.SmallParticulate = [0] * Section.NodeNumber
    Section.BigParticulate = [0] * Section.NodeNumber

    # Oxide thicknesses [g/cm^2]
    if Section in SGZones:
        Section.OuterFe3O4Thickness = [1.3E-4] * Section.NodeNumber
        Section.NiThickness = [1.3E-4] * Section.NodeNumber
        Section.OuterOxThickness = [1 * x + 1 * y for x, y in zip(Section.OuterFe3O4Thickness, Section.NiThickness)]

        Section.TubeThickness = 0.113
        Section.OuterDiameter = [x + 2 * Section.TubeThickness for x in Section.Diameter]

    if Section == Core:
        Section.OuterFe3O4Thickness = [0] * Section.NodeNumber
        Section.NiThickness = [0] * Section.NodeNumber
        Section.OuterOxThickness = [i * 1 for i in Section.OuterFe3O4Thickness]

    if Section == Outlet or Inlet:
        Section.OuterFe3O4Thickness = [2.5E-4] * Section.NodeNumber
        Section.NiThickness = [0] * Section.NodeNumber
        Section.OuterOxThickness = [i * 1 for i in Section.OuterFe3O4Thickness]

    Section.InnerOxThickness = [2.5E-4] * Section.NodeNumber
    Section.InnerIronOxThickness = [i * 1 for i in Section.InnerOxThickness]
    Section.CoThickness = [0] * Section.NodeNumber
    Section.OxThickness = [x + y for x, y in zip(Section.InnerOxThickness, Section.OuterOxThickness)]

    # Distance
    Section.Distance = np.cumsum(Section.Length.magnitude)


def ReynoldsNumber(Section, Diameter):
    # Diameter is an input due to difference in desired dimension (e.g., inner, outer, hydraulic, etc.)
    Reynolds = [x * y / (z * q) 
                for x, y, z, q in zip(Section.Velocity, Diameter, Section.ViscosityH2O, Section.DensityH2O)]    
    return Reynolds


def MassTransfer(Section):
    Schmidt = [x / (y * nc.FeDiffusivity) for x, y in zip(Section.ViscosityH2O, Section.DensityH2O)]
    Reynolds = ReynoldsNumber(Section, Section.Diameter)
    Sherwood = [0.0165 * (x ** 0.86) * (y ** 0.33) for x, y in zip(Reynolds, Schmidt)]  
    # Berger & Hau for straight pipe, single phase, fully developed (turbulent) flow
    return [nc.FeDiffusivity * x / y for x, y in zip(Sherwood, Section.Diameter)] 
