import lepreau_data as ld
import numpy as np
import thermochemistry_and_constants as nc
import random

NumberPluggedTubes = 8
TotalSGTubeNumber = 3542 - NumberPluggedTubes

YearStartup = 1983.25
YearCPP = 1988.25
YearSHTChemicalClean = 1996

YearDerates = [1999.75, 2000, 2000.25, 2000.5, 2000.75, 2001, 2001.25, 2001.5, 2001.75, 2002, 2002.25, 20002.5,
               2002.75, 2003, 2003.25, 2003.5, 2003.75, 2004, 2004.25, 2004.5, 2004.75, 2005, 2005.25, 2005.5, 2005.75,
               2006, 2006.25, 2006.5, 2006.75, 2007, 2007.25, 2007.5, 2007.75, 2008.25]

PercentDerates = [3, 4, 4.8, 5, 5.5, 4.07, 3.23, 3.43, 3.77, 4.27, 4.61, 4.74, 3.33, 3.64, 3.15, 4.08, 4, 3.88, 3.74,
                  4.18, 4.97, 5.17, 5.24, 6.30, 6.97, 7.41, 7.89, 8.29, 8.5, 8.75, 8.87, 9.39, 9, 9, 9]

PercentDerates = [i - 1 for i in PercentDerates]


h_i = nc.SGParameters()
h_o = nc.SGParameters()
U_h = nc.SGParameters()
U_c = nc.SGParameters()
U_total = nc.SGParameters()
R_F_primary = nc.SGParameters()
R_F_secondary = nc.SGParameters()
MassFlow_preheater = nc.SGParameters()
MassFlow_downcomer = nc.SGParameters()
MassFlow_c_total = nc.SGParameters()
ShellDiameter = nc.SGParameters()
MassFlux_c = nc.SGParameters()
MassFlux_h = nc.SGParameters()
EnthalpySaturatedSteam = nc.SGParameters()
MassFlow_h = nc.SGParameters()
MasssFlow_dividerplate = nc.SGParameters()
EquivalentDiameter = nc.SGParameters()
TubePitch = nc.SGParameters()

for i in [h_i, h_o]:
    i.unit = "K/W"

for i in [U_h, U_c, U_total]:
    i.unit = "W/cm^2 K"

for i in [R_F_primary, R_F_secondary]:
    i.unit = "cm^2 K/W"

T_sat_primary = 310 + 273.15
T_PreheaterIn = 186.5 + 273.15
RecirculationRatio = 5.3

MassFlow_h.magnitude = 1900 * 1000
# Steam flow for 4 steam generators in typical CANDU-6 = 1033.0 kg/s
# 240 kg/s pulled from AECL COG document and works well with the 1900 kg/s hot-side flow ?
MassFlow_preheater.magnitude = 239.25 * 1000
MassFlow_ReheaterDrains = 89 * 1000 / 4
MassFlow_downcomer.magnitude = (RecirculationRatio - 1) * MassFlow_preheater.magnitude
MassFlow_c_total.magnitude = MassFlow_downcomer.magnitude + MassFlow_preheater.magnitude + MassFlow_ReheaterDrains


def thermalplate_temp():
    LeakageRate = 0.1
    MassFlow_recirculation = (1 - LeakageRate) * \
    (MassFlow_preheater.magnitude * RecirculationRatio - MassFlow_preheater.magnitude)
    Enthalpy_recirculation = nc.enthalpy("SHT", 258.66 + 273.15, 4.593)
    Enthalpy_preheater = nc.enthalpy("SHT", T_PreheaterIn, 4.593)
    
    TotalMass = MassFlow_preheater.magnitude * LeakageRate + MassFlow_recirculation
    Enthalpy_under_thermalplate = ((Enthalpy_preheater * (MassFlow_preheater.magnitude * LeakageRate)
                                   + MassFlow_recirculation * Enthalpy_recirculation) / TotalMass)
    
    T_under_thermalplate = nc.temperature_from_enthalpy("SHT", Enthalpy_under_thermalplate, 4.593)
    return T_under_thermalplate                               


ShellDiameter.magnitude = 2.28 * 100

for i in [MassFlow_c_total, MassFlow_h, MasssFlow_dividerplate, MassFlow_preheater, MassFlow_downcomer]:
    i.unit = "g/s"

for i in [ShellDiameter, EquivalentDiameter, TubePitch]:
    i.unit = "cm"


def MassFlux_c(Section, i):
    if Section.Length.label[i] == "opposite preheater":
        PreheaterDiameter = ShellDiameter.magnitude / 2  # [cm]
        MassFlow = (RecirculationRatio - 1) * MassFlow_preheater.magnitude
        TotalTubes = TotalSGTubeNumber
    else:
        PreheaterDiameter = 0
        MassFlow = MassFlow_c_total.magnitude
        TotalTubes = TotalSGTubeNumber  # * 2
        
    ShellCrossSectionalArea = (np.pi / 4) * (
    (ShellDiameter.magnitude ** 2)
    - (ld.SteamGenerator[0].OuterDiameter[0] ** 2) * TotalTubes
    - PreheaterDiameter ** 2
    )
    Flux = MassFlow / ShellCrossSectionalArea
    return Flux


for i in [MassFlux_c, MassFlux_h]:
    i.unit = "g/cm^2 s"
    
EnthalpySaturatedSteam.magnitude = 2727.5  # 9.89 MPa
EnthalpySaturatedSteam.unit = "J/g"

TubePitch.magnitude = 2.413


def cleaned_tubes():
    # chooses tube bundles until 60% of total sg tube number reached, adds all chosen to "cleaned" tube list
    Cleaned = []
    NumberTubes = []
    
    for i in range(len(ld.SteamGenerator)):
        x = random.randint(0, 86)
        NumberTubes.append(ld.SteamGenerator[x].TubeNumber)
        
        # siva blast used on only 60% of tubes due to time/spacial constraints
        if sum(NumberTubes) <= (0.6 * TotalSGTubeNumber):
            Cleaned.append(ld.SteamGenerator[x])
        else:
            break
    
    return Cleaned


def closest_tubelength(TubeLength):
    difference = []
    for Zone in ld.SteamGenerator:
        difference.append(abs(TubeLength - Zone.Distance[Zone.NodeNumber - 1]))
    return difference.index(min(difference))


# searches through all u-bend arc lengths and chooses the one closest to input
def closest_ubend(UbendLength):
    difference = []
    for i in ld.u_bend_total:
        # calculates differences between input Number and all others in given list
        difference.append(abs(UbendLength - i))
    
    # returns index of value that has smallest difference with input Number
    return difference.index(min(difference))


SGFastMode = "yes"
# select all the SG tubes to be run (by arc length)
UBends = [1.52, 1.71, 2.31, 3.09]
UBends = [i * 100 for i in UBends]
TubeLengths = [1892, 1907]
# [1892, 1907, 1907, 1889, 1907.952, 1968.59, 2044.4]

def tube_picker(method):
    tubes = []
    tube_number = []
    
    # selects class initializations based on desired u-bend tube arc lengths
    if method == "arc length": 
        for k in UBends:  # want multiple u-bends
            x = closest_ubend(k)
            tube_number.append(x)

        for i in tube_number:
            tubes.append(ld.SteamGenerator[i])
    # selects class initializations based on desired overall tube lengths        
    elif method == "tube length":
        for k in TubeLengths:
            x = closest_tubelength(k)
            tube_number.append(x)
        
        for i in tube_number:
            tubes.append(ld.SteamGenerator[i])
    
    else:
        None
        
    return tubes, tube_number

    
if SGFastMode == "yes":
    selected_tubes = tube_picker("tube length")[0]
else:
    selected_tubes = ld.SteamGenerator  # [0::3]

tube_number = tube_picker("tube length")[1] 

Cleaned = cleaned_tubes()


def thermal_conductivity(Twall, material, SecondarySidePressure):
    if material == "Alloy-800" or material == "Alloy800" or material == "A800" or material == "Alloy 800":
        Twall_C = Twall - 273.15  # oC in alloy thermal conductivity equatin 
        conductivity = (11.450 + 0.0161 * Twall_C) / 100  # M.K.E thesis for Alloy-800 tubing [W/cm K] 
        return conductivity
    
    elif material == "water":
        return nc.thermal_conductivityH2O("PHT", Twall, SecondarySidePressure)
    
    elif material == "inner magnetite":
        return 2 / 100  # [W/cm K]
    
    elif material == "outer magnetite":
        return 1.3 / 100  # [W/cm K]
    else:
        return None


def sludge_fouling_resistance(Section, i, calendar_year):
    
    TubeGrowth = 0.00085  # [g/cm^2] /year = 6.5 um /year
    ReducedTubeGrowth = 0.0001  # [g/cm^2] /year = 3.25 um/year
        
    # between 1983 and 1986 (included)
    if YearStartup <= calendar_year <= YearCPP:
        # secondary side fouling slope based on 1/2 that for average primary side cold leg deposit
        Accumulation = (calendar_year - YearStartup) * TubeGrowth
     
    # CPP installation (late 1986) reduces secondary side crud by 50% 
    # between 1987 and 1995, not including, maybe some form of dissolution event to historic deposits
    elif YearCPP < calendar_year < YearSHTChemicalClean:
        Accumulation = (calendar_year - YearCPP) * ReducedTubeGrowth
    
    elif calendar_year >= YearSHTChemicalClean:
        Accumulation = (calendar_year - YearSHTChemicalClean) * ReducedTubeGrowth
        
    Thickness = Accumulation / nc.Fe3O4Density
    
    Fouling = Thickness / thermal_conductivity(None, "outer magnetite", None)
    
    return Fouling
    

def pht_fouling_resistance(Section, i, calendar_year, InnerAccumulation, OuterAccumulation):
    if i == 21:
        ThermalPlateReduction = 0.5
    else:
        ThermalPlateReduction = 1
    # [g/cm^2]/[g/cm^3] = [cm]
    # thickness/thermal conductivity [cm]/[W/cm K] = [cm^2 K/W]
    InnerThickness = InnerAccumulation / nc.Fe3O4Density
    OuterThickness = ThermalPlateReduction * OuterAccumulation / nc.Fe3O4Density 

    # [cm]/ [W/cm K] =[cm^2 K/W]
    # inner deposit is consolidated, predicted to have different thermal resistance than unconsolidated outer ox.
    InnerFouling = InnerThickness / thermal_conductivity(None, "inner magnetite", None)
    OuterFouling = OuterThickness / thermal_conductivity(None, "outer magnetite", None)
    # total pht thermal resistance
    
    return InnerFouling + OuterFouling  # [cm^2 K/W]


def conduction_resistance(Section, Twall, SecondarySidePressure, i):
    # R_conduction = L/kA (L = thickness)
    # A = area with shape factor correcrion factor for cylindrical pipe
    # R_conduction = ln(D_o/D_i)/(2pikl) [K/W]
    # l = length, k =thermal conductivity coefficient [W/cm K]
    
    k_w = thermal_conductivity(Twall, "Alloy-800", SecondarySidePressure)  # [W/cm K]

    # small conductivity = big resistance
    Rcyl_numerator = np.log(Section.OuterDiameter[i] / Section.Diameter[i])
    Rcyl_denominator = 2 * np.pi * k_w * Section.Length.magnitude[i]

    R_cond = Rcyl_numerator / Rcyl_denominator  # [K/W]
    return R_cond


def primary_convection_resistance(Section, T_PrimaryBulkIn, T_SecondaryBulkIn, T_PrimaryWall, T_SecondaryWall,
                                  SecondarySidePressure, x_pht, i, m_h_leakagecorrection, TotalAccumulation):
    # R_1,Convective = 1/(h_1,convectove *A)
    # A = inner area (based on inner diameter)
    # [W/cm K]/ [cm] = [W/K cm^2]
    
    # T_PrimaryBulkIn used in place of Tavg (of in and out for that node)
    Nu_D = nusseltnumber(
            Section, "PHT", T_PrimaryBulkIn, T_SecondaryBulkIn, T_SecondaryWall, None,
            i, x_pht, m_h_leakagecorrection
            )
    h_l =  Nu_D * thermal_conductivity(T_PrimaryBulkIn, "water", None) / Section.Diameter[i]
    
    # PHT boiling 
    if x_pht > 0:
        MassFlux_h.magnitude = m_h_leakagecorrection / (TotalSGTubeNumber * (np.pi / 4) * (Section.Diameter[i] ** 2))
        
        h_i = boiling_heat_transfer(
            x_pht, "PHT", MassFlux_h.magnitude, T_PrimaryWall, Section.Diameter[i],
            SecondarySidePressure, i, h_l 
            )
         
    else:
        h_i =  h_l

    return 1 / (h_i * inner_area(Section, TotalAccumulation)[i])  # [K/W]


def secondary_convection_resistance(
        Section, T_PrimaryBulkIn, T_SecondaryBulkIn, T_SecondaryWall, SecondaryT_film, x_in,SecondarySidePressure, i
        ):
    
    # split into boiling and non-boiling (preheater) sections
    if Section.Length.label[i] == "preheater":  # from Silpsrikul thesis
        # Based on PLGS data and Silpsrikul calculations
        f_b = 0.2119  # fraction of cross-sectional area of shell occupied by a baffle window
        N_b = f_b * TotalSGTubeNumber  # number of tubes in baffle window
        D_s = ((4 * 0.5 * np.pi / 4) * (ShellDiameter.magnitude / 100) ** 2) \
            / (0.5 * np.pi * (ShellDiameter.magnitude / 100) + ShellDiameter.magnitude / 100)  # [m]

        S_b = (f_b * np.pi * (D_s ** 2) / 4) - N_b * np.pi * ((Section.OuterDiameter[i] / 100) ** 2) / 4  # [cm^2]
        G_b = (MassFlow_preheater.magnitude / 1000) / S_b  # [kg/m^2s]
        S_c = 0.5373 * D_s * (1 - (Section.OuterDiameter[i] / 100) / 0.0219)
        G_c = MassFlow_preheater.magnitude / 1000 / S_c
        G_e = np.sqrt(G_b * G_c) * 1000 / (100 ** 2)  # [g/cm^2]

        # First two (raised to exponent) terms are unitless
        # [W/cm K]/[cm] = [W/cm^2 K]
        k = thermal_conductivity(SecondaryT_film, "water", SecondarySidePressure)
        h_o = (k * 0.2 / Section.OuterDiameter[i]) \
            * ((Section.OuterDiameter[i] * G_e / nc.viscosity("SHT", SecondaryT_film, SecondarySidePressure)) ** 0.6) \
            * (
                nc.HeatCapacity("SHT", SecondaryT_film, SecondarySidePressure) \
                * nc.viscosity("SHT", SecondaryT_film, SecondarySidePressure) \
                / thermal_conductivity(SecondaryT_film, "water", SecondarySidePressure)
                ) ** 0.33
    
    # outside preheater, no cross-flow/baffles
    else:
        # Zukasuskas correlation-based heat transfer coefficient for a bundle, non boiling, not in pre-heater
         # [cm]
       
        EquivalentDiameter.magnitude = 4 * ((TubePitch.magnitude ** 2) \
                                            - (np.pi / 4) * (Section.OuterDiameter[i] ** 2))\
                                            / (np.pi * Section.OuterDiameter[i])
        
        Nu_D_l = nusseltnumber(
            Section, "SHT", T_PrimaryBulkIn, T_SecondaryBulkIn, T_SecondaryWall, SecondarySidePressure, i, None,
            None
            )
    
        h_l = nc.thermal_conductivityH2O("SHT", T_SecondaryBulkIn, SecondarySidePressure) * Nu_D_l \
        / EquivalentDiameter.magnitude
        
        if x_in > 0:
            MassFlux_c.magnitude = MassFlux_c(Section, i)
            
            h_o = boiling_heat_transfer(
                x_in, "SHT", MassFlux_c.magnitude, T_SecondaryWall, Section.OuterDiameter[i],
                SecondarySidePressure, i, h_l)
       
        else:
            h_o = h_l
    
    return 1 / (h_o * outer_area(Section)[i])  # [K/W]


def boiling_heat_transfer(x, side, MassFlux, T_wall, Diameter, SecondarySidePressure, i, h_l):
    
    # x = steam quality in percent
    if side == "SHT":
        if SecondarySidePressure == 4.593:
            Density_v = 0.023187  # [g/cm^3]
            T_sat_secondary = 258.69 + 273.15
        elif SecondarySidePressure < 4.593:
            Density_v = 0.022529  # [g/cm^3]
            T_sat_secondary = 257 + 273.15
        T_sat = T_sat_secondary
        
        Density_l = nc.density(side, T_sat, SecondarySidePressure)
        Viscosity_sat = nc.viscosity(side, T_sat, SecondarySidePressure)
        Pressure = SecondarySidePressure
        F = 1.04
        
        k_f = nc.thermal_conductivityH2O(side, T_sat, SecondarySidePressure)
        Cp = nc.HeatCapacity(side, T_sat, SecondarySidePressure)
        SurfaceTension_H2O = 0.024 * 1000#[N/m]
        
#         h = 0.00122 * (
#             (k_f ** 0.79) * (Cp ** 0.45) * (Density_l ** 0.49) \
#             / ((SurfaceTension_H2O ** 0.79) * (Viscosity_sat ** 0.79) * (1128 ** 0.79) * (Density_v ** 0.79))
#             ) * ((T_wall - T_sat) ** 0.24) * ((2 / 100) * 1000 *1000000)**0.24 

        
#         # start with Rohsenow correlation for boiling over a single tube
#         SurfaceTension_H2O = 0.024 #[N/m]
#         #[g/cm^3] * [(100cm)^3/m^3] * [kg/1000g] = [kg/m^3]
#         GravitationalForce = 9.81 #[m/s^2]
#         deltaDensity = (Density_l - Density_v) * (100 ** 3) / 1000
        # [N/m] = [[kg m / s^2] /m]
        #[kg/m^3] * [m/s^2] = [kg/m^2 s^2] / [kg/s^2]
    
        # [1/m^2] --> [1/cm^2]
#         bubble_quotient = ((GravitationalForce * deltaDensity) / SurfaceTension_H2O) / (100 ** 2) 
#         deltaT = 1.3
# 
#         Cp = nc.HeatCapacity(side, T_sat, SecondarySidePressure)
#         # correlation-specific constant from a Table (for stainless steel surface)
#         C_f_constant = 0.0132
#         n = 1
#         Pr = prandtl(side, T_sat, SecondarySidePressure)
#         Enthalpy = 1128.4 #[J/g]
#         #[g / cm s] * [1/cm] * [J/g] = [J / cm^2 s] = [W/cm^2] (first part)
#         # [J/g K] * [K] / ([J/g] * [unitless] * [unitless]) = [unitless]
#         q = (Viscosity_sat * np.sqrt(bubble_quotient) * Enthalpy * \
#         (Cp * deltaT / (Enthalpy * C_f_constant * Pr ** n)) ** 3)
#     
#         HeatTransferCoefficient = q / deltaT
#         print (q, deltaT, HeatTransferCoefficient, side, i)
    
    elif side == "PHT":
        # quality only persists for first 2 nodes of PHT
        T_sat = T_sat_primary
        Density_v = 1000 * 55.462 / (100 ** 3)  # [g/cm^3]
        Density_l = nc.D2O_density(T_sat)
        Viscosity_sat = nc.D2O_viscosity(T_sat)
        Pressure = nc.PrimarySidePressure
        F = (1 + (x * prandtl(side, T_sat, Pressure) * ((Density_l / Density_v) - 1))) ** 0.35
        
    p_crit = 22.0640  # [MPa]
    
    MassFlux_liquid = MassFlux * (1 - x)

    Re_D = Diameter * MassFlux_liquid / Viscosity_sat
    Q_prime = F * h_l * (abs(T_wall - T_sat))  # [W/cm^2]
    
    S = (1 + 0.055 * (F ** 0.1) * (Re_D) ** 0.16) ** (-1)
    A_p = (55 * ((4.70 / p_crit) ** 0.12) * ((-np.log10(4.70 / p_crit)) ** (-0.55)) * (nc.H2MolarMass) ** (-0.5)) \
            / (100 ** 2)
    
    C = ((A_p * S) / (F * h_l) ** 2) * (Q_prime) ** (4 / 3)
    coeff = [1, -C, -1]
    cubic_solution = np.roots(coeff)
    q = cubic_solution[0]
    
    HeatTransferCoefficient =  F * (q ** (3 / 2)) * h_l  # [W/cm^2 K]

    return HeatTransferCoefficient


def nusseltnumber(
        Section, side, T_PrimaryBulkIn, T_SecondaryBulkIn, T_SecondaryWall, SecondarySidePressure, i, x_pht,
        m_h_leakagecorrection):
   
    if side == "PHT":        
        MassFlux_h.magnitude = m_h_leakagecorrection / (TotalSGTubeNumber * (np.pi / 4) * (Section.Diameter[i] ** 2))
        
        #Dittus-boelter evaluated at average T between inlet and outlet for section - approximated as T_in here
        Re_D = [(MassFlux_h.magnitude / nc.D2O_viscosity(T_PrimaryBulkIn)) * i for i in Section.Diameter]
        
        n = 1 / 3
        C = 0.023
        m = 4 / 5
        
        Nu_D = C * (Re_D[i] ** m) * (prandtl(side, T_PrimaryBulkIn, SecondarySidePressure)) ** n

    elif side == "SHT":
        # [g/s] / [g/cm^3] = [cm^3/s] / [cm^2]
        A_Cross = (np.pi / 4) * ((ShellDiameter.magnitude ** 2) - (TotalSGTubeNumber * Section.OuterDiameter[i] ** 2))
        
        Velocity = MassFlow_c_total.magnitude / (nc.density(side, T_SecondaryBulkIn, SecondarySidePressure) * A_Cross)
        
        V_max = (TubePitch.magnitude / (TubePitch.magnitude - Section.Diameter[i])) * Velocity
        
        Re_D_max = V_max * Section.Diameter[i] * nc.density(side, T_SecondaryBulkIn, SecondarySidePressure) \
        / nc.viscosity(side, T_SecondaryBulkIn, SecondarySidePressure)
        
        C1 = 0.40
        m = 0.60
        n = 0.36
        
        Prandtl_ratio = prandtl(side, T_SecondaryBulkIn, SecondarySidePressure) \
        / prandtl(side, T_SecondaryWall, SecondarySidePressure)
        
        Nu_D = C1 * (Re_D_max ** m) * ((prandtl(side, T_SecondaryBulkIn, SecondarySidePressure)) ** n) \
        * Prandtl_ratio ** 0.25
        
    else:
        None
    
    return Nu_D


def prandtl(side, Temperature, SecondarySidePressure):
    
    if side == "PHT":
        Viscosity = nc.D2O_viscosity(Temperature)
        HeatCapacity = nc.HeatCapacity(side, Temperature, None)
        ThermalConductivity = nc.thermal_conductivityH2O("PHT", Temperature, None)
    
    elif side == "SHT":
        Viscosity = nc.viscosity(side, Temperature, SecondarySidePressure)
        HeatCapacity = nc.HeatCapacity(side, Temperature, SecondarySidePressure)
        ThermalConductivity = nc.thermal_conductivityH2O(side, Temperature, SecondarySidePressure)
    
    else:
        None
        
    Pr = HeatCapacity * Viscosity / ThermalConductivity
    
    return Pr


def inner_area(Section, Oxide):
    return [np.pi * x * y for x, y in zip(Section.Diameter, Section.Length.magnitude)]
#     return [np.pi * (x - z) * y for x, y, z in zip(Section.Diameter, Section.Length.magnitude, Oxide)]


def outer_area(Section):
    return [np.pi * x * y for x, y in zip(Section.OuterDiameter, Section.Length.magnitude)]


def pht_steam_quality(Temperature, j):
    calendar_year = (j * nc.TIME_STEP / 8760) + YearStartup
         
    CoreMassFlow = (1925000 / 1000) * 4  # [kg /s]
    Delta_T = T_sat_primary - (259.25 + 273.15)  # [K]
    C_p_cold = nc.HeatCapacity("PHT", 259.25, SecondarySidePressure=None)
    C_p_hot = nc.HeatCapacity("PHT", T_sat_primary, SecondarySidePressure=None)
    C_p_avg = (C_p_cold + C_p_hot) / 2  # [kJ/kg K]
    
    Power = CoreMassFlow * C_p_avg * Delta_T  # [kW]
    YearDeratingStarts = 1999.75
    
    difference = []
    
    if YearDeratingStarts <= calendar_year <= 2008.25:
        for Yr in YearDerates:
            difference.append(abs(calendar_year - Yr))
        ClosestYear = difference.index(min(difference))
        
        Percent_derating = PercentDerates[ClosestYear]
   
    else:
        Percent_derating = 0
        
    Power = Power - Power * (Percent_derating / 100)
#     print (Power/1000)
    # core inlet enthalpy at RIHT + that added from fuel
    # H_fromfuel = Power / CoreMassFlow # [kJ/s /kg/s] = [kJ/kg]
    H_pht = Power / CoreMassFlow
    H_satliq_outlet = nc.enthalpy("PHT", T_sat_primary, None)
    # no quality in return primary flow
    H_current = nc.enthalpy("PHT", Temperature, None) + H_pht
    x = (H_current - H_satliq_outlet) / (EnthalpySaturatedSteam.magnitude - H_satliq_outlet)
    
    if x < 0:
        x = 0
    return x
# print (pht_steam_quality(262.13 +273.15, 876*0))

def sht_steam_quality(Q, T_sat_secondary, x, MassFlow_c, SecondarySidePressure):
        
    # [J/s] / [g/s] = [J /g]
    H_pht = Q / (MassFlow_c)  # [kJ/kg]
    H_satliq = nc.enthalpy("SHT", T_sat_secondary, SecondarySidePressure)
    if SecondarySidePressure == 4.593:
        H_SaturatedSteam = 2797.3  # [kJ/kg]
    elif SecondarySidePressure < 4.593:
        H_SaturatedSteam = 2798.2
    H_prev = H_satliq + x * (H_SaturatedSteam - H_satliq)
    H_current = H_pht + H_prev

    x = (H_current - H_satliq) / (H_SaturatedSteam - H_satliq)
    return x


def wall_temperature(
        Section, i, T_PrimaryBulkIn, T_SecondaryBulkIn, x_in, x_pht, InnerAccumulation, OuterAccumulation,
        calendar_year, SecondarySidePressure, m_h_leakagecorrection):
    
    TotalAccumulation = [x + y for x, y in zip (InnerAccumulation, OuterAccumulation)]
    
    # i = each node of SG tube
    T_PrimaryWall = T_PrimaryBulkIn - (1 / 3) * (T_PrimaryBulkIn - T_SecondaryBulkIn)
    T_SecondaryWall = T_PrimaryBulkIn - (1 / 3) * (T_PrimaryBulkIn - T_SecondaryBulkIn)

    for k in range(50):
        WT_h = T_PrimaryWall
        WT_c = T_SecondaryWall

        PrimaryT_film = (T_PrimaryBulkIn + T_PrimaryWall) / 2
        SecondaryT_film = (T_SecondaryBulkIn + T_SecondaryWall) / 2

        h_i.magnitude = 1 / (primary_convection_resistance(
            Section, T_PrimaryBulkIn, T_SecondaryBulkIn, T_PrimaryWall, T_SecondaryWall, SecondarySidePressure, x_pht,
            i, m_h_leakagecorrection, TotalAccumulation) * inner_area(Section, TotalAccumulation)[i])
        
        h_o.magnitude = 1 / (secondary_convection_resistance(
            Section, T_PrimaryBulkIn, T_SecondaryBulkIn, T_SecondaryWall, SecondaryT_film, x_in, SecondarySidePressure,
            i) * outer_area(Section)[i])
        # R = 1/hA --> h=A*1/R

        U_h1 = inner_area(Section, TotalAccumulation)[i] \
        * conduction_resistance(Section, T_PrimaryWall, SecondarySidePressure, i)
        
        U_h2 = inner_area(Section, TotalAccumulation)[i] * secondary_convection_resistance(
            Section, T_PrimaryBulkIn, T_SecondaryBulkIn, T_SecondaryWall, SecondaryT_film, x_in, SecondarySidePressure,
            i)
        U_h.magnitude = 1 / (U_h1 + U_h2)  # [W/ cm^2 K]

        U_c1 = outer_area(Section)[i] * conduction_resistance(Section, T_PrimaryWall, SecondarySidePressure, i)
        U_c2 = outer_area(Section)[i] * primary_convection_resistance(
            Section, T_PrimaryBulkIn, T_SecondaryBulkIn, T_PrimaryWall, T_SecondaryWall, SecondarySidePressure, x_pht,
            i, m_h_leakagecorrection, TotalAccumulation) 
        U_c.magnitude = 1 / (U_c1 + U_c2)  # [W/cm^2 K]

        T_PrimaryWall = (T_PrimaryBulkIn * h_i.magnitude + T_SecondaryBulkIn * U_h.magnitude)\
            / (h_i.magnitude + U_h.magnitude)
        
        T_SecondaryWall = (T_SecondaryBulkIn * h_o.magnitude + T_PrimaryBulkIn * U_c.magnitude)\
            / (h_o.magnitude + U_c.magnitude)

        RE1 = (T_PrimaryWall - WT_h)
        RE2 = (T_SecondaryWall - WT_c)

        # if converged
        if abs(RE1) <= 0.01 and abs(RE2) <= 0.01:
            # [cm^2 K/W]
            R_F_primary.magnitude = pht_fouling_resistance(
                Section, i, calendar_year, InnerAccumulation[i], OuterAccumulation[i]
                ) 
            
            R_F_secondary.magnitude = sludge_fouling_resistance(Section, i, calendar_year)
                        
            PCR = primary_convection_resistance(
                Section, T_PrimaryBulkIn, T_SecondaryBulkIn, T_PrimaryWall, T_SecondaryWall, SecondarySidePressure,
                x_pht, i, m_h_leakagecorrection, TotalAccumulation)
            SCR = secondary_convection_resistance(
                Section, T_PrimaryBulkIn, T_SecondaryBulkIn, T_SecondaryWall, SecondaryT_film, x_in,
                SecondarySidePressure, i)
#             if i ==0:
#                 print (PCR *100*100, "PCR", conduction_resistance(Section, T_PrimaryWall, SecondarySidePressure, i)*100*100,\
#                         "cond", SCR*100*100, "SCR", R_F_primary.magnitude*100*100, "primary fouling", \
#                         R_F_secondary.magnitude*100*100, "secondary fouling", i)
            
            # [cm^2 K/W] all resistances (convective, conductive, and fouling)
            inverseU_total = (PCR + conduction_resistance(Section, T_PrimaryWall, SecondarySidePressure, i) + SCR) \
            * outer_area(Section)[i] + R_F_primary.magnitude + R_F_secondary.magnitude

            U_total.magnitude = 1 / inverseU_total  # [W/ cm^2 K]

            return T_PrimaryWall, T_SecondaryWall, U_total.magnitude

def temperature_profile(
        Section, InnerOxide, OuterOxide, m_h_leakagecorrection, SecondarySidePressure, x_pht, calendar_year):
    
    # bulk temperatures guessed, wall temperatures and overall heat transfer coefficient calculated
    # U assumed to be constant over node's area, bulk temperatures at end of node calculated, repeat
    if SecondarySidePressure == 4.593: 
        T_sat_secondary = 258.69 + 273.15
    elif SecondarySidePressure < 4.593:
        T_sat_secondary = 257 + 273.15
    
    PrimaryWall = []
    PrimaryBulk = []
    SecondaryBulk = []
    SecondaryWall = []

    for i in range(Section.NodeNumber):
        if i == 0:
            # Temperatures entering SG (0 m of SG)--> not first "node" temp (several m into SG hot leg)
            T_PrimaryBulkIn = T_sat_primary  # [K]
            T_SecondaryBulkIn = T_sat_secondary
            x_in = 0
        
        Cp_h = nc.HeatCapacity("PHT", T_PrimaryBulkIn, None)
        Cp_c = nc.HeatCapacity("SHT", T_SecondaryBulkIn, SecondarySidePressure)
        
        T_wh, T_wc, U = wall_temperature(
            Section, i, T_PrimaryBulkIn, T_SecondaryBulkIn, x_in, x_pht, InnerOxide, OuterOxide, calendar_year,
            SecondarySidePressure, m_h_leakagecorrection)
        
        # [W/ cm^2 K] * [K] * [cm^2] = [W] = [J/s]
        Q = U * (T_PrimaryBulkIn - T_SecondaryBulkIn) * outer_area(Section)[i] * TotalSGTubeNumber
        # all nodes other than preheater
        
        if Section.Length.label[i] != "preheater start" and Section.Length.label[i] != "preheater" and \
            Section.Length.label[i] != "thermal plate":
            
            # for one tube --> multiply by NumberTubes = for all tubes
            # U based on one tube (heat transfer area x number tubes) --> would cancel out in U calc (1/h*NA) NA
            # [W/ cm^2 K] * [K] * [cm^2] = [W] = [J/s]
                        
            # Ti = Ti+1 = Tsat
            T_SecondaryBulkOut = T_sat_secondary

            if x_pht > 0:
                DeltaH = EnthalpySaturatedSteam.magnitude - nc.enthalpy("PHT", T_sat_primary, None)
                
                LatentHeatAvailable = x_pht * DeltaH * m_h_leakagecorrection
                # if x_pht = 0 (entire quality entering ith section) is not enough to transfer all of heat need to SHT
                if LatentHeatAvailable < Q:
                    # this much heat needs to be transferred from PHT fluid as sensible heat
                    Q_left = Q - LatentHeatAvailable
                    
                    T_PrimaryBulkOut = T_PrimaryBulkIn - Q_left / (Cp_h * m_h_leakagecorrection * (1 - x_pht))
                    # no more latent heat/quality remains
                    x_pht = 0
                else:    
                    # latent heat being used for heat transfer, T doesn't change
                    T_PrimaryBulkOut = T_sat_primary
                    x_pht = x_pht - Q / (m_h_leakagecorrection * DeltaH)
    
                    if x_pht < 0:
                        x_pht = 0

            else:
                # no latent heat available, only sensible heat used
                T_PrimaryBulkOut = T_PrimaryBulkIn - Q / (Cp_h * m_h_leakagecorrection)
#             Q = U * (T_PrimaryBulkOut - T_SecondaryBulkOut) * outer_area(Section)[i] * Section.TubeNumber
#             x_out = (x_in + Q / (
#                 MassFlow_c_total.magnitude * EnthalpySaturatedSteam.magnitude  * (Section.TubeNumber / TotalSGTubeNumber)
#                 ))               

            if Section.Length.label[i] == "opposite preheater":
                MassFlow_c = MassFlow_downcomer.magnitude
            else:
                MassFlow_c = MassFlow_c_total.magnitude
            
            x_in = sht_steam_quality(Q, T_sat_secondary, x_in, MassFlow_c, SecondarySidePressure)

            T_PrimaryBulkIn = T_PrimaryBulkOut
            T_SecondaryBulkIn = T_SecondaryBulkOut

            PrimaryBulk.append(T_PrimaryBulkOut)
            SecondaryBulk.append(T_SecondaryBulkOut)
            PrimaryWall.append(T_wh)
            SecondaryWall.append(T_wc)
            
#             if i > 11:
#                 x_in = (x_in * MassFlow_downcomer.magnitude) / MassFlow_c_total.magnitude
           
        # preheater sections 
        else:
            # just upstream of preheater
            if Section.Length.label[i] == "preheater start":
                Cp_h = nc.HeatCapacity("PHT", T_PrimaryBulkIn, None)
                Cp_c = nc.HeatCapacity("SHT", T_PreheaterIn, SecondarySidePressure)
                
                T_PrimaryBulkOut = T_PrimaryBulkIn - Q / (Cp_h * m_h_leakagecorrection)
                T_PrimaryBulkIn = T_PrimaryBulkOut
                
                C_min = Cp_c * MassFlow_preheater.magnitude  # [J/g K]*[g/s] = [J/Ks] ] [W/K]
                C_max = Cp_h * m_h_leakagecorrection
                C_r = C_min / C_max
                
                Q_max = C_min * (T_PrimaryBulkIn - T_PreheaterIn)
                
                # total area of only preheater area
                TotalArea = sum(outer_area(Section)[(i + 1):(Section.NodeNumber - 1)])
                NTU = U * TotalArea * TotalSGTubeNumber / C_min

#                 eta = (1 - np.exp(-NTU * (1 - C_r))) / (1 - C_r * np.exp(-NTU * (1 - C_r)))
                eta = 1 - np.exp(((NTU ** 0.28) / C_r) * (np.exp(-C_r * (NTU ** 0.78)) - 1))
                Q_NTU = eta * Q_max

                # 2 known temperatures at opposite ends (preheater in and PHT in at the top)
                # solving for other two ends (PHT out, and end of preheater before mixing with recirculating downcomer)
                T_PrimaryBulkOutEnd = T_PrimaryBulkIn - Q_NTU / (m_h_leakagecorrection * Cp_h)
                T_SecondaryBulkOutEnd = T_PreheaterIn + Q_NTU / (MassFlow_preheater.magnitude * Cp_c)
                # at end of preheater - about to mix with downcomer flow
                T_SecondaryBulkOut = T_sat_secondary#T_SecondaryBulkOutEnd
                T_SecondaryBulkIn = T_sat_secondary#T_SecondaryBulkOut
                
#                 print (T_SecondaryBulkOutEnd-273.15, i)

#                 x_in = (x_in * MassFlow_downcomer.magnitude) / MassFlow_c_total.magnitude
            # inside preheater
            else:
                x_in = 0
  
                # Guessing cold-side temperatures for remaining nodes inside preheater
                T_SecondaryBulkOut = T_SecondaryBulkIn - (T_sat_secondary - T_PreheaterIn) / 4
                
           
                T_PrimaryBulkOut = T_PrimaryBulkIn - Q / (Cp_h * m_h_leakagecorrection)
#                 print (i, T_SecondaryBulkIn-273.15, T_SecondaryBulkOut-273.15)
     
                # end of preheater 
                if i == Section.NodeNumber - 2:
                    T_SecondaryBulkOut = T_PreheaterIn
                    T_PrimaryBulkOut = T_PrimaryBulkOutEnd
                    # recirculating downcomer flow entering area under thermal plate
                    T_SecondaryBulkIn = T_sat_secondary

                elif i == Section.NodeNumber - 1:
                    T_SecondaryBulkIn = T_sat_secondary
                    T_SecondaryBulkOut = T_sat_secondary
                          
                else:
                    T_SecondaryBulkIn = T_SecondaryBulkOut
                   
                T_PrimaryBulkIn = T_PrimaryBulkOut
                        
            PrimaryBulk.append(T_PrimaryBulkOut)
            SecondaryBulk.append(T_SecondaryBulkOut)
            PrimaryWall.append(T_wh)
            SecondaryWall.append(T_wc)

#     print (calendar_year, T_sat_secondary-273.15)
#     print ([j - 273.15 for j in PrimaryBulk])
#     print ([j-273.15 for j in PrimaryWall])
#     print ()
#     print ([j-273.15 for j in SecondaryBulk])
#     print ([j-273.15 for j in SecondaryWall])
#     print()
    return PrimaryBulk

def pht_massflow():
    None

def station_events(calendar_year, x_pht):
    # Divider plate leakage, mechanical cleaning, and pressure changes 

    # Pressure changes only:
    if calendar_year < 1993:
        SecondarySidePressure = 4.593  # MPa
    
    # PLNGS pressure reduction in 1992 (september) by 125 kPa
    elif 1993 <= calendar_year < 1996.25:
        SecondarySidePressure = 4.593 - (125 / 1000)  # MPa
    
    # return to full boiler secondary side pressure, 4.593 MPa
    # pressure restored shortly after reactor back online from refurb.
    elif 1996.25 <= calendar_year < 1999:
        SecondarySidePressure = 4.593
    
    elif calendar_year >= 1999:
        SecondarySidePressure = 4.593 - (125 / 1000)  # MPa
    else:
        None
        
    # divider plate raplacement in 1995, assumed to stop increase in leak (2% constant going forward)
    if calendar_year < 1996:
        # divider plate leakage rates estimated based on AECL work
        # InitialLeakage = 0.035 # fraction of total SG inlet mass flow
        # YearlyRateLeakage = 0.0065 # yearly increase to fraction of total SG inlet mass flow
        InitialLeakage = 0.02
        YearlyRateLeakage = 0.005  # yearly increase to fraction of total SG inlet mass flow
         
    elif calendar_year == 1996:
        InitialLeakage = 0.015  # controls where first post-outage point is
        YearlyRateLeakage = 0
    elif calendar_year > 1996:
        InitialLeakage = 0.037  # helps second post-outage point rise
        YearlyRateLeakage = 0.00035  # either this or some SHT deposits help out Phase 4 from dipping so much
    else:
        None 
    
    if calendar_year < 1996:
        Leakage = InitialLeakage + (calendar_year - YearStartup) * YearlyRateLeakage
    else:
        Leakage = InitialLeakage + (calendar_year - YearSHTChemicalClean) * YearlyRateLeakage
        
    DividerPlateMassFlow = MassFlow_h.magnitude * Leakage
    # decreases as divider (bypass) flow increases
    m_h_leakagecorrection = MassFlow_h.magnitude - DividerPlateMassFlow

    return SecondarySidePressure, m_h_leakagecorrection, DividerPlateMassFlow


def energy_balance(
        SteamGeneratorOutputNode, InnerOxide, OuterOxide, CleanedInnerOxide, CleanedOuterOxide, x_pht, j, SGFastMode
        ):
    year = (j * nc.TIME_STEP / 8760) 
    calendar_year = year + YearStartup
    Energy = []
    
    [SecondarySidePressure, RemainingPHTMassFlow, MasssFlow_dividerplate.magnitude] = station_events(
        calendar_year, x_pht
        )
    
    for Zone in ld.SteamGenerator:
        if SGFastMode == "yes":
            if Zone in selected_tubes:
                # tracks oxide growth for these tubes specifically
                InnerOx = Zone.InnerOxThickness
                OuterOx = Zone.OuterOxThickness
  
            else:  # assumes same growth as in default passed tube for remaining tubes
                # pass through default cleaned and not cleaned tubes
                if Zone in Cleaned:   
                    InnerOx = CleanedInnerOxide
                    OuterOx = CleanedOuterOxide
                  
                else:
                    InnerOx = InnerOxide
                    OuterOx = OuterOxide
        else:
            InnerOx = Zone.InnerOxThickness
            OuterOx = Zone.OuterOxThickness
                       
         
        # [g/cm^2] / [g/cm^3] = [cm]
#         TotalIDDeposit = [(x + y) / nc.Fe3O4Density for x, y in zip (InnerOx, OuterOx)]
#         # insert adjusted mass flow here - eventually will do with pressure 
#         AverageDeposit = sum(TotalIDDeposit) / Zone.NodeNumber
#         RemainingPHTMassFlow_fouling = (
#             RemainingPHTMassFlow * ((Zone.Diameter[0] - 2 * AverageDeposit) ** 2) / (Zone.Diameter[0] ** 2))
        
        Zone.PrimaryBulkTemperature = temperature_profile(
            Zone, InnerOx, OuterOx, RemainingPHTMassFlow, SecondarySidePressure, x_pht, calendar_year
            ) 
        
        m_timesH = (Zone.TubeNumber / TotalSGTubeNumber) * RemainingPHTMassFlow \
            * nc.enthalpy("PHT", Zone.PrimaryBulkTemperature[SteamGeneratorOutputNode], None)

        Energy.append(m_timesH)
        
        Enthalpy_dp_sat_liq = nc.enthalpy("PHT", T_sat_primary, None)
        Enthalpy_dividerplate = Enthalpy_dp_sat_liq + x_pht * (EnthalpySaturatedSteam.magnitude - Enthalpy_dp_sat_liq)
        
        Enthalpy = (sum(Energy) + MasssFlow_dividerplate.magnitude * Enthalpy_dividerplate) / MassFlow_h.magnitude 

    RIHT = nc.temperature_from_enthalpy("PHT", Enthalpy, None)
#     print (calendar_year, x_pht, RIHT-273.15)
    return RIHT
    
UncleanedInner = ld.SteamGenerator[12].InnerOxThickness
UncleanedOuter = ld.SteamGenerator[12].InnerOxThickness
CleanedInner = [i * 0.67 for i in UncleanedInner]
CleanedOuter = [i * 0.67 for i in UncleanedOuter]
#       
print (energy_balance(21, UncleanedInner, UncleanedOuter, CleanedInner, CleanedOuter, 0.002, 876 * 0, SGFastMode="yes")
- 273.15)