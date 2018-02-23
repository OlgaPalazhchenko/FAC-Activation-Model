import lepreau_data as ld
import numpy as np
import thermochemistry_and_constants as nc
import random


NumberPluggedTubes = 8
TotalSGTubeNumber = 3550 - NumberPluggedTubes

YearStartup = 1983
YearCPP = 1987
YearSHTChemicalClean = 1995.5

# select all the SG tubes to be run (by arc length)
UBends = [1.52]#, 1.52, 2.31, 3.09]
UBends = [i * 100 for i in UBends]
TubeLengths = [1980]

h_i = nc.SGParameters()
h_o = nc.SGParameters()
U_h = nc.SGParameters()
U_c = nc.SGParameters()
U_total = nc.SGParameters()
R_F_primary = nc.SGParameters()
R_F_secondary = nc.SGParameters()
MassFlow_preheater = nc.SGParameters()
MassFlow_c = nc.SGParameters()
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

# T_sat_secondary = 260.1 + 273.15
T_sat_primary = 310 + 273.15
T_PreheaterIn = 187 + 273.15
RecirculationRatio = 5.3

MassFlow_h.magnitude = 1900 * 1000
# Steam flow for 4 steam generators in typical CANDU-6 = 1033.0 kg/s
# 240 kg/s pulled from AECL COG document and works well with the 1900 kg/s hot-side flow ?
MassFlow_preheater.magnitude = 258 * 1000
MassFlow_ReheaterDrains = 0#89 * 1000 / 4
MassFlow_c.magnitude = MassFlow_preheater.magnitude * RecirculationRatio + MassFlow_ReheaterDrains

ShellDiameter.magnitude = 2.311 * 100

for i in [MassFlow_c, MassFlow_h, MasssFlow_dividerplate, MassFlow_preheater]:
    i.unit = "g/s"

for i in [ShellDiameter, EquivalentDiameter, TubePitch]:
    i.unit = "cm"

ShellCrossSectionalArea = (np.pi / 4) * (
    (ShellDiameter.magnitude ** 2) - (ld.SteamGenerator[0].OuterDiameter[0]**2) * TotalSGTubeNumber
    )

MassFlux_c.magnitude = MassFlow_c.magnitude / ShellCrossSectionalArea

for i in [MassFlux_c, MassFlux_h]:
    i.unit = "g/cm^2 s"
    
EnthalpySaturatedSteam.magnitude = 2727.5 # 9.89 MPa
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
        if sum(NumberTubes) < (0.6 * TotalSGTubeNumber):
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
        difference.append(abs(UbendLength-i))
    
    # returns index of value that has smallest difference with input Number
    return difference.index(min(difference))


def tube_picker(method):
    tubes = []
    tube_number = []
    
    # selects class initializations based on desired u-bend tube arc lengths
    if method == "arc length": 
        for k in UBends: # want multiple u-bends
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

selected_tubes = tube_picker("tube length")[0]
tube_number = tube_picker("tube length")[1] 
Cleaned = cleaned_tubes()


def thermal_conductivity(Twall, material, SecondarySidePressure):
    if material == "Alloy-800" or material == "Alloy800" or material == "A800" or material == "Alloy 800":
        Twall_C = Twall - 273.15 # oC in alloy thermal conductivity equatin 
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
    
    TubeGrowth = 0.0026 # [g/cm^2] /year = 6.5 um /year
    ReducedTubeGrowth = 0.00156 # [g/cm^2] /year = 3.25 um/year
    
    if i == 0 or i == 21:
        SludgePileGrowth = 4 * TubeGrowth
    
    elif 0 < i < 21:
        SludgePileGrowth = 0
        
    else:
        None
    # between 1983 and 1986 (included)
    if YearStartup <= calendar_year <= YearCPP:
        # secondary side fouling slope based on 1/2 that for average primary side cold leg deposit
        Accumulation = (calendar_year - YearStartup) * (TubeGrowth + SludgePileGrowth)
     
    # CPP installation (late 1986) reduces secondary side crud by 50% 
    # (assumed proportional red. in deposit formation)
    # between 1987 and 1995, not including, maybe some form of dissolution event to historic deposits
    elif YearCPP < calendar_year < YearSHTChemicalClean: 
        Accumulation = ((YearCPP - YearStartup) * (TubeGrowth + SludgePileGrowth) * (0.6)
                        + (calendar_year - YearCPP) * (ReducedTubeGrowth + SludgePileGrowth * (1 - 0.7)))
    # tubes + sludge cleaned 70%, along with reduction in sludge growth rate (by 70%)
    # after and including 1995

    elif calendar_year >= YearSHTChemicalClean:
        Accumulation = ((calendar_year - YearSHTChemicalClean) * ReducedTubeGrowth
                        + (calendar_year - YearSHTChemicalClean) * SludgePileGrowth * (1 - 0.7)
                        + (YearSHTChemicalClean - YearStartup) * (SludgePileGrowth + ReducedTubeGrowth) * (1- 0.7))
        
    Thickness = Accumulation / nc.Fe3O4Density
    
    Fouling = Thickness / thermal_conductivity(None, "outer magnetite", None)
    
    return Fouling
    

def pht_fouling_resistance(Section, i, calendar_year, InnerAccumulation, OuterAccumulation):

    # [g/cm^2]/[g/cm^3] = [cm]
    # thickness/thermal conductivity [cm]/[W/cm K] = [cm^2 K/W]
    InnerThickness = InnerAccumulation / nc.Fe3O4Density
    OuterThickness = OuterAccumulation / nc.Fe3O4Density 

    # [cm]/ [W/cm K] =[cm^2 K/W]
    # inner deposit is consolidated, predicted to have different thermal resistance than unconsolidated outer ox.
    InnerFouling = InnerThickness / thermal_conductivity(None, "inner magnetite", None)
    OuterFouling = OuterThickness / thermal_conductivity(None, "outer magnetite", None)
    # total pht thermal resistance
    
    return InnerFouling + OuterFouling # [cm^2 K/W]


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


def primary_convection_resistance(Section, correlation, T_film, T_wall, SecondarySidePressure, x_pht, i):
    # R_1,Convective = 1/(h_1,convectove *A)
    # A = inner area (based on inner diameter)
    # [W/cm K]/ [cm] = [W/K cm^2]
    
    if Section.Length.label[i] == "PHT boiling":
        MassFlux_h.magnitude = MassFlow_h.magnitude / (TotalSGTubeNumber * (np.pi / 4) * (Section.Diameter[i] ** 2))
        
        h_i = boiling_heat_transfer(
            x_pht, "PHT", T_sat_primary, MassFlux_h.magnitude, T_wall, Section.Diameter[i], SecondarySidePressure, i
            )
         
    else:
        h_i = nusseltnumber(Section, correlation, T_film, SecondarySidePressure, i, x_pht) \
        * thermal_conductivity(T_film, "water", SecondarySidePressure) / Section.Diameter[i]

    return 1 / (h_i * inner_area(Section)[i])  # [K/W]


def secondary_convection_resistance(Section, T_film, T_wall, x_in, SecondarySidePressure, i):
    if SecondarySidePressure == 4.593: 
        T_sat_secondary = 258.69 + 273.15
    elif SecondarySidePressure < 4.593:
        T_sat_secondary = 257 + 273.15
    
    # split into boiling and non-boiling (preheater) sections
    if Section.Length.label[i] == "preheater" or Section.Length.label[i] == "preheater start":  # from Silpsrikul thesis
        # Based on PLGS data and Silpsrikul calculations
        f_b = 0.2119  # fraction of cross-sectional area of shell occupied by a baffle window
        N_b = f_b * TotalSGTubeNumber  # number of tubes in baffle window
        D_s = ((4 * 0.5 * np.pi / 4) * (ShellDiameter.magnitude / 100) ** 2) \
            / (0.5 * np.pi * (ShellDiameter.magnitude / 100) + ShellDiameter.magnitude / 100)  # [m]

        S_b = (f_b * np.pi * (D_s ** 2) / 4) - N_b * np.pi * ((Section.OuterDiameter[i] / 100) ** 2) / 4  # [cm^2]
        G_b = (MassFlow_preheater.magnitude / 1000) / S_b  # [kg/m^2s]
        S_c = 0.5373 * D_s * (1 - (Section.OuterDiameter[i] / 100) / 0.0219)
        G_c = MassFlow_preheater.magnitude / 1000 / S_c
        G_e = np.sqrt(G_b * G_c) * 1000 / (100 ** 2) #[g/cm^2]

        # First two (raised to exponent) terms are unitless
        # [W/cm K]/[cm] = [W/cm^2 K]
        k = thermal_conductivity(T_film, "water", SecondarySidePressure)
        h_o = (k * 0.2 / Section.OuterDiameter[i]) \
            * ((Section.OuterDiameter[i] * G_e / nc.viscosity("SHT", T_film, SecondarySidePressure)) ** 0.6) \
            * (
                nc.HeatCapacity("SHT", T_film, SecondarySidePressure) \
                * nc.viscosity("SHT", T_film, SecondarySidePressure) \
                / thermal_conductivity(T_film, "water", SecondarySidePressure)
                ) ** 0.33

    else:
        # [cm]
        EquivalentDiameter.magnitude = 4 * \
        ((TubePitch.magnitude ** 2) - (np.pi * (Section.OuterDiameter[i] ** 2)) / 4) \
        / (np.pi * Section.OuterDiameter[i])
        
        x = x_in #* 100
#         if i <= 12:
#             x = x_in #* 100
#         else:
#             x = .08

        h_o = boiling_heat_transfer(
            x, "SHT", T_sat_secondary, MassFlux_c.magnitude, T_wall, EquivalentDiameter.magnitude, 
            SecondarySidePressure, i)

    return 1 / (h_o * outer_area(Section)[i])  # K/W


def boiling_heat_transfer(x, side, T_sat, MassFlux, T_wall, Diameter, SecondarySidePressure, i):
    # x = steam quality in percent
      
    if side == "SHT":
        rho_v = 1000 * 23.187 / (100 ** 3)  # [g/cm^3]
        Density = nc.density(side, T_sat, SecondarySidePressure)
        Viscosity = nc.viscosity(side, T_sat, SecondarySidePressure)
        Pressure = SecondarySidePressure
    else:
        rho_v = 1000 * 55.462/ (100 ** 3)  # [g/cm^3]
        Density = nc.D2O_density(T_sat)
        Viscosity = nc.D2O_viscosity(T_sat)
        Pressure = nc.PrimarySidePressure
    
    p_crit = 22.0640  # [MPa]
   
    F = (1 + (x * prandtl(side, T_sat, Pressure, i) * ((Density / rho_v) - 1))) ** 0.35

    MassFlux_liquid = MassFlux * (1 - x) 

    Re_D = Diameter * MassFlux_liquid / Viscosity

    # 4*MassFlow/((nc.viscosity("SHT", T_sat_secondary)/1000)*np.pi*Shell_ID)
    h_l = 0.023 * nc.thermal_conductivityH2O(side, T_sat, Pressure) \
    * ((Re_D) ** 0.8) * (prandtl(side, T_sat, Pressure, i) ** 0.4) / Diameter  # [W/cm^2 K]
    
    Q_prime = F * h_l * (abs(T_wall - T_sat))  # [W/cm^2]
    
    S = (1 + 0.055 * (F ** 0.1) * (Re_D) ** 0.16) ** (-1)
    A_p = (55 * ((4.70 / p_crit) ** 0.12) * ((-np.log10(4.70 / p_crit)) ** (-0.55)) * (nc.H2MolarMass) ** (-0.5)) \
            / (100 ** 2)
    
    C = ((A_p * S) / (F * h_l) ** 2) * (Q_prime) ** (4 / 3)
    coeff = [1, -C, -1]
    cubic_solution = np.roots(coeff)
    q = cubic_solution[0]
    
    HeatTransferCoefficient = F * (q ** (3 / 2)) * h_l  # [W/cm^2 K]
    
    return HeatTransferCoefficient


def nusseltnumber(Section, correlation, Temperature, SecondarySidePressure, i, x_pht):
    # DensityH2O = nc.density("PHT", Temperature)

    # Re_D = Section.Velocity[i]*Section.Diameter[i]/(ViscosityH2O*DensityH2O)
    MassFlux_h.magnitude = MassFlow_h.magnitude * (1 - x_pht) \
    / (TotalSGTubeNumber * (np.pi / 4) * (Section.Diameter[i] ** 2))
    
    #only primary side heat transfer correlation uses Nu number
    Re_D = [(MassFlux_h.magnitude / nc.D2O_viscosity(Temperature)) * i for i in Section.Diameter]

    if correlation == "Dittus-Boelter":
        n = 1 / 3
    elif correlation == "Colburn":
        n = 0.4
    else:
        None
    C = 0.023
    m = 4 / 5

    Nu = C * (Re_D[i] ** m) * ((prandtl("PHT", Temperature, SecondarySidePressure, i)) ** n)
    return Nu


def prandtl(side, Temperature, SecondarySidePressure, i):
    if side == "PHT":
        Viscosity = nc.D2O_viscosity(Temperature)
    else:
        Viscosity = nc.viscosity("SHT", Temperature, SecondarySidePressure)
        
    Pr = nc.HeatCapacity("PHT", Temperature, SecondarySidePressure) \
    * Viscosity / nc.thermal_conductivityH2O("PHT", Temperature, SecondarySidePressure)
    return Pr


def inner_area(Section):
    return [np.pi * x * y for x, y in zip(Section.Diameter, Section.Length.magnitude)]


def outer_area(Section):
    return [np.pi * x * y for x, y in zip(Section.OuterDiameter, Section.Length.magnitude)]


def pht_steam_quality(Temperature):
    CoreMassFlow = (MassFlow_h.magnitude / 1000) * 4 # [kg /s]
    Delta_T = T_sat_primary - (262.5 + 273.15) # [K]
    C_p_cold = nc.HeatCapacity("PHT", Temperature, SecondarySidePressure = None)
    C_p_hot = nc.HeatCapacity("PHT", T_sat_primary, SecondarySidePressure = None)
    C_p_avg = (C_p_cold + C_p_hot) / 2 # [kJ/kg K]
    Power = CoreMassFlow * C_p_avg * Delta_T # [kW]
    
    # core inlet enthalpy at RIHT + that added from fuel
    # H_fromfuel = Power / CoreMassFlow # [kJ/s /kg/s] = [kJ/kg]
    H_pht = Power / CoreMassFlow
    H_satliq_outlet = nc.enthalpy("PHT", T_sat_primary, None)
    
    # no quality in return primary flow
    H_current = nc.enthalpy("PHT", Temperature, None) + H_pht
    x = (H_current - H_satliq_outlet) / (EnthalpySaturatedSteam.magnitude - H_satliq_outlet)
    return x


def sht_steam_quality(Q, T_sat_secondary, x):
    
    # [J/s] / [g/s] = [J /g]
    H_pht = Q / (MassFlow_c.magnitude) # [kJ/kg]
    H_satliq = nc.enthalpy("PHT", T_sat_secondary, None)
    H_SaturatedSteam = 2797.3 # [kJ/kg]
    
    H_prev = H_satliq + x * (H_SaturatedSteam - H_satliq)
    H_current = H_pht + H_prev

    x = (H_current - H_satliq) / (H_SaturatedSteam - H_satliq)
    return x
    

def wall_temperature(
        Section, i, T_PrimaryBulkIn, T_SecondaryBulkIn, x_in, x_pht, InnerAccumulation, OuterAccumulation,
        calendar_year, SecondarySidePressure):
    
    # i = each node of SG tube
    T_PrimaryWall = T_PrimaryBulkIn - (1 / 3) * (T_PrimaryBulkIn - T_SecondaryBulkIn)
    T_SecondaryWall = T_PrimaryBulkIn - (1 / 3) * (T_PrimaryBulkIn - T_SecondaryBulkIn)

    for k in range(50):
        WT_h = T_PrimaryWall
        WT_c = T_SecondaryWall

        PrimaryT_film = (T_PrimaryBulkIn + T_PrimaryWall) / 2
        SecondaryT_film = (T_SecondaryBulkIn + T_SecondaryWall) / 2

        h_i.magnitude = 1 / (primary_convection_resistance(
            Section, "Dittus-Boelter", PrimaryT_film, T_PrimaryWall, SecondarySidePressure, x_pht, i
            ) * inner_area(Section)[i])
        
        h_o.magnitude = 1 / (secondary_convection_resistance(
            Section, SecondaryT_film, T_SecondaryWall, x_in, SecondarySidePressure, i
            ) * outer_area(Section)[i])
        # R = 1/hA --> h=A*1/R

        U_h1 = inner_area(Section)[i] * conduction_resistance(Section, T_PrimaryWall, SecondarySidePressure, i)
        
        U_h2 = inner_area(Section)[i] * secondary_convection_resistance(
            Section, SecondaryT_film, T_SecondaryWall, x_in, SecondarySidePressure, i
            )
        U_h.magnitude = 1 / (U_h1 + U_h2)  # [W/ cm^2 K]

        U_c1 = outer_area(Section)[i] * conduction_resistance(Section, T_PrimaryWall, SecondarySidePressure, i)
        U_c2 = outer_area(Section)[i] * primary_convection_resistance(
            Section, "Dittus-Boelter", PrimaryT_film, T_PrimaryWall, SecondarySidePressure, x_pht, i
            ) 
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
                Section, "Dittus-Boelter", PrimaryT_film, T_PrimaryWall, SecondarySidePressure, x_pht, i
                )
            SCR = secondary_convection_resistance(
                Section, SecondaryT_film, T_SecondaryWall, x_in, SecondarySidePressure, i
                )
# #             if calendar_year == 1983 or calendar_year == 1984:
#             print (PCR *100*100, "PCR", conduction_resistance(Section, T_PrimaryWall, SecondarySidePressure, i)*100*100,\
#                     "cond", SCR*100*100, "SCR", R_F_primary.magnitude*100*100, "primary fouling", \
#                     R_F_secondary.magnitude*100*100, "secondary fouling", i)
            
            #[cm^2 K/W] all resistances (convective, conductive, and fouling)
            inverseU_total = (PCR + conduction_resistance(Section, T_PrimaryWall, SecondarySidePressure, i) + SCR) \
            * outer_area(Section)[i] + R_F_primary.magnitude + R_F_secondary.magnitude

            U_total.magnitude = 1 / inverseU_total  # [W/ cm^2 K]

            return T_PrimaryWall, T_SecondaryWall, U_total.magnitude


def temperature_profile(
        Section, InnerAccumulation, OuterAccumulation, m_h_leakagecorrection, SecondarySidePressure, T_primary_in,
        x_pht, calendar_year):
    
    # bulk temperatures guessed, wall temperatures and overall heat transfer coefficient calculated
    # U assumed to be constant over node's area, bulk temperatures at end of node calculated, repeat
    if SecondarySidePressure == 4.593: 
        T_sat_secondary = 258.69 + 273.15 # 258.69
    elif SecondarySidePressure < 4.593:
        T_sat_secondary = 257 + 273.15  # 257
    
    PrimaryWall = []
    PrimaryBulk = []
    SecondaryBulk = []
    SecondaryWall = []

    for i in range(Section.NodeNumber - 1):
        if i == 0:
            # Temperatures entering SG (0 m of SG)--> not first "node" temp (several m into SG hot leg)
            T_PrimaryBulkIn = T_primary_in  # [K]
            T_SecondaryBulkIn = T_sat_secondary
            x_in = 0

        Cp_h = nc.HeatCapacity("PHT", T_PrimaryBulkIn, SecondarySidePressure)
        Cp_c = nc.HeatCapacity("SHT", T_SecondaryBulkIn, SecondarySidePressure)
        T_wh, T_wc, U = wall_temperature(
            Section, i, T_PrimaryBulkIn, T_SecondaryBulkIn, x_in, x_pht, InnerAccumulation, OuterAccumulation,
            calendar_year, SecondarySidePressure)
        
        # all nodes other than preheater
        if Section.Length.label[i] != "preheater start" and Section.Length.label[i] != "preheater" and \
                Section.Length.label[i] != "thermal plate":

            T_PrimaryBulkOut = T_PrimaryBulkIn - (U * (T_PrimaryBulkIn - T_SecondaryBulkIn) * outer_area(Section)[i] \
                                                  * TotalSGTubeNumber) / (Cp_h * m_h_leakagecorrection)
            # Ti = Ti+1 = Tsat
            T_SecondaryBulkOut = T_sat_secondary

            # for one tube --> multiply by NumberTubes = for all tubes
            # Q = m_h_leakagecorrection*Cp_h*(T_PrimaryBulkIn-T_PrimaryBulkOut)

            # U based on one tube (heat transfer area x number tubes) --> would cancel out in U calc (1/h*NA) NA
            
#             Q = U * (T_PrimaryBulkOut - T_SecondaryBulkOut) * outer_area(Section)[i] * Section.TubeNumber
# 
#             x_out = (x_in + Q / (
#                 MassFlow_c.magnitude * EnthalpySaturatedSteam.magnitude  * (Section.TubeNumber / TotalSGTubeNumber)
#                 ))
   
            # [W/ cm^2 K] * [K] * [cm^2] = [W] = [J/s]
            Q = U * (T_PrimaryBulkOut - T_SecondaryBulkOut) * outer_area(Section)[i] * TotalSGTubeNumber
 
#             x_out = x_in + Q / (MassFlow_c.magnitude * EnthalpySaturatedSteam.magnitude)
            x_in = sht_steam_quality(Q, T_SecondaryBulkOut, x_in)
            T_PrimaryBulkIn = T_PrimaryBulkOut
            T_SecondaryBulkIn = T_SecondaryBulkOut

            PrimaryBulk.append(T_PrimaryBulkOut)
            SecondaryBulk.append(T_SecondaryBulkOut)
            PrimaryWall.append(T_wh)
            SecondaryWall.append(T_wc)
        else:
            if Section.Length.label[i] == "preheater start":
                C_min = Cp_c * MassFlow_preheater.magnitude  # [J/g K]*[g/s] = [J/Ks] ] [W/K]
                C_max = Cp_h * m_h_leakagecorrection
                C_r = C_min / C_max
                Q_max = C_min * (T_PrimaryBulkIn - T_PreheaterIn)
                TotalArea = sum(outer_area(Section)[i:(Section.NodeNumber - 1)])
                NTU = U * TotalArea * TotalSGTubeNumber / C_min

                eta = (1 - np.exp(-NTU * (1 - C_r))) / (1 - C_r * np.exp(-NTU * (1 - C_r)))
                # eta = 1-np.exp(((NTU**0.28)/C_r)*(np.exp(-C_r*(NTU**0.78))-1))
                Q_NTU = eta * Q_max

                T_PrimaryBulkOutEnd = T_PrimaryBulkIn - Q_NTU / (m_h_leakagecorrection * Cp_h)
                T_SecondaryBulkOutEnd = T_PreheaterIn + Q_NTU / (MassFlow_preheater.magnitude * Cp_c)
                T_SecondaryBulkIn = T_SecondaryBulkOutEnd

            # Guessing cold-side temperatures for remaining nodes inside preheater
            T_SecondaryBulkOut = T_SecondaryBulkOut - i
            T_PrimaryBulkOut = T_PrimaryBulkIn \
                - (U * outer_area(Section)[i] * TotalSGTubeNumber / (Cp_h * m_h_leakagecorrection)) \
                * (T_PrimaryBulkIn - T_SecondaryBulkIn)

            T_PrimaryBulkIn = T_PrimaryBulkOut
            T_SecondaryBulkIn = T_SecondaryBulkOut
            
            if i == Section.NodeNumber - 2:
                T_SecondaryBulkOut = T_PreheaterIn
                T_PrimaryBulkOut = T_PrimaryBulkOutEnd

            SecondaryBulk[16] = T_SecondaryBulkOutEnd
            PrimaryBulk.append(T_PrimaryBulkOut)
            SecondaryBulk.append(T_SecondaryBulkOut)
            PrimaryWall.append(T_wh)
            SecondaryWall.append(T_wc)
    
    # thermal plate below preheater
    PrimaryBulk.append(PrimaryBulk[20])
    PrimaryWall.append(PrimaryWall[20])

#     print ([j - 273.15 for j in PrimaryBulk])
#     print ([j-273.15 for j in PrimaryWall])
#     print ()
#     print ([j-273.15 for j in SecondaryBulk])
#     print ([j-273.15 for j in SecondaryWall])
#     print()
    return PrimaryBulk


def station_events(calendar_year, x_pht):
    # Divider plate leakage, mechanical cleaning, and pressure changes 

    # Pressure changes only:
    if calendar_year < 1992.5:
        SecondarySidePressure = 4.593   # MPa
    
    # PLNGS pressure reduction in 1992 (september) by 125 kPa
    elif 1992.5 <= calendar_year <= 1995.5:
        SecondarySidePressure = 4.593 - (125 / 1000) # MPa
    
    # return to full boiler secondary side pressure, 4.593 MPa
    # pressure restored shortly after reactor back online from refurb.
    elif 1995.5 < calendar_year < 1998.5:
        SecondarySidePressure = 4.593
    
    elif calendar_year >= 1998.5:
        SecondarySidePressure = 4.593 - (125 / 1000) # MPa
    else:
        None

    
    # divider plate raplacement in 1995, assumed to stop increase in leak (2% constant going forward)
    if calendar_year < 1995.5:
        # divider plate leakage rates estimated based on AECL work
        # InitialLeakage = 0.035 # fraction of total SG inlet mass flow
        # YearlyRateLeakage = 0.0065 # yearly increase to fraction of total SG inlet mass flow
        InitialLeakage = 0.03 
        YearlyRateLeakage = 0.0065 # yearly increase to fraction of total SG inlet mass flow
         
    elif calendar_year >= 1995.5:
        InitialLeakage = 0.02 
        YearlyRateLeakage = 0
    else:
        None
    
#     # partial mechanical clean of boiler tube primary side (67 % efficiency for the 60 % of tubes accessed)        
#     if calendar_year == 1995.5:
#         for Zone in ld.SteamGenerator:
#             if Zone in Cleaned:
#                 for i in range(Zone.NodeNumber - 1):
#                     if Zone.OuterOxThickness[i] > 0:
#                         Zone.OuterOxThickness[i] = Zone.OuterOxThickness[i] * (1 - 0.67)
#                     else:
#                         Zone.InnerOxThickness[i] = Zone.InnerOxThickness[i] * (1 - 0.67)
    MassFlow_h_liquid = MassFlow_h.magnitude * (1 - x_pht)    
    
    Leakage = InitialLeakage + (calendar_year - YearStartup) * YearlyRateLeakage
    DividerPlateMassFlow = MassFlow_h_liquid * Leakage
    # decreases as divider (bypass) flow increases
    m_h_leakagecorrection = MassFlow_h_liquid - DividerPlateMassFlow

    return SecondarySidePressure, m_h_leakagecorrection, DividerPlateMassFlow


def energy_balance(SteamGeneratorOutputNode, InnerAccumulation, OuterAccumulation, T_primary_in, x_pht, j):
    year = (j / 8760) 
    calendar_year = year + YearStartup
    
    [SecondarySidePressure, RemainingPHTMassFlow, MasssFlow_dividerplate.magnitude] = station_events(
        calendar_year, x_pht
        )
    
    Energy = []
    for Zone in ld.SteamGenerator:
        
#         # partial mechanical clean of boiler tube primary side (67 % efficiency for the 60 % of tubes accessed)
#         if calendar_year == 1995.5:        
#             if Zone in SGHX.Cleaned:
#                 for i in range(Zone.NodeNumber - 1):
#                     if Zone.OuterOxThickness[i] > 0:
#                         Zone.OuterOxThickness[i] = Zone.OuterOxThickness[i] * (1 - 0.67)
#                     else:
#                         Zone.InnerOxThickness[i] = Zone.InnerOxThickness[i] * (1 - 0.67)
#                 InnerOx = Zone.InnerOxThickness
#                 OuterOx = Zone.OuterOxThickness
    
        if Zone in selected_tubes:
            # tracks oxide growth for these tubes specifically
            InnerOx = Zone.InnerOxThickness
            OuterOx = Zone.OuterOxThickness
        else: # assumes same growth as in default passed tube for remaining tubes
            # pass through default cleaned and not cleaned tube? 
            # work cleaning into here
            InnerOx = InnerAccumulation
            OuterOx = OuterAccumulation
        
        Zone.PrimaryBulkTemperature = temperature_profile(
            Zone, InnerOx, OuterOx, RemainingPHTMassFlow, SecondarySidePressure, T_primary_in, x_pht, calendar_year
            ) 
        
        m_H = (Zone.TubeNumber / TotalSGTubeNumber) * RemainingPHTMassFlow \
            * nc.enthalpy("PHT", Zone.PrimaryBulkTemperature[SteamGeneratorOutputNode], SecondarySidePressure)
            
        MassFlow_h_liquid = MassFlow_h.magnitude * (1 - x_pht)

        Energy.append(m_H)
        Enthalpy = (
            sum(Energy) + MasssFlow_dividerplate.magnitude * nc.enthalpy("PHT", T_primary_in, SecondarySidePressure)
            ) / MassFlow_h_liquid

    RIHT = nc.temperature_from_enthalpy("PHT", Enthalpy, SecondarySidePressure)
    return RIHT
  
# print (energy_balance(21, ld.SteamGenerator[12].InnerOxThickness, ld.SteamGenerator[12].OuterOxThickness,
#                       583, 0, 0) - 273.15)