import lepreau_data as ld
import numpy as np
import constants as nc

# T_sat_secondary = 260.1 + 273.15
T_sat_primary = 310 + 273.15
T_PreheaterIn = 186.7 + 273.15
T_PrimaryIn = 310 + 273.15

h_i = ld.SGParameters()
h_o = ld.SGParameters()
U_h = ld.SGParameters()
U_c = ld.SGParameters()
U_total = ld.SGParameters()
R_F_primary = ld.SGParameters()
R_F_secondary = ld.SGParameters()
MassFlow_c = ld.SGParameters()
ShellDiameter = ld.SGParameters()
MassFlux_c = ld.SGParameters()
MassFlux_h = ld.SGParameters()
EnthalpySaturatedSteam = ld.SGParameters()
MassFlow_h = ld.SGParameters()
MasssFlow_dividerplate = ld.SGParameters()
EquivalentDiameter = ld.SGParameters()
TubePitch = ld.SGParameters()

for i in [h_i, h_o]:
    i.unit = "K/W"

for i in [U_h, U_c, U_total]:
    i.unit = "W/cm^2 K"

for i in [R_F_primary, R_F_secondary]:
    i.unit = "cm^2 K/W"
# Steam flow for 4 steam generators in typical CANDU-6 = 1033.0 kg/s
MassFlow_c.magnitude = 240 * 1000
MassFlow_h.magnitude = 1900 * 1000

ShellDiameter.magnitude = 2.28 * 100

for i in [MassFlow_c, MassFlow_h, MasssFlow_dividerplate]:
    i.unit = "g/s"

for i in [ShellDiameter, EquivalentDiameter, TubePitch]:
    i.unit = "cm"

MassFlux_c.magnitude = MassFlow_c.magnitude / ((np.pi / 4) * (ShellDiameter.magnitude ** 2))
for i in [MassFlux_c, MassFlux_h]:
    i.unit = "g/cm^2 s"
    
EnthalpySaturatedSteam.magnitude = 2800.4
EnthalpySaturatedSteam.unit = "J/g"

TubePitch.magnitude = 2.413


def thermal_conductivity(Twall, material, SecondarySidePressure):
    if material == "Alloy-800" or material == "Alloy800" or material == "A800" or material == "Alloy 800":
        return (11.450 + 0.0161 * Twall) / 100  # M.K.E thesis for Alloy-800 tubing [W/cm K]
    elif material == "water":
        return ld.thermal_conductivityH2O("PHT", Twall, SecondarySidePressure)
    elif material == "magnetite":
        return 1.4 / 100  # [W/cm K]
    else:
        print ("Error: material not specified")
        return None


def fouling_resistance(InnerAccumulation, OuterAccumulation, SecondarySidePressure):
    # thickness/thermal conductivity [cm]/[W/cm K] = [cm^2 K/W]
    # [g/cm^2]/[g/cm^3] = [cm]

    InnerThickness = [i / nc.Fe3O4Density for i in InnerAccumulation]
    OuterThickness = [i / nc.Fe3O4Density for i in OuterAccumulation]

    # [cm]/ [W/cm K] =[cm^2 K/W]
    InnerFouling = [i / 2 * thermal_conductivity(None, "magnetite", SecondarySidePressure) for i in InnerThickness]
    OuterFouling = [i / thermal_conductivity(None, "magnetite", SecondarySidePressure) for i in OuterThickness]
    return [x + y for x, y in zip(InnerFouling, OuterFouling)]


def conduction_resistance(Section, Twall, SecondarySidePressure, i):
    # R_conduction = L/kA (L = thickness)
    # A = area with shape factor correcrion factor for cylindrical pipe
    # R_conduction = ln(D_o/D_i)/(2pikl) [K/W]
    # l = length, k =thermal conductivity coefficient [W/cm K]

    k_w = thermal_conductivity(Twall, "Alloy-800", SecondarySidePressure)  # [W/cm K]

    Rcyl_numerator = np.log(Section.OuterDiameter[i] / Section.Diameter[i])
    Rcyl_denominator = 2 * np.pi * k_w * Section.Length.magnitude[i]

    R_cond = Rcyl_numerator / Rcyl_denominator  # [K/W]
    return R_cond


def primary_convection_resistance(Section, correlation, T_film, T_wall, SecondarySidePressure, i):
    # R_1,Convective = 1/(h_1,convectove *A)
    # A = inner area (based on inner diameter)
    # [W/cm K]/ [cm] = [W/K cm^2]
    
    if Section.Length.label[i] == "PHT boiling":
        x = 4.4
        MassFlux_h.magnitude = MassFlow_h.magnitude / (nc.TotalSGTubeNumber * (np.pi / 4) * (Section.Diameter[i] ** 2))
 
        h_i = boiling_heat_transfer(
            x, "PHT", T_sat_primary, MassFlux_h.magnitude, T_wall, Section.Diameter[i], SecondarySidePressure, i
            )
         
    else:
        h_i = nusseltnumber(Section, correlation, T_film, SecondarySidePressure, i) \
        * thermal_conductivity(T_film, "water", SecondarySidePressure) / Section.Diameter[i]

    return 1 / (h_i * inner_area(Section)[i])  # [K/W]


def secondary_convection_resistance(Section, T_film, T_wall, x_in, SecondarySidePressure, i):
    if SecondarySidePressure == 4.593: 
        T_sat_secondary = 260.1 + 273.15
    elif SecondarySidePressure < 4.593:
        T_sat_secondary = 258.1 + 273.15  
    
    # split into boiling and non-boiling (preheater) sections
    if Section.Length.label[i] == "preheater" or Section.Length.label[i] == "preheater start":  # from Silpsrikul thesis
        # Based on PLGS data and Silpsrikul calculations
        f_b = 0.2119  # fraction of cross-sectional area of shell occupied by a baffle window
        N_b = f_b * nc.TotalSGTubeNumber  # number of tubes in baffle window
        D_s = ((4 * 0.5 * np.pi / 4) * (ShellDiameter.magnitude / 100) ** 2) \
            / (0.5 * np.pi * (ShellDiameter.magnitude / 100) + ShellDiameter.magnitude / 100)  # [m]

        S_b = (f_b * np.pi * (D_s ** 2) / 4) - N_b * np.pi * ((Section.OuterDiameter[i] / 100) ** 2) / 4  # [cm^2]
        G_b = (MassFlow_c.magnitude / 1000) / S_b  # [kg/m^2s]
        S_c = 0.5373 * D_s * (1 - (Section.OuterDiameter[i] / 100) / 0.0219)
        G_c = MassFlow_c.magnitude / 1000 / S_c
        G_e = np.sqrt(G_b * G_c) * 1000 / (100 ** 2)

        # First two (raised to exponent) terms are unitless
        # [W/cm K]/[cm] = [W/cm^2 K]
        h_o = (thermal_conductivity(T_film, "water", SecondarySidePressure) * 0.25 / Section.OuterDiameter[i]) \
            * ((Section.OuterDiameter[i] * G_e / ld.Viscosity("water", "SHT", T_film)) ** 0.6) \
            * (
                ld.HeatCapacity("SHT", T_film, SecondarySidePressure) \
                * ld.Viscosity("water", "SHT", T_film) / thermal_conductivity(T_film, "water", SecondarySidePressure)
                ) ** 0.33

    else:
        EquivalentDiameter.magnitude = 4 * \
        ((TubePitch.magnitude ** 2) - (np.pi * (Section.OuterDiameter[i] ** 2)) / 4) \
        / (np.pi * Section.OuterDiameter[i])
        
        if i <= 9:
            x = 0.5 * x_in * 100  
        else:
            x = 24 - (i/1)
      
        h_o = boiling_heat_transfer(
            x, "SHT", T_sat_secondary, MassFlux_c.magnitude, T_wall, EquivalentDiameter.magnitude, 
            SecondarySidePressure, i)

    return 1 / (h_o * outer_area(Section)[i])  # K/W


def boiling_heat_transfer(x, side, T_sat, MassFlux, T_wall, Diameter, SecondarySidePressure, i):
    # x = steam quality in percent
      
    if side == "SHT":
        rho_v = 1000 * 23.187 / (100 ** 3)  # [g/cm^3]
    else:
        rho_v = 1000 * 55.462/ (100 ** 3)  # [g/cm^3]
    
    p_crit = 22.0640  # [MPa]
   
    F = (
        1 + (x * prandtl(T_sat, SecondarySidePressure, i) \
             * ((ld.Density("water", side, T_sat, SecondarySidePressure) / rho_v) - 1))
         ) ** 0.35
    
    Re_D = Diameter * MassFlux / ld.Viscosity("water", side, T_sat)
    
    # 4*MassFlow/((ld.Viscosity("water", "SHT", T_sat_secondary)/1000)*np.pi*Shell_ID)
    h_l = 0.023 * ld.thermal_conductivityH2O(side, T_sat, SecondarySidePressure) \
    * ((Re_D) ** 0.8) * (prandtl(T_sat, SecondarySidePressure, i) ** 0.4) / Diameter  # [W/cm^2 K]
    
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


def nusseltnumber(Section, correlation, Temperature, SecondarySidePressure, i):
    # ViscosityH2O = ld.Viscosity("water", "PHT", Temperature)
    # DensityH2O = ld.Density("water", "PHT", Temperature)

    # Re_D = Section.Velocity[i]*Section.Diameter[i]/(ViscosityH2O*DensityH2O)
    MassFlux_h.magnitude = MassFlow_h.magnitude / (nc.TotalSGTubeNumber * (np.pi / 4) * (Section.Diameter[i] ** 2))
    Re_D = [(MassFlux_h.magnitude / ld.Viscosity("water", "PHT", Temperature)) * i for i in Section.Diameter]
    # print (Re_D[i], Re_[i])

    if correlation == "Dittus-Boelter":
        n = 1 / 3
    elif correlation == "Colburn":
        n = 0.4
    else:
        print ("Error: Correlation not recognized in nusseltnumber function")
    C = 0.023
    m = 4 / 5

    Nu = C * (Re_D[i] ** m) * ((prandtl(Temperature, SecondarySidePressure, i)) ** n)
    return Nu


def prandtl(Temperature, SecondarySidePressure, i):
    # Need a better reference for Cp/viscosity data
    ViscosityH2O = ld.Viscosity("water", "PHT", Temperature)
    Pr = ld.HeatCapacity("PHT", Temperature, SecondarySidePressure) \
    * ViscosityH2O / ld.thermal_conductivityH2O("PHT", Temperature, SecondarySidePressure)
    return Pr


def inner_area(Section):
    return [np.pi * x * y for x, y in zip(Section.Diameter, Section.Length.magnitude)]


def outer_area(Section):
    return [np.pi * x * y for x, y in zip(Section.OuterDiameter, Section.Length.magnitude)]


def wall_temperature(
        Section, i, T_PrimaryBulkIn, T_SecondaryBulkIn, x_in, InnerAccumulation, OuterAccumulation, calendar_year, 
        SecondarySidePressure):
    
    # i = each node of SG tube
    T_PrimaryWall = T_PrimaryBulkIn - (1 / 3) * (T_PrimaryBulkIn - T_SecondaryBulkIn)
    T_SecondaryWall = T_PrimaryBulkIn - (1 / 3) * (T_PrimaryBulkIn - T_SecondaryBulkIn)

    for k in range(50):
        WT_h = T_PrimaryWall
        WT_c = T_SecondaryWall

        PrimaryT_film = (T_PrimaryBulkIn + T_PrimaryWall) / 2
        SecondaryT_film = (T_SecondaryBulkIn + T_SecondaryWall) / 2

        h_i.magnitude = 1 / (primary_convection_resistance(
            Section, "Dittus-Boelter", PrimaryT_film, T_PrimaryWall, SecondarySidePressure, i
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
            Section, "Dittus-Boelter", PrimaryT_film, T_PrimaryWall, SecondarySidePressure, i
            ) 
        U_c.magnitude = 1 / (U_c1 + U_c2)  # [W/cm^2 K]

        T_PrimaryWall = (T_PrimaryBulkIn * h_i.magnitude + T_SecondaryBulkIn * U_h.magnitude)\
            / (h_i.magnitude + U_h.magnitude)
        
        T_SecondaryWall = (T_SecondaryBulkIn * h_o.magnitude + T_PrimaryBulkIn * U_c.magnitude)\
            / (h_o.magnitude + U_c.magnitude)

        RE1 = (T_PrimaryWall - WT_h)
        RE2 = (T_SecondaryWall - WT_c)

        if abs(RE1) <= 0.01 and abs(RE2) <= 0.01:
            # [cm^2 K/W]
            R_F_primary.magnitude = fouling_resistance(InnerAccumulation, OuterAccumulation, SecondarySidePressure)[i]  
            
            Initial_SHTSFouling = 0.0005 # [g/cm2] same as for PHT (~ 1 um/a)
            
            if 1983 <= calendar_year < 1987:
                SHTSFouling = Initial_SHTSFouling + (calendar_year - 1983) * 0.0024
                # secondary side fouling slope based on 1/2 that for average primary side cold leg deposit
            else:
                # CPP installation (late 1986) reduces secondary side crud by 50% 
                # (assumed proportional red. in deposit formation)
                SHTSFouling = Initial_SHTSFouling + (calendar_year - 1983) * 0.0012
                
            R_F_secondary.magnitude = fouling_resistance(
                [0]*Section.NodeNumber, [SHTSFouling]*Section.NodeNumber, SecondarySidePressure
                )[i]
            
            PCR = primary_convection_resistance(
                Section, "Dittus-Boelter", PrimaryT_film, T_PrimaryWall, SecondarySidePressure, i
                )
            SCR = secondary_convection_resistance(
                Section, SecondaryT_film, T_SecondaryWall, x_in, SecondarySidePressure, i
                )
             
            # [cm^2 K/W]
            inverseU_total = (PCR + conduction_resistance(Section, T_PrimaryWall, SecondarySidePressure, i) + SCR) \
            * outer_area(Section)[i] + R_F_primary.magnitude + R_F_secondary.magnitude

            U_total.magnitude = 1 / inverseU_total  # [W/ cm^2 K]

            return T_PrimaryWall, T_SecondaryWall, U_total.magnitude


def temperature_profile(
        Section, InnerAccumulation, OuterAccumulation, m_h_leakagecorrection, SecondarySidePressure, calendar_year
        ):
    if SecondarySidePressure == 4.593: 
        T_sat_secondary = 260.1 + 273.15
    elif SecondarySidePressure < 4.593:
        T_sat_secondary = 258.1 + 273.15  
    
    PrimaryWall = []
    PrimaryBulk = []
    SecondaryBulk = []
    SecondaryWall = []

    for i in range(Section.NodeNumber - 1):
        if i == 0:
            # Temperatures entering SG --> not first node temps.
            T_PrimaryBulkIn = T_PrimaryIn  # [K]
            T_SecondaryBulkIn = T_sat_secondary
            x_in = 0

        Cp_h = ld.HeatCapacity("PHT", T_PrimaryBulkIn, SecondarySidePressure)
        Cp_c = ld.HeatCapacity("SHT", T_SecondaryBulkIn, SecondarySidePressure)
        T_wh, T_wc, U = wall_temperature(
            Section, i, T_PrimaryBulkIn, T_SecondaryBulkIn, x_in, InnerAccumulation, OuterAccumulation, calendar_year,
            SecondarySidePressure)

        if Section.Length.label[i] != "preheater start" and Section.Length.label[i] != "preheater" and \
                Section.Length.label[i] != "thermal plate":

            T_PrimaryBulkOut = T_PrimaryBulkIn - (U * (T_PrimaryBulkIn - T_SecondaryBulkIn) * outer_area(Section)[i] \
                                                  * nc.TotalSGTubeNumber) / (Cp_h * m_h_leakagecorrection)
            # Ti = Ti+1 = Tsat
            T_SecondaryBulkOut = T_sat_secondary

            # for one tube --> multiply by NumberTubes = for all tubes
            # Q = m_h_leakagecorrection*Cp_h*(T_PrimaryBulkIn-T_PrimaryBulkOut)

            # U based on one tube (heat transfer area x number tubes) --> would cancel out in U calc (1/h*NA) NA
            Q = U * (T_PrimaryBulkOut - T_SecondaryBulkOut) * outer_area(Section)[i] * nc.TotalSGTubeNumber

            x_out = x_in + Q / (MassFlow_c.magnitude * EnthalpySaturatedSteam.magnitude)

            x_in = x_out
            T_PrimaryBulkIn = T_PrimaryBulkOut
            T_SecondaryBulkIn = T_SecondaryBulkOut

            PrimaryBulk.append(T_PrimaryBulkOut)
            SecondaryBulk.append(T_SecondaryBulkOut)
            PrimaryWall.append(T_wh)
            SecondaryWall.append(T_wc)
        else:
            if Section.Length.label[i] == "preheater start":
                C_min = Cp_c * MassFlow_c.magnitude  # [J/g K]*[g/s] = [J/Ks] ] [W/K]
                C_max = Cp_h * m_h_leakagecorrection
                C_r = C_min / C_max
                Q_max = C_min * (T_PrimaryBulkIn - T_PreheaterIn)
                TotalArea = sum(outer_area(Section)[i:(Section.NodeNumber - 1)])
                NTU = U * TotalArea * nc.TotalSGTubeNumber / C_min

                eta = (1 - np.exp(-NTU * (1 - C_r))) / (1 - C_r * np.exp(-NTU * (1 - C_r)))
                # eta = 1-np.exp(((NTU**0.28)/C_r)*(np.exp(-C_r*(NTU**0.78))-1))
                Q_NTU = eta * Q_max

                T_PrimaryBulkOutEnd = T_PrimaryBulkIn - Q_NTU / (m_h_leakagecorrection * Cp_h)
                T_SecondaryBulkOutEnd = T_PreheaterIn + Q_NTU / (MassFlow_c.magnitude * Cp_c)
                T_SecondaryBulkIn = T_SecondaryBulkOutEnd

            # Guessing cold-side temperatures for remaining nodes inside preheater
            T_SecondaryBulkOut = T_SecondaryBulkOut - i
            T_PrimaryBulkOut = T_PrimaryBulkIn \
                - (U * outer_area(Section)[i] * nc.TotalSGTubeNumber / (Cp_h * m_h_leakagecorrection)) \
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

    PrimaryBulk.append(PrimaryBulk[20])
    PrimaryWall.append(PrimaryWall[20])

#     print ([j - 273.15 for j in PrimaryBulk])
#     print ([j-273.15 for j in PrimaryWall])
#     print ()
#     print ([j-273.15 for j in SecondaryBulk])
#     print ([j-273.15 for j in SecondaryWall])
#     print()
    return PrimaryBulk


def station_events(calendar_year):
    if calendar_year <= 1992:
        # divider plate leakage rates estimated based on AECL work
        InitialLeakage = 0.035 # fraction of total SG inlet mass flow
        YearlyRateLeakage = 0.0065 # yearly increase to fraction of total SG inlet mass flow
        SHTPressure = 4.593  # MPa
    else:
        # PLNGS pressure reduction in 1992 + divider plate replacement
        InitialLeakage = 0.02 
        YearlyRateLeakage = 0
        SHTPressure = 4.593 - (125 / 1000) # MPa
    
    Leakage = InitialLeakage + (calendar_year - 1983) * YearlyRateLeakage
    DividerPlateMassFlow = MassFlow_h.magnitude * Leakage
    m_h_leakagecorrection = MassFlow_h.magnitude - DividerPlateMassFlow
    
    return SHTPressure, m_h_leakagecorrection, DividerPlateMassFlow


def energy_balance(SteamGeneratorOutputNode, j):
    year = (j / 8760) 
    calendar_year = year + 1983
    
    [SecondarySidePressure, RemainingPHTMassFlow, MasssFlow_dividerplate.magnitude] = station_events(calendar_year)

    Balance = []
    for Zone in ld.SGZones:
        Zone.PrimaryBulkTemperature = temperature_profile(
            Zone, ld.SGZones[12].InnerOxThickness, ld.SGZones[12].OuterOxThickness, RemainingPHTMassFlow,
            SecondarySidePressure, calendar_year)
        
        x = (Zone.TubeNumber / nc.TotalSGTubeNumber) * RemainingPHTMassFlow \
            * ld.Enthalpy("PHT", Zone.PrimaryBulkTemperature[SteamGeneratorOutputNode], SecondarySidePressure)

        Balance.append(x)
        Enthalpy = (
            sum(Balance) + MasssFlow_dividerplate.magnitude * ld.Enthalpy("PHT", T_PrimaryIn, SecondarySidePressure)
            ) / MassFlow_h.magnitude

    RIHT = ld.TemperaturefromEnthalpy("PHT", Enthalpy, SecondarySidePressure)
    return RIHT
# print (energy_balance(21, 1)-273.15)

