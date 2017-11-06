import lepreau_data as ld
import numpy as np
import constants as nc

T_sat = 260.1 + 273.15
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

MassFlow_c.magnitude = 275 * 1000
MassFlow_h.magnitude = 1900 * 1000

ShellDiameter.magnitude = 2.28 * 100

for i in [MassFlow_c, MassFlow_h, MasssFlow_dividerplate]:
    i.unit = "g/s"

for i in [ShellDiameter, EquivalentDiameter, TubePitch]:
    i.unit = "cm"

MassFlux_c.magnitude = MassFlow_c.magnitude / ((np.pi / 4) * (ShellDiameter.magnitude ** 2))
MassFlux_c.unit = "g/cm^2 s"
EnthalpySaturatedSteam.magnitude = 2800.4
EnthalpySaturatedSteam.unit = "J/g"

TubePitch.magnitude = 2.413


def thermal_conductivity(Twall, material):
    if material == "Alloy-800" or material == "Alloy800" or material == "A800" or material == "Alloy 800":
        return (11.450 + 0.0161 * Twall) / 100  # M.K.E thesis for Alloy-800 tubing [W/cm K]
    elif material == "water":
        return ld.thermal_conductivityH2O("PHT", Twall)
    elif material == "magnetite":
        return 1.4 / 100  # [W/cm K]
    else:
        print ("Error: material not specified")
        return None


def fouling_resistance(InnerAccumulation, OuterAccumulation):
    # thickness/thermal conductivity [cm]/[W/cm K] = [cm^2 K/W]
    # [g/cm^2]/[g/cm^3] = [cm]

    InnerThickness = [i / nc.Fe3O4Density for i in InnerAccumulation]
    OuterThickness = [i / nc.Fe3O4Density for i in OuterAccumulation]

    # [cm]/ [W/cm K] =[cm^2 K/W]
    InnerFouling = [i / 2 * thermal_conductivity(None, "magnetite") for i in InnerThickness]
    OuterFouling = [i / thermal_conductivity(None, "magnetite") for i in OuterThickness]
    return [x + y for x, y in zip(InnerFouling, OuterFouling)]


def conduction_resistance(Section, Twall, i):
    # R_conduction = L/kA (L = thickness)
    # A = area with shape factor correcrion factor for cylindrical pipe
    # R_conduction = ln(D_o/D_i)/(2pikl) [K/W]
    # l = length, k =thermal conductivity coefficient [W/cm K]

    k_w = thermal_conductivity(Twall, "Alloy-800")  # [W/cm K]

    Rcyl_numerator = np.log(Section.OuterDiameter[i] / Section.Diameter[i])
    Rcyl_denominator = 2 * np.pi * k_w * Section.Length.magnitude[i]

    R_cond = Rcyl_numerator / Rcyl_denominator  # [K/W]
    return R_cond


def primary_convection_resistance(Section, correlation, Tfilm, i):
    # R_1,Convective = 1/(h_1,convectove *A)
    # A = inner area (based on inner diameter)
    # [W/cm K]/ [cm] = [W/K cm^2]
    h_i = nusseltnumber(Section, correlation, Tfilm, i) * thermal_conductivity(Tfilm, "water") / Section.Diameter[i]
    return 1 / (h_i * inner_area(Section)[i])  # [K/W]


def secondary_convection_resistance(Section, Tfilm, Twall, x_in, i):
    # split into boiling and non-boiling (preheater) sections
    EquivalentDiameter.magnitude = 4 * ((TubePitch.magnitude ** 2) - (np.pi * (Section.OuterDiameter[i] ** 2)) / 4) \
                                        / (np.pi * Section.OuterDiameter[i])

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
        h_o = (thermal_conductivity(Tfilm, "water") * 0.25 / Section.OuterDiameter[i]) \
            * ((Section.OuterDiameter[i] * G_e / ld.Viscosity("water", "SHT", Tfilm)) ** 0.6) \
            * (ld.HeatCapacity("SHT", Tfilm)
               * ld.Viscosity("water", "SHT", Tfilm) / thermal_conductivity(Tfilm, "water")) ** 0.33

    else:
        if i <= 10:
            x = x_in * 100
        else:
            x = 24 - 1.15 * i

        rho_v = 1000 * 23.753 / (100 ** 3)  # [g/cm^3]
        p_crit = 22.0640  # [MPa]
        # V*rho =[cm/s]([cm^2]*[g/cm^3] = [g/s]

        F = (1 + (x * prandtl(T_sat, i) * ((ld.Density("water", "SHT", T_sat) / rho_v) - 1))) ** 0.35

        Re_D = EquivalentDiameter.magnitude * MassFlux_c.magnitude / ld.Viscosity("water", "SHT", T_sat)
        # 4*MassFlow/((ld.Viscosity("water", "SHT", T_sat)/1000)*np.pi*Shell_ID)

        h_l = 0.023 * ld.thermal_conductivityH2O("SHT", T_sat) * ((Re_D) ** 0.8) * (prandtl(T_sat, i) ** 0.4) \
            / EquivalentDiameter.magnitude  # [W/cm^2 K]

        Q_prime = F * h_l * (Twall - T_sat)  # [W/cm^2]
        S = (1 + 0.055 * (F ** 0.1) * (Re_D) ** 0.16) ** (-1)
        A_p = (55 * ((4.70 / p_crit) ** 0.12) * ((-np.log10(4.70 / p_crit)) ** (-0.55)) * (nc.H2MolarMass) ** (-0.5)) \
            / (100 ** 2)
        C = ((A_p * S) / (F * h_l) ** 2) * (Q_prime) ** (4 / 3)
        coeff = [1, -C, -1]
        cubic_solution = np.roots(coeff)
        q = cubic_solution[0]

        h_o = F * (q ** (3 / 2)) * h_l  # [W/cm^2 K

    return 1 / (h_o * outer_area(Section)[i])  # K/W


def nusseltnumber(Section, correlation, Temperature, i):
    # ViscosityH2O = ld.Viscosity("water", "PHT", Temperature)
    # DensityH2O = ld.Density("water", "PHT", Temperature)

    # Re_D = Section.Velocity[i]*Section.Diameter[i]/(ViscosityH2O*DensityH2O)
    MassFlux_h = MassFlow_h.magnitude / (nc.TotalSGTubeNumber * ((np.pi / 4) * (Section.Diameter[i] ** 2)))
    Re_D = [(MassFlux_h / ld.Viscosity("water", "PHT", Temperature)) * i for i in Section.Diameter]
    # print (Re_D[i], Re_[i])

    if correlation == "Dittus-Boelter":
        n = 1 / 3
    elif correlation == "Colburn":
        n = 0.4
    else:
        print ("Error: Correlation not recognized in nusseltnumber function")
    C = 0.023
    m = 4 / 5

    Nu = C * (Re_D[i] ** m) * ((prandtl(Temperature, i)) ** n)
    return Nu


def prandtl(Temperature, i):
    # Need a better reference for Cp/viscosity data
    ViscosityH2O = ld.Viscosity("water", "PHT", Temperature)
    Pr = ld.HeatCapacity("PHT", Temperature) * ViscosityH2O / ld.thermal_conductivityH2O("PHT", Temperature)
    return Pr


def inner_area(Section):
    return [np.pi * x * y for x, y in zip(Section.Diameter, Section.Length.magnitude)]


def outer_area(Section):
    return [np.pi * x * y for x, y in zip(Section.OuterDiameter, Section.Length.magnitude)]


def wall_temperature(Section, i, T_PrimaryBulkIn, T_SecondaryBulkIn, x_in, InnerAccumulation, OuterAccumulation, j):
    # i = each node of SG tube
    T_PrimaryWall = T_PrimaryBulkIn - (1 / 3) * (T_PrimaryBulkIn - T_SecondaryBulkIn)
    T_SecondaryWall = T_PrimaryBulkIn - (1 / 3) * (T_PrimaryBulkIn - T_SecondaryBulkIn)

    for k in range(50):
        WT_h = T_PrimaryWall
        WT_c = T_SecondaryWall

        PrimaryT_film = (T_PrimaryBulkIn + T_PrimaryWall) / 2
        SecondaryT_film = (T_SecondaryBulkIn + T_SecondaryWall) / 2

        h_i.magnitude = 1 / (primary_convection_resistance(Section, "Dittus-Boelter", PrimaryT_film, i)
                             * inner_area(Section)[i])
        h_o.magnitude = 1 / (secondary_convection_resistance(Section, SecondaryT_film, T_SecondaryWall, x_in, i)
                             * outer_area(Section)[i])
        # R = 1/hA --> h=A*1/R

        U_h1 = conduction_resistance(Section, T_PrimaryWall, i) * inner_area(Section)[i]
        U_h2 = secondary_convection_resistance(
            Section, SecondaryT_film, T_SecondaryWall, x_in, i
            ) * inner_area(Section)[i]
        U_h.magnitude = 1 / (U_h1 + U_h2)  # [W/ cm^2 K]

        U_c1 = conduction_resistance(Section, T_PrimaryWall, i) * outer_area(Section)[i]
        U_c2 = primary_convection_resistance(Section, "Dittus-Boelter", PrimaryT_film, i) * outer_area(Section)[i]
        U_c.magnitude = 1 / (U_c1 + U_c2)  # [W/cm^2 K]

        T_PrimaryWall = (T_PrimaryBulkIn * h_i.magnitude + T_SecondaryBulkIn * U_h.magnitude)\
            / (h_i.magnitude + U_h.magnitude)
        T_SecondaryWall = (T_SecondaryBulkIn * h_o.magnitude + T_PrimaryBulkIn * U_c.magnitude)\
            / (h_o.magnitude + U_c.magnitude)

        RE1 = (T_PrimaryWall - WT_h)
        RE2 = (T_SecondaryWall - WT_c)

        if abs(RE1) <= 0.01 and abs(RE2) <= 0.01:
            R_F_primary.magnitude = fouling_resistance(InnerAccumulation, OuterAccumulation)[i]  # [cm^2 K/W]
            
            year = (j/8760) 
            calendaryear = year + 1983
            Initial_SHTSFouling = 0.0005 # [g/cm2] same as for PHT (~ 1 um/a)
            
            if 1983 <= calendaryear < 1987:
                SHTSFouling = Initial_SHTSFouling + year*0.0016
                # secondary side fouling slope based on 1/2 that for average primary side cold leg deposit
            else:
                # CPP installation (late 1986) reduces secondary side crud by 50% 
                # (assumed proportional red. in deposit formation)
                SHTSFouling = Initial_SHTSFouling + year*0.0008
                
            R_F_secondary.magnitude = fouling_resistance([0]*Section.NodeNumber, [SHTSFouling]*Section.NodeNumber)[i]
            
            # [cm^2 K/W]
            inverseU_total = (primary_convection_resistance(Section, "Dittus-Boelter", PrimaryT_film, i)
                              + conduction_resistance(Section, T_PrimaryWall, i)
                              + secondary_convection_resistance(Section, SecondaryT_film, T_SecondaryWall, x_in, i)) \
                * outer_area(Section)[i] + R_F_primary.magnitude + R_F_secondary.magnitude

            U_total.magnitude = 1 / inverseU_total  # [W/ cm^2 K]

            return T_PrimaryWall, T_SecondaryWall, U_total.magnitude


def temperature_profile(Section, InnerAccumulation, OuterAccumulation, m_h_leakagecorrection, j):
    PrimaryWall = []
    PrimaryBulk = []
    SecondaryBulk = []
    SecondaryWall = []

    for i in range(Section.NodeNumber - 1):
        if i == 0:
            # Temperatures entering SG --> not first node temps.
            T_PrimaryBulkIn = 583.15  # [K]
            T_SecondaryBulkIn = T_sat
            x_in = 0

        Cp_h = ld.HeatCapacity("PHT", T_PrimaryBulkIn)
        Cp_c = ld.HeatCapacity("SHT", T_SecondaryBulkIn)
        T_wh, T_wc, U = wall_temperature(
            Section, i, T_PrimaryBulkIn, T_SecondaryBulkIn, x_in, InnerAccumulation, OuterAccumulation, j
            )

        if Section.Length.label[i] != "preheater start" and Section.Length.label[i] != "preheater" and \
                Section.Length.label[i] != "thermal plate":

            T_PrimaryBulkOut = T_PrimaryBulkIn\
                             - (U * (T_PrimaryBulkIn - T_SecondaryBulkIn) * outer_area(Section)[i] * nc.TotalSGTubeNumber)\
                             / (Cp_h * m_h_leakagecorrection)
            # Ti = Ti+1 = Tsat
            T_SecondaryBulkOut = T_sat

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
                # print (T_PrimaryBulkOutEnd-273.15, T_SecondaryBulkOutEnd-273.15)

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


def energy_balance(SteamGeneratorOutputNode, j):
    InitialLeakage = 0.02
    YearlyRateLeakage = 0.006
    Leakage = InitialLeakage + (j / 8760) * YearlyRateLeakage
    MasssFlow_dividerplate.magnitude = MassFlow_h.magnitude * Leakage
    m_h_leakagecorrection = MassFlow_h.magnitude - MasssFlow_dividerplate.magnitude

    Balance = []
    for Zone in ld.SGZones:
        Zone.PrimaryBulkTemperature = temperature_profile(
            Zone, ld.SGZones[0].InnerOxThickness, ld.SGZones[0].OuterOxThickness, m_h_leakagecorrection, j
            )
        
        x = (Zone.TubeNumber / nc.TotalSGTubeNumber) * m_h_leakagecorrection \
            * ld.Enthalpy("PHT", Zone.PrimaryBulkTemperature[SteamGeneratorOutputNode])

        Balance.append(x)
        Enthalpy = (sum(Balance) + MasssFlow_dividerplate.magnitude * ld.Enthalpy("PHT", T_PrimaryIn)) \
            / MassFlow_h.magnitude

    RIHT = ld.TemperaturefromEnthalpy("PHT", Enthalpy)
    return RIHT
# print (energy_balance(21, 1)-273.15)