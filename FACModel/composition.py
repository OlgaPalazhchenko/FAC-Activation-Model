import thermochemistry_and_constants as nc
import lepreau_data as ld
import numpy as np
import sg_heattransfer as SGHX

CONCENTRATION_LITHIUM = 0.000225095734  # [mol/L] #0.00063095734

# Temperature thermochemistry_and_constants for equilibrium/hydrolysis thermochemistry_and_constants
DEBYE_HUCKEL_POLYNOMIAL = [3.29E-10, -1.709E-07, 0.00003315, -0.0009028, 0.5027]
KwPOLYNOMIAL = [-7.226E-10, 1.32661E-06, -0.000959311, 0.32765297, -55.86334915]
KLiPOLYNOMIAL = [0.00000675, -0.0048, -0.7532]

# KFe2SchikorrPOLYNOMIAL = [6.66667E-11, -2.128E-07, 2.52E-04, -0.142254891, 35.94096367]
KFeOHPOLYNOMIAL = [-4E-10, 8.1013E-07, -0.00062963, 0.230763, -41.5545]
KFeOH2POLYNOMIAL = [-4E-10, 8.63467E-07, -0.00073731, 0.311499, -67.8248]
KFeOH3POLYNOMIAL = [-4.667E-10, 1.0496E-06, -0.000935775, 0.413186, -97.4709]
# Fe3+
# KFeOH3SchikorrPOLYNOMIAL = [-1.35525E-20, 5.33333E-08, -0.00010368, 0.074951307, -27.34968224]
# KFeOH4SchikorrPOLYNOMIAL = [2.46667E-09, -4.72693E-06, 0.003334163, -1.007239881, 89.37956761]
 
# KNiOHPOLYNOMIAL = [3.10167E-09, -5.8294E-06, 0.004039757, -1.201661204, 119.296086]
# KNiOH2POLYNOMIAL = [-3.04309E-10, 6.66757E-07, -0.000580137, 0.255160282, -60.03240922]
# KNiOH3POLYNOMIAL = [8.49674E-10, -1.5126E-06, 0.000945292, -0.211217025, -16.71117746]

# def purification_system():
#     
#     
#     return None
# 
# 
# def arrhenius_activaton_energy():
#     Section = ld.SteamGenerator[57]
# # #     year = (j / 8760) 
# # #     calendar_year = year + 1983
# #     
# #     Tube = SGHX.closest_tubelength(1980)
# 
#         # based on chemical clean (AECL) of plngs pulled tubes in 1994
#         # assumed this average growth rate applies to all tubes (data only available for 2 x 19.8 m long tubes)
#         
#         # [g/m^2 /yr] --> [g/cm^2 /yr]
#     Growth_coldavg = 34.5 / (100 ** 2)  # preheater
#     Growth_hotavg = 16.15 / (100 ** 2)
#     
#     Saturation = ld.UnitConverter(
#         Section, "Mol per Kg", "Grams per Cm Cubed", Section.SolutionOxide.FeSatFe3O4, None,
#         None, None, nc.FeMolarMass, None
#         )
#     
#     Concentration = ld.UnitConverter(
#         Section, "Mol per Kg", "Grams per Cm Cubed", Section.SolutionOxide.FeTotal, None,
#         None, None, nc.FeMolarMass, None
#         )
#     
#     ConcentrationGradient_cold = Concentration[19] - Saturation[19]
#     ConcentrationGradient_hot = Concentration[7] - Saturation[7]
#     
#     if ConcentrationGradient_cold == 0 or ConcentrationGradient_hot == 0:
#         kp_cold = nc.KpFe3O4
#         kp_hot = nc.KpFe3O4
#     
#     else:
#         kp_cold = Growth_coldavg / ConcentrationGradient_cold
#         kp_hot = Growth_hotavg / ConcentrationGradient_hot
#         
#         kp_cold = kp_cold / (365 * 24 * 3600)  # [cm/yr] * 365 d/yr * 24 h/d * 3600 s/h = [cm/s]
#         kp_hot = kp_hot / (365 * 24 * 3600)
#     # assume first ~4 m of hot leg have minimal deposition 
#     
#     if kp_cold and kp_hot < 0:
#         kp_cold = nc.KpFe3O4
#         kp_hot = nc.KpFe3O4
#         
#     T_cold = []
#     T_hot = []
#         
# #     for i in range(SteamGenerator[Tube].NodeNumber - 1):
# #         if (SteamGenerator[Tube].Length.label[i] == "preheater"
# #             or SteamGenerator[Tube].Length.label[i] == "preheater start"):
# #             # preheater area
# #             T_cold.append(SteamGenerator[Tube].PrimaryBulkTemperature[i])
# #             
# #         if 3 < i < 12:
# #             T_hot.append(SteamGenerator[Tube].PrimaryBulkTemperature[i])
# 
#     T_coldavg = Section.PrimaryBulkTemperature[19]  # preheater
#     T_hotavg = Section.PrimaryBulkTemperature[7]
#     
# #     T_coldavg = sum(T_cold) / len(T_cold)  # preheater
# #     T_hotavg = sum(T_hot) / len(T_hot)  # hot leg not including preheater or first 4 m
# #     print (kp_cold, kp_hot)
#     # natural logarithm, ln, [J/K mol][K] = [J/mol]
#     ActivationEnergy = (nc.R * T_coldavg * T_hotavg / (T_coldavg - T_hotavg)) * np.log(kp_cold / kp_hot)
# 
#     A = kp_cold / np.exp(-ActivationEnergy / (nc.R * T_coldavg)) 
#     
#     return ActivationEnergy, A
 

def plngs_precipitation_kinetics(Section, j):
    
    ActivationEnergy = 4624 * nc.R  
    A = np.exp(5.369)    
    kp = [A * np.exp(-ActivationEnergy / (nc.R * i)) for i in Section.PrimaryBulkTemperature]
  
    return kp


def temperature_dependent_parameters(Section):
    Celsius = [i - 273.15 for i in Section.PrimaryBulkTemperature]
    # Equilibrium and Debye-Huckel thermochemistry and_ onstants - polynomials as a function of temperature
    # Coeff1*x^4 + Coeff2*x^3 + Coeff3*x^2 + Coeff4*x + Coeff5, depends on # elements in coeff list
    Section.DebyeHuckelConstant = (np.polyval(DEBYE_HUCKEL_POLYNOMIAL, Celsius))
    Section.k_W = 10 ** (np.polyval(KwPOLYNOMIAL, Section.PrimaryBulkTemperature))
    Section.k_Li = 10 ** (np.polyval(KLiPOLYNOMIAL, Section.PrimaryBulkTemperature))
        
    Section.k_FeOH_hydrolysis = 10 ** (np.polyval(KFeOHPOLYNOMIAL, Section.PrimaryBulkTemperature))
    Section.k_FeOH2_hydrolysis = 10 ** (np.polyval(KFeOH2POLYNOMIAL, Section.PrimaryBulkTemperature))
    Section.k_FeOH3_hydrolysis = 10 ** (np.polyval(KFeOH3POLYNOMIAL, Section.PrimaryBulkTemperature))
#     Section.k_NiOH = 10 ** (np.polyval(KNiOHPOLYNOMIAL, Section.PrimaryBulkTemperature))
#     Section.k_NiOH2 = 10 ** (np.polyval(KNiOH2POLYNOMIAL, Section.PrimaryBulkTemperature))
#     Section.k_NiOH3 = 10 ** (np.polyval(KNiOH3POLYNOMIAL, Section.PrimaryBulkTemperature))

    return None


def bulkpH_calculator(Section):  # Bulk pH calculation
    temperature_dependent_parameters(Section)
    H = []
    ConcentrationH = 0.0000000001  # initial guess [mol/kg]
    gamma_1 = 1  # initial guess
    for i in range(Section.NodeNumber):
        ConcentrationOH = Section.k_W[i] / ((gamma_1 ** 2) * (ConcentrationH))

        # At high temp, LiOH doesn't dissociate 100% - eq'm established: LiOH(aq) <-> Li+(aq) + OH-(aq)
        # KLi = ([Li+]gamma_1 *[OH-]gamma_1)/[LiOH]; 
        # CONCENTRATION_LITHIUM = ConcentrationLi  + LiConcentrationOH (sub in eq'm expression for LiConcentrationOH
        ConcentrationLi = Section.k_Li[i] * CONCENTRATION_LITHIUM / (
            Section.k_Li[i] + ((gamma_1 ** 2) * ConcentrationOH)
            )
        # k_W = ConcentrationH*gamma_1*ConcentrationOH*gamma_1;   ConcentrationOH = k_W/ConcentrationH*gamma_1**2
        # H+ + Li+ = OH-;                    ConcentrationH = ConcentrationLi + k_W/ConcentrationH*gamma_1**2
        for k in range(50):
            # no more than 10 iterations should really be necessary for convergence at the provided error    
            # function of H
            FH = ConcentrationH ** 2 + (ConcentrationH * ConcentrationLi) - (Section.k_W[i] / (gamma_1 ** 2))  
            DFH = 2 * ConcentrationH + ConcentrationLi  # derivative of function

            NewConcentrationH = ConcentrationH - (FH / DFH)

            RE = abs((NewConcentrationH - ConcentrationH) / NewConcentrationH)

            ConcentrationH = NewConcentrationH

            ConcentrationOH = Section.k_W[i] / ((gamma_1 ** 2) / ConcentrationH)
            ConcentrationLi = (Section.k_Li[i] * CONCENTRATION_LITHIUM) / (
                Section.k_Li[i] + ((gamma_1 ** 2) * ConcentrationOH)
                )
            IonicStrength = ((1 ** 2) * ConcentrationH + (1 ** 2) * ConcentrationLi + (1 ** 2) * ConcentrationOH) / 2
            # Davies equation loggamma_1 = -DebyeHuckConst*(z**2)*[(sqrt(I)/(1+sqrt(I)))-beta*I]
            gamma_1 = 10 ** (
                Section.DebyeHuckelConstant[i] * (1 ** 2) * (
                    ((IonicStrength ** 0.5) / (1 + (IonicStrength ** 0.5))) - 0.2 * IonicStrength
                    )
                ) 
            # All entries in ConcentrationH list must meet error convergence requirement
            if RE < 0.0000001:
                break 
        H.append(ConcentrationH)
    return H


def hydrolysis(Section, FeTotal, ConcentrationH):
    ActivityCoefficient1 = [1.0] * Section.NodeNumber  # initial estimate for activities (+/- 1 charged ions)
    ActivityCoefficient2 = [1.0] * Section.NodeNumber  # (+/- 2 charged ions)
    # Concentration NiTotal <<FeTotal, so Ni species left out of ionic strength calcs
    for i in range(50):
        gamma_1itr = ActivityCoefficient1
        gamma_2itr = ActivityCoefficient2

        # ConcentrationFe/NiTotal and ConcentrationH are not re-evaluated here, thus  not used for them
        ConcentrationOH = [x / (y * (z ** 2)) for x, y, z in zip(Section.k_W, ConcentrationH, ActivityCoefficient1)]
        ConcentrationLi = [CONCENTRATION_LITHIUM * x / (x + (y * (z ** 2))) for x, y, z in zip(
            Section.k_Li, ConcentrationOH, ActivityCoefficient1
            )]
        
        ConcentrationFe2 = [x / (1
                                 + (y * z / (r * (s ** 2)))
                                 + (t * z / ((r ** 2) * (s ** 2)))
                                 + (q * z / ((r ** 3) * (s ** 4)))) 
            for x, y, z, r, s, t, q in zip(
                FeTotal, Section.k_FeOH_hydrolysis, ActivityCoefficient2, ConcentrationH, ActivityCoefficient1,
                Section.k_FeOH2_hydrolysis, Section.k_FeOH3_hydrolysis
                )]
        ConcentrationFeOH = [x * y * z / (q * (r ** 2)) for x, y, z, q, r in zip(
            Section.k_FeOH_hydrolysis, ActivityCoefficient2, ConcentrationFe2, ConcentrationH, ActivityCoefficient1
            )]
        ConcentrationFeOH2 = [x * y * z / ((q ** 2) * (r ** 2)) for x, y, z, q, r in zip(
            Section.k_FeOH2_hydrolysis, ActivityCoefficient2, ConcentrationFe2, ConcentrationH, ActivityCoefficient1
            )] 
        ConcentrationFeOH3 = [x * y * z / ((q ** 3) * (r ** 4)) for x, y, z, q, r, in zip(
            Section.k_FeOH3_hydrolysis, ActivityCoefficient2, ConcentrationFe2, ConcentrationH, ActivityCoefficient1
            )]

#         ConcentrationNi2 = [
#             x / (1 + (y * z / (r * (s ** 2))) + (t * z / ((r ** 2) * (s ** 2))) + (q * z / ((r ** 3) * (s ** 4)))) 
#             for x, y, z, r, s, t, q in zip(NiTotal, Section.k_NiOH, ActivityCoefficient2, ConcentrationH, 
#                                            ActivityCoefficient1, Section.k_NiOH2, Section.k_NiOH3)
#             ]
#         ConcentrationNiOH = [x * y * z / (q * (r ** 2)) for x, y, z, q, r in zip(
#             Section.k_NiOH, ActivityCoefficient2, ConcentrationNi2, ConcentrationH, ActivityCoefficient1
#             )]
#         ConcentrationNiOH2 = [x * y * z / ((q ** 2) * (r ** 2)) for x, y, z, q, r in zip(
#             Section.k_NiOH2, ActivityCoefficient2, ConcentrationNi2, ConcentrationH, ActivityCoefficient1
#             )] 
#         ConcentrationNiOH3 = [x * y * z / ((q ** 3) * (r ** 4)) for x, y, z, q, r, in zip(
#             Section.k_NiOH3, ActivityCoefficient2, ConcentrationNi2, ConcentrationH, ActivityCoefficient1
#             )]

        # FeOH2 and NiOH2 left out of ionic strength calc due to electroneutrality principle
        IonicStrength = [((1 ** 2) * x
                          + ((-1) ** 2) * y
                          + (2 ** 2) * z
                          + (1 ** 2) * q
                          + ((-1) ** 2) * r
                          + (1 ** 2) * v) / 2 
                         for x, y, z, q, r, v in zip(
                             ConcentrationH, ConcentrationOH, ConcentrationFe2, ConcentrationFeOH, ConcentrationFeOH3,
                             ConcentrationLi
                             )
                         ]
        ActivityCoefficient1 = [10 ** (
            x * (1 ** 2) * (((y ** 0.5) / (1 + (y ** 0.5))) - 0.2 * y))
                                for x, y in zip(-Section.DebyeHuckelConstant, IonicStrength
            )] 
        ActivityCoefficient2 = [10 ** (
            -x * ((-2) ** 2) * (((y ** 0.5) / (1 + (y ** 0.5))) - 0.2 * y))
                                for x, y in zip(Section.DebyeHuckelConstant, IonicStrength
            )]

        RE = [((x - y) / x) for x, y in zip(gamma_1itr, ActivityCoefficient1)]

        RE2 = [((x - y) / x) for x, y in zip(gamma_2itr, ActivityCoefficient2)]
        if RE < [0.00000001] * Section.NodeNumber and RE2 < [0.00000001] * Section.NodeNumber:

            return ConcentrationFe2, ConcentrationFeOH2, ActivityCoefficient1, ActivityCoefficient2


def iron_solubility_TL(Section, Condition):
    # As temperature changes, thermochemistry and constants (k_Li, k_w, etc. change), pH changes as well, 
    # both affecting solubility --> Bulk FeSatFe3O4
    Section.SolutionOxide.ConcentrationH = bulkpH_calculator(Section)
        
    Gibbs_Energies = [-62700, -4300, 56000, 127840, 74000, 128800]  # [J/mol]
    Entropies = [-51.63, -108.96, -102.35, -188.71, -66.60, -187.76]  # [J/ K mol]
    Delta_T = [i - 298 for i in Section.PrimaryBulkTemperature]  # [K]
    
    Fe3O4_SolubilityConstants_TL = []
    
    for x, y in zip(Gibbs_Energies, Entropies):
        q = [np.exp(-(x - y * k) / (nc.R * i)) for i, k in zip(Section.PrimaryBulkTemperature, Delta_T)]
        Fe3O4_SolubilityConstants_TL.append(q)

    k_Fe2 = Fe3O4_SolubilityConstants_TL[0]
    k_FeOH = Fe3O4_SolubilityConstants_TL[1]
    k_FeOH2 = Fe3O4_SolubilityConstants_TL[2]
    k_FeOH3 = Fe3O4_SolubilityConstants_TL[3]
    k_FeOH3_Fe3 = Fe3O4_SolubilityConstants_TL[4]
    k_FeOH4_Fe3 = Fe3O4_SolubilityConstants_TL[5]

    # based on high temperature henry's law constant data for H2 (some taken from T&L ref. some from sol. program DHL
    k_H2 = [-4.1991 * i + 2633.2 for i in Section.PrimaryBulkTemperature]

    # convert H2 conc. from cm^3/kg to mol/kg then concentration to P ("fugacity") using Henry's constant
    P_H2 = [(nc.H2 * nc.H2Density / nc.H2MolarMass) * i for i in k_H2]
    
    # need to put these in line with hydrolysis activity coeffs    
    if Condition == "initial":
        gamma_1 = [0.95] * Section.NodeNumber
        gamma_2 = [0.95] * Section.NodeNumber
    else:
        gamma_1 = hydrolysis(Section, Section.SolutionOxide.FeTotal, Section.SolutionOxide.ConcentrationH)[2]
        gamma_2 = hydrolysis(Section, Section.SolutionOxide.FeTotal, Section.SolutionOxide.ConcentrationH)[3]
           
    b = [0, 1, 2, 3, 3, 4]
    
    Activity = []
    Constants = [k_Fe2, k_FeOH, k_FeOH2, k_FeOH3, k_FeOH3_Fe3, k_FeOH4_Fe3]
    
    for hydroxyls, k in zip (b, Constants):
        if hydroxyls == 0:
            gamma_oxidation = gamma_2
        elif hydroxyls == 1 or hydroxyls == 3 or hydroxyls == 4:  # -1 for Fe2+, +1 for Fe2+, -1 for Fe3+ species
            gamma_oxidation = gamma_1
        elif hydroxyls == 2 or hydroxyls == 4:
            # 0 charge on Fe(OH)2^0 and Fe(OH)3^0 = electroneutrality, no gamma
            gamma_oxidation = [1] * Section.NodeNumber
        else: 
            None
        
        if k == Constants[4] or k == Constants[5]:
            charge = 3  # ferric, Fe3+ species
        else:
            charge = 2  # ferrous, Fe2+ species

        q = [x * ((y * g_1) ** (charge - hydroxyls)) * (z ** ((4 / 3) - (charge / 2))) / (g_o) for
             x, y, z, g_1, g_o in zip(k, Section.SolutionOxide.ConcentrationH, P_H2, gamma_1, gamma_oxidation)]
        
        Activity.append(q)
    
    [ActivityFe2, ActivityFeOH, ActivityFeOH2, ActivityFeOH3, ActivityFeOH3_Fe3, ActivityFeOH4_Fe3] = Activity
    
    FeTotalActivity = [
        x + y + z + q + w + t for x, y, z, q, t, w in zip (
            ActivityFe2, ActivityFeOH, ActivityFeOH2, ActivityFeOH3, ActivityFeOH3_Fe3, ActivityFeOH4_Fe3
            )
        ]
    
    return FeTotalActivity


def iron_solubility_SB(Section):
    
    Section.SolutionOxide.ConcentrationH = bulkpH_calculator(Section)
    
    # based on high temperature henry's law constant data for H2 (some taken from T&L ref. some from sol. program DHL
    k_H2 = [-4.1991 * i + 2633.2 for i in Section.PrimaryBulkTemperature]

    # convert H2 conc. from cm^3/kg to mol/kg then concentration to P ("fugacity") using Henry's constant
    P_H2 = [(nc.H2 * nc.H2Density / nc.H2MolarMass) * i for i in k_H2]

    # no ferric (Fe3+) species considered by Sweeton and Baes
    # [cal/mol K], but constants are unitless
    A = [-26876, -11733, 4615, 9045] 
    B = [9.81, 3.35, 0, 0]
    D = [-81.21, -42.72, -23.57, -49.37]
    
    Fe3O4_SolubilityConstants_SB = []
    
    R = nc.R / 4.184 # [J/k mol to cal/K mol]

    for a, b, d in zip(A, B, D):
        k = [np.exp(((-a / i) + b * (np.log(i) - 1) + d) / R) for i in Section.PrimaryBulkTemperature]
        Fe3O4_SolubilityConstants_SB.append(k)
        
    k_Fe2 = Fe3O4_SolubilityConstants_SB[0]
    k_FeOH = Fe3O4_SolubilityConstants_SB[1]
    k_FeOH2 = Fe3O4_SolubilityConstants_SB[2]
    k_FeOH3 = Fe3O4_SolubilityConstants_SB[3]
    
    Solubility = []
    
    b = [0, 1, 2, 3]
    Constants = [k_Fe2, k_FeOH, k_FeOH2, k_FeOH3]
    
    for hydroxyls, k in zip(b, Constants):
        q = [
            x * (y ** (1 / 3)) * z ** (2 - hydroxyls) for x, y, z in zip(k, P_H2, Section.SolutionOxide.ConcentrationH)
            ]
        
        Solubility.append(q)
    
    [SolubilityFe2, SolubilityFeOH, SolubilityFeOH2, SolubilityFeOH3] = Solubility
    
    FeTotalSolubility = [
        x + y + z + w for x, y, z, w in zip(SolubilityFe2, SolubilityFeOH, SolubilityFeOH2, SolubilityFeOH3)
        ]
    # mol/kg
    return FeTotalSolubility


def cobalt_composition(Section):
    if Section in ld.SteamGenerator or Section in ld.SteamGenerator_2:
        CompositionCr_Alloy = nc.FractionCr_Alloy800
        # 5/3 term comes from lattice energy preference for Co chromite retention
        MolesCobalt = (nc.FractionCo_Alloy800 / (5 / 3)) / nc.CoMolarMass  
        # x + y = 0.0015;  x/y = 2/3 --> y =(3/2)x
        # 3/2x + x = 0.00015 --> 5/2x = 0.00015

    elif Section in ld.InletSections or Section in ld.OutletSections:
        CompositionCr_Alloy = nc.FractionCr_CS
        MolesCobalt = (nc.FractionCo_CS / (5 / 2)) / nc.CoMolarMass

    MolesChromium = CompositionCr_Alloy / nc.CrMolarMass
    return 2 / (MolesChromium / MolesCobalt)


def fraction_chromite(Section):
    if Section in ld.SteamGenerator or Section in ld.SteamGenerator_2:
        AlloyDensity = nc.Alloy800Density
        CompositionCr_Alloy = nc.FractionCr_Alloy800
        # MolesCobalt = (nc.FractionCo_Alloy800 / (5 / 3)) / nc.CoMolarMass #5/3 term comes from lattice energy 
        # preference for Co chromite retention
        # x + y = 0.00006;  x/y = 2/3 --> y =(3/2)x
        # 3/2x + x = 0.00006 --> 5/2x = 0.00006
    elif Section in ld.InletSections or Section in ld.OutletSections:
        AlloyDensity = nc.FeDensity
        CompositionCr_Alloy = nc.FractionCr_CS
        # MolesCobalt = (nc.FractionCo_CS / (5 / 2)) / nc.CoMolarMass

        # Assumptions:
        # 1. volume inner oxide replaces volume CS lost due to corrosion -used to solve for mass oxide via substitution
        # 2. all Cr is retained inside the inner oxide layer as FeCr2O4 (for both CS and Alloy-800)
        # 3. 1 g basis of metal lost

    if Section not in ld.FuelSections:
        VolumeAlloyCorroded = 1 / AlloyDensity
        MolesChromium = CompositionCr_Alloy / nc.CrMolarMass
        MolesChromite = MolesChromium / 2  # 1:2 stoichiometry in FeCr2O4 b/w compound and Cr
        MassChromite = MolesChromite * (
            cobalt_composition(Section) * nc.CoMolarMass
            + 2 * nc.CrMolarMass
            + 4 * nc.OMolarMass
            + (1 - cobalt_composition(Section)) * nc.FeMolarMass
            ) 
        VolumeChromite = MassChromite / nc.FeCr2O4Density

    return VolumeChromite / VolumeAlloyCorroded


def fraction_metal_inner_oxide(Section, Element):
    # "Second Oxide" refers to the secondary species comprising the inner layer in addition to the FeCr2O4 iron chromite
    # E.g., for carbon steel this is magnetite, and for alloy-800 it is a non-stoichiometric nickel ferrite
    if Section in ld.SteamGenerator or Section in ld.SteamGenerator_2:
        # Second oxide = Ni0.6Fe2.4O4
        FractionFeSecondOxide = 2.4 * nc.FeMolarMass / (2.4 * nc.FeMolarMass + 0.6 * nc.NiMolarMass + 4 * nc.OMolarMass)
        FractionNiSecondOxide = 0.6 * nc.NiMolarMass / (2.4 * nc.FeMolarMass + 0.6 * nc.NiMolarMass + 4 * nc.OMolarMass)
        FractionCrSecondOxide = 0
        # assumed that amount based on lattice energies retained in inner layer, while rest diffuses
        FractionCoSecondOxide = 0

    elif Section in ld.InletSections or Section in ld.OutletSections:  # Second oxide = Fe3O4
        FractionFeSecondOxide = 3 * nc.FeMolarMass / (3 * nc.FeMolarMass + 4 * nc.OMolarMass)
        FractionCrSecondOxide = 0
        FractionCoSecondOxide = 0
        FractionNiSecondOxide = 0

    if Section not in ld.FuelSections:
        FractionSecondOxide = 1 - fraction_chromite(Section)
        FractionFeChromite = (1 - cobalt_composition(Section)) * nc.FeMolarMass / (
            cobalt_composition(Section) * nc.CoMolarMass
            + 2 * nc.CrMolarMass + 4 * nc.OMolarMass
            + (1 - cobalt_composition(Section)) * nc.FeMolarMass
            )
        
        FractionCrChromite = 2 * nc.CrMolarMass / (
            cobalt_composition(Section) * nc.CoMolarMass
            + 2 * nc.CrMolarMass + 4 * nc.OMolarMass
            + (1 - cobalt_composition(Section)) * nc.FeMolarMass
            )
        
        FractionCoChromite = cobalt_composition(Section) * nc.CoMolarMass / (
            cobalt_composition(Section) * nc.CoMolarMass
            + 2 * nc.CrMolarMass + 4 * nc.OMolarMass
            + (1 - cobalt_composition(Section)) * nc.FeMolarMass
            )
        FractionNiChromite = 0 

        if Element == "Fe":
            return FractionFeChromite * fraction_chromite(Section) + FractionSecondOxide * FractionFeSecondOxide
        elif Element == "Ni":
            return FractionNiChromite * fraction_chromite(Section) + FractionSecondOxide * FractionNiSecondOxide
        elif Element == "Cr":
            return FractionCrChromite * fraction_chromite(Section)
        elif Element == "Co":
            return FractionCoChromite * fraction_chromite(Section)
        else:
            return None

def fraction_metal_in_oxide(Section, Element, Oxide):
    if Oxide == "Magnetite":
        if Element == "Fe":
            return (3 * nc.FeMolarMass) / (3 * nc.FeMolarMass + 4 * nc.OMolarMass)
    elif Oxide == "Nickel Ferrite":  # Ni0.6Fe2.4O4
        if Element == "Fe":
            return (2.4 * nc.FeMolarMass) / (2.4 * nc.FeMolarMass + 0.6 * nc.NiMolarMass + 4 * nc.OMolarMass)
        if Element == "Ni":
            return (0.6 * nc.NiMolarMass) / (2.4 * nc.FeMolarMass + 0.6 * nc.NiMolarMass + 4 * nc.OMolarMass)
    elif Oxide == "Cobalt Nickel Ferrite":  # Co0.24Ni0.22Fe2.54O4
        if Element == "Fe":
            return (2.54 * nc.FeMolarMass) / (0.24 * nc.CoMolarMass
                                              + 0.22 * nc.NiMolarMass
                                              + 2.54 * nc.FeMolarMass
                                              + 4 * nc.OMolarMass)
        if Element == "Ni":
            return  (0.22 * nc.NiMolarMass) / (0.24 * nc.CoMolarMass
                                               + 0.22 * nc.NiMolarMass
                                               + 2.54 * nc.FeMolarMass
                                               + 4 * nc.OMolarMass)
        if Element == "Co":
            return (0.24 * nc.CoMolarMass) / (0.24 * nc.CoMolarMass
                                              + 0.22 * nc.NiMolarMass
                                              + 2.54 * nc.FeMolarMass
                                              + 4 * nc.OMolarMass)
    elif Oxide == "Cobalt Ferrite":  # CoFe2O4
        if Element == "Fe":
            return (2 * nc.FeMolarMass) / (nc.CoMolarMass + 2 * nc.FeMolarMass + 4 * nc.OMolarMass)
        if Element == "Co":
            return (nc.CoMolarMass) / (nc.CoMolarMass + 2 * nc.FeMolarMass + 4 * nc.OMolarMass)
    
    elif Oxide == "Iron Chromite":  # FeCr2O4
        if Element == "Fe":
            return (nc.FeMolarMass * (1 - cobalt_composition(Section))) / (
                cobalt_composition(Section) * nc.CoMolarMass
                + 2 * nc.CrMolarMass
                + 4 * nc.OMolarMass
                + (1 - cobalt_composition(Section)) * nc.FeMolarMass
                )
        if Element == "Cr":
            return (2 * nc.CrMolarMass) / (cobalt_composition(Section) * nc.CoMolarMass
                                           + 2 * nc.CrMolarMass
                                           + 4 * nc.OMolarMass
                                           + (1 - cobalt_composition(Section)) * nc.FeMolarMass)
        if Element == "Co":
            return cobalt_composition(Section) * nc.CoMolarMass / (cobalt_composition(Section) * nc.CoMolarMass
                                                                  + 2 * nc.CrMolarMass
                                                                  + 4 * nc.OMolarMass
                                                                  + (1 - cobalt_composition(Section)) * nc.FeMolarMass)
    elif Oxide == "Nickel Chromite":  # NiCr2O4
        if Element == "Ni":
            return (nc.NiMolarMass) / (nc.CrMolarMass * 2 + nc.NiMolarMass + 4 * nc.OMolarMass)
        if Element == "Cr":
            return (2 * nc.CrMolarMass) / (nc.CrMolarMass * 2 + nc.NiMolarMass + 4 * nc.OMolarMass)
    elif Oxide == "Nickel":
        return 1
    else:
        return None
