import constants as nc
import lepreau_data as ld
import numpy as np


def temperature_dependent_parameters(Section):
    Celsius = [i - 273.15 for i in Section.PrimaryBulkTemperature]
    # Equilibrium and Debye-Huckel constants - polynomials as a function of temperature
    # Coeff1*x^4 + Coeff2*x^3 + Coeff3*x^2 + Coeff4*x + Coeff5, depends on # elements in coeff list
    Section.DebyeHuckelConstant = (np.polyval(nc.DebyeHuckPolynomial, Celsius))
    Section.k_W = 10 ** (np.polyval(nc.KwPolynomial, Section.PrimaryBulkTemperature))
    Section.k_Li = 10 ** (np.polyval(nc.KLiPolynomial, Section.PrimaryBulkTemperature))
    
    Gibbs_Energies = [-62700, -4300, 56000, 127840, 74000, 128800] # [J/mol]
    Entropies = [-51.63, -108.96, -102.35, -188.71, -66.60, -187.76] # [J/ K mol]
    Delta_T = [i - 298 for i in Section.PrimaryBulkTemperature] # [K]
    
    Fe3O4_SolubilityConstants = []
     
    for x, y, z in zip(Gibbs_Energies, Entropies, Delta_T):
        q = [np.exp(-(x - y * z) / (nc.R * i)) for i in Section.PrimaryBulkTemperature]
        Fe3O4_SolubilityConstants.append(q)
        
    Section.k_Fe2 = Fe3O4_SolubilityConstants[0]
    Section.k_FeOH = Fe3O4_SolubilityConstants[1]
    Section.k_FeOH2 = Fe3O4_SolubilityConstants[2]
    Section.k_FeOH3 = Fe3O4_SolubilityConstants[3]
    Section.k_FeOH3_Fe3 = Fe3O4_SolubilityConstants[4]
    Section.k_FeOH4_Fe3 = Fe3O4_SolubilityConstants[5]
    
    return None
    
#     Section.k_FeOH = 10 ** (np.polyval(nc.KFeOHPolynomial, Section.PrimaryBulkTemperature))
#     Section.k_FeOH2 = 10 ** (np.polyval(nc.KFeOH2Polynomial, Section.PrimaryBulkTemperature))
#     Section.k_FeOH3 = 10 ** (np.polyval(nc.KFeOH3Polynomial, Section.PrimaryBulkTemperature))
#     Section.k_NiOH = 10 ** (np.polyval(nc.KNiOHPolynomial, Section.PrimaryBulkTemperature))
#     Section.k_NiOH2 = 10 ** (np.polyval(nc.KNiOH2Polynomial, Section.PrimaryBulkTemperature))
#     Section.k_NiOH3 = 10 ** (np.polyval(nc.KNiOH3Polynomial, Section.PrimaryBulkTemperature))


def bulkpH_calculator(Section):  # Bulk pH calculation
    temperature_dependent_parameters(Section)
    H = []
    ConcentrationH = 0.0000000001  # initial guess [mol/kg]
    gamma_1 = 1  # initial guess
    for i in range(Section.NodeNumber):
        ConcentrationOH = Section.k_W[i] / ((gamma_1 ** 2) * (ConcentrationH))

        # At high temp, LiOH doesn't dissociate 100% - eq'm established: LiOH(aq) <-> Li+(aq) + OH-(aq)
        # KLi = ([Li+]gamma_1 *[OH-]gamma_1)/[LiOH]; 
        # ConcentrationLiTotal = ConcentrationLi  + LiConcentrationOH (sub in eq'm expression for LiConcentrationOH
        ConcentrationLi = Section.k_Li[i] * nc.ConcentrationLiTotal / (
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
            ConcentrationLi = (Section.k_Li[i] * nc.ConcentrationLiTotal) / (
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


def iron_solubility(Section):
    # As temperature changes, constants (k_Li, k_w, etc. change), pH changes as well, 
    # both effecting solubility --> Bulk FeSatFe3O4
    Section.SolutionOxide.ConcentrationH = bulkpH_calculator(Section)
    
    # based on high temperature henry's law constant data for H2 (some taken from T&L ref. some from sol. program DHL
    k_H2 = [-4.1991 * i + 2633.2 for i in Section.PrimaryBulkTemperature]
    
    # convert H2 conc. from cm^3/kg to mol/kg then concentration to P ("fugacity") using Henry's constant
    P_H2 = [(nc.H2 * nc.H2Density / nc.H2MolarMass) * i for i in k_H2]
    
    gamma_1 = 0.97
    gamma_2 = 0.85
    
    b = [0, 1, 2, 3, 3, 4]
    
    Activity = []
    Constants = [
        Section.k_Fe2, Section.k_FeOH, Section.k_FeOH2, Section.k_FeOH3, Section.k_FeOH3_Fe3, Section.k_FeOH4_Fe3
        ]
    
    for oxidation, k in zip (b, Constants):
        if oxidation == 0:
            gamma_oxidation = gamma_2
        elif oxidation == 1 or oxidation == 3 or oxidation == 4: # -1 for Fe2+, +1 for Fe2+, -1 for Fe3+ species
            gamma_oxidation = gamma_1
        elif oxidation == 2 or oxidation == 4:
            gamma_oxidation = 1 # 0 charge on Fe(OH)2^0 and Fe(OH)3^0 = electroneutrality, no gamma
        else: 
            None
        
        if k == Constants[4] or k == Constants[5]:
            charge = 3 # ferric, Fe3+ species
        else:
            charge = 2 # ferrous, Fe2+ species

        q = [x * ((y * gamma_1)**(charge - oxidation)) * (z**((4/3) - (charge / 2))) / (gamma_oxidation) for x, y, z in zip(
            k, Section.SolutionOxide.ConcentrationH, P_H2)]
        
        Activity.append(q)
    
    [ActivityFe2, ActivityFeOH, ActivityFeOH2, ActivityFeOH3, ActivityFeOH3_Fe3, ActivityFeOH4_Fe3] = Activity

    FeTotalActivity = [
        x + y + z + q + w + t for x, y, z, q, t, w in zip (
            ActivityFe2, ActivityFeOH, ActivityFeOH2, ActivityFeOH3, ActivityFeOH3_Fe3, ActivityFeOH4_Fe3
            )
        ]
    
    if Section == ld.FuelChannel:
        print (Section.PrimaryBulkTemperature, -np.log10(Section.SolutionOxide.ConcentrationH), FeTotalActivity)
    if Section == ld.OutletFeeder:
        print (Section.PrimaryBulkTemperature, -np.log10(Section.SolutionOxide.ConcentrationH), FeTotalActivity, "outlet")
    return FeTotalActivity


def hydrolysis(Section, FeTotal, NiTotal, ConcentrationH):
    ActivityCoefficient1 = [1.0] * Section.NodeNumber  # initial estimate for activities (+/- 1 charged ions)
    ActivityCoefficient2 = [1.0] * Section.NodeNumber  # (+/- 2 charged ions)
    # Concentration NiTotal <<FeTotal, so Ni species left out of ionic strength calcs
    for i in range(50):
        gamma_1itr = ActivityCoefficient1
        gamma_2itr = ActivityCoefficient2

        # ConcentrationFe/NiTotal and ConcentrationH are not re-evaluated here, thus  not used for them
        ConcentrationOH = [x / (y * (z ** 2)) for x, y, z in zip(Section.k_W, ConcentrationH, ActivityCoefficient1)]
        ConcentrationLi = [nc.ConcentrationLiTotal * x / (x + (y * (z ** 2))) for x, y, z in zip(
            Section.k_Li, ConcentrationOH, ActivityCoefficient1
            )]
        
        ConcentrationFe2 = [x / (1
                                 +(y * z / (r * (s ** 2)))
                                 + (t * z / ((r ** 2) * (s ** 2)))
                                 + (q * z / ((r ** 3) * (s ** 4)))) 
            for x, y, z, r, s, t, q in zip(
                FeTotal, Section.k_FeOH, ActivityCoefficient2, ConcentrationH, ActivityCoefficient1, Section.k_FeOH2, 
                Section.k_FeOH3
                )]
        ConcentrationFeOH = [x * y * z / (q * (r ** 2)) for x, y, z, q, r in zip(
            Section.k_FeOH, ActivityCoefficient2, ConcentrationFe2, ConcentrationH, ActivityCoefficient1
            )]
        ConcentrationFeOH2 = [x * y * z / ((q ** 2) * (r ** 2)) for x, y, z, q, r in zip(
            Section.k_FeOH2, ActivityCoefficient2, ConcentrationFe2, ConcentrationH, ActivityCoefficient1
            )] 
        ConcentrationFeOH3 = [x * y * z / ((q ** 3) * (r ** 4)) for x, y, z, q, r, in zip(
            Section.k_FeOH3, ActivityCoefficient2, ConcentrationFe2, ConcentrationH, ActivityCoefficient1
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
        # print (gamma_1itr, ActivityCoefficient1,"lol")
        RE2 = [((x - y) / x) for x, y in zip(gamma_2itr, ActivityCoefficient2)]
        if RE < [0.00000001] * Section.NodeNumber and RE2 < [0.00000001] * Section.NodeNumber:
            # print (i)
            # break
            return ConcentrationFe2, ConcentrationFeOH2, ActivityCoefficient1, ActivityCoefficient2


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
