'''
Created on Feb 27, 2018

@author: opalazhc
'''
import numpy as np

class Section():
    def __init__(self):
        self.PrimaryBulkTemperature = None
        self.NodeNumber = None
        self.ConcentrationH = None
        self.FeTotal = None

SG = Section()

# user input
# based on pH (298 K) of 10.3
CONCENTRATION_LITHIUM = 0.000225095734 # [mol/kg]
H2 = 10  # [cm^3/kg] dissolved hydrogen gas concentration

SG.PrimaryBulkTemperature = [265 + 273.15] # [K]
SG.NodeNumber = len(SG.PrimaryBulkTemperature)

#literature constants
R = 8.314 # Ideal gas constant [J/K *mol]
H2Density = 8.988e-5
H2MolarMass = 2.016

# Temperature thermochemistry_and_constants for equilibrium/hydrolysis thermochemistry_and_constants
DEBYE_HUCKEL_POLYNOMIAL = [3.29E-10, -1.709E-07, 0.00003315, -0.0009028, 0.5027]
KwPOLYNOMIAL = [-7.226E-10, 1.32661E-06, -0.000959311, 0.32765297, -55.86334915]
KLiPOLYNOMIAL = [0.00000675, -0.0048, -0.7532]

# KFe2SchikorrPOLYNOMIAL = [6.66667E-11, -2.128E-07, 2.52E-04, -0.142254891, 35.94096367]
KFeOHPOLYNOMIAL = [-4E-10, 8.1013E-07, -0.00062963, 0.230763, -41.5545]
KFeOH2POLYNOMIAL = [-4E-10, 8.63467E-07, -0.00073731, 0.311499, -67.8248]
KFeOH3POLYNOMIAL = [-4.667E-10, 1.0496E-06, -0.000935775, 0.413186, -97.4709]


def temperature_dependent_parameters(Section):
    Celsius = [i - 273.15 for i in Section.PrimaryBulkTemperature]
    # Equilibrium and Debye-Huckel thermochemistry and_ onstants - polynomials as a function of temperature
    # Coeff1*x^4 + Coeff2*x^3 + Coeff3*x^2 + Coeff4*x + Coeff5, depends on # elements in coeff list
    Section.DebyeHuckelConstant = (np.polyval(DEBYE_HUCKEL_POLYNOMIAL, Celsius))
    Section.k_W = 10 ** (np.polyval(KwPOLYNOMIAL, Section.PrimaryBulkTemperature))
    Section.k_Li = 10 ** (np.polyval(KLiPOLYNOMIAL, Section.PrimaryBulkTemperature))
    
    Gibbs_Energies = [-62700, -4300, 56000, 127840, 74000, 128800] # [J/mol]
    Entropies = [-51.63, -108.96, -102.35, -188.71, -66.60, -187.76] # [J/ K mol]
    Delta_T = [i - 298 for i in Section.PrimaryBulkTemperature] # [K]
    
    Fe3O4_SolubilityConstants = []
    
    for x, y in zip(Gibbs_Energies, Entropies):
        q = [np.exp(-(x - y * k) / (R * i)) for i, k in zip(Section.PrimaryBulkTemperature, Delta_T)]
        Fe3O4_SolubilityConstants.append(q)

    Section.k_Fe2 = Fe3O4_SolubilityConstants[0]
    Section.k_FeOH = Fe3O4_SolubilityConstants[1]
    Section.k_FeOH2 = Fe3O4_SolubilityConstants[2]
    Section.k_FeOH3 = Fe3O4_SolubilityConstants[3]
    Section.k_FeOH3_Fe3 = Fe3O4_SolubilityConstants[4]
    Section.k_FeOH4_Fe3 = Fe3O4_SolubilityConstants[5]
    
    Section.k_FeOH_hydrolysis = 10 ** (np.polyval(KFeOHPOLYNOMIAL, Section.PrimaryBulkTemperature))
    Section.k_FeOH2_hydrolysis = 10 ** (np.polyval(KFeOH2POLYNOMIAL, Section.PrimaryBulkTemperature))
    Section.k_FeOH3_hydrolysis = 10 ** (np.polyval(KFeOH3POLYNOMIAL, Section.PrimaryBulkTemperature))

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


def hydrolysis(Section, FeTotal):
    
    Section.ConcentrationH = bulkpH_calculator(Section)
    
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
                                 +(y * z / (r * (s ** 2)))
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
        
        # FeOH2 left out of ionic strength calc due to electroneutrality principle
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


def iron_solubility(Section, Condition):
    # As temperature changes, thermochemistry and constants (k_Li, k_w, etc. change), pH changes as well, 
    # both affecting solubility --> Bulk FeSatFe3O4
    Section.ConcentrationH = bulkpH_calculator(Section)
    
    # based on high temperature henry's law constant data for H2 (some taken from T&L ref. some from sol. program DHL
    k_H2 = [-4.1991 * i + 2633.2 for i in Section.PrimaryBulkTemperature]

    # convert H2 conc. from cm^3/kg to mol/kg then concentration to P ("fugacity") using Henry's constant
    P_H2 = [(H2 * H2Density / H2MolarMass) * i for i in k_H2]
    
    # if total iron concentration in solution unknown, use these activity coeffs    
    if Condition == "initial":
        gamma_1 = [0.95] * Section.NodeNumber
        gamma_2 = [0.78] * Section.NodeNumber
    else:
        gamma_1 = hydrolysis(Section, Section.FeTotal, Section.ConcentrationH)[2]
        gamma_2 = hydrolysis(Section, Section.FeTotal, Section.ConcentrationH)[3]
           
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
            # 0 charge on Fe(OH)2^0 and Fe(OH)3^0 = electroneutrality, no gamma
            gamma_oxidation = [1] * Section.NodeNumber
        else: 
            None
        
        if k == Constants[4] or k == Constants[5]:
            charge = 3 # ferric, Fe3+ species
        else:
            charge = 2 # ferrous, Fe2+ species

        q = [x * ((y * g_1)**(charge - oxidation)) * (z**((4/3) - (charge / 2))) / (g_o) for
             x, y, z, g_1, g_o in zip(k, Section.ConcentrationH, P_H2, gamma_1, gamma_oxidation)]
        
        Activity.append(q)
    
    [ActivityFe2, ActivityFeOH, ActivityFeOH2, ActivityFeOH3, ActivityFeOH3_Fe3, ActivityFeOH4_Fe3] = Activity
    
    FeTotalActivity = [
        x + y + z + q + w + t for x, y, z, q, t, w in zip (
            ActivityFe2, ActivityFeOH, ActivityFeOH2, ActivityFeOH3, ActivityFeOH3_Fe3, ActivityFeOH4_Fe3
            )
        ]
    
    return FeTotalActivity

print (iron_solubility(SG, "initial"))