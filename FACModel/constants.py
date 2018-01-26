# General electrochemistry constants
Beta = 0.5  # Symmetry coefficient
kH2 = 0.00078  # Henry's law constant for H2 @ 298.15 K [mol/L*atm]
F = 96485.3  # Faraday's constant [C/mol]
n = 2  # number of electrons transferred
R = 8.314  # Ideal gas constant [J/K *mol]
kb = 1.38E-23  # Boltzman constant
hp = 6.626E-34  # Planck's constant [m^2*kg/s]

# Molecular weights of compounds and elements [g/mol]
FeMolarMass = 55.847
OMolarMass = 15.99
Fe3O4MolarMass = 231.541
NiMolarMass = 58.694
CrMolarMass = 51.996
CoMolarMass = 58.9
NiFe2O4MolarMass = 234.388
FeCr2O4MoladMass = 223.779
H2MolarMass = 2.016
H2OMolarMass = H2MolarMass + OMolarMass

FractionCo_Alloy800 = 0.00015
FractionCr_Alloy800 = 0.21
FractionCo_CS = 0.00006
FractionCr_CS = 0.0002
FractionCr_Stellite = 0.27
FractionCo_Stellite = 0.65

# Density [g/cm^3]
Fe3O4Density = 5.2
NiFe2O4Density = 5.368
FeCr2O4Density = 4.8
H2Density = 8.988e-5
FeDensity = 7.86
NiDensity = 8.908
Alloy800Density = 7.94
CoDensity = 8.9

Fe3O4Porosity_inner = 0.1
Fe3O4Porosity_outer = 0.3
# FeCr2O4Porosity = 0.15
Fe3O4Tortuosity = 1.8
FeCr2O4Tortuosity = 1.2

# Diffusivity coefficients [cm^2/s]
FeDiffusivity = 0.00041
NiDiffusivity = 0.00041
CoDiffusivity = 0.00041

# Kinetic precipitation/dissolution constants [cm/s]
KpFe3O4 = 0.1
KdFe3O4 = 0.1
# KpFe_Ferrite = 0.014 #same as magnetite precipitation
# KdFe_Ferrite  = 0.044 'same as magnetite dissolution

FracNi_NiFe2O4 = 0.25
FracFe_Fe3O4 = 0.723
FracCr_Fe3O4 = 0.00019  # 0.0033
FracCo_Fe3O4 = 0.000049

ConcentrationLiTotal = 0.00063095734 # [mol/L]

# Temperature constants for equilibrium/hydrolysis constants
DebyeHuckPolynomial = [3.29E-10, -1.709E-07, 0.00003315, -0.0009028, 0.5027]
KwPolynomial = [-7.226E-10, 1.32661E-06, -0.000959311, 0.32765297, -55.86334915]
KLiPolynomial = [0.00000675, -0.0048, -0.7532]

KFe2SchikorrPolynomial = [6.66667E-11, -2.128E-07, 2.52E-04, -0.142254891, 35.94096367]
KFeOHPolynomial = [-4E-10, 8.1013E-07, -0.00062963, 0.230763, -41.5545]
KFeOH2Polynomial = [-4E-10, 8.63467E-07, -0.00073731, 0.311499, -67.8248]
KFeOH3Polynomial = [-4.667E-10, 1.0496E-06, -0.000935775, 0.413186, -97.4709]
# Fe3+
KFeOH3SchikorrPolynomial = [-1.35525E-20, 5.33333E-08, -0.00010368, 0.074951307, -27.34968224]
KFeOH4SchikorrPolynomial = [2.46667E-09, -4.72693E-06, 0.003334163, -1.007239881, 89.37956761]

KNiOHPolynomial = [3.10167E-09, -5.8294E-06, 0.004039757, -1.201661204, 119.296086]
KNiOH2Polynomial = [-3.04309E-10, 6.66757E-07, -0.000580137, 0.255160282, -60.03240922]
KNiOH3Polynomial = [8.49674E-10, -1.5126E-06, 0.000945292, -0.211217025, -16.71117746]

H2 = 10  # [cm^3/kg]
SAFactor = 1.73  # surface area factor
CobaltWear = 0.0000000015  # mol/kg (spike input term)

# Activation Energies [J/mol]
ActivationEnergyFe = 263509.1314  # at 80 um/a = 266254.4854
ActivationEnergyH2onFe = 263164.1035  # at 80 um/a = 265246.0692

DissolutionActivationEnergyFe3O4 = 131838.7907
DissolutionActivationEnergyH2onFe3O4 = 285989.8292

PrecipitationActivationEnergyFe3O4 = 179876.4692
PrecipitationActivationEnergyH2onFe3O4 = 281688.6659

ActivationEnergyAlloy800 = 295843.8269
ActivationEnergyH2onAlloy800 = 284892.0487

SimulationDuration = 10  # Total runtime [h]
TimeIncrement = 3600  # s  Based on desired time step (3600 s/h for 1h time step)

# Spalling constants depend heavily on particle size distribution
OutletOuterSpallConstant = 7000
OutletInnerSpallConstant = 1000
InletOuterSpallConstant = 1.00E+18  # Different units for inlet versus outlet (different functions)
InletInnerSpallConstant = 1.00E+5
ErosionConstant = 2.4e-11  # [g/cm^2 s]

OUTCORE_DEPOSITION = 0.0045
INCORE_DEPOSITION = 0.01

NumberPluggedTubes = 8
TotalSGTubeNumber = 3542 - NumberPluggedTubes
