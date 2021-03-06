#General electrochemistry constants
Beta = 0.5 #Symmetry coefficient
kH2 = 0.00078 #Henry's law constant for H2 @ 298.15 K [mol/L*atm]
F = 96485.3   #Faraday's constant [C/mol]
n = 2 #number of electrons transferred 
R = 8.314 #Ideal gas constant [J/K *mol]
kb = 1.38E-23 #Boltzman constant
hp = 6.626E-34 #Planck's constant [m^2*kg/s]

#Molecular weights of compounds and elements [g/mol]
FeMolarMass = 55.847 
OMolarMass = 15.99 
Fe3O4MolarMass = 231.541 
NiMolarMass = 58.694 
CrMolarMass = 51.996 
CoMolarMass = 58.9
NiFe2O4MolarMass = 234.388 
FeCr2O4MoladMass = 223.779  
H2MolarMass = 2.016 

FractionCo_Alloy800 =0.00015
FractionCr_Alloy800 = 0.21
FractionCo_CS = 0.00006
FractionCr_CS  = 0.0002

#Density [g/cm^3]
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
#FeCr2O4Porosity = 0.15
Fe3O4Tortuosity = 1.8
FeCr2O4Tortuosity = 1.2

#Diffusivity coefficients [cm^2/s]
FeDiffusivity = 0.00041
NiDiffusivity = 0.00041
CoDiffusivity = 0.00041
    
#Molecular weights of specific isotopes [g/mol]
MolarMassFe59 = 58.9348755
MolarMassFe55 = 54.9382934
MolarMassCo58 = 57.9357528
MolarMassNi63 = 62.9296694
MolarMassMn54 = 53.9403589
MolarMassCr51 = 50.9447674
MolarMassCo60 = 59.9338171

#Abundances (as fractions) [unitless]
AbundanceFe58 = 0.00282
AbundanceFe54 = 0.0584
AbundanceNi58 = 0.681
AbundanceNi62 = 0.0363
AbundanceCr50 = 0.0435
AbundanceCo59 = 1

#Decay constants [h^-1]
LambdaFe59 = 0.000650476
LambdaFe55 = 0.0000288782
LambdaCo58 = 0.00040735
LambdaNi63 = 0.000000790473
LambdaMn54 = 0.0000924788
LambdaCr51 = 0.00104264
LambdaCo60 = 0.000015011

#Cross Sections [cm^2]
CrossSectionFe58 = 1.14E-24
CrossSectionFe54 = 2.5E-24
CrossSectionNi58 = 1.46E-25
CrossSectionNi62 = 1.46E-23
CrossSectionCr50 = 1.59E-23
CrossSectionCo59 = 3.7E-23
 
#General nuclear constants
NeutronFlux = 50000000000000        #[neutrons/cm^2 s]
Avagadro = 6.022E+23                 #[atoms/mol]
CoreSurfaceArea = 2430                 #[cm^2]

#Kinetic precipitation/dissolution constants [cm/s]
KpFe3O4  = 0.01#0.004
KdFe3O4  = 0.01#0.008
#KpFe_Ferrite = 0.014 #same as magnetite precipitation
#KdFe_Ferrite  = 0.044 'same as magnetite dissolution
    
#Kinetic deposition/release constants [kg_coolant/m^2 s] (All from Burrill paper)
Kdeposition_OutCore = 0.0045
Kdeposition_InCore = 0.01
Krelease = 0.0000018 #[s^-1]

FracNi_NiFe2O4 = 0.25
FracFe_Fe3O4 = 0.723
FracCr_Fe3O4 = 0.00019 #0.0033
FracCo_Fe3O4 = 0.000049

ConcentrationLiTotal=0.00022585 #[mol/L]

#Temperature constants for equilibrium/hydrolysis constants
DebyeHuckPolynomial=[3.29E-10, -1.709E-07, 0.00003315, -0.0009028, 0.5027]
KwPolynomial= [-7.226E-10, 1.32661E-06, -0.000959311, 0.32765297,-55.86334915]
KLiPolynomial= [0.00000675, -0.0048, -0.7532]

KFeOHPolynomial= [-4E-10, 8.1013E-07, -0.00062963, 0.230763,-41.5545]
KFeOH2Polynomial= [-4E-10,  8.63467E-07, -0.00073731,  0.311499,-67.8248]
KFeOH3Polynomial= [-4.667E-10, 1.0496E-06, -0.000935775,  0.413186,-97.4709]
  
KNiOHPolynomial= [3.10167E-09, -5.8294E-06, 0.004039757, -1.201661204, 119.296086]
KNiOH2Polynomial = [-3.04309E-10, 6.66757E-07, -0.000580137, 0.255160282, -60.03240922]          
KNiOH3Polynomial= [8.49674E-10, -1.5126E-06, 0.000945292, -0.211217025, -16.71117746]               

H2 = 10 #[cm^3/kg]
SAFactor = 1.73 #surface area factor 
CobaltWear = 0.0000000015 #mol/kg (spike input term)

#Activation Energies [J/mol]
ActivationEnergyFe = 264578.5708
ActivationEnergyH2onFe = 263220.6351
DissolutionActivationEnergyFe3O4 = 144987.3286
DissolutionActivationEnergyH2onFe3O4 = 297150.6019
PrecipitationActivationEnergyFe3O4 = 191037.2419
PrecipitationActivationEnergyH2onFe3O4 = 291214.6173
ActivationEnergyAlloy800 = 294051.6401
ActivationEnergyH2onAlloy800 = 283099.8619

SimulationDuration = 10 #Total runtime
TimeIncrement = 3600 #s  Based on desired time step (3600 s/h for 1h time step)
#Spalling constants depend heavily on particle size distribution 
OutletOuterSpallConstant = 7000#10000
OutletInnerSpallConstant = 3000#100
InletOuterSpallConstant = 1.00E+18 #Different units for 
InletInnerSpallConstant = 1.00E+5
ErosionConstant = 2.4e-10 #[g/cm^2 s]