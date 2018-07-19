import lepreau_data as ld
import iteration as it
import numpy as np
import thermochemistry_and_constants as nc
import rk_4
    
# Molecular weights of specific isotopes [g/mol]
MOLARMASSFe59 = 58.9348755
MOLARMASSFe55 = 54.9382934
MOLARMASSCo58 = 57.9357528
MOLARMASSNi63 = 62.9296694
MOLARMASSMn54 = 53.9403589
MOLARMASSCr51 = 50.9447674
MOLARMASSCo60 = 59.9338171

# Abundances (as fractions) [unitless]
ABUNDANCEFe58 = 0.00282
ABUNDANCEFe54 = 0.0584
ABUNDANCENi58 = 0.681
ABUNDANCENi62 = 0.0363
ABUNDANCECr50 = 0.0435
ABUNDANCECo59 = 1

# Decay thermochemistry_and_constants [h^-1]
LAMBDAFe59 = 0.000650476
LAMBDAFe55 = 0.0000288782
LAMBDACo58 = 0.00040735
LAMBDANi63 = 0.000000790473
LAMBDAMn54 = 0.0000924788
LAMBDACr51 = 0.00104264
LAMBDACo60 = 0.000015011

# Cross Sections [cm^2]
CROSS_SECTIONFe58 = 1.14E-24
CROSS_SECTIONFe54 = 2.5E-24
CROSS_SECTIONNi58 = 1.46E-25
CROSS_SECTIONNi62 = 1.46E-23
CROSS_SECTIONCr50 = 1.59E-23
CROSS_SECTIONCo59 = 3.7E-23

# General nuclear thermochemistry_and_constants
NEUTRON_FLUX = 50000000000000  # [neutrons/cm^2 s]
AVOGADRO = 6.022E+23  # [atoms/mol]

# Kinetic deposition/release thermochemistry_and_constants [kg_coolant/m^2 s] (All from Burrill paper)
PARTICULATE_DISSOLUTION = 0.0000018  # [s^-1]

# In-core deposition thermochemistry_and_constants
VAPORIZATION_COEFFICIENT = 0.1
ENRICHMENT = 1

def eta(Section):
    # eta term is for out-of-core activation only
    Diffusion_innerox = it.Diffusion(Section, "Fe")
    MassInnerOxide = [Diffusion_innerox * x / Section.FractionFeInnerOxide for x in Section.CorrRate]
    
    ConvertedSaturation = ld.UnitConverter(
            Section, "Mol per Kg", "Grams per Cm Cubed", Section.SolutionOxide.FeSatFe3O4, None, None, None,
            nc.FeMolarMass, None
            ) 
    
    if Section in ld.SteamGenerator or Section in ld.SteamGenerator_2: 
    
        DeltaTemperature = [abs(x - y) for x, y in zip(
            Section.PrimaryBulkTemperature[1:], Section.PrimaryBulkTemperature
            )]
       
        # One less Eta than total number of nodes --> deltas = 21, so add small diff b/w first node and last node outlet
        # Last two nodes (preheater plate) are isothermal, same addition as above
        DeltaTemperature.insert(0, 0.1)
        DeltaTemperature[Section.NodeNumber - 1] = 0.1 
    
        DeltaTemperature_Length = [x / y for x, y in zip(DeltaTemperature, Section.Length.magnitude)]        
       
        DeltaSolubility = [abs(x - y) for x, y in zip(
            Section.SolutionOxide.FeSatFe3O4[1:], Section.SolutionOxide.FeSatFe3O4
            )]
        
        # One less Eta than total number of nodes --> deltas = 21, so add small diff b/w first node and last node outlet
        DeltaSolubility.insert(0, 1e-9)
        ConvertedDeltaSolubility = ld.UnitConverter(
            Section, "Mol per Kg", "Grams per Cm Cubed", DeltaSolubility, None, None, None, nc.FeMolarMass, None
            )
        
        DeltaSolubility_Temp = [x / y for x, y in zip(ConvertedDeltaSolubility, DeltaTemperature)]
        
        CRYST_o = [x / (y * z * t * u * v) for x, y, z, t, u, v in zip(
            ConvertedSaturation, [0.18] * Section.NodeNumber, Section.Diameter, Section.Velocity, DeltaSolubility_Temp,
            DeltaTemperature_Length)]
     
    DIFF_i = [
        i / (2 * nc.Fe3O4Density * nc.FeDiffusivity * nc.Fe3O4Porosity_inner * (1 - nc.Fe3O4Porosity_inner)) for i in
        Section.InnerOxLoading
        ]
          
    CRYST_i = [x / (y * z) for x, y, z in zip(ConvertedSaturation, [0.35] * Section.NodeNumber, MassInnerOxide)]
    
    if Section in ld.InletSections or Section in ld.OutletSections:
        CRYST_o = CRYST_i
            
    ActivationCoefficient = [(1 / (x + y)) + 1 / z for x, y, z in zip(DIFF_i, CRYST_i, CRYST_o)]  # [cm/s]
    
    return ActivationCoefficient


def particulate(Section, BulkCrud_0):
    if Section in ld.FuelSections:  # Section based, not node-specific 
        # Same deposition constant down length of core (assumes no temperature dependence) 
        DepositionConstant = nc.INCORE_DEPOSITION 
    else:  # isothermal (with exception of SG)
        DepositionConstant = nc.OUTCORE_DEPOSITION 
    
    # density converted from [g/cm^3] to [kg/m^3] 
    DensityCoolant = [(i / 1000) * (100 ** 3) for i in Section.DensityH2O]
    
    # Deposition constant converted from kg_coolant/m^2 s to cm/s, dividing by density coolant (kg/m^3) converting to cm
    # [kg_coolant/m^2 s]/[kg_coolant/m^3] = [m/s]*[100 cm/m] * [1/cm] = [s^-1]
    Deposition_sec = [(DepositionConstant / x) * (4 / y) * 100 for x, y in zip(DensityCoolant, Section.Diameter)] 
    
    # Erosion constant: [g/cm^2 s] --> [g/m^2 s]
    Erosion = nc.ErosionConstant * (100 ** 2)
    
    # [s^-1] * [cm] / [cm/s] = [unitless] 
    ExponentialTerm = [np.exp(-x * y / z) for x, y, z in zip(Deposition_sec, Section.Distance, Section.Velocity)]
    
    # [g/m^2 s] / [kg_coolant/m^2 s] = [g / kg_coolant]
    Concentration = [(Erosion / DepositionConstant) * (1 - x) + BulkCrud_0 * x for x in ExponentialTerm]
    
    return Concentration # [g / kg_coolant]
                    

def surface_activity(Section, BulkActivity, j, i, Isotope):
    if Section in ld.FuelSections:
        return None
    else:
        time = [j] * Section.NodeNumber   
        
        if Isotope == "Co60":
            Lambda = LAMBDACo60
        elif Isotope == "Co58":
            Lambda = LAMBDACo58  
        elif Isotope == "Fe55":
            Lambda = LAMBDAFe55
        elif Isotope == "Fe59":
            Lambda = LAMBDAFe59
        elif Isotope == "Mn54":
            Lambda = LAMBDAMn54
        elif Isotope == "Cr51":
            Lambda = LAMBDACr51
        elif Isotope == "Ni63":
            Lambda = LAMBDANi63
        else:
            None
               
        Lambda_sec = Lambda / 3600  # [s^-1], Converts from h^-1 to s^-1
        EtaTerm = eta(Section)
        
        ActivityTerm = EtaTerm[i] * BulkActivity / Lambda_sec
        ExponentialTerm = 1 - np.exp(-Lambda * j) 
        
        ActivityConcentration = ActivityTerm * ExponentialTerm # [Bq/cm^2]
        ActivityConcentration_metersquared = ActivityConcentration * (1000 * (100 ** 2))  # [mBq/m^2]
        
        # [mCi/m^2]
        CurieSurfaceActivity = ActivityConcentration / (3.7 * 10 ** 10)
         
        return CurieSurfaceActivity


def core_active_deposit(Section, j, Element, ParentAbundance, DecayConstant, CrossSection, MolarMass):
    # Section1 = Core, Section2 = Outlet
    # Only outlet is used for averaging the composition of the crud (majority of spalling here due to high velocity)
    if Section == ld.FuelChannel:
        Section2 = ld.OutletFeeder
    elif Section == ld.FuelChannel_2:
        Section2 = ld.OutletFeeder_2
#     Section2 = ld.OutletFeeder

    Composition = []
    
    # Checks which oxide layer is uppermost (i.e., if outer layer removed due to spalling)
    # calculates elemental composition, adds that node's composition   
    # to the list below until all nodes in outlet are computed
    for i in range(Section2.NodeNumber): 
        if Section2.OuterOxLoading[i] > 0:
            Outer = "yes"
            OxideType = Section2.OuterOxLoading[i]
        else:
            Outer = "no"
            if Section2.InnerIronOxLoading[i] > 0:
                OxideType = Section2.InnerIronOxLoading[i]
        
        x = rk_4.oxide_composition(
            Section2, Element, OxideType, Outer, Section2.OuterFe3O4Loading[i], Section2.NiLoading[i],
            Section2.CoLoading[i], Section2.InnerIronOxLoading[i]
            )
        
        Composition.append(x)
        
    # Average composition element in the outlet feeder oxide    
    ElementComposition = sum(Composition) / Section2.NodeNumber
    
    CoreDepositThickness = deposition(Section, j)
    # DecayConstant = [s^-1]
    PreExponentialConstant = DecayConstant * AVOGADRO * NEUTRON_FLUX * ElementComposition * ParentAbundance \
    / (MolarMass * (DecayConstant + PARTICULATE_DISSOLUTION))
    
    ActiveDeposit = [
        PreExponentialConstant * i * (1 - np.exp(-j * (DecayConstant + PARTICULATE_DISSOLUTION) * 3600)) for i in
        CoreDepositThickness
        ]
    # Activity of specific element in the in-core deposit [Bq/cm^2]
    return ActiveDeposit

 
def deposition(Section, j):
    # Can also be manually set to desired steady-state value, e.g., ~1 ppb    
    TotalParticulate = [x + y for x, y in zip(Section.BigParticulate, Section.SmallParticulate)] # [g / kg_coolant]
    
    Time_sec = j * 3600 # s
    DepositThickness = []

    # 10th node in-core is where coolant temperature > saturation = boiling 
    # Enthalpy and heat flux vary in the boiling region  
    for i in range(Section.NodeNumber):
        if Section in ld.FuelSections and i >= 10:
            if i == 10:
                HeatFlux = 0.0762  # [kJ/cm^2 s]
                Enthalpy = 1415.8  # [kJ/kg]
            elif i == 11:
                HeatFlux = 0.05692  
                Enthalpy = 1430.4  
            elif i == 12:
                HeatFlux = 0.02487  
                Enthalpy = 1437.3
            else:
                None 
            # Only time dependence (not positional)
            # [kJ / cm^2 s] / [kJ / kg] = [kg/ cm^2 s]
            BoilingDeposition = VAPORIZATION_COEFFICIENT * ENRICHMENT * HeatFlux / Enthalpy 
            # [kg/ cm^2 s] * [g/kg_coolant] / [s^-1] = [g/cm^2]
            PreExponentiaDepositionTerm = BoilingDeposition * TotalParticulate[i] / PARTICULATE_DISSOLUTION # [g/cm^2]
        
        else: # all other PHT sections or in-core, but in a non-boiling region 
            # dW_deposit/dt = DepositionConstant*C_particulate (total) - PARTICULATE_DISSOLUTION*W_deposit
            if Section in ld.FuelSections:
                # [kg_coolant/m^2 s]*[1m^2/(100cm^2)] = [kg/cm^2 s]
                DepositionConstant = nc.INCORE_DEPOSITION / (100**2)
            else:
                DepositionConstant = nc.OUTCORE_DEPOSITION / (100**2)
            # [kg_coolant/cm^2 s]*[g/kg_coolant]/[s^-1] = [g/cm^2]
            PreExponentiaDepositionTerm = DepositionConstant * TotalParticulate[i] / PARTICULATE_DISSOLUTION
        
        # all in-core nodes set to 0 deposition (fuel bundle average residence time of ~1 year) 
        if Section in ld.FuelSections and (j + 1) % 7000 == 0:
            x = 0
        else:
            x = PreExponentiaDepositionTerm * (1 - np.exp(-PARTICULATE_DISSOLUTION * Time_sec ))
        
        DepositThickness.append(x) # [g/cm^2]
    
    return DepositThickness
    

def bulk_activity(Section, BulkConcentration_o, Isotope, j, i):
    
    if Isotope == "Co60":
        Lambda = LAMBDACo60
        MolarMass = MOLARMASSCo60
        ParentAbundance = ABUNDANCECo59
        CrossSection = CROSS_SECTIONCo59
        Element = "Co"
        
    elif Isotope == "Co58":
        Lambda = LAMBDACo58
        MolarMass = MOLARMASSCo58
        ParentAbundance = ABUNDANCENi58
        CrossSection = CROSS_SECTIONNi58 
        Element = "Ni"
        
    elif Isotope == "Fe55":
        Lambda = LAMBDAFe55
        MolarMass = MOLARMASSFe55
        ParentAbundance = ABUNDANCEFe54
        CrossSection = CROSS_SECTIONFe54
        Element = "Fe"
        
    elif Isotope == "Fe59":
        Lambda = LAMBDAFe59
        MolarMass = MOLARMASSFe59
        ParentAbundance = ABUNDANCEFe58
        CrossSection = CROSS_SECTIONFe58
        Element = "Fe"
        
    elif Isotope == "Mn54":
        Lambda = LAMBDAMn54
        MolarMass = MOLARMASSMn54
        ParentAbundance = ABUNDANCEFe54
        CrossSection = CROSS_SECTIONFe54
        Element = "Fe"
        
    elif Isotope == "Cr51":
        Lambda = LAMBDACr51
        MolarMass = MOLARMASSCr51
        ParentAbundance = ABUNDANCECr50
        CrossSection = CROSS_SECTIONCr50
        Element = "Cr"
        
    elif Isotope == "Ni63":
        Lambda = LAMBDANi63
        MolarMass = MOLARMASSNi63
        ParentAbundance = ABUNDANCENi62
        CrossSection = CROSS_SECTIONNi62
        Element = "Ni"
    else:
        None
    
    Lambda_sec = Lambda / 3600 # [s ^-1]
    
    if Section in ld.FuelSections:  
        ActiveCoreDeposit = core_active_deposit(
            Section, j, Element, ParentAbundance, Lambda_sec, CrossSection, MolarMass
            )
       
        # variation with time and distance
        ExponentialTerm = np.exp(-Lambda_sec * Section.Distance[i] / Section.Velocity[i])
        
        PreExponentialReleaseTerm = 4 * PARTICULATE_DISSOLUTION * ActiveCoreDeposit[i] /\
        (Lambda_sec * Section.Diameter[i])
        
        BulkActivity = BulkConcentration_o * ExponentialTerm + PreExponentialReleaseTerm * (1 - ExponentialTerm)
        
    # out-of-core activity                
    else:
        EtaTerm = eta(Section)
        ExponentialTerm = np.exp(
            (-Lambda_sec * Section.Distance[i] / Section.Velocity[i])
                                 - (4 * Section.Distance[i] * EtaTerm[i] / (Section.Diameter[i] * Section.Velocity[i]))
                                )
#         if Section == ld.OutletFeeder:
#             print (BulkConcentration_o, Isotope)
        # Distance is not for full PHT, just from start of current PHT section (decay from BulkConcentration_o
        # of that section's first node)
        BulkActivity = BulkConcentration_o * ExponentialTerm
       
    # convert from Bq/cm^3 to Curie/cm^3 
    #CurieBulkActivity = BulkActivity  / (3.7 * (10 ** 10)) 
    
    return BulkActivity # [Ci/cm^3]

# The flow through the purification system is provided by the HT pumps.  It is taken from one inlet header on each 
# loop of the HT system and is passed through one side of a heat interchanger, a cooler, a filter and an ion exchange
# column before being returned to the inlet of an HT pump through the other side of the interchanger.
