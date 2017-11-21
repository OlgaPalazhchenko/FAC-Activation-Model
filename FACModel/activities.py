import lepreau_data as ld
import iteration as it
import numpy as np
import constants as nc
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

# Decay constants [h^-1]
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

# General nuclear constants
NEUTRON_FLUX = 50000000000000  # [neutrons/cm^2 s]
AVOGADRO = 6.022E+23  # [atoms/mol]
CoreSurfaceArea = 2430  # [cm^2]


def eta(Section):
    
    DeltaMassInnerOxide = [x * y / z for x, y, z in zip(
        [it.Diffusion(Section, "Fe")] * Section.NodeNumber, Section.CorrRate, [Section.FractionFeInnerOxide] *\
        Section.NodeNumber
        )]
    if Section in ld.SGZones:  # SG (this function not called for core)    
    
        DeltaTemperature = [abs(x - y) for x, y in zip(
            Section.PrimaryBulkTemperature[1:], Section.PrimaryBulkTemperature
            )]
        # One less Eta than total number of nodes --> deltas = 21, so add small diff b/w first node and last node outlet
        DeltaTemperature.insert(0, 0.0001)  
        
        DeltaTemperature_Length = [x / y for x, y in zip(DeltaTemperature, Section.Length.magnitude)]        
        
        DeltaSolubility = [abs(x - y) for x, y in zip(
            Section.SolutionOxide.FeSatFe3O4[1:], Section.SolutionOxide.FeSatFe3O4
            )]
        ConvertedDeltaSolubility = ld.UnitConverter(
            Section, "Mol per Kg", "Grams per Cm Cubed", DeltaSolubility, None, None, None, FeMolarMass, None
            )
        # One less Eta than total number of nodes --> deltas = 21, so add small diff b/w first node and last node outlet
        DeltaSolubility.insert(0, 1e-12)  
     
        DeltaSolubility_Temp = [x / y for x, y in zip(ConvertedDeltaSolubility, DeltaTemperature)]
        
        ConvertedSaturation = ld.UnitConverter(
            Section, "Mol per Kg", "Grams per Cm Cubed", Section.SolutionOxide.FeSatFe3O4, None, None, None,
            FeMolarMass, None
            ) 
        CRYST_o = [x / (y * z * t * u * v) for x, y, z, t, u, v in zip(
            ConvertedSaturation, [0.18] * Section.NodeNumber, Section.Diameter, Section.Velocity, DeltaSolubility_Temp,
            DeltaTemperature_Length)]
     
    DIFF_i = [
        i / (2 * nc.Fe3O4Density * nc.FeDiffusivity * nc.Fe3O4Porosity_inner * (1 - nc.Fe3O4Porosity_inner)) for i in
        Section.InnerOxThickness
        ]
          
    CRYST_i = [x / (y * z) for x, y, z in zip(ConvertedSaturation, [0.35] * Section.NodeNumber, DeltaMassInnerOxide)]
    
    if Section == ld.Inlet or Section == ld.Outlet:
        CRYST_o = CRYST_i
         
    ActivationCoefficient = [(1 / (x + y)) + 1 / z for x, y, z in zip(DIFF_i, CRYST_i, CRYST_o)]  # [cm/s]
        
    return ActivationCoefficient


def particulate(Section, BulkCrud):
    if Section == ld.Core:  # Section based, not node-specific 
        # Same deposition constant down length of core (assumes no temperature dependence) 
        DepositionConstant = [nc.Kdeposition_InCore] * Section.NodeNumber
    else:  # isothermal (with exception of SG)
        DepositionConstant = [nc.Kdeposition_OutCore] * Section.NodeNumber
    
    # Deposition constant converted from kg_coolant units to cm/s by dividing by density coolant(in kg/m^3) and converting to cm
    # [kg_coolant/m^2 s]/[kg_coolant/m^3] = [m/s]*[100 cm/m] = [cm/s]*[1/cm] = [s^-1]    
    Deposition = [100 * x * (4 / y) / (1000 * z) for x, y, z in
                  zip(DepositionConstant, Section.Diameter, Section.DensityH2O)]
    
    # Erosion constant: [g/cm^2 s] --> [g/m^2 s]/[kg_coolant/m^2 s]
    Erosion = nc.ErosionConstant * (100 ** 2)
    
    Concentration = [(Erosion / x) * (1 - np.exp(-y * z / u)) + q * (np.exp(-y * z / u)) for x, y, z, u, q in
                     zip(DepositionConstant, Deposition, Section.Distance, Section.Velocity, BulkCrud)]
    
    return Concentration
                    

def surface_activity(Section, BulkActivity, j, Isotope):
    
    t = [j] * Section.NodeNumber   
    
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
            
    Lambda_sec = [Lambda / 3600] * Section.NodeNumber  # [s^-1], Converts from h^-1 to s^-1
    EtaTerm = eta(Section)
    
    ActivityTerm = [x * y / z for x, y, z in zip(EtaTerm, BulkActivity, Lambda_sec)] 
    ExponentialTerm = [1 - np.exp(-Lambda * j)] * Section.NodeNumber 
    
    ActivityConcentration = [x * y for x, y in zip(ActivityTerm, ExponentialTerm)]  # [Bq/cm^2]
    ActivityConcentration_metersquared = [i * (1000 * (100 ** 2)) for i in ActivityConcentration]  # [mBq/m^2]
    
    # [mCi/m^2]
    CurieSurfaceActivity = ld.UnitConverter(
        Section, "Bq", "Ci", ActivityConcentration_metersquared, None, None, None, None, None
        ) 
    return CurieSurfaceActivity


def core_active_deposit(Section, j, Element, ParentAbundance, DecayConstant, CrossSection, MolarMass):
    # Section1 = Core, Section2 = Outlet
    # Only outlet is used for averaging the composition of the crud (majority of spalling here due to high velocity)
    Section2 = ld.Outlet
    
    Composition = []
    
    for i in range(Section2.NodeNumber): 
        if Section2.OuterOxThickness[i] > 0:
            Outer = "yes"
            OxideType = Section2.OuterOxThickness[i]
        else:
            Outer = "no"
            if Section2.InnerIronOxThickness[i] > 0:
                OxideType = Section2.InnerIronOxThickness[i]
        
        if Element == "Fe":
            x = rk_4.oxide_composition(
                Section2, "Fe", OxideType, Outer, Section2.OuterFe3O4Thickness[i], Section2.NiThickness[i],
                Section2.CoThickness[i], Section2.InnerIronOxThickness[i]
                )
        elif Element == "Ni":
            x = rk_4.oxide_composition(
                Section2, "Ni", OxideType, Outer, Section2.OuterFe3O4Thickness[i], Section2.NiThickness[i],
                Section2.CoThickness[i], Section2.InnerIronOxThickness[i]
                )
        elif Element == "Co":
            x = rk_4.oxide_composition(
                Section2, "Co", OxideType, Outer, Section2.OuterFe3O4Thickness[i], Section2.NiThickness[i],
                Section2.CoThickness[i], Section2.InnerIronOxThickness[i]
                )
        elif Element == "Cr":
            x = rk_4.oxide_composition(
                Section2, "Cr", OxideType, Outer, Section2.OuterFe3O4Thickness[i], Section2.NiThickness[i],
                Section2.CoThickness[i], Section2.InnerIronOxThickness[i]
                )
        else:
            None
        
        # Checks which oxide layer is uppermost (i.e., if outer layer removed due to spalling)
        # calculates elemental composition, adds that node's composition   
        # to the list below until all nodes in outlet are computed
        Composition.append(x)
    

    # Average composition of each element in the outlet feeder oxide    
    ElementComposition = sum(Composition) / Section2.NodeNumber
    
    CoreDepositThickness = deposition(Section, Section.BigParticulate, Section.SmallParticulate, j)
      
    ActiveDeposit = [(DecayConstant / 3600) * (AVOGADRO * (i / MolarMass) * NEUTRON_FLUX * CrossSection * ElementComposition * ParentAbundance / \
    ((DecayConstant / 3600) + nc.Krelease)) * (1 - np.exp(-j * (nc.TimeIncrement / 3600) * (DecayConstant + nc.Krelease * 3600))) for i in CoreDepositThickness]
    
    return ActiveDeposit
    
    
def deposition(Section, BigParticulate, SmallParticulate, j):
    # Can also be manually set to desired steady-state value, e.g., ~1 ppb    
    TotalParticulate = [x + y for x, y in zip(BigParticulate, SmallParticulate)]
   
    if Section == ld.Core:
        DepositionConstant = nc.Kdeposition_InCore
        incr = nc.TimeIncrement / 3600 
         
        for i in range(Section.NodeNumber):
            # 10th node in-core is where coolant temperature > saturation = boiling 
            if i >= 10:
                if i == 10:
                    HeatFlux = 0.0762  # [kJ/cm^2 s]
                    Enthalpy = 1415.8  # [kJ/kg]
                if i == 11:
                    HeatFlux = 0.05692  
                    Enthalpy = 1430.4  
                if i == 12:
                    HeatFlux = 0.02487  
                    Enthalpy = 1437.3  
        
        # Only time dependence (not positional)
        BoilingDeposition = [
            nc.KVap * nc.Enrichment * HeatFlux * x / (Enthalpy * nc.Krelease) for x in TotalParticulate
            ]
        Time_sec = j * 3600 # s
        
        DepositThickness = [x * (1 - np.exp(-nc.Krelease * Time_sec )) for x in BoilingDeposition]  # [g/cm^2]    
        # Fuel bundle average residence time of ~1 year 
        if (j + 1) % (7000 / incr) == 0:
            DepositThickness = 0
                    
    else:
        # dW_deposit/dt = DepositionConstant*C_particulate (total) - Krelease*W_deposit
        # [kg_coolant/m^2 s]*[g/kg_coolant]/[s^-1]*[1m^2/(100cm^2)] = [g/cm^2]
        DepositionConstant = nc.Kdeposition_OutCore
        DepositThickness = [((DepositionConstant * i) / nc.Krelease) * (1 - np.exp(-nc.Krelease * (j * 3600 * (nc.TimeIncrement / 3600)))) * (1 / (100 ** 2)) \
                for i in TotalParticulate]  # [g/cm^2]
    
    return DepositThickness
    

def bulk_activity(Section, BulkConcentration_o, Isotope, j, i):
    
    # [s^-1], Converts from h^-1 to s^-1
    if Isotope == "Co60":
        Lambda = LAMBDACo60
        MolarMass = MOLARMASSCo60
        ParentAbundance = ABUNDANCECo59
        CrossSection = CROSS_SECTIONCo59
        Element = "Co"
        
    if Isotope == "Co58":
        Lambda = LAMBDACo58
        MolarMass = MOLARMASSCo58
        ParentAbundance = ABUNDANCENi58
        CrossSection = CROSS_SECTIONNi58 
        Element = "Ni"
        
    if Isotope == "Fe55":
        Lambda = LAMBDAFe55
        MolarMass = MOLARMASSFe55
        ParentAbundance = ABUNDANCEFe54
        CrossSection = CROSS_SECTIONFe54
        Element = "Fe"
        
    if Isotope == "Fe59":
        Lambda = LAMBDAFe59
        MolarMass = MOLARMASSFe59
        ParentAbundance = ABUNDANCEFe58
        CrossSection = CROSS_SECTIONFe58
        Element = "Fe"
        
    if Isotope == "Mn54":
        Lambda = LAMBDAMn54
        MolarMass = MOLARMASSMn54
        ParentAbundance = ABUNDANCEFe54
        CrossSection = CROSS_SECTIONFe54
        Element = "Fe"
        
    if Isotope == "Cr51":
        Lambda = LAMBDACr51
        MolarMass = MOLARMASSCr51
        ParentAbundance = ABUNDANCECr50
        CrossSection = CROSS_SECTIONCr50
        Element = "Cr"
        
    if Isotope == "Ni63":
        Lambda = LAMBDANi63
        MolarMass = MOLARMASSNi63
        ParentAbundance = ABUNDANCENi62
        CrossSection = CROSS_SECTIONNi62
        Element = "Ni"
    
    Lambda_sec = Lambda / 3600 # s ^-1
    if Section == ld.Core:  
        CoreDeposit = core_active_deposit(Section, j, Element, ParentAbundance, Lambda, CrossSection, MolarMass)
        
        # print (CoreDeposit)
        Activity = (4 / (Section.Diameter[i] * Lambda_sec)) * nc.Krelease * CoreDeposit[i] * \
                    (1 - np.exp(-Lambda_sec * Section.Distance[i] / Section.Velocity[i])) + \
                    BulkConcentration_o * np.exp(-Lambda_sec * Section.Distance[i] / Section.Velocity[i]) 
                    
    else:
        EtaTerm = eta(Section)
        
        Activity = BulkConcentration_o * np.exp(
            (-Section.Distance[i] * Lambda_sec / Section.Velocity[i])
            - (4 * Section.Distance[i] * EtaTerm[i] / (Section.Diameter[i] * Section.Velocity[i]))
            ) 
    
    return Activity


#     t = [j]*Section1.NodeNumber   
#     #Array of decay constants for all isotopes 
#     LAMBDA = [nc.LAMBDACo60, nc.LAMBDACo58, nc.LAMBDAFe55, nc.LAMBDAFe59, nc.LAMBDAMn54, nc.LAMBDACr51, nc.LAMBDANi63]
#     #Takes each constant in the above array and turns that into an array based on number of nodes in Section
#     #e.g., Inlet: [[Co60, Co60,Co60,Co60,Co60,Co60,Co60], [Co58,Co58,Co58,Co58,Co58,Co58,Co58],...] - array of arrays (matrix)
#     LAMBDA_array = np.array([[i]*Section1.NodeNumber for i in LAMBDA])
#     LAMBDA_sec = LAMBDA_array/3600 #[s^-1], Converts from h^-1 to s^-1
#     #Bulk concentration in first node of respective PHTS section 
#     BulkConcentration_o_array = np.array([[i]*Section1.NodeNumber for i in BulkConcentration_o])
#     Diameter_Velocity = [x*y for x,y in zip(Section1.Diameter,Section1.Velocity)]
    
#     if Section1 == ld.Core:
#         None
#     else:
#         EtaTerm = np.array(Eta(Section1, CorrosionRate, FeSatFe3O4, InnerOxThickness)) #Don't want this to be called from the core
#         return BulkConcentration_o_array*np.exp(-(Section1.Distance*LAMBDA_sec)/Section1.Velocity)-(4*Section1.Distance*EtaTerm/(Diameter_Velocity))#[Bq/cm^3]
# #         if Section1 == ld.Core:
#             [Fe,Ni,Co,Cr] =InCoreActiveDeposit(Section1, Section2, InnerOxThickness, OuterOxThickness, OuterFe3O4Thickness,OuterNiThickness,OuterCoThickness,InnerIronOxThickness,j)
#             
#             return (4/(Section1.Diameter*LAMBDA_sec))*nc.Krelease  
#       
# f =[1,2,3]
# g =np.array([[1,1,1],[2,2,2],[3,3,3]])
# q = np.array([2,2,2])
#    
# h = g*f
# print (h)
# hh = h*(1000)
# print (hh)

    
