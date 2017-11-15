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
LambdaFe59 = 0.000650476
LambdaFe55 = 0.0000288782
LambdaCo58 = 0.00040735
LambdaNi63 = 0.000000790473
LambdaMn54 = 0.0000924788
LambdaCr51 = 0.00104264
LambdaCo60 = 0.000015011

# Cross Sections [cm^2]
CrossSectionFe58 = 1.14E-24
CrossSectionFe54 = 2.5E-24
CrossSectionNi58 = 1.46E-25
CrossSectionNi62 = 1.46E-23
CrossSectionCr50 = 1.59E-23
CrossSectionCo59 = 3.7E-23

# General nuclear constants
NeutronFlux = 50000000000000  # [neutrons/cm^2 s]
Avagadro = 6.022E+23  # [atoms/mol]
CoreSurfaceArea = 2430  # [cm^2]


def Eta(Section,CorrosionRate,FeSatFe3O4,InnerOxThickness):
    
    DeltaMassInnerOxide = [x*y/z for x,y,z in zip([it.Diffusion(Section, "Fe")]*Section.NodeNumber, CorrosionRate, [Section.FractionFeInnerOxide]*Section.NodeNumber)]
    if Section == ld.SteamGenerator: #SG (this function not called for core)    
    
        DeltaTemperature = [abs(x-y) for x,y in zip(Section.PrimaryBulkTemperature[1:],Section.PrimaryBulkTemperature)]
        DeltaTemperature.insert(0, 0.0001) #One less Eta than total number of nodes --> deltas = 21, so add small diff b/w first node and last node outlet
        
        DeltaTemperature_Length = [x/y for x,y in zip(DeltaTemperature,Section.Length.magnitude)]        
        
        DeltaSolubility = ld.UnitConverter(Section, "Mol per Kg", "Grams per Cm Cubed", [abs(x-y) for x,y in zip(FeSatFe3O4[1:],FeSatFe3O4)], \
                                            None, None, None,  FeMolarMass, None)
        DeltaSolubility.insert(0,1e-12) #One less Eta than total number of nodes --> deltas = 21, so add small diff b/w first node and last node outlet
     

        DeltaSolubility_Temp = [x/y for x,y in zip(DeltaSolubility,DeltaTemperature)]
         
        CRYST_o = [x/(y*z*t*u*v) for x,y,z,t,u,v in zip(ld.UnitConverter(Section, "Mol per Kg", "Grams per Cm Cubed", FeSatFe3O4, None, None, None,  FeMolarMass, None),     \
                                                 [0.18]*Section.NodeNumber, Section.Diameter, Section.Velocity, DeltaSolubility_Temp, DeltaTemperature_Length)]
     
    DIFF_i = [i/(2*nc.Fe3O4Density*nc.FeDiffusivity*nc.Fe3O4Porosity_inner*(1-nc.Fe3O4Porosity_inner)) for i in InnerOxThickness]
          
    CRYST_i = [x/(y*z) for x,y,z in zip(ld.UnitConverter(Section, "Mol per Kg", "Grams per Cm Cubed", FeSatFe3O4, None, None, None, nc.FeMolarMass, None),     \
                                                 [0.35]*Section.NodeNumber, DeltaMassInnerOxide)]
    if Section == ld.Inlet or Section == ld.Outlet:
        CRYST_o = CRYST_i
         
    Eta = [(1/(x+y)) + 1/z for x,y,z in zip(DIFF_i, CRYST_i, CRYST_o)] # [cm/s]
        
    return Eta


def Particulate(Section, BulkCrud):
    if Section == ld.Core: #Section based, not node-specific 
        #Same deposition constant down length of core (assumes no temperature dependence) 
        DepositionConstant = [nc.Kdeposition_InCore]*Section.NodeNumber
    else: #isothermal (with exception of SG)
        DepositionConstant = [nc.Kdeposition_OutCore]*Section.NodeNumber
    
    #Deposition constant converted from kg_coolant units to cm/s by dividing by density coolant(in kg/m^3) and converting to cm
    #[kg_coolant/m^2 s]/[kg_coolant/m^3] = [m/s]*[100 cm/m] = [cm/s]*[1/cm] = [s^-1]    
    Deposition = [100*x*(4/y)/(1000*z) for x,y,z in zip(DepositionConstant, Section.Diameter, Section.DensityH2O)]
    
    #Erosion constant: [g/cm^2 s] --> [g/m^2 s]/[kg_coolant/m^2 s]
    x = [(nc.ErosionConstant*(100**2)/x)*(1-np.exp(-y*z/u)) + q*(np.exp(-y*z/u)) for x,y,z,u,q in \
         zip(DepositionConstant, Deposition, Section.Distance, Section.Velocity, BulkCrud)]
    
    return x
                    

def SurfaceActivity(Section, CorrosionRate, FeSatFe3O4, InnerOxThickness, BulkActivity, j, Isotope):
    t = [j]*Section.NodeNumber   
    
    if Isotope == "Co60":
        Lambda =  LambdaCo60
    elif Isotope == "Co58":
        Lambda =  LambdaCo58  
    elif Isotope == "Fe55":
        Lambda =  LambdaFe55
    elif Isotope == "Fe59":
        Lambda =  LambdaFe59
    elif Isotope == "Mn54":
        Lambda =  LambdaMn54
    elif Isotope == "Cr51":
        Lambda =  LambdaCr51
    elif Isotope == "Ni63":
        Lambda =  LambdaNi63
    else:
        print ("Wrong isotope input")
            
    Lambda_sec = [Lambda/3600]*Section.NodeNumber #[s^-1], Converts from h^-1 to s^-1
    ActivityTerm = [x*y/z for x,y,z in zip(Eta(Section, CorrosionRate, FeSatFe3O4, InnerOxThickness),BulkActivity, Lambda_sec)] 
    ExponentialTerm = [1-np.exp(-Lambda*j)]*Section.NodeNumber 
    ActivityConcentration = [x*y for x,y in zip(ActivityTerm,ExponentialTerm)] #[Bq/cm^2]
    ActivityConcentration_metersquared = [i*(1000*(100**2)) for i in ActivityConcentration] #[mBq/m^2]
    
    return (ld.UnitConverter(Section, "Bq", "Ci", ActivityConcentration_metersquared, None, None, None, None, None)) #[mCi/m^2]


def InCoreActiveDeposit(Section1, Section2, InnerOxThickness, OuterOxThickness, OuterFe3O4Thickness, NiThickness, CoThickness, \
                        InnerIronOxThickness, j, Element, ParentAbundance, DecayConstant, CrossSection, MolarMass, BigParticulate, SmallParticulate):
    #Section1 = Core, Section2= Outlet
    #Only outlet is used for averaging the composition of the crud (majority of spalling here due to high velocity)
    FeComposition = []
    NiComposition = []
    CoComposition = []
    CrComposition = []
    
    for i in range(Section2.NodeNumber): 
        if OuterOxThickness[i] >0:
            Outer = "yes"
            OxideType = OuterOxThickness[i]
        else:
            if InnerOxThickness[i] >0:
                Outer = "no"
                OxideType = InnerOxThickness[i]
        
        Fe = rk_4.OxideComposition(Section2, "Fe", OxideType, Outer, OuterFe3O4Thickness[i], NiThickness[i], CoThickness[i], InnerIronOxThickness[i])
        Ni = rk_4.OxideComposition(Section2, "Ni", OxideType, Outer, OuterFe3O4Thickness[i], NiThickness[i], CoThickness[i], InnerIronOxThickness[i])
        Co = rk_4.OxideComposition(Section2, "Co", OxideType, Outer, OuterFe3O4Thickness[i], NiThickness[i], CoThickness[i], InnerIronOxThickness[i])
        Cr = rk_4.OxideComposition(Section2, "Cr", OxideType, Outer, OuterFe3O4Thickness[i], NiThickness[i], CoThickness[i], InnerIronOxThickness[i])
    
    #Checks which oxide layer is uppermost (i.e., if outer layer removed due to spalling) and calculates elemental composition, adds that node's composition   
    #to the list below until all nodes in outlet are computed
        FeComposition.append(Fe)
        NiComposition.append(Ni)
        CoComposition.append(Co)
        CrComposition.append(Cr)
    #Adds together all elements in the list (e.g., [Fe, Fe, Fe+...+] (at each node and same) for Ni, Co, etc.) and divides by total number of nodes
    #Average composition of each element in the outlet feeder oxide (outer or inner at each location used, depending on growth/spalling)      
    
    if Element == "Fe":
        ElementComposition = sum(FeComposition)/Section2.NodeNumber
    if Element == "Ni":
        ElementComposition = sum(NiComposition)/Section2.NodeNumber
    if Element =="Co":
        ElementComposition = sum(CoComposition)/Section2.NodeNumber
    if Element =="Cr":
        ElementComposition = sum(CrComposition)/Section2.NodeNumber
    
    CoreDepositThickness = Deposition(Section1, BigParticulate, SmallParticulate, j)
      
    ActiveDeposit = [(DecayConstant/3600)*(nc.Avagadro*(i/MolarMass)* NeutronFlux*CrossSection*ElementComposition*ParentAbundance/ \
    ((DecayConstant/3600)+nc.Krelease))*(1-np.exp(-j*(nc.TimeIncrement/3600)* (DecayConstant + nc.Krelease*3600))) for i in CoreDepositThickness]
    
    return ActiveDeposit
    
    
def Deposition(Section, BigParticulate, SmallParticulate, j):
    #Can also be manually set to desired steady-state value, e.g., ~1 ppb    
    TotalParticulate = [x+y for x,y in zip(BigParticulate, SmallParticulate)]
    #print (TotalParticulate)
    if Section == ld.Core:
        DepositionConstant = nc.Kdeposition_InCore
        incr = nc.TimeIncrement/3600 
         
        for i in range(Section.NodeNumber):
            #10th node in-core is where coolant temperature > saturation = boiling 
            if i >=10:
                if i == 10:
                    HeatFlux = 0.0762 #[kJ/cm^2 s]
                    Enthalpy = 1415.8  #[kJ/kg]
                if i == 11:
                    HeatFlux = 0.05692 #[kJ/cm^2 s]
                    Enthalpy = 1430.4  #[kJ/kg]
                if i == 12:
                    HeatFlux =0.02487 #[kJ/cm^2 s]
                    Enthalpy = 1437.3 #[kJ/kg]
        
        #Only time dependence (not positional)
        DepositThickness = [(nc.KVap*(HeatFlux/Enthalpy)*nc.Enrichment*x/nc.Krelease)*(1-np.exp(-nc.Krelease*(j*3600*(nc.TimeIncrement/3600)))) \
                for x in TotalParticulate] #[g/cm^2]    
        #Fuel bundle average residence time of ~1 year 
        if (j+1)%(7000/incr)==0:
            DepositThickness = 0
                    
    else:
        #dW_deposit/dt = DepositionConstant*C_particulate (total) - Krelease*W_deposit
        #[kg_coolant/m^2 s]*[g/kg_coolant]/[s^-1]*[1m^2/(100cm^2)] = [g/cm^2]
        DepositionConstant = nc.Kdeposition_OutCore
        DepositThickness = [((DepositionConstant*i)/nc.Krelease)*(1-np.exp(-nc.Krelease*(j*3600*(nc.TimeIncrement/3600))))*(1/(100**2)) \
                for i in TotalParticulate] #[g/cm^2]
    
    return DepositThickness
    

def BulkActivity(Section1, Section2, BulkConcentration_o, CorrosionRate, FeSatFe3O4, InnerOxThickness, OuterOxThickness, OuterFe3O4Thickness,\
                 NiThickness, CoThickness, InnerIronOxThickness, Isotope, BigParticulate, SmallParticulate, j, i):
    
    #[s^-1], Converts from h^-1 to s^-1
    if Isotope == "Co60":
        Lambda =  LambdaCo60
        MolarMass =  MOLARMASSCo60
        ParentAbundance =  ABUNDANCECo59
        CrossSection =  CrossSectionCo59
        Element = "Co"
        
    if Isotope == "Co58":
        Lambda =  LambdaCo58
        MolarMass =  MOLARMASSCo58
        ParentAbundance =  ABUNDANCENi58
        CrossSection =  CrossSectionNi58 
        Element = "Ni"
        
    if Isotope == "Fe55":
        Lambda =  LambdaFe55
        MolarMass =  MOLARMASSFe55
        ParentAbundance =  ABUNDANCEFe54
        CrossSection =  CrossSectionFe54
        Element = "Fe"
        
    if Isotope == "Fe59":
        Lambda =  LambdaFe59
        MolarMass =  MOLARMASSFe59
        ParentAbundance =  ABUNDANCEFe58
        CrossSection =  CrossSectionFe58
        Element = "Fe"
        
    if Isotope == "Mn54":
        Lambda =  LambdaMn54
        MolarMass =  MOLARMASSMn54
        ParentAbundance =  ABUNDANCEFe54
        CrossSection =  CrossSectionFe54
        Element = "Fe"
        
    if Isotope == "Cr51":
        Lambda =  LambdaCr51
        MolarMass =  MOLARMASSCr51
        ParentAbundance =  ABUNDANCECr50
        CrossSection =  CrossSectionCr50
        Element = "Cr"
        
    if Isotope == "Ni63":
        Lambda =  LambdaNi63
        MolarMass =  MOLARMASSNi63
        ParentAbundance =  ABUNDANCENi62
        CrossSection =  CrossSectionNi62
        Element = "Ni"
    
    if Section1 == ld.Core: #Section2 = Outlet
        CoreDeposit = InCoreActiveDeposit(Section1, Section2, Section2.InnerOxThickness, Section2.OuterOxThickness, \
                        Section2.OuterFe3O4Thickness, Section2.NiThickness, Section2.CoThickness, Section2.InnerIronOxThickness, j, Element, \
                        ParentAbundance, Lambda, CrossSection, MolarMass, BigParticulate, SmallParticulate)
        
        #print (CoreDeposit)
        Activity = (4/(Section1.Diameter[i]*Lambda/3600))*nc.Krelease*CoreDeposit[i] * \
                    (1-np.exp(-(Lambda/3600)*Section1.Distance[i]/Section1.Velocity[i])) + \
                    BulkConcentration_o*np.exp(-(Lambda/3600)*Section1.Distance[i]/Section1.Velocity[i]) 
                    
    else:
        EtaTerm = Eta(Section1, CorrosionRate, FeSatFe3O4, InnerOxThickness)
        Activity = BulkConcentration_o*np.exp((-Section1.Distance[i]*Lambda/Section1.Velocity[i]) - \
                                              (4*Section1.Distance[i]*EtaTerm[i]/(Section1.Diameter[i]*Section1.Velocity[i]))) 
    
    return Activity


#     t = [j]*Section1.NodeNumber   
#     #Array of decay constants for all isotopes 
#     Lambda = [nc.LambdaCo60, nc.LambdaCo58, nc.LambdaFe55, nc.LambdaFe59, nc.LambdaMn54, nc.LambdaCr51, nc.LambdaNi63]
#     #Takes each constant in the above array and turns that into an array based on number of nodes in Section
#     #e.g., Inlet: [[Co60, Co60,Co60,Co60,Co60,Co60,Co60], [Co58,Co58,Co58,Co58,Co58,Co58,Co58],...] - array of arrays (matrix)
#     Lambda_array = np.array([[i]*Section1.NodeNumber for i in Lambda])
#     Lambda_sec = Lambda_array/3600 #[s^-1], Converts from h^-1 to s^-1
#     #Bulk concentration in first node of respective PHTS section 
#     BulkConcentration_o_array = np.array([[i]*Section1.NodeNumber for i in BulkConcentration_o])
#     Diameter_Velocity = [x*y for x,y in zip(Section1.Diameter,Section1.Velocity)]
    
#     if Section1 == ld.Core:
#         None
#     else:
#         EtaTerm = np.array(Eta(Section1, CorrosionRate, FeSatFe3O4, InnerOxThickness)) #Don't want this to be called from the core
#         return BulkConcentration_o_array*np.exp(-(Section1.Distance*Lambda_sec)/Section1.Velocity)-(4*Section1.Distance*EtaTerm/(Diameter_Velocity))#[Bq/cm^3]
# #         if Section1 == ld.Core:
#             [Fe,Ni,Co,Cr] =InCoreActiveDeposit(Section1, Section2, InnerOxThickness, OuterOxThickness, OuterFe3O4Thickness,OuterNiThickness,OuterCoThickness,InnerIronOxThickness,j)
#             
#             return (4/(Section1.Diameter*Lambda_sec))*nc.Krelease  
#       
# f =[1,2,3]
# g =np.array([[1,1,1],[2,2,2],[3,3,3]])
# q = np.array([2,2,2])
#    
# h = g*f
# print (h)
# hh = h*(1000)
# print (hh)

    