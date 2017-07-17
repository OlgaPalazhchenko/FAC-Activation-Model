import LepreauData as ld
import Iteration as it
import numpy as np
import NumericConstants as nc
import RK4
    
 
def Eta(Section,CorrosionRate,FeSatFe3O4,InnerOxThickness):
    
    DeltaMassInnerOxide = [x*y/z for x,y,z in zip([it.Diffusion(Section, "Fe")]*Section.NodeNumber, CorrosionRate, [Section.FractionFeInnerOxide]*Section.NodeNumber)]
    if Section == ld.SteamGenerator: #SG (this function not called for core)    
    
        DeltaTemperature = [abs(x-y) for x,y in zip(Section.Kelvin[1:],Section.Kelvin)]
        DeltaTemperature.insert(0, 0.0001) #One less Eta than total number of nodes --> deltas = 21, so add small diff b/w first node and last node outlet
        
        DeltaTemperature_Length = [x/y for x,y in zip(DeltaTemperature,Section.Length)]        
        
        DeltaSolubility = ld.UnitConverter(Section, "Mol per Kg", "Grams per Cm Cubed", [abs(x-y) for x,y in zip(FeSatFe3O4[1:],FeSatFe3O4)], \
                                            None, None, None, nc.FeMolarMass)
        DeltaSolubility.insert(0,1e-12) #One less Eta than total number of nodes --> deltas = 21, so add small diff b/w first node and last node outlet
     

        DeltaSolubility_Temp = [x/y for x,y in zip(DeltaSolubility,DeltaTemperature)]
         
        CRYST_o = [x/(y*z*t*u*v) for x,y,z,t,u,v in zip(ld.UnitConverter(Section, "Mol per Kg", "Grams per Cm Cubed", FeSatFe3O4, None, None, None, nc.FeMolarMass),     \
                                                 [0.18]*Section.NodeNumber, Section.Diameter, Section.Velocity, DeltaSolubility_Temp, DeltaTemperature_Length)]
     
    DIFF_i = [i/(2*nc.Fe3O4Density*nc.FeDiffusivity*nc.Fe3O4Porosity_inner*(1-nc.Fe3O4Porosity_inner)) for i in InnerOxThickness]
          
    CRYST_i = [x/(y*z) for x,y,z in zip(ld.UnitConverter(Section, "Mol per Kg", "Grams per Cm Cubed", FeSatFe3O4, None, None, None, nc.FeMolarMass),     \
                                                 [0.35]*Section.NodeNumber, DeltaMassInnerOxide)]
    if Section == ld.Inlet or Section == ld.Outlet:
        CRYST_o = CRYST_i
         
    Eta = [(1/(x+y)) + 1/z for x,y,z in zip(DIFF_i, CRYST_i, CRYST_o)] # [cm/s]
        
    return Eta


def Particulate(Section, BulkCrud):
    if Section == ld.Core: #Section based, not node-specific 
        DepositionConstant = [nc.Kdeposition_InCore]*Section.NodeNumber
    else:
        DepositionConstant = [nc.Kdeposition_OutCore]*Section.NodeNumber
    
    #Deposition constant converted from kg_coolant units to cm/s by dividing by density coolant(in kg/m^3) and converting to cm
    #[kg_coolant/m^2 s]/[kg_coolant/m^3] = [m/s]*[100 cm/m] = [cm/s]*[1/cm] = [s^-1]    
    Deposition = [100*x*(4/y)/(1000*z) for x,y,z in zip(DepositionConstant, Section.Diameter, Section.DensityH2O)]
    
    #Erosion constant: [g/cm^2 s] --> [g/m^2 s]/[kg_coolant/m^2 s]
    x = [(nc.ErosionConstant*(100**2)/x)*(1-np.exp(-y*z/u)) + q*(np.exp(-y*z/u)) for x,y,z,u,q in zip(DepositionConstant, Deposition, Section.Distance, Section.Velocity, BulkCrud)]
    
    return x
                    

def SurfaceActivity(Section, CorrosionRate, FeSatFe3O4, InnerOxThickness, BulkActivity, j):
    t = [j]*Section.NodeNumber   
    #Array of decay constants for all isotopes 
    Lambda = [nc.LambdaCo60, nc.LambdaCo58, nc.LambdaFe55, nc.LambdaFe59, nc.LambdaMn54, nc.LambdaCr51, nc.LambdaNi63]
    #Takes each constant in the above array and turns that into an array based on number of nodes in Section
    #e.g., Inlet: [[Co60, Co60,Co60,Co60,Co60,Co60,Co60], [Co58,Co58,Co58,Co58,Co58,Co58,Co58],...] - array of arrays (matrix)
    Lambda_array = np.array([[i]*Section.NodeNumber for i in Lambda])
    Lambda_sec = Lambda_array/3600 #[s^-1], Converts from h^-1 to s^-1
    EtaTerm = np.array(Eta(Section, CorrosionRate, FeSatFe3O4, InnerOxThickness))
    
    ActivityTerm = EtaTerm*BulkActivity/Lambda_sec  
    ExponentialTerm = 1-np.exp(-np.array(Lambda_array)*t) 
    ActivityConcentration = ActivityTerm*ExponentialTerm #[Bq/cm^2]
    ActivityConcentration_metersquared =ActivityConcentration*1000*(100**2) #[mBq/m^2]
     
    #Lambda_sec = [Lambda/3600]*Section.NodeNumber #[s^-1], Converts from h^-1 to s^-1
    #ActivityTerm = [x*y/z for x,y,z in zip(Eta(Section, CorrosionRate, FeSatFe3O4, InnerOxThickness),BulkActivity, Lambda_sec)] 
    #ExponentialTerm = [1-np.exp(-Lambda*j)]*Section.NodeNumber 
    #ActivityConcentration = [x*y for x,y in zip(ActivityTerm,ExponentialTerm)] #[Bq/cm^2]
    #ActivityConcentration_metersquared = [i*(1000*(100**2)) for i in ActivityConcentration] #[mBq/m^2]
    return (ld.UnitConverter(Section, "Bq", "Ci", ActivityConcentration_metersquared, None, None, None, None)) #[mCi/m^2]

def InCoreActiveDeposit(Section1, Section2, InnerOxThickness, OuterOxThickness, OuterFe3O4Thickness, NiThickness, CoThickness, InnerIronOxThickness, \
                        j):
    #Creates a blank list 
    FeComposition = []
    NiComposition = []
    CoComposition = []
    CrComposition = []
    #All Section 2 sent through here?
    for i in range(Section1.NodeNumber): #Only outlet is used for averaging the composition of the crud (majority of spalling here due to high velocity)
        if OuterOxThickness[i] >0:
            Outer = "yes"
            OxideType = OuterOxThickness[i]
        else:
            if OuterOxThickness[i] == 0:
                if InnerOxThickness[i] >0:
                    Outer = "no"
                    OxideType = InnerOxThickness[i]
        
        Fe = RK4.OxideComposition(Section2, "Fe", OxideType, Outer, OuterFe3O4Thickness[i], NiThickness[i], CoThickness[i], InnerIronOxThickness[i])
        Ni = RK4.OxideComposition(Section2, "Ni", OxideType, Outer, OuterFe3O4Thickness[i], NiThickness[i], CoThickness[i], InnerIronOxThickness[i])
        Co = RK4.OxideComposition(Section2, "Co", OxideType, Outer, OuterFe3O4Thickness[i], NiThickness[i], CoThickness[i], InnerIronOxThickness[i])
        Cr = RK4.OxideComposition(Section2, "Cr", OxideType, Outer, OuterFe3O4Thickness[i], NiThickness[i], CoThickness[i], InnerIronOxThickness[i])
    #Checks which oxide layer is uppermost (i.e., if outer layer removed due to spalling) and calculates elemental composition, adds that node's composition   
    #to the list below until all nodes in outlet are computed
        FeComposition.append(Fe)
        NiComposition.append(Ni)
        CoComposition.append(Co)
        CrComposition.append(Cr)
    #Adds together all elements in the list (e.g., [Fe, Fe, Fe+...+] (at each node and same) for Ni, Co, etc.) and divides by total number of nodes
    #Average composition of each element in the outlet feeder oxide (outer or inner at each location used, depending on growth/spalling)      
    FeAverage = sum(FeComposition)/Section2.NodeNumber
    NiAverage = sum(NiComposition)/Section2.NodeNumber
    CoAverage = sum(CoComposition)/Section2.NodeNumber
    CrAverage = sum(CrComposition)/Section2.NodeNumber
        
    #Maps the inputed isotope data: [Co60, Co58, Fe55, Fe59, Mn54, Cr51, Ni63]    
    ElementComposition = [CoAverage, CoAverage, FeAverage, FeAverage, FeAverage, CrAverage, NiAverage]
    
    
def Deposition(Section):
    KVap  = 0.1
    Enrichment = 1
    
    if Section == ld.Core:
        DepositionConstant = nc.Kdeposition_InCore
        for i in range(Section.NodeNumber):
            if i >=10:
                if i == 10:
                    HeatFlux = 0.0762 #[kJ/cm^2 s]
                    Enthalpy = 1415.8  #[kJ/kg]
                if i == 11:
                    HeatFlux = 0.05692
                    Enthalpy = 1430.4
                if i == 12:
                    HeatFlux =0.02487
                    Enthalpy = 1437.3
            
            return KVap*(HeatFlux/Enthalpy)*Enrichment
                
#             else
    else:
        DepositionConstant = nc.Kdeposition_OutCore
    
    
    
    
          
# def LOL(rofl):
#     b= []
#     for i in rofl:
#         if i == 0:
#             i=0.0001
#         b.append(i)
#         
#     c = sum(b)/4
#     print (c)
#     return b
#  
# a = LOL([1,2,0,2])
# print (a)
        
# def LOL(rofl):
#     list1= []
#     list2= []
#     for i in rofl:
#         if i >0:
#             a= 1
#             b= 3
#         else:
#             a=0
#             b=2
#         list1.append(a)
#         list2.append(b)
#         #print (list1,list2)
#     #print (list1,list2)
#     return [list1,list2]
# LOL([0,1,2,0])
# [q,r] = LOL([0,1,2,0])
# print (q)

def BulkActivity(Section1, Section2, BulkConcentration_o, CorrosionRate, FeSatFe3O4, InnerOxThickness, OuterOxThickness, OuterFe3O4Thickness,\
                 OuterNiThickness, OuterCoThickness, InnerIronOxThickness, j):
    t = [j]*Section1.NodeNumber   
    #Array of decay constants for all isotopes 
    Lambda = [nc.LambdaCo60, nc.LambdaCo58, nc.LambdaFe55, nc.LambdaFe59, nc.LambdaMn54, nc.LambdaCr51, nc.LambdaNi63]
    #Takes each constant in the above array and turns that into an array based on number of nodes in Section
    #e.g., Inlet: [[Co60, Co60,Co60,Co60,Co60,Co60,Co60], [Co58,Co58,Co58,Co58,Co58,Co58,Co58],...] - array of arrays (matrix)
    Lambda_array = np.array([[i]*Section1.NodeNumber for i in Lambda])
    Lambda_sec = Lambda_array/3600 #[s^-1], Converts from h^-1 to s^-1
    #Bulk concentration in first node of respective PHTS section 
    BulkConcentration_o_array = np.array([[i]*Section1.NodeNumber for i in BulkConcentration_o])
    Diameter_Velocity = [x*y for x,y in zip(Section1.Diameter,Section1.Velocity)]
    
    if Section1 == ld.Core:
        None
                    
    else:
        EtaTerm = np.array(Eta(Section1, CorrosionRate, FeSatFe3O4, InnerOxThickness)) #Don't want this to be called from the core
        return BulkConcentration_o_array*np.exp(-(Section1.Distance*Lambda_sec)/Section1.Velocity)-(4*Section1.Distance*EtaTerm/(Diameter_Velocity))#[Bq/cm^3]
#         if Section1 == ld.Core:
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

    