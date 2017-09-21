import LepreauData as ld
import numpy as np
import NumericConstants as nc


def ThermalConductivity(Twall, material):
    if material == "Alloy-800" or material == "Alloy800" or material == "A800" or material == "Alloy 800":
        return (11.450 + 0.0161*Twall)/100 #M.K.E thesis for Alloy-800 tubing [W/cm K]
    elif material == "water":
        return 666.4/100 #W/cm K
    else: print ("Error: material not specified")


def ConductionResistance(Section, Twall, i):
    #R_conduction = L/kA (L = thickness)
    #A = area with shape factor correcrion factor for cylindrical pipe
    #R_conduction = ln(D_o/D_i)/(2pikl) [K/W]
    #l = length, k =thermal conductivity coefficient [W/cm K]
    
    k_w = ThermalConductivity(Twall, "Alloy-800")  #[W/cm K]  
    
    Rcyl_numerator = np.log(Section.OuterDiameter[i]/Section.Diameter[i]) 
    Rcyl_denominator = 2*np.pi*k_w*Section.Length[i] 
    return Rcyl_numerator/Rcyl_denominator #[K/W]


def PrimaryConvectionResistance(Section, correlation, Twall, preheater, i):
    #R_1,Convective = 1/(h_1,convectove *A)
    #A = inner area (based on inner diameter)
    #[W/cm K]/ [cm] = [W/K cm^2]
    h_i = NusseltNumber(Section, correlation, Twall, i)*ThermalConductivity(Twall, "Alloy-800")/HeatedEquivalentDiameter(Section)[i]
    return 1/(h_i*InnerArea(Section)[i]) #[K/W]
    

def SecondaryConvectionResistance(Section, Twall, preheater, i):
    if preheater == "yes": #from Silpsrikul thesis
        f_b = 0.1783 #fraction of cross-sectional area of shell occupied by a baffle window
        NumberTubes = 3550
        N_b = f_b*NumberTubes #number of tubes in baffle window 
        S_b = 0.1271*(100**2) #[cm^2]  
        G_e = 1512*1000/(100**2) #[g/cm^2 s] weighted average mass velocity in preheater 
        
        #First two (raised to exponent) terms are unitless 
        h_o = (ThermalConductivity(Twall, "water")*0.2/Section.OuterDiameter[i])*\
        ((Section.OuterDiameter[i]*G_e/Section.ViscosityH2O[i])**0.6)*\
        (Cp(Twall)*Section.ViscosityH2O[i]/ThermalConductivity(Twall, "water"))**0.33
        
    elif preheater == "no":
        h_o = 1 #boiling heat transfer 
    else:
        print ("Error: preheater not specified")
    
    return 1/(h_o*OuterArea(Section)[i])
    #split into boiling and non-boiling (preheater) sections 


def HeatedEquivalentDiameter(Section):
    #for channels heated only on one side 
    #D_H = A_Cross/P_Heated = (Pi/4)*(D_o-D_i)^2/(Pi*Di)
    return  [(x**2 -y**2)/y for x,y in zip(Section.OuterDiameter, Section.Diameter)] #[cm]
    
    
def NusseltNumber(Section, correlation, Twall,i):
    Re_D = ld.ReynoldsNumber(Section, HeatedEquivalentDiameter(Section))
    
    if correlation == "Dittus-Boetler":
        n= 1/3
    elif correlation == "Colburn":
        n = 0.4
    else:
        print ("Error: Correlation not recognized in NusseltNumber function")
    
    C = 0.023
    m = 4/5
    n = 0.4     
    return C*(Re_D[i]**m)*((Prandtl(Section, Twall, i))**n)


def Cp(Twall):
    A = 92.053
    B = -3.9953E-02
    C = -2.1103E-04
    D = 5.3469E-07        
    #[J/mol K]/[g/mol] = [J/g K]
    return (A + B*Twall + C*(Twall**2) +D*(Twall**3))/(nc.H2OMolarMass) #[J/g K] 


def Prandtl(Section, Twall, i):
    #Need a better reference for Cp/viscosity data 
    return Cp(Twall)*Section.ViscosityH2O[i]/ConductionResistance(Section, Twall, i) 


def InnerArea(Section):
    return [np.pi*x*y for x,y in zip(Section.Diameter,Section.Length)]


def OuterArea(Section):
    return [np.pi*x*y for x,y in zip(Section.OuterDiameter,Section.Length)]

#Can either do this as per node or for all nodes at once (depends on if output from Node 1 = needed for Node 2, i.e., U_overall for
#initial guess of bulk T
#Looks like output of one section = input of next for bulk temperatures 

def WallTemperature(Section, i, PrimaryBulkTemp, SecondaryBulkTemp):
    #i = each node of SG tube 
    end = Section.NodeNumber-1
    
    for k in range(2):
        if k ==0: #PrimaryBulkTemp will update for next time iteration for WallTemperature called 
            #assumes initial linear profile across bulk hot and bulk cold temperatures
            PrimaryWallTemperature = PrimaryBulkTemp[i]-(1/3)*(PrimaryBulkTemp[i]-SecondaryBulkTemp[i]) #perhaps call from outside function
            SecondaryWallTemperature = PrimaryWallTemperature*1
            
        WT_h = PrimaryWallTemperature
        WT_c = SecondaryWallTemperature
        
        
        
        PrimaryWallTemperature = 188+273.15
        
        if Section.Distance[i] >= Section.Distance[end]-263.5: #(2.635 meters up cold-side leg) length at which pre-heater ends (cross-flow ends) 
            preheater = "yes"
        else:
            preheater = "no"
        h_i = PrimaryConvectionResistance(Section, "Dittus-Boetler", PrimaryWallTemperature, preheater, i)
        k_w = ConductionResistance(Section, PrimaryWallTemperature, i)
        h_o = None
    
        #return (h_i)


WallTemperature(ld.SG_Zone1, 20, ld.SG_Zone1.PrimaryBulkTemperature, [187+273.15]*ld.SG_Zone1.NodeNumber)

        
        
        
        
        
        
        
