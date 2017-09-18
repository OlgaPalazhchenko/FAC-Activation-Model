import LepreauData as ld
import numpy as np

def ThermalConductivity(Twall, material):
    if material == "Alloy-800":
        return (11.450 + 0.0161*Twall)/100 #M.K.E thesis for Alloy-800 tubing [W/cm K]
    elif material == "water":
        return None
    else: print ("Error: material not specified")


def ConductionResistance(Section, Twall, i):
    #R_conduction = L/kA (L = thickness)
    #A = area with shape factor correcrion factor for cylindrical pipe
    #R_conduction = ln(D_o/D_i)/(2pikl) [K/W]
    #l = length, k =thermal conductivity coefficient [W/cm K]
    
    k_w = ThermalConductivity(Twall)  #[W/cm K]  
    Rcyl_numerator = np.log(Section.OuterDiameter[i]/Section.Diameter[i]) 
    Rcyl_denominator = 2*np.pi*k_w*Section.Length[i] 
    return Rcyl_numerator/Rcyl_denominator #[K/W]


def PrimaryConvectionResistance(Section, correlation, Twall, preheater, i):
    #R_1,Convective = 1/(h_1,convectove *A)
    #A = inner area (based on inner diameter)
    #[W/cm K]/ [cm] = [W/K cm^2]
    if preheater == "yes": #from Silpsrikul thesis
        f_b = 0.1282 #fractional 
        
    elif preheater == "no":
        h_i = NusseltNumber(Section, correlation, Twall, i)*ThermalConductivity(Twall)/HeatedEquivalentDiameter(Section)[i]
        return 1/(h_i*InnerArea(Section)[i]) #[K/W]
    else: print ("Error: pre-heater location not specified")


def SecondaryConvectionResistance(Section, Twall, i):
    None
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
        print ("Error: Correlation not recognized")
    
    C = 0.023
    m = 4/5
    n = 0.4     
    return C*(Re_D[i]**m)*((Prandtl(Section, Twall, i))**n)


def Prandtl(Section, Twall, i):
    #Need a better reference for Cp/viscosity data 
    A = 92.053
    B = -3.9953E-02
    C = -2.1103E-04
    D = 5.3469E-07    
    Cp = A + B*Twall + C*(Twall**2) +D*(Twall**3) 
    return Cp*Section.ViscosityH2O[i]/ConductionResistance(Section, Twall, i) 


def InnerArea(Section):
    return [np.pi*x*y for x,y in zip(Section.Diameter,Section.Length)]


def OuterArea(Section):
    return [np.pi*x*y for x,y in zip(Section.OuterDiameter,Section.Length)]

#Can either do this as per node or for all nodes at once (depends on if output from Node 1 = needed for Node 2, i.e., U_overall for
#initial guess of bulk T
#Looks like output of one section = input of next for bulk temperatures 

def WallTemperature(Section, i, PrimaryBulkTemp, SecondaryBulkTemp):
    #i = each node of SG tube 
    for k in range(2):
        if k ==0: #PrimaryBulkTemp will update for next time iteration for WallTemperature called 
            #assumes initial linear profile across bulk hot and bulk cold temperatures
            PrimaryWallTemperature = (PrimaryBulkTemp[i]-(1/3)*(PrimaryBulkTemp[i]-SecondaryBulkTemp[i]))-273.15 #perhaps call from outside function
        
        WT = PrimaryWallTemperature
        PrimaryWallTemperature = 1
        print (WT, PrimaryWallTemperature, Section.Distance)
        
        h_i = PrimaryConvectionResistance(Section, "Dittus-Boetler", PrimaryWallTemperature, i)
        
        if Section.Length[i] <= 263.5: #(2.635 meters up cold-side leg) length at which pre-heater ends (cross-flow ends) 
            h_i = 1 #'preheater h_i calculated with baffles/cross flow - data from Silpsrikul thesis 
        
        k_w = ConductionResistance(Section, PrimaryWallTemperature, i)



WallTemperature(ld.SG_Zone1, 0, ld.SG_Zone2.PrimaryBulkTemperature, [260+273.15]*ld.SG_Zone1.NodeNumber)

        
        
        
        
        
        
        
