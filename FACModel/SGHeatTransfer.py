import LepreauData as ld
import numpy as np


def ConductionResistance(Section, Twall, i):
    #R_conduction
    ThermalConductivity = 11.450 + 0.0161*Twall #M.K.E thesis for Alloy-800 tubing [W/m K]   
    
    Thickness = 0.125 #cm
    Do = Section.Diameter[i] + 2*Thickness 
    Rcyl_numerator = np.log(Do/Section.Diameter) 
    
    Rcyl_denominator = 2*np.pi*ThermalConductivity*Section.Length[i] 
    
    return Rcyl_numerator/Rcyl_denominator #[K/W]


def PrimaryConvectionResistance(Section, correlation, Twall, i):
    return NusseltNumber(Section, correlation, Twall, i)*ConductionResistance(Section, Twall, i)*HeatedEquivalentDiameter(Section)[i]


def HeatedEquivalentDiameter(Section):
    #for channels heated only on one side 
    #D_H = A_Cross/P_Heated = (Pi/4)*(D_o-D_i)^2/(Pi*Di)
    Thickness = 0.125 #cm
    Do = [i + 2*Thickness for i in Section.Diameter] 
    return  [(x**2 -y**2)/y for x,y,z in zip(Do, Section.Diameter)] #[cm]
    
    
def NusseltNumber(Section, correlation, Twall,i):
    Re_D = ld.ReynoldsNumber(Section, HeatedEquivalentDiameter(Section))
    
    if correlation == "Dittus-Boetler":
        n= 1/3
    elif correlation == "Colburn":
        n = 0.4
    else:
        print ("Correlation not recognized")
    
    C = 0.023
    m = 4/5
    n = 0.4     
    return C*(Re_D**m)*((Prandtl(Section, Twall, i))**n)


def Prandtl(Section, Twall, i):
    #Need a better reference for Cp/viscosity data 
    A = 92.053
    B = -3.9953E-02
    C = -2.1103E-04
    D = 5.3469E-07
    
    Cp = A + B*Twall + C*(Twall**2) +D*(Twall**3) 
    
    return Cp*Section.ViscosityH2O[i]/ConductionResistance(Section, Twall) 

