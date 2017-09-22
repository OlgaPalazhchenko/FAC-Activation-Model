import LepreauData as ld
import numpy as np
import NumericConstants as nc


def ThermalConductivity(Twall, material):
    if material == "Alloy-800" or material == "Alloy800" or material == "A800" or material == "Alloy 800":
        return (11.450 + 0.0161*Twall)/100 #M.K.E thesis for Alloy-800 tubing [W/cm K]
    elif material == "water":
        return ld.ThermalConductivityH2O("PHT", Twall)
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


def PrimaryConvectionResistance(Section, correlation, Twall, i):
    #R_1,Convective = 1/(h_1,convectove *A)
    #A = inner area (based on inner diameter)
    #[W/cm K]/ [cm] = [W/K cm^2]
    h_i = NusseltNumber(Section, correlation, Twall, i)*ThermalConductivity(Twall, "water")/HeatedEquivalentDiameter(Section)[i]
    return 1/(h_i*InnerArea(Section)[i]) #[K/W]
    

def SecondaryConvectionResistance(Section, Twall, preheater, i):
    T_sat = 260 + 273.15 #[K]
    if preheater == "yes": #from Silpsrikul thesis
        f_b = 0.1783 #fraction of cross-sectional area of shell occupied by a baffle window
        NumberTubes = 3550
        N_b = f_b*NumberTubes #number of tubes in baffle window 
        S_b = 0.1271*(100**2) #[cm^2]  
        G_e = 1512*1000/(100**2) #[g/cm^2 s] weighted average mass velocity in preheater 
        
        #First two (raised to exponent) terms are unitless 
        #[W/cm K]/[cm] = [W/cm^2 K]
        h_o = (ThermalConductivity(Twall, "water")*0.2/Section.OuterDiameter[i])*\
        ((Section.OuterDiameter[i]*G_e/ld.Viscosity("water", "SHT", Twall))**0.6)*\
        (ld.HeatCapacity("PHT", Twall)*ld.Viscosity("water", "PHT", Twall)/ThermalConductivity(Twall, "water"))**0.33
       
    elif preheater == "no":
        h_o = 2.54*(Twall-T_sat)*np.exp(4.96/1.551) #[W/cm^2 K] #boiling heat transfer 
    else:
        print ("Error: preheater not specified")
    
    return 1/(h_o*OuterArea(Section)[i]) #K/W
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


def Prandtl(Section, Twall, i):
    #Need a better reference for Cp/viscosity data 
    return ld.HeatCapacity("PHT", Twall)*Section.ViscosityH2O[i]/ConductionResistance(Section, Twall, i) 


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
    PrimaryWallTemperature = PrimaryBulkTemp[i]-(1/3)*(PrimaryBulkTemp[i]-SecondaryBulkTemp[i]) #perhaps call from outside function
    SecondaryWallTemperature = PrimaryWallTemperature*1
    print (PrimaryWallTemperature-273.15)
    
    for k in range(100):        
        WT_h = PrimaryWallTemperature
        WT_c = SecondaryWallTemperature
        
        PrimaryT_film = (PrimaryBulkTemp[i] + PrimaryWallTemperature)/2
        SecondaryT_film = (SecondaryBulkTemp[i] + SecondaryWallTemperature)/2
        
        if Section.Distance[i] >= Section.Distance[end]-263.5: #(2.635 meters up cold-side leg) length at which pre-heater ends (cross-flow ends) 
            preheater = "yes"
        else:
            preheater = "no"
       
        R_i = PrimaryConvectionResistance(Section, "Dittus-Boetler", PrimaryT_film, i)
        R_cond = ConductionResistance(Section, PrimaryWallTemperature, i)
        R_o = SecondaryConvectionResistance(Section, SecondaryT_film, preheater, i)
        
        U_h = 1/((R_o + R_cond)*InnerArea(Section)[i])
        U_c = 1/((R_i + R_cond)*OuterArea(Section)[i])
        
        h_i = 1/(R_i*InnerArea(Section)[i])
        h_o = 1/(R_o*OuterArea(Section)[i])
        

        PrimaryWallTemperature = (PrimaryBulkTemp[i]*h_i + SecondaryBulkTemp[i]*(U_h))/(h_i + U_h)
        SecondaryWallTemperature = (SecondaryBulkTemp[i]*h_o + PrimaryBulkTemp[i]*(U_c))/(h_o + U_c)
        
        RE = (PrimaryWallTemperature-WT_h)
        print (PrimaryWallTemperature-273.15, k, RE)
        
        if abs(RE) < 0.1:
            
            return PrimaryWallTemperature-273.15
        
        
        
WallTemperature(ld.SG_Zone1, 21, ld.SG_Zone1.PrimaryBulkTemperature, [260+273.15]*ld.SG_Zone1.NodeNumber)

        

        
        
        
        
