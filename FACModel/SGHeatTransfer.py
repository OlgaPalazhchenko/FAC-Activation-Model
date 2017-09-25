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


def PrimaryConvectionResistance(Section, correlation, Tfilm, i):
    #R_1,Convective = 1/(h_1,convectove *A)
    #A = inner area (based on inner diameter)
    #[W/cm K]/ [cm] = [W/K cm^2]
    h_i = NusseltNumber(Section, correlation, Tfilm, i)*ThermalConductivity(Tfilm, "water")/HeatedEquivalentDiameter(Section)[i]
    return 1/(h_i*InnerArea(Section)[i]) #[K/W]
    

def SecondaryConvectionResistance(Section, Tfilm, Twall, preheater, i):
    T_sat = 260 + 273.15 #[K]
    if preheater == "yes": #from Silpsrikul thesis
        f_b = 0.1783 #fraction of cross-sectional area of shell occupied by a baffle window
        NumberTubes = 3550
        N_b = f_b*NumberTubes #number of tubes in baffle window 
        S_b = 0.1271*(100**2) #[cm^2]  
        G_e = 1512*1000/(100**2) #[g/cm^2 s] weighted average mass velocity in preheater 
        
        #First two (raised to exponent) terms are unitless 
        #[W/cm K]/[cm] = [W/cm^2 K]
        h_o = (ThermalConductivity(Tfilm, "water")*0.2/Section.OuterDiameter[i])*\
        ((Section.OuterDiameter[i]*G_e/ld.Viscosity("water", "SHT", Tfilm))**0.6)*\
        (ld.HeatCapacity("PHT", Twall)*ld.Viscosity("water", "PHT", Tfilm)/ThermalConductivity(Tfilm, "water"))**0.33
    
    elif preheater == "no":
        x=25
        rho_v = 1000*23.753/(100**3) #[g/cm^3]
        p_crit = 22.0640 #[MPa]
       
        F = (1+(x*Prandtl(Section, T_sat, i)*((ld.Density("water", "SHT", T_sat)/rho_v)-1)))**0.35
        
        MassFlow = 239.53 #kg/s
        Shell_ID = 2000 #[cm]
        Re_D = 4*MassFlow/((ld.Viscosity("water", "SHT", T_sat)/1000)*np.pi*Shell_ID)
       
        h_l = 0.023*ld.ThermalConductivityH2O("SHT", T_sat)*((Re_D)**0.8)*(Prandtl(Section, T_sat, i)**0.4)\
        /Shell_ID #[W/cm^2 K]
        
        Q_prime =F*h_l*(Twall-T_sat) #[W/cm^2]
                
        S = (1+0.055*(F**0.1)*(Re_D)**0.16)**(-1) 
        A_p = (55*((4.70/p_crit)**0.12)*((-np.log10(4.70/p_crit))**(-0.55))*(nc.H2MolarMass)**(-0.5))/(100**2)
        C = ((A_p*S)/(F*h_l)**2)*(Q_prime)**(4/3)
        coeff = [1, -C, -1]
        cubic_solution = np.roots(coeff)
        q = cubic_solution[0]
        
        h_o = F*(q**(3/2))*h_l #[W/cm^2 K
        
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

def TemperatureProfile(Section):
    PrimaryBulkTemp = []
    SecondaryBulkTemp = []
    for i in range(Section.NodeNumber):
        if i == 0:
            PrimaryBulkTemp[i] = Section.PrimaryBulkTemperature[i]
            SecondaryBulkTemp[i] = 260 + 273 #[K]
        
            

def WallTemperature(Section, i, PrimaryBulkTemp, SecondaryBulkTemp):
    #i = each node of SG tube 
    end = Section.NodeNumber-1
    PrimaryWallTemperature = PrimaryBulkTemp-(1/3)*(PrimaryBulkTemp-SecondaryBulkTemp) #perhaps call from outside function
    SecondaryWallTemperature =PrimaryBulkTemp*1
    
    for k in range(100):        
        WT_h = PrimaryWallTemperature
        WT_c = SecondaryWallTemperature
        
        PrimaryT_film = (PrimaryBulkTemp + PrimaryWallTemperature)/2
        SecondaryT_film = (SecondaryBulkTemp + SecondaryWallTemperature)/2

        if Section.Distance[i] >= Section.Distance[end]-263.5: #(2.635 meters up cold-side leg) length at which pre-heater ends (cross-flow ends) 
            preheater = "yes"
        else:
            preheater = "no"
       
        
        U_h = 1/((ConductionResistance(Section, PrimaryWallTemperature, i) + \
               SecondaryConvectionResistance(Section, SecondaryT_film, SecondaryWallTemperature, preheater, i))*\
               InnerArea(Section)[i])
        
        U_c = 1/((ConductionResistance(Section, PrimaryWallTemperature, i) +\
               PrimaryConvectionResistance(Section, "Dittus-Boetler", PrimaryT_film, i))*OuterArea(Section)[i])
        
    
        
        #U_total = 1/(R_total*(InnerArea(Section)[i]+OuterArea(Section)[i])/2)
        
        h_i = 1/(PrimaryConvectionResistance(Section, "Dittus-Boetler", PrimaryT_film, i)*InnerArea(Section)[i])
        h_o = 1/(SecondaryConvectionResistance(Section, SecondaryT_film, SecondaryWallTemperature, preheater, i)*OuterArea(Section)[i])
        
        PrimaryWallTemperature = (PrimaryBulkTemp*h_i + SecondaryBulkTemp*(U_h))/(h_i + U_h)
        SecondaryWallTemperature = ((SecondaryBulkTemp*h_o) + PrimaryBulkTemp*U_c)/(h_o + U_c) 
        
        RE1 = (PrimaryWallTemperature-WT_h)
        RE2 = (SecondaryWallTemperature - WT_c)
    
        if abs(RE1) < 0.01 and abs(RE2) <0.01:
            print (PrimaryWallTemperature-273.15, SecondaryWallTemperature-273.15, k, h_o)
            return SecondaryWallTemperature-273.15, SecondaryBulkTemp-273.15
        
        
        
WallTemperature(ld.SG_Zone1, 1, 305+273.15, 260+273)

        

        
        
        
        
