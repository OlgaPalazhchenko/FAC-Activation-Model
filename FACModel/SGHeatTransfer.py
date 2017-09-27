import LepreauData as ld
import numpy as np
import NumericConstants as nc
import scipy.optimize
import math

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
        x=2
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
   
    
def NusseltNumber(Section, correlation, Temperature,i):
    Re_D = ld.ReynoldsNumber(Section, HeatedEquivalentDiameter(Section))
    
    if correlation == "Dittus-Boetler":
        n= 1/3
    elif correlation == "Colburn":
        n = 0.4
    else:
        print ("Error: Correlation not recognized in NusseltNumber function")
    
    C = 0.023
    m = 4/5
    
    return C*(Re_D[i]**m)*((Prandtl(Section, Temperature, i))**n)


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
        
        h_i = 1/(PrimaryConvectionResistance(Section, "Dittus-Boetler", PrimaryT_film, i)*InnerArea(Section)[i])
        h_o = 1/(SecondaryConvectionResistance(Section, SecondaryT_film, SecondaryWallTemperature, preheater, i)*OuterArea(Section)[i])
        
        U_h1 = ConductionResistance(Section, PrimaryWallTemperature, i)*InnerArea(Section)[i]
        U_h2 = SecondaryConvectionResistance(Section, SecondaryT_film, SecondaryWallTemperature, preheater, i)*InnerArea(Section)[i]
        U_h = 1/(U_h1+ U_h2) #[W/ cm^2 K]
        
        U_c1 = ConductionResistance(Section, PrimaryWallTemperature, i)*OuterArea(Section)[i]
        U_c2 = PrimaryConvectionResistance(Section, "Dittus-Boetler", PrimaryT_film, i)*OuterArea(Section)[i]
        U_c = 1/(U_c1+U_c2) #[W/ cm^2 K]

        inverseU_total = (PrimaryConvectionResistance(Section, "Dittus-Boetler", PrimaryT_film, i) + \
        ConductionResistance(Section, PrimaryWallTemperature, i) + \
        SecondaryConvectionResistance(Section, SecondaryT_film, SecondaryWallTemperature, preheater, i))*InnerArea(Section)[i]\
        + (100*100*(6E-5))
        
        U_total = 1/inverseU_total #[W/ cm^2 K]
        
        #[W/K cm^2] * [K] * [cm^2] = [W] = [J/s]
        #Q = U_h*InnerArea(Section)[i]*(PrimaryWallTemperature-SecondaryBulkTemp)    
        Q = h_i*InnerArea(Section)[i]*(PrimaryBulkTemp-PrimaryWallTemperature)
        #[W] = [W/cm^2 K][cm^2]
        #print (Q)
        T_logmean = Q/(U_total*OuterArea(Section)[i])
        

            
#         U_totes = 1/(((OuterArea(Section)[i]/InnerArea(Section)[i])/h_i) + \
#         (0.5*Section.OuterDiameter[i]*np.log(Section.OuterDiameter[i]/Section.Diameter[i])/ThermalConductivity(PrimaryWallTemperature, "Alloy-800"))+\
#         (1/h_o))
        
        #R = 1/h*A, 1/R = h*A --> h =1/R*A   
#         U_z1 =  Section.Diameter[i]*np.log(Section.OuterDiameter[i]/Section.Diameter[i])/(2*ThermalConductivity(PrimaryWallTemperature, "Alloy-800"))
#         U_z2 = Section.Diameter[i]/(Section.OuterDiameter[i]*h_o)
#         U_z = 1/(U_z1+U_z2)
#           
#         U_y1 =  Section.OuterDiameter[i]*np.log(Section.OuterDiameter[i]/Section.Diameter[i])/(2*ThermalConductivity(PrimaryWallTemperature, "Alloy-800"))
#         U_y2 = Section.OuterDiameter[i]/(Section.Diameter[i]*h_i)
#         U_y = 1/(U_y1+U_y2)
        
        #U_total = 1/(R_total*(InnerArea(Section)[i]+OuterArea(Section)[i])/2)
        
        PrimaryWallTemperature = (PrimaryBulkTemp*h_i + SecondaryBulkTemp*U_h)/(h_i + U_h)
        SecondaryWallTemperature = (SecondaryBulkTemp*h_o + PrimaryBulkTemp*U_c)/(h_o + U_c) 
        
        RE1 = (PrimaryWallTemperature-WT_h)
        RE2 = (SecondaryWallTemperature - WT_c)
    
        if abs(RE1) < 0.01 and abs(RE2) <0.01:
            T_bh =scipy.optimize.fsolve(lambda x: (x-PrimaryBulkTemp)/\
            (math.log((x-SecondaryBulkTemp)/(PrimaryBulkTemp-SecondaryBulkTemp)))-T_logmean, (PrimaryBulkTemp-5))[0] # [0] is guess value

            #print (100*100*h_i, 100*100*h_o)
            
            return PrimaryWallTemperature, SecondaryWallTemperature, T_bh
        
        
WallTemperature(ld.SG_Zone1, 1, 583.15, 533.15)



    
def TemperatureProfile(Section):
    PrimaryWall = []
    PrimaryBulk = []
    SecondaryBulk = []
    SecondaryWall = []
    T_bc = 533.15
    for i in range(18):
        
        if i == 0:
            T_bh = 583.15 #[K]
            
             
        PrimaryBulk.append(T_bh)    
        #SecondaryBulk.append(T_bc)
        T_wh, T_wc, T_bh = WallTemperature(Section, i, T_bh, T_bc)
        print (T_bh)
        #m = V*rho = u*A*rho = [cm/s]*[cm^2]*[g/cm^3] = [g/s] 
         
#         T_bh = ((m_h*ld.HeatCapacity("PHT", PrimaryBulkTemp)-0.5*U_total*OuterArea(Section)[i])*PrimaryBulkTemp + 0.5*U_total*OuterArea(Section)[i]*\
#         (533.15 + 533.15))/(m_h*ld.HeatCapacity("PHT", PrimaryBulkTemp)-0.5*U_total*OuterArea(Section)[i])
#         
        
        PrimaryWall.append(T_wh)
        SecondaryWall.append(T_wc)
        
        
    print ([j-273 for j in PrimaryBulk])
    print ([j-273 for j in PrimaryWall])
TemperatureProfile(ld.SteamGenerator)
# rofl = []
# copter = []
# for i in range(3):
#     if i == 0:
#         a=4
#     else:
#         a=c
#     b=a*3
#     c=3
#     rofl.append(a)
#     copter.append(b)
# 
# print (rofl,copter)
    
    