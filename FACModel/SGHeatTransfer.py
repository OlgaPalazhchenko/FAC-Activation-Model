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
    h_i = NusseltNumber(Section, correlation, Tfilm, i)*ThermalConductivity(Tfilm, "water")/Section.Diameter[i]
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
        #V*rho =[cm/s]([cm^2]*[g/cm^3] = [g/s]
        
        F = (1+(x*Prandtl(Section, T_sat, i)*((ld.Density("water", "SHT", T_sat)/rho_v)-1)))**0.35
       
        MassFlow = 239.53*1000 #[g/s]
        Shell_ID = 2.28*100 #[cm]
        MassFlux = MassFlow/((np.pi/4)*(Shell_ID**2)) #[g/cm^2 s]
        
        #[kg/s]*[kg/cm s]
        Re_D = Section.OuterDiameter[i]*MassFlux/ld.Viscosity("water", "SHT", T_sat)#4*MassFlow/((ld.Viscosity("water", "SHT", T_sat)/1000)*np.pi*Shell_ID)
        h_l = 0.023*ld.ThermalConductivityH2O("SHT", T_sat)*((Re_D)**0.8)*(Prandtl(Section, T_sat, i)**0.4)/Section.OuterDiameter[i] #[W/cm^2 K]
        
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
    Re_D = ld.ReynoldsNumber(Section, Section.Diameter)
    
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


def WallTemperature(Section, i, T_PrimaryBulkIn, T_SecondaryBulkIn):
    #i = each node of SG tube 
    end = Section.NodeNumber-1

    T_PrimaryWall = T_PrimaryBulkIn-(1/3)*(T_PrimaryBulkIn-T_SecondaryBulkIn) #perhaps call from outside function
    T_SecondaryWall =T_PrimaryBulkIn*1
    
    for k in range(50):        
        WT_h = T_PrimaryWall
        WT_c = T_SecondaryWall
        
        PrimaryT_film = (T_PrimaryBulkIn + T_PrimaryWall)/2
        SecondaryT_film = (T_SecondaryBulkIn + T_SecondaryWall)/2

        if Section.Distance[i] >= Section.Distance[end]-263.5: #(2.635 meters up cold-side leg) length at which pre-heater ends (cross-flow ends) 
            preheater = "yes"
        else:
            preheater = "no"
        
        h_i = 1/(PrimaryConvectionResistance(Section, "Dittus-Boetler", PrimaryT_film, i)*InnerArea(Section)[i])
        h_o = 1/(SecondaryConvectionResistance(Section, SecondaryT_film, T_SecondaryWall, preheater, i)*OuterArea(Section)[i])
        
        U_h1 = ConductionResistance(Section, T_PrimaryWall, i)*InnerArea(Section)[i]
        U_h2 = SecondaryConvectionResistance(Section, SecondaryT_film, T_SecondaryWall, preheater, i)*InnerArea(Section)[i]
        U_h = 1/(U_h1+ U_h2) #[W/ cm^2 K]
        
        U_c1 = ConductionResistance(Section, T_PrimaryWall, i)*OuterArea(Section)[i]
        U_c2 = PrimaryConvectionResistance(Section, "Dittus-Boetler", PrimaryT_film, i)*OuterArea(Section)[i]
        U_c = 1/(U_c1+U_c2) #[W/ cm^2 K]
        
        
        
        
#         U_totes = 1/(((OuterArea(Section)[i]/InnerArea(Section)[i])/h_i) + \
#         (0.5*Section.OuterDiameter[i]*np.log(Section.OuterDiameter[i]/Section.Diameter[i])/ThermalConductivity(T_PrimaryWall, "Alloy-800"))+\
#         (1/h_o))
        
        #R = 1/h*A, 1/R = h*A --> h =1/R*A   
#         U_z1 =  Section.Diameter[i]*np.log(Section.OuterDiameter[i]/Section.Diameter[i])/(2*ThermalConductivity(T_PrimaryWall, "Alloy-800"))
#         U_z2 = Section.Diameter[i]/(Section.OuterDiameter[i]*h_o)
#         U_z = 1/(U_z1+U_z2)
#           
#         U_y1 =  Section.OuterDiameter[i]*np.log(Section.OuterDiameter[i]/Section.Diameter[i])/(2*ThermalConductivity(T_PrimaryWall, "Alloy-800"))
#         U_y2 = Section.OuterDiameter[i]/(Section.Diameter[i]*h_i)
#         U_y = 1/(U_y1+U_y2)
        
        #U_total = 1/(R_total*(InnerArea(Section)[i]+OuterArea(Section)[i])/2)
        
        T_PrimaryWall = (T_PrimaryBulkIn*h_i + T_SecondaryBulkIn*U_h)/(h_i + U_h)
        T_SecondaryWall = (T_SecondaryBulkIn*h_o + T_PrimaryBulkIn*U_c)/(h_o + U_c) 
        
        RE1 = (T_PrimaryWall-WT_h)
        RE2 = (T_SecondaryWall - WT_c)
        
        if abs(RE1) <= 0.01 and abs(RE2) <= 0.01:
        
            foulingresistance = (100*100*(1.15E-5))
        
            inverseU_total = (PrimaryConvectionResistance(Section, "Dittus-Boetler", PrimaryT_film, i) + ConductionResistance(Section, T_PrimaryWall, i) + \
                    SecondaryConvectionResistance(Section, SecondaryT_film, T_SecondaryWall, preheater, i))*InnerArea(Section)[i]\
                    + foulingresistance
                    
            U_total = 1/inverseU_total #[W/ cm^2 K]
        
            #[W/K cm^2] * [K] * [cm^2] = [W] = [J/s]
            #Q = U_h*InnerArea(Section)[i]*(T_PrimaryWall-T_SecondaryBulkIn)    
            Q = h_i*InnerArea(Section)[i]*(T_PrimaryBulkIn-T_PrimaryWall)
            #[W] = [W/cm^2 K][cm^2]
        
            T_logmean = Q/(U_total*OuterArea(Section)[i])
            
            if preheater == "no":
                GuessValue = T_PrimaryBulkIn-5
                
                T_PrimaryBulkOut =scipy.optimize.fsolve(lambda x: ((x-T_PrimaryBulkIn)/\
            math.log((x-T_SecondaryBulkIn)/(T_PrimaryBulkIn-T_SecondaryBulkIn)))-T_logmean, (GuessValue))[0] # [0] is guess value
                
                T_SecondaryBulkOut = 533.15
                tt = T_SecondaryBulkOut
            elif preheater == "yes":
                
                
                T_SecondaryBulkOut = T_SecondaryBulkIn
                T_SecondaryBulkIn = T_SecondaryBulkOut-4 #Guessed temperature 
                
                 
                T_PrimaryBulkOut = scipy.optimize.fsolve(lambda x: (((x-T_SecondaryBulkIn)-(T_PrimaryBulkIn-T_SecondaryBulkOut))/\
                                 math.log((x-T_SecondaryBulkIn)/(T_PrimaryBulkIn-T_SecondaryBulkOut)))-T_logmean, (T_PrimaryBulkIn-2))[0]
                 
            #SOLVING for BulkIN 
                T_SecondaryBulkIn = scipy.optimize.fsolve(lambda x: (((T_PrimaryBulkOut-x)-(T_PrimaryBulkIn-T_SecondaryBulkOut))/\
                                math.log((T_PrimaryBulkOut-x)/(T_PrimaryBulkIn-T_SecondaryBulkOut)))-T_logmean, (T_SecondaryBulkOut-1))[0]
                    
                tt = T_SecondaryBulkIn
                    
            else:
                print ("Error: must specify if in preheater")
            
            
            #print (T_PrimaryBulkIn-273.15, T_SecondaryBulkIn-273.15, T_SecondaryBulkOut-273.15)
            return T_PrimaryWall, T_SecondaryWall, T_PrimaryBulkOut, tt
        
        
#WallTemperature(ld.SG_Zone1, 1, 583.15, 533.15)

# def Preheater(Section):
#     end = Section.NodeNumber -1
#     PreheaterNodes = []
#     for i in range(Section.NodeNumber):
#         if Section.Distance[i] >= Section.Distance[end]-263.5: #(2.635 meters up cold-side leg) length at which pre-heater ends (cross-flow ends)
#             PreheaterNodes.append(Section.Length[i])
#     PreheaterLengths = list(reversed(PreheaterNodes))
#     
#     for i in PreheaterLengths:
#         
#     
# Preheater(ld.SteamGenerator)

#Need to model preheater separately
    
def TemperatureProfile(Section):
    PrimaryWall = []
    PrimaryBulk = []
    SecondaryBulk = []
    SecondaryWall = []
    
    for i in range(15):    
        if i == 0:
            T_bh = 583.15 #[K]
            T_bc = 533.15
    
        PrimaryBulk.append(T_bh)    
        SecondaryBulk.append(T_bc)
        
        
        T_wh, T_wc, T_bh, T_bc = WallTemperature(Section, i, T_bh, T_bc)
#         
        PrimaryWall.append(T_wh)
        SecondaryWall.append(T_wc)
        
        
    print ([j-273.15 for j in PrimaryBulk])
    print ([j-273.15 for j in PrimaryWall])
    print ()
    print ([j-273.15 for j in SecondaryBulk])
    print ([j-273.15 for j in SecondaryWall])
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
    
    