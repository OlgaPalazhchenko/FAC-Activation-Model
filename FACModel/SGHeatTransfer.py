import LepreauData as ld
import numpy as np
import NumericConstants as nc
import scipy.optimize
import math

h_i = ld.SGParameters()
h_o = ld.SGParameters()
U_h = ld.SGParameters()
U_c = ld.SGParameters()
U_total = ld.SGParameters()
R_F = ld.SGParameters()
MassFlow_c = ld.SGParameters()
ShellDiameter = ld.SGParameters()
MassFlux = ld.SGParameters() 
EnthalpySaturatedSteam = ld.SGParameters()

for i in [h_i, h_o]:
    i.unit = "K/W"
    
for i in [U_h, U_c, U_total]:
    i.unit = "W/cm^2 K"

R_F.unit = "cm^2 K/W"
MassFlow_c.magnitude = 240*1000 
MassFlow_c.unit = "g/s"
ShellDiameter.magnitude = 2.28*100 #[cm]
ShellDiameter.unit = "cm"
MassFlux.magnitude = MassFlow_c.magnitude/((np.pi/4)*(ShellDiameter.magnitude**2)) #[g/cm^2 s]
MassFlux.unit = "g/cm^2 s"
EnthalpySaturatedSteam.magnitude = 2800.4
EnthalpySaturatedSteam.unit = "J/g"


def ThermalConductivity(Twall, material):
    if material == "Alloy-800" or material == "Alloy800" or material == "A800" or material == "Alloy 800":
        return (11.450 + 0.0161*Twall)/100 #M.K.E thesis for Alloy-800 tubing [W/cm K]
    elif material == "water":
        return ld.ThermalConductivityH2O("PHT", Twall)
    elif material == "magnetite":
        return 1.4/100 #[W/cm K]
    else: print ("Error: material not specified")


def FoulingResistance(Section):
    #thickness/thermal conductivity [cm]/[W/cm K] = [cm^2 K/W]
    #[g/cm^2]/[g/cm^3] = [cm]
    OxideThickness = [i/nc.Fe3O4Density for i in  Section.OuterOxThickness]
    return [i/ThermalConductivity(None, "magnetite") for i in OxideThickness]


def ConductionResistance(Section, Twall, i):
    #R_conduction = L/kA (L = thickness)
    #A = area with shape factor correcrion factor for cylindrical pipe
    #R_conduction = ln(D_o/D_i)/(2pikl) [K/W]
    #l = length, k =thermal conductivity coefficient [W/cm K]
    
    k_w = ThermalConductivity(Twall, "Alloy-800")  #[W/cm K]  
    
    Rcyl_numerator = np.log(Section.OuterDiameter[i]/Section.Diameter[i]) 
    Rcyl_denominator = 2*np.pi*k_w*Section.Length.magnitude[i] 
    
    R_cond= Rcyl_numerator/Rcyl_denominator #[K/W]
    return R_cond

def PrimaryConvectionResistance(Section, correlation, Tfilm, i):
    #R_1,Convective = 1/(h_1,convectove *A)
    #A = inner area (based on inner diameter)
    #[W/cm K]/ [cm] = [W/K cm^2]
    h_i = NusseltNumber(Section, correlation, Tfilm, i)*ThermalConductivity(Tfilm, "water")/Section.Diameter[i]
    return 1/(h_i*InnerArea(Section)[i]) #[K/W]


def SecondaryConvectionResistance(Section, Tfilm, Twall, i):
    T_sat = 261 + 273.15 #[K]
    
    if Section.Length.label[i] == "preheater": #from Silpsrikul thesis
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
    
    else:
        x=2
        rho_v = 1000*23.753/(100**3) #[g/cm^3]
        p_crit = 22.0640 #[MPa]
        #V*rho =[cm/s]([cm^2]*[g/cm^3] = [g/s]
        
        F = (1+(x*Prandtl(Section, T_sat, i)*((ld.Density("water", "SHT", T_sat)/rho_v)-1)))**0.35
       
        Re_D = Section.OuterDiameter[i]*MassFlux.magnitude/ld.Viscosity("water", "SHT", T_sat)#4*MassFlow/((ld.Viscosity("water", "SHT", T_sat)/1000)*np.pi*Shell_ID)
        
        h_l = 0.023*ld.ThermalConductivityH2O("SHT", T_sat)*((Re_D)**0.8)*(Prandtl(Section, T_sat, i)**0.4)/Section.OuterDiameter[i] #[W/cm^2 K]
        
        Q_prime =F*h_l*(Twall-T_sat) #[W/cm^2]
        S = (1+0.055*(F**0.1)*(Re_D)**0.16)**(-1) 
        A_p = (55*((4.70/p_crit)**0.12)*((-np.log10(4.70/p_crit))**(-0.55))*(nc.H2MolarMass)**(-0.5))/(100**2)
        C = ((A_p*S)/(F*h_l)**2)*(Q_prime)**(4/3)
        coeff = [1, -C, -1]
        cubic_solution = np.roots(coeff)
        q = cubic_solution[0]
        
        h_o = F*(q**(3/2))*h_l #[W/cm^2 K
    
    return 1/(h_o*OuterArea(Section)[i]) #K/W
    #split into boiling and non-boiling (preheater) sections 


def HeatedEquivalentDiameter(Section):
    #for channels heated only on one side 
    #D_H = A_Cross/P_Heated = (Pi/4)*(D_o-D_i)^2/(Pi*Di)
    return  [(x**2 -y**2)/y for x,y in zip(Section.OuterDiameter, Section.Diameter)] #[cm]
   
MassFlow_h = ld.SGParameters()

def NusseltNumber(Section, correlation, Temperature,i):
    Re_D = ld.ReynoldsNumber(Section, Section.Diameter)
    
#     MassFlow_h.magnitude = 1791*1000/3500 #[g/s]
#     ID = Section.Diameter[i]#[cm]
#     MassFlux.magnitude = MassFlow_h.magnitude/((np.pi/4)*(ID**2)) #[g/cm^2 s]
#      
#      
#     Re_ =[(MassFlow_h.magnitude/ld.Viscosity("water", "PHT", Temperature))*i for i in Section.Diameter] 
#     print (Re_D[i], Re_[i])
    if correlation == "Dittus-Boetler":
        n= 1/3
    elif correlation == "Colburn":
        n = 0.4
    else:
        print ("Error: Correlation not recognized in NusseltNumber function")

    C = 0.023
    m = 4/5
    
    Nu = C*(Re_D[i]**m)*((Prandtl(Section, Temperature, i))**n)
    
    return Nu

def Prandtl(Section, Temperature, i):
    #Need a better reference for Cp/viscosity data 
    Pr =  ld.HeatCapacity("PHT", Temperature)*Section.ViscosityH2O[i]/ld.ThermalConductivityH2O("PHT", Temperature)
    #print (Pr)
    return Pr

def InnerArea(Section):
    return [np.pi*x*y for x,y in zip(Section.Diameter,Section.Length.magnitude)]


def OuterArea(Section):
    return [np.pi*x*y for x,y in zip(Section.OuterDiameter,Section.Length.magnitude)]


def PreheaterArea(Section):
        return [i*1 for i in [np.pi*x*y for x,y in zip(Section.OuterDiameter, [i*1 for i in Section.Length.magnitude])]]#[1*1 for i in [np.pi*x*y for x,y in zip(Section.OuterDiameter,Section.Length)]]
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

        h_i.magnitude = 1/(PrimaryConvectionResistance(Section, "Dittus-Boetler", PrimaryT_film, i)*InnerArea(Section)[i])
        h_o.magnitude = 1/(SecondaryConvectionResistance(Section, SecondaryT_film, T_SecondaryWall, i)*OuterArea(Section)[i])
        #R = 1/hA --> h=A*1/R
        
        U_h1 = ConductionResistance(Section, T_PrimaryWall, i)*InnerArea(Section)[i]
        U_h2 = SecondaryConvectionResistance(Section, SecondaryT_film, T_SecondaryWall, i)*InnerArea(Section)[i]
        U_h.magnitude = 1/(U_h1+ U_h2) #[W/ cm^2 K]
        
        U_c1 = ConductionResistance(Section, T_PrimaryWall, i)*OuterArea(Section)[i]
        U_c2 = PrimaryConvectionResistance(Section, "Dittus-Boetler", PrimaryT_film, i)*OuterArea(Section)[i]
        U_c.magnitude = 1/(U_c1+U_c2) #[W/cm^2 K]
        
#         if i==0: print (100*100*h_i.magnitude, "h_ih", 100*100*h_o.magnitude, "h_oc", \
#                (100*100)/(InnerArea(Section)[i]*ConductionResistance(Section, T_PrimaryWall, i)), "k_cond", "unit for all = [m^2 K/W]")


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
        
        T_PrimaryWall = (T_PrimaryBulkIn*h_i.magnitude + T_SecondaryBulkIn*U_h.magnitude)/(h_i.magnitude + U_h.magnitude)
        
        T_SecondaryWall = (T_SecondaryBulkIn*h_o.magnitude + T_PrimaryBulkIn*U_c.magnitude)/(h_o.magnitude + U_c.magnitude) 
        
        RE1 = (T_PrimaryWall-WT_h)
        RE2 = (T_SecondaryWall - WT_c)
        
        if abs(RE1) <= 0.01 and abs(RE2) <= 0.01:
            R_F.magnitude = FoulingResistance(Section)[i] #[cm^2 K/W] (100*100)*(1e-5)
             
            inverseU_total = (PrimaryConvectionResistance(Section, "Dittus-Boetler", PrimaryT_film, i) + ConductionResistance(Section, T_PrimaryWall, i) + \
                    SecondaryConvectionResistance(Section, SecondaryT_film, T_SecondaryWall, i))*OuterArea(Section)[i]\
                    + R_F.magnitude
                    
            U_total.magnitude = 1/inverseU_total #[W/ cm^2 K]
           
            #[W/K cm^2] * [K] * [cm^2] = [W] = [J/s]
            #Q = U_c.magnitude*OuterArea(Section)[i]*(T_PrimaryBulkIn-T_SecondaryWall)    
            Q = h_i.magnitude*InnerArea(Section)[i]*(T_PrimaryBulkIn-T_PrimaryWall)
            
            #Guess Tbh out --> Q hot = Q cold --> calc mass flow steam --> steam frac
            #Guess steam frac --> Qcold --> Qhot --> Tbh out
            SteamFraction = 0.03
            MassFlowSteam = MassFlow_c.magnitude*SteamFraction
            
            Q_c = MassFlowSteam*EnthalpySaturatedSteam.magnitude
            # Q_h = 1718.25*1000*ld.HeatCapacity("PHT", T_PrimaryBulkIn)*(T_PrimaryBulkIn)
                
            T_logmean = Q_c/(3542*U_total.magnitude*(OuterArea(Section)[i]+InnerArea(Section)[i])/2)
        
            #[W] = [W/cm^2 K][cm^2]
            if Section.Length.label[i] == "preheater":
                Area = PreheaterArea(Section)[i]
            else:
                Area = OuterArea(Section)[i]
                
            #T_logmean = Q/(U_total.magnitude*Area)
            
            #print (T_logmean, T_logmeanprime)
            
            return T_PrimaryWall, T_SecondaryWall, T_logmean#T_PrimaryBulkOut, tt
        

    
def TemperatureProfile(Section):
    PrimaryWall = []
    PrimaryBulk = []
    SecondaryBulk = []
    SecondaryWall = []
    
    for i in range(2):    
        if i == 0:
            T_PrimaryBulkIn = 583.15 #[K]
            T_SecondaryBulkIn = 534.15
        
        else: 
            T_PrimaryBulkIn = T_PrimaryBulkOut
            T_SecondaryBulkIn = T_SecondaryBulkOut
            
        PrimaryBulk.append(T_PrimaryBulkIn)    
        SecondaryBulk.append(T_SecondaryBulkIn)
        
        T_wh, T_wc, T_logmean = WallTemperature(Section, i, T_PrimaryBulkIn, T_SecondaryBulkIn)
#       
        PrimaryWall.append(T_wh)
        SecondaryWall.append(T_wc)
        
        if Section.Length.label[i] != "preheater":
                GuessValue = T_PrimaryBulkIn-0.5
                
                T_PrimaryBulkOut =scipy.optimize.fsolve(lambda x: ((x-T_PrimaryBulkIn)/\
            math.log((x-T_SecondaryBulkIn)/(T_PrimaryBulkIn-T_SecondaryBulkIn)))-T_logmean, (GuessValue))[0] # [0] is guess value

                T_SecondaryBulkOut = 534.15
                
        elif Section.Length.label[i] == "preheater":
            T_SecondaryBulkIn = T_SecondaryBulkOut-4.5  
            
            T_PrimaryBulkOut = scipy.optimize.fsolve(lambda x: (((x-T_SecondaryBulkIn)-(T_PrimaryBulkIn-T_SecondaryBulkOut))/\
                                 math.log((x-T_SecondaryBulkIn)/(T_PrimaryBulkIn-T_SecondaryBulkOut)))-T_logmean, (T_PrimaryBulkIn-6))[0]

            
            for k in range(2):
                T_PrimaryBulkIn=T_PrimaryBulkOut    
                
                T_SecondaryBulkIn = scipy.optimize.fsolve(lambda x: (((T_PrimaryBulkOut-x)-(T_PrimaryBulkIn-T_SecondaryBulkOut))/\
                                math.log((T_PrimaryBulkOut-x)/(T_PrimaryBulkIn-T_SecondaryBulkOut)))-T_logmean, (T_SecondaryBulkOut-4.5))[0]
                     
            
                T_wh, T_wc, T_logmean = WallTemperature(Section, i, T_PrimaryBulkIn, T_SecondaryBulkIn)
                #print (T_logmean, T_SecondaryBulkIn-273, i)
                
            
            #print (T_logmean, T_PrimaryBulkIn-273, T_PrimaryBulkOut-273, T_SecondaryBulkIn-273, T_SecondaryBulkOut-273)
                                
            T_SecondaryBulkOut = T_SecondaryBulkIn
        
        
        
        
    print ([j-273.15 for j in PrimaryBulk])
    print ([j-273.15 for j in PrimaryWall])
    print ()
    print ([j-273.15 for j in SecondaryBulk])
    print ([j-273.15 for j in SecondaryWall])

TemperatureProfile(ld.SG_Zone1)

#[g/s]


# m_c = 240*1000 
# Q_c = (1134.78)*240*1000 
# Q_c = ld.HeatCapacity("SHT", 533)*m_c*(260-253)
# 
# m_h = 1719*1000
#             m_c = 240*1000 
#             
#             Q_c = (1134.78)*240*1000 
#             
#             for k in range(10):
#                 Cp_h = ld.HeatCapacity("PHT", (T_PrimaryBulkIn+T_PrimaryBulkOut)/2)
#                 Cp_c = ld.HeatCapacity("SHT", (T_SecondaryBulkIn+T_SecondaryBulkOut)/2)
#                 
#                 Tbh = T_PrimaryBulkOut
#                 TbC = T_SecondaryBulkIn
#                 
#                 T_PrimaryBulkOut = (T_PrimaryBulkIn*(m_h*Cp_h-0.5*U_total*InnerArea(Section)[i])+\
#                                 0.5*U_total*InnerArea(Section)[i]*(T_SecondaryBulkIn+T_SecondaryBulkOut))/(m_h*Cp_h-0.5*U_total*InnerArea(Section)[i])
#         
#                 T_SecondaryBulkIn = ((m_c*Cp_c-0.5*U_total*InnerArea(Section)[i])*T_SecondaryBulkOut-\
#                                 0.5*U_total*InnerArea(Section)[i]*(T_PrimaryBulkIn+T_PrimaryBulkOut))/(m_c*Cp_c-0.5*U_total*InnerArea(Section)[i])
#             
#                 [g/s]*[J/g K]*[K] = [J/s] = [W]
#                 Q_h = ld.HeatCapacity("PHT", T_PrimaryBulkIn)*m_h*(T_PrimaryBulkIn-T_PrimaryBulkOut)
#                 [J/g]*[g/s] = [W]
#                 Q_c = ld.HeatCapacity("SHT", T_SecondaryBulkIn)*m_c*(T_SecondaryBulkIn-T_SecondaryBulkOut)
#               
#                 print (T_PrimaryBulkIn-273, T_PrimaryBulkOut-273)
#                 delta1 = abs(Tbh-T_PrimaryBulkOut)
#                 if delta1<=0.01:
#                     break

