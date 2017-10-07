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
MassFlow_h = ld.SGParameters()

for i in [h_i, h_o]:
    i.unit = "K/W"
    
for i in [U_h, U_c, U_total]:
    i.unit = "W/cm^2 K"

R_F.unit = "cm^2 K/W"
MassFlow_c.magnitude = 240*1000 
MassFlow_c.unit = "g/s"
MassFlow_h.magnitude = 1719.25*1000/3542
MassFlow_h.unit = "g/s"
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
    T_sat = 260.1 + 273.15 #[K]
    
    if Section.Length.label[i] == "preheater" or Section.Length.label[i] == "preheater start": #from Silpsrikul thesis
        f_b = 0.1783 #fraction of cross-sectional area of shell occupied by a baffle window
        NumberTubes = 3542
        N_b = f_b*NumberTubes #number of tubes in baffle window 
        S_b = 0.1271*(100**2) #[cm^2]  
        G_e = 1512*1000/(100**2) #[g/cm^2 s] weighted average mass velocity in preheater 
        
        #First two (raised to exponent) terms are unitless 
        #[W/cm K]/[cm] = [W/cm^2 K]
        h_o = (ThermalConductivity(Tfilm, "water")*0.36/Section.OuterDiameter[i])*\
        ((Section.OuterDiameter[i]*G_e/ld.Viscosity("water", "SHT", Tfilm))**0.6)*\
        (ld.HeatCapacity("SHT", Tfilm)*ld.Viscosity("water", "SHT", Tfilm)/ThermalConductivity(Tfilm, "water"))**0.36
    
    else:
        x=Section.Length.steam_quality[i]*100
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
   


def NusseltNumber(Section, correlation, Temperature,i):
    Re_D = ld.ReynoldsNumber(Section, Section.Diameter)
    
#     MassFlow_h.magnitude = 1719.25*1000/3542 #[g/s]
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
        return [i*1 for i in [np.pi*x*y for x,y in zip([2.89]*Section.NodeNumber, [i*1 for i in Section.Length.magnitude])]]#[1*1 for i in [np.pi*x*y for x,y in zip(Section.OuterDiameter,Section.Length)]]
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
            #print (U_total.magnitude)
#             print (100*100*h_i.magnitude, "h_ih", 100*100*h_o.magnitude, "h_oc", \
#                (100*100)/(InnerArea(Section)[i]*ConductionResistance(Section, T_PrimaryWall, i)), i,"k_cond", "unit for all = [m^2 K/W]")
            
            return T_PrimaryWall, T_SecondaryWall, U_total.magnitude#T_PrimaryBulkOut, tt
            

    
def TemperatureProfile(Section):
    PrimaryWall = []
    PrimaryBulk = []
    SecondaryBulk = []
    SecondaryWall = []
    
    for i in range(Section.NodeNumber):    
        if i == 0:
            T_PrimaryBulkIn = 583.15 #[K]
            T_SecondaryBulkIn = 533.25
            m_in = 0
        
        PrimaryBulk.append(T_PrimaryBulkIn)    
        SecondaryBulk.append(T_SecondaryBulkIn)
        
        T_wh, T_wc, U = WallTemperature(Section, i, T_PrimaryBulkIn, T_SecondaryBulkIn)
        
        PrimaryWall.append(T_wh)
        SecondaryWall.append(T_wc)
        
        Cp_h=ld.HeatCapacity("PHT", T_PrimaryBulkIn)
        
        if Section.Length.label[i] != "preheater" or Section.Length.label != "preheater start":
            T_PrimaryBulkOut = T_PrimaryBulkIn-(U*(T_PrimaryBulkIn-T_SecondaryBulkIn)*OuterArea(Section)[i])/(Cp_h*MassFlow_h.magnitude)
                
            T_SecondaryBulkOut = 533.25
            #for one tube x 3542 = for all tubes
            Q1 = 3542*MassFlow_h.magnitude*Cp_h*(T_PrimaryBulkIn-T_PrimaryBulkOut)
                
            #for one tube (U based on one tube) (heat transfer area x number tubes) --> would cancel out in U calc (1/h*NA) NA
            Q = U*(T_PrimaryBulkIn-T_SecondaryBulkIn)*OuterArea(Section)[i]*3542
                
            Q3 = MassFlow_h.magnitude*3500*ld.HeatCapacity("PHT", T_PrimaryBulkIn)
                
            m_out = m_in + Q/(MassFlow_c.magnitude*EnthalpySaturatedSteam.magnitude)
            
        if Section.Length.label[i] == "preheater" or Section.Length.label[i] == "preheater start":
            if Section.Length.label[i] == "preheater start":
                T_SecondaryBulkOut = 260.1+273.15
                
            Cp_c = ld.HeatCapacity("SHT", T_SecondaryBulkIn)
            inverseC_c = 1/(Cp_c*MassFlow_c.magnitude)
            inverseC_h = 1/(Cp_h*MassFlow_h.magnitude*3452)
            
            T_PrimaryBulkOut = T_PrimaryBulkIn-(U*OuterArea(Section)[i]*3542/(Cp_h*MassFlow_h.magnitude*3542))*(T_PrimaryBulkIn-T_SecondaryBulkIn)
                
                
            T_SecondaryBulkOut = T_SecondaryBulkIn-(U*OuterArea(Section)[i]*3542/(Cp_c*MassFlow_c.magnitude))*(T_PrimaryBulkIn-T_SecondaryBulkIn)
            print (T_SecondaryBulkOut-273,i)   
            #deltaT_out = (U*OuterArea(Section)[i]*3542*(inverseC_c-inverseC_h)+1)*(T_PrimaryBulkIn-T_SecondaryBulkIn)
            
        m_in = m_out
        T_PrimaryBulkIn = T_PrimaryBulkOut
        T_SecondaryBulkIn = T_SecondaryBulkOut
        
        #print (T_SecondaryBulkOut-273.15,T_PrimaryBulkIn-273.15, i)    

        
    print ([j-273.15 for j in PrimaryBulk])
    print ([j-273.15 for j in PrimaryWall])
    print ()
    print ([j-273.15 for j in SecondaryBulk])
    print ([j-273.15 for j in SecondaryWall])

TemperatureProfile(ld.SG_Zone1)


