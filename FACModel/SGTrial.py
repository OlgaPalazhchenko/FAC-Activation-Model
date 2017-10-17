import LepreauData as ld
import numpy as np
import NumericConstants as nc

T_sat = 260.1 + 273.15
T_PreheaterIn = 186.5 + 273.15

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
EquivalentDiameter = ld.SGParameters()
TubePitch = ld.SGParameters()

for i in [h_i, h_o]:
    i.unit = "K/W"
    
for i in [U_h, U_c, U_total]:
    i.unit = "W/cm^2 K"

R_F.unit = "cm^2 K/W"

MassFlow_c.magnitude = 239.25*1000 
MassFlow_h.magnitude = 1719.25*1000

ShellDiameter.magnitude = 2.28*100 #[cm]
for i in [MassFlow_c, MassFlow_h]:
    i.unit = "g/s"

for i in [ShellDiameter, EquivalentDiameter,TubePitch]:
    i.unit = "cm"
    
MassFlux.magnitude = MassFlow_c.magnitude/((np.pi/4)*(ShellDiameter.magnitude**2)) #[g/cm^2 s]
MassFlux.unit = "g/cm^2 s"
EnthalpySaturatedSteam.magnitude = 2800.4
EnthalpySaturatedSteam.unit = "J/g"

TubePitch.magnitude = 2.413

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

def SecondaryConvectionResistance(Section, Tfilm, Twall, x_in, i):
    #split into boiling and non-boiling (preheater) sections 

    EquivalentDiameter.magnitude = 4*((TubePitch.magnitude**2)-(np.pi*(Section.OuterDiameter[i]**2))/4)/(np.pi*Section.OuterDiameter[i])
    
    
    if Section.Length.label[i] == "preheater" or Section.Length.label[i] == "preheater start": #from Silpsrikul thesis
        #Based on PLGS data and Silpsrikul calculations 
        f_b = 0.2119 #fraction of cross-sectional area of shell occupied by a baffle window
        NumberTubes = 3542
        N_b = f_b*NumberTubes #number of tubes in baffle window 
        D_s = ((4*0.5*np.pi/4)*(ShellDiameter.magnitude/100)**2)/(0.5*np.pi*(ShellDiameter.magnitude/100)+ShellDiameter.magnitude/100)#[m]
         
        S_b = (f_b*np.pi*(D_s**2)/4)-N_b*np.pi*((Section.OuterDiameter[i]/100)**2)/4 #[cm^2]  
        G_b = (MassFlow_c.magnitude/1000)/S_b#[kg/m^2s]
        S_c = 0.5373*D_s*(1-(Section.OuterDiameter[i]/100)/0.0219)
        G_c = MassFlow_c.magnitude/1000/S_c
        G_e = np.sqrt(G_b*G_c)*1000/(100**2)       
        #G_e = 1512*1000/(100**2) #[g/cm^2 s] weighted average mass velocity in preheater 
        
        #First two (raised to exponent) terms are unitless 
        #[W/cm K]/[cm] = [W/cm^2 K]
        h_o = (ThermalConductivity(Tfilm, "water")*0.25/Section.OuterDiameter[i])*\
        ((Section.OuterDiameter[i]*G_e/ld.Viscosity("water", "SHT", Tfilm))**0.6)*\
        (ld.HeatCapacity("SHT", Tfilm)*ld.Viscosity("water", "SHT", Tfilm)/ThermalConductivity(Tfilm, "water"))**0.33
    
    else:
        if i<=10:
            x= x_in*100
        else:
            x=2
        
        rho_v = 1000*23.753/(100**3) #[g/cm^3]
        p_crit = 22.0640 #[MPa]
        #V*rho =[cm/s]([cm^2]*[g/cm^3] = [g/s]
        
        F = (1+(x*Prandtl(Section, T_sat, i)*((ld.Density("water", "SHT", T_sat)/rho_v)-1)))**0.35
       
        Re_D = EquivalentDiameter.magnitude*MassFlux.magnitude/ld.Viscosity("water", "SHT", T_sat)#4*MassFlow/((ld.Viscosity("water", "SHT", T_sat)/1000)*np.pi*Shell_ID)
        
        h_l = 0.023*ld.ThermalConductivityH2O("SHT", T_sat)*((Re_D)**0.8)*(Prandtl(Section, T_sat, i)**0.4)/EquivalentDiameter.magnitude #[W/cm^2 K]
        
        Q_prime =F*h_l*(Twall-T_sat) #[W/cm^2]
        S = (1+0.055*(F**0.1)*(Re_D)**0.16)**(-1) 
        A_p = (55*((4.70/p_crit)**0.12)*((-np.log10(4.70/p_crit))**(-0.55))*(nc.H2MolarMass)**(-0.5))/(100**2)
        C = ((A_p*S)/(F*h_l)**2)*(Q_prime)**(4/3)
        coeff = [1, -C, -1]
        cubic_solution = np.roots(coeff)
        q = cubic_solution[0]
        
        h_o = F*(q**(3/2))*h_l #[W/cm^2 K
       
    return 1/(h_o*OuterArea(Section)[i]) #K/W
   

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
    return Pr


def InnerArea(Section):
    return [np.pi*x*y for x,y in zip(Section.Diameter,Section.Length.magnitude)]


def OuterArea(Section):
    return [np.pi*x*y for x,y in zip(Section.OuterDiameter,Section.Length.magnitude)]


def WallTemperature(Section, i, T_PrimaryBulkIn, T_SecondaryBulkIn, x_in):
    #i = each node of SG tube 
    T_PrimaryWall = T_PrimaryBulkIn-(1/3)*(T_PrimaryBulkIn-T_SecondaryBulkIn) #perhaps call from outside function
    T_SecondaryWall =T_PrimaryBulkIn
    
    for k in range(50):        
        WT_h = T_PrimaryWall
        WT_c = T_SecondaryWall
        
        PrimaryT_film = (T_PrimaryBulkIn + T_PrimaryWall)/2
        SecondaryT_film = (T_SecondaryBulkIn + T_SecondaryWall)/2

        h_i.magnitude = 1/(PrimaryConvectionResistance(Section, "Dittus-Boetler", PrimaryT_film, i)*InnerArea(Section)[i])
        h_o.magnitude = 1/(SecondaryConvectionResistance(Section, SecondaryT_film, T_SecondaryWall, x_in, i)*OuterArea(Section)[i])
        #R = 1/hA --> h=A*1/R
        
        U_h1 = ConductionResistance(Section, T_PrimaryWall, i)*InnerArea(Section)[i]
        U_h2 = SecondaryConvectionResistance(Section, SecondaryT_film, T_SecondaryWall, x_in, i)*InnerArea(Section)[i]
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
                    SecondaryConvectionResistance(Section, SecondaryT_film, T_SecondaryWall, x_in, i))*OuterArea(Section)[i]\
                    
                    
            U_total.magnitude = 1/inverseU_total #[W/ cm^2 K]
#             if i ==19:print (U_total.magnitude, PrimaryConvectionResistance(Section, "Dittus-Boetler", PrimaryT_film, i), \
#                    SecondaryConvectionResistance(Section, SecondaryT_film, T_SecondaryWall, x_in, i), ConductionResistance(Section, T_PrimaryWall, i))
            
            return T_PrimaryWall, T_SecondaryWall, U_total.magnitude
            

    
def TemperatureProfile(Section):
    PrimaryWall = []
    PrimaryBulk = []
    SecondaryBulk = []
    SecondaryWall = []
    
    for i in range(0,19):
        if i == 0:
            #Temperatures entering SG --> not first node temps.
            T_PrimaryBulkIn = 583.15 #[K]
            T_SecondaryBulkIn = 533.25
            x_in = 0
        
        if Section.Length.label[i] != "preheater start" and Section.Length.label[i] != "preheater":    
        
            Cp_h=ld.HeatCapacity("PHT", T_PrimaryBulkIn)
            Cp_c = ld.HeatCapacity("SHT", T_SecondaryBulkIn)
            T_wh, T_wc, U = WallTemperature(Section, i, T_PrimaryBulkIn, T_SecondaryBulkIn, x_in)
            
            T_PrimaryBulkOut = T_PrimaryBulkIn-(U*(T_PrimaryBulkIn-T_SecondaryBulkIn)*OuterArea(Section)[i])/(Cp_h*MassFlow_h.magnitude/3542)
            #Ti = Ti+1 = Tsat    
            T_SecondaryBulkOut = T_sat
            
            #for one tube --> multiply by 3542 = for all tubes
            #Q = MassFlow_h.magnitude*Cp_h*(T_PrimaryBulkIn-T_PrimaryBulkOut)

            #for one tube (U based on one tube) (heat transfer area x number tubes) --> would cancel out in U calc (1/h*NA) NA
            Q = U*(T_PrimaryBulkIn-T_SecondaryBulkIn)*OuterArea(Section)[i]*3542    
            #Q = MassFlow_h.magnitude*ld.HeatCapacity("PHT", T_PrimaryBulkIn)
           
            x_out = x_in + Q/(MassFlow_c.magnitude*EnthalpySaturatedSteam.magnitude)
            
            x_in = x_out
            T_PrimaryBulkIn = T_PrimaryBulkOut
            T_SecondaryBulkIn = T_SecondaryBulkOut
            
            PrimaryBulk.append(T_PrimaryBulkIn)   
            SecondaryBulk.append(T_SecondaryBulkIn)
            PrimaryWall.append(T_wh)
            SecondaryWall.append(T_wc)
            
        else:
            if Section.Length.label[i] == "preheater start":
                T_wh, T_wc, U = WallTemperature(Section, i, T_PrimaryBulkIn, T_SecondaryBulkIn, x_in)
                #Tcold,out = Ti = Tsat = 260.1 oC (533.25 K)
                T_SecondaryBulkIn = 264+273 
                #TotalArea= sum(OuterArea(Section)[i:Section.NodeNumber]) 
    
                #Tcold,i+1 - Tcold,i
                dQ_c = MassFlow_c.magnitude*Cp_c*(T_SecondaryBulkIn-T_PreheaterIn)
                #Thot,i+1 = Thot,i 
                T_PrimaryBulkOutEnd = T_PrimaryBulkIn - dQ_c/(Cp_h*MassFlow_h.magnitude) #estimate based on heat balance only  
                
                #T_LMTD = ((T_PrimaryBulkIn-T_SecondaryBulkIn)- (T_PrimaryBulkOut-T_PreheaterIn))/(np.log(T_PrimaryBulkIn-T_SecondaryBulkIn)-np.log(T_PrimaryBulkOut-T_PreheaterIn))
                #F = 0.8725 #https://checalc.com/solved/LMTD_Chart.html (cross-flow LMTD correction factor)
                print (T_PrimaryBulkOutEnd-273.15)
             
    print ([j-273.15 for j in PrimaryBulk])
    print ([j-273.15 for j in PrimaryWall])
    print ()
    print ([j-273.15 for j in SecondaryBulk])
    print ([j-273.15 for j in SecondaryWall])
    #18 element arrays
        
    for j in range(1):
#         PrimaryBulk2=[]
#         SecondaryBulk2=[]
        for i in range(18, Section.NodeNumber):                 
            Cp_h=ld.HeatCapacity("PHT", T_PrimaryBulkIn)
            Cp_c = ld.HeatCapacity("SHT", T_SecondaryBulkIn)
            if i == 18: 
                T_PrimaryBulkIn =PrimaryBulk[17]
                 
            T_PrimaryBulkOut = T_PrimaryBulkIn-(U*OuterArea(Section)[i]*3542/(Cp_h*MassFlow_h.magnitude))*(T_PrimaryBulkIn-T_SecondaryBulkIn)
            T_SecondaryBulkOut = T_SecondaryBulkIn-(U*OuterArea(Section)[i]*3542/(Cp_c*MassFlow_c.magnitude))*(T_PrimaryBulkIn-T_SecondaryBulkIn)
            
            print (T_PrimaryBulkOut-273.15,T_SecondaryBulkOut-273.15,i,j)   
            T_PrimaryBulkIn = T_PrimaryBulkOut
            T_SecondaryBulkIn = T_SecondaryBulkOut
            
#             PrimaryBulk2.append(T_PrimaryBulkOut)
#             SecondaryBulk2.append(T_SecondaryBulkOut)
            
            
            if i == Section.NodeNumber-1: #21
                dQ_h = MassFlow_h.magnitude*Cp_h*(PrimaryBulk[17]-T_PrimaryBulkOut)
                T_SecondaryBulkOutEnd = T_PreheaterIn + dQ_h/(Cp_c*MassFlow_c.magnitude)
                print (T_SecondaryBulkOutEnd-273.15,i)
            
                T_SecondaryBulkIn = T_SecondaryBulkOutEnd
        
            #T_wh, T_wc, U = WallTemperature(Section, i, T_PrimaryBulkIn, T_SecondaryBulkIn, x_in)
            #Tcold, i+1 starts at 186.5 (preheater entrance temperature end) and work backwards to 
            #Thot,i --> starts at 268-272 hot in end and works down to outlet (requires 260.1 Tcold,i assumption) 
            #don't call wall temp calcs, assume basic U stays the same 
            #U only adjusted in these nodes via oxide thicknesses 
        
#            T_PrimaryBulkOut = T_PrimaryBulkIn-(U*OuterArea(Section)[i]*3542/(Cp_h*MassFlow_h.magnitude))*(T_LMTD)*F
#            T_SecondaryBulkOut = T_SecondaryBulkIn-(U*OuterArea(Section)[i]*3542/(Cp_c*MassFlow_c.magnitude))*(T_LMTD)*F
                      
        

TemperatureProfile(ld.SG_Zone1)

