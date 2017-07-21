import LepreauData as ld
import NumericConstants as nc
import Composition as c
import numpy as np
import Iteration as it 
import Electrochemistry as e
import random


def Particulate(Section, ErosionConstant, BulkCrud):
    if Section == ld.Core:
        DepositionConstant = nc.Kdeposition_InCore
    else:
        DepositionConstant = nc.Kdeposition_OutCore

    Deposition = [DepositionConstant*(4/x)*100/(y*1000) for x,y in zip(Section.Diameter,Section.DensityH2O)]
    
    Crud = [((ErosionConstant*100*100)/DepositionConstant)*(1-np.exp((-x/y)*z))+q*np.exp((-x/y)*z) for x,y,z,q in zip(Deposition,Section.Velocity, Section.Distance, BulkCrud)]
    return Crud


def OxideComposition(Section, Element, OxideType, Outer, OuterFe3O4Thickness, NiThickness, CoThickness, InnerIronOxThickness):     
#Note: this function does not take list input 
    
    if Outer == "yes": #for determining outer layer only composition
        #Outer Fe3O4 layer 
        if OxideType >0: #For the total oxide value in each node
                #Fraction of iron in Fe3O4 * fraction of magnetite in outer oxide
            if Element == "Fe":
                return OuterFe3O4Thickness*c.FractionMetalInOxide(Section, Element, "Magnetite")/OxideType
                #Outer ox thickness includes Fe3O4 and any incorporated Ni and Co
            if Element == "Ni":
                return NiThickness/OxideType 
            if Element == "Co":
                return CoThickness/OxideType
            if Element == "Cr":
                #Cr assumed to be 100% retained in inner layer
                return 0
        else:
            return 0#composition within outer layer only 
         
    if Outer == "no": #inner oxide layer composition only
        if OxideType >0: #For the InnerIronOxThickness value in each node
            if Element == "Fe":
                return 4#c.FractionMetalInnerOxide(Section, Element)             
            if OuterFe3O4Thickness > 0:
                if Element == "Ni":    
                    return c.FractionMetalInnerOxide(Section, Element) #no Ni or Co incorporation via precip - only if present in inner layer through diffusion
                if Element == "Co":
                    return c.FractionMetalInnerOxide(Section, Element)
            else: #InnerOxThickness includes Ni incorporation from precipitation 
                if Element == "Ni":
                    return (NiThickness +c.FractionMetalInnerOxide(Section, Element)*InnerIronOxThickness)/OxideType  
                if Element == "Co":
                    return (CoThickness+c.FractionMetalInnerOxide(Section, Element)*InnerIronOxThickness)/OxideType
            if Element == "Cr":
                return c.FractionMetalInnerOxide(Section, "Cr")
     
            if OxideType == 0: #InnerIronOxide = 0 
                return 0


def Spatial(Section, Solution, Bulk, km, Diameter, Velocity, Length):
    #Cylindrical pipe: Surface Area / Volume = Pi*Diameter*Length/(Pi*(Diameter^2)/4) = 4/Diameter
    #Allows conversion between amount/area (assuming uniform distribution across the area) to amount/volume
    Solution_Bulk = Solution - Bulk
    Delta = (km*4/(Velocity*Diameter))*Solution_Bulk#[(4*x/(y*z))*q for x,y,z,q in zip(km, Section.Diameter, Section.Velocity, Solution_Bulk)]
    #[cm]*[1/cm]*[mol/kg] + [mol/kg] = [mol/kg]
    return Bulk+Delta*Length#[x + y*Length for x,y in zip(Bulk, Delta)] 

 
def OxideGrowth(Section, BulkFeTotal, BulkNiTotal, BulkCoTotal, BulkFeSat, SolutionOxideFeSat, SolutionOxideNiSat, SolutionOxideCoSat, \
                SolutionOxideFeTotal, SolutionOxideNiTotal, SolutionOxideCoTotal, OuterOxThickness, InnerOxThickness, InnerIronOxThickness,\
                 OuterFe3O4Thickness, NiThickness, CoThickness, KdFe3O4electrochem, KpFe3O4electrochem, ConcentrationH, CorrRate):    
    
    
    RK4_InnerIronOxThickness = InnerIronOxThickness
    RK4_OuterFe3O4Thickness = OuterFe3O4Thickness
    RK4_NiThickness = NiThickness
    RK4_CoThickness = CoThickness
#     if Section==ld.Inlet:
#         import matplotlib.pyplot as plt
#         from matplotlib import rc
#         rc('mathtext', default='regular')
#         fig, ax1 = plt.subplots()
#         ax1.plot(InnerOxThickness, linestyle =None, marker='o', color='0.50', label='Inner Oxide')
#         ax1.plot(OuterOxThickness, linestyle =None, marker='o', color='k', label='Outer Oxide')
#         ax1.set_xlabel('Distance (m)')
#         ax1.set_ylabel('Oxide Layer Loadings (${g/m^2}$)')
#          
#         ax2 = ax1.twinx()
#         ax2.plot(NiThickness, linestyle =None, marker='o', color='c', label='Nickel')
#         ax2.plot(CoThickness, linestyle =None, marker='o', color='m', label='Cobalt')
#         ax2.set_ylabel('Ni, Co, and Cr Loadings (${g/m^2}$)', rotation=270, labelpad=20)
#          
#         lines, labels = ax1.get_legend_handles_labels()
#         lines2, labels2 = ax2.get_legend_handles_labels()
#         ax2.legend(lines + lines2, labels + labels2, loc=0)        
#         fig.tight_layout()
#         # plt.axis([4, 69, 0, 10]) 
#         plt.show()
    
    L = []
    M = []
    N = []
    O = []
    #if Section ==ld.Outlet: print ()
    for approximation in range(4): #4 approximations in the "RK4" method       
        ##Solves S/O elemental concentrations at current approximation of oxide thickness(es)
        #At start of new time step, input S/O concentrations based on previous time step's evaluation (needed inside SolutionOxideBalance function to 
        #determine if <> saturation)
                        
        SolutionOxideFeTotal = it.SolutionOxide(Section, CorrRate, BulkFeTotal, SolutionOxideFeSat, SolutionOxideFeTotal, RK4_InnerIronOxThickness, \
                                                RK4_OuterFe3O4Thickness, NiThickness, CoThickness, KdFe3O4electrochem, KpFe3O4electrochem, "Fe")
        
        SolutionOxideNiTotal = it.SolutionOxide(Section, CorrRate, BulkNiTotal, SolutionOxideNiSat, SolutionOxideNiTotal, RK4_InnerIronOxThickness, \
                                                RK4_OuterFe3O4Thickness, NiThickness, CoThickness, KdFe3O4electrochem, KpFe3O4electrochem, "Ni")
            
        SolutionOxideCoTotal = it.SolutionOxide(Section, CorrRate, BulkCoTotal, SolutionOxideCoSat, SolutionOxideCoTotal, RK4_InnerIronOxThickness, \
                                                RK4_OuterFe3O4Thickness, NiThickness, CoThickness, KdFe3O4electrochem, KpFe3O4electrochem, "Co")
        ##
        #print (SolutionOxideFeTotal)
        MixedECP, EqmPotentialFe3O4 = e.ECP(Section, SolutionOxideFeTotal, SolutionOxideFeSat, SolutionOxideNiTotal, ConcentrationH)
        
        if Section == ld.SteamGenerator:
            MetalOxideNi = it.MetalOxideInterfaceConcentration(Section, "Ni", SolutionOxideNiTotal, InnerOxThickness, OuterOxThickness, CorrRate)
        else:
            MetalOxideNi = [0]*Section.NodeNumber
            #MetalOxideCo = it.MetalOxideInterfaceConcentration(Section, "Co", SolutionOxideCoTotal, InnerOxThickness, OuterOxThickness, CorrRate)
    
        if Section == ld.Core:
            CorrRate = [0]*Section.NodeNumber
            MetalOxideFe = [0]*Section.NodeNumber
            MOMixedECP = [0]*Section.NodeNumber
        else:
            
            MetalOxideFe = it.MetalOxideInterfaceConcentration(Section, "Fe", SolutionOxideFeTotal, InnerOxThickness, OuterOxThickness, CorrRate)
            
            CorrRate, MOMixedECP = it.CorrosionRate(Section, MetalOxideFe, MetalOxideNi, ConcentrationH)
            #if Section==ld.Outlet: print (MetalOxideFe, ld.UnitConverter(Section, "Corrosion Rate Grams", "Corrosion Rate Micrometers", None, CorrRate, None, None, None))
        
        KpFe3O4electrochem, KdFe3O4electrochem, SolutionOxideFeSat, MOConcentrationH = e.ElectrochemicalAdjustment(Section, EqmPotentialFe3O4, \
                                            MixedECP, MOMixedECP, SolutionOxideFeTotal, SolutionOxideFeSat, BulkFeSat, ConcentrationH)
        
        
        ##Converts all concentration units to be used within RK4 from [mol/kg] to [g/cm^3]
        Concentrations = [SolutionOxideFeTotal, SolutionOxideFeSat, SolutionOxideNiTotal, SolutionOxideNiSat, SolutionOxideCoTotal, SolutionOxideCoSat]
        MolarMasses = [nc.FeMolarMass, nc.FeMolarMass, nc.NiMolarMass, nc.NiMolarMass, nc.CoMolarMass, nc.CoMolarMass]
        ConvertedConcentrations = []
        for i,h in zip (Concentrations, MolarMasses): #Concentrations has 6 lists in it
            x = [(t/1000)*(h*y) for t, y in zip(i, Section.DensityH2O)]
            ConvertedConcentrations.append(x)
        FeTotal, FeSat, NiTotal, NiSat, CoTotal, CoSat = ConvertedConcentrations
        ##
        
        ##Determines dy/dt (growth functions) for each type of solid corrosion product at each node based on respective element's saturation behaviour
        GrowthOuterMagnetite = []
        GrowthCobalt = []  
        GrowthNickel = []
        GrowthInnerIronOxide = []
        
        ##Iterate through this using previously solved RK4 thickness: re-evaluates growth functions based on S/O + M/O concentrations using new thickness     
        for i in range(Section.NodeNumber):                 
            #Magnetite (outer and inner)
            if RK4_OuterFe3O4Thickness[i] > 0:
                #With outer layer present, dFe3O4_i/dt only depends on corrosion rate/diffusion (function of Fe3O4_i thickness)
                if Section == ld.Core:
                    q = 0
                else:
                    #While outer layer present, inner Fe3O4 not affected by oxide kinetics (only diffusion)
                    q = CorrRate[i]*it.Diffusion(Section, "Fe")/Section.FractionFeInnerOxide #dy/dt = f(y)
                if FeTotal[i] >= FeSat[i]: #Precipitation of outer layer magnetite 
                #dy/dt = f(y) via Diff term in SO concentration
                    x = KpFe3O4electrochem[i]*(FeTotal[i] - FeSat[i])
                else:#if FeTotal[i] < FeSat[i]:
                    x = KdFe3O4electrochem[i]*(FeTotal[i] - FeSat[i])
            else: #OuterFe3O4Thickness[i] ==0
                if FeTotal[i] >= FeSat[i]:
                    if Section == ld.Core:
                        q = 0
                    else:
                        q = CorrRate[i]*it.Diffusion(Section, "Fe")/Section.FractionFeInnerOxide 
                    x = KpFe3O4electrochem[i]*(FeTotal[i] - FeSat[i])
                else: #if FeTotal[i] < FeSat[i]:
                    x = 0 #nothing to dissolve
                    if Section == ld.Core:
                        q = 0
                    else:
                        q = CorrRate[i]*it.Diffusion(Section, "Fe")/Section.FractionFeInnerOxide + KdFe3O4electrochem[i]*(FeTotal[i] - FeSat[i])
                        
            GrowthOuterMagnetite.append(x)            
            GrowthInnerIronOxide.append(q)
            
            #Cobalt incorporation
            if RK4_InnerIronOxThickness[i] > 0:#Oxide layer present for Co to incorporate into
                if CoTotal[i] >= CoSat[i]:
                    y = KpFe3O4electrochem[i]*(CoTotal[i] - CoSat[i]) 
                else: #CoTotal[i] < CoSat[i]
                    if RK4_CoThickness[i] >0:
                        y = KdFe3O4electrochem[i]*(CoTotal[i] - CoSat[i])
                    else: #CoThickness == 0
                        y= 0 #Nothing to dissolve 
            else: #InnerIronOxThickness[i] == 0 (Core)
                y = 0 #Nothing to incorporate into 
           
            GrowthCobalt.append(y)           
           
            #Nickel incorporation or Ni(s) deposition (independent of presence of an oxide layer backbone)
            if NiTotal[i] >= NiSat[i]:
                z= KpFe3O4electrochem[i]*(NiTotal[i] - NiSat[i]) 
            else: #NiTotal[i] < NiSat[i]:
                if RK4_NiThickness[i] > 0: #Ni must be present for its dissolution to occur
                    z= KdFe3O4electrochem[i]*(NiTotal[i] - NiSat[i])
                else:
                    z = 0 #Nothing to dissolve
            GrowthNickel.append(z)
        ##  
        
        [RK4_InnerIronOxThickness, a] = RK4(Section, InnerIronOxThickness, GrowthInnerIronOxide, approximation)
        L.append(a)
        
        [RK4_OuterFe3O4Thickness, b] = RK4(Section, OuterFe3O4Thickness, GrowthOuterMagnetite, approximation)
        M.append(b)
        
        [RK4_CoThickness, c] = RK4(Section, CoThickness, GrowthCobalt, approximation)
        N.append(c)
        
        [RK4_NiThickness, d] = RK4(Section, NiThickness, GrowthNickel, approximation)
        O.append(d)
        ##
        
    InnerIronOxThickness = [a+(b+2*c+2*d+e)/6 for a, b, c, d, e in zip(InnerIronOxThickness, L[0], L[1], L[2], L[3])]
    OuterFe3O4Thickness = [a+(b+2*c+2*d+e)/6 for a, b, c, d, e in zip(OuterFe3O4Thickness, M[0], M[1], M[2], M[3])]
    CoThickness = [a+(b+2*c+2*d+e)/6 for a, b, c, d, e in zip(CoThickness, N[0], N[1], N[2], N[3])]
    NiThickness = [a+(b+2*c+2*d+e)/6 for a, b, c, d, e in zip(NiThickness, O[0], O[1], O[2], O[3])]
    
    Layers = [InnerIronOxThickness, OuterFe3O4Thickness, CoThickness, NiThickness]
    #4 different layers (potentially) at each node. If any oxide thicknesses are negative due to dissolution of respective layer, thickness = 0
    for i in range(4):
        for x in range(Section.NodeNumber):
            if Layers[i][x] < 0:
                Layers[i][x] = 0       

    
    return InnerOxThickness, OuterOxThickness, InnerIronOxThickness, OuterFe3O4Thickness, CoThickness, NiThickness, CorrRate, SolutionOxideFeTotal, \
        SolutionOxideNiTotal, SolutionOxideCoTotal, MixedECP, EqmPotentialFe3O4, MetalOxideFe, KdFe3O4electrochem

    
def RK4(Section, InitialThickness, GrowthFunction, approximation):
    x = [i*nc.TimeIncrement  for i in GrowthFunction] #f(t)*delta t
    if approximation <2:
        thickness= [a+b/2 for a,b in zip(InitialThickness,x)] #y = 
    elif approximation==2:
        thickness = [a+b for a,b in zip(InitialThickness,x)]
    if approximation > 2:
        thickness = InitialThickness
    
    for i in range(Section.NodeNumber):
        if thickness[i] <0:
            thickness[i]= 0
        
    return thickness, x


def ParticleSize():
    r = random.random()
    if 0<r<0.005:
        Size = 0.15
    elif 0.005<r<0.009:
        Size = 0.5
    elif 0.009<r<0.17:
        Size = 1
    elif 0.17<r<0.25:
        Size = 1.5
    elif 0.25<r<0.31:
        Size = 2
    elif 0.31<r<0.4:
        Size = 2.5
    elif 0.4<r<0.5:
        Size = 3.2
    elif 0.5<r<0.52:
        Size = 3.2
    elif 0.52<r<0.55:
        Size = 4
    elif 0.55<r<0.6:
        Size = 4.5
    elif 0.6<r<0.62:
        Size = 5
    elif 0.62<r<0.67:
        Size = 5.5
    elif 0.67<r<0.69:
        Size = 6
    elif 0.69<r<0.72:
        Size = 7
    elif 0.72<r<0.75:
        Size = 8
    elif 0.75<r<0.8:
        Size = 9
    elif 0.8<r<0.83:
        Size = 10
    elif 0.83<r<0.9:
        Size = 15
    elif 0.9<r<0.98:
        Size = 40
    elif 0.98<r<1:
        Size = 100
    else:
        print (r)
    
    return Size*0.0001*nc.Fe3O4Density #[um to cm to g/cm^2]


def SpallingTime(Section, Particle, SolutionOxideFeSat, SolutionOxideFeTotal, KdFe3O4electrochem, OuterOxThickness,Velocity):
    
    if SolutionOxideFeSat > SolutionOxideFeTotal: #Outlet 
        if OuterOxThickness > 0:
            SpTime = (nc.OutletOuterSpallConstant*Particle/3600)/((Velocity**2)*nc.Fe3O4Porosity_outer*KdFe3O4electrochem*(SolutionOxideFeSat-SolutionOxideFeTotal))
        else: #no outer oxide layer 
            SpTime = (nc.OutletInnerSpallConstant*Particle/3600)/((Velocity**2)*nc.Fe3O4Porosity_inner*KdFe3O4electrochem*(SolutionOxideFeSat-SolutionOxideFeTotal))
    else: #Inlet 
        if OuterOxThickness > 0:
            SpTime = (nc.InletOuterSpallConstant * Particle / 3600) / (Velocity**2)
        else: #no outer oxide layer 
            SpTime = (nc.InletInnerSpallConstant * Particle / 3600) / (Velocity**2)
        
    return SpTime


def Oxide(Layer, TotalOxideThickness, Spalling):
    #inner/outer oxthickness & spalling = [g/cm^2]
    #This assumes that solid corrosion product species (Fe3O4, Ni0.6Fe2.4O4, etc.), if present simultaneously, are uniformally distributed with depth/layer thickness
    if TotalOxideThickness > 0:
        Ox = Layer-(Layer/TotalOxideThickness)*Spalling 
        if Ox < 0:
            Ox = 0
        return Ox


def Spall(Section,j, Particle, SolutionOxideFeSat, SolutionOxideFeTotal, KdFe3O4electrochem, OuterOxThickness, InnerOxThickness, OuterFe3O4Thickness, \
          CoThickness, NiThickness, InnerIronOxThickness, ElapsedTime, SpallTime, SolutionOxideNiTotal, SolutionOxideCoTotal):
    ##Current time step's RK4 input from OxideGrowth function for each of InnerOx, InnerIronOx, OuterFe3O4, Co, and Ni at each node of current section 
    
    #Silences spalling for desired sections
    if Section == ld.SteamGenerator or Section == ld.Core or Section==ld.Inlet:
        Particle = [0]*Section.NodeNumber 
    
    ConvertedConcentrations = []
    Concentrations = [SolutionOxideFeSat, SolutionOxideFeTotal]
    for i in range(2):
        x = ld.UnitConverter(Section, "Mol per Kg", "Grams per Cm Cubed", Concentrations[i], None, None, None, nc.FeMolarMass)
        ConvertedConcentrations.append(x)
    
    FeSat, FeTotal = ConvertedConcentrations
    
#     for i in range(Section.NodeNumber):    
#         if Section== ld.Core:
#             InnerIronOxThickness[i] =0
#                 
#         if OuterFe3O4Thickness[i] > 0: #from previous time step
#             #With outer magnetite layer present, Ni and Co incorporate into overall "outer" oxide layer  
#             OuterOxThickness[i] = OuterFe3O4Thickness[i] + CoThickness[i] + NiThickness[i]
#             InnerOxThickness[i] = InnerIronOxThickness[i]
#         else: #OuterFe3O4Thickness == 0
#             InnerOxThickness[i] = InnerIronOxThickness[i] + CoThickness[i] + NiThickness[i]
    
    for i in range(Section.NodeNumber):    
#         if Section==ld.Outlet or Section==ld.Inlet or Section==ld.SteamGenerator:
#             if InnerIronOxThickness[i] <= 1e-6:
#                 InnerIronOxThickness[i] = 0.00025 #Resets to original thickness (resetting corrosion rate here to ~75 um/a
#                 #if Section==ld.SteamGenerator:
#                     #print (i, Section, "too much dissolution")
            
        ## Oxide totals for RK4 iterations (M/O Concentration depends on total oxide thickness) and before spalling function
        if OuterFe3O4Thickness[i]> 0: #from previous time step
            #With outer magnetite layer present, Ni and Co incorporate into overall "outer" oxide layer  
            OuterOxThickness[i] = OuterFe3O4Thickness[i] + CoThickness[i] + NiThickness[i]
            InnerOxThickness[i] = InnerIronOxThickness[i]
        else: #OuterFe3O4Thickness == 0
            InnerOxThickness[i] = InnerIronOxThickness[i] + CoThickness[i] + NiThickness[i]
        
        
        
        
        
        if j == 0:
        #First time step call generate particle sizes and calc spalling times, respectively
            x = ParticleSize()
            #Amount of time it will take each node's particle to spall 
            y = SpallingTime(Section, x, FeSat[i], FeTotal[i], KdFe3O4electrochem[i], OuterFe3O4Thickness[i], Section.Velocity[i])
            
            Particle.append(x) #list of random particle sizes based on input distribution in ParticleSize function
            SpallTime.append(y)
            
#         if Section == ld.Outlet: #Outlet header spalling corrector 
#             if SpallTime[7] > 4000:
#                 SpallTime[7] = 3000 
            ElapsedTime = [0]*Section.NodeNumber #No time has elapsed yet at first time step for all nodes    
        
        else: #after first time step
            if ElapsedTime[i] >= SpallTime[i]: #enough time elapsed for particle of that size (with respective spall time) to come off
                if OuterOxThickness[i] >0: 
                    #[cm]*[g/cm^3] = [g/cm^2] of layer removed due to spalling
                    OuterFe3O4Thickness[i] = Oxide(OuterFe3O4Thickness[i], OuterOxThickness[i], Particle[i]) 
                    CoThickness[i] =  Oxide(CoThickness[i], OuterOxThickness[i], Particle[i])
                    NiThickness[i] = Oxide(NiThickness[i], OuterOxThickness[i], Particle[i])
                    OuterOxThickness[i]  = Oxide(OuterOxThickness[i], OuterOxThickness[i], Particle[i])#new total outer ox
                else: #OuterOXThickness == 0, spalling takes place at inner layer  
                    InnerIronOxThickness[i]= Oxide(InnerIronOxThickness[i], InnerOxThickness[i], Particle[i])
                    CoThickness[i] = Oxide(CoThickness[i], InnerOxThickness[i], Particle[i])
                    NiThickness[i] = Oxide(NiThickness[i], InnerOxThickness[i], Particle[i])
                    InnerOxThickness[i] = Oxide(InnerOxThickness[i], InnerOxThickness[i], Particle[i])#new total inner ox
                    
                    if InnerIronOxThickness[i] <=1e-7: 
                        InnerIronOxThickness[i] = 0.00025
                        #print (i, Section, "help")
                        
                Particle[i]= ParticleSize() #once particle spalled off, another particle size (just at that node) is randomly generated (w/i distribution)
            
                SpallTime[i] = SpallingTime(Section, Particle[i], FeSat[i], FeTotal[i], KdFe3O4electrochem[i], OuterFe3O4Thickness[i], Section.Velocity[i])
                if Section == ld.Outlet: #Outlet header spalling corrector 
                    if SpallTime[7] >= 4000:
                        SpallTime[7] = 3000              
                ElapsedTime[i] = 0 #once a particle has "spalled" off, elapsed time since spalling resets to zero and counter restarts at that node
        
            else: #not enough time has passed 
                ElapsedTime[i] = ElapsedTime[i] + 1*(nc.TimeIncrement/3600)   
    
    
    for i in range(Section.NodeNumber):    
        if Section==ld.Outlet or Section==ld.Inlet or Section==ld.SteamGenerator:
            if InnerIronOxThickness[i] <= 1e-6:
                InnerIronOxThickness[i] = 0.00025 #Resets to original thickness (resetting corrosion rate here to ~75 um/a
                #if Section==ld.SteamGenerator:
                    #print (i, Section, "too much dissolution")
            
            ## Oxide totals for RK4 iterations (M/O Concentration depends on total oxide thickness) and before spalling function
            if OuterFe3O4Thickness[i]> 0: #from previous time step
                #With outer magnetite layer present, Ni and Co incorporate into overall "outer" oxide layer  
                OuterOxThickness[i] = OuterFe3O4Thickness[i] + CoThickness[i] + NiThickness[i]
                InnerOxThickness[i] = InnerIronOxThickness[i]
            else: #OuterFe3O4Thickness == 0
                InnerOxThickness[i] = InnerIronOxThickness[i] + CoThickness[i] + NiThickness[i]
    
    return InnerIronOxThickness, OuterFe3O4Thickness, CoThickness, NiThickness, InnerOxThickness, OuterOxThickness, Particle, ElapsedTime, SpallTime


