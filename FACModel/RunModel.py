import numpy as np
import NumericConstants as nc
import LepreauData as ld
import Composition as c
import RK4
#import Activities as a
import Electrochemistry as e
import Iteration as it
import matplotlib.pyplot as plt
from matplotlib import rc
rc('mathtext', default='regular')

##Initial Concentrations 
#All concentrations in mol/kg
Sections = [ld.Inlet, ld.Core, ld.Outlet, ld.SteamGenerator] 
for Section in Sections:
    
    Section.FractionFeInnerOxide = c.FractionMetalInnerOxide(Section, "Fe")
    Section.FractionNiInnerOxide = c.FractionMetalInnerOxide(Section, "Ni")
    
    #[Section.Bulk.Co60, Section.Bulk.Co58, Section.Bulk.Fe55, Section.Bulk.Fe59, Section.Bulk.Mn54, Section.Bulk.Cr51, Section.Bulk.Ni63] = [[0]*Section.NodeNumber]*7
    #Section.CorrRate = [10e-9]*Section.NodeNumber
    
    Interfaces = [Section.MetalOxide, Section.SolutionOxide, Section.Bulk]
    for Interface in Interfaces:
        Interface.ConcentrationH2 =[nc.H2*nc.H2Density/nc.H2MolarMass]*Section.NodeNumber #makes into array with appropriate # of nodes for that PHTS section
        Interface.ConcentrationH = np.array(c.BulkpHCalculator(Section)).tolist() #from system pH calculation

        if (Section == ld.Inlet and Interface == Section.SolutionOxide) or (Section == ld.SteamGenerator and Interface == Section.SolutionOxide):
            FeTotal = np.array(Section.SolubilityFe)*0.8
        
        if Section == ld.Core and Interface == Section.MetalOxide:
            Interface.FeTotal = [0]*Section.NodeNumber #makes into array with appropriate # of nodes for that PHTS section
        
        if Section == ld.Outlet and Interface == Section.MetalOxide:
            Interface.FeTotal = [0.00000026]*Section.NodeNumber #From Cook's thesis - experimental corrosion rate measurements and calcs
    
        if Section == ld.Inlet or Section == ld.Core or Section == ld.Outlet:
            Section.MetalOxide.NiTotal = [0]*Section.NodeNumber
        

class RunModel():
    def __init__(self, Section1, Section2, j):
        #j = overall time step
        self.Section1 = Section1
        self.Section2 = Section2
        
        x = []
        y= []
        z = []
        for i in range(self.Section1.NodeNumber):
            if i ==0: 
                self.Section1.Bulk.FeTotal[i] = self.Section1.Bulk.FeTotal[i]
                self.Section1.Bulk.NiTotal[i] = self.Section1.Bulk.NiTotal[i]
                self.Section1.Bulk.CoTotal[i] = self.Section1.Bulk.CoTotal[i]
                
            else:
                self.Section1.Bulk.FeTotal[i] = RK4.Spatial(self.Section1, self.Section1.SolutionOxide.FeTotal[i-1], self.Section1.Bulk.FeTotal[i-1], \
                                                    ld.MassTransfer(self.Section1)[i-1], self.Section1.Diameter[i-1], self.Section1.Velocity[i-1],\
                                                    self.Section1.Length[i-1])
            
                self.Section1.Bulk.NiTotal[i] = RK4.Spatial(self.Section1, self.Section1.SolutionOxide.NiTotal[i-1], self.Section1.Bulk.NiTotal[i-1], \
                                                    ld.MassTransfer(self.Section1)[i-1], self.Section1.Diameter[i-1], self.Section1.Velocity[i-1],\
                                                    self.Section1.Length[i-1])
                 
                self.Section1.Bulk.CoTotal[i] = RK4.Spatial(self.Section1, self.Section1.SolutionOxide.CoTotal[i-1], self.Section1.Bulk.CoTotal[i-1], \
                                                    ld.MassTransfer(self.Section1)[i-1], self.Section1.Diameter[i-1], self.Section1.Velocity[i-1],\
                                                    self.Section1.Length[i-1])
                
                if self.Section1 == ld.Inlet and i ==3:
                    self.Section1.Bulk.FeTotal[i] = 0.59*self.Section1.Bulk.FeTotal[i-1]
                    self.Section1.Bulk.NiTotal[i] = 0.59*self.Section1.Bulk.NiTotal[i-1]
                    self.Section1.Bulk.CoTotal[i] = 0.59*self.Section1.Bulk.CoTotal[i-1]
                
            x.append(self.Section1.Bulk.FeTotal[i])
            y.append(self.Section1.Bulk.NiTotal[i])
            z.append(self.Section1.Bulk.CoTotal[i])
            
        self.Section1.Bulk.FeTotal = x
        self.Section1.Bulk.NiTotal = y
        self.Section1.Bulk.CoTotal = z
        
        end = self.Section1.NodeNumber-1
        [self.Section2.Bulk.FeTotal[0], self.Section2.Bulk.NiTotal[0], self.Section2.Bulk.CoTotal[0], self.Section2.Bulk.CrTotal[0]] = \
        [self.Section1.Bulk.FeTotal[end], self.Section1.Bulk.NiTotal[end], self.Section1.Bulk.CoTotal[end], self.Section1.Bulk.CrTotal[end]]       
        
        if self.Section1 == ld.Core:
            self.Section2.Bulk.CoTotal[0]=self.Section1.Bulk.CoTotal[end]+nc.CobaltWear
            self.Section1.CorrRate, self.Section1.MetalOxide.MixedPotential = [0]*self.Section1.NodeNumber, [0]*self.Section1.NodeNumber
        else:
            #Calculates CS and Alloy-800 corrosion rates based on MO concentrations - for j=0, these are initial concentration values
            #Not called for core- no "corrosion" here
            self.Section1.CorrRate, self.Section1.MetalOxide.MixedPotential =it.CorrosionRate(self.Section1, self.Section1.MetalOxide.FeTotal, \
                                                                self.Section1.MetalOxide.NiTotal, self.Section1.MetalOxide.ConcentrationH)
            
        
        self.Section1.SolutionOxide.MixedPotential, self.Section1.SolutionOxide.EqmPotentialFe3O4 = \
        e.ECP(self.Section1, self.Section1.SolutionOxide.FeTotal, self.Section1.SolutionOxide.FeSatFe3O4, self.Section1.SolutionOxide.NiTotal, \
              self.Section1.SolutionOxide.ConcentrationH)
        
        ##Electrochemical adjustement to Kp, Kd, FeSat and Concentration H+ at M/O interface --> uses M/O and S/O interface mixted potentials
        self.Section1.KpFe3O4electrochem, self.Section1.KdFe3O4electrochem, self.Section1.SolutionOxide.FeSatFe3O4, self.Section1.MetalOxide.ConcentrationH \
        = e.ElectrochemicalAdjustment(self.Section1, self.Section1.SolutionOxide.EqmPotentialFe3O4, self.Section1.SolutionOxide.MixedPotential, \
            self.Section1.MetalOxide.MixedPotential, self.Section1.SolutionOxide.FeTotal, self.Section1.SolutionOxide.FeSatFe3O4, \
            self.Section1.Bulk.FeSatFe3O4, self.Section1.SolutionOxide.ConcentrationH)
        ##
          
        ##Calculates S/O elemental concentrations based on updated oxide thicknesses at each time step
        self.Section1.SolutionOxide.FeTotal = it.SolutionOxide(self.Section1, self.Section1.CorrRate, self.Section1.Bulk.FeTotal, \
        self.Section1.SolutionOxide.FeSatFe3O4, self.Section1.SolutionOxide.FeTotal, self.Section1.InnerIronOxThickness, \
        self.Section1.OuterFe3O4Thickness, self.Section1.NiThickness, self.Section1.CoThickness, self.Section1.KdFe3O4electrochem, \
        self.Section1.KpFe3O4electrochem, "Fe")
            
        self.Section1.SolutionOxide.NiTotal = it.SolutionOxide(self.Section1, self.Section1.CorrRate, self.Section1.Bulk.NiTotal, \
        self.Section1.SolutionOxide.NiSatFerrite, self.Section1.SolutionOxide.NiTotal, self.Section1.InnerIronOxThickness, \
        self.Section1.OuterFe3O4Thickness, self.Section1.NiThickness, self.Section1.CoThickness, self.Section1.KdFe3O4electrochem, \
        self.Section1.KpFe3O4electrochem, "Ni")
            
        self.Section1.SolutionOxide.FeTotal = it.SolutionOxide(self.Section1, self.Section1.CorrRate, self.Section1.Bulk.CoTotal, \
        self.Section1.SolutionOxide.CoSatFerrite, self.Section1.SolutionOxide.CoTotal, self.Section1.InnerIronOxThickness, \
        self.Section1.OuterFe3O4Thickness, self.Section1.NiThickness, self.Section1.CoThickness, self.Section1.KdFe3O4electrochem, \
        self.Section1.KpFe3O4electrochem, "Co")
        ##     
                
        ##RK4 oxide thickness calculation (no spalling)
        self.Section1.InnerOxThickness, self.Section1.OuterOxThickness, self.Section1.InnerIronOxThickness, self.Section1.OuterFe3O4Thickness, \
        self.Section1.CoThickness, self.Section1.NiThickness, self.Section1.CorrRate, self.Section1.SolutionOxide.FeTotal, \
        self.Section1.SolutionOxide.NiTotal, self.Section1.SolutionOxide.CoTotal,  \
        self.Section1.SolutionOxide.MixedPotential, self.Section1.SolutionOxide.EqmPotentialFe3O4, self.Section1.MetalOxide.FeTotal, \
        self.Section1.KdFe3O4electrochem =\
            RK4.OxideGrowth(self.Section1, self.Section1.Bulk.FeTotal, self.Section1.Bulk.NiTotal, self.Section1.Bulk.CoTotal, \
            self.Section1.Bulk.FeSatFe3O4, self.Section1.SolutionOxide.FeSatFe3O4, self.Section1.SolutionOxide.NiSatFerrite, \
            self.Section1.SolutionOxide.CoSatFerrite, self.Section1.SolutionOxide.FeTotal, self.Section1.SolutionOxide.NiTotal, \
            self.Section1.SolutionOxide.CoTotal, self.Section1.OuterOxThickness, self.Section1.InnerOxThickness, self.Section1.InnerIronOxThickness, \
            self.Section1.OuterFe3O4Thickness, self.Section1.NiThickness, self.Section1.CoThickness, self.Section1.KdFe3O4electrochem, \
            self.Section1.KpFe3O4electrochem, self.Section1.SolutionOxide.ConcentrationH, self.Section1.CorrRate)
            
               
            #if self.Section1==ld.Outlet: print (self.Section1.OuterOxThickness)
            #Spalling    
        self.Section1.InnerIronOxThickness, self.Section1.OuterFe3O4Thickness, self.Section1.CoThickness, self.Section1.NiThickness, \
        self.Section1.InnerOxThickness, self.Section1.OuterOxThickness, self.Section1.Particle, self.Section1.ElapsedTime, \
        self.Section1.SpallTime = \
             RK4.Spall(self.Section1, j, self.Section1.Particle, self.Section1.SolutionOxide.FeSatFe3O4, self.Section1.SolutionOxide.FeTotal, \
            self.Section1.KdFe3O4electrochem, self.Section1.OuterOxThickness, self.Section1.InnerOxThickness, self.Section1.OuterFe3O4Thickness, \
            self.Section1.CoThickness, self.Section1.NiThickness, self.Section1.InnerIronOxThickness, self.Section1.ElapsedTime, \
            self.Section1.SpallTime, self.Section1.SolutionOxide.NiTotal, self.Section1.SolutionOxide.CoTotal)
            
                     
        
#         a.InCoreActiveDeposit(self.Section1, self.Section2,self.Section1.InnerOxThickness, self.Section1.OuterOxThickness, self.Section1.OuterFe3O4Thickness, \
#                               self.Section1.NiThickness, self.Section1.CoThickness, self.Section1.InnerOxThickness,j)
#        


##Plot spalltime and maybe elapsed time alongside (or on different plot) with oxide thicknesses? 
##Echem plot...mixed potential versus eqm, especially at S/O interface

InnerOx= []         
import time

start_time = time.time()
for j in range(250):#nc.SimulationDuration
    I = RunModel(ld.Inlet, ld.Core, j)
    C = RunModel(ld.Core, ld.Outlet, j)
    O = RunModel(ld.Outlet, ld.SteamGenerator, j)
    S = RunModel(ld.SteamGenerator, ld.Inlet, j)

#     if (j+1)%100==0:
#         print (O.Section1.ElapsedTime, O.Section1.SpallTime, O.Section1.OuterOxThickness, O.Section1.InnerOxThickness)
#         fig, ax1 = plt.subplots()
#         ax1.plot(S.Section1.InnerOxThickness, linestyle =None, marker='o', color='0.50', label='Inner Oxide')
#         ax1.plot(S.Section1.OuterOxThickness, linestyle =None, marker='o', color='k', label='Outer Oxide')
#         ax1.set_xlabel('Distance (m)')
#         ax1.set_ylabel('Oxide Layer Loadings (${g/m^2}$)')
#           
#         ax2 = ax1.twinx()
#         ax2.plot(S.Section1.NiThickness, linestyle =None, marker='o', color='c', label='Nickel')
#         ax2.plot(S.Section1.CoThickness, linestyle =None, marker='o', color='m', label='Cobalt')
#         ax2.set_ylabel('Ni, Co, and Cr Loadings (${g/m^2}$)', rotation=270, labelpad=20)
#           
#         lines, labels = ax1.get_legend_handles_labels()
#         lines2, labels2 = ax2.get_legend_handles_labels()
#         ax2.legend(lines + lines2, labels + labels2, loc=4)        
#         fig.tight_layout()
#         # plt.axis([4, 69, 0, 10]) 
#         plt.show()

end_time = time.time()
delta_time = end_time - start_time

hours = delta_time//3600
temp = delta_time - 3600*hours
minutes = delta_time//60
seconds = delta_time - 60*minutes
print('%d:%d:%d' %(hours,minutes,seconds))

#     InnerOx.append(O.Section1.InnerOxThickness[2])
# print (InnerOx)

#if self.Section1 == ld.SteamGenerator: print (self.Section1.CorrRate) 
#             [self.Section1.MetalOxide.Co60, self.Section1.MetalOxide.Co58, self.Section1.MetalOxide.Fe55, self.Section1.MetalOxide.Fe59, \
#              self.Section1.MetalOxide.Mn54, self.Section1.MetalOxide.Cr51, self.Section1.MetalOxide.Ni63] = \
#              a.SurfaceActivity(self.Section1, self.Section1.CorrRate, self.Section1.SolutionOxide.FeSatFe3O4, self.Section1.InnerOxThickness, \
#                                BulkActivities, j)
#    

       #Initial activity at the start (first node) of each section 
#         BulkConcentration_o = [self.Section1.Bulk.Co60[1], self.Section1.Bulk.Co58[1], self.Section1.Bulk.Fe55[1], self.Section1.Bulk.Fe59[1], \
#                                              self.Section1.Bulk.Mn54[1], self.Section1.Bulk.Cr51[1], self.Section1.Bulk.Ni63[1]]
#             
#         
        
                 
#         BulkActivities = np.array([self.Section1.Bulk.Co60, self.Section1.Bulk.Co58, self.Section1.Bulk.Fe55, self.Section1.Bulk.Fe59, self.Section1.Bulk.Mn54, self.Section1.Bulk.Cr51, self.Section1.Bulk.Ni63])
#    










LoopDistance = I.Section1.Length + C.Section1.Length + O.Section1.Length + S.Section1.Length
TotalLoopDistance = [i/100 for i in np.cumsum(LoopDistance)] #Distance down length of PHTS [m]

BulkFeSat = [np.log10(i) for i in I.Section1.Bulk.FeSatFe3O4 + C.Section1.Bulk.FeSatFe3O4 + O.Section1.Bulk.FeSatFe3O4 + S.Section1.Bulk.FeSatFe3O4]
BulkFeTotal = [np.log10(i) for i in I.Section1.Bulk.FeTotal + C.Section1.Bulk.FeTotal + O.Section1.Bulk.FeTotal + S.Section1.Bulk.FeTotal]
SolutionOxideFeTotal = [np.log10(i) for i in I.Section1.SolutionOxide.FeTotal + C.Section1.SolutionOxide.FeTotal + \
                        O.Section1.SolutionOxide.FeTotal + S.Section1.SolutionOxide.FeTotal]
SolutionOxideFeSat = [np.log10(i) for i in I.Section1.SolutionOxide.FeSatFe3O4 + C.Section1.SolutionOxide.FeSatFe3O4 + \
                        O.Section1.SolutionOxide.FeSatFe3O4 + S.Section1.SolutionOxide.FeSatFe3O4]


BulkNiTotal = [np.log10(i) for i in I.Section1.Bulk.NiTotal + C.Section1.Bulk.NiTotal + O.Section1.Bulk.NiTotal + S.Section1.Bulk.NiTotal]
SolutionOxideNiTotal = [np.log10(i) for i in I.Section1.SolutionOxide.NiTotal + C.Section1.SolutionOxide.NiTotal + \
                        O.Section1.SolutionOxide.NiTotal + S.Section1.SolutionOxide.NiTotal]
SolutionOxideNiSat = [np.log10(i) for i in I.Section1.SolutionOxide.NiSatFerrite + C.Section1.SolutionOxide.NiSatFerrite + \
                        O.Section1.SolutionOxide.NiSatFerrite + S.Section1.SolutionOxide.NiSatFerrite]


BulkCoTotal = [np.log10(i) for i in I.Section1.Bulk.CoTotal + C.Section1.Bulk.CoTotal + O.Section1.Bulk.CoTotal + S.Section1.Bulk.CoTotal]
SolutionOxideCoTotal = [np.log10(i) for i in I.Section1.SolutionOxide.CoTotal + C.Section1.SolutionOxide.CoTotal + \
                        O.Section1.SolutionOxide.CoTotal + S.Section1.SolutionOxide.CoTotal]
SolutionOxideCoSat = [np.log10(i) for i in I.Section1.SolutionOxide.CoSatFerrite + C.Section1.SolutionOxide.CoSatFerrite + \
                        O.Section1.SolutionOxide.CoSatFerrite + S.Section1.SolutionOxide.CoSatFerrite]


Corrosion= []
Rates = [I.Section1.CorrRate, C.Section1.CorrRate, O.Section1.CorrRate, S.Section1.CorrRate]
for Rate, Section in zip (Rates, Sections):
    x = ld.UnitConverter(Section, "Corrosion Rate Grams", "Corrosion Rate Micrometers", None, Rate, None, None, None)
    Corrosion.append(x)
CorrosionRate = Corrosion[0] + Corrosion[1] + Corrosion[2] + Corrosion[3]


InnerOxide = []
InnerOxideThicknesses = [I.Section1.InnerOxThickness, C.Section1.InnerOxThickness, O.Section1.InnerOxThickness, S.Section1.InnerOxThickness]
for Thickness, Sect in zip (InnerOxideThicknesses, Sections):
    z = ld.UnitConverter(Sect, "Grams per Cm Squared", "Grams per M Squared", None, None, Thickness, None, None)
    InnerOxide.append(z)
InnerOxideThickness = InnerOxide[0]+InnerOxide[1]+InnerOxide[2]+InnerOxide[3]


OuterOxide = []
OuterOxideThicknesses = [I.Section1.OuterOxThickness, C.Section1.OuterOxThickness, O.Section1.OuterOxThickness, S.Section1.OuterOxThickness]
for Thickness, Sect in zip (OuterOxideThicknesses, Sections):
    q = ld.UnitConverter(Sect, "Grams per Cm Squared", "Grams per M Squared", None, None, Thickness, None, None)
    OuterOxide.append(q)
OuterOxideThickness = OuterOxide[0]+OuterOxide[1]+OuterOxide[2]+OuterOxide[3]

#LogCorrosion= [np.log10(i) for i in Corrosion]


# plt.plot(TotalLoopDistance, CorrosionRate, 'ro')
# #plt.axis([4, 69, -2, 3]) 
# plt.show()


plt.plot(TotalLoopDistance, SolutionOxideNiTotal, 'ro')
plt.plot(TotalLoopDistance, SolutionOxideNiSat, '-')
plt.plot(TotalLoopDistance, BulkNiTotal, 'g-')
 #plt.axis([4, 69, -2, 3]) 
plt.show()

plt.plot(TotalLoopDistance, SolutionOxideCoTotal, 'ro')
plt.plot(TotalLoopDistance, SolutionOxideCoSat, '-')
plt.plot(TotalLoopDistance, BulkCoTotal, 'g-')
#plt.axis([4, 69, -2, 3]) 
plt.show()

 
plt.plot(TotalLoopDistance, SolutionOxideFeTotal, 'ro')
plt.plot(TotalLoopDistance, SolutionOxideFeSat, '-')
plt.plot(TotalLoopDistance, BulkFeTotal, 'g-')
 #plt.axis([4, 69, -2, 3]) 
plt.show()
print (CorrosionRate)
  
plt.plot(TotalLoopDistance, InnerOxideThickness, 'ro')
plt.plot(TotalLoopDistance, OuterOxideThickness, 'p-')
#plt.axis([4, 69, 0, 200]) 
plt.show()




# def OxideThicknessConverter(OxideType, OxideThicknesses, Sects):
#     InnerOxide = []
#     OxideDensity = []
#    
#     for Thickness, Sect, Section in zip(OxideThicknesses, Sects, [ld.Inlet, ld.Core, ld.Outlet, ld.SteamGenerator]):
#         if OxideType == "Inner":
#             for i in range(Sect.NodeNumber):
#                 if Sect.OuterOxThickness == 0 and (Sect.CoThickness[i] + Sect.NiThickness[i] > 0):
#                     x = nc.NiFe2O4Density
#                 elif Sect == S.Section1: 
#                     x = nc.NiFe2O4Density
#                 else:
#                     x = nc.Fe3O4Density
#         
#             OxideDensity.append(x)
#             
#             if Sect == C.Section1:
#                 y = 0
#             else:
#                 y1 = [c.FractionChromite(Section)*i for i in ld.UnitConverter(Sect, "Oxide Thickness Grams", "Oxide Thickness Micrometers", None, None, \
#                                                                           Thickness, [nc.FeCr2O4Density]*Sect.NodeNumber, None)]
#                 y2 = [(1-c.FractionChromite(Section))*i for i in ld.UnitConverter(Sect, "Oxide Thickness Grams", "Oxide Thickness Micrometers", None, None, \
#                                                                           Thickness, OxideDensity, None)]
#             
#                 y = [j+k for j,k in zip(y1, y2)]
#             
#             InnerOxide.append(y)
#         
#         
#         else:
#             if Sect.CoThickness[i] + Sect.NiThickness[i] > 0:
#                 q = nc.NiFe2O4Density
#             else:
#                 q=nc.Fe3O4Density
#             
#             OxideDensity.append(q)
#             
#             z= ld.UnitConverter(Sect, "Oxide Thickness Grams", "Oxide Thickness Micrometers", None, None, Thickness, OxideDensity, None)
#             
#             InnerOxide.append(z)
#                
#     return InnerOxide
# 
# print (OxideThicknessConverter("Inner", InnerOxideThicknesses, [I.Section1, C.Section1, O.Section1, S.Section1]))





#             
#             self.Section1.KpFe3O4electrochem, self.Section1.KdFe3O4electrochem, self.Section1.SolutionOxide.FeSatFe3O4, self.Section1.MetalOxide.ConcentrationH \
#             = e.ElectrochemicalAdjustment(self.Section1, self.Section1.SolutionOxide.EqmPotentialFe3O4, self.Section1.SolutionOxide.MixedPotential, \
#                                     self.Section1.MetalOxide.MixedPotential, self.Section1.SolutionOxide.FeTotal, self.Section1.SolutionOxide.FeSatFe3O4, \
#                                     self.Section1.Bulk.FeSatFe3O4, self.Section1.SolutionOxide.ConcentrationH)
#             
#             if self.Section1 == ld.SteamGenerator or self.Section1 == ld.Inlet or self.Section1 == ld.Outlet:
#                 self.Section1.MetalOxide.FeTotal = it.MetalOxideInterfaceConcentration(self.Section1, "Fe", self.Section1.SolutionOxide.FeTotal, \
#                     self.Section1.InnerOxThickness, self.Section1.OuterOxThickness, self.Section1.CorrRate)
#         
#             
#                 #Calculates CS and Alloy-800 corrosion rates based on MO concentrations - for j=1, these are initial concentration values
#                 #Not called for core- no "corrosion" here
#                 self.Section1.CorrRate, self.Section1.MetalOxide.MixedPotential =it.CorrosionRate(self.Section1, self.Section1.MetalOxide.FeTotal, \
#                                                                 self.Section1.MetalOxide.NiTotal, self.Section1.MetalOxide.ConcentrationH)
#             
#             else:
#                 self.Section1.CorrRate, self.Section1.MetalOxide.MixedPotential = [0]*self.Section1.NodeNumber, [0]*self.Section1.NodeNumber  
#             
#         
#             
#         ##