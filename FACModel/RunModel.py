import numpy as np
import NumericConstants as nc
import LepreauData as ld
import Composition as c
import RK4
import Activities as a
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
    
    [Section.Bulk.Co60, Section.Bulk.Co58, Section.Bulk.Fe55, Section.Bulk.Fe59, Section.Bulk.Mn54, Section.Bulk.Cr51, Section.Bulk.Ni63] \
    = [[1e-10]*Section.NodeNumber]*7
    
    Interfaces = [Section.MetalOxide, Section.SolutionOxide, Section.Bulk]
    for Interface in Interfaces:
        Interface.ConcentrationH2 =[nc.H2*nc.H2Density/nc.H2MolarMass]*Section.NodeNumber #makes into array with appropriate # of nodes for that PHTS section
        Interface.ConcentrationH = c.BulkpHCalculator(Section) #from system pH calculation
            
    Section.SolutionOxide.MixedPotential, Section.SolutionOxide.EqmPotentialFe3O4 = \
    e.ECP(Section)
    
    if Section==ld.Core:
        Section.CorrRate = [0]*Section.NodeNumber
    else:
        Section.CorrRate, Section.MetalOxide.MixedPotential = it.CorrosionRate(Section)
        
    Section.KpFe3O4electrochem, Section.KdFe3O4electrochem, Section.SolutionOxide.FeSatFe3O4, Section.MetalOxide.ConcentrationH = \
    e.ElectrochemicalAdjustment(Section, Section.SolutionOxide.EqmPotentialFe3O4, Section.SolutionOxide.MixedPotential, \
                                Section.MetalOxide.MixedPotential, Section.SolutionOxide.FeTotal, Section.SolutionOxide.FeSatFe3O4, \
                                Section.Bulk.FeSatFe3O4, Section.SolutionOxide.ConcentrationH)
##


class RunModel():
    def __init__(self, Section1, Section2, j): #j = overall time step
        self.Section1 = Section1
        self.Section2 = Section2
        
        #BulkActivities = [self.Section1.Bulk.Co60, self.Section1.Bulk.Co58, self.Section1.Bulk.Fe59, self.Section1.Bulk.Fe55, \
        #                  self.Section1.Bulk.Mn54, self.Section1.Bulk.Cr51, self.Section1.Bulk.Ni63]
        
        #BulkActivities2 = [self.Section2.Bulk.Co60, self.Section2.Bulk.Co58, self.Section2.Bulk.Fe59, self.Section2.Bulk.Fe55, \
        #                  self.Section2.Bulk.Mn54, self.Section2.Bulk.Cr51, self.Section2.Bulk.Ni63]
        
        #SurfaceActivities = [self.Section1.Bulk.Co60, self.Section1.Bulk.Co58, self.Section1.Bulk.Fe59, self.Section1.Bulk.Fe55, \
        #                  self.Section1.Bulk.Mn54, self.Section1.Bulk.Cr51, self.Section1.Bulk.Ni63]
        
        #Tags = ["Co60", "Co58", "Fe59", "Fe55", "Mn54", "Cr51", "Ni63"]
        
        BulkConcentrations = [self.Section1.Bulk.FeTotal, self.Section1.Bulk.NiTotal, self.Section1.Bulk.CoTotal, self.Section1.Bulk.CrTotal]
        BulkConcentrations2 = [self.Section2.Bulk.FeTotal, self.Section2.Bulk.NiTotal, self.Section2.Bulk.CoTotal, self.Section2.Bulk.CrTotal]
        
        SolutionOxideConcentrations = [self.Section1.SolutionOxide.FeTotal, self.Section1.SolutionOxide.NiTotal, \
                                       self.Section1.SolutionOxide.CoTotal, self.Section1.SolutionOxide.CrTotal]
        Saturations = [self.Section1.SolutionOxide.FeSatFe3O4, self.Section1.SolutionOxide.NiSatFerrite, self.Section1.SolutionOxide.CoSatFerrite]
        
        for i in range(self.Section1.NodeNumber):
            if i >0: 
                for r,q in zip (BulkConcentrations, SolutionOxideConcentrations):
                    r[i] = RK4.Spatial(self.Section1, q[i-1], r[i-1], ld.MassTransfer(self.Section1)[i-1], self.Section1.Diameter[i-1], \
                                       self.Section1.Velocity[i-1], self.Section1.Length.magnitude[i-1])

                #Exponential decay of bulk particulate at start of section as function of distance + removal due to deposition and erosion source
                self.Section1.BigParticulate[i] = RK4.Particulate(self.Section1, self.Section1.BigParticulate[0], self.Section1.Diameter[i],\
                                                self.Section1.DensityH2O[i], self.Section1.Velocity[i], self.Section1.Distance[i])
                self.Section1.SmallParticulate[i] = RK4.Particulate(Section, self.Section1.SmallParticulate[0], self.Section1.Diameter[i],\
                                                self.Section1.DensityH2O[i], self.Section1.Velocity[i], self.Section1.Distance[i])
                                
#                 for x,y in zip (BulkActivities, Tags):
#                     x[i] = a.BulkActivity(self.Section1, self.Section2, x[0], self.Section1.CorrRate, self.Section1.Bulk.FeSatFe3O4, \
#                                 self.Section1.InnerOxThickness, self.Section1.OuterOxThickness, self.Section1.OuterFe3O4Thickness, \
#                                 self.Section1.NiThickness, self.Section1.CoThickness, self.Section1.InnerIronOxThickness, \
#                                 y, self.Section1.BigParticulate, self.Section1.SmallParticulate, j, i)
                    
                ##Inlet header purification system
                if self.Section1 == ld.Inlet:      
                    if i ==3:
                        #for x,y in zip(BulkConcentrations, BulkActivities):
                        for x in BulkConcentrations: 
                            x[i] = 0.59*x[i-1]
                            #y[i] = 0.59*y[i-1]
                    
                        self.Section1.SmallParticulate[i] = 0.59*self.Section1.SmallParticulate[i-1]
                        self.Section1.BigParticulate[i] = 0 #0.45 um filter removes everything over this size
                ##    
                    
                    if i >3: #decay for the rest of the inlet section depends on purification system 
                        self.Section1.BigParticulate[i] = RK4.Particulate(Section, self.Section1.BigParticulate[3], self.Section1.Diameter[i], \
                                                        self.Section1.DensityH2O[i], self.Section1.Velocity[i], self.Section1.Distance[i])
                        self.Section1.SmallParticulate[i] = RK4.Particulate(Section, self.Section1.SmallParticulate[3], self.Section1.Diameter[i], \
                                                        self.Section1.DensityH2O[i], self.Section1.Velocity[i], self.Section1.Distance[i]) 
                        
#                         for x,y in zip (BulkActivities, Tags):
#                             #print (y, x[i], i, self.Section1.NodeNumber)
#                             x[i] = a.BulkActivity(self.Section1, self.Section2, x[3], self.Section1.CorrRate, self.Section1.Bulk.FeSatFe3O4, \
#                                     self.Section1.InnerOxThickness, self.Section1.OuterOxThickness, self.Section1.OuterFe3O4Thickness, \
#                                     self.Section1.NiThickness, self.Section1.CoThickness, self.Section1.InnerIronOxThickness, y, \
#                                     self.Section1.BigParticulate, self.Section1.SmallParticulate, j, i)
                        
                    
        end = self.Section1.NodeNumber-1
        for x,y in zip(BulkConcentrations,BulkConcentrations2):
        #for x,y,z,q in zip(BulkConcentrations,BulkConcentrations2, BulkActivities, BulkActivities2):
            y[0] = x[end]
            #q[0] = z[end]
        
        self.Section2.BigParticulate[0], self.Section2.SmallParticulate[0] = self.Section1.BigParticulate[end], self.Section1.SmallParticulate[end]       
        
        ##Stellite wear bulk input term for cobalt and chromium    
        if self.Section1 == ld.Core:
            self.Section2.Bulk.CoTotal[0]=self.Section1.Bulk.CoTotal[end]+nc.CobaltWear
            self.Section2.Bulk.CrTotal[0]=self.Section1.Bulk.CrTotal[end]+nc.CobaltWear*(nc.FractionCr_Stellite/nc.FractionCo_Stellite)
        ##
        
#         else:
#             ##Surface activities 
#             self.Section1.MetalOxide.Co60 = a.SurfaceActivity(self.Section1, self.Section1.CorrRate, self.Section1.SolutionOxide.FeSatFe3O4, \
#                                                           self.Section1.InnerOxThickness, self.Section1.Bulk.Co60, j, "Co60")
#             
#             self.Section1.MetalOxide.Co58 = a.SurfaceActivity(self.Section1, self.Section1.CorrRate, self.Section1.SolutionOxide.FeSatFe3O4, \
#                                                           self.Section1.InnerOxThickness, self.Section1.Bulk.Co58, j, "Co58")
#             
#             self.Section1.MetalOxide.Fe59 = a.SurfaceActivity(self.Section1, self.Section1.CorrRate, self.Section1.SolutionOxide.FeSatFe3O4, \
#                                                           self.Section1.InnerOxThickness, self.Section1.Bulk.Fe59, j, "Fe59")
#             
#             self.Section1.MetalOxide.Fe55 = a.SurfaceActivity(self.Section1, self.Section1.CorrRate, self.Section1.SolutionOxide.FeSatFe3O4, \
#                                                           self.Section1.InnerOxThickness, self.Section1.Bulk.Fe55, j, "Fe55")
#             
#             self.Section1.MetalOxide.Mn54 = a.SurfaceActivity(self.Section1, self.Section1.CorrRate, self.Section1.SolutionOxide.FeSatFe3O4, \
#                                                           self.Section1.InnerOxThickness, self.Section1.Bulk.Mn54, j, "Mn54")
#             
#             self.Section1.MetalOxide.Cr51 = a.SurfaceActivity(self.Section1, self.Section1.CorrRate, self.Section1.SolutionOxide.FeSatFe3O4, \
#                                                           self.Section1.InnerOxThickness, self.Section1.Bulk.Cr51, j, "Cr51")
#             
#             self.Section1.MetalOxide.Ni63 = a.SurfaceActivity(self.Section1, self.Section1.CorrRate, self.Section1.SolutionOxide.FeSatFe3O4, \
#                                                           self.Section1.InnerOxThickness, self.Section1.Bulk.Ni63, j, "Ni63")
            ##    
            
        ##Deposit thickness around PHTS
        self.Section1.DepositThickness = a.Deposition(Section1, self.Section1.BigParticulate, self.Section1.SmallParticulate, j)
        ##
                        
        ##RK4 oxide thickness calculation (no spalling)
        RK4.OxideGrowth(self.Section1, Saturations, BulkConcentrations)
        ##
                
        ##Spalling    
        self.Section1.ElapsedTime, self.Section1.SpallTime = RK4.Spall(self.Section1, j, self.Section1.ElapsedTime, self.Section1.SpallTime)
        ##
        

##Echem plot...mixed potential versus eqm, especially at S/O interface
##Reduce default oxide thickness from 0.00025 g/cm^2 to less than that? Average corrosion rate? 
#What is the minimum allowable in outlet? SG?
    
TotalTime = []
SACo60 = [] 
SACo58 = [] 
SAFe59 = []
SAFe55 = []
SAMn54 = []
SACr51 = []
SANi63 = []       

import time
start_time = time.time()
for j in range(70):#nc.SimulationDuration
    I = RunModel(ld.Inlet, ld.Core, j)
    C = RunModel(ld.Core, ld.Outlet, j)
    O = RunModel(ld.Outlet, ld.SteamGenerator, j)
    S = RunModel(ld.SteamGenerator, ld.Inlet, j)
    
#     TotalTime.append(j/7000)# converts to EFPY
#     SACo60.append(sum(S.Section1.MetalOxide.Co60)/S.Section1.NodeNumber)
#     SACo58.append(sum(S.Section1.MetalOxide.Co58)/S.Section1.NodeNumber)
#     SAFe59.append(sum(S.Section1.MetalOxide.Fe59)/S.Section1.NodeNumber)
#     SAFe55.append(sum(S.Section1.MetalOxide.Fe55)/S.Section1.NodeNumber)
#     SAMn54.append(sum(S.Section1.MetalOxide.Mn54)/S.Section1.NodeNumber)
#     SACr51.append(sum(S.Section1.MetalOxide.Cr51)/S.Section1.NodeNumber)
#     SANi63.append(sum(S.Section1.MetalOxide.Ni63)/S.Section1.NodeNumber)
    
end_time = time.time()
delta_time = end_time - start_time

hours = delta_time//3600
temp = delta_time - 3600*hours
minutes = delta_time//60
seconds = delta_time - 60*minutes
print('%d:%d:%d' %(hours,minutes,seconds))

ActivityTime = TotalTime[0::3]
Co60SGActivity = SACo60[0::3]
Co58SGActivity = SACo58[0::3]
Fe59SGActivity = SAFe59[0::3]
Fe55SGActivity = SAFe55[0::3]
Mn54SGActivity = SAMn54[0::3]
Cr51SGActivity = SACr51[0::3]
Ni63SGActivity = SANi63[0::3]

LoopDistance = I.Section1.Length.magnitude + C.Section1.Length.magnitude + O.Section1.Length.magnitude + S.Section1.Length.magnitude
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


BulkCrTotal = [np.log10(i) for i in I.Section1.Bulk.CrTotal + C.Section1.Bulk.CrTotal + O.Section1.Bulk.CrTotal + S.Section1.Bulk.CrTotal]

SolutionOxideCrSat = [np.log10(i) for i in I.Section1.SolutionOxide.CrSat + C.Section1.SolutionOxide.CrSat + O.Section1.SolutionOxide.CrSat + \
                      S.Section1.SolutionOxide.CrSat]


Corrosion= []
Rates = [I.Section1.CorrRate, C.Section1.CorrRate, O.Section1.CorrRate, S.Section1.CorrRate]
for Rate, Section in zip (Rates, Sections):
    x = ld.UnitConverter(Section, "Corrosion Rate Grams", "Corrosion Rate Micrometers", None, Rate, None, None, None, None)
    Corrosion.append(x)
CorrosionRate = Corrosion[0] + Corrosion[1] + Corrosion[2] + Corrosion[3]

InnerOxide = []
InnerOxideThicknesses = [I.Section1.InnerOxThickness, C.Section1.InnerOxThickness, O.Section1.InnerOxThickness, S.Section1.InnerOxThickness]
for Thickness, Sect in zip (InnerOxideThicknesses, Sections):
    z = ld.UnitConverter(Sect, "Grams per Cm Squared", "Grams per M Squared", None, None, Thickness, None, None, None)
    InnerOxide.append(z)
InnerOxideThickness = InnerOxide[0]+InnerOxide[1]+InnerOxide[2]+InnerOxide[3]

OuterOxide = []
OuterOxideThicknesses = [I.Section1.OuterOxThickness, C.Section1.OuterOxThickness, O.Section1.OuterOxThickness, S.Section1.OuterOxThickness]
for Thickness, Sect in zip (OuterOxideThicknesses, Sections):
    q = ld.UnitConverter(Sect, "Grams per Cm Squared", "Grams per M Squared", None, None, Thickness, None, None, None)
    OuterOxide.append(q)
OuterOxideThickness = OuterOxide[0]+OuterOxide[1]+OuterOxide[2]+OuterOxide[3]

Cobalt = []
CoThicknesses = [I.Section1.CoThickness, C.Section1.CoThickness, O.Section1.CoThickness, S.Section1.CoThickness]
for Thickness, Sect in zip (CoThicknesses, Sections):
    c = ld.UnitConverter(Sect, "Grams per Cm Squared", "Grams per M Squared", None, None, Thickness, None, None, None)
    Cobalt.append(c)
CobaltThickness = Cobalt[0]+Cobalt[1]+Cobalt[2]+Cobalt[3]

Nickel = []
NiThicknesses = [I.Section1.NiThickness, C.Section1.NiThickness, O.Section1.NiThickness, S.Section1.NiThickness]
for Thickness, Sect in zip (NiThicknesses, Sections):
    n = ld.UnitConverter(Sect, "Grams per Cm Squared", "Grams per M Squared", None, None, Thickness, None, None, None)
    Nickel.append(n)
NickelThickness = Nickel[0]+Nickel[1]+Nickel[2]+Nickel[3]


#LogCorrosion= [np.log10(i) for i in Corrosion]


fig1 = plt.figure()
ax1 = fig1.add_subplot(221)
ax1.plot(TotalLoopDistance, SolutionOxideFeTotal, marker='o', color='r', label = 'S/O Interface')
ax1.plot(TotalLoopDistance, SolutionOxideFeSat, marker='o', color='w', label = 'Magnetite Solubility')
ax1.plot(TotalLoopDistance, BulkFeTotal, marker='o', color='0.50', label = 'Bulk Coolant')
ax1.set_xlabel('Distance (m)')
ax1.set_ylabel('$Log_{10}$[Fe Concentration (mol/kg)]')
#ax1.legend()

ax2 = fig1.add_subplot(222)
ax2.plot(TotalLoopDistance, SolutionOxideNiTotal, marker='o', color='b', label = 'S/O Interface')
ax2.plot(TotalLoopDistance, SolutionOxideNiSat, marker='o', color='w', label = 'Ni-Ferrite Solubility')
ax2.plot(TotalLoopDistance, BulkNiTotal, marker='o', color='0.50', label = 'Bulk Coolant') 
ax2.set_xlabel('Distance (m)')
ax2.set_ylabel('$Log_{10}$[Ni Concentration (mol/kg)]')
#ax2.legend()

ax3 = fig1.add_subplot(223)
ax3.plot(TotalLoopDistance, SolutionOxideCoTotal, marker='^', color='m', label = 'S/O Interface')
ax3.plot(TotalLoopDistance, SolutionOxideCoSat, marker='^', color='w', label = 'Co-Ferrite Solubility')
ax3.plot(TotalLoopDistance, BulkCoTotal, marker='^', color='0.50', label = 'Bulk Coolant') 
ax3.set_xlabel('Distance (m)')
ax3.set_ylabel('$Log_{10}$[Co Concentration (mol/kg)]')
#ax3.legend()


ax4 = fig1.add_subplot(224)
ax4.plot(TotalLoopDistance, SolutionOxideCrSat, marker='*', color='0', label = 'Chromite Solubility')
ax4.plot(TotalLoopDistance, BulkCrTotal, marker='*', color='0.50', label = 'Bulk Coolant') 
ax4.set_xlabel('Distance (m)')
ax4.set_ylabel('$Log_{10}$[Cr Concentration (mol/kg)]')

plt.tight_layout()
plt.show()

print (CorrosionRate)

 
fig2, ax1 = plt.subplots()
#ax1 = fig2.add_subplot(221)
ax1.plot(TotalLoopDistance, InnerOxideThickness, linestyle =None, marker='o', color='0.50', label='Inner Oxide')
ax1.plot(TotalLoopDistance, OuterOxideThickness, linestyle =None, marker='o', color='k', label='Outer Oxide')
#ax1.axis([51,69, 0, 30])
ax1.set_xlabel('Distance (m)')
ax1.set_ylabel('Oxide Layer Loadings (${g/m^2}$)')
           
ax2 = ax1.twinx()
ax2.plot(TotalLoopDistance, NickelThickness, linestyle =None, marker='o', color='c', label='Nickel')
ax2.plot(TotalLoopDistance, CobaltThickness, linestyle =None, marker='o', color='m', label='Cobalt')
ax2.set_ylabel('Ni, Co, and Cr Loadings (${g/m^2}$)', rotation=270, labelpad=20)
lines, labels = ax1.get_legend_handles_labels()
lines2, labels2 = ax2.get_legend_handles_labels()
ax2.legend(lines + lines2, labels + labels2, loc=0)        
#plt.axis([51, 69, 0, 30])
plt.tight_layout()
plt.show()

# fig3, ax1 = plt.subplots()
# #ax1 = fig2.add_subplot(221)
# ax1.plot(ActivityTime, Co60SGActivity, linestyle ="-", marker='o', color='0.50', label='Co-60')
# ax1.plot(ActivityTime, Co58SGActivity, linestyle ="-", marker='o', color='b', label='Co-58')
# ax1.plot(ActivityTime, Fe59SGActivity, linestyle ="-", marker='o', color='r', label='Fe-59')
# ax1.plot(ActivityTime, Fe55SGActivity, linestyle ="-", marker='o', color='g', label='Fe-55')
# ax1.plot(ActivityTime, Mn54SGActivity, linestyle ="-", marker='o', color='m', label='Mn-54')
# ax1.plot(ActivityTime, Cr51SGActivity, linestyle ="-", marker='o', color='y', label='Cr-51')
# ax1.plot(ActivityTime, Ni63SGActivity, linestyle ="-", marker='o', color='k', label='Ni-63')
#  
# #ax1.axis([51,69, 0, 30])
# ax1.set_xlabel('EFPY')
# ax1.set_ylabel('Surface Activity (${mCi/m^2}$)')
#             
# plt.tight_layout()
# plt.show()






# fig3, ax1 = plt.subplots()
# #ax3 = fig2.add_subplot(222)
# ax1.plot(TotalLoopDistance, SpallTime, linestyle = '-', marker = 'o', color = 'b', label='Spall Time')
# ax1.plot(TotalLoopDistance, ElapsedTime, linestyle = '--', marker = 'o', color = 'r', label = 'Elapsed Time')
# ax1.set_xlabel('Distance (m)')
# ax1.set_ylabel('Spalling and Elapsed Time (h)')
# ax1.legend()
# plt.tight_layout()
# plt.show()


# plt.plot(TotalLoopDistance, CorrosionRate, 'ro')
# #plt.axis([4, 69, -2, 3]) 
# plt.show()








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

