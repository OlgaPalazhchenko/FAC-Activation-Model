import numpy as np
import csv 
import NumericConstants as nc


SizingParameters = open('SizingParameters.txt','r')      
SizingParametersReader = list(csv.reader(SizingParameters, delimiter = ',')) #Assign file data to list, Reader[row][column], delimiter = comma 


def UnitConverter(Section, UnitInput, UnitOutput, Concentration, Rate, Oxide, OxideDensity, MolarMass, Temperature):
    #Converts temperature from Celsius to Kelvin
    if UnitInput == "Celsius" and UnitOutput == "Kelvin":
        return [i+273.15 for i in Temperature]
    
    #Converts surface area from /cm^2 to /m^2
    elif UnitInput == "Grams per Cm Squared" and UnitOutput == "Grams per M Squared":
        return [i*(100**2) for i in Oxide]
    
    #Converts coolant activity from Bq/cm^3 to microCi/m^3    
    elif UnitInput == "Bq" and UnitOutput == "Ci":
        #return Concentration/(3.7*10**10)
        return [i/(3.7*10**10) for i in Concentration]
    #Converts concentrations between mol/kg to g/cm^3 and vice versa    
    #[mol_solute/kg_coolant] *[1 kg/1000 g coolant]* [g/mol solute] = [g_solute]/[g_coolant]* [g/cm^3 coolant] = [g_solute/cm^3 coolant]
    elif UnitInput == "Mol per Kg" and UnitOutput == "Grams per Cm Cubed":
        return [(x/1000)*(MolarMass*y) for x,y in zip(Concentration, Section.DensityH2O)]
    
    elif UnitInput == "Grams per Cm Cubed" and UnitOutput == "Mol per Kg":
    #[g_solute/cm^3 coolant]/([g/cm^3 coolant]*[g/mol solute]) = [mol_solute/g coolant]*[1000 g/1 kg coolant] = [mol solute/kg coolant]
        return [x*1000/(MolarMass*y) for x,y in zip(Concentration, Section.DensityH2O)]
        
    #Converts corrosion rates from [g/cm^2 s] to micrometers per annum [um/yr]    
    elif UnitInput =="Corrosion Rate Grams" and UnitOutput == "Corrosion Rate Micrometers":
    #([g/cm^2 s] / [g/cm^3] *[10000 um/cm]) = ([cm/s] * [3600 s/h] * [24 h/d] * [365 d/yr]) = [um/yr]
        if Section == Inlet or Section == Outlet:
            return [i*((1/nc.FeDensity)*3600*24*365*10000) for i in Rate]   
        elif Section ==SteamGenerator:
            return [x*y for x,y in zip(Rate, [(1/nc.Alloy800Density)*3600*24*365*10000]*Section.NodeNumber)]
        else: 
            return [0]*Section.NodeNumber
     
    #Converts from oxide mass per surface area to oxide thickness [g/cm^2] to [um]
    #assuming uniform distribution on surface   
    elif UnitInput =="Oxide Thickness Grams" and UnitOutput == "Oxide Thickness Micrometers":
        return [x*10000/y for x,y in zip(Oxide, OxideDensity)]

    else: print ("Error: incorrect unit input/output combination")


class Interface():
    def __init__(self):
        self.ConcentrationH2 = None
        self.ConcentrationH = None
        self.ConcentrationOH = None
        self.FeSatFe3O4 = None
        self.FeTotal = None
        self.ConcentrationFe2 = None
        self.ConcentrationFeOH = None
        self.ConcentrationFeOH2 = None
        self.ConcentrationFeOH3 = None
        
        self.NiSatFerrite = None
        self.NiSatMetallicNi = None
        self.NiTotal = None
        self.ConcentrationNi2 = None
        self.ConcentrationNiOH = None
        self.ConcentrationNiOH2 = None
        self.ConcentrationNiOH3 = None
        self.CrSat = None
        self.CoSatFerrite = None
        self.CrTotal = None
        self.CoTotal = None
        self.ActivityCoefficient1 =None
        self.ActivityCoefficient2 =None
        
        self.Co60 = None
        self.Co58 = None
        self.Fe55 = None
        self.Fe59 = None
        self.Mn54 = None
        self.Cr51= None
        self.Ni63 = None
        
        self.MixedPotential = None
        self.EqmPotentialFe3O4 = None

class Section(): #Defining each primary heat transport section as a class 
    def __init__(self, j, RowStart, RowEnd): 
        self.RowStart = RowStart
        self.RowEnd = RowEnd
        self.j =j
        
        self.MetalOxide = Interface()
        self.SolutionOxide = Interface()
        self.Bulk = Interface()
        
        #General section parameter input (interface specification not necessary)
        self.Diameter = None#[float(SizingParametersReader[j][i]) for i in range(self.RowStart,self.RowEnd)] #Reader([row][column]) #[cm]
        self.OuterDiameter = None
        self.TubeThickness = None
        self.Velocity = None#[float(SizingParametersReader[j+1][i]) for i in range(self.RowStart,self.RowEnd)] #[cm/s]
        self.Length = None #[float(SizingParametersReader[j+2][i]) for i in range(self.RowStart,self.RowEnd)] #[cm]
        self.Distance = None #np.cumsum(self.Length) #[cm]
        self.PrimaryBulkTemperature = None#[x + 273.15 for x in self.Celsius] # [K]
        self.PrimaryWallTemperature = None
        self.SecondaryWallTemperature = None
        self.SecondaryBulkTemperature = None
        self.NernstConstant = None#[x*(2.303*nc.R/(2*nc.F)) for x in self.Kelvin] #2.303RT/nF
        
        self.DensityH2O = [float(SizingParametersReader[j+4][i]) for i in range(self.RowStart,self.RowEnd)] #[g/cm^3]
        self.ViscosityH2O = [float(SizingParametersReader[j+5][i]) for i in range(self.RowStart,self.RowEnd)] #[g/cm s]
        self.SolubilityFe = None#[float(SizingParametersReader[j+6][i]) for i in range(self.RowStart,self.RowEnd)] #[mol/kg] (Tremaine & LeBlanc)
        self.SolubilityCr = None#[float(SizingParametersReader[j+7][i]) for i in range(self.RowStart,self.RowEnd)] #[mol/kg] 
        self.SolubilityCo = None#[float(SizingParametersReader[j+8][i]) for i in range(self.RowStart,self.RowEnd)] #[mol/kg] (Sandler & Kunig) 
        self.SolubilityNi = None#[float(SizingParametersReader[j+9][i]) for i in range(self.RowStart,self.RowEnd)] #[mol/kg] (Sandler & Kunig)
        
        self.InnerOxThickness = None#[float(SizingParametersReader[j+10][i]) for i in range(self.RowStart,self.RowEnd)] #[g/cm^2] 
        self.InnerIronOxThickness = None#self.InnerOxThickness
        self.OuterFe3O4Thickness = None#[float(SizingParametersReader[j+11][i]) for i in range(self.RowStart,self.RowEnd)] #[g/cm^2]
        self.NiThickness = None#[float(SizingParametersReader[j+12][i]) for i in range(self.RowStart,self.RowEnd)] #[g/cm^2]
        self.OuterOxThickness = None#self.OuterFe3O4Thickness #[g/cm^2]
        self.CoThickness = None#[i*0 for i in self.OuterFe3O4Thickness]
        self.OxThickness = None#[x+y for x,y in zip(self.OuterFe3O4Thickness,self.InnerOxThickness)] #adds the lists together element by element
        
        self.StandardEqmPotentialFe = [float(SizingParametersReader[j+13][i]) for i in range(self.RowStart,self.RowEnd)]
        self.StandardEqmPotentialFe3O4red = [float(SizingParametersReader[j+14][i]) for i in range(self.RowStart,self.RowEnd)]
        self.StandardEqmPotentialFe3O4oxid = [float(SizingParametersReader[j+15][i]) for i in range(self.RowStart,self.RowEnd)]
        self.StandardEqmPotentialH2 = [float(SizingParametersReader[j+16][i]) for i in range(self.RowStart,self.RowEnd)]
        self.StandardEqmPotentialNi = [float(SizingParametersReader[j+17][i]) for i in range(self.RowStart,self.RowEnd)]
        self.NodeNumber = None#int(SizingParametersReader[j+18][0])      
        
        self.km = None
        self.KpFe3O4electrochem = None
        self.KdFe3O4electrochem = None
        self.CorrRate = None
        
        self.DepositThickness = None
        self.Particle = []
        self.FractionFeInnerOxide = None
        self.FractionNiInnerOxide = None
        
        self.ElapsedTime = None
        self.SpallTime = []
        
        #[mol/kg] initial estimate for iron saturation based on Fe3O4 solubility (Tremaine & LeBlanc)
        self.Bulk.FeTotal = None#[float(SizingParametersReader[j+6][i]) for i in range(self.RowStart,self.RowEnd)]
        self.SolutionOxide.FeTotal = None#[float(SizingParametersReader[j+6][i]) for i in range(self.RowStart,self.RowEnd)]
        self.MetalOxide.FeTotal = None#[float(SizingParametersReader[j+6][i]) for i in range(self.RowStart,self.RowEnd)]   
        self.Bulk.FeSatFe3O4 = None#[float(SizingParametersReader[j+6][i]) for i in range(self.RowStart,self.RowEnd)]
        self.SolutionOxide.FeSatFe3O4 = None#[float(SizingParametersReader[j+6][i]) for i in range(self.RowStart,self.RowEnd)]
        self.MetalOxide.FeSatFe3O4 = None#[float(SizingParametersReader[j+6][i]) for i in range(self.RowStart,self.RowEnd)]
        
        #[mol/kg] (Sandler & Kunig)
        self.Bulk.NiTotal =None #[float(SizingParametersReader[j+9][i]) for i in range(self.RowStart,self.RowEnd)]
        self.SolutionOxide.NiTotal = None#[float(SizingParametersReader[j+9][i]) for i in range(self.RowStart,self.RowEnd)]
        self.MetalOxide.NiTotal = None#[float(SizingParametersReader[j+9][i]) for i in range(self.RowStart,self.RowEnd)]
        self.Bulk.NiSatFerrite = None#[float(SizingParametersReader[j+9][i]) for i in range(self.RowStart,self.RowEnd)]
        self.SolutionOxide.NiSatFerrite = None#[float(SizingParametersReader[j+9][i]) for i in range(self.RowStart,self.RowEnd)]
        self.MetalOxide.NiSatFerrite = None#[float(SizingParametersReader[j+9][i]) for i in range(self.RowStart,self.RowEnd)] 
        
        #[mol/kg] (Sandler & Kunig)
        self.Bulk.CoTotal = None#[float(SizingParametersReader[j+8][i]) for i in range(self.RowStart,self.RowEnd)]
        self.SolutionOxide.CoTotal = None#[float(SizingParametersReader[j+8][i]) for i in range(self.RowStart,self.RowEnd)] 
        self.MetalOxide.CoTotal = None#[float(SizingParametersReader[j+8][i]) for i in range(self.RowStart,self.RowEnd)] 
        self.Bulk.CoSatFerrite = None#[float(SizingParametersReader[j+8][i]) for i in range(self.RowStart,self.RowEnd)] 
        self.SolutionOxide.CoSatFerrite = None#[float(SizingParametersReader[j+8][i]) for i in range(self.RowStart,self.RowEnd)] 
        self.MetalOxide.CoSatFerrite = None#[float(SizingParametersReader[j+8][i]) for i in range(self.RowStart,self.RowEnd)]  
        
        self.Bulk.CrTotal = None#[float(SizingParametersReader[j+7][i]) for i in range(self.RowStart,self.RowEnd)]
        self.SolutionOxide.CrTotal = None#[float(SizingParametersReader[j+7][i]) for i in range(self.RowStart,self.RowEnd)]
        self.MetalOxide.CrTotal = None#[float(SizingParametersReader[j+7][i]) for i in range(self.RowStart,self.RowEnd)]
        self.Bulk.CrSat = None#[float(SizingParametersReader[j+7][i]) for i in range(self.RowStart,self.RowEnd)]
        self.SolutionOxide.CrSat = None#[float(SizingParametersReader[j+7][i]) for i in range(self.RowStart,self.RowEnd)]
        self.MetalOxide.CrSat = None#[float(SizingParametersReader[j+7][i]) for i in range(self.RowStart,self.RowEnd)] 
        
#Creates the 4 PHTS sections and their methods based on the Sections class template/blueprint
Inlet= Section(21,0,7)
Core=Section(42,0,12)
Outlet=Section(63,0,9)
SteamGenerator= Section(84,0,22)

## Steam generator split into multiple sections/zones
# Temporary 2 "zones" added 
SG_Zone1= Section(84,0,22)
SG_Zone2 =Section(84,0,22)
##

#Node Number 
Inlet.NodeNumber = 7
Core.NodeNumber = 12
Outlet.NodeNumber = 9
SteamGenerator.NodeNumber = 22
SG_Zone1.NodeNumber = 22
SG_Zone2.NodeNumber = 22 


#Diameter [cm]
Inlet.Diameter = [44.3, 50, 106, 5.68, 5.68, 5.68, 5.68] 
Core.Diameter = [1.3]*Core.NodeNumber 
Outlet.Diameter = [6.4, 6.4, 6.4, 6.4, 8.9, 8.9, 8.9, 116, 40.8] 
SteamGenerator.Diameter = [1.59]*SteamGenerator.NodeNumber 
SG_Zone1.Diameter = [1.59]*SG_Zone1.NodeNumber
SG_Zone2.Diameter = [1.59]*SG_Zone2.NodeNumber


#Velocity [cm/s]
Inlet.Velocity = [1530, 1200, 270, 985, 985, 985, 985]
Core.Velocity =[883.08, 890.66, 900.3, 910.64, 920.97, 932.68, 945.08, 958.17, 973.32, 989.16, 1073.89, 1250.92]
Outlet.Velocity =[1619, 1619, 1619, 1619, 857, 857, 857, 306, 1250]
SteamGenerator.Velocity =[533.002, 533.001, 533, 531, 524, 517, 511, 506, 506, 502, 498, 494, 491, 489, 487, 484, 483, 481, 480, 479, 476, 474]
SG_Zone1.Velocity =[533.002, 533.001, 533, 531, 524, 517, 511, 506, 506, 502, 498, 494, 491, 489, 487, 484, 483, 481, 480, 479, 476, 474]
SG_Zone2.Velocity =[533.002, 533.001, 533, 531, 524, 517, 511, 506, 506, 502, 498, 494, 491, 489, 487, 484, 483, 481, 480, 479, 476, 474]


#Length [cm]
Inlet.Length = [477.6, 281.8, 78.6, 350, 350, 350, 350] 
Core.Length = [49.5]*Core.NodeNumber 
Outlet.Length = [17, 3.5, 139.5, 432, 225.5, 460.3, 460.3, 400, 100]
SteamGenerator.Length = [60.96, 88.9, 88.9, 88.9, 91.44, 91.44, 91.44, 91.44, 91.44, 103.05, 103.5, 91.44, 91.44, 91.44, 91.44, 91.44, 76.2, 58.42, 58.42, 58.42, 45.72, 30.48]
SG_Zone1.Length = [60.96, 88.9, 88.9, 88.9, 91.44, 91.44, 91.44, 91.44, 91.44, 103.05, 103.5, 91.44, 91.44, 91.44, 91.44, 91.44, 76.2, 58.42, 58.42, 58.42, 45.72, 30.48]
SG_Zone2.Length = [60.96, 88.9, 88.9, 88.9, 91.44, 91.44, 91.44, 91.44, 91.44, 103.05, 103.5, 91.44, 91.44, 91.44, 91.44, 91.44, 76.2, 58.42, 58.42, 58.42, 45.72, 30.48]


#Solubility (mol/kg)
Inlet.SolubilityFe = [1.28494E-08]*Inlet.NodeNumber
Inlet.SolubilityNi = [2.71098E-09]*Inlet.NodeNumber
Inlet.SolubilityCo = [3.77742E-09]*Inlet.NodeNumber
Inlet.SolubilityCr = [8.81E-11]*Inlet.NodeNumber

Core.SolubilityFe = [1.28521E-08, 1.35852E-08, 1.449E-08, 1.54787E-08, 1.65654E-08, 1.77722E-08, 1.91309E-08, 2.06833E-08, 2.25034E-08, 2.46889E-08, 2.64104E-08, 2.64104E-08]
Core.SolubilityNi = [2.71098E-09, 2.66196E-09, 2.60372E-09, 2.54535E-09, 2.29445E-09, 2.2137E-09, 1.91294E-09, 1.82788E-09, 1.54081E-09, 1.54384E-09, 1.54584E-09, 1.54584E-09]
Core.SolubilityCo = [3.77742E-09, 3.51525E-09, 3.20371E-09, 2.8915E-09, 2.60041E-09, 2.3011E-09, 2.08369E-09, 1.86963E-09, 1.67578E-09, 1.53187E-09, 1.43654E-09, 1.43654E-09]
Core.SolubilityCr = [8.81E-11, 9.61E-11, 1.01E-10, 9.40E-11, 8.69E-11, 7.98E-11, 7.28E-11, 6.57E-11, 5.86E-11, 5.16E-11, 4.69E-11, 4.69E-11]

Outlet.SolubilityFe = [2.58268E-08]*Outlet.NodeNumber
Outlet.SolubilityNi = [1.54584E-09]*Outlet.NodeNumber
Outlet.SolubilityCo = [1.44E-09]*Outlet.NodeNumber
Outlet.SolubilityCr = [4.84E-11]*Outlet.NodeNumber

#Replicate solubilities and temperatures here for now 
SGInterfaces = [SteamGenerator, SG_Zone1, SG_Zone2]
for SGInterface in SGInterfaces:
    SGInterface.SolubilityFe = [2.58E-08, 2.56E-08, 2.55E-08, 2.52546E-08, 2.32347E-08, 2.16094E-08, 2.03107E-08, 1.9246E-08, 1.83579E-08, 1.75642E-08, 1.68473E-08, 1.62692E-08, 1.58017E-08, 1.53906E-08, 1.50304E-08, 1.47083E-08, 1.44316E-08, 1.42108E-08, 1.39481E-08, 1.35926E-08, 1.31784E-08, 1.27883E-08]
    SGInterface.SolubilityNi = [1.5452E-09, 1.5452E-09, 1.5452E-09, 1.54453E-09, 1.54189E-09, 1.78271E-09, 1.84719E-09, 1.9062E-09, 1.96011E-09, 2.22698E-09, 2.27478E-09, 2.31567E-09, 2.35035E-09, 2.3821E-09, 2.41091E-09, 2.59037E-09, 2.60733E-09, 2.62118E-09, 2.63802E-09, 2.66147E-09, 2.68978E-09, 2.71747E-09]
    SGInterface.SolubilityCo = [1.46729E-09, 1.46729E-09, 1.46729E-09, 1.49896E-09, 1.62443E-09, 1.75595E-09, 1.91821E-09, 2.06673E-09, 2.2024E-09, 2.35035E-09, 2.5275E-09, 2.67907E-09, 2.80762E-09, 2.92529E-09, 3.0321E-09, 3.13232E-09, 3.22305E-09, 3.2971E-09, 3.38716E-09, 3.51258E-09, 3.66401E-09, 3.81211E-09]
    SGInterface.SolubilityCr = [4.84E-11, 4.84E-11, 4.84E-11, 4.99E-11, 5.61E-11, 6.20E-11, 6.73E-11, 7.22E-11, 7.67E-11, 8.10E-11, 8.52E-11, 8.88E-11, 9.18E-11, 9.46E-11, 9.71E-11, 9.94E-11, 1.01E-10, 1.03E-10, 1.00E-10, 9.62E-11, 9.15E-11, 8.70E-11]
    SGInterface.PrimaryBulkTemperature = \
    UnitConverter(SGInterface, "Celsius", "Kelvin", None, None, None, None, None, \
                  [310.002, 310.001, 310, 308.97, 304.89, 301.02, 297.48, 294.24, 291.28, 288.42, 285.65, 283.28, 281.27, 279.43, 277.76, 276.22, 274.86, 273.75, 272.4, 270.52, 268.25, 266.03])
    
#Temperature [Celsius]
Inlet.PrimaryBulkTemperature = UnitConverter(Inlet, "Celsius", "Kelvin", None, None, None, None, None, [266]*Inlet.NodeNumber)
Core.PrimaryBulkTemperature = \
UnitConverter(Inlet, "Celsius", "Kelvin", None, None, None, None, None, [266.55, 270.48, 275.15, 279.83, 284.51, 289.19, 293.87, 298.54, 303.22, 307.9, 311, 311])
Outlet.PrimaryBulkTemperature = UnitConverter(Inlet, "Celsius", "Kelvin", None, None, None, None, None, [310]*Outlet.NodeNumber)


# InletInterfaces = [Inlet.Bulk, Inlet.SolutionOxide, Inlet.MetalOxide]
# for i in InletInterfaces:
#     i.FeTotal = [1.28494E-08]*Inlet.NodeNumber
#     i.FeSatFe3O4 = [1.28494E-08]*Inlet.NodeNumber
#     i.NiTotal = [2.71098E-09]*Inlet.NodeNumber
#     i.NiSatFerrite = [2.71098E-09]*Inlet.NodeNumber
#     i.CoTotal = [3.77742E-09]*Inlet.NodeNumber
#     i.CoSatFerrite = [3.77742E-09]*Inlet.NodeNumber
#     i.CrTotal = [8.81E-11]*Inlet.NodeNumber
#     i.CrSat= [8.81E-11]*Inlet.NodeNumber

# CoreInterfaces = [Core.Bulk, Core.SolutionOxide, Core.MetalOxide]
# for i in CoreInterfaces:
#     i.FeTotal = [1.28521E-08, 1.35852E-08, 1.449E-08, 1.54787E-08, 1.65654E-08, 1.77722E-08, 1.91309E-08, 2.06833E-08, 2.25034E-08, 2.46889E-08, 2.64104E-08, 2.64104E-08]
#     i.FeSatFe3O4 = [i*1 for i in i.FeTotal]#[1.28521E-08, 1.35852E-08, 1.449E-08, 1.54787E-08, 1.65654E-08, 1.77722E-08, 1.91309E-08, 2.06833E-08, 2.25034E-08, 2.46889E-08, 2.64104E-08, 2.64104E-08]
#     i.NiTotal = [2.71098E-09, 2.66196E-09, 2.60372E-09, 2.54535E-09, 2.29445E-09, 2.2137E-09, 1.91294E-09, 1.82788E-09, 1.54081E-09, 1.54384E-09, 1.54584E-09, 1.54584E-09]
#     i.NiSatFerrite = [i*1 for i in i.NiTotal]
#     #[2.71098E-09, 2.66196E-09, 2.60372E-09, 2.54535E-09, 2.29445E-09, 2.2137E-09, 1.91294E-09, 1.82788E-09, 1.54081E-09, 1.54384E-09, 1.54584E-09, 1.54584E-09]
#     i.CoTotal = [3.77742E-09, 3.51525E-09, 3.20371E-09, 2.8915E-09, 2.60041E-09, 2.3011E-09, 2.08369E-09, 1.86963E-09, 1.67578E-09, 1.53187E-09, 1.43654E-09, 1.43654E-09]
#     i.CoSatFerrite = [3.77742E-09, 3.51525E-09, 3.20371E-09, 2.8915E-09, 2.60041E-09, 2.3011E-09, 2.08369E-09, 1.86963E-09, 1.67578E-09, 1.53187E-09, 1.43654E-09, 1.43654E-09]
#     i.CrTotal = [8.81E-11, 9.61E-11, 1.01E-10, 9.40E-11, 8.69E-11, 7.98E-11, 7.28E-11, 6.57E-11, 5.86E-11, 5.16E-11, 4.69E-11, 4.69E-11] 
#     i.CrSat = [8.81E-11, 9.61E-11, 1.01E-10, 9.40E-11, 8.69E-11, 7.98E-11, 7.28E-11, 6.57E-11, 5.86E-11, 5.16E-11, 4.69E-11, 4.69E-11]
# 
# OutletInterfaces = [Outlet.Bulk, Outlet.SolutionOxide, Outlet.MetalOxide]
# for i in OutletInterfaces:
#     i.FeTotal = [2.58268E-08]*Outlet.NodeNumber
#     i.FeSatFe3O4 = [2.58268E-08]*Outlet.NodeNumber
#     i.NiTotal = [1.54584E-09]*Outlet.NodeNumber
#     i.NiSatFerrite = [1.54584E-09]*Outlet.NodeNumber
#     i.CoTotal = [1.44E-09]*Outlet.NodeNumber
#     i.CoSatFerrite = [1.44E-09]*Outlet.NodeNumber
#     i.CrTotal = [4.84E-11]*Outlet.NodeNumber
#     i.CrSat = [4.84E-11]*Outlet.NodeNumber
# 
# SGInterfaces = [SteamGenerator.Bulk, SteamGenerator.SolutionOxide, SteamGenerator.MetalOxide]
# for i in SGInterfaces:
#     i.FeTotal = [2.58E-08, 2.56E-08, 2.55E-08, 2.52546E-08, 2.32347E-08, 2.16094E-08, 2.03107E-08, 1.9246E-08, 1.83579E-08, 1.75642E-08, 1.68473E-08, 1.62692E-08, 1.58017E-08, 1.53906E-08, 1.50304E-08, 1.47083E-08, 1.44316E-08, 1.42108E-08, 1.39481E-08, 1.35926E-08, 1.31784E-08, 1.27883E-08]
#     i.FeSatFe3O4 = [2.58E-08, 2.56E-08, 2.55E-08, 2.52546E-08, 2.32347E-08, 2.16094E-08, 2.03107E-08, 1.9246E-08, 1.83579E-08, 1.75642E-08, 1.68473E-08, 1.62692E-08, 1.58017E-08, 1.53906E-08, 1.50304E-08, 1.47083E-08, 1.44316E-08, 1.42108E-08, 1.39481E-08, 1.35926E-08, 1.31784E-08, 1.27883E-08]
#     i.NiTotal = [1.5452E-09, 1.5452E-09, 1.5452E-09, 1.54453E-09, 1.54189E-09, 1.78271E-09, 1.84719E-09, 1.9062E-09, 1.96011E-09, 2.22698E-09, 2.27478E-09, 2.31567E-09, 2.35035E-09, 2.3821E-09, 2.41091E-09, 2.59037E-09, 2.60733E-09, 2.62118E-09, 2.63802E-09, 2.66147E-09, 2.68978E-09, 2.71747E-09] 
#     i.NiSatFerrite = [1.5452E-09, 1.5452E-09, 1.5452E-09, 1.54453E-09, 1.54189E-09, 1.78271E-09, 1.84719E-09, 1.9062E-09, 1.96011E-09, 2.22698E-09, 2.27478E-09, 2.31567E-09, 2.35035E-09, 2.3821E-09, 2.41091E-09, 2.59037E-09, 2.60733E-09, 2.62118E-09, 2.63802E-09, 2.66147E-09, 2.68978E-09, 2.71747E-09]
#     i.CoTotal = [1.46729E-09, 1.46729E-09, 1.46729E-09, 1.49896E-09, 1.62443E-09, 1.75595E-09, 1.91821E-09, 2.06673E-09, 2.2024E-09, 2.35035E-09, 2.5275E-09, 2.67907E-09, 2.80762E-09, 2.92529E-09, 3.0321E-09, 3.13232E-09, 3.22305E-09, 3.2971E-09, 3.38716E-09, 3.51258E-09, 3.66401E-09, 3.81211E-09]
#     i.CoSatFerrite = [1.46729E-09, 1.46729E-09, 1.46729E-09, 1.49896E-09, 1.62443E-09, 1.75595E-09, 1.91821E-09, 2.06673E-09, 2.2024E-09, 2.35035E-09, 2.5275E-09, 2.67907E-09, 2.80762E-09, 2.92529E-09, 3.0321E-09, 3.13232E-09, 3.22305E-09, 3.2971E-09, 3.38716E-09, 3.51258E-09, 3.66401E-09, 3.81211E-09]
#     i.CrTotal = [4.84E-11, 4.84E-11, 4.84E-11, 4.99E-11, 5.61E-11, 6.20E-11, 6.73E-11, 7.22E-11, 7.67E-11, 8.10E-11, 8.52E-11, 8.88E-11, 9.18E-11, 9.46E-11, 9.71E-11, 9.94E-11, 1.01E-10, 1.03E-10, 1.00E-10, 9.62E-11, 9.15E-11, 8.70E-11]
#     i.CrSat = [4.84E-11, 4.84E-11, 4.84E-11, 4.99E-11, 5.61E-11, 6.20E-11, 6.73E-11, 7.22E-11, 7.67E-11, 8.10E-11, 8.52E-11, 8.88E-11, 9.18E-11, 9.46E-11, 9.71E-11, 9.94E-11, 1.01E-10, 1.03E-10, 1.00E-10, 9.62E-11, 9.15E-11, 8.70E-11]



#Section-specific parameters
Sections = [Inlet, Core, Outlet, SteamGenerator, SG_Zone1, SG_Zone2]
for Section in Sections:
    
    Interfaces = [Section.MetalOxide, Section.SolutionOxide, Section.Bulk]
    
    for Interface in Interfaces:
        ##Concentration/Saturation Input [mol/kg]
        Interface.FeTotal = [i*1 for i in Section.SolubilityFe]
        Interface.FeSatFe3O4 = [i*1 for i in Section.SolubilityFe]
        Interface.NiTotal = [i*1 for i in Section.SolubilityNi]
        Interface.NiSatFerrite = [i*1 for i in Section.SolubilityNi]
        Interface.NiSatMetallicNi = [i*1 for i in Section.SolubilityNi]
        
        Interface.CoTotal = [i*1 for i in Section.SolubilityCo]
        Interface.CoSatFerrite = [i*1 for i in Section.SolubilityCo]
        Interface.CrTotal = [i*1 for i in Section.SolubilityCr]
        Interface.CrSat = [i*1 for i in Section.SolubilityCr]
        
        if Section == Inlet or Section == SteamGenerator or Section == SG_Zone1 or Section == SG_Zone2:
            if Interface == Section.SolutionOxide:
                Interface.FeTotal = [i*0.8 for i in Section.SolubilityFe]
        
        if Section == Core and Interface == Section.MetalOxide:
            Interface.FeTotal = [0]*Section.NodeNumber
        
        if Section == Outlet and Interface == Section.MetalOxide:
            Interface.FeTotal = [0.00000026]*Section.NodeNumber #From Cook's thesis - experimental corrosion rate measurements and calcs 
        
        if Section == SteamGenerator or Section == SG_Zone1 or Section == SG_Zone2:
            if Interface == Section.MetalOxide:
                Interface.NiTotal = [0]*Section.NodeNumber
        else:
            Interface.NiTotal = [i*1 for i in Section.SolubilityNi]
        ## 
           
    ##Particulate #[mg/kg] (ppm)
    Section.SmallParticulate = [0]*Section.NodeNumber
    Section.BigParticulate = [0]*Section.NodeNumber
    ##
    
    ##Oxide thicknesses [g/cm^2]
    if Section == SteamGenerator or Section == SG_Zone1 or Section== SG_Zone2:
        Section.OuterFe3O4Thickness = [1.3E-4]*Section.NodeNumber
        Section.NiThickness = [1.3E-4]*Section.NodeNumber
        
        Section.TubeThickness = 0.125
        Section.OuterDiameter = [x+2*Section.TubeThickness for x in Section.Diameter]
        
        
    elif Section == Core:
        Section.OuterFe3O4Thickness = [0]*Section.NodeNumber
    else:
        Section.OuterFe3O4Thickness = [2.54E-4]*Section.NodeNumber
    
    Section.OuterOxThickness = [i*1 for i in Section.OuterFe3O4Thickness]    
    Section.InnerOxThickness = [2.5E-4]*Section.NodeNumber
    Section.InnerIronOxThickness = [i*1 for i in Section.InnerOxThickness]
    Section.NiThickness = [0]*Section.NodeNumber
    Section.CoThickness = [0]*Section.NodeNumber
    Section.OxThickness = [x+y for x,y in zip(Section.InnerOxThickness, Section.OuterOxThickness)]
    ##
    
    ##Distance
    Section.Distance = np.cumsum(Section.Length) 
    ##
    
    ##Temperature-dependent parameters            
    Kelvin = UnitConverter(Section, "Celsius", "Kelvin", None, None, None, None, None, Section.PrimaryBulkTemperature) # [K]
    Section.NernstConstant = [x*(2.303*nc.R/(2*nc.F)) for x in Kelvin]
    #Equilibrium and Debye-Huckel constants - polynomials as a function of temperature 
    #Coeff1*x^4 + Coeff2*x^3 + Coeff3*x^2 + Coeff4*x + Coeff5, depends on # elements in coeff list
    Section.DebyeHuckelConstant=(np.polyval(nc.DebyeHuckPolynomial, Section.PrimaryBulkTemperature)) 
    Section.k_W=10**(np.polyval(nc.KwPolynomial, Kelvin))
    Section.k_Li=10**(np.polyval(nc.KLiPolynomial, Kelvin)) 
    Section.k_FeOH=10**(np.polyval(nc.KFeOHPolynomial, Kelvin))
    Section.k_FeOH2=10**(np.polyval(nc.KFeOH2Polynomial, Kelvin))
    Section.k_FeOH3=10**(np.polyval(nc.KFeOH3Polynomial, Kelvin))
    Section.k_NiOH=10**(np.polyval(nc.KNiOHPolynomial, Kelvin))
    Section.k_NiOH2=10**(np.polyval(nc.KNiOH2Polynomial, Kelvin))
    Section.k_NiOH3=10**(np.polyval(nc.KNiOH3Polynomial, Kelvin))
    ##

def ReynoldsNumber(Section, Diameter):
    #Diameter is an input due to difference in desired dimension (e.g., inner, outer, hydraulic, etc.)
    Reynolds = [x*y/(z*q) for x,y,z,q in zip(Section.Velocity, Diameter, Section.ViscosityH2O, Section.DensityH2O)]    
    return Reynolds


def MassTransfer(Section):
    Schmidt = [x/(y*nc.FeDiffusivity) for x,y in zip(Section.ViscosityH2O, Section.DensityH2O)]
    Reynolds = ReynoldsNumber(Section, Section.Diameter)
    Sherwood = [0.0165*(x**0.86)*(y**0.33) for x,y in zip(Reynolds, Schmidt)] #Berger & Hau for straight pipe, single phase, fully developed (turbulent) flow
    return [nc.FeDiffusivity*x/y for x,y in zip(Sherwood,Section.Diameter)] 

