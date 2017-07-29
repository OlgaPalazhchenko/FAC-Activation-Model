import numpy as np
import csv 
import NumericConstants as nc

SizingParameters = open('SizingParameters.txt','r')      
SizingParametersReader = list(csv.reader(SizingParameters, delimiter = ',')) #Assign file data to list, Reader[row][column], delimiter = comma 
#CoolantParameters = open('CoolantParameters.txt','r')
#CoolantParametersReader = list(csv.reader(CoolantParameters, delimiter = ','))

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
    def __init__(self,j,RowStart,RowEnd): 
        self.RowStart = RowStart
        self.RowEnd = RowEnd
        self.j =j
        
        self.MetalOxide = Interface()
        self.SolutionOxide = Interface()
        self.Bulk = Interface()
        
        #General section parameter input (interface specification not necessary)
        self.Diameter = [float(SizingParametersReader[j][i]) for i in range(self.RowStart,self.RowEnd)] #Reader([row][column]) #[cm]
        self.Velocity = [float(SizingParametersReader[j+1][i]) for i in range(self.RowStart,self.RowEnd)] #[cm/s]
        self.Length = [float(SizingParametersReader[j+2][i]) for i in range(self.RowStart,self.RowEnd)] #[cm]
        self.Distance = np.cumsum(self.Length) #[cm]
        self.Celsius = [float(SizingParametersReader[j+3][i]) for i in range(self.RowStart,self.RowEnd)] #[oC]
        self.Kelvin = [x + 273.15 for x in self.Celsius] # [K]
        self.NernstConstant = [x*(2.303*nc.R/(2*nc.F)) for x in self.Kelvin] #2.303RT/nF
        
        self.DensityH2O = [float(SizingParametersReader[j+4][i]) for i in range(self.RowStart,self.RowEnd)] #[g/cm^3]
        self.ViscosityH2O = [float(SizingParametersReader[j+5][i]) for i in range(self.RowStart,self.RowEnd)] #[g/cm s]
        self.SolubilityFe = [float(SizingParametersReader[j+6][i]) for i in range(self.RowStart,self.RowEnd)] #[mol/kg] (Tremaine & LeBlanc)
        self.SolubilityCr = [float(SizingParametersReader[j+7][i]) for i in range(self.RowStart,self.RowEnd)] #[mol/kg] 
        self.SolubilityCo = [float(SizingParametersReader[j+8][i]) for i in range(self.RowStart,self.RowEnd)] #[mol/kg] (Sandler & Kunig) 
        self.SolubilityNi = [float(SizingParametersReader[j+9][i]) for i in range(self.RowStart,self.RowEnd)] #[mol/kg] (Sandler & Kunig)
        
        self.InnerOxThickness = [float(SizingParametersReader[j+10][i]) for i in range(self.RowStart,self.RowEnd)] #[g/cm^2] 
        self.InnerIronOxThickness = self.InnerOxThickness
        self.OuterFe3O4Thickness = [float(SizingParametersReader[j+11][i]) for i in range(self.RowStart,self.RowEnd)] #[g/cm^2]
        self.NiThickness = [float(SizingParametersReader[j+12][i]) for i in range(self.RowStart,self.RowEnd)] #[g/cm^2]
        self.OuterOxThickness = self.OuterFe3O4Thickness #[g/cm^2]
        self.CoThickness = [i*0 for i in self.OuterFe3O4Thickness]
        self.OxThickness = [x+y for x,y in zip(self.OuterFe3O4Thickness,self.InnerOxThickness)] #adds the lists together element by element
        
        self.StandardEqmPotentialFe = [float(SizingParametersReader[j+13][i]) for i in range(self.RowStart,self.RowEnd)]
        self.StandardEqmPotentialFe3O4red = [float(SizingParametersReader[j+14][i]) for i in range(self.RowStart,self.RowEnd)]
        self.StandardEqmPotentialFe3O4oxid = [float(SizingParametersReader[j+15][i]) for i in range(self.RowStart,self.RowEnd)]
        self.StandardEqmPotentialH2 = [float(SizingParametersReader[j+16][i]) for i in range(self.RowStart,self.RowEnd)]
        self.StandardEqmPotentialNi = [float(SizingParametersReader[j+17][i]) for i in range(self.RowStart,self.RowEnd)]
        self.NodeNumber = int(SizingParametersReader[j+18][0])      
        
        #Equilibrium and Debye-Huckel constants - polynomials as a function of temperature 
        #Coeff1*x^4 + Coeff2*x^3 + Coeff3*x^2 + Coeff4*x + Coeff5, depends on # elements in coeff list
        self.DebyeHuckelConstant=(np.polyval(nc.DebyeHuckPolynomial,self.Celsius)) 
        self.k_W=10**(np.polyval(nc.KwPolynomial,self.Kelvin))
        self.k_Li=10**(np.polyval(nc.KLiPolynomial,self.Kelvin)) 
        self.k_FeOH=10**(np.polyval(nc.KFeOHPolynomial,self.Kelvin))
        self.k_FeOH2=10**(np.polyval(nc.KFeOH2Polynomial,self.Kelvin))
        self.k_FeOH3=10**(np.polyval(nc.KFeOH3Polynomial,self.Kelvin))
        self.k_NiOH=10**(np.polyval(nc.KNiOHPolynomial,self.Kelvin))
        self.k_NiOH2=10**(np.polyval(nc.KNiOH2Polynomial,self.Kelvin))
        self.k_NiOH3=10**(np.polyval(nc.KNiOH3Polynomial,self.Kelvin))
        
        self.km = None
        self.KpFe3O4electrochem = None
        self.KdFe3O4electrochem = None
        self.CorrRate = None
        self.Co60Source = None
        self.Co58Source = None
        self.Mn54Source = None
        self.Fe55Source = None
        self.Fe59Source = None
        self.Cr51Source = None
        self.Ni63Source = None
        
        self.DepositThickness = None
        self.Particle = []
        self.FractionFeInnerOxide = None
        self.FractionNiInnerOxide = None
        
        self.BigParticulate =[0.001]*self.NodeNumber #[mg/kg] (ppm)
        self.SmallParticulate = [0.001]*self.NodeNumber #[mg/kg] (ppm)
        self.ElapsedTime = None
        self.SpallTime = []
        
        #[mol/kg] initial estimate for iron saturation based on Fe3O4 solubility (Tremaine & LeBlanc)
        self.Bulk.FeTotal = [float(SizingParametersReader[j+6][i]) for i in range(self.RowStart,self.RowEnd)]
        self.SolutionOxide.FeTotal = [float(SizingParametersReader[j+6][i]) for i in range(self.RowStart,self.RowEnd)]
        self.MetalOxide.FeTotal = [float(SizingParametersReader[j+6][i]) for i in range(self.RowStart,self.RowEnd)]   
        self.Bulk.FeSatFe3O4 = [float(SizingParametersReader[j+6][i]) for i in range(self.RowStart,self.RowEnd)]
        self.SolutionOxide.FeSatFe3O4 = [float(SizingParametersReader[j+6][i]) for i in range(self.RowStart,self.RowEnd)]
        self.MetalOxide.FeSatFe3O4 = [float(SizingParametersReader[j+6][i]) for i in range(self.RowStart,self.RowEnd)]
        
        #[mol/kg] (Sandler & Kunig)
        self.Bulk.NiTotal = [float(SizingParametersReader[j+9][i]) for i in range(self.RowStart,self.RowEnd)]
        self.SolutionOxide.NiTotal = [float(SizingParametersReader[j+9][i]) for i in range(self.RowStart,self.RowEnd)]
        self.MetalOxide.NiTotal = [float(SizingParametersReader[j+9][i]) for i in range(self.RowStart,self.RowEnd)]
        self.Bulk.NiSatFerrite = [float(SizingParametersReader[j+9][i]) for i in range(self.RowStart,self.RowEnd)]
        self.SolutionOxide.NiSatFerrite = [float(SizingParametersReader[j+9][i]) for i in range(self.RowStart,self.RowEnd)]
        self.MetalOxide.NiSatFerrite = [float(SizingParametersReader[j+9][i]) for i in range(self.RowStart,self.RowEnd)] 
        
        #[mol/kg] (Sandler & Kunig)
        self.Bulk.CoTotal = [float(SizingParametersReader[j+8][i]) for i in range(self.RowStart,self.RowEnd)]
        self.SolutionOxide.CoTotal = [float(SizingParametersReader[j+8][i]) for i in range(self.RowStart,self.RowEnd)] 
        self.MetalOxide.CoTotal = [float(SizingParametersReader[j+8][i]) for i in range(self.RowStart,self.RowEnd)] 
        self.Bulk.CoSatFerrite = [float(SizingParametersReader[j+8][i]) for i in range(self.RowStart,self.RowEnd)] 
        self.SolutionOxide.CoSatFerrite = [float(SizingParametersReader[j+8][i]) for i in range(self.RowStart,self.RowEnd)] 
        self.MetalOxide.CoSatFerrite = [float(SizingParametersReader[j+8][i]) for i in range(self.RowStart,self.RowEnd)]  
        
        self.Bulk.CrTotal = [float(SizingParametersReader[j+7][i]) for i in range(self.RowStart,self.RowEnd)]
        self.SolutionOxide.CrTotal = [float(SizingParametersReader[j+7][i]) for i in range(self.RowStart,self.RowEnd)]
        self.MetalOxide.CrTotal = [float(SizingParametersReader[j+7][i]) for i in range(self.RowStart,self.RowEnd)]
        self.Bulk.CrSatFerrite = [float(SizingParametersReader[j+7][i]) for i in range(self.RowStart,self.RowEnd)]
        self.SolutionOxide.CrSatFerrite = [float(SizingParametersReader[j+7][i]) for i in range(self.RowStart,self.RowEnd)]
        self.MetalOxide.CrSatFerrite = [float(SizingParametersReader[j+7][i]) for i in range(self.RowStart,self.RowEnd)] 
        
#Creates the 4 PHTS sections and their methods based on the Sections class template/blueprint
Inlet= Section(21,0,7)
Core=Section(42,0,12)
Outlet=Section(63,0,9)
SteamGenerator=Section(84,0,22)        


def MassTransfer(Section):
    Schmidt = [x/(y*nc.FeDiffusivity) for x,y in zip(Section.ViscosityH2O, Section.DensityH2O)]
    Reynolds = [x*y/(z*q) for x,y,z,q in zip(Section.Velocity, Section.Diameter, Section.ViscosityH2O, Section.DensityH2O)]
    Sherwood = [0.0165*(x**0.86)*(y**0.33) for x,y in zip(Reynolds, Schmidt)] #Berger & Hau for straight pipe, single phase, fully developed (turbulent) flow
    return [nc.FeDiffusivity*x/y for x,y in zip(Sherwood,Section.Diameter)] 


def UnitConverter(Section, UnitInput, UnitOutput, Concentration, Rate, Oxide, OxideDensity, MolarMass):
    #Converts surface area from /cm^2 to /m^2
    if UnitInput == "Grams per Cm Squared" and UnitOutput == "Grams per M Squared":
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


