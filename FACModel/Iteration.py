import LepreauData as ld
import numpy as np
import NumericConstants as nc
import Composition as c
import Electrochemistry as e

def Diffusion(Section, Element):
    #Inner oxide density is determined as a weighted average of the densities of the 2 comprising oxides 
    if Section == ld.SteamGenerator:
        AlloyDensity = nc.Alloy800Density
        OxideDensity = c.FractionChromite(Section)*nc.FeCr2O4Density + (1-c.FractionChromite(Section))*nc.NiFe2O4Density
        FractionCo = nc.FractionCo_Alloy800
    elif Section == ld.Inlet or Section == ld.Outlet:
        AlloyDensity = nc.FeDensity
        OxideDensity = c.FractionChromite(Section)*nc.FeCr2O4Density + (1-c.FractionChromite(Section))*nc.Fe3O4Density
        FractionCo = nc.FractionCo_CS
        
    if Section == ld.Inlet or Section == ld.Outlet or Section == ld.SteamGenerator:
        if Element == "Fe":
            MetalRetainedInsideInnerOxide = c.FractionMetalInnerOxide(Section, "Fe")*(OxideDensity / AlloyDensity) * (1 -nc.Fe3O4Porosity_inner) #[gFe/g CSlost] 
        elif Element == "Ni":
            MetalRetainedInsideInnerOxide = c.FractionMetalInnerOxide(Section, "Ni")*(OxideDensity / AlloyDensity) * (1 -nc.Fe3O4Porosity_inner)
        if Element == "Co":
            MetalRetainedInsideInnerOxide = c.FractionMetalInnerOxide(Section, "Co")*(OxideDensity / AlloyDensity) * (1 -nc.Fe3O4Porosity_inner)
        
    if Element == "Co":
        return FractionCo - MetalRetainedInsideInnerOxide
    elif Element == "Ni" or Element == "Fe":
        return 1-MetalRetainedInsideInnerOxide #[g metal/ g alloy lost]


def MetalOxideInterfaceConcentration(Section, Element, SolutionOxideInterfaceConcentration, InnerOxThickness, OuterOxThickness, Corrosion):
    
    if Section == ld.SteamGenerator:
        OxideDensity =nc.NiFe2O4Density
    else:
        OxideDensity = nc.Fe3O4Density
        
    if Element == "Fe":
        Diffusivity = [nc.FeDiffusivity]*Section.NodeNumber
        MolarMass = nc.FeMolarMass
    elif Element == "Ni":
        Diffusivity = [nc.NiDiffusivity]*Section.NodeNumber
        MolarMass = nc.NiMolarMass
    elif Element == "Co":
        Diffusivity = [nc.CoDiffusivity]*Section.NodeNumber
        MolarMass = nc.CoMolarMass
    
    PathLength = [(x*nc.Fe3O4Tortuosity/(OxideDensity*(1-nc.Fe3O4Porosity_inner))) + (y*nc.Fe3O4Tortuosity/(OxideDensity*(1-nc.Fe3O4Porosity_outer))) \
                  for x,y in zip(InnerOxThickness, OuterOxThickness)] 
    
    
#     for i in range(Section.NodeNumber):
#         if Section==ld.Outlet and InnerOxThickness[i] <= 0:
#             InnerOxThickness[i] = 0.00025 #Resets to original thickness (resetting corrosion rate here to ~75 um/a
    
    #inner ox decrease = path length decreases --> diffusivity increases --> metal oxide concentration decreases --> corr rate incr
    #minimum oxide thickness = 1e-8 g/cm^2 --> corr rate = 2.41e-9 g/cm^2 s
    
    DiffusivityTerm = [(x*nc.Fe3O4Porosity_inner)/y for x,y in zip(PathLength, Diffusivity)]
    
    SolutionConcentration = ld.UnitConverter(Section, "Mol per Kg", "Grams per Cm Cubed", SolutionOxideInterfaceConcentration, None, None, None, MolarMass, None)
    
    MOConcentration = [(x*Diffusion(Section, Element)/ y) +z for x,y,z in zip(Corrosion, DiffusivityTerm, SolutionConcentration)]
    

#     ##Pre-set corrosion rate for outlet (lower MOConc = higher FAC rate (e-chem effect))
    #inner ox decrease = path length decreases --> diffusivity term increases --> metal oxide concentration decreases --> corr rate incr
    # Can also reset the oxide thickness such that the minimum results in acceptable MO Conc and corrrate
    
    #if Section == ld.Outlet: 
        #MOConcentration =[1.22e-8]*Section.NodeNumber #110 um/a rate at this fixed M/O Fe concentration
    
    ##       
    return ld.UnitConverter(Section, "Grams per Cm Cubed", "Mol per Kg", MOConcentration, None, None, None, MolarMass, None)
    

def TransfertoBulk(Section,BulkConcentration,MolarMass,km): #[g/cm^2 s]
    Concentration = ld.UnitConverter(Section, "Mol per Kg", "Grams per Cm Cubed", BulkConcentration, None, None, None, MolarMass, None)
    return [x*y for x,y in zip(km,Concentration)] 
    

def InterfaceOxideKinetics (Section, KineticConstant, SaturationConcentration, MolarMass):    
#print (MetalOxideInterfaceConcentration(ld.Inlet, "Fe", [1e-8]*ld.Inlet.NodeNumber, [1e-3]*ld.Inlet.NodeNumber, [1e-3]*ld.Inlet.NodeNumber, [1e-10]*ld.Inlet.NodeNumber, [1e-9]*ld.Inlet.NodeNumber))
    Concentration = ld.UnitConverter(Section, "Mol per Kg", "Grams per Cm Cubed", SaturationConcentration, None, None, None, MolarMass, None) 
    return [x*y for x,y in zip(KineticConstant, Concentration)] 
    

def SolutionOxide(Section, BulkConcentration, SaturationConcentration, SolutionOxideConcentration, InnerIronOxThickness, \
                  OuterFe3O4Thickness, NiThickness, CoThickness, Element):
    
    km = ld.MassTransfer(Section)
    if Section == ld.Core:
        Diff = [0]*Section.NodeNumber
    else:
        Diff = [i* Diffusion(Section, Element) for i in Section.CorrRate]
        if Section==ld.Inlet or Section==ld.Outlet:
            if Element == "Ni": 
                Diff = [0]*Section.NodeNumber
    
    #Corrosion rate becomes negative --> Diff --> S/O Concentration 
    
    if Element == "Fe": 
        MolarMass = nc.FeMolarMass
    elif Element == "Ni": 
        MolarMass = nc.NiMolarMass
    elif Element == "Co": 
        MolarMass = nc.CoMolarMass
    else:
        print ("Error: unknown element specified")            
     
    BTrans = TransfertoBulk(Section, BulkConcentration, MolarMass, km) #Same operation applied within function to every element in list (all nodes)
    KineticConstant = []
    Concentration = []
    
    for i in range(Section.NodeNumber): 
        if SaturationConcentration[i] > SolutionOxideConcentration[i]:
            x = Section.KdFe3O4electrochem[i]
        else:
            x = Section.KpFe3O4electrochem[i]
        KineticConstant.append(x) #Based on saturation behaviour at each node, e.g., [kd, kd, kp, kd, kp, etc.]
        #Generally, all nodes will be kd or all kp w/i a section, but this allows for inidividual calculation 
     
    OxideKinetics = InterfaceOxideKinetics(Section, KineticConstant, SaturationConcentration, MolarMass) #Same operation applied to all nodes

    Oxide = [x+y for x,y in zip(InnerIronOxThickness, OuterFe3O4Thickness)]
    
    for i in range(Section.NodeNumber):  
        if Element == "Fe": #same expression for dissolution (Oxide >0) or precipitation, subbing out appropriate kinetic constant 
            if SolutionOxideConcentration[i] >= SaturationConcentration[i]:
                y =  (Diff[i] + OxideKinetics[i] + BTrans[i])/(Section.KpFe3O4electrochem[i] + km[i])
            else: #SolutionOxideConcentration[i] < SaturationConcentration[i]:
                y = (Diff[i] + OxideKinetics[i] + BTrans[i])/(Section.KdFe3O4electrochem[i] + km[i])
            if Oxide[i] == 0:# Core 
                y = (Diff[i] + BTrans[i])/km[i] #No contribution from oxide kinetics
 
     
        if Element == "Co":               
            if Oxide[i] >0:# Precipitation within existing inner/outer layer or dissolution of Co from layer(s)
                if SolutionOxideConcentration[i] >= SaturationConcentration[i]:
                    y = (Diff[i] + OxideKinetics[i] + BTrans[i])/(Section.KpFe3O4electrochem[i] + km[i])
                else: #SolutionOxideConcentration[i] < SaturationConcentration[i]:
                    if CoThickness[i] > 0:
                        y = (Diff[i] + OxideKinetics[i] + BTrans[i])/(Section.KdFe3O4electrochem[i] + km[i])
                    else: #CoThickness == 0
                        y = (Diff[i] + BTrans[i])/km[i] #No contribution from oxide kinetics
            else: #Oxide[i] == 0 (Core) 
                y = (Diff[i] + BTrans[i])/km[i]
            
        if Element == "Ni": #Doesn't require existing Oxide layer for Ni precip (Ni(s))
            if SolutionOxideConcentration[i] >= SaturationConcentration[i]:
                y = (Diff[i] + OxideKinetics[i] + BTrans[i])/(Section.KpFe3O4electrochem[i] +km[i])
            else: #SolutionOxideConcentration < SaturationConcentration:
                if NiThickness[i] >0:
                    y = (Diff[i] + OxideKinetics[i] + BTrans[i])/(Section.KdFe3O4electrochem[i] +km[i])
                else:
                    y = (Diff[i] + BTrans[i])/km[i] 
        Concentration.append(y)
        
    Conc = ld.UnitConverter(Section, "Grams per Cm Cubed", "Mol per Kg", Concentration, None, None, None, MolarMass, None)
    return Conc    
        
            
def CorrosionRate(Section):
    ConcentrationFe2, ConcentrationFeOH2, ActivityCoefficient1, ActivityCoefficient2 = c.Hydrolysis(Section, Section.MetalOxide.FeTotal, \
                                                                            Section.MetalOxide.NiTotal, Section.MetalOxide.ConcentrationH)
    ProductConcentration = Section.MetalOxide.ConcentrationH2
    
    EqmPotentialH2 = [] #x
    EqmPotentialFe = [] #y
    ExchangeCurrentFe = [] #z
    ExchangeCurrentH2onFe = [] #w
    
    for i in range(Section.NodeNumber):
        
        x = e.Potential(Section, Section.StandardEqmPotentialH2[i], ProductConcentration[i], 1, 1, Section.MetalOxide.ConcentrationH[i], \
                        ActivityCoefficient1[i], 2, Section.DensityH2O[i], Section.NernstConstant[i], "gas")
        
        y = e.Potential(Section, Section.StandardEqmPotentialFe[i], 1, 1, 1, ConcentrationFe2[i], ActivityCoefficient2[i], 1, \
                        Section.DensityH2O[i], Section.NernstConstant[i], "aqueous")
        
        if Section == ld.Inlet or Section == ld.Outlet:
            w = e.ExchangeCurrentDensity(Section, nc.ActivationEnergyH2onFe, Section.MetalOxide.ConcentrationH[i], x, Section.DensityH2O[i], \
                                         Section.PrimaryBulkTemperature[i], "Acceptor")
            z = e.ExchangeCurrentDensity(Section, nc.ActivationEnergyFe, ConcentrationFe2[i], y, Section.DensityH2O[i], \
                                         Section.PrimaryBulkTemperature[i], "Acceptor")
            
        elif Section == ld.SteamGenerator:
            #EqmPotentialNi = e.Potential(Section, Section.StandardEqmPotentialNi, [1]*Section.NodeNumber, 1, 1, Composition.ConcentrationNi2, Composition.ActivityCoefficient2, 1, "aqueous")
            w = e.ExchangeCurrentDensity(Section, nc.ActivationEnergyH2onAlloy800, Section.MetalOxide.ConcentrationH[i], x, Section.DensityH2O[i], \
                                         Section.PrimaryBulkTemperature[i], "Acceptor")
            z = e.ExchangeCurrentDensity(Section, nc.ActivationEnergyAlloy800, ConcentrationFe2[i], y, Section.DensityH2O[i], \
                                         Section.PrimaryBulkTemperature[i], "Acceptor")
            
        EqmPotentialH2.append(x)
        EqmPotentialFe.append(y)
        ExchangeCurrentFe.append(z)
        ExchangeCurrentH2onFe.append(w)
    
        
    MixedECP = e.MixedPotential(Section, ExchangeCurrentH2onFe, EqmPotentialH2, ExchangeCurrentFe, EqmPotentialFe)
    CorrosionCurrent = [x*(np.exp((nc.Beta*nc.n*nc.F*(y-z))/(nc.R*q))-np.exp((-(1-nc.Beta)*nc.n*nc.F*(y-z))/(nc.R*q))) for x,y,z,q in zip(ExchangeCurrentFe,MixedECP,EqmPotentialFe,Section.PrimaryBulkTemperature)]
    
        
    if Section == ld.SteamGenerator: #icorr = ia = -ic   equivalent molar mass for Alloy 800 (0.46 Fe, 0.33 Ni, -> 2e- processes; 0.21 Cr -> 3 e- process)   
        Constant = [(1/nc.F)*((nc.FeMolarMass*0.46/nc.n) + (nc.NiMolarMass*0.33/nc.n) + (0.21*nc.CrMolarMass/3))]*Section.NodeNumber
        
    else:
        Constant = [nc.FeMolarMass/(nc.n*nc.F)]*Section.NodeNumber
        
    rate = [x*y for x,y in zip(CorrosionCurrent, Constant)] #[g/cm^2*s] 
    if Section == ld.Core:
        rate = [0]*Section.NodeNumber
    
    return rate, MixedECP 
    