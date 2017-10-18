import LepreauData as ld
import NumericConstants as nc
import numpy as np
import Composition as c

def Potential(Section, StandardPotential, ProductConcentration, ProductGamma, ProductStoich, ReactantConcentration, ReactantGamma, \
              ReactantStoich, DensityH2O, NernstConstant, Phase):
    
    if ProductConcentration <0:
        print ("Error in electrochemistry calculation: concentration < 0", ProductConcentration, Section.InnerOxThickness, "Inner Oxide g/cm2", Section.OuterOxThickness, "Outer Oxide g/cm2")
        ProductConcentration = 5e-10
        print ("Concentration reset to", ProductConcentration)
    
    if Phase == "gas":
        Product = (ProductConcentration*DensityH2O/nc.kH2)**ProductStoich
        Reactant = (ReactantConcentration*DensityH2O*ReactantGamma)**ReactantStoich 
    
    else:
        Product = (ProductConcentration*DensityH2O*ProductGamma)**ProductStoich 
        Reactant = (ReactantConcentration*DensityH2O*ReactantGamma)**ReactantStoich 
             
    E = StandardPotential- NernstConstant*np.log10(Product/Reactant)
    return E

def ElectrochemicalKineticConstant(Section, KineticConstant, EquilibriumPotential, Dissolution, MixedECP):
    if Dissolution == "yes":
        Beta_prime = nc.Beta*(-1)
    else:
        Beta_prime = nc.Beta
    
    return [nc.SAFactor*x*np.exp(Beta_prime*nc.n*nc.F*(y-z)/(nc.R*q)) for x,y,z,q in zip(KineticConstant,MixedECP,EquilibriumPotential,Section.PrimaryBulkTemperature)] 


def ElectrochemicalSaturation(Section, BulkSatFe3O4, EqmPotentialFe3O4, MixedPotential, Kelvin, Dissolution):
    if Dissolution == "yes":
        Beta_prime = (-1)*nc.Beta
    else:
        Beta_prime = nc.Beta
    
    AdjustedSaturation = BulkSatFe3O4*np.exp(Beta_prime*nc.n*nc.F*(MixedPotential-EqmPotentialFe3O4)/(nc.R*Kelvin)) 

    return AdjustedSaturation

def ElectrochemicalAdjustment(Section, EqmPotentialFe3O4, SOMixedPotential, MOMixedPotential, FeTotal, FeSatFe3O4, BulkSatFe3O4, SOConcentrationH):
    
    KpFe3O4electrochem = ElectrochemicalKineticConstant(Section, [nc.KpFe3O4]*Section.NodeNumber, EqmPotentialFe3O4, "no", SOMixedPotential)
    KdFe3O4electrochem = ElectrochemicalKineticConstant(Section, [nc.KdFe3O4]*Section.NodeNumber, EqmPotentialFe3O4, "yes", SOMixedPotential)
    
    AdjustedSaturation = []
    for i in range(Section.NodeNumber):
        if FeTotal[i] >= FeSatFe3O4[i]:
            
            x = ElectrochemicalSaturation(Section, BulkSatFe3O4[i], EqmPotentialFe3O4[i], SOMixedPotential[i], Section.PrimaryBulkTemperature[i], "no")
        else:
            x =ElectrochemicalSaturation(Section, BulkSatFe3O4[i], EqmPotentialFe3O4[i], SOMixedPotential[i], Section.PrimaryBulkTemperature[i], "yes")
        AdjustedSaturation.append(x)    

    
    if Section != ld.Core:
        MOConcentrationH = SOConcentrationH#[x*np.exp(-(y-z)*nc.F/(nc.R*t)) for x,y,z,t in zip(SOConcentrationH, MOMixedPotential, SOMixedPotential, Section.PrimaryBulkTemperature)]
    else: MOConcentrationH = None
    
    return KpFe3O4electrochem, KdFe3O4electrochem, AdjustedSaturation, MOConcentrationH                                                                      
    

def MixedPotential(Section, CathodeCurrent, CathodePotential, AnodeCurrent, AnodePotential):
    Numerator = [x*np.exp((nc.Beta*nc.n*nc.F*y)/(nc.R*z))+q*np.exp((nc.Beta*nc.n*nc.F*r)/(nc.R*z)) for x,y,z,q,r, in zip(CathodeCurrent,CathodePotential,Section.PrimaryBulkTemperature,AnodeCurrent,AnodePotential)]
    
    Denominator = [x*np.exp((-(1-nc.Beta)*nc.n*nc.F*y)/(nc.R*z))+q*np.exp((-(1-nc.Beta)*nc.n*nc.F*r)/(nc.R*z)) for x,y,z,q,r, in zip(CathodeCurrent,CathodePotential,Section.PrimaryBulkTemperature,AnodeCurrent,AnodePotential)]
    
    return  [(nc.R*x/(nc.n*nc.F))*np.log(y/z) for x,y,z in zip(Section.PrimaryBulkTemperature, Numerator, Denominator)]
    

def ExchangeCurrentDensity(Section, ActivationE, Concentration, Potential, DensityH2O, Kelvin, Species):
    #A^z+ + ze- -> D,   where A = e- acceptor and D = e- donor (Bockris & Reddy - Modern Electrochemistry)
    #If in terms of acceptor:  io = (F(C_A)kT/h)*exp(-ActivationEnergy/RT)*exp((-BnF/RT)Eeqm), C donor/acceptor is [mol/cm^2] (raise to the 2/3)
    #If in terms of donor: io =(F(C_D)kT/h)*exp(-ActivationEnergy/RT)*exp(((1-B)nF/RT)Eeqm), the B term is positive in the second exponential
    
    #examples:
    #Fe -> Fe2+ + 2e-,  io = (F[Fe2]^(2/3)kT/h)exp(-ActivationEnergyFe/RT)exp(-BnFEeqm/RT) (acceptor)
    #2H+ + 2e- -> H2,   io =(F[H]^(2/3)kT/h)exp(-ActivationEnergyFe/RT)exp(-BnFEeqm/RT) (acceptor)
    #Fe3O4(s) + 2H+(aq) + 2H2O(l) + 2e- -> 3Fe(OH)2(aq), io = (F[Fe(OH)2]^(2/3))kT/h)exp(-ActivationEnergyFe3O4/RT)exp((1-B)nFEeqm/RT), or just use +B (donor)
    #3Fe(OH)2(s) -> Fe3O4(s) + 2H+ + 2H2O + 2e-, io = (FkT/h)exp(-ActivationEnergyFe3O4/RT)exp((1-B)nFEeqm/RT), or just use +B (donor)
    if Species == "Acceptor":
        Beta_prime = nc.Beta*(-1)
    else:
        Beta_prime = nc.Beta
    if Concentration == 1:
        io= ((nc.F*nc.kb*Kelvin)/nc.hp)*np.exp(-ActivationE/(nc.R*Kelvin))*np.exp((Beta_prime*nc.n*nc.F*Potential)/(nc.R*Kelvin)) 
    else:  
        #print ([((nc.F*nc.kb*x)/nc.hp)*np.exp(-y/(nc.R*x))*((w*p/1000)**(2/3))*np.exp((Beta_prime*nc.n*nc.F*z)/(nc.R*x)) for x,y,z,w,p in zip(Section.PrimaryBulkTemperature, ActivationEnergy, Potential, Concentration, Section.DensityH2O)])
        io= ((nc.F*nc.kb*Kelvin)/nc.hp)*np.exp(-ActivationE/(nc.R*Kelvin))*\
            ((Concentration*DensityH2O/1000)**(2/3))*np.exp((Beta_prime*nc.n*nc.F*Potential)/(nc.R*Kelvin))    
    
    return io
    #Convert concentrations from mol/kg to mol/cm^2 for io eq'n: [mol/kg]/1000 -> [mol/g]*[g/cm^3] -> [mol/cm^3]^(2/3) -> mol^(2/3)/cm^2    


def Energy(Section,ExchangeCurrent, Concentration, Acceptor, EquilibriumPotential):
    if Acceptor == "yes":
        return [nc.R*x*-np.log10(y*nc.hp/(nc.F*((z*q/1000)**(2/3))*nc.kb*x*np.exp(-nc.Beta*nc.n*nc.F*w/(nc.R*x)))) for x,y,z,q,w in (Section.PrimaryBulkTemperature,ExchangeCurrent,Concentration,Section.DensityH2O, EquilibriumPotential)]
    else: #(1-Beta) vs. -Beta (sign change)
        return [nc.R*x*-np.log10(y*nc.hp/(nc.F*((z*q/1000)**(2/3))*nc.kb*x*np.exp(nc.Beta*nc.n*nc.F*w/(nc.R*x)))) for x,y,z,q,w in (Section.PrimaryBulkTemperature,ExchangeCurrent,Concentration,Section.DensityH2O, EquilibriumPotential)]
        
    if Concentration == 1 and Acceptor == "yes":
        return [nc.R*x*-np.log10(y*nc.hp/(nc.F*nc.kb*x*np.exp(-nc.Beta*nc.n*nc.F*w/(nc.R*x)))) for x,y,w in (Section.PrimaryBulkTemperature,ExchangeCurrent, EquilibriumPotential)]
    elif Concentration == 1 and Acceptor == "no":
        return [nc.R*x*-np.log10(y*nc.hp/(nc.F*nc.kb*x*np.exp(nc.Beta*nc.n*nc.F*w/(nc.R*x)))) for x,y,w in (Section.PrimaryBulkTemperature,ExchangeCurrent, EquilibriumPotential)]


def ECP(Section):
    for i in range(Section.NodeNumber):
        if Section.SolutionOxide.FeTotal[i] < 0:
            print ("Error: Fe concentration <0 in ECP function")
            Section.SolutionOxide.FeTotal[i]= 5e-10
            print ("Corrected, concentration reset to:", Section.SolutionOxide.FeTotal[i])
             
    ConcentrationFe2, ConcentrationFeOH2, ActivityCoefficient1, ActivityCoefficient2 = c.Hydrolysis(Section, Section.SolutionOxide.FeTotal, \
                        Section.SolutionOxide.NiTotal, Section.SolutionOxide.ConcentrationH)
    #ProductConcentration = Section.MetalOxide.ConcentrationH2
    
    ProductConcentration = Section.SolutionOxide.ConcentrationH2

    EqmPotentialFe3O4 = [] #x
    ExchangeCurrentFe3O4 = [] #y
    ExchangeCurrentH2onFe3O4 = [] #z
    EqmPotentialH2 = [] #q
    
    for i in range (Section.NodeNumber):
        
        q = Potential(Section, Section.StandardEqmPotentialH2[i], ProductConcentration[i], 1, 1, Section.SolutionOxide.ConcentrationH[i], \
                      ActivityCoefficient1[i], 2, Section.DensityH2O[i], Section.NernstConstant[i], "gas")
        
        if Section.SolutionOxide.FeTotal[i] >= Section.SolutionOxide.FeSatFe3O4[i]:
            x = Potential(Section, Section.StandardEqmPotentialFe3O4oxid[i], 1, 1, 1, Section.SolutionOxide.ConcentrationH[i], \
                          ActivityCoefficient1[i], 2, Section.DensityH2O[i], Section.NernstConstant[i], "aqueous")
            
            y = ExchangeCurrentDensity(Section, nc.PrecipitationActivationEnergyFe3O4, 1, x, Section.DensityH2O[i], Section.PrimaryBulkTemperature[i], "Donor")
            
            z = ExchangeCurrentDensity(Section, nc.PrecipitationActivationEnergyH2onFe3O4, Section.SolutionOxide.ConcentrationH[i], q, \
                                       Section.DensityH2O[i], Section.PrimaryBulkTemperature[i], "Acceptor")
        
        if Section.SolutionOxide.FeSatFe3O4[i] > Section.SolutionOxide.FeTotal[i]:
            x = Potential(Section, Section.StandardEqmPotentialFe3O4red[i], ConcentrationFeOH2[i], 1, 3, Section.SolutionOxide.ConcentrationH[i], \
                          ActivityCoefficient1[i], 2, Section.DensityH2O[i], Section.NernstConstant[i], "aqueous")
            y = ExchangeCurrentDensity(Section, nc.DissolutionActivationEnergyFe3O4, ConcentrationFeOH2[i], x, Section.DensityH2O[i], \
                                       Section.PrimaryBulkTemperature[i], "Donor") 
            z = ExchangeCurrentDensity(Section, nc.DissolutionActivationEnergyH2onFe3O4, Section.SolutionOxide.ConcentrationH[i], q, Section.DensityH2O[i], \
                                       Section.PrimaryBulkTemperature[i], "Acceptor")
        
        EqmPotentialFe3O4.append(x)
        ExchangeCurrentFe3O4.append(y)
        ExchangeCurrentH2onFe3O4.append(z) 
        EqmPotentialH2.append(q)
         
    MixedECP = MixedPotential(Section, ExchangeCurrentFe3O4, EqmPotentialFe3O4, ExchangeCurrentH2onFe3O4, EqmPotentialH2)    
    
    return MixedECP, EqmPotentialFe3O4


