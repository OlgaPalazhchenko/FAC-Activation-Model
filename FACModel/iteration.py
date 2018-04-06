import lepreau_data as ld
import numpy as np
import thermochemistry_and_constants as nc
import composition as c
import electrochemistry as e

# Activation Energies [J/mol]
ACTIVATION_ENERGY_Fe = 265813.5395 #269850.7017
ACTIVATION_ENERGY_H2onFe = 263239.8834 #267277.0456

ACTIVATION_ENERGY_ALLOY800 = 321276.4779
ACTIVATION_ENERGYH2onALLOY800 = 310324.6997


def Diffusion(Section, Element):
    # Inner oxide density is determined as a weighted average of the densities of the 2 comprising oxides 
    if Section in ld.SteamGenerator or Section in ld.SteamGenerator_2:
        AlloyDensity = nc.Alloy800Density
        OxideDensity = (c.fraction_chromite(Section) * nc.FeCr2O4Density
        + (1 - c.fraction_chromite(Section)) * nc.NiFe2O4Density)
        FractionCo = nc.FractionCo_Alloy800
    
    if Section in ld.InletSections or Section in ld.OutletSections:
        AlloyDensity = nc.FeDensity
        OxideDensity = (c.fraction_chromite(Section) * nc.FeCr2O4Density
                        + (1 - c.fraction_chromite(Section)) * nc.Fe3O4Density)
        FractionCo = nc.FractionCo_CS
        
    if Section not in ld.FuelSections:
        if Element == "Fe":
            MetalRetainedInsideInnerOxide = c.fraction_metal_inner_oxide(Section, "Fe") \
            * (OxideDensity / AlloyDensity) * (1 - nc.Fe3O4Porosity_inner)  # [gFe/g CSlost]
            
        elif Element == "Ni":
            MetalRetainedInsideInnerOxide = c.fraction_metal_inner_oxide(Section, "Ni") \
            * (OxideDensity / AlloyDensity) * (1 - nc.Fe3O4Porosity_inner)
        if Element == "Co":
            MetalRetainedInsideInnerOxide = c.fraction_metal_inner_oxide(Section, "Co") \
            * (OxideDensity / AlloyDensity) * (1 - nc.Fe3O4Porosity_inner)
        
    if Element == "Co":
        return FractionCo - MetalRetainedInsideInnerOxide
    elif Element == "Ni" or Element == "Fe":
        
        return 1 - MetalRetainedInsideInnerOxide  # [g metal/ g alloy lost]


def MetalOxideInterfaceConcentration(
        Section, Element, SolutionOxideInterfaceConcentration, InnerOxThickness, OuterOxThickness, Corrosion
        ):
    if Section in ld.SteamGenerator or Section in ld.SteamGenerator_2:
        OxideDensity = nc.NiFe2O4Density
    
    if Section in ld.InletSections or Section in ld.OutletSections:
        OxideDensity = nc.Fe3O4Density
        
    if Element == "Fe":
        Diffusivity = [nc.FeDiffusivity] * Section.NodeNumber
        MolarMass = nc.FeMolarMass
    elif Element == "Ni":
        Diffusivity = [nc.NiDiffusivity] * Section.NodeNumber
        MolarMass = nc.NiMolarMass
    elif Element == "Co":
        Diffusivity = [nc.CoDiffusivity] * Section.NodeNumber
        MolarMass = nc.CoMolarMass
    
    PathLength = [(x * nc.Fe3O4Tortuosity / (OxideDensity * (1 - nc.Fe3O4Porosity_inner)))
                  + (y * nc.Fe3O4Tortuosity / (OxideDensity * (1 - nc.Fe3O4Porosity_outer)))
                  for x, y in zip(InnerOxThickness, OuterOxThickness)] 
    
    DiffusivityTerm = [x * nc.Fe3O4Porosity_inner / y for x, y in zip(Diffusivity, PathLength)]
    
    SolutionConcentration = ld.UnitConverter(
        Section, "Mol per Kg", "Grams per Cm Cubed", SolutionOxideInterfaceConcentration, None, None, None, MolarMass,
        None
        )
    
    MOConcentration = [(x * Diffusion(Section, Element) / y) + z for x, y, z in
                       zip(Corrosion, DiffusivityTerm, SolutionConcentration)]

    # #Pre-set corrosion rate for outlet (lower MOConc = higher FAC rate (e-chem effect))
    # inner ox decrease = path length decreases --> diffusivity term increases --> metal oxide concentration decreases --> corr rate incr
    # Can also reset the oxide thickness such that the minimum results in acceptable MO Conc and corrrate
    
    # if Section in ld.OutletSections: 
    # MOConcentration =[1.22e-8]*Section.NodeNumber #110 um/a rate at this fixed M/O Fe concentration
      
    return ld.UnitConverter(
        Section, "Grams per Cm Cubed", "Mol per Kg", MOConcentration, None, None, None, MolarMass, None
        )
    

def TransfertoBulk(Section, BulkConcentration, MolarMass, km):  # [g/cm^2 s]
    Concentration = ld.UnitConverter(
        Section, "Mol per Kg", "Grams per Cm Cubed", BulkConcentration, None, None, None, MolarMass, None
        )
    
    return [x * y for x, y in zip(km, Concentration)] 
    

def InterfaceOxideKinetics (Section, KineticConstant, SaturationConcentration, MolarMass):    
    Concentration = ld.UnitConverter(
        Section, "Mol per Kg", "Grams per Cm Cubed", SaturationConcentration, None, None, None, MolarMass, None
        )
    
    SOTerm = [x * y for x, y in zip(KineticConstant, Concentration)]
    return SOTerm
    

def SolutionOxide(
        Section, BulkConcentration, SaturationConcentration, SolutionOxideConcentration, InnerIronOxThickness,
        OuterFe3O4Thickness, NiThickness, CoThickness, Element
        ):

    km = ld.MassTransfer(Section)
    
    if Section in ld.FuelSections:
        Diff = [0] * Section.NodeNumber
    else:
        Diff = [i * Diffusion(Section, Element) for i in Section.CorrRate]
        if Section in ld.InletSections or Section in ld.OutletSections:
            if Element == "Ni": 
                Diff = [0] * Section.NodeNumber
    
    # Corrosion rate becomes negative --> Diff --> S/O Concentration 
    
    if Element == "Fe": 
        MolarMass = nc.FeMolarMass
    elif Element == "Ni": 
        MolarMass = nc.NiMolarMass
    elif Element == "Co": 
        MolarMass = nc.CoMolarMass
    else:
        print ("Error: unknown element specified")            
    # Same operation applied within function to every element in list (all nodes) 
    
    # Bulk concentration in mol/kg here, need to convert to g/cm^3

#     if Element == "Fe": print (BulkConcentration[0], BulkConcentration_gm_cm_3[0])
    BTrans = TransfertoBulk(Section, BulkConcentration, MolarMass, km)
    
    KineticConstant = []
    Concentration = []
    
    for i in range(Section.NodeNumber): 
        if SaturationConcentration[i] > SolutionOxideConcentration[i]:
            x = Section.KdFe3O4electrochem[i]
        else:
            x = Section.KpFe3O4electrochem[i]
        KineticConstant.append(x)  # Based on saturation behaviour at each node, e.g., [kd, kd, kp, kd, kp, etc.]
        # Generally, all nodes will be kd or all kp w/i a section, but this allows for inidividual calculation 
    
    # Same operation applied to all nodes
    OxideKinetics = InterfaceOxideKinetics(Section, KineticConstant, SaturationConcentration, MolarMass)

    Oxide = [x + y for x, y in zip(InnerIronOxThickness, OuterFe3O4Thickness)]
    
    for i in range(Section.NodeNumber):  
        if Element == "Fe":
            # same expression for dissolution (Oxide >0) or precipitation, subbing out appropriate kinetic constant 
            if SolutionOxideConcentration[i] >= SaturationConcentration[i]:
                y = (Diff[i] + OxideKinetics[i] + BTrans[i]) / (Section.KpFe3O4electrochem[i] + km[i])
#                 if Section == ld.OutletFeeder: print (Diff[i], OxideKinetics[i], BTrans[i], Section.KpFe3O4electrochem[i] + km[i], i)
            
            else:  # SolutionOxideConcentration[i] < SaturationConcentration[i]:
                y = (Diff[i] + OxideKinetics[i] + BTrans[i]) / (Section.KdFe3O4electrochem[i] + km[i])
#                 if Section == ld.OutletFeeder: print (Diff[i], OxideKinetics[i], BTrans[i])
            
            if Oxide[i] == 0:  
                y = (Diff[i] + BTrans[i]) / km[i]  # No contribution from oxide kinetics
                
            if Section in ld.FuelSections:  
                y = BTrans[i] / km[i]  # No contribution from oxide kinetics
 
        if Element == "Co":               
            if Oxide[i] > 0:  # Precipitation within existing inner/outer layer or dissolution of Co from layer(s)
                if SolutionOxideConcentration[i] >= SaturationConcentration[i]:
                    y = (Diff[i] + OxideKinetics[i] + BTrans[i]) / (Section.KpFe3O4electrochem[i] + km[i])
                else:  # SolutionOxideConcentration[i] < SaturationConcentration[i]:
                    if CoThickness[i] > 0:
                        y = (Diff[i] + OxideKinetics[i] + BTrans[i]) / (Section.KdFe3O4electrochem[i] + km[i])
                    else:  # CoThickness == 0
                        y = (Diff[i] + BTrans[i]) / km[i]  # No contribution from oxide kinetics
            else:  # Oxide[i] == 0 (Core) 
                y = (Diff[i] + BTrans[i]) / km[i]
            
        if Element == "Ni":  # Doesn't require existing Oxide layer for Ni precip (Ni(s))
            if SolutionOxideConcentration[i] >= SaturationConcentration[i]:
                y = (Diff[i] + OxideKinetics[i] + BTrans[i]) / (Section.KpFe3O4electrochem[i] + km[i])
            else:  # SolutionOxideConcentration < SaturationConcentration:
                if NiThickness[i] > 0:
                    y = (Diff[i] + OxideKinetics[i] + BTrans[i]) / (Section.KdFe3O4electrochem[i] + km[i])
                else:
                    y = (Diff[i] + BTrans[i]) / km[i] 
        Concentration.append(y)
      
    Conc = ld.UnitConverter(
        Section, "Grams per Cm Cubed", "Mol per Kg", Concentration, None, None, None, MolarMass, None
        )
    return Conc    
        
            
def FAC_solver(Section, ConstantRate):
    
    # updates hydrolysis distribution of Fe species. Even if FAC rate kept constant, oxide thickness changes
    # M/O Fe total concentration changes w.r.t. thickness, so species cncentrations change too
    [ConcentrationFe2, ConcentrationFeOH2, ActivityCoefficient1, ActivityCoefficient2] = c.hydrolysis(
        Section, Section.MetalOxide.FeTotal, Section.MetalOxide.ConcentrationH
        )
    ProductConcentration = Section.MetalOxide.ConcentrationH2
    
    EqmPotentialH2 = []  # x
    EqmPotentialFe = []  # y
    ExchangeCurrentFe = []  # z
    ExchangeCurrentH2onFe = []  # w
    
    for i in range(Section.NodeNumber):
        
        x = e.potential(
            Section, Section.StandardEqmPotentialH2[i], ProductConcentration[i], 1, 1,
            Section.MetalOxide.ConcentrationH[i], ActivityCoefficient1[i], 2, Section.DensityH2O[i],
            Section.NernstConstant[i], "gas"
            )
        
        y = e.potential(
            Section, Section.StandardEqmPotentialFe[i], 1, 1, 1, ConcentrationFe2[i], ActivityCoefficient2[i], 1,
            Section.DensityH2O[i], Section.NernstConstant[i], "aqueous"
            )
        
        if Section in ld.InletSections or Section in ld.OutletSections:
            w = e.exchangecurrentdensity(
                Section, ACTIVATION_ENERGY_H2onFe, Section.MetalOxide.ConcentrationH[i], x, Section.DensityH2O[i],
                Section.PrimaryBulkTemperature[i], "Acceptor"
                )
            z = e.exchangecurrentdensity(
                Section, ACTIVATION_ENERGY_Fe, ConcentrationFe2[i], y, Section.DensityH2O[i],
                Section.PrimaryBulkTemperature[i], "Acceptor"
                )
            
        if Section in ld.SteamGenerator or Section in ld.SteamGenerator_2:
            # EqmPotentialNi = e.potential(
            # Section, Section.StandardEqmPotentialNi, [1]*Section.NodeNumber, 1, 1, composition.ConcentrationNi2,
            # composition.ActivityCoefficient2, 1, "aqueous"
            # )
            
            # assumed that different activation energies for half-cells if redox occurs on Alloy-800 vs. carbon steel
            w = e.exchangecurrentdensity(
                Section, ACTIVATION_ENERGYH2onALLOY800, Section.MetalOxide.ConcentrationH[i], x,
                Section.DensityH2O[i], Section.PrimaryBulkTemperature[i], "Acceptor"
                )
            
            z = e.exchangecurrentdensity(
                Section, ACTIVATION_ENERGY_ALLOY800, ConcentrationFe2[i], y, Section.DensityH2O[i],
                Section.PrimaryBulkTemperature[i], "Acceptor"
                )
            
        EqmPotentialH2.append(x)
        EqmPotentialFe.append(y)
        ExchangeCurrentFe.append(z)
        ExchangeCurrentH2onFe.append(w)
        
    MixedECP = e.mixed_potential(Section, ExchangeCurrentH2onFe, EqmPotentialH2, ExchangeCurrentFe, EqmPotentialFe)
    
    if Section in ld.FuelSections:
        rate = [0] * Section.NodeNumber
    # preset desired FAC rate
#     elif Section in ld.OutletSections and ConstantRate == "yes":
#         rate = [1.8e-09, 2.6e-09, 1.8e-09, 1.60e-09, 1.50e-09, 1.60e-09, 1.60e-09, 9.00e-10, 1.6e-09]
#         rate = [i * 0.4 for i in rate]
    
    # corrosion current calculation not required of rate has been set as constant
    else:
        CorrosionCurrent = [x * (np.exp((nc.Beta * nc.n * nc.F * (y - z)) / (nc.R * q))
                               - np.exp((-(1 - nc.Beta) * nc.n * nc.F * (y - z)) / (nc.R * q))) for x, y, z, q in
                            zip(ExchangeCurrentFe, MixedECP, EqmPotentialFe, Section.PrimaryBulkTemperature)]
      
        if Section in ld.SteamGenerator or Section in ld.SteamGenerator_2:
            # icorr = ia = -ic  equivalent MW for Alloy 800 (0.46 Fe, 0.33 Ni, -> 2e- processes; 0.21 Cr -> 3 e- process)   
            Constant = [(1 / nc.F) * ((nc.FeMolarMass * 0.46 / nc.n)
                                  + (nc.NiMolarMass * 0.33 / nc.n)
                                  + (0.21 * nc.CrMolarMass / 3))] * Section.NodeNumber
            
        elif Section in ld.OutletSections or Section in ld.InletSections:
            Constant = [nc.FeMolarMass / (nc.n * nc.F)] * Section.NodeNumber
        else:
            Constant = 0
            
        rate = [x * y for x, y in zip(CorrosionCurrent, Constant)]  # [g/cm^2*s] 
    
    return rate, MixedECP 


def interface_concentrations(Section, ConstantRate, BulkConcentrations, Saturations, RK4_InnerIronOxThickness,
                                  RK4_OuterFe3O4Thickness, RK4_NiThickness, RK4_CoThickness, j):
        
    #Solves S/O elemental concentrations at current approximation of oxide thickness(es)
    # (needed inside SolutionOxideBalance function to determine if <> saturation)

    # Calculates S/O elemental concentrations based on updated oxide thicknesses at each time step
    SolutionOxideConcentrations = [
        Section.SolutionOxide.FeTotal, Section.SolutionOxide.NiTotal, Section.SolutionOxide.CoTotal
        ]
    SOConc = []
    # Excludes Cr concentrations --> purely based on stellite transport (all Cr-oxides assumed to be insoluble)
    for x, y, z, w in zip (SolutionOxideConcentrations, BulkConcentrations[0:3], Saturations, ["Fe", "Ni", "Co"]):
        q = SolutionOxide(
            Section, y, z, x, RK4_InnerIronOxThickness, RK4_OuterFe3O4Thickness, RK4_NiThickness, RK4_CoThickness,
            w
            )

        SOConc.append(q)
    Section.SolutionOxide.FeTotal, Section.SolutionOxide.NiTotal, Section.SolutionOxide.CoTotal = SOConc

    Section.SolutionOxide.MixedPotential, Section.SolutionOxide.EqmPotentialFe3O4 = e.ECP(Section)

    if Section in ld.FuelSections:
        Section.CorrRate = [0] * Section.NodeNumber
        Section.MetalOxide.FeTotalFe = [0] * Section.NodeNumber
        Section.MetalOxide.MixedPotential = [0] * Section.NodeNumber
    else:

        Section.MetalOxide.FeTotal = MetalOxideInterfaceConcentration(
            Section, "Fe", Section.SolutionOxide.FeTotal, RK4_InnerIronOxThickness, RK4_OuterFe3O4Thickness,
            Section.CorrRate
            )
        Section.CorrRate, Section.MetalOxide.MixedPotential = FAC_solver(Section, ConstantRate)

    if Section in ld.SteamGenerator or Section in ld.SteamGenerator_2:
        Section.MetalOxide.NiTotal = MetalOxideInterfaceConcentration(
            Section, "Ni", Section.SolutionOxide.NiTotal, Section.InnerOxThickness, Section.OuterOxThickness,
            Section.CorrRate
            )

    elif Section in ld.InletSections or Section in ld.OutletSections:
        Section.MetalOxide.NiTotal = [0] * Section.NodeNumber
        
    # Cobalt M/O interface concentration not currently tracked but can be enabled
#         MetalOxideCo = MetalOxideInterfaceConcentration(
#         Section, "Co", SolutionOxideCoTotal, InnerOxThickness, OuterOxThickness, CorrRate
#         )
    
    
    [
        Section.KpFe3O4electrochem, Section.KdFe3O4electrochem, Section.SolutionOxide.FeSatFe3O4,
        Section.MetalOxide.ConcentrationH
        ] = e.electrochemical_adjustment(
        Section, Section.SolutionOxide.EqmPotentialFe3O4, Section.SolutionOxide.MixedPotential,
        Section.MetalOxide.MixedPotential, Section.SolutionOxide.FeTotal, Section.SolutionOxide.FeSatFe3O4,
        Section.Bulk.FeSatFe3O4, Section.SolutionOxide.ConcentrationH, j
        )
    
    return None
