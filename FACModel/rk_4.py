import lepreau_data as ld
import thermochemistry_and_constants as nc
import composition as c
import numpy as np
import iteration as it
import electrochemistry as e
import random
import sg_heattransfer as SGHX
from lepreau_data import SteamGeneratorSections
from datetime import date, timedelta

# Spalling thermochemistry_and_constants depend heavily on particle size distribution
OUTLET_OUTER_SPALL_CONSTANT = 5
OUTLET_INNER_SPALL_CONSTANT = 5
INLET_OUTER_SPALL_CONSTANT = 5.00E+10  # Different units for inlet versus outlet (different functions)
INLET_INNER_SPALL_CONSTANT = 1.00E+4


def oxide_composition(
        Section, Element, OxideType, Outer, OuterFe3O4Loading, NiLoading, CoLoading, InnerIronOxLoading
        ):
    # Note: this function does not take list input
   
    if Outer == "yes":  # for determining outer layer only composition
        # Outer Fe3O4 layer
        if OxideType > 0:  # For the total oxide value in each node
                # Fraction of iron in Fe3O4 * fraction of magnetite in outer oxide
            if Element == "Fe":
                return OuterFe3O4Loading * c.fraction_metal_in_oxide(Section, Element, "Magnetite") / OxideType
                # Outer ox thickness includes Fe3O4 and any incorporated Ni and Co
            elif Element == "Ni":
                return NiLoading / OxideType
            elif Element == "Co":
                return CoLoading / OxideType
            elif Element == "Cr":
                # Cr assumed to be 100% retained in inner layer
                return 0
            else:
                None
        else:
            return 0  # composition within outer layer only

    if Outer == "no":  # inner oxide layer composition only
        if OxideType > 0:  # For the InnerIronOxLoading value in each node
            if Element == "Fe":
                return c.fraction_metal_inner_oxide(Section, Element)
            if OuterFe3O4Loading > 0:
                if Element == "Ni":
                    # no Ni or Co incorporation via precip - only if present in inner layer through diffusion
                    return c.fraction_metal_inner_oxide(Section, Element)
                if Element == "Co":
                    return c.fraction_metal_inner_oxide(Section, Element)
            else:  # InnerOxLoading includes Ni incorporation from precipitation
                if Element == "Ni":
                    return (NiLoading
                            + c.fraction_metal_inner_oxide(Section, Element) * InnerIronOxLoading) / OxideType
                if Element == "Co":
                    return (CoLoading
                            + c.fraction_metal_inner_oxide(Section, Element) * InnerIronOxLoading) / OxideType
            if Element == "Cr":
                return c.fraction_metal_inner_oxide(Section, "Cr")

            if OxideType == 0:  # InnerIronOxide = 0
                return "lol"
        else:
            return None


def spatial(Section, Solution, Bulk, km, Diameter, Velocity, Length, i):
    # Cylindrical pipe: Surface Area / Volume = Pi*Diameter*Length/(Pi*(Diameter^2)/4) = 4/Diameter
    # Allows conversion between amount/area (assuming uniform distribution across the area) to amount/volume
    
    Solution_Bulk = Solution - Bulk
    Delta = (km * 4 / (Velocity * Diameter)) * Solution_Bulk
    # [cm]*[1/cm]*[mol/kg] + [mol/kg] = [mol/kg]
    BulkConcentration = Bulk + Delta * Length  # [x + y*Length for x,y in zip(Bulk, Delta)] 
    
    return BulkConcentration 


def oxide_growth(
        Section, ElementTracking,  RK4_InnerIronOxLoading, RK4_OuterFe3O4Loading, RK4_NiLoading, RK4_CoLoading,
        j):

    MolarMasses = [nc.FeMolarMass, nc.FeMolarMass, nc.NiMolarMass, nc.NiMolarMass, nc.CoMolarMass, nc.CoMolarMass]
    
    # Converts all concentration units to be used within RK4 from [mol/kg] to [g/cm^3]
    Concentrations = [
        Section.SolutionOxide.FeTotal, Section.SolutionOxide.NiTotal, Section.SolutionOxide.CoTotal,
        Section.SolutionOxide.FeSatFe3O4, Section.SolutionOxide.NiSatFerrite, Section.SolutionOxide.CoSatFerrite
        ]
    ConvertedConcentrations = []
    for Conc, MW in zip(Concentrations, MolarMasses):  # Concentrations has 6 lists in it
        x = [(z / 1000) * (MW * y) for z, y in zip(Conc, Section.DensityH2O)]
        ConvertedConcentrations.append(x)
    FeTotal, NiTotal, CoTotal, FeSat, NiSat, CoSat = ConvertedConcentrations

    # growth functions for each solid corrosion product at each node based on respective element's saturation behaviour
    GrowthOuterMagnetite = []
    GrowthInnerIronOxide = []
    
    if ElementTracking == "yes":
        GrowthCobalt = []
        GrowthNickel = []
    
    
    for i in range(Section.NodeNumber):
        
        # magnetite (outer and inner)
        if RK4_OuterFe3O4Loading[i] > 0:
            # with outer layer present, inner layer growth only depends on corrosion rate/diffusion
            # not affected by oxide kinetics (only diffusion)
            if Section in ld.FuelSections:
                q = 0
            else:
                q = Section.CorrRate[i] * it.Diffusion(Section, "Fe") / Section.FractionFeInnerOxide  # dy/dt = f(y)
            
            if FeTotal[i] >= FeSat[i]:  # precipitation of outer layer magnetite 
                
                if Section in ld.SteamGenerator or Section in ld.SteamGenerator_2:
                
                    if Section.Length.label[i] == "preheater":
                        Q = Section.HeatFlux[i]
                        HFE = (1.0866 / 100) * (np.exp((9.4157e-3) * Q - (3.8084e-6) * (Q ** 2)) - 1)
                       
                    else:
                        HFE = 0
                else:
                    HFE = 0 
                    
                x = Section.KpFe3O4electrochem[i] * (FeTotal[i] - FeSat[i]) * (1 + 4 * HFE)
                
            else:  # FeTotal[i] < FeSat[i]:
                x = Section.KdFe3O4electrochem[i] * (FeTotal[i] - FeSat[i])
        
        else:  # OuterFe3O4Loading[i] ==0
            if FeTotal[i] >= FeSat[i]: # precipitation of outer layer magnetite
                if Section in ld.FuelSections:
                    q = 0
                else:
                    q = Section.CorrRate[i] * it.Diffusion(Section, "Fe") / Section.FractionFeInnerOxide 
        
                x = Section.KpFe3O4electrochem[i] * (FeTotal[i] - FeSat[i]) 
            
            else:  # if FeTotal[i] < FeSat[i]:
                if Section in ld.FuelSections:
                    q = 0
                else: # no outer layer and dissolution conditions - inner layer begins to dissolve 
                    OxideGrowth = (Section.CorrRate[i] * it.Diffusion(Section, "Fe") / Section.FractionFeInnerOxide)
                    Dissolution = Section.KdFe3O4electrochem[i] * (FeTotal[i] - FeSat[i])
                    
                    q = OxideGrowth + Dissolution

                x = 0  # nothing to dissolve (no outer layer and dissolution conditions)
            
        GrowthOuterMagnetite.append(x)
        GrowthInnerIronOxide.append(q)
        
        if ElementTracking == "yes":
            # cobalt incorporation
            if RK4_InnerIronOxLoading[i] > 0:  # oxide layer present for Co to incorporate into
                
                if CoTotal[i] >= CoSat[i]: # Co incorporation 
                    y = Section.KpFe3O4electrochem[i] * (CoTotal[i] - CoSat[i]) 
                
                else:  # CoTotal[i] < CoSat[i]
                    if RK4_CoLoading[i] > 0: # if any cobalt present in current layer 
                        y = Section.KdFe3O4electrochem[i] * (CoTotal[i] - CoSat[i])
                    else:  # CoLoading == 0
                        y = 0  # Nothing to dissolve (no cobalt in layer + dissolution conditions)
            else:  # InnerIronOxLoading[i] == 0 (Core)
                y = 0  # Nothing to incorporate into 

            GrowthCobalt.append(y)           

            # nickel incorporation or Ni(s) deposition (independent of presence of an oxide layer backbone)
            if NiTotal[i] >= NiSat[i]:
                z = Section.KpFe3O4electrochem[i] * (NiTotal[i] - NiSat[i])
            else:  # NiTotal[i] < NiSat[i]:
                if RK4_NiLoading[i] > 0:  # Ni must be present for its dissolution to occur
                    z = Section.KdFe3O4electrochem[i] * (NiTotal[i] - NiSat[i])
                else:
                    z = 0  # Nothing to dissolve
            GrowthNickel.append(z)
    
    
        else:
            GrowthNickel = [0] * Section.NodeNumber
            GrowthCobalt = [0] * Section.NodeNumber
    
    return GrowthInnerIronOxide, GrowthOuterMagnetite, GrowthNickel, GrowthCobalt
    

def pht_cleaning(Section, InnerOxide, OuterOxide, Year_Month_Day_Hour):
    
    #need to change cleaning efficiency for each individual clean 
    if Year_Month_Day_Hour == (1995, 5, 8, 8):
        CleaningEfficiency = 0.67
    elif Year_Month_Day_Hour == (2008, 3, 1, 14):
        CleaningEfficiency = 0.30 #(IR-33110-0039-001-A)
    else:
        CleaningEfficiency = 0 # no cleaning, oxide layers returned without reduction in thickness
        
    Inner = []
    Outer = []
    for i in range(Section.NodeNumber):
        if OuterOxide[i] > 0:
            OuterOxide[i] = OuterOxide[i] * (1 - CleaningEfficiency)
            InnerOxide[i] = InnerOxide[i]
        else:
            InnerOxide[i] = InnerOxide[i] * (1 - CleaningEfficiency) #inner layer reduced in thickness instead
            OuterOxide[i] = OuterOxide[i] # remains zero
            
    Inner.append(InnerOxide)
    Outer.append(OuterOxide)
    
    return InnerOxide, OuterOxide


def oxide_layers(Section, ConstantRate, Saturations, BulkConcentrations, ElementTracking, j, SGFastMode):  
    
    start = SGHX.YearStartup
    delta = timedelta(hours = j * nc.TIME_STEP)
    CalendarDate = start + delta
        
    Year_Month_Day_Hour = (CalendarDate.year, CalendarDate.month, CalendarDate.day, CalendarDate.hour) 
    Year_Month = (CalendarDate.year, CalendarDate.month)
    
    if Section in ld.SteamGenerator or Section in ld.SteamGenerator_2:
        
        if Section in ld.SteamGenerator:
            SteamGenerator = ld.SteamGenerator
        elif Section in ld.SteamGenerator_2:
            SteamGenerator = ld.SteamGenerator_2
        
        TotalSGTubeNumber = SGHX.total_tubes_plugged(SteamGenerator, Year_Month_Day_Hour)
        
        # here Section = a bundle from the SG, but tube picker needs entire SG to be passed through
        if Section in ld.SteamGenerator: 
            SteamGenerator = ld.SteamGenerator
        else:
            SteamGenerator = ld.SteamGenerator_2

        if SGHX.SGFastMode == "yes":
            
            SelectedTubes = SGHX.tube_picker(SGHX.Method, SteamGenerator)[0]
    
            if Section == SelectedTubes[0]: #first tube from list of desired tubes per each SG
                
                [Section.InnerIronOxLoading, Section.OuterFe3O4Loading] = pht_cleaning(
                Section, Section.InnerIronOxLoading, Section.OuterFe3O4Loading, Year_Month_Day_Hour)
            else:
                Section.InnerIronOxLoading = Section.InnerIronOxLoading
                Section.OuterFe3O4Loading = Section.OuterFe3O4Loading
        
        else: #all sg tubes run (each bundle passed through, not entire SG) 
            Cleaned = primaryside_cleaned_tubes(Section, TotalSGTubeNumber, Year_Month)
            
            if Section in Cleaned:
                [Section.InnerIronOxLoading, Section.OuterFe3O4Loading] = pht_cleaning(
                Section, Section.InnerIronOxLoading, Section.OuterFe3O4Loading, Year_Month_Day_Hour)
            else:
               None 
    else:
        None
    
    
    if Year_Month in SGHX.OutageYearsMonths:
         # no oxide growth anywhere in the PHT system
        Section.InnerIronOxLoading = Section.InnerIronOxLoading
        Section.OuterFe3O4Loading = Section.OuterFe3O4Loading
        Section.NiLoading = Section.NiLoading
        Section.CoLoading = Section.CoLoading
        
    else:
       #regular oxide growth
   
        if ConstantRate == "yes":
            # updates M/O and S/O concentrations based on oxide thickness
            it.interface_concentrations(
                Section, ConstantRate, BulkConcentrations, Saturations, Section.InnerIronOxLoading,
                Section.OuterFe3O4Loading, Section.NiLoading, Section.CoLoading, j
                )
            # oxide growth functions based on S/O and M/O concentrations
            [GrowthInnerIronOxide, GrowthOuterMagnetite, GrowthNickel, GrowthCobalt] = oxide_growth(
                    Section, ElementTracking,  Section.InnerIronOxLoading, Section.OuterFe3O4Loading,
                    Section.NiLoading, Section.CoLoading, j
                    )
            
            # uniform oxide growth (preset corrosion rate)
            Section.InnerIronOxLoading = [
                x + y * nc.TIME_INCREMENT for x, y in zip(Section.InnerIronOxLoading, GrowthInnerIronOxide)
                ]
            Section.OuterFe3O4Loading = [
                x + y * nc.TIME_INCREMENT for x, y in zip(Section.OuterFe3O4Loading, GrowthOuterMagnetite)
                ]
            
            if ElementTracking == "yes":
                Section.NiLoading = [x + y * nc.TIME_INCREMENT for x, y in zip(Section.NiLoading, GrowthNickel)]
                Section.CoLoading = [x + y * nc.TIME_INCREMENT for x, y in zip(Section.CoLoading, GrowthCobalt)]
            
            else:
                None
        
        # layers continue to grow  
        else:
        # FAC solver for rate   
            RK4_InnerIronOxLoading = Section.InnerIronOxLoading
            RK4_OuterFe3O4Loading = Section.OuterFe3O4Loading
            # If element tracking is off, Co/Ni thickness remains same as initial loadings
            RK4_NiLoading = Section.NiLoading
            RK4_CoLoading = Section.CoLoading
            
            
            L = []
            M = []
            
            if ElementTracking == "yes": # lists for Ni and Co growth
                N = []
                P = []
            
            for approximation in range(4):  # 4 approximations in the "RK4" method
                
                it.interface_concentrations(
                    Section, ConstantRate, BulkConcentrations, Saturations, RK4_InnerIronOxLoading,
                    RK4_OuterFe3O4Loading, RK4_NiLoading, RK4_CoLoading, j
                    )
                
                GrowthInnerIronOxide, GrowthOuterMagnetite, GrowthNickel, GrowthCobalt = oxide_growth(
                    Section, ElementTracking,  RK4_InnerIronOxLoading, RK4_OuterFe3O4Loading, RK4_NiLoading,
                    RK4_CoLoading, j
                    )
                
                # iterate using previously solved RK4 thickness: re-evaluates growth functions based on S/O + M/O 
                # concentrations using new thickness     
               
                [RK4_InnerIronOxLoading, a] = RK4(
                    Section, Section.InnerIronOxLoading, GrowthInnerIronOxide, approximation
                    )
                L.append(a)
        
                [RK4_OuterFe3O4Loading, b] = RK4(
                    Section, Section.OuterFe3O4Loading, GrowthOuterMagnetite, approximation
                    )
                M.append(b)
        
                if ElementTracking == "yes": # otherwise, no growth 
                    [RK4_CoLoading, c] = RK4(Section, Section.CoLoading, GrowthCobalt, approximation)
                    N.append(c)
            
                    [RK4_NiLoading, d] = RK4(Section, Section.NiLoading, GrowthNickel, approximation)
                    P.append(d)
        
            Section.InnerIronOxLoading = [
                x + (y + 2 * z + 2 * q + e) / 6 for x, y, z, q, e in zip(
                    Section.InnerIronOxLoading, L[0], L[1], L[2], L[3]
                    )
                ]
            Section.OuterFe3O4Loading = [
                x + (y + 2 * z + 2 * q + e) / 6 for x, y, z, q, e in zip(
                    Section.OuterFe3O4Loading, M[0], M[1], M[2], M[3]
                    )
                ]
            
            if ElementTracking == "yes": # otherwise, not updated 
                Section.CoLoading = [
                    x + (y + 2 * z + 2 * q + e) / 6 for x, y, z, q, e in zip(Section.CoLoading, N[0], N[1], N[2], N[3])
                    ]
                Section.NiLoading = [
                    x + (y + 2 * z + 2 * q + e) / 6 for x, y, z, q, e in zip(Section.NiLoading, P[0], P[1], P[2], P[3])
                    ]
                
        Layers = [Section.InnerIronOxLoading, Section.OuterFe3O4Loading, Section.CoLoading, Section.NiLoading]
        # 4 different layers at each node. If any thicknesses are negative due to dissolution of respective layer, 
        # thickness = 0
        for i in range(4):
            for x in range(Section.NodeNumber):
                if Layers[i][x] < 0:
                    Layers[i][x] = 0
    

    return None
    

def RK4(Section, InitialThickness, GrowthFunction, approximation):
    x = [i * nc.TIME_INCREMENT for i in GrowthFunction]  # f(t)*delta t
    if approximation < 2:
        thickness = [a + b / 2 for a, b in zip(InitialThickness, x)]
    elif approximation == 2:
        thickness = [a + b for a, b in zip(InitialThickness, x)]
    if approximation > 2:
        thickness = InitialThickness

    for i in range(Section.NodeNumber):
        if thickness[i] < 0:
            thickness[i] = 0

    return thickness, x


def particle_size():
    r = random.random()
    if 0 < r < 0.005:
        Size = 0.15
    elif 0.005 < r < 0.009:
        Size = 0.5
    elif 0.009 < r < 0.17:
        Size = 1
    elif 0.17 < r < 0.25:
        Size = 1.5
    elif 0.25 < r < 0.31:
        Size = 2
    elif 0.31 < r < 0.4:
        Size = 2.5
    elif 0.4 < r < 0.5:
        Size = 3.2
    elif 0.5 < r < 0.52:
        Size = 3.2
    elif 0.52 < r < 0.55:
        Size = 4
    elif 0.55 < r < 0.6:
        Size = 4.5
    elif 0.6 < r < 0.62:
        Size = 5
    elif 0.62 < r < 0.67:
        Size = 5.5
    elif 0.67 < r < 0.69:
        Size = 6
    elif 0.69 < r < 0.72:
        Size = 7
    elif 0.72 < r < 0.75:
        Size = 8
    elif 0.75 < r < 0.8:
        Size = 9
    elif 0.8 < r < 0.83:
        Size = 10
    elif 0.83 < r < 0.9:
        Size = 15
    elif 0.9 < r < 0.98:
        Size = 40
    elif 0.98 < r < 1:
        Size = 100
    else:
        print (r)

    return Size * 0.0001 * nc.Fe3O4Density  # [um to cm to g/cm^2]


def spalling_time(Section, Particle, SolutionOxideFeSat, SolutionOxideFeTotal, KdFe3O4, OuterOxLoading, Velocity):

    if SolutionOxideFeSat > SolutionOxideFeTotal:  # Outlet 
        if OuterOxLoading > 0:
            SpTime = (OUTLET_OUTER_SPALL_CONSTANT * Particle) / \
            ((Velocity ** 2) * nc.Fe3O4Porosity_outer * KdFe3O4 * (SolutionOxideFeSat - SolutionOxideFeTotal))
        
        else:  # no outer oxide layer 
            SpTime = (OUTLET_INNER_SPALL_CONSTANT * Particle) / \
            ((Velocity ** 2) * nc.Fe3O4Porosity_inner * KdFe3O4 * (SolutionOxideFeSat - SolutionOxideFeTotal))
    
    else:  # Inlet 
        if OuterOxLoading > 0:
            SpTime = (INLET_OUTER_SPALL_CONSTANT * Particle) / (Velocity ** 2)
        else:  # no outer oxide layer 
            SpTime = (INLET_INNER_SPALL_CONSTANT * Particle) / (Velocity ** 2)
    
    # [hr]        
    return SpTime


def layer_spalling(Layer, TotalOxideThickness, Spalling):
    # inner/outer oxthickness & spalling = [g/cm^2]
    # This assumes that solid corrosion product species (Fe3O4, Ni0.6Fe2.4O4, etc.), if present simultaneously, 
    # are uniformally distributed with depth/layer thickness
    if TotalOxideThickness > 0:
        Ox = Layer - (Layer / TotalOxideThickness) * Spalling 
        if Ox < 0:
            Ox = 0
        return Ox


def spall(Section, j, SimulationStart, ElapsedTime, SpallTime, ElementTracking):
    # Current time step's RK4 input from oxidegrowth function for each of InnerOx, InnerIronOx, OuterFe3O4, Co, and 
    # Ni at each node of current section 

    # Silences spalling for desired sections
    if Section not in ld.OutletSections:
        Section.Particle = [0] * Section.NodeNumber 

    ConvertedConcentrations = []
    Concentrations = [Section.SolutionOxide.FeSatFe3O4, Section.SolutionOxide.FeTotal]
    for i in range(2):
        x = ld.UnitConverter(
            Section, "Mol per Kg", "Grams per Cm Cubed", Concentrations[i], None, None, None, nc.FeMolarMass, None
            )
        ConvertedConcentrations.append(x)

    FeSat, FeTotal = ConvertedConcentrations

    for i in range(Section.NodeNumber):
        # Oxide totals for RK4 iterations (M/O Concentration depends on total oxide thickness)
        if Section.OuterFe3O4Loading[i] > 0:  # from previous time step
            # With outer magnetite layer present, Ni and Co incorporate into overall "outer" oxide layer
            if ElementTracking == "yes":
                Section.OuterOxLoading[i] = (Section.OuterFe3O4Loading[i]
                + Section.CoLoading[i]
                + Section.NiLoading[i])
            else:
                Section.OuterOxLoading[i] = Section.OuterFe3O4Loading[i]
            
            Section.InnerOxLoading[i] = Section.InnerIronOxLoading[i]
        
        else:  
            if ElementTracking == "yes":
                Section.InnerOxLoading[i] = (Section.InnerIronOxLoading[i]
                + Section.CoLoading[i]
                + Section.NiLoading[i])
            else:
                Section.InnerOxLoading[i] = Section.InnerIronOxLoading[i]
        
        
        if j == SimulationStart:
            # First time step call generate particle sizes and calc spalling times, respectively
            x = particle_size()
            # Amount of time it will take each node's particle to spall 
            y = spalling_time(
                Section, x, FeSat[i], FeTotal[i], Section.KdFe3O4electrochem[i], Section.OuterFe3O4Loading[i],
                Section.Velocity[i]
                )
            
            # list of random particle sizes based on input distribution in particle_size function
            Section.Particle.append(x)
            SpallTime.append(y)
            ElapsedTime = [0] * Section.NodeNumber  # No time has elapsed yet at first time step for all nodes

        else:  # after first time step
#             if Section ==ld.OutletFeeder: print (Section.Particle[i], SpallTime[i], ElapsedTime[i], j, i)
            if ElapsedTime[i] >= SpallTime[i]:
                # enough time elapsed for particle of that size (with respective spall time) to come off
                if Section.OuterOxLoading[i] > 0:
                    # [cm]*[g/cm^3] = [g/cm^2] of layer removed due to spalling
                    Section.OuterFe3O4Loading[i] = layer_spalling(
                        Section.OuterFe3O4Loading[i], Section.OuterOxLoading[i], Section.Particle[i]
                        ) 
                    if ElementTracking == "yes":
                        Section.CoLoading[i] = layer_spalling(
                            Section.CoLoading[i], Section.OuterOxLoading[i], Section.Particle[i]
                            )
                        Section.NiLoading[i] = layer_spalling(
                            Section.NiLoading[i], Section.OuterOxLoading[i], Section.Particle[i]
                        )
                    Section.OuterOxLoading[i] = layer_spalling(
                        Section.OuterOxLoading[i], Section.OuterOxLoading[i], Section.Particle[i]
                        )  # new total outer ox
                
                else:  # OuterOxLoading == 0, spalling takes place at inner layer  
                    Section.InnerIronOxLoading[i] = layer_spalling(
                        Section.InnerIronOxLoading[i], Section.InnerOxLoading[i], Section.Particle[i]
                        )
                    if ElementTracking == "yes":
                        
                        Section.CoLoading[i] = layer_spalling(
                            Section.CoLoading[i], Section.InnerOxLoading[i], Section.Particle[i]
                            )
                        Section.NiLoading[i] = layer_spalling(
                            Section.NiLoading[i], Section.InnerOxLoading[i], Section.Particle[i]
                            )
                    Section.InnerOxLoading[i] = layer_spalling(
                        Section.InnerOxLoading[i], Section.InnerOxLoading[i], Section.Particle[i]
                        )  # new total inner ox

                Section.Particle[i] = particle_size()
                # once particle spalled off, another size (just at that node) is randomly generated w/i distribution

                SpallTime[i] = spalling_time(
                    Section, Section.Particle[i], FeSat[i], FeTotal[i], Section.KdFe3O4electrochem[i],
                    Section.OuterFe3O4Loading[i], Section.Velocity[i]
                    )
                
                # Elapsed time needs to reset to zero after spalling event occurs
                ElapsedTime[i] = 0 

            else:  # not enough time has passed
                ElapsedTime[i] = ElapsedTime[i] + 1 * (nc.TIME_INCREMENT / 3600)

    for i in range(Section.NodeNumber):
        
        if Section in ld.OutletSections:  # Outlet header spalling corrector 

#             print (SpallTime, i)
            
            
            if SpallTime[i] > 4000:
                SpallTime[i] = 3000 
        
        if Section not in ld.FuelSections:
            if Section.InnerIronOxLoading[i] <= 5e-6:
                Section.InnerIronOxLoading[i] = 0.000025  # Resets to original thickness

    return ElapsedTime, SpallTime
