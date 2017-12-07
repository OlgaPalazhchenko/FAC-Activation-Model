import lepreau_data as ld
import constants as nc
import composition as c
import numpy as np
import iteration as it
import electrochemistry as e
import random


def oxide_composition(
        Section, Element, OxideType, Outer, OuterFe3O4Thickness, NiThickness, CoThickness, InnerIronOxThickness
        ):
    # Note: this function does not take list input
   
    if Outer == "yes":  # for determining outer layer only composition
        # Outer Fe3O4 layer
        if OxideType > 0:  # For the total oxide value in each node
                # Fraction of iron in Fe3O4 * fraction of magnetite in outer oxide
            if Element == "Fe":
                return OuterFe3O4Thickness * c.fraction_metal_in_oxide(Section, Element, "Magnetite") / OxideType
                # Outer ox thickness includes Fe3O4 and any incorporated Ni and Co
            elif Element == "Ni":
                return NiThickness / OxideType
            elif Element == "Co":
                return CoThickness / OxideType
            elif Element == "Cr":
                # Cr assumed to be 100% retained in inner layer
                return 0
            else:
                None
        else:
            return 0  # composition within outer layer only

    if Outer == "no":  # inner oxide layer composition only
        if OxideType > 0:  # For the InnerIronOxThickness value in each node
            if Element == "Fe":
                return c.fraction_metal_inner_oxide(Section, Element)
            if OuterFe3O4Thickness > 0:
                if Element == "Ni":
                    # no Ni or Co incorporation via precip - only if present in inner layer through diffusion
                    return c.fraction_metal_inner_oxide(Section, Element)
                if Element == "Co":
                    return c.fraction_metal_inner_oxide(Section, Element)
            else:  # InnerOxThickness includes Ni incorporation from precipitation
                if Element == "Ni":
                    return (NiThickness
                            + c.fraction_metal_inner_oxide(Section, Element) * InnerIronOxThickness) / OxideType
                if Element == "Co":
                    return (CoThickness
                            + c.fraction_metal_inner_oxide(Section, Element) * InnerIronOxThickness) / OxideType
            if Element == "Cr":
                return c.fraction_metal_inner_oxide(Section, "Cr")

            if OxideType == 0:  # InnerIronOxide = 0
                return 0
        else:
            return None


def spatial(Solution, Bulk, km, Diameter, Velocity, Length):
    # Cylindrical pipe: Surface Area / Volume = Pi*Diameter*Length/(Pi*(Diameter^2)/4) = 4/Diameter
    # Allows conversion between amount/area (assuming uniform distribution across the area) to amount/volume
    Solution_Bulk = Solution - Bulk
    Delta = (km * 4 / (Velocity * Diameter)) * Solution_Bulk
    # [cm]*[1/cm]*[mol/kg] + [mol/kg] = [mol/kg]
    BulkConccentration = Bulk + Delta * Length  # [x + y*Length for x,y in zip(Bulk, Delta)] 
    return BulkConccentration 


def constant_oxide_growth():
    # if the corrosion rate is manually set, no need for RK4 method to solve for oxides (no longer ODE)
    None 


def oxidegrowth(Section, Saturations, BulkConcentrations, ElementTracking):
    RK4_InnerIronOxThickness = Section.InnerIronOxThickness
    RK4_OuterFe3O4Thickness = Section.OuterFe3O4Thickness
    # If element tracking is off, Co/Ni thickness remains same as initial loadings
    RK4_NiThickness = Section.NiThickness
    RK4_CoThickness = Section.CoThickness
    
    L = []
    M = []
    
    if ElementTracking == "yes": # lists for Ni and Co growth
        N = []
        P = []
      
    MolarMasses = [nc.FeMolarMass, nc.FeMolarMass, nc.NiMolarMass, nc.NiMolarMass, nc.CoMolarMass, nc.CoMolarMass]

    for approximation in range(4):  # 4 approximations in the "RK4" method
        # #Solves S/O elemental concentrations at current approximation of oxide thickness(es)
        # (needed inside SolutionOxideBalance function to determine if <> saturation)

        # Calculates S/O elemental concentrations based on updated oxide thicknesses at each time step
        SolutionOxideConcentrations = [
            Section.SolutionOxide.FeTotal, Section.SolutionOxide.NiTotal, Section.SolutionOxide.CoTotal
            ]
        SOConc = []
        # Excludes Cr concentrations --> purely based on stellite transport (all Cr-oxides assumed to be insoluble)
        for x, y, z, w in zip (SolutionOxideConcentrations, BulkConcentrations[0:3], Saturations, ["Fe", "Ni", "Co"]):
            q = it.SolutionOxide(
                Section, y, z, x, RK4_InnerIronOxThickness, RK4_OuterFe3O4Thickness, RK4_NiThickness, RK4_CoThickness,
                w
                )

            SOConc.append(q)
        Section.SolutionOxide.FeTotal, Section.SolutionOxide.NiTotal, Section.SolutionOxide.CoTotal = SOConc

        Section.SolutionOxide.MixedPotential, Section.SolutionOxide.EqmPotentialFe3O4 = e.ECP(Section)

        if Section in ld.FuelChannels:
            Section.CorrRate = [0] * Section.NodeNumber
            Section.MetalOxide.FeTotalFe = [0] * Section.NodeNumber
            Section.MetalOxide.MixedPotential = [0] * Section.NodeNumber
        else:

            Section.MetalOxide.FeTotal = it.MetalOxideInterfaceConcentration(
                Section, "Fe", Section.SolutionOxide.FeTotal, RK4_InnerIronOxThickness, RK4_OuterFe3O4Thickness,
                Section.CorrRate
                )
            Section.CorrRate, Section.MetalOxide.MixedPotential = it.CorrosionRate(Section)

        if Section in ld.SGZones:
            Section.MetalOxide.NiTotal = it.MetalOxideInterfaceConcentration(
                Section, "Ni", Section.SolutionOxide.NiTotal, Section.InnerOxThickness, Section.OuterOxThickness,
                Section.CorrRate
                )

        elif Section in ld.InletSections or Section in ld.OutletSections:
            Section.MetalOxide.NiTotal = [0] * Section.NodeNumber
#             MetalOxideCo = it.MetalOxideInterfaceConcentration(
#             Section, "Co", SolutionOxideCoTotal, InnerOxThickness, OuterOxThickness, CorrRate
#             )

        [
            Section.KpFe3O4electrochem, Section.KdFe3O4electrochem, Section.SolutionOxide.FeSatFe3O4,
            Section.MetalOxide.ConcentrationH
            ] = \
            e.ElectrochemicalAdjustment(
            Section, Section.SolutionOxide.EqmPotentialFe3O4, Section.SolutionOxide.MixedPotential,
            Section.MetalOxide.MixedPotential, Section.SolutionOxide.FeTotal, Section.SolutionOxide.FeSatFe3O4,
            Section.Bulk.FeSatFe3O4, Section.SolutionOxide.ConcentrationH
            )

        # Converts all concentration units to be used within RK4 from [mol/kg] to [g/cm^3]
        Concentrations = [
            Section.SolutionOxide.FeTotal, Section.SolutionOxide.NiTotal, Section.SolutionOxide.CoTotal,
            Section.SolutionOxide.FeSatFe3O4, Section.SolutionOxide.NiSatFerrite, Section.SolutionOxide.CoSatFerrite
            ]
        ConvertedConcentrations = []
        for i, k in zip(Concentrations, MolarMasses):  # Concentrations has 6 lists in it
            x = [(z / 1000) * (k * y) for z, y in zip(i, Section.DensityH2O)]
            ConvertedConcentrations.append(x)
        FeTotal, NiTotal, CoTotal, FeSat, NiSat, CoSat = ConvertedConcentrations

        # Determines dy/dt (growth functions) for each type of solid corrosion product at each node based on 
        # respective element's saturation behaviour
        GrowthOuterMagnetite = []
        GrowthInnerIronOxide = []
        
        if ElementTracking == "yes":
            GrowthCobalt = []
            GrowthNickel = []

        # Iterate using previously solved RK4 thickness: re-evaluates growth functions based on S/O + M/O 
        # concentrations using new thickness     
        for i in range(Section.NodeNumber):
            # Magnetite (outer and inner)
            if RK4_OuterFe3O4Thickness[i] > 0:
                # With outer layer present, dFe3O4_i/dt only depends on corrosion rate/diffusion
                if Section in ld.FuelChannels:
                    q = 0
                else:
                    # While outer layer present, inner Fe3O4 not affected by oxide kinetics (only diffusion)
                    q = Section.CorrRate[i] * it.Diffusion(Section, "Fe") / Section.FractionFeInnerOxide  # dy/dt = f(y)
                if FeTotal[i] >= FeSat[i]:  # Precipitation of outer layer magnetite 
                # dy/dt = f(y) via Diff term in SO concentration
                    x = Section.KpFe3O4electrochem[i] * (FeTotal[i] - FeSat[i])
                else:  # if FeTotal[i] < FeSat[i]:
                    x = Section.KdFe3O4electrochem[i] * (FeTotal[i] - FeSat[i])
            
            else:  # OuterFe3O4Thickness[i] ==0
                if FeTotal[i] >= FeSat[i]:
                    if Section in ld.FuelChannels:
                        q = 0
                    else:
                        q = Section.CorrRate[i] * it.Diffusion(Section, "Fe") / Section.FractionFeInnerOxide 
                    
                    x = Section.KpFe3O4electrochem[i] * (FeTotal[i] - FeSat[i])
                else:  # if FeTotal[i] < FeSat[i]:
                    x = 0  # nothing to dissolve
                    if Section in ld.FuelChannels:
                        q = 0
                    else:
                        q = ((Section.CorrRate[i] * it.Diffusion(Section, "Fe") / Section.FractionFeInnerOxide)
                             + Section.KdFe3O4electrochem[i] * (FeTotal[i] - FeSat[i]))

            GrowthOuterMagnetite.append(x)
            GrowthInnerIronOxide.append(q)
            
            if ElementTracking == "yes":
                # Cobalt incorporation
                if RK4_InnerIronOxThickness[i] > 0:  # Oxide layer present for Co to incorporate into
                    if CoTotal[i] >= CoSat[i]:
                        y = Section.KpFe3O4electrochem[i] * (CoTotal[i] - CoSat[i]) 
                    else:  # CoTotal[i] < CoSat[i]
                        if RK4_CoThickness[i] > 0:
                            y = Section.KdFe3O4electrochem[i] * (CoTotal[i] - CoSat[i])
                        else:  # CoThickness == 0
                            y = 0  # Nothing to dissolve 
                else:  # InnerIronOxThickness[i] == 0 (Core)
                    y = 0  # Nothing to incorporate into 
    
                GrowthCobalt.append(y)           
    
                # Nickel incorporation or Ni(s) deposition (independent of presence of an oxide layer backbone)
                if NiTotal[i] >= NiSat[i]:
                    z = Section.KpFe3O4electrochem[i] * (NiTotal[i] - NiSat[i])
                else:  # NiTotal[i] < NiSat[i]:
                    if RK4_NiThickness[i] > 0:  # Ni must be present for its dissolution to occur
                        z = Section.KdFe3O4electrochem[i] * (NiTotal[i] - NiSat[i])
                    else:
                        z = 0  # Nothing to dissolve
                GrowthNickel.append(z)

        # Section.Oxidelayer = initial input from previous time step's RK4 solution 
        [RK4_InnerIronOxThickness, a] = RK4(Section, Section.InnerIronOxThickness, GrowthInnerIronOxide, approximation)
        L.append(a)

        [RK4_OuterFe3O4Thickness, b] = RK4(Section, Section.OuterFe3O4Thickness, GrowthOuterMagnetite, approximation)
        M.append(b)

        if ElementTracking == "yes": # otherwise, no growth 
            [RK4_CoThickness, c] = RK4(Section, Section.CoThickness, GrowthCobalt, approximation)
            N.append(c)
    
            [RK4_NiThickness, d] = RK4(Section, Section.NiThickness, GrowthNickel, approximation)
            P.append(d)

#         for i in range(Section.NodeNumber):
#             # Need the overall inner and outer oxides to be updated for M/O concentration
#             if RK4_OuterFe3O4Thickness[i] > 0:  # from previous time step
#                 # With outer magnetite layer present, Ni and Co incorporate into overall "outer" oxide layer
#                 Section.OuterOxThickness[i] = RK4_OuterFe3O4Thickness[i] + RK4_CoThickness[i] + RK4_NiThickness[i]
#                 Section.InnerOxThickness[i] = RK4_InnerIronOxThickness[i]
#             else:  # OuterFe3O4Thickness == 0
#                 Section.InnerOxThickness[i] = RK4_InnerIronOxThickness[i] + RK4_CoThickness[i] + RK4_NiThickness[i]

    Section.InnerIronOxThickness = [
        x + (y + 2 * z + 2 * q + e) / 6 for x, y, z, q, e in zip(Section.InnerIronOxThickness, L[0], L[1], L[2], L[3])
        ]
    Section.OuterFe3O4Thickness = [
        x + (y + 2 * z + 2 * q + e) / 6 for x, y, z, q, e in zip(Section.OuterFe3O4Thickness, M[0], M[1], M[2], M[3])
        ]
    
    if ElementTracking == "yes": # otherwise, not updated 
        Section.CoThickness = [
            x + (y + 2 * z + 2 * q + e) / 6 for x, y, z, q, e in zip(Section.CoThickness, N[0], N[1], N[2], N[3])
            ]
        Section.NiThickness = [
            x + (y + 2 * z + 2 * q + e) / 6 for x, y, z, q, e in zip(Section.NiThickness, P[0], P[1], P[2], P[3])
            ]

#     Layers = [Section.InnerIronOxThickness, Section.OuterFe3O4Thickness, Section.CoThickness, Section.NiThickness]
    # 4 different layers at each node. If any thicknesses are negative due to dissolution of respective layer, 
    # thickness = 0
#     for i in range(4):
#         for x in range(Section.NodeNumber):
#             if Layers[i][x] < 0:
#                 Layers[i][x] = 0
    return None


def RK4(Section, InitialThickness, GrowthFunction, approximation):
    x = [i * nc.TimeIncrement for i in GrowthFunction]  # f(t)*delta t
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


def spalling_time(
        Section, Particle, SolutionOxideFeSat, SolutionOxideFeTotal, KdFe3O4electrochem, OuterOxThickness, Velocity
        ):

    if SolutionOxideFeSat > SolutionOxideFeTotal:  # Outlet 
        if OuterOxThickness > 0:
            SpTime = (nc.OutletOuterSpallConstant * Particle / 3600) / (
                (Velocity ** 2) * nc.Fe3O4Porosity_outer * KdFe3O4electrochem * (SolutionOxideFeSat - SolutionOxideFeTotal)
                )
        else:  # no outer oxide layer 
            SpTime = (nc.OutletInnerSpallConstant * Particle / 3600) / (
                (Velocity ** 2) * nc.Fe3O4Porosity_inner * KdFe3O4electrochem * (SolutionOxideFeSat - SolutionOxideFeTotal)
                )
    else:  # Inlet 
        if OuterOxThickness > 0:
            SpTime = (nc.InletOuterSpallConstant * Particle / 3600) / (Velocity ** 2)
        else:  # no outer oxide layer 
            SpTime = (nc.InletInnerSpallConstant * Particle / 3600) / (Velocity ** 2)

    return SpTime


def oxide(Layer, TotalOxideThickness, Spalling):
    # inner/outer oxthickness & spalling = [g/cm^2]
    # This assumes that solid corrosion product species (Fe3O4, Ni0.6Fe2.4O4, etc.), if present simultaneously, 
    # are uniformally distributed with depth/layer thickness
    if TotalOxideThickness > 0:
        Ox = Layer - (Layer / TotalOxideThickness) * Spalling 
        if Ox < 0:
            Ox = 0
        return Ox


def spall(Section, j, ElapsedTime, SpallTime, ElementTracking):
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

        # Oxide totals for RK4 iterations (M/O Concentration depends on total ox thickness) and before spalling function
        if Section.OuterFe3O4Thickness[i] > 0:  # from previous time step
            # With outer magnetite layer present, Ni and Co incorporate into overall "outer" oxide layer
            Section.OuterOxThickness[i] = (Section.OuterFe3O4Thickness[i]
                + Section.CoThickness[i]
                + Section.NiThickness[i])  
#             if ElementTracking == "yes":
#                 Section.OuterOxThickness[i] = (Section.OuterFe3O4Thickness[i]
#                 + Section.CoThickness[i]
#                 + Section.NiThickness[i])
#             else:
#                 Section.OuterOxThickness[i] = Section.OuterFe3O4Thickness[i]
            
            
            Section.InnerOxThickness[i] = Section.InnerIronOxThickness[i]
        else:  # OuterFe3O4Thickness == 0
            Section.InnerOxThickness[i] = (Section.InnerIronOxThickness[i]
            + Section.CoThickness[i]
            + Section.NiThickness[i])

        if j == 0:
            # First time step call generate particle sizes and calc spalling times, respectively
            x = particle_size()
            # Amount of time it will take each node's particle to spall 
            y = spalling_time(
                Section, x, FeSat[i], FeTotal[i], Section.KdFe3O4electrochem[i], Section.OuterFe3O4Thickness[i],
                Section.Velocity[i]
                )
            
            # list of random particle sizes based on input distribution in particle_size function
            Section.Particle.append(x)
            SpallTime.append(y)

            if Section in ld.OutletSections:  # Outlet header spalling corrector 
                if SpallTime[i] >= 4000:
                    SpallTime[i] = 2000 

            ElapsedTime = [0] * Section.NodeNumber  # No time has elapsed yet at first time step for all nodes

        else:  # after first time step
            if ElapsedTime[i] >= SpallTime[i]:
                # enough time elapsed for particle of that size (with respective spall time) to come off
                if Section.OuterOxThickness[i] > 0:
                    # [cm]*[g/cm^3] = [g/cm^2] of layer removed due to spalling
                    Section.OuterFe3O4Thickness[i] = oxide(
                        Section.OuterFe3O4Thickness[i], Section.OuterOxThickness[i], Section.Particle[i]
                        ) 
                    if ElementTracking == "yes":
                        Section.CoThickness[i] = oxide(
                            Section.CoThickness[i], Section.OuterOxThickness[i], Section.Particle[i]
                            )
                        Section.NiThickness[i] = oxide(
                            Section.NiThickness[i], Section.OuterOxThickness[i], Section.Particle[i]
                        )
                    Section.OuterOxThickness[i] = oxide(
                        Section.OuterOxThickness[i], Section.OuterOxThickness[i], Section.Particle[i]
                        )  # new total outer ox
                
                else:  # OuterOXThickness == 0, spalling takes place at inner layer  
                    Section.InnerIronOxThickness[i] = oxide(
                        Section.InnerIronOxThickness[i], Section.InnerOxThickness[i], Section.Particle[i]
                        )
                    if ElementTracking == "yes":
                        
                        Section.CoThickness[i] = oxide(
                            Section.CoThickness[i], Section.InnerOxThickness[i], Section.Particle[i]
                            )
                        Section.NiThickness[i] = oxide(
                            Section.NiThickness[i], Section.InnerOxThickness[i], Section.Particle[i]
                            )
                    Section.InnerOxThickness[i] = oxide(
                        Section.InnerOxThickness[i], Section.InnerOxThickness[i], Section.Particle[i]
                        )  # new total inner ox

                Section.Particle[i] = particle_size()
                # once particle spalled off, another size (just at that node) is randomly generated w/i distribution

                SpallTime[i] = spalling_time(
                    Section, Section.Particle[i], FeSat[i], FeTotal[i], Section.KdFe3O4electrochem[i],
                    Section.OuterFe3O4Thickness[i], Section.Velocity[i]
                    )

                if Section in ld.OutletSections:  # Outlet header spalling corrector
                    if SpallTime[i] >= 4000:
                        SpallTime[i] = 2000
                ElapsedTime[i] = 0
                # once a particle has "spalled" off, elapsed time since spalling resets to zero and counter
                # restarts at that node

            else:  # not enough time has passed
                ElapsedTime[i] = ElapsedTime[i] + 1 * (nc.TimeIncrement / 3600)

    for i in range(Section.NodeNumber):
        if Section not in ld.FuelChannels:
            if Section.InnerIronOxThickness[i] <= 8e-6:
                Section.InnerIronOxThickness[i] = 0.00025  # Resets to original thickness

            # Oxide totals for RK4 iterations (M/O Concentration depends on total oxide thickness)
        if Section.OuterFe3O4Thickness[i] > 0:  # from previous time step
            # With outer magnetite layer present, Ni and Co incorporate into overall "outer" oxide layer
#             Section.OuterOxThickness[i] = (Section.OuterFe3O4Thickness[i]
#                 + Section.CoThickness[i]
#                 + Section.NiThickness[i])  
            if ElementTracking == "yes":
                Section.OuterOxThickness[i] = (Section.OuterFe3O4Thickness[i]
                + Section.CoThickness[i]
                + Section.NiThickness[i])
            else:
                Section.OuterOxThickness[i] = Section.OuterFe3O4Thickness[i]
            
            Section.InnerOxThickness[i] = Section.InnerIronOxThickness[i]
        else:  # OuterFe3O4Thickness == 0
#             Section.InnerOxThickness[i] = (Section.InnerIronOxThickness[i]
#             + Section.CoThickness[i]
#             + Section.NiThickness[i])
            
            if ElementTracking == "yes":
                 Section.InnerOxThickness[i] = (Section.InnerIronOxThickness[i]
                + Section.CoThickness[i]
                + Section.NiThickness[i])
            else:
                Section.InnerOxThickness[i] = Section.InnerIronOxThickness[i]

    return ElapsedTime, SpallTime
