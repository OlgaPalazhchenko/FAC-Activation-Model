import lepreau_data as ld
import numpy as np
import thermochemistry_and_constants as nc
import random
import matplotlib.pyplot as plt
import csv


def RIHT_plots():
    RIHT_data = open('RIHTOutputSG2.csv', 'r')
    RIHTReader = list(csv.reader(RIHT_data, delimiter=',')) 
    
    RIHT = [float(RIHTReader[1][i]) for i in range((2016-1983)*4)]
    Year = [float(RIHTReader[2][i]) for i in range((2016-1983)*4)]
   
    ax1 = plt.subplot()
    
    ax1.plot(
        Year, RIHT, linestyle=None, marker='o', color='0.50',
        label='Inner Oxide'
        )
    ax1.set_xlabel('Year')
    ax1.set_ylabel('Reactor Inlet Header Temperature (${^o C}$)')
    plt.tight_layout()
    plt.show()
   
    # converts from string object to time series  
    
    
# RIHT_plots()

SGFastMode = "yes"
Method = "tube length"

FullTubeComplement = 3541

YearStartup = 1983.25
YearCPP = 1988.25
YearOutage = 1995.25
YearOutageRestart = 1996
YearRefurbishment = 2008.25
YearRefurbRestart = 2014

T_sat_primary = 310 + 273.15
T_PreheaterIn = 186.5 + 273.15
RecirculationRatio = 5.3


h_i = nc.SGParameters()
h_o = nc.SGParameters()
U_h = nc.SGParameters()
U_c = nc.SGParameters()
U_total = nc.SGParameters()
R_F_primary = nc.SGParameters()
R_F_secondary = nc.SGParameters()
MassFlow_preheater = nc.SGParameters()
MassFlow_downcomer = nc.SGParameters()
MassFlow_c_total = nc.SGParameters()
ShellDiameter = nc.SGParameters()
MassFlux_c = nc.SGParameters()
MassFlux_h = nc.SGParameters()
EnthalpySaturatedSteam = nc.SGParameters()
MassFlow_h = nc.SGParameters()
MasssFlow_dividerplate = nc.SGParameters()
EquivalentDiameter = nc.SGParameters()
TubePitch = nc.SGParameters()

for i in [h_i, h_o]:
    i.unit = "K/W"

for i in [U_h, U_c, U_total]:
    i.unit = "W/cm^2 K"

for i in [R_F_primary, R_F_secondary]:
    i.unit = "cm^2 K/W"

MassFlow_h.magnitude = 1900 * 1000
# Steam flow for 4 steam generators in typical CANDU-6 = 1033.0 kg/s
# 240 kg/s pulled from AECL COG document and works well with the 1900 kg/s hot-side flow ?
MassFlow_preheater.magnitude = 239 * 1000
MassFlow_ReheaterDrains = 0#72 * 1000 / 4
MassFlow_downcomer.magnitude = (RecirculationRatio - 1) * MassFlow_preheater.magnitude
MassFlow_c_total.magnitude = MassFlow_downcomer.magnitude + MassFlow_preheater.magnitude + MassFlow_ReheaterDrains


ShellDiameter.magnitude = 2.28 * 100

for i in [MassFlow_c_total, MassFlow_h, MasssFlow_dividerplate, MassFlow_preheater, MassFlow_downcomer]:
    i.unit = "g/s"

for i in [ShellDiameter, EquivalentDiameter, TubePitch]:
    i.unit = "cm"

for i in [MassFlux_c, MassFlux_h]:
    i.unit = "g/cm^2 s"
    
EnthalpySaturatedSteam.magnitude = 2116.3  # 9.89 MPa
EnthalpySaturatedSteam.unit = "J/g"

TubePitch.magnitude = 2.413

# select desired SG tubes to be run by arc length
UBends = [1.49, 0.685, 2.31, 3.09]
UBends = [i * 100 for i in UBends]
# select desired tubes to be run by total tube length
TubeLengths = [1887, 1807, 1970, 2046]


def total_tubes_plugged(Bundle, CalendarYear):
    #total tubes called sometimes by bundle inside an SG or for a specific boiler itself
    
    if Bundle in ld.SteamGenerator:
        SteamGenerator = ld.SteamGenerator
        YearPlugged = [YearStartup, 1988, 1992, 1993, 1994, 1996, 1999, 2002, 2009]
        AmountPlugged = [1, 1, 10, 6, 6, 1, 8, 1, 2]
    
    
    elif Bundle in ld.SteamGenerator_2:
        SteamGenerator = ld.SteamGenerator_2
        
        YearPlugged = [YearStartup, 1996, 1999]
        AmountPlugged = [6, 5, 11]

    
#     TotalAmountPlugged = []
#     for i in range(len(AmountPlugged)):
#         x = sum(AmountPlugged[0:(i+1)])
#         TotalAmountPlugged.append(x)    
    
    if CalendarYear < YearPlugged[0]:
        NumberPluggedTubes = 0
    
    elif CalendarYear >= YearPlugged[len(YearPlugged) - 1]:
        NumberPluggedTubes = sum(AmountPlugged)
    else:
        for i in range(len(YearPlugged) - 1):

            if YearPlugged[i] <= CalendarYear < YearPlugged[i + 1]:
                NumberPluggedTubes = sum(AmountPlugged[0:(i+1)])
            
 
    TotalSGTubeNumber = 3542 - NumberPluggedTubes
    
    return TotalSGTubeNumber


def bundle_sizes(SteamGenerator, TotalSGTubeNumber):
    if TotalSGTubeNumber < 3542:
        for i in range(3542 - TotalSGTubeNumber):
            x = random.randint(0, 86)
            SteamGenerator[x].TubeNumber = SteamGenerator[x].TubeNumber - 1
            
    
def primaryside_cleaned_tubes(Bundle, CalendarYear):
    #amount of tubes cleaned per each cleaning will have to be custom
    if CalendarYear == YearOutageRestart:
        # siva blast in 1995 used on only 60% of tubes due to time/spacial constraints
        PercentTubesCleaned = 0.6
    elif CalendarYear == YearRefurbRestart:
        PercentTubesCleaned = 1
    else:
        # none cleaned outside of outages, cleaned tubes list below has no appended bundles
        PercentTubesCleaned = 0
    
    TotalSGTubeNumber = total_tubes_plugged(Bundle, CalendarYear)
    
    # chooses tube bundles until % of total sg tube number reached, adds all chosen to "cleaned" tube list
    Cleaned = []
    NumberTubes = []
    
    for i in range(len(ld.SteamGenerator)):
        x = random.randint(0, 86)
        NumberTubes.append(ld.SteamGenerator[x].TubeNumber)
        
        if sum(NumberTubes) <= (PercentTubesCleaned * TotalSGTubeNumber):
            Cleaned.append(ld.SteamGenerator[x])
        else:
            break
    
    return Cleaned, TotalSGTubeNumber


def closest_tubelength(TubeLength):
    # input = desired tube length
    # function scans through all available 87 lengths and outputs the number of the tube that has the closest length
    
    difference = []
    for Zone in ld.SteamGenerator:
        difference.append(abs(TubeLength - Zone.Distance[Zone.NodeNumber - 1]))
    
    return difference.index(min(difference))


def closest_ubend(UbendLength):
    # searches through all 87 u-bend arc lengths and chooses the one closest to input
    
    difference = []
    for i in ld.u_bend_total:
        # calculates differences between input Number and all others in given list
        difference.append(abs(UbendLength - i))
    
    # returns index of value that has smallest difference with input Number
    return difference.index(min(difference))


def tube_picker(Method, SteamGenerator):
    tubes = []
    tube_number = [] #number (from 0 to 86) of the tube bundle
    
    # selects class initializations (for desired steam generator) based on desired u-bend tube arc lengths
    if Method == "arc length": 
        for k in UBends:  # want multiple u-bends
            x = closest_ubend(k)
            tube_number.append(x)

        for i in tube_number:
            tubes.append(SteamGenerator[i])
    # selects class initializations based on desired overall tube lengths        
    elif Method == "tube length":
        #if different tubes needed to be run in each steam generator, customize line below
        if SteamGenerator == ld.SteamGenerator or SteamGenerator == ld.SteamGenerator_2:
            Lengths = TubeLengths

        for k in Lengths:
            x = closest_tubelength(k)
            tube_number.append(x)
        
        for i in tube_number:
            tubes.append(SteamGenerator[i])
    
    else:
        None
        
    return tubes, tube_number


Default_Tube = closest_ubend(1.52 * 100)



def MassFlux_c(Bundle, i):
    if Bundle.Length.label[i] == "opposite preheater":
        PreheaterDiameter = 1.3 * 100  # [cm]
        MassFlow = (RecirculationRatio - 1) * MassFlow_preheater.magnitude
        TotalTubes = FullTubeComplement
    else:
        PreheaterDiameter = 0
        MassFlow = MassFlow_c_total.magnitude
        TotalTubes = FullTubeComplement  # * 2
        
    ShellCrossBundlealArea = (np.pi / 4) * (
    (ShellDiameter.magnitude ** 2)
    - (ld.Bundle[0].OuterDiameter[0] ** 2) * TotalTubes
    - PreheaterDiameter ** 2
    )
    Flux = MassFlow / ShellCrossSectionalArea
    return Flux


def thermal_conductivity(Twall, material, SecondarySidePressure):
    if material == "Alloy-800" or material == "Alloy800" or material == "A800" or material == "Alloy 800":
        Twall_C = Twall - 273.15  # oC in alloy thermal conductivity equation 
        conductivity = (11.450 + 0.0161 * Twall_C) / 100  # M.K.E thesis for Alloy-800 tubing [W/cm K] 
        return conductivity
    
    elif material == "inner magnetite":
        return 2 / 100  # [W/cm K]
    
    elif material == "outer magnetite":
        return 1.3 / 100  # [W/cm K]
    else:
        return None


def sludge_fouling_resistance(Bundle, i, calendar_year):
    
    TubeGrowth = 0.00135  # [g/cm^2] /year = 6.5 um /year
    ReducedTubeGrowth = 0.0001  # [g/cm^2] /year = 3.25 um/year
    
    InitialAccumulation = 0.0001 #g/cm^2
      
    # between 1983 and 1986 (included)
    if YearStartup <= calendar_year <= YearCPP:
        # secondary side fouling slope based on 1/2 that for average primary side cold leg deposit
        Accumulation = InitialAccumulation + (calendar_year - YearStartup) * TubeGrowth
     
    # CPP installation (late 1986) reduces secondary side crud by 50% 
    # between 1987 and 1995, not including, maybe some form of dissolution event to historic deposits
    elif YearCPP < calendar_year < YearOutageRestart:
        Accumulation = InitialAccumulation + (calendar_year - YearCPP) * ReducedTubeGrowth
    
    elif calendar_year >= YearOutageRestart:
        Accumulation = (calendar_year - YearOutageRestart) * ReducedTubeGrowth
        
    elif YearRefurbishment < calendar_year < YearRefurbRestart:
        #no more SHT fouling growth during refurbishment outage - stays at this while those years are run blank
        Accumulation = (YearRefurbishment - YearOutageRestart) * ReducedTubeGrowth
    #
    #
    #
    #
    #
    # did all of SHT get cleaned?
    elif calendar_year >= YearRefurbRestart:
        Accumulation = (calendar_year - YearRefurbRestart) * ReducedTubeGrowth
    
    Thickness = Accumulation / nc.Fe3O4Density
    
    Fouling = Thickness / thermal_conductivity(None, "outer magnetite", None)
    
    return Fouling
    

def pht_fouling_resistance(i, calendar_year, InnerAccumulation, OuterAccumulation):
    
    # [g/cm^2]/[g/cm^3] = [cm]
    # thickness/thermal conductivity [cm]/[W/cm K] = [cm^2 K/W]
    InnerThickness = InnerAccumulation / nc.Fe3O4Density
    OuterThickness = OuterAccumulation / nc.Fe3O4Density 

    # [cm]/ [W/cm K] =[cm^2 K/W]
    # inner deposit is consolidated, predicted to have different thermal resistance than unconsolidated outer ox.
    InnerFouling = InnerThickness / thermal_conductivity(None, "inner magnetite", None)
    OuterFouling = OuterThickness / thermal_conductivity(None, "outer magnetite", None)
    # total pht thermal resistance
    
    return InnerFouling + OuterFouling  # [cm^2 K/W]


def conduction_resistance(Bundle, Twall, SecondarySidePressure, i):
    # R_conduction = L/kA (L = thickness)
    # A = area with shape factor correcrion factor for cylindrical pipe
    # R_conduction = ln(D_o/D_i)/(2pikl) [K/W]
    # l = length, k =thermal conductivity coefficient [W/cm K]
    
    k_w = thermal_conductivity(Twall, "Alloy-800", SecondarySidePressure)  # [W/cm K]

    # small conductivity = big resistance
    Rcyl_numerator = np.log(Bundle.OuterDiameter[i] / Bundle.Diameter[i])
    Rcyl_denominator = 2 * np.pi * k_w * Bundle.Length.magnitude[i]

    R_cond = Rcyl_numerator / Rcyl_denominator  # [K/W]
    return R_cond


def primary_convection_resistance(
        Bundle, T_PrimaryBulkIn, T_PrimaryWall, x_pht, m_h_leakagecorrection, i, CalendarYear
        ):
    
    TotalSGTubeNumber = total_tubes_plugged(Bundle, CalendarYear)
    
    # R_1,Convective = 1/(h_1,convectove *A)
    # A = inner area (based on inner diameter)
    # [W/cm K]/ [cm] = [W/K cm^2]
    
    Diameter = Bundle.Diameter[i]
    # T_PrimaryBulkIn used in place of Tavg (of in and out for that node)
    
    MassFlux_h.magnitude = m_h_leakagecorrection / (TotalSGTubeNumber * (np.pi / 4) * (Diameter ** 2))
   
    Nu_D = nusseltnumber(Bundle, "PHT", T_PrimaryBulkIn, MassFlux_h.magnitude, i)
    
    h_l =  Nu_D * nc.thermal_conductivityD2O_liquid(T_PrimaryBulkIn) / Bundle.Diameter[i]
    
     # quality only persists for first 2 nodes of PHT
#     if x_pht > 0:
#         h_i = in_tube_boiling_LiuWinterton(x_pht, MassFlux_h.magnitude, Diameter, T_PrimaryWall, h_l, i)
#     else:
    h_i = h_l
    
    return 1 / (h_i * inner_area(Bundle)[i])  # [K/W]


def secondary_convection_resistance(
        Bundle, T_SecondaryBulkIn, T_SecondaryWall, T_SecondaryFilm, SecondarySidePressure, x_sht, i
        ):
    
#     print (T_SecondaryBulkIn-273.15, T_SecondaryWall-273.15, T_SecondaryFilm-273.15, SecondarySidePressure, x_sht,i)
    # split into boiling and non-boiling (preheater) Bundles
    if Bundle.Length.label[i] == "preheater":  # from Silpsrikul thesis
        
        # Based on PLGS data and Silpsrikul calculations
        f_b = 0.2119  # fraction of cross-Bundleal area of shell occupied by a baffle window
        N_b = f_b * FullTubeComplement  # number of tubes in baffle window
        D_s = ((4 * 0.5 * np.pi / 4) * (ShellDiameter.magnitude / 100) ** 2) \
            / (0.5 * np.pi * (ShellDiameter.magnitude / 100) + ShellDiameter.magnitude / 100)  # [m]

        S_b = (f_b * np.pi * (D_s ** 2) / 4) - N_b * np.pi * ((Bundle.OuterDiameter[i] / 100) ** 2) / 4  # [cm^2]
        G_b = (MassFlow_preheater.magnitude / 1000) / S_b  # [kg/m^2s]
        S_c = 0.5373 * D_s * (1 - (Bundle.OuterDiameter[i] / 100) / 0.0219)
        G_c = MassFlow_preheater.magnitude / 1000 / S_c
        G_e = np.sqrt(G_b * G_c) * 1000 / (100 ** 2)  # [g/cm^2]

        # First two (raised to exponent) terms are unitless
        # [W/cm K]/[cm] = [W/cm^2 K]
        k = nc.thermal_conductivityH2O_liquid(T_SecondaryFilm, SecondarySidePressure)
#         k = thermal_conductivity(T_SecondaryFilm, "water", SecondarySidePressure)
        h_o = (k * 0.2 / Bundle.OuterDiameter[i]) \
            * ((Bundle.OuterDiameter[i] * G_e / nc.viscosityH2O_liquid(T_SecondaryFilm, SecondarySidePressure)) ** 0.6) \
            * (
                nc.heatcapacityH2O_liquid(T_SecondaryFilm, SecondarySidePressure) \
                * nc.viscosityH2O_liquid(T_SecondaryFilm, SecondarySidePressure) / k
                ) ** 0.33
    
    # outside preheater, no cross-flow/baffles, but two-phase flow
    else:
        T_sat_secondary = nc.saturation_temperatureH2O(SecondarySidePressure)
        # [W/cm^2 K]
        h_o = outside_bundle_pool_boiling(
            "FZ", Bundle, x_sht, T_sat_secondary, T_SecondaryWall, T_SecondaryBulkIn, SecondarySidePressure, i
            )
#         print (h_o, i)
#         h_1 = outside_bundle_pool_boiling(
#             "lol", Bundle, x_sht, T_sat_secondary, T_SecondaryWall, T_SecondaryBulkIn, SecondarySidePressure, i
#             )
#         h_2 = outside_bundle_pool_boiling(
#             "None", Bundle, x_sht, T_sat_secondary, T_SecondaryWall, T_SecondaryBulkIn, SecondarySidePressure, i
#             )
        
#         print (h_o, h_1, h_2, i)
        # outside_bundle_pool_boiling(None, T_sat_secondary, T_SecondaryWall, SecondarySidePressure)
        
        # Zukasuskas correlation-based heat transfer coefficient for a bundle, non boiling, not in pre-heater
        # [cm]
        
#         EquivalentDiameter.magnitude = 4 * ((TubePitch.magnitude ** 2) \
#                                             - (np.pi / 4) * (Bundle.OuterDiameter[i] ** 2))\
#                                             / (np.pi * Bundle.OuterDiameter[i])
#          
#         Nu_D_l = nusseltnumber(
#             Bundle, "SHT", T_PrimaryBulkIn, T_SecondaryBulkIn, T_SecondaryWall, SecondarySidePressure, i, None,
#             None
#             )
#      
#         h_l = nc.thermal_conductivityH2O("SHT", T_SecondaryBulkIn, SecondarySidePressure) * Nu_D_l \
#         / EquivalentDiameter.magnitude
#          
#         if x_in > 0:
#             MassFlux_c.magnitude = MassFlux_c(Bundle, i)
#              
#             h_o = boiling_heat_transfer(
#                 x_in, "SHT", MassFlux_c.magnitude, T_SecondaryWall, Bundle.OuterDiameter[i],
#                 SecondarySidePressure, i, h_l)
#         
#         else:
#             h_o = h_l
    
    return 1 / (h_o * outer_area(Bundle)[i])  # [K/W]


# def PHT_pressure_drop(m_h_leakagecorrection):
#     Diameter = SteamGenerator.Diameter[i]
#     # T_PrimaryBulkIn used in place of Tavg (of in and out for that node)
#     
#     MassFlux_h.magnitude = m_h_leakagecorrection / (TotalSGTubeNumber * (np.pi / 4) * (Diameter ** 2))
#     return None
    

def nusseltnumber(Bundle, side, T_PrimaryBulkIn, MassFlux, i):
   
    if side == "PHT":
        #Dittus-boelter evaluated at average T between inlet and outlet for Bundle - approximated as T_in here
        Re_D = [(MassFlux / nc.viscosityD2O_liquid(T_PrimaryBulkIn)) * i for i in Bundle.Diameter]
        
        n = 1 / 3
        C = 0.023
        m = 4 / 5
        
        Nu_D = C * (Re_D[i] ** m) * (prandtl(side, T_PrimaryBulkIn, SecondarySidePressure = None)) ** n

#     elif side == "SHT":
#         # [g/s] / [g/cm^3] = [cm^3/s] / [cm^2]
#         A_Cross = (np.pi / 4) * ((ShellDiameter.magnitude ** 2) - (TotalSGTubeNumber * Bundle.OuterDiameter[i] ** 2))
#         
#         Velocity = (
#             MassFlow_c_total.magnitude / (nc.densityH2O_liquid(T_SecondaryBulkIn, SecondarySidePressure) * A_Cross))
#         
#         V_max = (TubePitch.magnitude / (TubePitch.magnitude - Bundle.Diameter[i])) * Velocity
#         
#         Re_D_max = V_max * Bundle.Diameter[i] * nc.densityH2O_liquid(T_SecondaryBulkIn, SecondarySidePressure) \
#         / nc.viscosityH2O_liquid(T_SecondaryBulkIn, SecondarySidePressure)
#         
#         C1 = 0.40
#         m = 0.60
#         n = 0.36
#         C2 = 0.955 # lower value = more supression = higher RIHT growth
#         
#         Prandtl_ratio = prandtl(side, T_SecondaryBulkIn, SecondarySidePressure) \
#         / prandtl(side, T_SecondaryWall, SecondarySidePressure)
#         
#         Nu_D = C2 * C1 * (Re_D_max ** m) * ((prandtl(side, T_SecondaryBulkIn, SecondarySidePressure)) ** n) \
#         * Prandtl_ratio ** 0.25
#         
#     else:
#         None
    
    return Nu_D


def prandtl(side, Temperature, SecondarySidePressure):
    
    if side == "PHT":
        Viscosity = nc.viscosityD2O_liquid(Temperature)
        HeatCapacity = nc.heatcapacityD2O_liquid(Temperature)
        ThermalConductivity = nc.thermal_conductivityD2O_liquid(Temperature)
    
    elif side == "SHT":
        Viscosity = nc.viscosityH2O_liquid(Temperature, SecondarySidePressure)
        HeatCapacity = nc.heatcapacityH2O_liquid(Temperature, SecondarySidePressure)
        ThermalConductivity = nc.thermal_conductivityH2O_liquid(Temperature, SecondarySidePressure)
    
    else:
        None
        
    Pr = HeatCapacity * Viscosity / ThermalConductivity
    
    return Pr


def inner_area(SteamGenerator):
    return [np.pi * x * y for x, y in zip(SteamGenerator.Diameter, SteamGenerator.Length.magnitude)]
#     return [np.pi * (x - z) * y for x, y, z in zip(SteamGenerator.Diameter, SteamGenerator.Length.magnitude, Oxide)]


def outer_area(SteamGenerator):
    return [np.pi * x * y for x, y in zip(SteamGenerator.OuterDiameter, SteamGenerator.Length.magnitude)]


def reactor_power(CalendarYear):
    from matplotlib import pyplot as plt
    import pandas as pd
    
    
    if CalendarYear == YearOutage:
        Filename = 'C:\\Users\\opalazhc\\git\\FAC-Activation-Model\\FACModel\\OutagePower.csv'
    else:
        FileName = Filename = 'C:\\Users\\opalazhc\\git\\FAC-Activation-Model\\FACModel\\RefurbPower.csv'
    
    df = pd.read_csv(Filename, header=None)
    df.columns = ['date','power']

    # converts from string object to time series  
    df.date = pd.to_datetime(df.date)
    df.set_index('date', inplace=True)

#     plt.plot(df)
#     plt.xlabel('Year')
#     plt.ylabel('Power')
#     plt.show()
##     prints first 5 rows to check everything is assigned properly
#     print (df.head())

    df['year'] = [df.index[i].year for i in range(len(df))]
    df['quarter'] = [df.index[i].quarter for i in range(len(df))]
    # df['day'] = [df.index[i].day for i in range(len(df))]
    
    quarterly_grouping = df.groupby(["year","quarter"])
    
    year_output = df.as_matrix(['year'])
    quarter_output = df.as_matrix(['quarter'])
    
    #first and last year and quarter of averaged power data
    initial_quart = quarter_output[0]
    initial_yr = year_output[0]
  
    final_quart = quarter_output[len(quarter_output) - 1]
    final_yr = year_output[len(year_output) - 1]
   
    # year/mo converted to year.quarter (i.e., 1999.75) based on which quarter first data point is in
    if initial_quart == 1:
        year_start = initial_yr
    elif initial_quart == 2:
        year_start = initial_yr + 0.25
    elif initial_quart == 3:
        year_start = initial_yr + 0.5
    else:
        year_start = initial_yr + 0.75
    
    
    if final_quart == 1:
        year_end = final_yr
       
    elif final_quart == 2:
        year_end = final_yr + 0.25
    elif final_quart == 3:
        year_end = final_yr + 0.5
    else:
        year_end = final_yr + 0.75
    
    year_end = year_end + 0.25 # arange function does not include last element
    YearDerates = np.arange(year_start, year_end, 0.25)
     
    # For each group, calculate the average of only the Power column
    quarterly_average = quarterly_grouping.aggregate({"power":np.mean})
    
    power_output = quarterly_average.as_matrix(['power'])
    # converts matrix notation to 1D array
    power_output = power_output.flatten()

    return power_output, YearDerates


OutagePower, OutageYearDerate = reactor_power(YearOutage)
RefurbPower, RefurbYearDerate = reactor_power(YearRefurbishment)


def pht_steam_quality(Temperature, j):
    
    CalendarYear = (j * nc.TIME_STEP / 8760) + YearStartup
    
    if CalendarYear in OutageYearDerate:
        PercentDerates = OutagePower
        YearDerates = OutageYearDerate
    
    elif CalendarYear in RefurbYearDerate:
        PercentDerates = RefurbPower
        YearDerates = RefurbYearDerate
    
    else:
        PercentDerates = [0]
        YearDerates = [0]
         
    CoreMassFlow =  7600  # [kg /s]
#     Delta_T = T_sat_primary - (262 + 273.15)  # [K]
#     C_p_cold = nc.heatcapacityD2O_liquid(262 + 273.15)
#     C_p_hot = nc.heatcapacityD2O_liquid(T_sat_primary)
#     C_p_avg = (C_p_cold + C_p_hot) / 2  # [kJ/kg K]

    Power = 2061.4 * 1000 #CoreMassFlow * C_p_avg * Delta_T  # [kW]
#     print (CalendarYear, PercentDerates)
    difference = []
    if CalendarYear in YearDerates:
        for Yr in YearDerates:
            difference.append(abs(CalendarYear - Yr))
            ClosestYear = difference.index(min(difference))
        Percent_derating = PercentDerates[ClosestYear]
        print (Percent_derating)
    else:
        Percent_derating = 1
       
    P = Power * Percent_derating

    # core inlet enthalpy at RIHT + that added from fuel
    # H_fromfuel = Power / CoreMassFlow # [kJ/s /kg/s] = [kJ/kg]
    H_pht = P / CoreMassFlow

    H_satliq_outlet = nc.enthalpyD2O_liquid(T_sat_primary)
#     H_satliq_outlet = nc.enthalpyH2O_liquid(nc.PrimarySidePressure, T_sat_primary)
    
    # no quality in return primary flow
    H_current = nc.enthalpyD2O_liquid(Temperature) + H_pht
#     H_current = nc.enthalpyH2O_liquid(nc.PrimarySidePressure, Temperature) + H_pht
    x = (H_current - H_satliq_outlet) / (EnthalpySaturatedSteam.magnitude - H_satliq_outlet)
    
#     print (CalendarYear, P, Percent_derating)
    if x < 0:
        x = 0
    return x

# print (pht_steam_quality(264 +273.15, (25.5) * 876))


def sht_steam_quality(Q, T_sat_secondary, x, MassFlow_c, SecondarySidePressure):
        
    # [J/s] / [g/s] = [J /g]
    H_sht = Q / (MassFlow_c)  # [kJ/kg]
    H_satliq = nc.enthalpyH2O_liquid(SecondarySidePressure, T_sat_secondary)
    H_SaturatedSteam = nc.enthalpyH2O_vapour(SecondarySidePressure, T_sat_secondary)  # [kJ/kg]

    H_prev = H_satliq + x * (H_SaturatedSteam - H_satliq)
    H_current = H_sht + H_prev

    x = (H_current - H_satliq) / (H_SaturatedSteam - H_satliq)
    return x


def Mostinski_outside_tube_boiling(T_sat_secondary, T_SecondaryWall, SecondarySidePressure):
    
    P_crit = 22.064 * 10 **3 #kPa
    
    if T_SecondaryWall <= T_sat_secondary:
        T_SecondaryWall = T_sat_secondary + 0.5
    
    T_e = T_SecondaryWall - T_sat_secondary # K
    P_r = (SecondarySidePressure * 10 **3) / P_crit
    
    F_p = 1.8 * (P_r ** 0.17) + 4 * (P_r ** 1.2) + 10 * (P_r ** 10)
    
    h_nb = (1.167 * 10 ** (-8)) * (P_crit ** 2.3) * (T_e ** 2.333) * F_p ** 3.333
    
    h_nb = h_nb / (100 ** 2) # [W/cm^2 K]
    
    return h_nb


def Zukauskas_outside_tube_boiling(
        Bundle, side, T_sat_secondary, T_SecondaryWall, T_SecondaryBulkIn, SecondarySidePressure, i
        ):
    
    A_Cross = (np.pi / 4) * (
        (ShellDiameter.magnitude ** 2) - (TotalSGTubeNumber * Bundle.OuterDiameter[i] ** 2)
        )

    Density_l = nc.densityH2O_liquid(T_sat_secondary, SecondarySidePressure)
    Velocity = MassFlow_c_total.magnitude / (Density_l * A_Cross)
    
    V_max = (TubePitch.magnitude / (TubePitch.magnitude - Bundle.Diameter[i])) * Velocity
    
    Re_D_max = V_max * Bundle.OuterDiameter[i] * Density_l \
    / nc.viscosityH2O_liquid(T_sat_secondary, SecondarySidePressure)
    
    C1 = 0.40
    m = 0.60
    n = 0.36
    C2 = .75 # lower value = more supression = higher RIHT growth
    
    Prandtl_ratio = prandtl(side, T_SecondaryBulkIn, SecondarySidePressure) \
    / prandtl(side, T_SecondaryWall, SecondarySidePressure)
    
    Nu_D_l = C2 * C1 * (Re_D_max ** m) * ((prandtl(side, T_SecondaryBulkIn, SecondarySidePressure)) ** n) \
    * Prandtl_ratio ** 0.25
    
    OD = Bundle.OuterDiameter[i]
    
    EquivalentDiameter.magnitude = 4 * ((TubePitch.magnitude ** 2) - (np.pi / 4) * (OD ** 2)) / (np.pi * OD)
    
    ED = EquivalentDiameter.magnitude
    h_l = nc.thermal_conductivityH2O_liquid(T_SecondaryBulkIn, SecondarySidePressure) * Nu_D_l / ED
    
    return h_l
    

# def LiuWinterton_outside_bundle(
#         Bundle, x_sht, T_sat_secondary, T_SecondaryWall, T_SecondaryBulkIn, SecondarySidePressure, i
#         ):
#     
#     h_l = Zukauskas_outside_tube_boiling(
#         Bundle, "SHT", T_sat_secondary, T_SecondaryWall, T_SecondaryBulkIn, SecondarySidePressure, i
#         )
#     
#     h_l = h_l * (100 ** 2) 
#     Pressure = SecondarySidePressure
#     T_sat_secondary = nc.saturation_temperatureH2O(Pressure)
#     rho_v = nc.densityH2O_vapour(Pressure, T_sat_secondary)  # [g/cm^3]
#     rho_l = nc.densityH2O_liquid(T_sat_secondary, SecondarySidePressure)
#     Viscosity_sat = nc.viscosityH2O_liquid(T_sat_secondary, SecondarySidePressure)
#      
#     F = (1 + x_sht * prandtl("SHT", T_sat_secondary, Pressure) * ((rho_l / rho_v) - 1)) ** 0.35
#    
#     MassFlux_c.magnitude = MassFlux_c(Bundle, i)
#     MassFlux_liquid = MassFlux_c.magnitude * (1 - x_sht)
#      
#     Re_D = Bundle.OuterDiameter[i] * MassFlux_liquid / Viscosity_sat
#      
#     S = (1 + 0.055 * (F ** 0.1) * (Re_D) ** 0.16) ** (-1)
#     
#     p_crit = 22.0640  # [MPa]
#     P_r = Pressure / p_crit   
#    
#     q_l = F * h_l * (abs(T_SecondaryWall - T_sat_primary))  # [K W/cm^2]
#     
#     A_p = 55 * (P_r ** 0.12) * ((-np.log10(P_r)) ** (-0.55)) * (nc.H2OMolarMass) ** (-0.5) 
#     
#     C = ((A_p * S / (F * h_l)) ** 2) * q_l ** (4 / 3)
#     coeff = [1, -C, 0, -1]
#     cubic_solution = np.roots(coeff)
#     
#     # only real roots usable in correlation
# #     q = []
# #     for x in cubic_solution:
# #         if np.isreal(x):
# #             q = x
#     q = cubic_solution[0]
#     h_twophase =  F * (q ** (3 / 2)) * h_l / (100 ** 2)  # [W/cm^2 K]
#     return h_twophase

 
def Cooper_outside_tube_boiling(T_sat_secondary, T_SecondaryWall, SecondarySidePressure):
    P_crit = 22.064 # MPa
    
    if T_SecondaryWall <= T_sat_secondary:
        T_SecondaryWall = T_sat_secondary + 0.5
        
    deltaT_sat = T_SecondaryWall - T_sat_secondary # K
    
    P_r = SecondarySidePressure / P_crit
     
    h = (55 * P_r ** 0.12) * ((-np.log10(P_r)) ** (-0.55)) * (nc.H2OMolarMass ** (-0.5)) * (deltaT_sat ** 0.67)
     
    h_nb = (h  ** (1 / 0.33)) / (100 ** 2)
    
    return h_nb     


# def ForsterZuber_outside_tube_boiling(T_sat_secondary, T_SecondaryWall, SecondarySidePressure):
# 
#     # W/cm K / 1000 -- > kW/cm K * 100 --> kW/m K thermal conductivity
#     lambda_l = (nc.thermal_conductivityH2O_liquid(T_sat_secondary, SecondarySidePressure) / 10) ** 0.79 
#     #J / g K --> liquid specific heat
#     cp_l = (nc.heatcapacityH2O_liquid(T_sat_secondary, SecondarySidePressure)) ** 0.45 
#     rho_l = (nc.densityH2O_liquid(T_sat_secondary, SecondarySidePressure) * (100 ** 3) / 1000) ** 0.49 # kg/m^3
#     rho_v = (nc.densityH2O_vapour(SecondarySidePressure, T_sat_secondary) * (100 ** 3) / 1000) ** 0.24 # kg/m^3
#     sigma = (0.024) ** 0.5 #surface tension N/m
#     mu = (nc.viscosityH2O_liquid(T_sat_secondary, SecondarySidePressure) / 10) ** 0.29 # liquid viscosity #kg/m s 
#     
#     enthalpy_l = nc.enthalpyH2O_liquid(SecondarySidePressure, T_sat_secondary)  # J / g (kJ/kg)
#     enthalpy_v = nc.enthalpyH2O_vapour(SecondarySidePressure, T_sat_secondary) 
#     h_lg = (enthalpy_v - enthalpy_l) ** 0.24 
#     
#     if T_SecondaryWall <= T_sat_secondary:
#         T_SecondaryWall = T_sat_secondary + 0.5
#         
#     deltaT_sat = (T_SecondaryWall - T_sat_secondary) ** 0.24 # K
#     P_w = nc.saturation_pressureH2O(T_SecondaryWall) * 10 ** 6
#     P_s = nc.saturation_pressureH2O(T_sat_secondary) * 10 ** 6
#     
#     deltaP_sat = (P_w - P_s) ** 0.75 # MPa
#         
#     h = 0.00122 * lambda_l * cp_l * rho_l * deltaT_sat * deltaP_sat / (sigma * mu * h_lg * rho_v) 
#     h_nb = 1000 * h / (100 ** 2) # [W/cm^2 K]
#     print (h, h_nb)
# 
#     return h_nb


def ForsterZuber_outside_tube_boiling(T_sat_secondary, T_SecondaryWall, SecondarySidePressure):
 
    # W/cm K  * 100 --> W/m K thermal conductivity
    lambda_l = (nc.thermal_conductivityH2O_liquid(T_sat_secondary, SecondarySidePressure) * 100) ** 0.79 
    #J / g K --> liquid specific heat
    cp_l = (nc.heatcapacityH2O_liquid(T_sat_secondary, SecondarySidePressure) * 1000) ** 0.45 
    rho_l = (nc.densityH2O_liquid(T_sat_secondary, SecondarySidePressure) * (100 ** 3) / 1000) ** 0.49 # kg/m^3
    rho_v = (nc.densityH2O_vapour(SecondarySidePressure, T_sat_secondary) * (100 ** 3) / 1000) ** 0.24 # kg/m^3
    sigma = (0.024) ** 0.5 #surface tension N/m
    mu = (nc.viscosityH2O_liquid(T_sat_secondary, SecondarySidePressure) / 10) ** 0.29 # liquid viscosity #kg/m s 
     
    enthalpy_l = nc.enthalpyH2O_liquid(SecondarySidePressure, T_sat_secondary) * 1000  # J / g (kJ/kg)
    enthalpy_v = nc.enthalpyH2O_vapour(SecondarySidePressure, T_sat_secondary) * 1000
    h_lg = (enthalpy_v - enthalpy_l) ** 0.24 
     
    if T_SecondaryWall <= T_sat_secondary:
        T_SecondaryWall = T_sat_secondary + 1
         
    deltaT_sat = (T_SecondaryWall - T_sat_secondary) ** 0.24 # K
    P_w = nc.saturation_pressureH2O(T_SecondaryWall) * 10 ** 6
    P_s = nc.saturation_pressureH2O(T_sat_secondary) * 10 ** 6
     
    deltaP_sat = (P_w - P_s) ** 0.75 # MPa
         
    h = 0.00122 * lambda_l * cp_l * rho_l * deltaT_sat * deltaP_sat / (sigma * mu * h_lg * rho_v) # W / m^2 K
     
    h_nb = h / (100 ** 2) # [W/cm^2 K]
#     print (h, h_nb)
 
    return h_nb


def outside_bundle_pool_boiling(
        Correlation, Bundle, x_sht, T_sat_secondary, T_SecondaryWall, T_SecondaryBulkIn, SecondarySidePressure,
        i
        ):
    
    F = .11 # bundle boiling factor (empirical)
    if Correlation == "FZ":
        
        h_nb = ForsterZuber_outside_tube_boiling(T_sat_secondary, T_SecondaryWall, SecondarySidePressure)
        #h_nb = Cooper_outside_tube_boiling(T_sat_secondary, T_SecondaryWall, SecondarySidePressure)
    elif Correlation == "other":
        h_nb = Cooper_outside_tube_boiling(T_sat_secondary, T_SecondaryWall, SecondarySidePressure)

    else:
        h_nb = Mostinski_outside_tube_boiling(T_sat_secondary, T_SecondaryWall, SecondarySidePressure)
        #h_nb = LiuWinterton_outside_bundle(
#             Bundle, x_sht, T_sat_secondary, T_SecondaryWall, T_SecondaryBulkIn,
#                                            SecondarySidePressure, i
#                                            )
        #h_nb = outside_tube_boiling_Cooper(T_sat_secondary, T_SecondaryWall, SecondarySidePressure)
        # outside_tube_boiling_Mostinski(T_sat_secondary, T_SecondaryWall, SecondarySidePressure)

    h_twophase = F * h_nb
#     print (h_twophase, h_nb)
    return h_twophase
    

# click to open

# def in_tube_boiling_LiuWinterton(x_pht, MassFlux, Diameter, T_PrimaryWall, h_l, i):
#     h_l = h_l * (100 ** 2)
#     
#     Pressure = nc.PrimarySidePressure
#     Density_v = nc.densityD2O_vapour(Pressure, T_sat_primary)  # [g/cm^3]
#     Density_l = nc.densityD2O_liquid(T_sat_primary)
#     Viscosity_sat = nc.D2O_viscosity(T_sat_primary)
#     
#     F = (1 + (x_pht * prandtl("PHT", T_sat_primary, Pressure) * ((Density_l / Density_v) - 1))) ** 0.35
#     
#     MassFlux_liquid = MassFlux * (1 - x_pht)
#     Re_D = Diameter * MassFlux_liquid / Viscosity_sat
#     S = (1 + 0.055 * (F ** 0.1) * (Re_D ** 0.16)) ** (-1)
#     
#     p_crit = 22.0640  # [MPa]
#     P_r = Pressure / p_crit
#     
#     deltaT_sat = 3 # T_PrimaryWall - T_sat_primary
#     q_l = F * h_l * (abs(deltaT_sat))  # [K W/m^2]
#     A_p = 55 * (P_r ** 0.12) * ((-np.log10(P_r)) ** (-0.55)) * (nc.H2OMolarMass) ** (-0.5)
#     
#     C = ((A_p * S / (F * h_l)) ** 2) * q_l ** (4 / 3)
#     
#     coeff = [1, -C, 0, -1]
#     cubic_solution = np.roots(coeff)
# 
#     # only real roots usable in correlation
#     q = []
#     for x in cubic_solution:
#         if np.isreal(x):
#             q.append(x)
#         else:
#             q.append(1)
#     h_twophase =  F * (q[0] ** (3 / 2)) * h_l / (100 ** 2)  # [W/m^2 K]
# 
#     return h_twophase


def wall_temperature(
        Bundle, i, T_PrimaryBulkIn, T_SecondaryBulkIn, x_in, x_pht, InnerAccumulation, OuterAccumulation,
        CalendarYear, SecondarySidePressure, m_h_leakagecorrection):
     
    # i = each node of SG tube
    T_PrimaryWall = T_PrimaryBulkIn - (1 / 3) * (T_PrimaryBulkIn - T_SecondaryBulkIn)
    T_SecondaryWall = T_PrimaryBulkIn - (2 / 3) * (T_PrimaryBulkIn - T_SecondaryBulkIn)
    
    A_i = inner_area(Bundle)[i]
    A_o = outer_area(Bundle)[i]
    
    for k in range(50):
        
        WT_h = T_PrimaryWall
        WT_c = T_SecondaryWall

#         PrimaryT_film = (T_PrimaryBulkIn + T_PrimaryWall) / 2
        T_SecondaryFilm = (T_SecondaryBulkIn + T_SecondaryWall) / 2
               
        h_i.magnitude = (1 / (A_i * primary_convection_resistance(
            Bundle, T_PrimaryBulkIn, T_PrimaryWall, x_pht, m_h_leakagecorrection, i, CalendarYear
            )))
        
        h_o.magnitude = (1 / (A_o * secondary_convection_resistance(
                Bundle, T_SecondaryBulkIn, T_SecondaryWall, T_SecondaryFilm, SecondarySidePressure, x_in, i
            )))
        # R = 1/hA --> h=A*1/R
#         print (A_i, h_i.magnitude, h_o.magnitude, A_o, i) 
        U_h1 = A_i * conduction_resistance(Bundle, T_PrimaryWall, SecondarySidePressure, i)
        
        U_h2 = (A_i * secondary_convection_resistance(
            Bundle, T_SecondaryBulkIn, T_SecondaryWall, T_SecondaryFilm, SecondarySidePressure, x_in, i
            ))
        U_h.magnitude = 1 / (U_h1 + U_h2)  # [W/ cm^2 K]

        U_c1 = A_o * conduction_resistance(Bundle, T_PrimaryWall, SecondarySidePressure, i)
        U_c2 = A_o * primary_convection_resistance(
            Bundle, T_PrimaryBulkIn, T_PrimaryWall, x_pht, m_h_leakagecorrection, i, CalendarYear) 
        
        U_c.magnitude = 1 / (U_c1 + U_c2)  # [W/cm^2 K]

        T_PrimaryWall = (T_PrimaryBulkIn * h_i.magnitude + T_SecondaryBulkIn * U_h.magnitude)\
            / (h_i.magnitude + U_h.magnitude)
        
        T_SecondaryWall = (T_SecondaryBulkIn * h_o.magnitude + T_PrimaryBulkIn * U_c.magnitude)\
            / (h_o.magnitude + U_c.magnitude)

        RE1 = (T_PrimaryWall - WT_h)
        RE2 = (T_SecondaryWall - WT_c)

        # if converged
        if abs(RE1) <= .05 and abs(RE2) <= .05:
#             print (h_i.magnitude, h_o.magnitude, i)
            # [cm^2 K/W]
            R_F_primary.magnitude = pht_fouling_resistance(i, CalendarYear, InnerAccumulation[i], OuterAccumulation[i]) 
            
            R_F_secondary.magnitude = sludge_fouling_resistance(Bundle, i, CalendarYear)
                        
            PCR = primary_convection_resistance(
                Bundle, T_PrimaryBulkIn, T_PrimaryWall, x_pht, m_h_leakagecorrection, i, CalendarYear
                )
            
            SCR = secondary_convection_resistance(
                Bundle, T_SecondaryBulkIn, T_SecondaryWall, T_SecondaryFilm, SecondarySidePressure, x_in, i
               )
            
#             if i ==0:
#                 print (PCR *100*100, "PCR", conduction_resistance(Bundle, T_PrimaryWall, SecondarySidePressure, i)*100*100,\
#                         "cond", SCR*100*100, "SCR", R_F_primary.magnitude*100*100, "primary fouling", \
#                         R_F_secondary.magnitude*100*100, "secondary fouling", i)
              
            # [cm^2 K/W] all resistances (convective, conductive, and fouling)
            
            ConductionResistance = conduction_resistance(Bundle, T_PrimaryWall, SecondarySidePressure, i)
            
            inverseU_total = (PCR + ConductionResistance + SCR) * A_o + R_F_primary.magnitude + R_F_secondary.magnitude

            U_total.magnitude = 1 / inverseU_total  # [W/ cm^2 K]
            
            return T_PrimaryWall, T_SecondaryWall, U_total.magnitude


def temperature_profile(
        Bundle, InnerOxide, OuterOxide, m_h_leakagecorrection, SecondarySidePressure, x_pht, CalendarYear):
    
    TotalSGTubeNumber = total_tubes_plugged(Bundle, CalendarYear)
    
    # bulk temperatures guessed, wall temperatures and overall heat transfer coefficient calculated
    # U assumed to be constant over node's area, bulk temperatures at end of node calculated, repeat
    
    T_sat_secondary = nc.saturation_temperatureH2O(SecondarySidePressure)
    
    PrimaryWall = []
    PrimaryBulk = []
    SecondaryBulk = []
    SecondaryWall = []
    HeatFlux = []

    for i in range(Bundle.NodeNumber):
        if i == 0:
            # Temperatures entering SG (0 m of SG)--> not first "node" temp (several m into SG hot leg)
            T_PrimaryBulkIn = T_sat_primary  # [K]
            T_SecondaryBulkIn = T_sat_secondary
            x_in = 0
        
        Cp_h = nc.heatcapacityD2O_liquid(T_PrimaryBulkIn)
        Cp_c = nc.heatcapacityH2O_liquid(T_SecondaryBulkIn, SecondarySidePressure)
        
        T_wh, T_wc, U = wall_temperature(
            Bundle, i, T_PrimaryBulkIn, T_SecondaryBulkIn, x_in, x_pht, InnerOxide, OuterOxide, CalendarYear,
            SecondarySidePressure, m_h_leakagecorrection)
        
        # [W/ cm^2 K] * [K] * [cm^2] = [W] = [J/s]
        Q = U * (T_PrimaryBulkIn - T_SecondaryBulkIn) * outer_area(Bundle)[i] * TotalSGTubeNumber
        # all nodes other than preheater
        
        q_Flux = (T_PrimaryBulkIn - T_SecondaryBulkIn) * U * 100 * 100 / 1000 #[kW/m^2]
#         if Bundle == selected_tubes[0]: print (HeatFlux, T_PrimaryBulkIn - 273.15, T_SecondaryBulkIn - 273.15, i)
        
        if Bundle.Length.label[i] != "preheater start" and Bundle.Length.label[i] != "preheater" and \
            Bundle.Length.label[i] != "thermal plate":
            
            # for one tube --> multiply by NumberTubes = for all tubes
            # U based on one tube (heat transfer area x number tubes) --> would cancel out in U calc (1/h*NA) NA
            # [W/ cm^2 K] * [K] * [cm^2] = [W] = [J/s]
                        
            # Ti = Ti+1 = Tsat
            T_SecondaryBulkOut = T_sat_secondary

            if x_pht > 0:
                DeltaH = EnthalpySaturatedSteam.magnitude - nc.enthalpyD2O_liquid(T_sat_primary)
#                 DeltaH = EnthalpySaturatedSteam.magnitude - nc.enthalpyH2O_liquid(nc.PrimarySidePressure, T_sat_primary)
                
                LatentHeatAvailable = x_pht * DeltaH * m_h_leakagecorrection
                # if x_pht = 0 (entire quality entering ith Bundle) is not enough to transfer all of heat need to SHT
                if LatentHeatAvailable < Q:
                    # this much heat needs to be transferred from PHT fluid as sensible heat
                    Q_left = Q - LatentHeatAvailable
                    
                    T_PrimaryBulkOut = T_PrimaryBulkIn - Q_left / (Cp_h * m_h_leakagecorrection * (1 - x_pht))
                    # no more latent heat/quality remains
                    x_pht = 0
                else:    
                    # latent heat being used for heat transfer, T doesn't change
                    T_PrimaryBulkOut = T_sat_primary
                    x_pht = x_pht - Q / (m_h_leakagecorrection * DeltaH)
    
                    if x_pht < 0:
                        x_pht = 0

            else:
                # no latent heat available, only sensible heat used
                T_PrimaryBulkOut = T_PrimaryBulkIn - Q / (Cp_h * m_h_leakagecorrection)
#             Q = U * (T_PrimaryBulkOut - T_SecondaryBulkOut) * outer_area(Bundle)[i] * Bundle.TubeNumber
#             x_out = (x_in + Q / (
#                 MassFlow_c_total.magnitude * EnthalpySaturatedSteam.magnitude  * (Bundle.TubeNumber / TotalSGTubeNumber)
#                 ))               

            if Bundle.Length.label[i] == "opposite preheater":
                MassFlow_c = MassFlow_c_total.magnitude
            else:
                MassFlow_c = MassFlow_c_total.magnitude
            
            x_in = sht_steam_quality(Q, T_sat_secondary, x_in, MassFlow_c, SecondarySidePressure)

            T_PrimaryBulkIn = T_PrimaryBulkOut
            T_SecondaryBulkIn = T_SecondaryBulkOut

            PrimaryBulk.append(T_PrimaryBulkOut)
            SecondaryBulk.append(T_SecondaryBulkOut)
            PrimaryWall.append(T_wh)
            SecondaryWall.append(T_wc)
            HeatFlux.append(q_Flux)
            
            if i > 11:
                x_in = (x_in * MassFlow_downcomer.magnitude) / MassFlow_c_total.magnitude
           
        # preheater Bundles 
        else:
            # just upstream of preheater
            if Bundle.Length.label[i] == "preheater start":
                Cp_h = nc.heatcapacityD2O_liquid(T_PrimaryBulkIn)
                Cp_c = nc.heatcapacityH2O_liquid(T_PreheaterIn, SecondarySidePressure)
                
                T_PrimaryBulkOut = T_PrimaryBulkIn - Q / (Cp_h * m_h_leakagecorrection)
                T_PrimaryBulkIn = T_PrimaryBulkOut
                
                C_min = Cp_c * MassFlow_preheater.magnitude  # [J/g K]*[g/s] = [J/Ks] ] [W/K]
                C_max = Cp_h * m_h_leakagecorrection
                C_r = C_min / C_max
                
                Q_max = C_min * (T_PrimaryBulkIn - T_PreheaterIn)
                
                # total area of only preheater area
                TotalArea = sum(outer_area(Bundle)[(i + 1):(Bundle.NodeNumber - 1)])
                NTU = U * TotalArea * TotalSGTubeNumber / C_min

#                 eta = (1 - np.exp(-NTU * (1 - C_r))) / (1 - C_r * np.exp(-NTU * (1 - C_r)))
                eta = 1 - np.exp(((NTU ** 0.28) / C_r) * (np.exp(-C_r * (NTU ** 0.78)) - 1))
                Q_NTU = eta * Q_max

                # 2 known temperatures at opposite ends (preheater in and PHT in at the top)
                # solving for other two ends (PHT out, and end of preheater before mixing with recirculating downcomer)
                T_PrimaryBulkOutEnd = T_PrimaryBulkIn - Q_NTU / (m_h_leakagecorrection * Cp_h)
                T_SecondaryBulkOutEnd = T_PreheaterIn + Q_NTU / (MassFlow_preheater.magnitude * Cp_c)
                # at end of preheater - about to mix with downcomer flow
                T_SecondaryBulkOut = T_sat_secondary#T_SecondaryBulkOutEnd
                T_SecondaryBulkIn = T_sat_secondary#T_SecondaryBulkOut
                
#                 x_in = (x_in * MassFlow_downcomer.magnitude) / MassFlow_c_total.magnitude
            # inside preheater
            else:
                x_in = 0
  
                # Guessing cold-side temperatures for remaining nodes inside preheater
                T_SecondaryBulkOut = T_SecondaryBulkIn - (T_sat_secondary - T_PreheaterIn) / 4
                
           
                T_PrimaryBulkOut = T_PrimaryBulkIn - Q / (Cp_h * m_h_leakagecorrection)
#                 print (i, T_SecondaryBulkIn-273.15, T_SecondaryBulkOut-273.15)
     
                # end of preheater 
                if i == Bundle.NodeNumber - 2:
                    T_SecondaryBulkOut = T_PreheaterIn
                    T_PrimaryBulkOut = T_PrimaryBulkOutEnd
                    # recirculating downcomer flow entering area under thermal plate
                    T_SecondaryBulkIn = T_sat_secondary

                elif i == Bundle.NodeNumber - 1:
                    T_SecondaryBulkIn = T_sat_secondary
                    T_SecondaryBulkOut = T_sat_secondary
                          
                else:
                    T_SecondaryBulkIn = T_SecondaryBulkOut
                   
                T_PrimaryBulkIn = T_PrimaryBulkOut
                        
            PrimaryBulk.append(T_PrimaryBulkOut)
            SecondaryBulk.append(T_SecondaryBulkOut)
            PrimaryWall.append(T_wh)
            SecondaryWall.append(T_wc)
            HeatFlux.append(q_Flux)

#     print ([j - 273.15 for j in PrimaryBulk])
#     print ([j-273.15 for j in PrimaryWall])
#     print ()
#     print ([j-273.15 for j in SecondaryWall])
#     print ([j-273.15 for j in SecondaryBulk])
#     print()
    return PrimaryBulk, HeatFlux


# separate function for divider plate leakage
def divider_plate(CalendarYear):
    None


def station_events(calendar_year):

    # Pressure changes only:
    
    if calendar_year < 1993:
        SecondarySidePressure = 4.593  # MPa
    
    # PLNGS pressure reduction in 1992 (september) by 125 kPa
    elif 1993 <= calendar_year < 1996.25:
        SecondarySidePressure = 4.593 - (125 / 1000)  # MPa
    
    # return to full boiler secondary side pressure, 4.593 MPa
    # pressure restored shortly after reactor back online from refurb.
    elif 1996.25 <= calendar_year < 1999:
        SecondarySidePressure = 4.593
    
    elif calendar_year >= 1999:
        SecondarySidePressure = 4.593 - (125 / 1000)  # MPa
    else:
        None
    
    
    PostOutageYearlyLeakage = 0.00035
    PostOutageInitialLeakage = 0.02   
    
    # Divider plate raplacement only:
    
    if calendar_year < YearOutageRestart: # (1983.25 - 1996, not including)
        # divider plate leakage rates estimated based on AECL work (2-3.5 % range estimated)
        # InitialLeakage = 0.035 # fraction of total SG inlet mass flow
        # YearlyRateLeakage = 0.0065 # yearly increase to fraction of total SG inlet mass flow
        InitialLeakage = 0.025
        YearlyRateLeakage = 0.0065  # yearly increase to fraction of total SG inlet mass flow
        InitialYear = YearStartup
         
    elif calendar_year == YearOutageRestart: #1996
        InitialLeakage = 0.0075  # controls where first post-outage point is (underpredicted w/o this)
        YearlyRateLeakage = 0
        InitialYear = YearOutageRestart
    
    elif calendar_year > YearOutageRestart: # after 1996
        # helps second post-outage point rise (second leak of 2.5% magnitude initialized)
        InitialLeakage = PostOutageInitialLeakage
        # either this or some SHT deposits help out Phase 4 from dipping so much
        YearlyRateLeakage = PostOutageYearlyLeakage
        InitialYear = YearOutageRestart    
    
    elif YearRefurbishment < calendar_year < YearRefurbRestart: # (2008.25 - 2014)
        InitialLeakage = (YearRefurbishment - YearOutageRestart) * PostOutageYearlyRateLeakage + PostOutageYearlyLeakage
        # plateau during outage 
        YearlyRateLeakage = 0
        InitialYear = YearRefurbishment #doesn't matter, growth rate is zero throughout refurb outage
    
    elif calendar_year >= YearRefurbRestart: # after and incl. 2014
        InitialLeakage = (YearRefurbishment - YearOutageRestart) * PostOutageYearlyLeakage + PostOutageInitialLeakage
        YearlyRateLeakage = 0.0015
        InitialYear = YearRefurbRestart
        
    else:
        None
        
    
    Leakage = InitialLeakage + (calendar_year - InitialYear) * YearlyRateLeakage
    
    DividerPlateMassFlow = MassFlow_h.magnitude * Leakage
    # decreases as divider (bypass) flow increases
    m_h_leakagecorrection = MassFlow_h.magnitude - DividerPlateMassFlow

    return SecondarySidePressure, m_h_leakagecorrection, DividerPlateMassFlow


def energy_balance(SteamGenerator, x_pht, j, SGFastMode):
    
    year = (j * nc.TIME_STEP / 8760) 
    CalendarYear = year + YearStartup

    if j == 0:
        x_pht = 0.01
    
    Energy = []
    
    Cleaned, TotalSGTubeNumber = primaryside_cleaned_tubes(SteamGenerator[0], CalendarYear)
    # adusts how many tubes per bundle to account for tubes plugged
    bundle_sizes(SteamGenerator, TotalSGTubeNumber)
    
    [SecondarySidePressure, RemainingPHTMassFlow, MasssFlow_dividerplate.magnitude] = station_events(CalendarYear)
    
    for Bundle in SteamGenerator:
        
        if SGFastMode == "yes":
            
            #not all tubes run (only those selected)
            Selected_Tubes = tube_picker(Method, SteamGenerator)[0]
            
            if Bundle in Selected_Tubes:
                # tracks oxide growth for these tubes specifically
                InnerOx = Bundle.InnerOxThickness
                OuterOx = Bundle.OuterOxThickness
            
            elif Bundle in Cleaned:
                # first tube in selected tube list (for that steam generator) simulates all cleaned tubes in "fast mode"
                CleanedTubeNumber = tube_picker(Method, SteamGenerator)[1][0]
                CleanedInnerOxide = SteamGenerator[CleanedTubeNumber].InnerOxThickness
                CleanedOuterOxide = SteamGenerator[CleanedTubeNumber].OuterOxThickness
                
                InnerOx = CleanedInnerOxide
                OuterOx = CleanedOuterOxide
  
            else:  # assumes same growth as in default passed tube for remaining tubes
                
                DefaultUncleanedInnerOxide = SteamGenerator[Default_Tube].InnerOxThickness
                DefaultUncleanedOuterOxide = SteamGenerator[Default_Tube].OuterOxThickness
                
                InnerOx = DefaultUncleanedInnerOxide
                OuterOx = DefaultUncleanedOuterOxide
                
        # in non-fast-mode all tubes are passed through (all cleaned/uncleaned) through oxide growth functions
        else:
            
            # prevents default (1.52 m, 57th index number) tube from being run twice through PHTS/oxide growth functions
            Selected_Tubes= SteamGenerator[0:57] + SteamGenerator[58:87]
            # if want to select all tubes as selected tubes
            
            if Bundle in Selected_Tubes:
                # tracks oxide growth for all tubes
                InnerOx = Bundle.InnerOxThickness
                OuterOx = Bundle.OuterOxThickness
            else:
                None
                       
        
        [Bundle.PrimaryBulkTemperature, Bundle.HeatFlux] = temperature_profile(
            Bundle, InnerOx, OuterOx, RemainingPHTMassFlow, SecondarySidePressure, x_pht, CalendarYear
            )
        
        SteamGeneratorOutputNode = SteamGenerator[Default_Tube].NodeNumber - 1
        
        m_timesH = (Bundle.TubeNumber / TotalSGTubeNumber) * RemainingPHTMassFlow \
        * nc.enthalpyD2O_liquid(Bundle.PrimaryBulkTemperature[SteamGeneratorOutputNode])
#         m_timesH = (Bundle.TubeNumber / TotalSGTubeNumber) * RemainingPHTMassFlow \
#             * nc.enthalpyH2O_liquid(nc.PrimarySidePressure, Bundle.PrimaryBulkTemperature[SteamGeneratorOutputNode])

        Energy.append(m_timesH)
        
        Enthalpy_dp_sat_liq = nc.enthalpyD2O_liquid(T_sat_primary)

        #nc.enthalpyH2O_liquid(nc.PrimarySidePressure, T_sat_primary)
        Enthalpy_dividerplate = Enthalpy_dp_sat_liq + x_pht * (EnthalpySaturatedSteam.magnitude - Enthalpy_dp_sat_liq)

        Enthalpy = (sum(Energy) + MasssFlow_dividerplate.magnitude * Enthalpy_dividerplate) / MassFlow_h.magnitude 

#     RIHT = nc.temperature_from_enthalpyH2O_liquid("PHT", Enthalpy, None)
   
    RIHT = nc.temperature_from_enthalpyD2O_liquid(Enthalpy)

    return RIHT

             
# print (energy_balance(ld.SteamGenerator_2, 0.02, 876 * 0, SGFastMode="yes")- 273.15)
