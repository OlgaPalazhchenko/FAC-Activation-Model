import lepreau_data as ld
import numpy as np
import thermochemistry_and_constants as nc
import random
import matplotlib.pyplot as plt 
import csv
import pandas as pd
from datetime import date, timedelta, datetime

    
def station_RIHT():
    Filename = 'C:\\Users\\opalazhc\\Dropbox\\PLNGS Modelling\\PLNGS_RIHT_raw.csv'
    df = pd.read_csv(Filename, header=None)
    
    Header1 = 'RIHT 2'
    Header2 = 'RIHT 4'
    Header3 = 'RIHT 6'
    Header4 = 'RIHT 8'
    
    df.columns = ['date', 'power', Header1, Header2, Header3, Header4]

    # converts from string object to time series  
    df.date = pd.to_datetime(df.date)
    df.set_index('date', inplace=True)
    
    # average of all 4 
    df['mean'] = df.iloc[:, 1:5].mean(axis = 1)
    
    #filters data to remove anything below 0.99 FP (except for phase 4, where power deratings took place to ~0.9 FP
    df = df[df['power'] > 0.90]
        
    #resamples all data by desired time period ('D' for daily) and performs operation, here taking the mean
    daily_average = df.resample('D').mean().copy()
    
    daily_average_no_missing = daily_average.dropna(axis = 0, how='any')
    
    
    # separate sub dataframes based on PLNGS 'Phases' of operation
    daily_average_phase3 = daily_average_no_missing['1996-1-1':'1998-10-1']
    daily_average_phase4 = daily_average_no_missing['1998-10-2':'2008-03-1']
    daily_average_phase5 = daily_average_no_missing['2013-1-29': '2013-12-14']
    daily_average_phase6 = daily_average_no_missing['2013-12-15' : '2017-08-02']
    daily_average_phase7 = daily_average_no_missing['2017-08-02':]
    
    daily_average_phase3 = daily_average_phase3[daily_average_phase3['power'] > 0.99]
    daily_average_phase4 = daily_average_phase4[daily_average_phase4['power'] > 0.90]
    daily_average_phase5 = daily_average_phase5[daily_average_phase5['power'] > 0.99]
    daily_average_phase6 = daily_average_phase6[daily_average_phase6['power'] > 0.99]
    daily_average_phase7 = daily_average_phase7[daily_average_phase7['power'] > 0.99]
    
    writer = pd.ExcelWriter('PLNGS_RIHT_by_Phase.xlsx', engine='xlsxwriter', datetime_format='mm-dd-yyyy')
    
    daily_average_phase3.to_excel(writer, sheet_name = 'Phase 3')
    daily_average_phase4.to_excel(writer, sheet_name = 'Phase 4')
    daily_average_phase5.to_excel(writer, sheet_name = 'Phase 5')
    daily_average_phase6.to_excel(writer, sheet_name = 'Phase 6')
    daily_average_phase7.to_excel(writer, sheet_name = 'Phase 7')
    
    # sets spacing between columns A and B so date column (A) is more clear
    workbook  = writer.book
    worksheet1 = writer.sheets['Phase 3']
    worksheet2 = writer.sheets['Phase 4']
    worksheet3 = writer.sheets['Phase 5']
    worksheet4 = writer.sheets['Phase 6']
    worksheet5 = writer.sheets['Phase 7']
    
    worksheets = [worksheet1, worksheet2, worksheet3, worksheet4, worksheet5]
    
    for sheet in worksheets:
        sheet.set_column('A:B', 12)
    
    writer.save()
    
station_RIHT()


def reactor_power():
    TrackedOutageDays = []
    EstimatedOutageDays = []
    AllStationDataDates = []
    
    Filename_stationpower = 'C:\\Users\\opalazhc\\Dropbox\\PLNGS Modelling\\PLNGS_Power_raw.csv'
    Filename_power_estimated = "C:\\Users\\opalazhc\\Dropbox\\PLNGS Modelling\\PLNGS_Power_raw_estimated (pre 1996).csv"
    
    df = pd.read_csv(Filename_stationpower, header=None)
    df.columns = ['date','power']

    # converts from string object to time series  
    df.date = pd.to_datetime(df.date)
    df.set_index('date', inplace=True)
    
    # frequency of power data changed from every few minutes to daily average
    daily_average = df.resample('D').mean().copy()

    # for days that are blank, e.g, 1996-8, forward filled, i.e., power from 1996-7 applied forward to 1996-8
    daily_average = daily_average.fillna(method='ffill')
#     monthly_average = monthly_average.dropna(how = "any")

    for date in daily_average.index:
        x = (date.year, date.month, date.day)
        AllStationDataDates.append(x)
    
    # secondary side pressure wil be added here eventually (maybe also primary side pressure)
#     writer = pd.ExcelWriter('PLNGS_Power.xlsx', engine='xlsxwriter', datetime_format='dd-mm-yyyy')
#     daily_average.to_excel(writer, sheet_name = 'Daily Averaged FP')
    
#     workbook = writer.book
#     worksheet1 = writer.sheets['Daily Averaged FP']
#     writer.save()
    
    TrackedOutages = daily_average.index[daily_average['power'] < 0.1]
    
    for date in TrackedOutages:
        x = (date.year, date.month, date.day)
        TrackedOutageDays.append(x)
   
   
    StationPower = daily_average.as_matrix(['power']).flatten()
    

    df2 = pd.read_csv(Filename_power_estimated, header=None)
    df2.columns = ['date','power']
    
    df2.date = pd.to_datetime(df2.date)
    df2.set_index('date', inplace=True)
    
    daily_average_2 = df2.resample('D').mean().copy()
    daily_average_2 = daily_average_2.fillna(method='ffill')
    
    writer = pd.ExcelWriter('PLNGS_Power_Estimated (Pre 1996).xlsx', engine='xlsxwriter', datetime_format='dd-mm-yyyy')
    daily_average_2.to_excel(writer, sheet_name = 'Daily Averaged FP')
    
    workbook = writer.book
    worksheet1 = writer.sheets['Daily Averaged FP']
    writer.save()
    
    EstimatedOutages = daily_average_2.index[daily_average_2['power'] < 0.1]
    
    
    for date in EstimatedOutages:
        x = (date.year, date.month, date.day)
        EstimatedOutageDays.append(x)
    
    # Estimated operating power (from estimates in spreadsheet)
    EstimatedPower = daily_average_2.as_matrix(['power']).flatten()
    
    # All outage days, estimated and from station data
    TrackedOutageDays = EstimatedOutageDays + TrackedOutageDays

#     plt.plot(df)
#     plt.xlabel('Year')
#     plt.ylabel('Power')
#     plt.show()
##     prints first 5 rows to check everything is assigned properly
#     print (df.head())

    # all dates for which station data available (1996-June 2018)
    # all outage days (estimated + station data), where power is < 0.1
    # only the estimated outage days pre 1996 station data
    # estimated operating power during the outage days
    # station power for dates where station data available (1996- June 2018)
    return AllStationDataDates, TrackedOutageDays, EstimatedOutageDays, EstimatedPower, StationPower


AllStationDataDates, TrackedOutageDays, EstimatedOutageDays, EstimatedPower, StationPower = reactor_power()
# print (len(OutageYearsMonths))

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

FullTubeComplement = 3542

DayStartup = (1983, 4, 8)
DayCPP = (1987, 3, 1)
DayOutage = (1995, 4, 13)
DayOutageRestart = (1996, 1, 1)
DayRefurbishment = (2008, 3, 29)
DayRefurbishmentRestart = (2012, 12, 10)

EventWeeks = []
EventDays = (DayStartup, DayCPP, DayOutage, DayOutageRestart, DayRefurbishment, DayRefurbishmentRestart) 
for x in EventDays:
    y = (datetime(*x).isocalendar()[0], datetime(*x).isocalendar()[1]) 
    EventWeeks.append(y)

WeekStartup, WeekCPP, WeekOutage, WeekOutageRestart, WeekRefurbishment, WeekRefurbishmentRestart = EventWeeks


FULLPOWER = 2061.4 * 1000

T_sat_primary = 310.77 + 273.15
T_PreheaterIn = 190 + 273.15
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

# # select desired SG tubes to be run by arc length
# UBends = [1.49, 0.685, 2.31, 3.09]
# UBends = [i * 100 for i in UBends]
# # select desired tubes to be run by total tube length
# TubeLengths = []
# for i in ld.SteamGenerator_2:
#     TubeLengths.append(i.Distance[21])


def total_tubes_plugged(SteamGenerator, Date):
    
    #plug events arbitrarily assumed to happen at the start of each year (no month/quarter granularity currently used) 
    if SteamGenerator == ld.SteamGenerator:
        DatePlugged = [
            (1983, 4, 8), (1988, 1, 8), (1992, 1, 8), (1993, 1, 8), (1994, 1, 8), (1996, 1, 8), (1999, 1, 8),
            (2002, 1, 8), (2009, 1, 8)
            ]
        AmountPlugged = [1, 1, 10, 6, 6, 1, 8, 1, 2]
        
    
    elif SteamGenerator == ld.SteamGenerator_2:
        DatePlugged = [(1983, 4, 8), (1996, 1, 8), (1999, 1, 8)]
        AmountPlugged = [6, 5, 11]


    if Date < DatePlugged[0]:
        NumberPluggedTubes = 0
    
    elif Date >= DatePlugged[len(DatePlugged) - 1]:
        NumberPluggedTubes = sum(AmountPlugged)
    
    else:
        for i in range(len(DatePlugged) - 1):

            if DatePlugged[i] <= Date < DatePlugged[i + 1]:
                # last tube in range not included, hence the + 1
                # e.g., SG 1: if date 1989-4, Date is <= i =1, so amount plugged summed from 0 to 1 --> 1 + 1 =2 
                NumberPluggedTubes = sum(AmountPlugged[0:(i + 1)])
        
    TotalSGTubeNumber = 3542 - NumberPluggedTubes
    
    return TotalSGTubeNumber


def bundle_sizes(SteamGenerator, TotalSGTubeNumber):
    
    TotalTubesinBundles = []
    for i in range(87):
        TotalTubesinBundles.append(SteamGenerator[i].TubeNumber)
        
    TotalTubesinBundles = sum(TotalTubesinBundles)

    if TotalSGTubeNumber < TotalTubesinBundles:
        
        for i in range(3542 - TotalSGTubeNumber):
            x = random.randint(0, 86)
            SteamGenerator[x].TubeNumber = SteamGenerator[x].TubeNumber - 1
    
    return None

   
def primaryside_cleaned_tubes(SteamGenerator, Date):
    
    CalendarDate = datetime(*Date).isocalendar() 
    WeeklyDate = (CalendarDate[0], CalendarDate[1])
    
    if SteamGenerator == ld.SteamGenerator:
        TubesCleaned = 2093
        Fraction_TubesCleaned = 0.969
    elif SteamGenerator == ld.SteamGenerator_2:
        TubesCleaned = 2105
        Fraction_TubesCleaned = 0.987
    else:
        None
    
    TotalSGTubeNumber = total_tubes_plugged(SteamGenerator, Date)
    
    #amount of tubes cleaned per each cleaning will have to be custom
    if Date == DayOutage or WeeklyDate == WeekOutage:
        # siva blast in 1995 used on only 60% of tubes due to time/spacial constraints

        PercentTubesCleaned = TubesCleaned / TotalSGTubeNumber
    
    elif Date == DayRefurbishment or WeeklyDate == WeekRefurbishment:

        PercentTubesCleaned = Fraction_TubesCleaned
    else:
        # none cleaned outside of outages, cleaned tubes list below has no appended bundles
        PercentTubesCleaned = 0
     
    # chooses tube bundles until % of total sg tube number reached, adds all chosen to "cleaned" tube list
    Cleaned = []
    NumberTubes = []
    AllTubes = list(range(0, 86))
    random.shuffle(AllTubes)
    
    for i in range(len(SteamGenerator)):
        x = AllTubes[i]
        
        NumberTubes.append(SteamGenerator[x].TubeNumber)
        
        if sum(NumberTubes) <= (PercentTubesCleaned * TotalSGTubeNumber):
            Cleaned.append(SteamGenerator[x])
        else:
            break

    return Cleaned


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


def sludge_fouling_resistance(Bundle, HeatTransferTimeStep, Date, i):
    
    CalendarDate = datetime(*Date).isocalendar() 
    WeeklyDate = (CalendarDate[0], CalendarDate[1])
    
    Time_Step = HeatTransferTimeStep / 24 / 365 #hr --> yr
    # default values
    ReducedTubeGrowth = 0.0003  # [g/cm^2] /year = 3.25 um/year
    
    # CPP installation (late 1986) reduces secondary side crud by 50% 
    
    if Date <= DayCPP:
        Growth = 0.0013#[g/cm^2]/yr
    
    elif Date in TrackedOutageDays:
        Growth = 0
    
    else:
        Growth = ReducedTubeGrowth 
        
    #estimated decrease in pre-existing sludge deposits on tubes due to CPP installation + draining + chemistry change
    if  Date == (1988, 4, 1) or WeeklyDate == (1988, 13):
        Bundle.SludgeLoading[i] = 0.75 * Bundle.SludgeLoading[i]
        
    elif Date == DayOutage or WeeklyDate == (1995, 15):
        Bundle.SludgeLoading[i] = Bundle.SludgeLoading[i] * 0.8
    
    elif Date == DayRefurbishment:
        Bundle.SludgeLoading[i] = Bundle.SludgeLoading[i] * 0.25
         
    else:
        Bundle.SludgeLoading[i] = Bundle.SludgeLoading[i] + Growth * Time_Step #  [g/cm^2] + [g/cm^2]/yr * 1/12th of a year
    
    Thickness = Bundle.SludgeLoading[i] / nc.Fe3O4Density
    
    Fouling = Thickness / thermal_conductivity(None, "outer magnetite", None)
    
    return Fouling # [cm^2 K/W]


def pht_fouling_resistance(i, InnerAccumulation, OuterAccumulation):
    
    # [g/cm^2]/[g/cm^3] = [cm]
    # thickness/thermal conductivity [cm]/[W/cm K] = [cm^2 K/W]
    InnerThickness = InnerAccumulation / nc.Fe3O4Density
    OuterThickness = OuterAccumulation / nc.Fe3O4Density 

    # [cm]/ [W/cm K] =[cm^2 K/W]
    # inner deposit is consolidated, predicted to have different thermal resistance than unconsolidated outer ox.
    InnerFouling = InnerThickness / thermal_conductivity(None, "inner magnetite", None)
    OuterFouling = OuterThickness / thermal_conductivity(None, "outer magnetite", None)
    # total pht thermal resistance
#     print (OuterAccumulation, OuterThickness, OuterFouling, InnerFouling + OuterFouling)
    return InnerFouling + OuterFouling  # [cm^2 K/W]


def conduction_resistance(Bundle, Twall, SecondarySidePressure, i):
    # R_conduction = L/kA (L = thickness)
    # A = area with shape factor correcrion factor for cylindrical pipe
    # R_conduction = ln(D_o/D_i)/(2pikl) [K/W]
    # l = length, k =thermal conductivity coefficient [W/cm K]
    
    k_w = thermal_conductivity(Twall, "Alloy-800", SecondarySidePressure)  # [W/cm K]
    
    r_o = Bundle.OuterDiameter[i] / 2
    r_i = Bundle.Diameter[i] / 2

    # small conductivity = big resistance
    Rcyl_numerator = np.log(r_o / r_i)
    Rcyl_denominator = 2 * np.pi * k_w * Bundle.Length.magnitude[i]

    R_cond = Rcyl_numerator / Rcyl_denominator  # [K/W]
    return R_cond


def primary_convection_resistance(
        Bundle, T_PrimaryBulkIn, T_PrimaryWall, x_pht, m_h_leakagecorrection, i, TotalSGTubeNumber
        ):
    
    # R_1,Convective = 1/(h_1,convectove *A)
    # A = inner area (based on inner diameter)
    # [W/cm K]/ [cm] = [W/K cm^2]
    
    # T_PrimaryBulkIn used in place of Tavg (of in and out for that node)
    
    Area_CrossSection = TotalSGTubeNumber * (np.pi / 4) * (Bundle.Diameter[i] ** 2)
    MassFlux_h.magnitude = m_h_leakagecorrection / Area_CrossSection
   
    Nu_D = nusseltnumber(Bundle, "PHT", T_PrimaryBulkIn, MassFlux_h.magnitude, i)
    
    h_l =  Nu_D * nc.thermal_conductivityD2O_liquid(T_PrimaryBulkIn) / Bundle.Diameter[i]
    
     # quality only persists for first 2 nodes of PHT
#     if x_pht > 0:
#         h_i = in_tube_boiling_LiuWinterton(x_pht, MassFlux_h.magnitude, Diameter, T_PrimaryWall, h_l, i)
#     else:
    h_i = h_l
    
    return 1 / (h_i * inner_area(Bundle)[i])  # [K/W]


def secondary_convection_resistance(
        Bundle, T_SecondaryBulkIn, T_SecondaryWall, T_SecondaryFilm, SecondarySidePressure, x_sht, i, TotalSGTubeNumber
        ):
    
#     print (T_SecondaryBulkIn-273.15, T_SecondaryWall-273.15, T_SecondaryFilm-273.15, SecondarySidePressure, x_sht,i)
    # split into boiling and non-boiling (preheater) Bundles
    if Bundle.Length.label[i] == "preheater":  # from Silpsrikul thesis
        
        # Based on PLGS data and Silpsrikul calculations
        f_b = 0.2119  # fraction of cross-Bundleal area of shell occupied by a baffle window
        N_b = f_b * TotalSGTubeNumber  # number of tubes in baffle window
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



def pht_steam_quality(Temperature, Date):
   
    CoreMassFlow =  7600  # [kg /s]
#     Delta_T = T_sat_primary - (262 + 273.15)  # [K]
#     C_p_cold = nc.heatcapacityD2O_liquid(262 + 273.15)
#     C_p_hot = nc.heatcapacityD2O_liquid(T_sat_primary)
#     C_p_avg = (C_p_cold + C_p_hot) / 2  # [kJ/kg K]

     #CoreMassFlow * C_p_avg * Delta_T  # [kW]
    
    difference = []

    # outages are accounted for in the 1996-01 to 2018-06 data via power derating (down to 0.01 or lower)
    if Date in AllStationDataDates:
        
        for x in AllStationDataDates:
            
            delta_time = [x1 - x2 for (x1, x2) in zip(Date, x)] 
            
            delta_time = [abs(number) for number in delta_time]
            
            difference.append(delta_time)

        # min function --> smallest difference between first index in all lists in tuple compared, followed by second     
        # index = location of list with this minimum difference relative to that of current (year, month) input 
        ClosestDay = difference.index(min(difference))
        
        Fraction_of_FP = StationPower[ClosestDay]
    
    # the outage in 1995 precedes the station log data available, so this is manually set to 0 % FP here in addition
    # to the other approximated outage dates
    
    # many of these estimated outages are not the entire month (so a recovery power is taken into account)
    elif Date in EstimatedOutageDays:
        Fraction_of_FP = 0
    
        
    # outside of data in the logs and the 1995 outage (1983-1995 data) 99 % full power is assumed
    else:
        Fraction_of_FP = 0.99
       
    P = FULLPOWER * Fraction_of_FP

    # core inlet enthalpy at RIHT + that added from fuel
    # H_fromfuel = Power / CoreMassFlow # [kJ/s /kg/s] = [kJ/kg]
    H_pht = P / CoreMassFlow

    H_satliq_outlet = nc.enthalpyD2O_liquid(T_sat_primary)
    
    # no quality in return primary flow
    H_current = nc.enthalpyD2O_liquid(Temperature) + H_pht
    x = (H_current - H_satliq_outlet) / (EnthalpySaturatedSteam.magnitude - H_satliq_outlet)
#     print (Date, x, P, Fraction_of_FP)
    if x < 0:
        x = 0
    return x, P

# pht_steam_quality(261.8 +273.15, (1991, 25, 10))

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
    print (MassFlow_c_total.magnitude, V_max)
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
     
#     if T_SecondaryWall <= T_sat_secondary:
#         T_SecondaryWall = T_sat_secondary #+ 1
         
    deltaT_sat = (T_SecondaryWall - T_sat_secondary) ** 0.24 # K
    P_w = nc.saturation_pressureH2O(T_SecondaryWall) * 10 ** 6
    P_s = nc.saturation_pressureH2O(T_sat_secondary) * 10 ** 6
     
    deltaP_sat = (P_w - P_s) ** 0.75 # MPa
         
    h = 0.00122 * lambda_l * cp_l * rho_l * deltaT_sat * deltaP_sat / (sigma * mu * h_lg * rho_v) # W / m^2 K

    
    h_nb = h / (100 ** 2) # [W/cm^2 K]
 
    return h_nb


def outside_bundle_pool_boiling(
        Correlation, Bundle, x_sht, T_sat_secondary, T_SecondaryWall, T_SecondaryBulkIn, SecondarySidePressure,
        i
        ):
    
    F = .19 # bundle boiling factor (empirical)
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

    return h_twophase
    

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
        Bundle, HeatTransferTimeStep, i, T_PrimaryBulkIn, T_SecondaryBulkIn, x_in, x_pht, InnerAccumulation,
        OuterAccumulation, Date, SecondarySidePressure, m_h_leakagecorrection, TotalSGTubeNumber):
    
    # i = each node of SG tube
    T_PrimaryWall = T_PrimaryBulkIn - (1 / 3) * (T_PrimaryBulkIn - T_SecondaryBulkIn)
    T_SecondaryWall = T_PrimaryBulkIn - (2 / 3) * (T_PrimaryBulkIn - T_SecondaryBulkIn)
    
    A_i = inner_area(Bundle)[i]
    A_o = outer_area(Bundle)[i]
    
    for k in range(100):
        
        WT_h = T_PrimaryWall
        WT_c = T_SecondaryWall

#         PrimaryT_film = (T_PrimaryBulkIn + T_PrimaryWall) / 2
        T_SecondaryFilm = (T_SecondaryBulkIn + T_SecondaryWall) / 2
               
        h_i.magnitude = (1 / (A_i * primary_convection_resistance(
            Bundle, T_PrimaryBulkIn, T_PrimaryWall, x_pht, m_h_leakagecorrection, i, TotalSGTubeNumber
            )))
        
        h_o.magnitude = (1 / (A_o * secondary_convection_resistance(
                Bundle, T_SecondaryBulkIn, T_SecondaryWall, T_SecondaryFilm, SecondarySidePressure, x_in, i,
                TotalSGTubeNumber
            )))
        # R = 1/hA --> h=A*1/R
        
        U_h1 = A_i * conduction_resistance(Bundle, T_PrimaryWall, SecondarySidePressure, i)
        
        U_h2 = (A_i * secondary_convection_resistance(
            Bundle, T_SecondaryBulkIn, T_SecondaryWall, T_SecondaryFilm, SecondarySidePressure, x_in, i,
            TotalSGTubeNumber
            ))
        U_h.magnitude = 1 / (U_h1 + U_h2)  # [W/ cm^2 K]

        U_c1 = A_o * conduction_resistance(Bundle, T_PrimaryWall, SecondarySidePressure, i)
        U_c2 = A_o * primary_convection_resistance(
            Bundle, T_PrimaryBulkIn, T_PrimaryWall, x_pht, m_h_leakagecorrection, i, TotalSGTubeNumber) 
        
        U_c.magnitude = 1 / (U_c1 + U_c2)  # [W/cm^2 K]

        T_PrimaryWall = ((T_PrimaryBulkIn * h_i.magnitude + T_SecondaryBulkIn * U_h.magnitude)\
            / (h_i.magnitude + U_h.magnitude))
        
        T_SecondaryWall = (T_SecondaryBulkIn * h_o.magnitude + T_PrimaryBulkIn * U_c.magnitude)\
            / (h_o.magnitude + U_c.magnitude)

        RE1 = (T_PrimaryWall - WT_h)
        RE2 = (T_SecondaryWall - WT_c)

        # if converged
        if abs(RE1) <= .05 and abs(RE2) <= .05:

            # [cm^2 K/W]
            R_F_primary.magnitude = pht_fouling_resistance(i, InnerAccumulation[i], OuterAccumulation[i]) 
            
            R_F_secondary.magnitude = sludge_fouling_resistance(Bundle, HeatTransferTimeStep, Date, i)
                        
            PCR = primary_convection_resistance(
                Bundle, T_PrimaryBulkIn, T_PrimaryWall, x_pht, m_h_leakagecorrection, i, TotalSGTubeNumber
                )
            
            SCR = secondary_convection_resistance(
                Bundle, T_SecondaryBulkIn, T_SecondaryWall, T_SecondaryFilm, SecondarySidePressure, x_in, i,
                TotalSGTubeNumber
               )
            
            # [cm^2 K/W] all resistances (convective, conductive, and fouling)
            CR = conduction_resistance(Bundle, T_PrimaryWall, SecondarySidePressure, i)
            
            inverseU_total = (PCR + CR + SCR) * A_o + R_F_primary.magnitude + R_F_secondary.magnitude

            U_total.magnitude = 1 / inverseU_total  # [W/ cm^2 K]
            
            
            return T_PrimaryWall, T_SecondaryWall, U_total.magnitude, h_i.magnitude


# def total_area_per_node(SteamGenerator):
#     TotalArea = []
#     for Bundle in SteamGenerator:
#         Areas = outer_area(Bundle)
#         A = [x * Bundle.TubeNumber for x in Areas]
#         TotalArea.append(A)
#     
#     TotalAreaBundles = []
#     for i in range(len(TotalArea[0])):
#         B = sum(x[i] for x in TotalArea)
#         TotalAreaBundles.append(B)
#     return TotalAreaBundles


def temperature_profile(
        Bundle, HeatTransferTimeStep, InnerOxide, OuterOxide, m_h_leakagecorrection, SecondarySidePressure, x_pht,
        Date, TotalSGTubeNumber
        ):
    
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
        
        T_wh, T_wc, U, h_i = wall_temperature(
            Bundle, HeatTransferTimeStep, i, T_PrimaryBulkIn, T_SecondaryBulkIn, x_in, x_pht, InnerOxide, OuterOxide,
            Date, SecondarySidePressure, m_h_leakagecorrection, TotalSGTubeNumber)
        
        
        # [W/ cm^2 K] * [K] * [cm^2] = [W] = [J/s]        
        Q = U * (T_PrimaryBulkIn - T_SecondaryBulkIn) * outer_area(Bundle)[i] * TotalSGTubeNumber
        
        q_Flux = (T_PrimaryBulkIn - T_SecondaryBulkIn) * U * 100 * 100 / 1000 #[kW/m^2]
        
        # all nodes other than preheater
        if Bundle.Length.label[i] != "preheater start" and Bundle.Length.label[i] != "preheater" and \
            Bundle.Length.label[i] != "thermal plate":
            
        
            # for one tube --> multiply by NumberTubes = for all tubes
            # U based on one tube (heat transfer area x number tubes) --> would cancel out in U calc (1/h*NA) NA
            # [W/ cm^2 K] * [K] * [cm^2] = [W] = [J/s]
                        
            # Ti = Ti+1 = Tsat
            T_SecondaryBulkOut = T_sat_secondary

            if x_pht > 0:
                DeltaH = EnthalpySaturatedSteam.magnitude - nc.enthalpyD2O_liquid(T_sat_primary)
                
                LatentHeatAvailable = x_pht * DeltaH * m_h_leakagecorrection

                # if x_pht = 0 (entire quality entering ith Bundle) is not enough to transfer all of heat need to SHT
                if LatentHeatAvailable < Q:
                    # this much heat needs to be transferred from PHT fluid as sensible heat
                    Q_left = Q - LatentHeatAvailable
                    
                    T_PrimaryBulkOut = T_PrimaryBulkIn - Q_left / (Cp_h * m_h_leakagecorrection * (1 - x_pht))
                    # no more latent heat/quality remains
                    x_pht = 0
                else:    
                    # latent heat being used for heat transfer, T doesn't change, PHT steam quality decreases
                    T_PrimaryBulkOut = T_sat_primary
                   
                    x_pht = x_pht - Q / (m_h_leakagecorrection * DeltaH)
    
                    if x_pht < 0:
                        x_pht = 0

            else:
            
                # no latent heat available, only sensible heat used
                T_PrimaryBulkOut = T_PrimaryBulkIn - Q / (Cp_h * m_h_leakagecorrection)
                    
#                 for k in range(20):
#                     if k > 0:
#                         Q = Q1
#                         
#                     # no latent heat available, only sensible heat used
#                     T_PrimaryBulkOut = T_PrimaryBulkIn - Q / (Cp_h * m_h_leakagecorrection)
#                     
#                     ln_term = np.log((T_PrimaryBulkOut - T_sat_secondary) / (T_PrimaryBulkIn - T_sat_secondary))
#                     LMTD = ((T_PrimaryBulkOut - T_sat_secondary) - (T_PrimaryBulkIn - T_sat_secondary)) / ln_term
#                     
#                     Q1 = U * outer_area(Bundle)[i] * TotalSGTubeNumber * LMTD
#                     RE = abs(Q1 - Q)
# #                     print (i, k, LMTD, Q, Q1, Q1/Q)
#                     if RE < 0.05:
#                         break 
            
            x_in = sht_steam_quality(Q, T_sat_secondary, x_in, MassFlow_c_total.magnitude, SecondarySidePressure)

            
            

            PrimaryBulk.append(T_PrimaryBulkOut)
            SecondaryBulk.append(T_SecondaryBulkOut)
            PrimaryWall.append(T_wh)
            SecondaryWall.append(T_wc)
            HeatFlux.append(q_Flux)
            
            T_PrimaryBulkIn = T_PrimaryBulkOut
            T_SecondaryBulkIn = T_SecondaryBulkOut
            
            if i > 11:
                x_in = 3#(x_in * MassFlow_downcomer.magnitude) / MassFlow_c_total.magnitude

           
        # preheater Bundles 
        else:
            # just upstream of preheater
            if Bundle.Length.label[i] == "preheater start":
                
                Cp_h = nc.heatcapacityD2O_liquid(T_PrimaryBulkIn)
                Cp_c = nc.heatcapacityH2O_liquid(T_PreheaterIn, SecondarySidePressure)
                
                T_PrimaryBulkOut = T_PrimaryBulkIn - Q / (Cp_h * m_h_leakagecorrection)
                
                C_min = Cp_c * MassFlow_preheater.magnitude  # [J/g K]*[g/s] = [J/Ks] ] [W/K]
                C_max = Cp_h * m_h_leakagecorrection
                C_r = C_min / C_max
                
                #needs to change to bulk out
                Q_max = C_min * (T_PrimaryBulkOut - T_PreheaterIn)
                
                # total area of only preheater area
                TotalArea = sum(outer_area(Bundle)[(i + 1):(Bundle.NodeNumber - 1)])
                NTU = U * TotalArea * TotalSGTubeNumber / C_min

#                 eta = 1 - np.exp(((NTU ** 0.22) / C_r) * (np.exp(-C_r * (NTU ** 0.78)) - 1))
                eta = 1 - np.exp(-(1 / C_r) * (1 - np.exp(-C_r * NTU)))
                
                Q_NTU = eta * Q_max

                # 2 known temperatures at opposite ends (preheater in and PHT in at the top)
                # solving for other two ends (PHT out, and end of preheater before mixing with recirculating downcomer)

                T_PrimaryBulkOutEnd = T_PrimaryBulkOut - Q_NTU / (m_h_leakagecorrection * Cp_h)
                T_SecondaryBulkOutEnd = T_PreheaterIn + Q_NTU / (MassFlow_preheater.magnitude * Cp_c)
                # at end of preheater - about to mix with downcomer flow

                T_SecondaryBulkOut = T_SecondaryBulkIn - (T_sat_secondary - T_PreheaterIn) / 5
                T_SecondaryBulkIn = T_SecondaryBulkOut
                
            # inside preheater
            else:
                x_in = 0
  
                #from end of first preheater node
                T_PrimaryBulkIn = T_PrimaryBulkOut
                Cp_h = nc.heatcapacityD2O_liquid(T_PrimaryBulkIn)
                # Guessing cold-side temperatures for remaining nodes inside preheater
                T_SecondaryBulkOut = T_SecondaryBulkIn - (T_sat_secondary - T_PreheaterIn) / 5
                
                T_PrimaryBulkOut = T_PrimaryBulkIn - Q / (Cp_h * m_h_leakagecorrection)
     
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

#     print ("Primary Bulk)", [j - 273.15 for j in PrimaryBulk])
#     print ("Primary Wall", [j-273.15 for j in PrimaryWall])
#     print ()
#     print ("Secondary Wall", [j-273.15 for j in SecondaryWall])
#     print ("Secondary Bulk", [j-273.15 for j in SecondaryBulk])
#     print()
    
    return PrimaryBulk, HeatFlux


def divider_plate(Date, HeatTransferTimeStep, DividerPlateLeakage):
    # time input (j) is in hours, converted to yearly
    Time_Step = HeatTransferTimeStep / 24 / 365 # hr --> yr
    
    PostOutageYearlyLeakage = 0.00035 # per year rate
    
    # (year, week number in the year, day number in the week) 
    CalendarDate = datetime(*Date).isocalendar() 
    WeeklyDate = (CalendarDate[0], CalendarDate[1])
    
    # changes to rate of leakage growth
    if Date in TrackedOutageDays:
        LeakageRate = 0 # per year rate
      
    elif Date == DayStartup:
        LeakageRate = 0
        
    elif Date < DayOutage:
        LeakageRate = 0.0018 # per year rate
    
    elif Date >= DayOutage:
        LeakageRate = PostOutageYearlyLeakage # per year rate
    
    else:
        None 
       
    
    # changes in overall leakage amount
    # replacement of divider plate
    
    # Leakage through the replaced welded divider plates immediately following replacement was estimated as 2% PHT flow
    # impulse input of leak site
    if Date == (1983, 5, 1) or WeeklyDate == (1983, 17):
        DividerPlateLeakage = 0.043
    
    elif Date == DayOutageRestart or WeeklyDate == WeekOutageRestart:
        DividerPlateLeakage = 0.025
    
    # Development of additional leak site
    # impulse input of leak site
    elif Date == (1996, 3, 1) or WeeklyDate == (1996, 9):
        DividerPlateLeakage = 0.0325
 
    DividerPlateLeakage = DividerPlateLeakage + Time_Step * LeakageRate
    
    # impulse input of leak site plateau
    if Date >= (1998, 11, 1):
        DividerPlateLeakage = 0.048

    return DividerPlateLeakage


def secondary_side_pressure(Date):
    # first drop estimated
    FirstPressureReduction = (1992, 11, 1)
    # second rediction date from plngs data
    SecondPressureReduction = (1998, 10, 25)
    
    if Date < FirstPressureReduction:
        SecondarySidePressure = 4.593  # MPa
    
    # PLNGS pressure reduction in 1992 (september) by 125 kPa
    elif FirstPressureReduction <= Date <= DayOutageRestart:
        SecondarySidePressure = 4.593 - (125 / 1000)  # MPa
    
    # return to full boiler secondary side pressure, 4.593 MPa
    # pressure restored shortly after reactor back online from refurb.
    elif DayOutageRestart < Date < SecondPressureReduction:
        SecondarySidePressure = 4.593
    
    elif Date >= SecondPressureReduction:
        SecondarySidePressure = 4.593 - (125 / 1000)  # MPa
    else:
        None

    return SecondarySidePressure


Cleaned_OutageSG1 = primaryside_cleaned_tubes(ld.SteamGenerator, DayOutage)
CleanedOutageSG2 = primaryside_cleaned_tubes(ld.SteamGenerator_2, DayOutage)

CleanedRefurbishmentSG1 = primaryside_cleaned_tubes(ld.SteamGenerator, DayRefurbishment)
CleanedRefurbishmentSG2 = primaryside_cleaned_tubes(ld.SteamGenerator_2, DayRefurbishment)


def energy_balance(SteamGenerator, x_pht, DividerPlateLeakage, Date, HeatTransferTimeStep, SGFastMode):
   
    Energy = []

    TotalSGTubeNumber = total_tubes_plugged(SteamGenerator, Date)
    
    # adjusts how many tubes per bundle to account for tubes plugged
    bundle_sizes(SteamGenerator, TotalSGTubeNumber)
    
    SecondarySidePressure = secondary_side_pressure(Date)
    
    MasssFlow_dividerplate.magnitude = MassFlow_h.magnitude * DividerPlateLeakage
    # decreases as divider (bypass) flow increases
    
    # g / s
    RemainingPHTMassFlow = MassFlow_h.magnitude - MasssFlow_dividerplate.magnitude

    for Bundle in SteamGenerator:
        
        [Bundle.PrimaryBulkTemperature, Bundle.HeatFlux] = temperature_profile(
            Bundle, HeatTransferTimeStep, Bundle.InnerOxLoading, Bundle.OuterOxLoading, RemainingPHTMassFlow,
            SecondarySidePressure, x_pht, Date, TotalSGTubeNumber
            )
        
        SteamGeneratorOutputNode = SteamGenerator[Default_Tube].NodeNumber - 1
        
        m_timesH = (Bundle.TubeNumber / TotalSGTubeNumber) * RemainingPHTMassFlow \
        * nc.enthalpyD2O_liquid(Bundle.PrimaryBulkTemperature[SteamGeneratorOutputNode])

        Energy.append(m_timesH)
        
        Enthalpy_dp_sat_liq = nc.enthalpyD2O_liquid(T_sat_primary)

        Enthalpy_dividerplate = Enthalpy_dp_sat_liq + x_pht * (EnthalpySaturatedSteam.magnitude - Enthalpy_dp_sat_liq)

        Enthalpy = (sum(Energy) + MasssFlow_dividerplate.magnitude * Enthalpy_dividerplate) / MassFlow_h.magnitude
    
    
    RIHT = nc.temperature_from_enthalpyD2O_liquid(Enthalpy)
       
    return RIHT

        
# print (energy_balance(
#     ld.SteamGenerator_2, 0.01, DividerPlateLeakage= 0.03, Date= (1984, 4, 8), HeatTransferTimeStep = nc.TIME_STEP * 1,
#     SGFastMode="yes"
#     )- 273.15)
