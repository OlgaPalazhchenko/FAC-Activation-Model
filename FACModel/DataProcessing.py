'''
Created on Oct 22, 2017

@author: opalazhc
'''
for j in range(7000):#nc.SimulationDuration
    I = RunModel(ld.Inlet, ld.Core, j)
    C = RunModel(ld.Core, ld.Outlet, j)
    O = RunModel(ld.Outlet, ld.SteamGenerator,  j)
    #Sets input concentrations for all SG zones to be same as input of first (Outlet output)
#     for Zone in ld.SGZones[1:len(ld.SGZones)]:
#         BulkOutletOutput= [O.Section1.Bulk.FeTotal, O.Section1.Bulk.NiTotal, O.Section1.Bulk.CoTotal, O.Section1.Bulk.CrTotal]
#         BulkSGInput = [Zone.Bulk.FeTotal, Zone.Bulk.NiTotal, Zone.Bulk.CoTotal, Zone.Bulk.CrTotal]
#         for x,y in zip(BulkSGInput,BulkOutletOutput):
#             x[0] = y[O.Section1.NodeNumber-1]
    
    SG = RunModel(ld.SteamGenerator, ld.Inlet, j)
#     SG1 = RunModel(ld.SG_Zone1, ld.Inlet, j)
#     SG2 = RunModel(ld.SG_Zone2, ld.Inlet, j)
#     SG3 = RunModel(ld.SG_Zone3, ld.Inlet, j)
     
#     I.Section1.PrimaryBulkTemperature = SGHX.EnergyBalance(21, I.Section1.NodeNumber)
#     I.Section1.Bulk.FeSatFe3O4 = c.IronSolubility(I.Section1)
#     print([i-273.15 for i in I.Section1.PrimaryBulkTemperature])
     
#     if j %699==0: #8759
#         RIHT.append(ld.Inlet.PrimaryBulkTemperature[0]-273.15)
#         T_SG.append(ld.SteamGenerator.PrimaryBulkTemperature[21]-273.15)
#         T_SG1.append(ld.SG_Zone1.PrimaryBulkTemperature[21]-273.15)
#         T_SG2.append(ld.SG_Zone2.PrimaryBulkTemperature[21]-273.15)
#         T_SG3.append(ld.SG_Zone3.PrimaryBulkTemperature[21]-273.15)
#           
#         for Zone in ld.SGZones:    
#             totalthicknessSG=[x+y for x,y in zip(Zone.OuterOxThickness, Zone.InnerOxThickness)]
#             convertedtotalthicknessSG = ld.UnitConverter(Zone, "Grams per Cm Squared", "Grams per M Squared", None, None, totalthicknessSG, None, None, None)
#             coldlegavg = sum(convertedtotalthicknessSG[11:22])/(22-11)
#             hotlegavg = sum(convertedtotalthicknessSG[0:10])/(22-11)
#             coldlegoxide.append(coldlegavg)
#             hotlegoxide.append(hotlegavg)
#           
#         CorrosionRate = ld.UnitConverter(ld.Outlet, "Corrosion Rate Grams", "Corrosion Rate Micrometers", None, O.Section1.CorrRate, None, None, None, None)   
#         FACrate.append(sum(CorrosionRate)/len(CorrosionRate))  

hours = delta_time//3600
temp = delta_time - 3600*hours
minutes = delta_time//60
seconds = delta_time - 60*minutes
print('%d:%d:%d' %(hours,minutes,seconds))

# print (FACrate, "FAC rate, [um/a]")
# print()
# print (RIHT, "RIHT, [oC]")
# print(T_SG, "SteamGenerator")
# print (T_SG1, "SG_Zone1")
# print (T_SG2, "SG_Zone2")
# print (T_SG3, "SG_Zone3")
# print()
# print(coldlegoxide, "cold leg averages, [g/m^2]")
# print(hotlegoxide, "hot leg averages, [g/m^2]")    