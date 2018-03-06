import csv
import numpy as np

# General electrochemistry constants
Beta = 0.5  # Symmetry coefficient
kH2 = 0.00078  # Henry's law constant for H2 @ 298.15 K [mol/L*atm]
F = 96485.3  # Faraday's constant [C/mol]
n = 2  # number of electrons transferred
R = 8.314  # Ideal gas constant [J/K *mol]
kb = 1.38E-23  # Boltzman constant
hp = 6.626E-34  # Planck's constant [m^2*kg/s]

# Molecular weights of compounds and elements [g/mol]
FeMolarMass = 55.847
OMolarMass = 15.99
Fe3O4MolarMass = 231.541
NiMolarMass = 58.694
CrMolarMass = 51.996
CoMolarMass = 58.9
NiFe2O4MolarMass = 234.388
FeCr2O4MoladMass = 223.779
H2MolarMass = 2.016
H2OMolarMass = H2MolarMass + OMolarMass

FractionCo_Alloy800 = 0.00015
FractionCr_Alloy800 = 0.21
FractionCo_CS = 0.00006
FractionCr_CS = 0.0002
FractionCr_Stellite = 0.27
FractionCo_Stellite = 0.65

# Density [g/cm^3]
Fe3O4Density = 5.2
NiFe2O4Density = 5.368
FeCr2O4Density = 4.8
H2Density = 8.988e-5
FeDensity = 7.86
NiDensity = 8.908
Alloy800Density = 7.94
CoDensity = 8.9

Fe3O4Porosity_inner = 0.1
Fe3O4Porosity_outer = 0.3
# FeCr2O4Porosity = 0.15
Fe3O4Tortuosity = 1.8
FeCr2O4Tortuosity = 1.2

# Diffusivity coefficients [cm^2/s]
FeDiffusivity = 0.00041
NiDiffusivity = 0.00041
CoDiffusivity = 0.00041

# Kinetic precipitation/dissolution constants [cm/s]
KpFe3O4 = .07
KdFe3O4 = .35
# KpFe_Ferrite = 0.014 #same as magnetite precipitation
# KdFe_Ferrite  = 0.044 'same as magnetite dissolution

FracNi_NiFe2O4 = 0.25
FracFe_Fe3O4 = 0.723
FracCr_Fe3O4 = 0.00019  # 0.0033
FracCo_Fe3O4 = 0.000049

H2 = 10  # [cm^3/kg]
SAFactor = 1.73  # surface area factor
CobaltWear = 0.0000000015  # mol/kg (spike input term)

ErosionConstant = 2.4e-11  # [g/cm^2 s]

OUTCORE_DEPOSITION = 0.0045
INCORE_DEPOSITION = 0.01


class SGParameters():
    def __init__(self):
        self.magnitude = None
        self.label = None
        self.unit = None
        self.steam_quality = None

PrimarySidePressure = 9.89  # MPa
# SecondarySidePressure = 4.593  # MPa


n_IAPWS = [
    0.14632971213167, -0.84548187169114, -0.37563603672040 * 10 ** 1, 0.33855169168385 * 10 ** 1, -0.95791963387872,
    0.15772038513228, -0.16616417199501 * 10 ** (-1), 0.81214629983568 * 10 ** (-3),
    0.28319080123804 * 10 ** (-3), -0.60706301565874 * 10 ** (-3), -0.18990068218419 * 10 ** (-1),
    - 0.32529748770505 * 10 ** (-1), -0.21841717175414 * 10 ** (-1), -0.52838357969930 * 10 ** (-4),
    - 0.47184321073267 * 10 ** (-3), -0.30001780793026 * 10 ** (-3), 0.47661393906987 * 10 ** (-4),
    - 0.44141845330846 * 10 ** (-5), -0.72694996297594 * 10 ** (-15), -0.31679644845054 * 10 ** (-4),
    - 0.28270797985312 * 10 ** (-5), -0.85205128120103 * 10 ** (-9), -0.22425281908000 * 10 ** (-5),
    - 0.65171222895601 * 10 ** (-6), -0.14341729937924 * 10 ** (-12), -0.40516996860117 * 10 ** (-6),
    - 0.12734301741641 * 10 ** (-8), -0.17424871230634 * 10 ** (-9), -0.68762131295531 * 10 ** (-18),
    0.14478307828521 * 10 ** (-19), 0.26335781662795 * 10 ** (-22), -0.11947622640071 * 10 ** (-22),
    0.18228094581404 * 10 ** (-23), -0.93537087292458 * 10 ** (-25)
    ]

I_IAPWS = [0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 3, 3, 3, 4, 4, 4, 5, 8, 8, 21, 23, 29, 30, 31, 32]
J_IAPWS = [
    - 2, -1, 0, 1, 2, 3, 4, 5, -9, -7, -1, 0, 1, 3 , -3, 0, 1, 3, 17, -4, 0, 6, -5, -2, 10, -8, -11, -6, -29, -31, -38,
    - 39, -40, -41
    ]

n_backwards_IAPWS = [
    - 0.23872489924521 * 10 ** 3, 0.40421188637945 * 10 ** 3, 0.11349746881718 * 10 ** 3, -0.58457616048039 * 10 ** 1,
    - 0.15285482413140 * 10 ** (-3), -0.10866707695377 * 10 ** (-5), -0.13391744872602 * 10 ** 2,
    0.43211039183559 * 10 ** 2, -0.54010067170506 * 10 ** 2, 0.30535892203916 * 10 ** 2, -0.65964749423638 * 10 ** 1,
    0.93965400878363 * 10 ** (-2), 0.11573647505340 * 10 ** (-6), -0.25858641282073 * 10 ** (-4),
    - 0.40644363084799 * 10 ** (-8), 0.66456186191635 * 10 ** (-7), 0.80670734103027 * 10 ** (-10),
    - 0.93477771213947 * 10 ** (-12), 0.58265442020601 * 10 ** (-14), -0.15020185953503 * 10 ** (-16)
    ]

I_backwards_IAPWS = [0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 2, 2, 3, 3, 4, 5, 6]

J_backwards_IAPWS = [0, 1, 2, 6, 22, 32, 0, 1, 2, 3, 4, 10, 32, 10, 32, 10, 32, 32, 32, 32] 

R_IAPWS = 0.461526  # [kJ/kg K] 


class ThermoDimensionlessInput():
    def __init__(self):
        self.p_ref = 16.53  # [MPa] reducing pressure [MPa]
        self.T_ref = 1386  # [K] reducing temperature [K]
        self.rho_ref = 0.3220  # [g/cm^3]
        self.mu_ref = 1  # [micro Pa s]
        self.lambda_ref = 1  # [mW/mK]
        self.enthalpy = 2500  # [J/g]


ReferenceValues = ThermoDimensionlessInput()


def enthalpy(side, Temperature, SecondarySidePressure):
    if side == "PHT" or side == "PHTS" or side == "phts" or side == "pht":
        p = PrimarySidePressure  # MPa
    else:
        p = SecondarySidePressure  # MPa secondary side

    ratio_pressures = p / ReferenceValues.p_ref  # p/p*, reduced pressure [unitless]
    ratio_temperatures = ReferenceValues.T_ref / Temperature  # T/T*, reduced temperature [unitless]

    Gibbs_T = sum(
        [x * ((7.1 - ratio_pressures) ** y) * z * ((ratio_temperatures - 1.222) ** (z - 1))
        for x, y, z in zip(n_IAPWS, I_IAPWS, J_IAPWS)]
        )
    return ratio_temperatures * Gibbs_T * R_IAPWS * Temperature  # [J/g mol]*[K] = [J/g]


def enthalpy_D2O(Temperature):
    T = Temperature - 273.15
    return 4.7307 * T - 546.55 # [kJ/kg]


def temperature_from_enthalpy_D2O(Enthalpy):
    T = (Enthalpy + 546.55) / 4.7307
    Temperature = T + 273.15
    return Temperature


def temperature_from_enthalpy(side, Enthalpy, SecondarySidePressure):
    if side == "PHT" or side == "PHTS" or side == "phts" or side == "pht":
        p = PrimarySidePressure
    else:
        p = SecondarySidePressure  # MPa secondary side
    p_ref = 1  # MPa
    T_ref = 1  # K
    
    ratio_pressures = p / p_ref  # p/p*, reduced pressure [unitless]
    ratio_enthalpies = Enthalpy / ReferenceValues.enthalpy
    
    ratio_temmperatures = sum(
        [x * ((ratio_pressures) ** y) * (ratio_enthalpies + 1) ** z
        for x, y, z in zip(n_backwards_IAPWS, I_backwards_IAPWS, J_backwards_IAPWS)]
        )
    
    return ratio_temmperatures * T_ref


def saturation_temperature():
    None

           
def density(side, Temperature, SecondarySidePressure):
    
    if side == "phts" or side == "PHTS" or side == "PHT":
        p = PrimarySidePressure
    
    elif side == "shts" or side == "SHTS" or side == "SHT":
        p = SecondarySidePressure 
    
    else:
        None
        
    ratio_pressures = p / ReferenceValues.p_ref
    ratio_temperatures = ReferenceValues.T_ref / Temperature  # ratio is reverse here from all other functions 
    
    Gibbs_p = sum(
        [-x * y * ((7.1 - ratio_pressures) ** (y - 1)) * (ratio_temperatures - 1.222) ** z
        for x, y, z in zip(n_IAPWS, I_IAPWS, J_IAPWS)]
        )
  
    return (1 / (R_IAPWS * Temperature * ratio_pressures * (Gibbs_p / (p * 1000)))) * 1000 / (100 ** 3)  # [g/cm^3] 

        
def D2O_density(Temperature):
    T = Temperature - 273
    return -0.0023 * T + 1.4646


def D2O_viscosity(Temperature):
    rho_ref = 0.3580  # [g/cm^3]
    T_ref = Temperature / 643.847
    mu_ref = 10 * 55.2651 / 1000000
    
    # D2O is in primary side only, don't need to carry through SHT pressure
    rho_rel = D2O_density(Temperature) / rho_ref

    no = [1.0, 0.940695, 0.578377, -0.202044]
    fi0 = T_ref ** 0.5 / sum([n / T_ref ** i for i, n in enumerate(no)])

    Li = [0, 1, 2, 3, 4, 5, 0, 1, 2, 3, 0, 1, 2, 5, 0, 1, 2, 3, 0, 1, 3, 5, 0, 1, 5, 3]
    Lj = [0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4, 5, 5, 5, 6]
    
    Lij = [0.4864192, -0.2448372, -0.8702035, 0.8716056, -1.051126, 0.3458395, 0.3509007, 1.315436, 1.297752, 1.353448,
           -0.2847572, -1.037026, -1.287846, -0.02148229, 0.07013759, 0.4660127, 0.2292075, -0.4857462, 0.01641220,
           -0.02884911, 0.1607171, -0.009603846, -0.01163815, -0.008239587, 0.004559914, -0.003886659]

    arr = [lij * (1 / T_ref - 1) ** i * (rho_rel - 1) ** j for i, j, lij in zip(Li, Lj, Lij)]
    fi1 = np.exp(rho_rel * sum(arr))
    
    return mu_ref * fi0 * fi1 #[g/cm s]
    
    
def viscosity(side, Temperature, SecondarySidePressure):
    
    T_ref = 647.096  # K
    rho_ref = ReferenceValues.rho_ref
    mu_ref = ReferenceValues.mu_ref
    
    i_1 = [0, 1, 2, 3, 0, 1, 2, 3, 5, 0, 1, 2, 3, 4, 0, 1, 0, 3, 4, 3, 5]
    j_1 = [0, 0, 0, 0, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 3, 3, 4, 4, 5, 6, 6]
    H_1 = [
    0.520094, 0.0850895, -1.08374, -0.289555, 0.222531, 0.999115, 1.88797, 1.26613, 0.120573, -0.281378, -0.906851,
    - 0.772479, -0.489837, -0.25704, 0.161913, 0.257399, -0.0325372, 0.0698452, 0.00872102, -0.00435673, -5.93E-04
    ]
        
    T_rel = Temperature / T_ref
    rho_rel = density(side, Temperature, SecondarySidePressure) / rho_ref
    
    H_0 = [1.67752, 2.20462, 0.6366564, -0.241605]
    terms = [0, 1, 2, 3]
    summation = sum([x / (T_rel) ** y for x, y in zip(H_0, terms)])
    
    mu_0 = 100 * np.sqrt(T_rel) / summation
    
    summation1 = sum([x * (((1 / T_rel) - 1) ** y) * ((rho_rel - 1) ** z) for x, y, z in zip(H_1, i_1, j_1)])
    mu_1 = np.exp(rho_rel * summation1)
    mu_2 = 1
    # [microPa s] --> [g/cm s]
    return (mu_ref * (mu_0 * mu_1 * mu_2) / 1000000) * 10  # [g/cm s]


def HeatCapacity(side, Temperature, SecondarySidePressure):
    if side == "PHT" or side == "PHTS" or side == "phts" or side == "pht":
        p = PrimarySidePressure
    else:
        p = SecondarySidePressure
    
    ratio_pressures = p / ReferenceValues.p_ref  # p/p*, reduced pressure [unitless]
    ratio_temperatures = ReferenceValues.T_ref / Temperature  # T/T*, reduced temperature [unitless]
    
    Gibbs_TT = sum(
        [x * ((7.1 - ratio_pressures) ** y) * z * (z - 1) * ((ratio_temperatures - 1.222) ** (z - 2))
        for x, y, z in zip(n_IAPWS, I_IAPWS, J_IAPWS)]
        )
    
    return (-ratio_temperatures ** 2) * Gibbs_TT * R_IAPWS  # [kJ/kg K] ([J/g K] or *H2OMolarMass for [J/mol K]
 
 
def thermal_conductivityH2O(side, Temperature, SecondarySidePressure):
#     if side == "PHT" or side == "PHTS" or side == "phts" or side=="pht":
#         p = PrimarySidePressure
#     else:
#         p = SecondarySidePressure
    T_ref = 647.096  # [K]
#     p_ref = 22.064 #[MPa]

    ratio_densities = density(side, Temperature, SecondarySidePressure) / ReferenceValues.rho_ref
    ratio_temperatures = Temperature / T_ref
    # ratio_pressures = p/p_ref

    L_IAPWS = [2.443221E-03, 1.323095E-02, 6.770357E-03, -3.454586E-03, 4.096266E-04]
    k_IAPWS = [0, 1, 2, 3, 4]
    summation = sum([x / (ratio_temperatures ** y) for x, y in zip(L_IAPWS, k_IAPWS)])

    lambda_0 = np.sqrt(ratio_temperatures) / summation

    i_IAPWS = [0, 1, 2, 3, 4]
    j_IAPWS = [0, 1, 2, 3, 4, 5]
    L1_IAPWS = [
        [1.60397357, -0.646013523, 0.111443906, 0.102997357, -0.050412363, 0.006098593],
        [2.33771842, -2.78843778, 1.53616167, -0.463045512, 0.083282702, -0.007192012],
        [2.19650529, -4.54580785, 3.55777244, -1.40944978, 0.275418278, -0.020593882],
        [-1.21051378, 1.60812989, -0.621178141, 0.071637322, 0, 0],
        [-2.720337, 4.57586331, -3.18369245, 1.1168348, -0.19268305, 0.012913842]
        ]

    summation1 = []
    for row in L1_IAPWS:
        y = sum([x * (ratio_densities - 1) ** y for x, y in zip(row, j_IAPWS)])
        summation1.append(y)

    lambda_1 = np.exp(
        ratio_densities * sum([(((1 / ratio_temperatures) - 1) ** x) * y for x, y in zip(i_IAPWS, summation1)])
        )

    return ((lambda_0 * lambda_1 * ReferenceValues.lambda_ref) / 1000) / 100  # [W/m K] 
