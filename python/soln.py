#! /Library/Frameworks/Python.framework/Versions/Current/bin/python


# Functions for calculating various properties of aqueous solutions
# including H2O/D2O mixtures
# 
#  Modified 8 December 2011 so that for urea solution calculations
# the parameter vfd refers to the volume fraction of D2O for the
# total solution volume, not the water component as used before.


########### constants  #################
k = 1.380658E-23 # Boltzmann constant in J/K
Mh2o = 18.0151528  # molecular weight of H2O in g/mol
Md2o = 20.027604  # molecular weight of D2O in g/mol  
Mu = 60.056     # molecular weight of urea
Na = 6.0221367E23 # Avogadro's number
be = 0.28179E-12 # X-ray scattering length of an electron, in cm
bH = -3.7409E-13 # neutron scattering length of natural H (99.985% 1H), in cm
bD = 6.674E-13   # neutron scattering length of 2H, in cm
bC = 6.6484E-13  # neutron scattering length of natural C (98.89% 12C), in cm
bN = 9.36E-13  # neutron scattering length of natural N (99.635% 14N), in cm
bO = 5.805E-13  # neutron scattering length of natural O (99.75% 16O), in cm
bh2o = 2.0*bH + bO # neutron scattering length of natural h2o
bd2o = 2.0*bD + bO # neutron scattering length of d2o
bu = 4.0*bH + 2.0*bN + bC + bO  # neutron scattering length of urea

def h2oDens(t):
    """Calculates mass density of H2O as a function of 
    temperature, which is used for calibration of absolute scattering
    intensities from reference watering scattering intensity.  
    Input is temperature, in Celsius.
    Output is density in g/cm^3
    Uses equation 16 from:	Kell, G. S. (1975). Density, thermal expansivity, 
    and compressibility of liquid water from 0 to 100 C: Correlations and 
    tables for atmospheric pressure and stauration reviewed and expressed 
    on 1968 temperature scale. J. Chem. Eng. Data, 20, 97-105. 
    http://dx.doi.org/10.1021/je60064a005"""
    
    a = 999.83952
    b = 16.945176
    c = -7.9870401E-3
    d = -46.170461E-6
    e = 105.56302E-9
    f = -280.54253E-12
    g = 16.879850E-3
    
    dens = 1.0E-3*(a + b*t + c*t**2 + d*t**3 + e*t**4 + f*t**5)/(1.0 + g*t)
    return dens

def h2oCompr(t):
    """Calculates isothermal compressiblity of H2O as a function of 
    temperature, which is used for calibration of absolute scattering
    intensities from reference watering scattering intensity.  
    Input is temperature, in Celsius.
    Output is compressibility in Pa-1
    Uses equation 20 from:	Kell, G. S. (1975). Density, thermal expansivity, 
    and compressibility of liquid water from 0 to 100 C: Correlations and 
    tables for atmospheric pressure and stauration reviewed and expressed 
    on 1968 temperature scale. J. Chem. Eng. Data, 20, 97-105. 
    http://dx.doi.org/10.1021/je60064a005"""
    
    a = 50.88496
    b = 0.6163813
    c = 1.459187E-3
    d = 20.08438E-6
    e = -58.47727E-9
    f = 410.4110E-12
    g = 19.67348E-3
    
    comp = 1.0E-11*(a + b*t + c*t**2 + d*t**3 + e*t**4 + f*t**5)/(1.0 + g*t)
    return comp
    

def d2oDens(t):

    """Calculates mass density of D2O as a function of temperature.
    intensities from reference watering scattering intensity.  
    Input is temperature, in Celsius.
    Output is density in g/cm^3
    Uses data from
        Kell, G. S. (1967). Precise representation of volume properties 
        of water at one atmosphere. J. Chem. Eng. Data, 12, 66-69. 
        http://dx.doi.org/10.1021/je60032a018"""
    
    a = 1.104690
    b = 20.09315E-3
    c = -9.24227E-6
    d = -55.9509E-9
    e = 79.9512E-12
    f = 0.0
    g = 17.96190E-3
    
    dens = (a + b*t + c*t**2 + d*t**3 + e*t**4 + f*t**5)/(1.0 + g*t)
    return dens
    

def ureaSolnDens(cUrea,vfd):
    """Calculates mass density of urea solution at 25 C in H2O/D2O mixtures.
    Inputs are molar concentration of urea and volume fraction of D2O
    volume fraction is fraction of *total* volume.
    Checks to make sure that volume fraction of D2O does not exceed
    total volume fraction available for water.  If it does, the function 
    returns a negative value, equal to -1 * the maximum vfd for the
    specified urea concentration.
    Uses data of:
        Gucker, Jr., F. T., Gage, F. W. & Moser, C. E. (1938). 
        The densities of aqueous solutions of urea at 25 and 30 and 
        the apparent molal volume of urea. J. Am. Chem. Soc., 60, 2582-2588. 
        http://dx.doi.org/10.1021/ja01278a008
    to calculate density of urea solutions in H2O at 25 C
    and data of:
        Kell, G. S. (1967). Precise representation of volume properties 
        of water at one atmosphere. J. Chem. Eng. Data, 12, 66-69. 
        http://dx.doi.org/10.1021/je60032a018
        
        Kell, G. S. (1975). Density, thermal expansivity, and compressibility
        of liquid water from 0 to 150 C: Correlations and tables for atmospheric 
        pressure and stauration reviewed and expressed on 1968 temperature scale. 
        J. Chem. Eng. Data, 20, 97-105. 
        http://dx.doi.org/10.1021/je60064a005
    For densities of H2O and D2O at 25 C.  Applies a small fudge factor to adjust
    the Gucker et al. data so that 0 M urea density matches the Kell data at 25 C.
    Assumes that volumes of H2O and D2O add ideally, with and without urea.

    """
    densH2O = 0.997045  # at 25 C
    densD2O= 1.10441    # at 25 C

    ##### Equation from fit of data of Gucker et al. ############
    a0 = 0.99715
    a1 = 0.015826
    a2 = -0.00010036
    
    fudge = densH2O/a0  # a small fudge factor that scales the published urea solution densities
                        # to published H2O and D2O densities
    
    dens = fudge*(a0 + a1*cUrea + a2*cUrea**2) # density of the urea solution if it were in pure H2O
    
    if vfd > 0.0:
        # mass of water per cm^3 of the hypothetical urea solution in H2O
        mWat = dens - cUrea*Mu/1000.0  
        # volume of water per cm^3 of the hypothetical urea solution in H2O
        volWat = mWat/densH2O         
        volD2O = vfd                  # volume D2O per cm^3
        if volD2O > volWat:
            return -1.0*volWat
        volH2O = volWat - volD2O      # volume H2O per cm^3
        
        mWater = volH2O*densH2O + volD2O*densD2O  # mass of H2O + D2O 
        mTotal = mWater + cUrea*Mu/1000.0   # Total mass per cm^3
        return mTotal

    return dens

def rhoXrayWat(t, vfd):
    densH2O = h2oDens(t)
    densD2O = d2oDens(t)
    densWat = mixDens(vfd,densH2O,densD2O)

    massD2O = mFracB(vfd,densH2O,densD2O)*densWat # mass of D2O per cm^3
    massH2O = densWat - massD2O       # mass of H2O per cm^3
    
    rhoH2O = 10.0*be*Na*massH2O/Mh2o   # contribution to scattering density from H2O
    rhoD2O = 10.0*be*Na*massD2O/Md2o   # contribution to scattering density from D2O
        
    return rhoH2O + rhoD2O
    
def rhoNeutWat(t, vfd):
    densH2O = h2oDens(t)
    densD2O = d2oDens(t)
    densWat = mixDens(vfd,densH2O,densD2O)

    massD2O = mFracB(vfd,densH2O,densD2O)*densWat # mass of D2O per cm^3
    massH2O = densWat - massD2O       # mass of H2O per cm^3
    
    rhoH2O = bh2o*Na*massH2O/Mh2o   # contribution to scattering density from H2O
    rhoD2O = bd2o*Na*massD2O/Md2o   # contribution to scattering density from D2O
        
    return rhoH2O + rhoD2O    

def h2oI0(t):
    """Calculates absolute scattering intensity of liquid H2O at 0 scattering angle.
    Input is temperature in Celsius.  Output is scattering intensity in cm-1
    Equations from:
    Orthaber, D., Bergmann, A. & Glatter, O. (2000). SAXS experiments on absolute scale
    with Kratky systems using water as a secondary standard. J. Appl. Cryst., 33, 218-225. 
    http://dx.doi.org/10.1107/S0021889899015216
    and from http://physchem.kfunigraz.ac.at/sm/
    """
    
    ePerWater = 10 # number of electrons per water molecule
    dens = h2oDens(t) # density of water in g/cm^3
    comp = h2oCompr(t) # compressibility of water in Pa-1
    T = t+273.15 # absolute temperature
    
    s0 = 1E6*(Na*dens/Mh2o)*k*T*comp
    
    i0 = (ePerWater*be)**2*s0*Na*dens/Mh2o
    return i0


    
def rhoXrayUrea(cUrea, vfd):
    """Calculates neutron scattering density of urea solution, with 
    specified volume fraction of D2O.
    vfd is volume fraction D2O of total solution"""

    densH2O = 0.997045  # at 25 C
    densD2O= 1.10441    # at 25 C

    dens = ureaSolnDens(cUrea,vfd)
    massUrea = cUrea*Mu/1000.0  # mass of urea per cm^3
    massWater = dens - massUrea # mass of H2O + D2O per cm^3
    massD2O = vfd*densD2O # mass of D2O per cm^3
    massH2O = massWater - massD2O       # mass of H2O per cm^3
    
    rhou = 32.0*be*Na*cUrea/1000.0       # contribution to scattering density from urea
    rhoH2O = 10.0*be*Na*massH2O/Mh2o   # contribution to scattering density from H2O
    rhoD2O = 10.0*be*Na*massD2O/Md2o   # contribution to scattering density from D2O
        
    return rhou + rhoH2O + rhoD2O

def rhoNeutUrea(cUrea, vfd):
    """Calculates neutron scattering density of urea solution, with 
    specified volume fraction of D2O
    vfd is volume fraction D2O of total solution"""
    
    densH2O = 0.997045  # at 25 C
    densD2O= 1.10441    # at 25 C
    
    dens = ureaSolnDens(cUrea,vfd)
    massUrea = cUrea*Mu/1000.0  # mass of urea per cm^3
    massWater = dens - massUrea # mass of H2O + D2O per cm^3
    massD2O = vfd*densD2O # mass of D2O per cm^3
    massH2O = massWater - massD2O       # mass of H2O per cm^3
    
    rhou = bu*Na*cUrea/1000.0       # contribution to scattering density from urea
    rhoH2O = bh2o*Na*massH2O/Mh2o   # contribution to scattering density from H2O
    rhoD2O = bd2o*Na*massD2O/Md2o   # contribution to scattering density from D2O
        
    return rhou + rhoH2O + rhoD2O
       
def mFracB(vFracB, densA, densB):
    mfd = vFracB*densB/(densA+vFracB*(densB-densA))
    return mfd
    
def vFracB(mFracB, densA, densB):
    vfd = mFracB*densA/(densB-mFracB*(densB-densA))
    return vfd

def mixDens(vFracB,densA,densB):
    dens = densA + vFracB*(densB - densA)
    return dens 
