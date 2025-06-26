#! /Library/Frameworks/Python.framework/Versions/Current/bin/python

# Functions for manipulating saxs data
# including file read and write,
# data plotting, slit-smearing functions
# and a general data smoothing function from 
# the SciPy Cookbook: http://www.scipy.org/Cookbook/SignalSmooth

#  (c) 2009-2011 by David P. Goldenberg
#  Please send feature requests, bug reports, or feedback to this address:
#           Department of Biology
#           University of Utah
#           257 South 1400 East
#           Salt Lake City, UT
#     
#           goldenberg@biology.utah.edu
#     
#  This software is distributed under the conditions of the BSD license.
#  Please see the documentation for further details.


import types
import numpy as np
import sys

from scipy import interpolate
from scipy import integrate
import matplotlib.pyplot as plt 


######## Data File Read and Write Functions    #############

def readPdh(pdhFileName):
    """Opens and reads a saxs data file in the "PDH" (Primary Data Handling) format
    used in the  PCG SAXS software suite developed by the Glatter group at the
    University of Graz, Austria This format is described in the appendix 
    of the PCG manual (page 123 of the 2005 version).
    In the PDH format, lines 1-5 contain header information, 
    followed by the SAXS data.
    Line 3 contains 8 integer constants, and lines 4 and 5 each contain five 
    floating-point constants each.
    The first integer constant in line 3 is reserved for the number of data points.
    The 4th floating point constant in line 4 is reserved for a 
    normalization  factor.
    The software for the SAXSess instrument uses two additional floating-point
    constants in line 4:
      2nd constant: sample to detector distance (mm)
      5th constant: x-ray wavelength
    The Utah SAXS tools also use one integer constant in line 3 and four 
    floating point constants in line 5:
    2nd integer constant in line 3: beam length profile type 
    defined by three possible values:
       0: no smearing applied to model functions
       1: sigmoidal beam-length profile
       2: trapezoidal beam-length profile
    Floating-point parameters in line 5: 
       1: beam-length profile parameter a 
       2: beam-length profile parameter b
       3: half-width of the beam-width profile (at half-height)
       4: detector slit length (equivalent to width of integration area
          used for image plate in SAXSess)
       5: q units scale: 1.0 for A-1, 10.0 for nm-1
    This function returns two lists of lists:
    [line1,line2,line3fields,line4fields,line5fields]
    [q, I, iErr]
    In the first list, line1 is a text string, line 2 is a list of text strings, and the other three items
    are lists of the numeric parameters, formatted as int (line3) 
    or floats (line4 and5).
    In the second list, q, I and iErr are the saxs data, in lists of floats
    """
    
    if pdhFileName == None:
        dataFile = sys.stdin
    else:
        dataFile = open(pdhFileName,'r')

    # first line is experiment description
    line1 = dataFile.readline().rstrip()
    # second line contains up to 16 keywords of up to 4 characters each,
    # in 4-space fields separated by one blank space.  
    # Usually used to specify experiment type. (eg "SAXS")
    line = dataFile.readline().rstrip()
    line2fields  = line.split()
    # third line is integer parameters
    line = dataFile.readline().rstrip()
    line3fields  = line.split()
    for i in range(8):
        line3fields[i] = int(line3fields[i])
    # fourth line is floating-point parameters
    line = dataFile.readline()[:-1]
    line4fields  = line.split()
    for i in range(5):
        line4fields[i] = float(line4fields[i])
    # fifth line is floating-point parameters
    line = dataFile.readline()
    line5fields  = line.split()
    for i in range(5):
        line5fields[i] = float(line5fields[i])

    # saxs data lines
    # Reads experimental data, selecting lines with exactly three fields,
    # And checking that all three are floating point numbers
    q=[]
    I=[]
    iErr=[]
    for line in dataFile.readlines():
        line = line[:-1]    # Chop newline character
        fields = line.split()     # Split line into fields
        if len(fields) ==3:
            try:
                thisQ = float(fields[0])
                thisI = float(fields[1])
                thisErr = float(fields[2])
                q.append(thisQ)
                I.append(thisI)
                iErr.append(thisErr)
            except ValueError:
                continue
            
    q = np.array(q)
    I = np.array(I)
    iErr = np.array(iErr)

    dataFile.close()
    
    return [line1,line2fields,line3fields,line4fields,line5fields],[q, I, iErr]

def writePdh(pdhFileName, header,data):
    
    if pdhFileName != None:
        outFile = open(pdhFileName,'w')
    else: 
        outFile = None
    
    print >> outFile, header[0]
    for field in header[1]:
        print >> outFile, field[:4],
    print >> outFile, ""
    for field in header[2]:
        print >> outFile, "%9i" % field,
    print >> outFile, ""  
    for field in header[3]:
        print >> outFile, "%14.6E" % field,
    print >> outFile, ""  
    for field in header[4]:
        print >> outFile, "%14.6E" % field,
    print >> outFile, ""
    
    for i in range(len(data[0])):
        print >> outFile, "%14.6E" % data[0][i],
        print >> outFile, "%14.6E" % data[1][i],
        print >> outFile, "%14.6E" % data[2][i],
        print >> outFile, ""    
    
    if outFile != None:
        outFile.close()
        
def readQIE(dataFileName, numbCols=3):
    """Reads generic SAXS data file for Q, I, and Ierr
    Reads experimental data, selecting lines that do not begin with #
    and contain exactly three fields with floating point numbers, or
    the number of fields specified by the optional numbCols argrument.
    number must be >=3
    Returns 3 lists, containing Q, I and Ierr
    """
    if dataFileName == None:
        dataFile = sys.stdin
    else:
        dataFile = open(dataFileName,'r')
    
    # saxs data lines
    q=[]
    I=[]
    iErr=[]
    for line in dataFile.readlines():
        line = line[:-1]    # Chop newline character
        fields = line.split()     # Split line into fields
        if ((len(fields) ==numbCols) & (fields[0][0] != "#")):  # check that there are  
       # the specified number of fields and that the first field does not begin with a #
            try:
                thisQ = float(fields[0])
                thisI = float(fields[1])
                thisErr = float(fields[2])
                q.append(thisQ)
                I.append(thisI)
                iErr.append(thisErr)
            except ValueError:
                continue
            
    q = np.array(q)
    I = np.array(I)
    iErr = np.array(iErr)

    dataFile.close()
    
    return [q, I, iErr]

def writeQIE(outFileName, data, noErr=False):
    if outFileName != None:
        outFile = open(outFileName,'w')
    else: 
        outFile = None


    if noErr == False:  
        print >> outFile, '#Q \tI \tIerror'
    else:
        print >> outFile, '#Q \tI'
        
    for i in range(len(data[0])):
        if noErr == False:
            print >> outFile, '%-8.5e \t %8.5e \t %8.5e' % (data[0][i], data[1][i], data[2][i])
        else:
            print >> outFile, '%-8.5e \t %8.5e' % (data[0][i], data[1][i])

    if outFile != None:
        outFile.close()

def readQI(dataFileName):
    """Reads generic SAXS data file for Q and I, no error field.
    This is convenient for reading simulated data files
    Selects lines that do not begin with #
    and contain either two or three fields with floating point numbers, 
     Returns 2 lists, containing Q and I
    """
    if dataFileName == None:
        dataFile = sys.stdin
    else:
        dataFile = open(dataFileName,'r')
    
    # saxs data lines
    q=[]
    I=[]

    for line in dataFile.readlines():
        line = line[:-1]    # Chop newline character
        fields = line.split()     # Split line into fields
        if (fields[0][0] != "#"):  # check that the first field does not begin with a #
            if ((len(fields) ==3 or len(fields) ==2)): # check that there are two or three fields            
                try:
                    thisQ = float(fields[0])
                    thisI = float(fields[1])
                    q.append(thisQ)
                    I.append(thisI)
                except ValueError:
                    continue
             
    q = np.array(q)
    I = np.array(I)

    dataFile.close()
    
    return [q, I]

def fileNamer(options, args, argNo, label):
    """ A function to handle file name options in the usToo scripts.
    inputs are:
    options: the options object returned by the OptionParser and is assumed to include
    stdin, stdout and raw.
    args: the argument list returned by the OptionParser
    If stdin is not used, the input file name is assumed to be the last element of args.
    argNo: the expected number of arguments, assuming that stdin is *not* used.
    label: a string that is appended to the root of the input file name 
    to make the output file name.
    Returns None if the number of arguments is not correct.
    Otherwise, returns rootname. input and output file names.
    If standard input or output are to be used, filenames are set to None
    If raw file option is used, output file name has .txt extension.
    Otherwise, output filename has .pdh extension.
    """
    rootName = ''
    if options.stdin:
        if len(args) !=argNo-1:
            return None
        else:
            inFileName = None
            outFileName = None           
    else:
        if (len(args) != argNo):        
            return None
        else:
            inFileName = args[-1]  # input file name is always the last argument!     
            #### construct output file name, unless standard output option is selected
            if options.stdo:
                outFileName = None
            else:
                rootName = inFileName.rpartition('.')[0]
                if options.rawFile:
                    outFileName = rootName + '_' + label + '.txt'
                else:
                    outFileName = rootName + '_' + label + '.pdh'
                    
    return rootName,inFileName,outFileName


############## Slit-smearing Functions  #######################

def setWparams(pdhHeader):
    """extract smear weighting parameters from pdh header, as read using
    readPdh(pdhFileName)"""

    line3fields=pdhHeader[2]
    line5fields= pdhHeader[4]
        
    type = line3fields[1]
    if type == 1:
        beamProfType = 'sig'
    elif type ==2:
        beamProfType = 'trap'
    else:
        beamProfType = 0
    # fifth line has the beam parameters we need
    a = line5fields[0]
    b = line5fields[1]
    bhw = line5fields[2]
    Ld = line5fields[3]
        
    return [beamProfType, a, b, bhw, Ld]


def w(y,wParams):
    """Calculates weighting factors for smearing due to beam length
    and detector slit width. First parameter is y value for which weight is to
    to be calculated. Second parameter is a list of parameters.
    #   wParam[0] = beam length profile type ('trap' or 'sig')
    wParam[1] = beam lenght profile parameter a
    wParam[2] = beam length profile parameter b
    wParam[3] = beam width parameter (Gaussian half-width at half height)
    wParam[4] = detector slit length (integration rectangle width for image
    plates)"""
    
    # all weighting functions are assumed to be symmetrical about y=0
    y = abs(y)
    
    a = wParams[1]
    b = wParams[2]
    Ld = wParams[4]

            
    if (wParams[0]=="trap"):
        # Trapezoidal beam profile function
        # Beam profile is defined by parameters a and b
        # a is the long side of the trapezoid, and b is short side
        # Ld is the detector slit length
        
        if (Ld==0):
            #if Ld=0, weights are simply given by the beam profile
            if (y > a/2):
                weight = 0
            elif (y>b/2):
                weight = 1-(y-b/2)/(a/2-b/2)
            else:
                weight = 1
        elif (a==0):
            # point collimation
            if (y>Ld/2):
                weight = 0
            else:
                weight = 1
        else:
            weight = trapBeamIntegral(y-Ld/2,y+Ld/2,a,b)
        
        
    if (wParams[0] =="sig"):
        # Sigmiodal model of beam profile
        # beam profile parameters are a and b
        # plus Ld (detector slit length)
        # Beam profile function: b(y)=(1+a)/(a+exp(abs(y)/b))
        
        if (Ld==0):
            #if Ld=0, weights are simply given by the beam profile 
            weight = (1+a)/(a+np.exp(y/b))
        
        elif ((y-Ld/2) < 0):
            weight = sigBeamIntegral(0,Ld/2-y,a,b) + sigBeamIntegral(0,y+Ld/2,a,b)
        else:
            weight = sigBeamIntegral(y-Ld/2,y+Ld/2,a,b)
        
    return weight

def trapBeamIntegral(y1,y2,a,b):
    """Calculates the definite integral of the trapezoidal beam profile
    function betwee y1 and y2. Used to calculate weighting function."""
    
    if (y1 <= -a/2):
        if (y2 <= -a/2):
            integral = 0
        elif (y2 <= -b/2):
            integral = (a+2*y2)**2/(4*(a-b))
        elif (y2 <= b/2):
            integral = (a-b)/4
            integral += b/2 + y2
        elif (y2 <= a/2):
            integral = (a-b)/4
            integral += b
            integral += -(2*a-b-2*y2)*(b-2*y2)/(4*(a-b))
        else:
            integral = (a+b)/2
    elif (y1 <= -b/2):
        if (y2 <= -b/2):
            integral = (y2-y1)*(a+y1+y2)/(a-b)
        elif (y2 <= b/2):
            integral = (b-2*y1)*(2*a+b+2*y1)/(4*(a-b))
            integral += b/2 + y2
        elif (y2 <= a/2):
            integral = (b-2*y1)*(2*a+b+2*y1)/(4*(a-b))
            integral += b
            integral += -(2*a-b-2*y2)*(b-2*y2)/(4*(a-b))
        else:
            integral = (b-2*y1)*(2*a+b+2*y1)/(4*(a-b))
            integral += b
            integral += (a-b)/4
    elif (y1 <= b/2):
        if (y2 <= b/2):
            integral = y2-y1
        elif (y2 <= a/2):
            integral = b/2-y1
            integral += -(2*a-b-2*y2)*(b-2*y2)/(4*(a-b))
        else:
            integral = b/2-y1
            integral += (a-b)/4
    elif (y1 <= a/2):
        if (y2 <= a/2):
            integral = (y1-y2)*(y1+y2-a)/(a-b)
        else:
            integral = (a-2*y1)**2/(4*(a-b))
    else:
        integral = 0
    
    return integral
            

def sigBeamIntegral(y1,y2,a,b):
    """Calculates the definite integral of the sigmoidal beam profile
    function betwee y1 and y2. Used to calculate weighting function. 
    This is only valid for y1>=0 and y2>=0"""

    upperInt = (1+a)*(y2-b*np.log(a+np.exp(y2/b)))/a
    lowerInt = (1+a)*(y1-b*np.log(a+np.exp(y1/b)))/a

    return upperInt-lowerInt


def wG(x,bhw):
    """Calculates weighting factor for beam width, using a Gaussian
    function to describe the width profile.  First parameter is the x value.
    Second parameter is the half-width (at half-height) that defines the Gaussian.
    returns normalized value"""
    
    if (bhw==0):
        return 1.0
    
    sigma = bhw/np.sqrt(np.log(2.0))
    normFact = np.sqrt(np.pi)*sigma
    
    return np.exp(-x**2/sigma**2)/normFact


def wNormFact(wParams,wym):
    """calculates integral of unnormalized weighting factor
    over range 0 to wym. Parameters are defined as in w(y,wParams)"""   
    
    a=wParams[1]
    Ld=wParams[4]
    
    if (a==0 and Ld==0):
        return 1

    
    deltaY=wym/40.0
    yArray  = np.arange(0,wym+deltaY,deltaY)
    integrand = []
    for y in yArray:
        integrand.append(w(y,wParams))
        
    integrand=np.array(integrand)
    wIntegral = integrate.simps(integrand,yArray)
    
    return wIntegral


def wYmax(wParams,eps):
    """Find the value for y such that w(y)/w(0)<= eps"""
    
    a=wParams[1]
    Ld=wParams[4]
    
    if (a==0 and Ld==0):
        return 0
        
    w0=w(0,wParams)
    deltaX=0.001
    y=0
    wY=w(y,wParams)
    while (wY > w0*eps):
        y += deltaX
        wY=w(y,wParams)
    
    return y
    
def wXmax(bhw,eps):
    """Find the value for x such that w(x)/w(0)= eps.
    Assumes Gaussian function"""

    if bhw==0:
        return 0
    
    sigma = bhw/np.sqrt(np.log(2.0))
    return sigma*np.sqrt(-np.log(eps))

def smearPt(q,tckUnsm,wParams,wym,wxm):
    """Calculates smeared intensity for a single value of q
       assumes as input a value q and the tck parameters
       defining a cubic spline representing the unsmeared scattering curve,
       the three parameters that define the smearing, the max y value for w(y),
       the beam half-width in x-direction and the max x value for w(x) """
    
    # This is the heart of things!
    # Integrates all of the contributions to the observed scattering
    # intensity at a nominal q-value
    # Contributions of scattering from other q-values are determined by
    # the beam geometry and the detector slit width.
    # x refers to the beam width (short dimension)
    # y refers to the beam length (long dimension)
    # integrate over y-axis up to wym
    # integrate over x-axis up to wxm

    bhw = wParams[3]
    
    if (wxm==0):
        # don't integrate over x
        if (wym ==0):
            # there is no smearing! Just return interpolated value
            integral = interpolate.splev(q,tckUnsm,der=0)
        else:
            # integrate over y
            yMax = wym
            deltaY = wym/20
            yArray = np.arange(0,yMax+deltaY,deltaY)
            yIntegrand=[]
            for y in yArray:    
                iVal = interpolate.splev(np.sqrt(q**2 + y**2),tckUnsm,der=0)
                wt = w(y,wParams)
                yIntegrand.append(wt*iVal)
            yIntegrand=np.array(yIntegrand)
            if type(q) == type(np.ndarray(shape=())):
                yIntegrand = np.transpose(yIntegrand)
                integral = []
                for i,qVal in enumerate(q):
                    integral.append(integrate.simps(yIntegrand[i],yArray))
                integral=np.array(integral)            
            else:
                integral = integrate.simps(yIntegrand,yArray)

                    
    else:
        # integrate over x
        if (wym ==0):
            # integrate only over x
            xMin=-wxm
            xMax= wxm
            deltaX = xMax/10
            xArray = np.arange(xMin,xMax+deltaX,deltaX)
            xIntegrand=[]
            for x in xArray:
                xWt = wG(x,bhw)
                iVal = interpolate.splev(q+x,tckUnsm,der=0)
                xIntegrand.append(xWt*iVal)
            xIntegrand = np.array(xIntegrand)
            if type(q) == type(np.ndarray(shape=())):
                xIntegrand = np.transpose(xIntegrand)
                integral = []
                for i,qVal in enumerate(q):
                    integral.append(integrate.simps(xIntegrand[i],yArray))
                integral=np.array(integral)            
            else:
                integral = integrate.simps(xIntegrand,xArray)
        else:
            # integrate over both x and y
            xMin=-wxm
            xMax= wxm
            deltaX = xMax/10
            xArray = np.arange(xMin,xMax+deltaX,deltaX)
            yMax = wym
            deltaY = wym/40
            yArray = np.arange(0,yMax+deltaY,deltaY)
            if type(q) == type(np.ndarray(shape=())):  
                integral = []
                for i,qVal in enumerate(q):
                    xIntegrand=[]
                    for x in xArray:
                        xWt = wG(x,bhw)
                        yIntegrand=[]
                        for y in yArray:
                            iVal = interpolate.splev(np.sqrt((qVal+x)**2 + y**2),tckUnsm,der=0)
                            wt = w(y,wParams)*xWt
                            yIntegrand.append(wt*iVal)
                        yIntegrand = np.array(yIntegrand)
                        yIntegral = integrate.simps(yIntegrand,yArray)  
                        xIntegrand.append(yIntegral)
                    xIntegrand = np.array(xIntegrand)
                    integral.append(integrate.simps(xIntegrand,xArray))
                    
                
            else:
                xIntegrand=[]
                for x in xArray:
                    xWt = wG(x,bhw)
                    yIntegrand=[]
                    for y in yArray:
                        iVal = interpolate.splev(np.sqrt((q+x)**2 + y**2),tckUnsm,der=0)
                        wt = w(y,wParams)*xWt
                        yIntegrand.append(wt*iVal)
                    yIntegrand = np.array(yIntegrand)
                    yIntegral = integrate.simps(yIntegrand,yArray)  
                    xIntegrand.append(yIntegral)
                xIntegrand = np.array(xIntegrand)
                integral = integrate.simps(xIntegrand,xArray)

    return integral                   


######### Data Plotting Functions  ########################


def saxsPlot(fig,plotPos, plotParam, scattData, legends=[],colors =[], qUnits=0, absI=False):
    """Creates MatPlotLib axis and plots one or more saxs profiles. 
    Input parameters are:
    fig: MatPlotLib figure, usually created with something like: fig = plt.figure()
    plotPos: A three digit number representing position of subplot in the figure.
        The first digit is the number of rows of subplots in the figure.
        The second digit is the number of columns of subplots in the figure.
        The third digit is position of the new subplot in the figure.
        The figure positions fill the first row, left to right, then the 
        second row, left to right and so on.
    plotParam: list of plot parameters: [plotType, qPlotLim,iPlotLim]
        plotType: one of 'linear', 'log', 'loglog', 'guinier' or 'kratky'
        qlim:  q plot limits - either list of two floats or other data type (eg 0)
        ilim: i plot limits  - either list of two floats or other data type (eg 0)
        qlim and ilim are type checked and only used if they are lists of two floats
    scattData: list of one or more data sets.  Each data set is a list of 
        two or three lists of floats: q, I and (optional) I error
    legends: optional list of legends for scattering data, as text strings
    colors: optional list of colors for individual curves
    qUnits: optional specification for q-units to be indicated in x-axis label, 
        default is no units specified. 1 = A-1, 10 = nm-1, 0 = none
    absolute intensity calibration, boolean. Defaut is False
    Returns axis object.
    """

    plotType= plotParam[0]
    qlim = plotParam[1]
    ilim = plotParam[2]
            

    ax = fig.add_subplot(plotPos)
    ax.autoscale(enable=True,tight=True)
    if len(colors)>0:    
        ax.set_color_cycle(colors)
    
    if qUnits ==1:
        qUnitsStr = " ($\AA^{-1}$)"
        qSqUnitsStr = " ($\AA^{-2}$)"
    elif qUnits ==10:
        qUnitsStr = " (nm$^{-1}$)"
        qSqUnitsStr = " (nm${-2}$)"
    else:
        qUnitsStr = ""
        qSqUnitsStr = ""



        
    for i,dataSet in enumerate(scattData):
        
        if i<len(legends):
            label = legends[i]
        else:
            label = ''
        if plotType == 'linear':
            if len(dataSet)== 3:
                ax.errorbar(dataSet[0],dataSet[1],yerr=dataSet[2],
                    linestyle='None', marker='.',capsize=0,ms=1,label=label)
            else:
                ax.plot(dataSet[0],dataSet[1], label=label)
            if type(ilim) == types.ListType and len(ilim) ==2:
                ax.set_ylim(ilim)
            if type(qlim) == types.ListType and len(qlim) ==2:
                ax.set_xlim(qlim)
            ax.set_xlabel('$q$' + qUnitsStr)
            ax.set_ylabel('$I$')

        if plotType == 'log':
            if len(dataSet)== 3:
                ax.errorbar(dataSet[0],dataSet[1],yerr=dataSet[2],
                    linestyle='None', marker='.',capsize=0,ms=1, label=label)
            else:
                ax.plot(dataSet[0],dataSet[1], label=label)
            ax.set_yscale('log')
            if type(ilim) == types.ListType and len(ilim) ==2:
                ax.set_ylim(ilim)
            if type(qlim) == types.ListType and len(qlim) ==2:
                ax.set_xlim(qlim)
            
            ax.set_xlabel('$q$' + qUnitsStr)
            ax.set_ylabel('$I$')
            
        if plotType == 'loglog':        
            if len(dataSet)== 3:
                ax.errorbar(dataSet[0],dataSet[1],yerr=dataSet[2],
                    linestyle='None', marker='.',capsize=0,ms=1, label=label)
            else:
                ax.plot(dataSet[0],dataSet[1], label=label)
            ax.set_xscale('log')
            ax.set_yscale('log')
            if type(ilim) == types.ListType and len(ilim) ==2:
                ax.set_ylim(ilim)
            if type(qlim) == types.ListType and len(qlim) ==2:
                ax.set_xlim(qlim)
            
            ax.set_xlabel('$q$' + qUnitsStr)
            ax.set_ylabel('$I$')
        
        if plotType =='guinier':
            plotQ=np.square(dataSet[0])
            plotI=np.log(dataSet[1])
            if len(dataSet)== 3:
                plotErr = np.log((dataSet[1]+dataSet[2])/dataSet[1])
                ax.errorbar(plotQ,plotI,yerr=plotErr,linestyle='None', 
                    marker='.',capsize=0,ms=1, label=label)
            else:
                ax.plot(plotQ,plotI, label=label)
            if type(ilim) == types.ListType and len(ilim) ==2:
                ax.set_ylim([np.log(ilim[0]),np.log(ilim[1])])
            if type(qlim) == types.ListType and len(qlim) ==2:
                if qlim[0]<0:
                    qlim[0]=0
                ax.set_xlim([np.square(qlim[0]),np.square(qlim[1])])
            ax.set_xlabel('$q^2$' + qSqUnitsStr)
            ax.set_ylabel('$\ln(I)$')
        
        if plotType == 'kratky':
            plotI=dataSet[1]*np.square(dataSet[0])
            plotQ = dataSet[0]
            if len(dataSet)== 3:
                plotErr = dataSet[2]*np.square(dataSet[0])
                ax.errorbar(plotQ,plotI,yerr=plotErr,linestyle='None',
                    marker='.',capsize=0,ms=1, label=label)
            else:
                ax.plot(plotQ,plotI, label=label)
            if type(qlim) == types.ListType and len(qlim) ==2:
                ax.set_xlim(qlim)
        
            ax.set_xlabel('$q$' + qUnitsStr)
            ax.set_ylabel('$Iq^2$')

    if len(legends)>0:
        leg = ax.legend(markerscale=10,scatterpoints=1, numpoints=1) 
        for t in leg.get_texts():
            t.set_fontsize('small')    # the legend text fontsize
            
    if absI:
        if plotType == 'kratky':
            ax.set_ylabel('$Iq^2$ (cm$^{-1} \AA^{-2}$)')
        else:        
            label = ax.get_ylabel()
            ax.set_ylabel(label + ' (cm$^{-1}$)')
    
    return ax
    


def plotBeamLengthProfile(fig,plotPos,wParams,wym):
    """Creates MatPlotLib axis and plots beam length profile specified by wParams
    and the weighting function.  The two curves are identical if the detector slit length
    is 0.  Otherwise, the detector length further smears the weighting function.
    Input parameters are:
    fig: MatPlotLib figure, usually created with something like: fig = plt.figure()
    plotPos: A three digit number representing position of subplot in the figure.
        The first digit is the number of rows of subplots in the figure.
        The second digit is the number of columns of subplots in the figure.
        The third digit is position of the new subplot in the figure.
        The figure positions fill the first row, left to right, then the 
        second row, left to right and so on.
    wParams: list of beam parameters, as defined in w(y,wParams)
    wym: maximum value of wym to plot
    """
    deltaY = wym/100
    yArray = np.arange(-wym,wym+deltaY,deltaY)
    wY=[]
    bY=[]
    bParams= [wParams[0],wParams[1],wParams[2],wParams[3],0]

    bnf = wNormFact(bParams,wym)
    wnf = wNormFact(wParams,wym)
    for y in yArray:
        bY.append(w(abs(y),bParams)/bnf)
        wY.append(w(abs(y),wParams)/wnf)     
    wY =np.array(wY)
    bY = np.array(bY)

    ax=fig.add_subplot(plotPos)
    
    ax.plot(yArray,bY,"b-")
    ax.plot(yArray,wY,"r-")
    
    ax.set_xlabel('y')
    ax.set_ylabel('w(y)')
    
    leg=ax.legend(('beam length profile','weighting function'),loc='lower center')
    ax.set_title("Beam Length Profile and Weighting Function")
    ax.set_xlim([-wym,wym])
    
    for t in leg.get_texts():
        t.set_fontsize('small')    # the legend text fontsize
        
    text = "Beam Profile: \nType ="
    if wParams[0]=='trap':
        text += 'trapezoidal \n'
    elif wParams[0] == 'sig':
        text += 'sigmoidal \n'
    text += 'a = ' + str(wParams[1]) + '\n'
    text += 'b = ' + str(wParams[2])
    
    plt.text(0.05, 0.9, text,
    horizontalalignment='left',
    verticalalignment='top',
    transform = ax.transAxes,size='small')
    
    text = 'Detector slit-length = ' + str(wParams[4])
    plt.text(0.95, 0.9, text,
    horizontalalignment='right',
    verticalalignment='top',
    transform = ax.transAxes,size='small')
    
    return ax

def plotBeamWidthProfile(fig,plotPos,wParams,wxm,qRange):
    """Creates MatPlotLib axis and plots width profile specified by wParams
    and the weighting function. 
    Input parameters are:
    fig: MatPlotLib figure, usually created with something like: fig = plt.figure()
    plotPos: A three digit number representing position of subplot in the figure.
        The first digit is the number of rows of subplots in the figure.
        The second digit is the number of columns of subplots in the figure.
        The third digit is position of the new subplot in the figure.
        The figure positions fill the first row, left to right, then the 
        second row, left to right and so on.
    wParams: list of beam parameters, as defined in w(y,wParams)
    wym: maximum value of wym to plot
    qRange: range of q values for x-axis.  
    """
    
    bhw = wParams[3]
    
    deltaX=wxm/100
    xArray = np.arange(-wxm,wxm+deltaX,deltaX)
    wX =[]
    for x in xArray:
        wX.append(wG(x,bhw))
    wX=np.array(wX)/wG(0,bhw)
    
    ax=fig.add_subplot(plotPos)    

    ax.plot(xArray,wX)
    ax.set_xlabel('x')
    ax.set_ylabel('w(x)')    
    
    ax.set_title("Beam Width Profile")
    
    ax.set_ylim([0.0,1.1])
    if qRange>wxm:
        ax.set_xlim([-qRange,qRange])
        
    text = "Beam half-width = " + str(bhw)
    
    plt.text(0.05, 0.85, text,
    transform = ax.transAxes,size='small')

    return ax


######################### Smoothing Function ########################

def smooth(x,window_len=11,window='hanning'):
    """smooth the data using a window with requested size.
    
    From SciPy Cookbook: http://www.scipy.org/Cookbook/SignalSmooth
    
    This method is based on the convolution of a scaled window with the signal.
    The signal is prepared by introducing reflected copies of the signal 
    (with the window size) in both ends so that transient parts are minimized
    in the begining and end part of the output signal.
    
    input:
        x: the input signal 
        window_len: the dimension of the smoothing window; should be an odd integer
        window: the type of window from 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'
            flat window will produce a moving average smoothing.

    output:
        the smoothed signal
        
    example:

    t=linspace(-2,2,0.1)
    x=sin(t)+randn(len(t))*0.1
    y=smooth(x)
    
    see also: 
    
    numpy.hanning, numpy.hamming, numpy.bartlett, numpy.blackman, numpy.convolve
    scipy.signal.lfilter
 
    TODO: the window parameter could be the window itself if an array instead of a string   
    """

    if x.ndim != 1:
        raise ValueError, "smooth only accepts 1 dimension arrays."

    if x.size < window_len:
        raise ValueError, "Input vector needs to be bigger than window size."


    if window_len<3:
        return x


    if not window in ['flat', 'hanning', 'hamming', 'bartlett', 'blackman']:
        raise ValueError, "Window is not one of 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'"


    s=np.r_[2*x[0]-x[window_len:1:-1],x,2*x[-1]-x[-1:-window_len:-1]]
    #print(len(s))
    if window == 'flat': #moving average
        w=np.ones(window_len,'d')
    else:
        w=eval('np.'+window+'(window_len)')

    y=np.convolve(w/w.sum(),s,mode='same')
    return y[window_len-1:-window_len+1]


############## Physical Properties of Water  ######################

def waterCompr(t):
    """Calculates isothermal compressiblity of water as a function of 
    temperature, which is used for calibration of absolute scattering
    intensities from reference watering scattering intensity.  
    Input is temperature, in Celsius.
    Output is compressibility in Pa-1
    Uses equation 20 from:  Krell, G. S. (1975). Density, thermal expansivity, 
    and compressibility of liquid water from 0 to 100 C: Correlations and 
    tables for atmospheric pressure and stauration reviewed and expressed 
    on 1968 temperature scale. J. Chem. Eng. Data, 20, 97-105. 
    http://dx.doi.org/10.1021/je60064a005"""
    
    a = 50.88496
    b = 0.616383
    c = 1.459187E-3
    d = 20.08438E-6
    e = -58.47727E-9
    f = 410.4110E-12
    g = 19.67348E-3
    
    comp = 1.0E-11*(a + b*t + c*t**2 + d*t**3 + e*t**4 + f*t**5)/(1.0 + g*t)
    return comp
    
    
def waterRho(t):
    """Calculates density of water as a function of 
    temperature, which is used for calibration of absolute scattering
    intensities from reference watering scattering intensity.  
    Input is temperature, in Celsius.
    Output is density in g/cm^3
    Uses equation 16 from:  Krell, G. S. (1975). Density, thermal expansivity, 
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
    
    rho = 1.0E-3*(a + b*t + c*t**2 + d*t**3 + e*t**4 + f*t**5)/(1.0 + g*t)
    return rho
    
    
def waterS0(t):
    """Calculates structure factor of liquid water at 0 scattering angle.
    Input is temperature in Celsius.  Output is dimensionless structure factor
    Equation from http://physchem.kfunigraz.ac.at/sm/"""
    
    k = 1.380658E-23 # Boltzmann constant in J/K
    Mr = 18.0151528  # molecular weight of water in g/mol
    Na = 6.0221367E23 # Avogadro's number
    
    rho = waterRho(t) # density of water in g/cm^3
    comp = waterCompr(t) # compressibility of water in Pa-1
    
    T = t+273.15 # absolute temperature
    
    s0 = 1E6*(Na*rho/Mr)*k*T*comp
    return s0
    
def waterI0(t):
    """Calculates absolute scattering intensity of liquid water at 0 scattering angle.
    Input is temperature in Celsius.  Output is scattering intensity in cm-1
    Equation from http://physchem.kfunigraz.ac.at/sm/"""
    
    Mr = 18.0151528  # molecular weight of water in g/mol
    Na = 6.0221367E23 # Avogadro's number
    ePerWater = 10 # number of electrons per water molecule
    eScattLen = 0.28179E-12 # Scattering length of an electron
    rho = waterRho(t) # density of water in g/cm^3
    s0 = waterS0(t)
    
    i0 = (ePerWater*eScattLen)**2*s0*Na*rho/Mr
    return i0
    
def rhoEwat(t):
    """Calculates scattering density of water from temp, using mass densitity"""
    Mrwat = 18.0151528  # molecular weight of water in g/mol
    Na = 6.0221367E23 # Avogadro's number
    ePerWater = 10 # number of electrons per water molecule
    eScattLen = 0.28179E-12 # Scattering length of an electron, cm
    
    rhoWat = waterRho(t) # density of water in g/cm^3
    rhoEwat = (rhoWat/Mrwat)*Na*ePerWater*eScattLen # water scattering density (cm-2)
    return rhoEwat
    
    
    
######## Calculation of molecular weight from absolute I0 and other stuff ###########
    
def mw(i0,c, vbar,deltaRho):
    """Calculates molecular weight from absolute X-ray scattering 
    intensity at 0 scattering angle.
    Inputs are scattering intensity (i0, cm-1), concentration in mg/cm^3
    partial specific volume of particle (vbar, cm^3/g)
    number of electrons per unit molecular weight.
    Outputs molecular weight in g/mol"""

    Na = 6.0221367E23 # Avogadro's number
    c = c*1E-3  # convert mg/mL conc to g/cm^3


    m = i0*Na/(c*vbar**2*deltaRho**2)
    return m

def conc(i0,mw, vbar,deltaRho):
    """Calculates concentration from absolute X-ray scattering 
    intensity at 0 scattering angle.
    Inputs are scattering intensity (i0, cm-1), relative molecular mass
    partial specific volume of particle (vbar, cm^3/g)
    scattering density difference of molecule and solvent"""

    Na = 6.0221367E23 # Avogadro's number
    
    c = i0*Na/(mw*vbar**2*deltaRho**2)
    
    c = c*1E3  # convert g/mL conc to mg/cm^3
    return c
    
    
def hsStructF(q,phi,r):
    """Function to calculate theoretical structure factor for solution of hard-sphere particles
    Based on Percus-Yevick approximation, as given by Boualem Hammouda in the SANS Toolbox
    http://www.ncnr.nist.gov/staff/hammouda/the_SANS_toolbox.pdf
    input parameters are scattering momentum (q), volume fraction density of particles (phi)
    and radius of particles"""
    qd = 2.0*r*q
    if (qd <= 0.0):
        return (phi-1.0)**4/(1.0+2.0*phi)**2

    sin = np.sin(qd)
    cos = np.cos(qd)

    alpha = (1+2.0*phi)**2/(1.0-phi)**4
    beta = -6.0*phi*(1.0+phi/2.0)**2/(1.0-phi)**4
    gamma = (phi/2.0)*alpha
    
    term1 = sin-qd*cos
    term2 = 2.0*qd*sin - (qd**2-2.0)*cos-2.0
    term3 = (4.0*qd**3-24.0*qd)*sin - (qd**4-12.0*qd**2 +24.0)*cos+24.0
    
    c = -(24.0*phi/qd**6)*(alpha*qd**3*term1 + beta*qd**2*term2 + gamma*term3)
    s = 1.0/(1.0-c)
    
    return s

# This vectorizes the structure factor function so that it works with numpy arrays
hsStructFact = np.vectorize(hsStructF)


def splitThousands(s, tSep, dSep=None):
    """A little function to format string representations of numbers with 
    thousands and decimal separators.  Found on http://code.activestate.com/recipes/498181/"""
    if s.rfind('.')>0:
        rhs=s[s.rfind('.')+1:]
        s=s[:s.rfind('.')-1]
        if len(s) <= 3: return s + dSep + rhs
        return splitThousands(s[:-3], tSep) + tSep + s[-3:] + dSep + rhs
    else:
        if len(s) <= 3: return s
        return splitThousands(s[:-3], tSep) + tSep + s[-3:]
