4 November 2011
Installation notes for USToo SAXS programs

There are two components to USToo:

1. Macros for ImageJ.  

These can be placed any directory to which the user has read access. The macros are then loaded into ImageJ with the Plugins -> Macros -> Install command. But, it is most convenient to place them in the ImageJ/macros folder.  If the name of the file is changed to StartupMacros.txt and placed in this location, then the macros will be automatically loaded when ImageJ is launched.  It is probably a good idea to rename the original StartupMacros.txt file provided with ImageJ, to keep from overwriting it.

2. Python scripts. 

The Python scripts include those for the individual processing commands, and two library files, saxs.py and soln.py, that are used by the command scripts.  It is easiest if all of these files are kept in a single directory.

The first line of each script must begin with the characters "#!", followed by the path to the python executable.
As distributed, the scripts specify:

#! /Library/Frameworks/Python.framework/Versions/Current/bin/python

This is the location for OS X with the Enthought Python Distribution (EPD) installed.  The location may be different for other python installations or other systems.

On Unix or Unix-like systems the Python files can be placed anywhere that the user has access to.  For a single user, the easiest thing may be to place them in a directory called bin in the login directory.   If multiple users of the same system will be using the programs, it may be helpful to put them in a common directory, for instance a directory called usToo in /usr/local.  In either case, the users' path variable should then be set to include this directory on the search path.  The file permissions must be set to allow execution by the users.

For Windows, I'm afraid that I can't be of as much help,except to make the following notes:
1. The #! line is ignored under Windows.
2. The command interpreter uses the extensions of file names to determine the program that executes them.  As a consequence, the ".py" extension must be added the names of the script files, and they are executed by typing, for instance, saxsPlot.py at the command prompt.
3. The programs will work if they are in the working directory, though that is generally inconvenient.  
4. The path variable, %PATH%, can be set so that the programs can be called from other directories, but I'm not sure exactly how this is done in Windows!

############# Other Required Software  #############################

The ImageJ program can be obtained from: http://rsbweb.nih.gov/ij/
The macros work with versions up to 1.44o, but I have encountered problems with version 1.45 and 1.46.

The Python scripts require that the Python interpreter and supporting files be installed.  Python version 2.7 is recommended.  Do not use version 3! 
Python can be obtained at http://www.python.org/

In addition to the standard Python installation, the following packages are required:

1. scipy and numpy - http://www.scipy.org
2. matplotlib - http://sourceforge.net/projects/matplotlib/files/
3. uncertainties - http://pypi.python.org/pypi/uncertainties/

The easiest way to install everything (except uncertainties) is by using the Enthought Python Distribution, a very robust installation specifically intended for scientific computing.  It is available free for academic use at: http://www.enthought.com/products/epd.php.  The uncertainties package must be installed separately after installing Python.


#  (c) 2009-2012 by David P. Goldenberg
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
