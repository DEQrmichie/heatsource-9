#!/usr/bin/env python

from distutils.core import setup, Extension
from sys import version_info as VI

DISTUTILS_DEBUG = True

if VI < (2,7):
    v = "%i.%i" %(VI[0],VI[1])
    raise Exception("Default Python version must be >2.7, not %s" % v)

setup(name='heatsource9',
      version='9.0.0b11',
      classifiers=[
          'Development Status :: 4 - Beta',
          'Environment :: Console',
          'Intended Audience :: End Users/Desktop',
          'Intended Audience :: Science/Research',
          'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
          'Natural Language :: English',
          'Operating System :: MacOS :: MacOS X',
          'Operating System :: Microsoft :: Windows',
          'Programming Language :: Python :: 2.7',
          'Topic :: Scientific/Engineering',
          ],
      long_description="""Heat Source is a computer model used by the Oregon Department of Environmental Quality to simulate
      stream thermodynamics and hydraulic routing. It was originally developed by Matt Boyd in 1996 as a Masters Thesis at Oregon
      State University in the Departments of Bioresource Engineering and Civil Engineering. Since then it has grown and changed 
      significantly. Oregon DEQ currently maintains the Heat Source methodology and computer programming. 
      Appropriate model use and application are the sole responsibility of the user.""",
      description='One-dimensional stream temperature modeling program',
      url='http://www.deq.state.or.us/wq/TMDLs/tools.htm',
      author='Matt Boyd, Brian Kasper, John Metta, Ryan Michie, Dan Turner',
      maintainer='Ryan Michie, Oregon DEQ',
      maintainer_email='michie.ryan@deq.state.or.us',
      platforms = ['darwin', 'win32'],
      license = ['GNU General Public License v3 (GPLv3)'],
      packages=['heatsource9',
                'heatsource9.CSV',
                'heatsource9.Dieties',
                'heatsource9.Stream',
                'heatsource9.Utils'])
