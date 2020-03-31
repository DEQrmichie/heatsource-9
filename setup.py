#!/usr/bin/env python3

from setuptools import setup
from setuptools import Extension
from Cython.Build import cythonize
from sys import version_info as vi

installed_version = (vi[0], vi[1])

if installed_version < (3, 0):
    raise Exception("The default Python version must be 3.0 or higher, not {0}.{1}".format(vi[0], vi[1]))

setup(name='heatsource9',
      version='9.0.0b25',
      classifiers=[
          'Development Status :: 4 - Beta',
          'Environment :: Console',
          'Intended Audience :: Science/Research',
          'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
          'Natural Language :: English',
          'Operating System :: MacOS :: MacOS X',
          'Operating System :: Microsoft :: Windows',
          'Programming Language :: Python :: 3',
          'Topic :: Scientific/Engineering'
      ],
      long_description="""Heat Source is a computer model used by the
      Oregon Department of Environmental Quality to simulate stream
      thermodynamics and hydraulic routing. It was originally developed
      by Matt Boyd in 1996 as a Masters Thesis at Oregon State University
      in the Departments of Bioresource Engineering and Civil
      Engineering. Since then it has grown and changed significantly.
      Oregon DEQ currently maintains the Heat Source methodology and
      computer programming. Appropriate model use and application are
      the sole responsibility of the user.""",
      description='One-dimensional stream temperature modeling program',
      url='http://www.deq.state.or.us/wq/TMDLs/tools.htm',
      project_urls={
          'Documentation': 'https://www.oregon.gov/deq/FilterDocs/heatsourcemanual.pdf',
          'Source': 'https://github.com/rmichie/heatsource-9/'},
      author='Matt Boyd, Brian Kasper, John Metta, Ryan Michie, Dan Turner',
      maintainer='Ryan Michie, Oregon DEQ',
      maintainer_email='michie.ryan@deq.state.or.us',
      platforms=['darwin', 'win32'],
      license=['GNU General Public License v3 (GPLv3)'],
      zip_safe=False,
      entry_points={'console_scripts': ['hs = heatsource9.BigRedButton:hs']},
      packages=['heatsource9',
                'heatsource9.ModelSetup',
                'heatsource9.Dieties',
                'heatsource9.Stream',
                'heatsource9.Utils'],
      package_dir={'': 'src'},
      ext_modules=cythonize('src/heatsource9/Stream/*.pyx', compiler_directives={'language_level': "3"})
      )
