#!/usr/bin/env python3

from setuptools import setup
from setuptools import Extension
from sys import version_info as vi

installed_version = (vi[0], vi[1])

if installed_version < (3, 8):
    raise Exception("The default Python version must be 3.8 or higher, not {0}.{1}".format(vi[0], vi[1]))

USE_CYTHON = True

if USE_CYTHON:
    from Cython.Build import cythonize
    extensions = cythonize(['src/heatsource9/Stream/StreamNode.pyx',
                            'src/heatsource9/Stream/PyHeatsource.pyx'],
                           compiler_directives={'language_level': "3"})
else:
    extensions = [Extension('heatsource9.Stream.PyHeatsource', ['src/heatsource9/Stream/PyHeatsource.c']),
                  Extension('heatsource9.Stream.StreamNode', ['src/heatsource9/Stream/StreamNode.c'])]

setup(name='heatsource9',
      version='9.0.0b28',
      classifiers=[
          'Development Status :: 4 - Beta',
          'Environment :: Console',
          'Intended Audience :: Science/Research',
          'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
          'Natural Language :: English',
          'Operating System :: MacOS :: MacOS X',
          'Operating System :: POSIX :: Linux',
          'Operating System :: Microsoft :: Windows',
          'Programming Language :: Python :: 3',
          'Programming Language :: Python :: 3 :: Only',
          'Programming Language :: Python :: 3.8',
          'Programming Language :: Python :: 3.9',
          'Programming Language :: Python :: 3.10',
          'Programming Language :: Python :: 3.11',
          'Programming Language :: Python :: 3.12',
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
      url='https://www.oregon.gov/deq/wq/tmdls/Pages/TMDLs-Tools.aspx',
      project_urls={
          'Documentation': 'https://www.oregon.gov/deq/FilterDocs/heatsourcemanual.pdf',
          'Source': 'https://github.com/DEQrmichie/heatsource-9/'},
      author='Matt Boyd, Brian Kasper, John Metta, Ryan Michie, Dan Turner',
      maintainer='Ryan Michie, Oregon DEQ',
      maintainer_email='ryan.michie@deq.state.or.us',
      platforms=['darwin', 'linux', 'win32'],
      license='GNU General Public License v3 (GPLv3)',
      zip_safe=False,
      entry_points={'console_scripts': ['hs = heatsource9.BigRedButton:hs']},
      packages=['heatsource9',
                'heatsource9.ModelSetup',
                'heatsource9.Dieties',
                'heatsource9.Stream',
                'heatsource9.Utils'],
      package_dir={'': 'src'},
      install_requires=['Cython==3.0.10', 'openpyxl <= 3.1.3'],
      ext_modules=extensions,
      python_requires='>=3.8, <4'
      )
