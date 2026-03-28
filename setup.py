#!/usr/bin/env python3

from setuptools import setup
from setuptools import Extension
from setuptools import find_namespace_packages
from sys import version_info as vi

installed_version = (vi[0], vi[1])

if installed_version < (3, 8):
    raise Exception(
        "The default Python version must be 3.8 or higher, not {0}.{1}".format(vi[0], vi[1]))

USE_CYTHON = True

STREAMNODE_PYX = "src/heatsource9/model/streamnode.pyx"
PYHEATSOURCE_PYX = "src/heatsource9/model/pyheatsource.pyx"

extensions = [
    Extension(
        name="heatsource9.model.streamnode",
        sources=[STREAMNODE_PYX],
    ),
    Extension(
        name="heatsource9.model.pyheatsource",
        sources=[PYHEATSOURCE_PYX],
    ),
]

if USE_CYTHON:
    from Cython.Build import cythonize

    extensions = cythonize(
        extensions,
        compiler_directives={"language_level": "3"},
    )

setup(
    name="heatsource9",
    version="9.0.0b29",
    classifiers=[
        "Development Status :: 4 - Beta",
        "Environment :: Console",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
        "Natural Language :: English",
        "Operating System :: MacOS :: MacOS X",
        "Operating System :: POSIX :: Linux",
        "Operating System :: Microsoft :: Windows",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3 :: Only",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
        "Programming Language :: Python :: 3.12",
        "Programming Language :: Python :: 3.13",
        "Programming Language :: Python :: 3.14",
        "Topic :: Scientific/Engineering",
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
    description="One-dimensional stream temperature modeling program",
    url="https://www.oregon.gov/deq/wq/tmdls/Pages/TMDLs-Tools.aspx",
    project_urls={
        "Documentation": "https://www.oregon.gov/deq/FilterDocs/heatsourcemanual.pdf",
        "Source": "https://github.com/DEQrmichie/heatsource-9/",
    },
    author="Matt Boyd, Brian Kasper, Terra Metta, Ryan Michie, Dan Turner",
    maintainer="Ryan Michie, Oregon DEQ",
    maintainer_email="ryan.michie@deq.oregon.gov",
    platforms=["darwin", "linux", "win32"],
    license="GNU General Public License v3 (GPLv3)",
    zip_safe=False,
    entry_points={"console_scripts": ["hs = heatsource9.__main__:main"]},
    packages=find_namespace_packages("src"),
    package_dir={"": "src"},
    install_requires=["Cython==3.1.4", "openpyxl <= 3.1.5"],
    ext_modules=extensions,
    python_requires=">=3.8, <4",
)
