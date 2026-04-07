#!/usr/bin/env python3

from setuptools import setup
from setuptools import Extension
from setuptools import find_namespace_packages

from Cython.Build import cythonize

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

setup(
    packages=find_namespace_packages("src"),
    package_dir={"": "src"},
    ext_modules=cythonize(extensions,
                          compiler_directives={"language_level": 3},
                          ),
)
