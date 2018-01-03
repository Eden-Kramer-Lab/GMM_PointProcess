#!/usr/bin/env python3

from setuptools import find_packages, setup

INSTALL_REQUIRES = ['numpy', 'scipy', 'matplotlib']
TESTS_REQUIRE = ['pytest >= 2.7.1']

setup(
    name='GMM_PointProcess',
    version='0.0.1.dev0',
    license='MIT',
    description=('Closed Form Solution for 2D Replay Decoder Using Gaussian '
                 'Mixture Models'),
    author='Ali Yousefi',
    author_email='AYOUSEFI@mgh.harvard.edu',
    url='https://github.com/Eden-Kramer-Lab/GMM_PointProcess',
    packages=find_packages(),
    install_requires=INSTALL_REQUIRES,
    tests_require=TESTS_REQUIRE,
)
