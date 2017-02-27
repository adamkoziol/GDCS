#!/usr/bin/env python
from setuptools import setup, find_packages
__author__ = 'adamkoziol'

setup(
    name="genesipprOOP",
    version="0.1",
    packages=find_packages(),
    include_package_data=True,
    license='MIT',
    author='Adam Koziol',
    author_email='adam.koziol@inspection.gc.ca',
    description='Object oriented raw read typing software',
    url='https://github.com/OLC-LOC-Bioinformatics/geneSipprV2/tree/dev',
    long_description=open('README.md').read(),
    install_requires=['biopython >= 1.65',
                      'pysam',
                      'pysamstats',
                      'xlsxwriter'],
)
