#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""The setup script."""

from setuptools import setup, find_packages

with open('README.rst') as readme_file:
    readme = readme_file.read()

requirements = ['numpy', 'scipy', 'sympy', 'python-libsbml', 'pyscipopt', 'pandas']

setup_requirements = []

test_requirements = requirements + ['cplex']

setup(
    author="Daniel Machado",
    author_email='cdanielmachado@gmail.com',
    python_requires='>=3.8',
    classifiers=[
        'License :: OSI Approved :: Apache Software License',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
        'Programming Language :: Python :: 3.10',
        'Programming Language :: Python :: 3.11',
        'Programming Language :: Python :: 3.12',
        'Programming Language :: Python :: 3.13',
    ],
    description="metabolic modeling package",
    install_requires=requirements,
    license="Apache Software License 2.0",
    long_description=readme,
    include_package_data=True,
    keywords='reframed',
    name='reframed',
    packages=find_packages(),
    setup_requires=setup_requirements,
    test_suite='tests',
    tests_require=test_requirements,
    url='https://github.com/cdanielmachado/reframed',
    version='1.6.0',
    zip_safe=False,
    entry_points={ 
        "console_scripts": [
            "fba=reframed.cli.cli:main",
        ],
    },
)
