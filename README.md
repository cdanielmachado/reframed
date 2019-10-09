[![License](https://img.shields.io/badge/License-Apache%202.0-blue.svg)](https://opensource.org/licenses/Apache-2.0)
[![PyPI version](https://badge.fury.io/py/reframed.svg)](https://badge.fury.io/py/reframed)
[![Build Status](https://travis-ci.org/cdanielmachado/reframed.svg?branch=master)](https://travis-ci.org/cdanielmachado/reframed)
[![Documentation Status](https://readthedocs.org/projects/reframed/badge/?version=latest)](https://reframed.readthedocs.io/en/latest/?badge=latest)
[![DOI](https://zenodo.org/badge/212059108.svg)](https://zenodo.org/badge/latestdoi/212059108)

![ReFramed](reframed_logo.png)

# ReFramed: metabolic modeling package

**ReFramed** implements many constraint-based simulation methods (see list below), and contains interfaces to other
libraries of the *COBRA* ecosystem including [**escher**](https://escher.github.io),
[**cobrapy**](https://opencobra.github.io/cobrapy/), and [**optlang**](https://github.com/biosustain/optlang).

### Available methods

Name | Long name | Reference
--- | --- | ---
FBA | Flux Balance Analysis | (Varma and Palsson, 1993)
FVA | Flux Variability Analysis | (Mahadevan and Schilling, 2003)
pFBA | Parsimonious FBA | (Lewis et al, 2010)
FBrAtio | FBA with flux ratios| (Yen et al, 2013)
CAFBA | Constrained Allocation FBA | (Mori et al, 2016)
MOMA | Minimization of Metabolic Adjustment | (Segre et al, 2002)
lMOMA | linear MOMA | (Segre et al, 2002)
ROOM | Regulatory On/Off Minimization | (Shlomi et al, 2005)
ll-FBA | loopless FBA | (Schellenberger et al, 2011)
TFA | Thermodynamic Flux Analysis | (Henry et al, 2007)
TVA | Thermodynamic Variability Analysis | (Henry et al, 2007)
NET | Network-embedded thermodynamic analysis |  (Kummel et al 2006)
GIMME | Gene Inactivity Moderated by Met and Exp | (Becker and Palsson, 2008)
E-Flux | E-Flux | (Colijn et al, 2009)
SteadyCom | Community simulation | (Chan et al, 2017)

### Credits and License

Developed by Daniel Machado at the European Molecular Biology Laboratory (2019).

**ReFramed** is a refactored version of the [**framed**](https://github.com/cdanielmachado/framed) library.

Released under an Apache License.

