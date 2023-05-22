.. figure:: reframed_logo.png
   :alt: ReFramed

   ReFramed

ReFramed: metabolic modeling package
====================================

**ReFramed** implements many constraint-based simulation methods (see
list below), and contains interfaces to other libraries of the *COBRA*
ecosystem including `escher <https://escher.github.io>`__,
`cobrapy <https://opencobra.github.io/cobrapy/>`__, and
`optlang <https://github.com/biosustain/optlang>`__.

Available methods
-----------------

+-----------------------+-----------------------+-----------------------+
| Name                  | Long name             | Reference             |
+=======================+=======================+=======================+
| FBA                   | Flux Balance Analysis | (Varma and Palsson,   |
|                       |                       | 1993)                 |
+-----------------------+-----------------------+-----------------------+
| FVA                   | Flux Variability      | (Mahadevan and        |
|                       | Analysis              | Schilling, 2003)      |
+-----------------------+-----------------------+-----------------------+
| pFBA                  | Parsimonious FBA      | (Lewis et al, 2010)   |
+-----------------------+-----------------------+-----------------------+
| FBrAtio               | FBA with flux ratios  | (Yen et al, 2013)     |
+-----------------------+-----------------------+-----------------------+
| CAFBA                 | Constrained           | (Mori et al, 2016)    |
|                       | Allocation FBA        |                       |
+-----------------------+-----------------------+-----------------------+
| MOMA                  | Minimization of       | (Segre et al, 2002)   |
|                       | Metabolic Adjustment  |                       |
+-----------------------+-----------------------+-----------------------+
| lMOMA                 | linear MOMA           | (Segre et al, 2002)   |
+-----------------------+-----------------------+-----------------------+
| ROOM                  | Regulatory On/Off     | (Shlomi et al, 2005)  |
|                       | Minimization          |                       |
+-----------------------+-----------------------+-----------------------+
| ll-FBA                | loopless FBA          | (Schellenberger et    |
|                       |                       | al, 2011)             |
+-----------------------+-----------------------+-----------------------+
| TFA                   | Thermodynamic Flux    | (Henry et al, 2007)   |
|                       | Analysis              |                       |
+-----------------------+-----------------------+-----------------------+
| TVA                   | Thermodynamic         | (Henry et al, 2007)   |
|                       | Variability Analysis  |                       |
+-----------------------+-----------------------+-----------------------+
| NET                   | Network-embedded      | (Kummel et al 2006)   |
|                       | thermodynamic         |                       |
|                       | analysis              |                       |
+-----------------------+-----------------------+-----------------------+
| GIMME                 | Gene Inactivity       | (Becker and Palsson,  |
|                       | Moderated by Met and  | 2008)                 |
|                       | Exp                   |                       |
+-----------------------+-----------------------+-----------------------+
| E-Flux                | E-Flux                | (Colijn et al, 2009)  |
+-----------------------+-----------------------+-----------------------+
| SteadyCom             | Community simulation  | (Chan et al, 2017)    |
+-----------------------+-----------------------+-----------------------+

Credits and License
-------------------

Developed by Daniel Machado at the Norwegian University of Science and Technology (NTNU).

**ReFramed** is a refactored version of the
`framed <https://github.com/cdanielmachado/framed>`__ library.

Released under an Apache License.
