==========
Interfaces
==========


CobraPy
_______

ReFramed supports conversion of its models from and to CobraPy_ models.
This allows users to take full advantage of all the features present in both packages.

.. _CobraPy: https://opencobra.github.io/cobrapy/

You can load a model in ReFramed and simulate it with CobraPy:

::

    > from reframed import load_cbmodel, to_cobrapy

    > rf_model = load_cbmodel("e_coli_core.xml")
    > cb_model = to_cobrapy(rf_model)
    > cb_model.optimize()

    Optimal solution with objective value 0.874

You can load a model with CobraPy and simulate it with ReFramed:

::

    > from cobra.io import read_sbml_model
    > from reframed import from_cobrapy, FBA

    > cb_model = read_sbml_model("e_coli_core.xml")
    > rf_model = from_cobrapy(cb_model)
    > FBA(rf_model)

    Objective: 0.8739215069684306
    Status: Optimal


Escher
______

If you are using Jupyter notebooks, you can easily visualize flux distributions using Escher_.

.. _Escher: https://escher.readthedocs.io

The following will start the Escher widget inside a Jupyter notebook cell:

::

    > from reframed import fluxes2escher
    > sol = pFBA(model)
    > fluxes2escher(sol.values)

You can also specify which Escher map to use and send additional arguments directly to the Escher API:

::

    > fluxes2escher(sol.values, map_name='e_coli_core.Core metabolism')

    > fluxes2escher(sol.values, hide_secondary_metabolites=True)

