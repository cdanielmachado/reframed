==========
Simulation
==========

Basics
______

You can run any of the available simulation methods by importing it and applying it to the model.

::

    > from reframed import FBA
    > FBA(model)

    Objective: 0.8739215069684306
    Status: Optimal

Ideally you will want to store the solution object for further inspection.

::

    > sol = FBA(model)
    > sol.show_values()

    R_ACONTa      6.00725
    R_ACONTb      6.00725
    R_AKGDH       5.06438
    R_ATPM        8.39
    R_ATPS4r      45.514
    (...)

You can perform some basic queries on the solution values using regular expressions.

::

    > sol.show_values(pattern="R_EX", sort=True)

    R_EX_o2_e    -21.7995
    R_EX_glc__D_e -10
    R_EX_nh4_e   -4.76532
    R_EX_pi_e    -3.2149
    R_EX_h_e      17.5309
    R_EX_co2_e    22.8098
    R_EX_h2o_e    29.1758

All methods support additional constraints. These are temporary constraints used for simulation that
will not be permanently stored in the model.

::

    > from reframed import MOMA
    > sol = MOMA(model, constraints={"R_PGI": 0, "R_EX_o2_e": (-15, 0)})

For simulating multiple knockouts you can use higher level methods:

::

    > from reframed import gene_knockout
    > sol = gene_knockout(model, ["G_b1723", "G_b3916"], constraints={"R_EX_o2_e: (-15, 0)"})

    > from reframed import reaction_knockout
    > sol = reaction_knockout(model, ["R_PGI", "R_PFK", "R_MDH"])

There are multiple methods for other tasks:

::

    > from reframed import essential_genes
    > essential = essential_genes(model)

    > from reframed import blocked_reactions
    > blocked = blocked_reactions(model)


Using Environments
__________________


The **Environment** class supports many helper functions to create environmental conditions for a model.

::

    > from reframed import Environment
    > Environment.empty(model)

    R_EX_ac_e   	0	inf
    R_EX_acald_e	0	inf
    R_EX_akg_e  	0	inf
    (...)

::

    > Environment.complete(model, max_uptake=10)

    R_EX_ac_e   	-10	inf
    R_EX_acald_e	-10	inf
    R_EX_akg_e  	-10	inf
    (...)

::

    > Environment.from_compounds(["glc", "o2", "nh4"])

    R_EX_glc_e	-10.0	inf
    R_EX_o2_e	-10.0	inf
    R_EX_nh4_e	-10.0	inf

You can create an environment and use it as temporary simulation constraints:

::

    > env = Environment.complete(model)
    > sol = FBA(model, constraints=env)

Or you can apply the environment as permanent constraints to the model:

::

    > Environment.complete(model, inplace=True)
    > sol = FBA(model)

