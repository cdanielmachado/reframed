===========
Basic usage
===========


Loading and saving
__________________


ReFramed can load and save models SBML format.
You optionally can use the *flavor* argument to indicate if the model uses the new *fbc2* format (default) or the legacy *cobra* format.

::

    from reframed import load_cbmodel
    model = load_cbmodel('e_coli_core.xml', flavor='fbc2')

If you make any changes to the model, you can choose to save it in a different flavor.

::

    from reframed import save_cbmodel
    save_cbmodel(model, 'e_coli_core_modified.xml', flavor='cobra')


Model manipulation
__________________

Once you load a model you can easily access the model's attributes and perform different kinds of operations using the *model* object.

You can quickly inspect your model by printing it:

::

    > model

    R_ACALD: M_acald_c + M_coa_c + M_nad_c <-> M_accoa_c + M_h_c + M_nadh_c
    R_ACALDt: M_acald_e <-> M_acald_c
    R_ACKr: M_ac_c + M_atp_c <-> M_actp_c + M_adp_c
    R_ACONTa: M_cit_c <-> M_acon_C_c + M_h2o_c
    R_ACONTb: M_acon_C_c + M_h2o_c <-> M_icit_c
    (...)

The summary method also gives a quick overview of metabolites, compartments, and reactions:

::

    > model.summary()

    Metabolites:
    c 52
    e 20

    Reactions:
    enzymatic 48
    transport 25
    exchange 20
    sink 0
    other 2

You can also inspect particular metabolites and reactions.

::

    > model.reactions.R_PGI

    R_PGI: M_g6p_c <-> M_f6p_c

::

    > model.metabolites.M_g6p_c

    D-Glucose-6-phosphate

::

    > model.reactions.R_TKT1.gpr

    (G_b2465 or G_b2935)

You can easily make changes to a model. For instance, let's change the maximum glucose uptake rate:

::

    model.reactions.R_EX_glc_e.lb = -10

But for programmatic access, you may want to use reactions as a dictionary instead:

::

    for reaction in model.reactions.values():
        reaction.lb = 0

You can easily add a new reactions to a model (optionally specifying flux bounds) :

::

    model.add_reaction_from_str('R_abc: A + 0.5 B --> 2 C')
    model.add_reaction_from_str('R_s2p: S <-> P [-10, 10]')

You can add ratio constraints:

::

    model.add_ratio_constraint("R_PGI", "R_G6PDH", 1.5)

You can also remove reactions, metabolites, genes or compartments:

::

    model.remove_reaction('R_PGI')
    model.remove_metabolites(['M_h2o_c', 'M_h_c'])


There is a lot more you can try! Just take a look into the API.

