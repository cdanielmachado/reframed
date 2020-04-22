from ..solvers import solver_instance
from ..solvers.solution import Status
from ..cobra.simulation import pFBA, lMOMA, FBA
from .GPRtransform import gpr_transform


def marge(model, rel_expression, transformed=False, constraints_a=None, constraints_b=None,
          growth_frac_a=1.0, growth_frac_b=0.0, activation=0.0, gene_prefix='G_', pseudo_genes=None):
    """ Metabolic Analysis with Relative Gene Expression 2.0 (MARGE2)

    This method integrates gene expression into flux balance analysis in two steps:

    - Compute enzyme usage fluxes for reference condition (condition A) using gene-pFBA [1].
    - Compute enzyme usage fluxes for perturbed condition (condition B) using gene-lMOMA [1]
       and the relative gene-expression data between the two conditions (B / A).

    [1] Machado et al, PLoS Computational Biology, 2016.

    Args:
        model (CBModel): organism model
        rel_expression (dict): relative gene expression (condition B / condition A)
        transformed (bool): True if the model is already in extended GPR format (default: False)
        constraints_a (dict): additional constraints to use for condition A (optional)
        constraints_b (dict): additional constraints to use for condition B (optional)
        growth_frac_a (float): minimum growth rate in condition A (default: 1.0)
        growth_frac_b (float): minimum growth rate in condition B (default: 0.0)
        gene_prefix (str): prefix used in gene identifiers (default: 'G_')
        pseudo_genes (list): pseudo-genes in model to ignore (e.g: 'spontaneous') (optional)

    Returns:
        dict: fluxes in condition A
        dict: fluxes in condition B
        Solution: solver solution for step 1 (pFBA)
        Solution: solver solution for step 2 (lMOMA)

    """

    if not transformed:
        model = gpr_transform(model, inplace=False, gene_prefix=gene_prefix, pseudo_genes=pseudo_genes)

    solver = solver_instance(model)

    if constraints_a is None:
        constraints_a = {}
    else:
        constraints_a = model.convert_constraints(constraints_a)

    if constraints_b is None:
        constraints_b = {}
    else:
        constraints_b = model.convert_constraints(constraints_b)

    sol1 = pFBA(model, obj_frac=growth_frac_a, constraints=constraints_a, solver=solver,
                reactions=model.u_reactions, cleanup=False)

    if sol1.status != Status.OPTIMAL:
        print('Failed to solve first problem.')
        return None, None, sol1, None

    u_ref = {}
    for g_id, val in rel_expression.items():
        if g_id not in model.genes:
            continue
        u_id = 'u_' + g_id[len(gene_prefix):]
        if val > 1:
            u_ref[u_id] = val * max(sol1.values[u_id], activation)
        else:
            u_ref[u_id] = val * sol1.values[u_id]

    if growth_frac_b > 0:
        tmp = FBA(model, constraints=constraints_b, solver=solver)
        constraints_b[model.biomass_reaction] = (growth_frac_b * tmp.fobj, tmp.fobj)

    sol2 = lMOMA(model, reference=u_ref, constraints=constraints_b, reactions=model.u_reactions, solver=solver)

    if sol2.status != Status.OPTIMAL:
        print('Failed to solve sol2 problem.')
        return None, None, sol1, sol2

    fluxes_a = {r_id: sol1.values[r_id] for r_id in model.reactions}
    fluxes_b = {r_id: sol2.values[r_id] for r_id in model.reactions}

    fluxes_a = model.convert_fluxes(fluxes_a)
    fluxes_b = model.convert_fluxes(fluxes_b)

    return fluxes_a, fluxes_b, sol1, sol2
