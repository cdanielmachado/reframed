from ..solvers import solver_instance
from .simulation import FBA, pFBA
from numpy import percentile
from math import inf


def gene2rxn(gpr, gene_exp, and_func=min, or_func=sum):

    def f_and(x):
        x2 = [xi for xi in x if xi is not None]
        return and_func(x2) if x2 else None

    def f_or(x):
        x2 = [xi for xi in x if xi is not None]
        return or_func(x2) if x2 else None

    level = f_or([f_and([gene_exp[gene]
                         for gene in protein.genes if gene in gene_exp])
                  for protein in gpr.proteins])

    return level


def gene_to_reaction_expression(model, gene_exp, and_func=min, or_func=sum):
    rxn_exp = {}
    for r_id, reaction in model.reactions.items():
        if reaction.gpr is not None:
            level = gene2rxn(reaction.gpr, gene_exp, and_func, or_func)
            if level is not None:
                rxn_exp[r_id] = level
    return rxn_exp


def GIMME(model, gene_exp, cutoff=25, growth_frac=0.9, constraints=None, parsimonious=False):
    """ Run a GIMME simulation (Becker and Palsson, 2008).

    Arguments:
        model (CBModel): model
        gene_exp (dict): transcriptomics data
        cutoff (int): percentile cuttof (default: 25)
        growth_frac (float): minimum growth requirement (default: 0.9)
        constraints (dict): additional constraints
        parsimonious (bool): compute a parsimonious solution (default: False)

    Returns:
        Solution: solution
    """

    rxn_exp = gene_to_reaction_expression(model, gene_exp, or_func=max)

    threshold = percentile(list(rxn_exp.values()), cutoff)
    coeffs = {r_id: threshold-val for r_id, val in rxn_exp.items() if val < threshold}

    solver = solver_instance(model)

    wt_solution = FBA(model, constraints=constraints, solver=solver)

    if not constraints:
        constraints = {}

    biomass = model.biomass_reaction
    constraints[biomass] = (growth_frac * wt_solution.values[biomass], inf)

    for r_id in model.reactions:
        if model.reactions[r_id].reversible:
            pos, neg = r_id + '+', r_id + '-'
            solver.add_variable(pos, 0, inf, update=False)
            solver.add_variable(neg, 0, inf, update=False)
    solver.update()

    for r_id in model.reactions:
        if model.reactions[r_id].reversible:
            pos, neg = r_id + '+', r_id + '-'
            solver.add_constraint('c' + pos, {r_id: -1, pos: 1}, '>', 0, update=False)
            solver.add_constraint('c' + neg, {r_id: 1, neg: 1}, '>', 0, update=False)
    solver.update()

    objective = dict()
    for r_id, val in coeffs.items():
        if model.reactions[r_id].reversible:
            pos, neg = r_id + '+', r_id + '-'
            objective[pos] = val
            objective[neg] = val
        else:
            objective[r_id] = val

    solution = solver.solve(objective, minimize=True, constraints=constraints)

    if parsimonious:
        pre_solution = solution

        solver.add_constraint('obj', objective, '=', pre_solution.fobj)
        objective = dict()

        for r_id in model.reactions:
            if model.reactions[r_id].reversible:
                pos, neg = r_id + '+', r_id + '-'
                objective[pos] = 1
                objective[neg] = 1
            else:
                objective[r_id] = 1

        solution = solver.solve(objective, minimize=True, constraints=constraints)
        solver.remove_constraint('obj')
        solution.pre_solution = pre_solution

    for r_id in model.reactions:
        if model.reactions[r_id].reversible:
            pos, neg = r_id + '+', r_id + '-'
            del solution.values[pos]
            del solution.values[neg]

    return solution


def eFlux(model, gene_exp, scale_rxn=None, scale_value=1, constraints=None, parsimonious=False):
    """ Run an E-Flux simulation (Colijn et al, 2009).

    Arguments:
        model (CBModel): model
        gene_exp (dict): transcriptomics data
        scale_rxn (str): reaction to scale flux vector (optional)
        scale_value (float): scaling factor (mandatory if scale_rxn is specified)
        constraints (dict): additional constraints (optional)
        parsimonious (bool): compute a parsimonious solution (default: False)

    Returns:
        Solution: solution
    """

    rxn_exp = gene_to_reaction_expression(model, gene_exp)
    max_exp = max(rxn_exp.values())
    bounds = {}

    for r_id, reaction in model.reactions.items():
        val = rxn_exp[r_id] / max_exp if r_id in rxn_exp else 1
        lb2 = -val if reaction.lb < 0 else 0
        ub2 = val if reaction.ub > 0 else 0
        bounds[r_id] = (lb2, ub2)

    if constraints:
        for r_id, x in constraints.items():
            lb, ub = x if isinstance(x, tuple) else (x, x)
            lb2 = -1 if lb < 0 else 0
            ub2 = 1 if ub > 0 else 0
            bounds[r_id] = (lb2, ub2)

    if parsimonious:
        sol = pFBA(model, constraints=bounds)
    else:
        sol = FBA(model, constraints=bounds)

    if scale_rxn is not None:

        if sol.values[scale_rxn] != 0:
            k = abs(scale_value / sol.values[scale_rxn])
        else:
            k = 0

        for r_id, val in sol.values.items():
            sol.values[r_id] = val * k

    return sol

