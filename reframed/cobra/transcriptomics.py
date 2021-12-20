from ..solvers import solver_instance
from ..solvers.solution import Status
from ..core.transformation import gpr_transform
from .simulation import FBA, pFBA
from math import inf
from numpy import percentile
from collections.abc import Iterable


def marge(model, expr_a=None, expr_b=None, rel_expr=None, constraints_a=None, constraints_b=None,
         growth_frac_a=1.0, growth_frac_b=1.0, step2_relax=0.1, get_ranges=False, gene_prefix='G_', pseudo_genes=None):
    """ Metabolic Analysis with Relative Gene Expression (MARGE)

    Args:
        model (CBModel): organism model
        expr_a (dict): gene expression for condition A (optional)
        expr_b (dict): gene expression for condition B (optional)
        rel_expr (dict): relative gene expression (B / A) (optional)
        constraints_a (dict): additional constraints to use for condition A (optional)
        constraints_b (dict): additional constraints to use for condition B (optional)
        growth_frac_a (float): minimum growth rate in condition A (default: 1.0)
        growth_frac_b (float): minimum growth rate in condition B (default: 1.0)
        step2_relax (float): relaxation from main objective during second step (default: 0.1)
        get_ranges (bool): return flux ranges instead of single flux distributions (default: False)
        gene_prefix (str): prefix used in gene identifiers (default: 'G_')
        pseudo_genes (list): pseudo-genes in model to ignore (e.g: 'spontaneous') (optional)

    Returns:
        dict: fluxes (or flux ranges) in conditions A and B
    """

    if rel_expr is None:
        if expr_a is None or expr_b is None:
            raise RuntimeError("Please provide either relative or absolute expression values.")

        if min(expr_a.values()) < 0 or min(expr_b.values()) < 0:
            raise RuntimeError("Expression values must be non-negative.")

        common_genes = set(expr_a) & set(expr_b) & set(model.genes)
        if len(common_genes) == 0:
            raise RuntimeError("Gene identifiers do not match any genes in the model.")

        expr_a = {x: expr_a[x] for x in common_genes}
        expr_b = {x: expr_b[x] for x in common_genes}

    else:
        if expr_a is not None or expr_b is not None:
            raise RuntimeError("Please provide either relative or absolute expression values.")

        if min(rel_expr.values()) <= 0:
            raise RuntimeError("Relative expression values must be positive.")

        common_genes = set(rel_expr) & set(model.genes)
        if len(common_genes) == 0:
            raise RuntimeError("Gene identifiers do not match any genes in the model.")

        # converting relative to absolute gene expression such that geometric mean(a, b) = 1
        expr_a = {x: rel_expr[x] ** -0.5 for x in common_genes}
        expr_b = {x: rel_expr[x] ** 0.5 for x in common_genes}

    model0 = model
    model = gpr_transform(model0, inplace=False, add_proteome=True,
                              gene_prefix=gene_prefix, pseudo_genes=pseudo_genes)

    if constraints_a is None:
        constraints_a = {}
    else:
        constraints_a = model.convert_constraints(constraints_a)

    if constraints_b is None:
        constraints_b = {}
    else:
        constraints_b = model.convert_constraints(constraints_b)

    pre_solver = solver_instance(model)

    sol_a = pFBA(model, solver=pre_solver, reactions=["proteome_synth"], constraints=constraints_a)
    if sol_a.status != Status.OPTIMAL:
        raise RuntimeError('Failed to solve reference model for condition A.')

    sol_b = pFBA(model, solver=pre_solver, reactions=["proteome_synth"], constraints=constraints_b)
    if sol_b.status != Status.OPTIMAL:
        raise RuntimeError('Failed to solve reference model for condition B.')

    max_growth_a = sol_a.values[model.biomass_reaction]
    constraints_a[model.biomass_reaction] = (max_growth_a * growth_frac_a, inf)

    max_growth_b = sol_b.values[model.biomass_reaction]
    constraints_b[model.biomass_reaction] = (max_growth_b * growth_frac_b, inf)

    proteome_a = sol_a.values["proteome_synth"]
    proteome_b = sol_b.values["proteome_synth"]

    solver = solver_instance()

    for r_id, reaction in model.reactions.items():
        lb_a, ub_a = constraints_a.get(r_id, (reaction.lb, reaction.ub))
        solver.add_variable(r_id + '_a', lb_a, ub_a, update=False)

        lb_b, ub_b = constraints_b.get(r_id, (reaction.lb, reaction.ub))
        solver.add_variable(r_id + '_b', lb_b, ub_b, update=False)

    for g_id in common_genes:
        solver.add_variable(g_id + '_+', 0, inf, update=False)
        solver.add_variable(g_id + '_-', 0, inf, update=False)

    solver.update()

    table = model.metabolite_reaction_lookup()

    for m_id in model.metabolites:
        stoich_a = {r_id + '_a': val for r_id, val in table[m_id].items()}
        solver.add_constraint(m_id + '_a', stoich_a, update=False)

        stoich_b = {r_id + '_b': val for r_id, val in table[m_id].items()}
        solver.add_constraint(m_id + '_b', stoich_b, update=False)

    objective = {}

    for g_id in common_genes:
        g_id_p, g_id_m = g_id + '_+', g_id + '_-'
        g_c_p, g_c_m = g_id + '_c+', g_id + '_c-'

        objective[g_id_p] = 1
        objective[g_id_m] = 1

        u_id_a = 'u_' + g_id[len(gene_prefix):] + '_a'
        u_id_b = 'u_' + g_id[len(gene_prefix):] + '_b'

        solver.add_constraint(g_c_p, {u_id_b: expr_a[g_id], u_id_a: -expr_b[g_id], g_id_m: 1}, '>', 0, update=False)
        solver.add_constraint(g_c_m, {u_id_b: expr_a[g_id], u_id_a: -expr_b[g_id], g_id_p: -1}, '<', 0, update=False)

    solver.add_constraint("proteome_a", {"proteome_synth_a": 1}, '<', proteome_a, update=False)
    solver.add_constraint("proteome_b", {"proteome_synth_b": 1}, '<', proteome_b, update=False)

    solver.update()

    sol = solver.solve(objective, minimize=True)

    if sol.status != Status.OPTIMAL:
        raise RuntimeError("Failed to solve first problem.")

    obj1_max = sol.fobj * (1 + step2_relax)
    solver.add_constraint("obj1", objective, '<', obj1_max, update=True)

    if not get_ranges:

        objective2 = {"proteome_synth_a": 1, "proteome_synth_b": 1}
        sol2 = solver.solve(objective2, minimize=True)

        if sol.status != Status.OPTIMAL:
            raise RuntimeError("Failed to solve second problem.")

        fluxes_a = {r_id: sol2.values[r_id + "_a"] for r_id in model.reactions}
        fluxes_a = model.convert_fluxes(fluxes_a)

        fluxes_b = {r_id: sol2.values[r_id + "_b"] for r_id in model.reactions}
        fluxes_b = model.convert_fluxes(fluxes_b)

        return fluxes_a, fluxes_b

    else:
        if isinstance(get_ranges, Iterable):
            reactions = list(get_ranges)
        else:
            reactions = model0.reactions

        ranges_a = {r_id: [-inf, inf] for r_id in reactions}
        ranges_b = {r_id: [-inf, inf] for r_id in reactions}

        for r_id in reactions:
            obj_a = {x + '_a': val for x, val in model.convert_id_to_expr(r_id).items()}
            tmp = solver.solve(obj_a, minimize=True, get_values=False)
            if tmp.status == Status.OPTIMAL:
                ranges_a[r_id][0] = tmp.fobj

            obj_b = {x + '_b': val for x, val in model.convert_id_to_expr(r_id).items()}
            tmp = solver.solve(obj_b, minimize=True, get_values=False)
            if tmp.status == Status.OPTIMAL:
                ranges_b[r_id][0] = tmp.fobj

        for r_id in reactions:
            obj_a = {x + '_a': val for x, val in model.convert_id_to_expr(r_id).items()}
            tmp = solver.solve(obj_a, minimize=False, get_values=False)
            if tmp.status == Status.OPTIMAL:
                ranges_a[r_id][1] = tmp.fobj

            obj_b = {x + '_b': val for x, val in model.convert_id_to_expr(r_id).items()}
            tmp = solver.solve(obj_b, minimize=False, get_values=False)
            if tmp.status == Status.OPTIMAL:
                ranges_b[r_id][1] = tmp.fobj

        return ranges_a, ranges_b


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

