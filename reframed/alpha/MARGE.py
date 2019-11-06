from ..solvers import solver_instance
from ..solvers.solution import Status
from ..cobra.simulation import FBA
from .GPRtransform import gpr_transform
from math import inf


def marge(model, rel_expression, transformed=False, constraints_a=None, constraints_b=None, rel_constraints=None,
          growth_frac_a=0.1, growth_frac_b=0.1, activation_frac=0.1, activation_max=1.0, step2_tol=0.1,
          gene_prefix='G_', pseudo_genes=None):
    """ Metabolic Analysis with Relative Gene Expression (MARGE)

    First step minimizes the objective (fluxes B - rel expr (B/A) * fluxes A)

    Second step adds a parsimony criteria (with a slight relaxation of the first objective)

    Args:
        model (CBModel): organism model
        rel_expression (dict): relative gene expression (condition B / condition A)
        transformed (bool): True if the model is already in extended GPR format (default: False)
        constraints_a (dict): additional constrants to use for condition A (optional)
        constraints_b (dict): additional constrants to use for condition B (optional)
        rel_constraints (dict): relative constraints between conditions (such as flux ratios) (default: False)
        growth_frac_a (float): minimum growth rate in condition A (default: 0.1)
        growth_frac_b (float): minimum growth rate in condition B (default: 0.1)
        activation_frac (float): flux activation threshold for expressed genes (fraction of max flux, default: 0.1)
        activation_max (float): max value for flux activation threshold (default: 1.0)
        step2_tol (float): relaxation from main objective during second step (default: 0.1)
        gene_prefix (str): prefix used in gene identifiers (default: 'G_')

        pseudo_genes (list): pseudo-genes in model to ignore (e.g: 'spontaneous') (optional)

    Returns:
        dict: fluxes in condition A
        dict: fluxes in condition B
        Solution: solver solution for step 1 (minimize flux changes)
        Solution: solver solution for step 2 (minimize absolute fluxes)

    """

    if not transformed:
        model = gpr_transform(model, inplace=False, gene_prefix=gene_prefix, pseudo_genes=pseudo_genes)

    if constraints_a is None:
        constraints_a = {}
    else:
        constraints_a = model.convert_constraints(constraints_a)

    if constraints_b is None:
        constraints_b = {}
    else:
        constraints_b = model.convert_constraints(constraints_b)

    if rel_constraints is None:
        rel_constraints = {}

    pre_solver = solver_instance(model)

    if growth_frac_a > 0:
        biomass = model.biomass_reaction
        sol_a = FBA(model, solver=pre_solver, constraints=constraints_a)
        if sol_a.status != Status.OPTIMAL:
            print('Failed to solve reference model for condition A.')
            return None, None, None, None
        constraints_a[biomass] = (sol_a.fobj * growth_frac_a, inf)

    if growth_frac_b > 0:
        biomass = model.biomass_reaction
        sol_b = FBA(model, solver=pre_solver, constraints=constraints_b)
        if sol_b.status != Status.OPTIMAL:
            print('Failed to solve reference model for condition B.')
            return None, None, None, None
        constraints_b[biomass] = (sol_b.fobj * growth_frac_b, inf)

    if activation_frac > 0:
        act_constraints_a = {}
        act_constraints_b = {}

        for g_id in rel_expression:
            r_id = 'u_' + g_id[len(gene_prefix):]

            if r_id not in constraints_a:
                max_flux = FBA(model, objective=r_id, solver=pre_solver, constraints=constraints_a)
                if max_flux.status == Status.OPTIMAL:
                    min_flux = min(max_flux.fobj * activation_frac, activation_max)
                    act_constraints_a[r_id] = (min_flux, model.reactions[r_id].ub)

            if r_id not in constraints_b:
                max_flux = FBA(model, objective=r_id, solver=pre_solver, constraints=constraints_b)
                if max_flux.status == Status.OPTIMAL:
                    min_flux = min(max_flux.fobj * activation_frac, activation_max)
                    act_constraints_b[r_id] = (min_flux, model.reactions[r_id].ub)

        constraints_a.update(act_constraints_a)
        constraints_b.update(act_constraints_b)

    solver = solver_instance()

    for r_id, reaction in model.reactions.items():
        lb_a, ub_a = constraints_a.get(r_id, (reaction.lb, reaction.ub))
        solver.add_variable(r_id + '_a', lb_a, ub_a, update=False)

        lb_b, ub_b = constraints_b.get(r_id, (reaction.lb, reaction.ub))
        solver.add_variable(r_id + '_b', lb_b, ub_b, update=False)

    for g_id, val in rel_expression.items():
        solver.add_variable(g_id + '_+', 0, inf, update=False)
        solver.add_variable(g_id + '_-', 0, inf, update=False)

    solver.update()

    table = model.metabolite_reaction_lookup()

    for m_id in model.metabolites:
        stoich_a = {r_id + '_a': val for r_id, val in table[m_id].items()}
        solver.add_constraint(m_id + '_a', stoich_a, update=False)

        stoich_b = {r_id + '_b': val for r_id, val in table[m_id].items()}
        solver.add_constraint(m_id + '_b', stoich_b, update=False)

    for r_id, ratio in rel_constraints.items():
        constr = {}
        expr_a = model.convert_id_to_expr(r_id, -ratio)
        expr_b = model.convert_id_to_expr(r_id, 1)
        constr.update({r_id2 + '_a': val for r_id2, val in expr_a.items()})
        constr.update({r_id2 + '_b': val for r_id2, val in expr_b.items()})
        solver.add_constraint(r_id + '_rel', constr, update=False)

    for g_id, val in rel_expression.items():
        u_id_a = 'u_' + g_id[len(gene_prefix):] + '_a'
        u_id_b = 'u_' + g_id[len(gene_prefix):] + '_b'
        solver.add_constraint(g_id + '_c+', {g_id + '_+': 1, u_id_b: -1, u_id_a: val}, '>', 0, update=False)
        solver.add_constraint(g_id + '_c-', {g_id + '_-': 1, u_id_b: 1, u_id_a: -val}, '>', 0, update=False)

    solver.update()

    objective1 = {}
    for g_id in rel_expression.keys():
        objective1[g_id + '_+'] = 1
        objective1[g_id + '_-'] = 1

    solution1 = solver.solve(objective1, minimize=True)

    if solution1.status != Status.OPTIMAL:
        print('Failed to solve first problem.')
        return None, None, solution1, None

    obj1_max = solution1.fobj * (1 + step2_tol)
    solver.add_constraint("obj1", objective1, '<', obj1_max, update=True)

    objective2 = {}
    for r_id in model.u_reactions:
        objective2[r_id + '_a'] = 1
        objective2[r_id + '_b'] = 1

    solution2 = solver.solve(objective2, minimize=True)

    if solution2.status != Status.OPTIMAL:
        print('Failed to solve second problem.')
        return None, None, solution1, solution2

    fluxes_a = {r_id: solution2.values[r_id + '_a'] for r_id in model.reactions}
    fluxes_b = {r_id: solution2.values[r_id + '_b'] for r_id in model.reactions}

    fluxes_a = model.convert_fluxes(fluxes_a)
    fluxes_b = model.convert_fluxes(fluxes_b)

    return fluxes_a, fluxes_b, solution1, solution2
