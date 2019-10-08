from ..solvers import solver_instance
from ..solvers.solver import VarType
from ..solvers.solution import Status
from ..core.model import ReactionType
from warnings import warn
from collections import Iterable
from math import inf


def FBA(model, objective=None, minimize=False, constraints=None, solver=None, get_values=True,
        shadow_prices=False, reduced_costs=False):
    """ Run a Flux Balance Analysis (FBA) simulation:

    Arguments:
        model (CBModel): a constraint-based model
        objective (dict: objective coefficients (optional)
        minimize (bool): minimize objective function (False by default)
        constraints (dict): environmental or additional constraints (optional)
        solver (Solver): solver instance instantiated with the model, for speed (optional)
        get_values (bool): set to false for speedup if you only care about the objective value (optional, default: True)
        shadow_prices (bool): retrieve shadow prices (default: False)
        reduced_costs (bool): retrieve reduced costs (default: False)

    Returns:
        Solution: solution
    """

    if not objective:
        objective = model.get_objective()

        if len(objective) == 0:
            warn('Model objective undefined.')

    if not solver:
        solver = solver_instance(model)

    solution = solver.solve(objective, minimize=minimize, constraints=constraints, get_values=get_values,
                            shadow_prices=shadow_prices, reduced_costs=reduced_costs)
    return solution


def pFBA(model, objective=None, obj_frac=None, minimize=False, constraints=None, reactions=None, solver=None,
         cleanup=True):
    """ Run a parsimonious Flux Balance Analysis (pFBA) simulation:

    Arguments:
        model (CBModel): a constraint-based model
        objective (dict): objective coefficients (optional)
        obj_frac (float): require only a fraction of the main objective during the flux minimization step (optional)
        minimize (bool): sense of optimization (maximize by default)
        constraints (dict: environmental or additional constraints (optional)
        reactions (list): list of reactions to be minimized (optional, default: all)
        solver (Solver): solver instance instantiated with the model, for speed (optional)
        cleanup (bool): remove temporary variables from solution (default: True)

    Returns:
        Solution: solution
    """

    if not solver:
        solver = solver_instance(model)

    if not objective:
        objective = model.get_objective()

    pre_solution = FBA(model, objective, minimize, constraints, solver)

    if pre_solution.status != Status.OPTIMAL:
        return pre_solution

    if obj_frac is None:
        solver.add_constraint('obj', objective, '=', pre_solution.fobj)
    else:
        solver.add_constraint('obj', objective, '>', obj_frac * pre_solution.fobj)

    if not reactions:
        reactions = model.reactions.keys()

    if not hasattr(solver, 'pFBA_flag'):
        solver.pFBA_flag = True
        for r_id in reactions:
            if model.reactions[r_id].reversible:
                pos, neg = r_id + '+', r_id + '-'
                solver.add_variable(pos, 0, inf, update=False)
                solver.add_variable(neg, 0, inf,  update=False)
        solver.update()
        for r_id in reactions:
            if model.reactions[r_id].reversible:
                pos, neg = r_id + '+', r_id + '-'
                solver.add_constraint('c' + pos, {r_id: -1, pos: 1}, '>', 0, update=False)
                solver.add_constraint('c' + neg, {r_id: 1, neg: 1}, '>', 0, update=False)
        solver.update()

    objective = dict()
    for r_id in reactions:
        if model.reactions[r_id].reversible:
            pos, neg = r_id + '+', r_id + '-'
            objective[pos] = 1
            objective[neg] = 1
        else:
            objective[r_id] = 1

    solution = solver.solve(objective, minimize=True, constraints=constraints)
    solver.remove_constraint('obj')
    solution.pre_solution = pre_solution

    if cleanup:
        for r_id in reactions:
            if model.reactions[r_id].reversible:
                pos, neg = r_id + '+', r_id + '-'
                del solution.values[pos]
                del solution.values[neg]

    return solution


def CAFBA(model, objective=None, minimize=False, wc=0, we=8.3e-4, wr=0.169, pmax=0.484, carbon_source="M_glc__D_e",
          constraints=None, solver=None, cleanup=True):
    """ Constrained Allocation Flux Balance Analysis (Mori et al, 2016).

    Arguments:
        model (CBModel): a constraint-based model
        objective (dict): objective coefficients (optional)
        minimize (bool): sense of optimization (maximize by default)
        wc (float): weight of carbon uptake sector (default: 0.0)
        we (float): weight of biosynthetic enzymes sector (default: 8.3e-4)
        wr (float): weight of ribosome-affiliated proteins (default: 0.169)
        pmax (float): maximum protein allocation excluding Q-sector (default: 0.484)
        constraints (dict: environmental or additional constraints (optional)
        solver (Solver): solver instance instantiated with the model, for speed (optional)
        cleanup (bool): remove temporary variables from solution (default: True)

    Notes:
        The default parameter values used here are the same as in the original publication.

    Returns:
        Solution: solution
    """
    if solver is None:
        solver = solver_instance(model)

    if objective is None:
        objective = model.get_objective()

    if hasattr(solver, 'CAFBA_flag'):
        uptake = solver.uptake
        enzymatic = solver.enzymatic
    else:
        solver.CAFBA_flag = True
        solver.uptake = []
        solver.enzymatic = []
        uptake = solver.uptake
        enzymatic = solver.enzymatic
        splits = {}

        for r_id in model.get_reactions_by_type(ReactionType.TRANSPORT):
            if carbon_source in model.reactions[r_id].stoichiometry:
                if model.reactions[r_id].reversible:
                    pos, neg = r_id + '+', r_id + '-'
                    solver.add_variable(pos, 0, inf, update=False)
                    solver.add_variable(neg, 0, inf, update=False)
                    uptake.append(r_id + '+')
                    uptake.append(r_id + '-')
                    splits[r_id] = pos, neg
                else:
                    uptake.append(r_id)

        for r_id in model.get_reactions_by_type(ReactionType.ENZYMATIC):
            if model.reactions[r_id].reversible:
                pos, neg = r_id + '+', r_id + '-'
                solver.add_variable(pos, 0, inf, update=False)
                solver.add_variable(neg, 0, inf, update=False)
                enzymatic.append(r_id + '+')
                enzymatic.append(r_id + '-')
                splits[r_id] = pos, neg
            else:
                enzymatic.append(r_id)

        solver.update()

        for r_id, (pos, neg) in splits.items():
            solver.add_constraint('c' + pos, {r_id: -1, pos: 1}, '>', 0, update=False)
            solver.add_constraint('c' + neg, {r_id: 1, neg: 1}, '>', 0, update=False)
        solver.update()

    main_constr = {}
    for r_id in enzymatic:
        main_constr[r_id] = we
    for r_id in uptake:
        main_constr[r_id] = wc
    main_constr[model.biomass_reaction] = wr

    solver.add_constraint('alloc', main_constr, '=', pmax, update=True)

    solution = solver.solve(objective, minimize=minimize, constraints=constraints)

    solver.remove_constraint('alloc')

    if cleanup:
        for (pos, neg) in splits.values():
            del solution.values[pos]
            del solution.values[neg]

    return solution


def MOMA(model, reference=None, constraints=None, reactions=None, solver=None):
    """ Run a Minimization Of Metabolic Adjustment (MOMA) simulation:

    Arguments:
        model (CBModel): a constraint-based model
        reference (dict): reference flux distribution (optional)
        constraints (dict): environmental or additional constraints (optional)
        reactions (list): list of reactions to include in the objective (optional, default: all)
        solver (Solver): solver instance instantiated with the model, for speed (optional)

    Returns:
        Solution: solution
    """

    if reference is None:
        wt_solution = pFBA(model)
        reference = wt_solution.values
    else:
        reactions = reference.keys()

    if reactions is None:
        reactions = model.reactions.keys()

    quad_obj = {(r_id, r_id): 1 for r_id in reactions}
    lin_obj = {r_id: -2 * reference[r_id] for r_id in reactions}

    if not solver:
        solver = solver_instance(model)

    solution = solver.solve(lin_obj, quadratic=quad_obj, minimize=True, constraints=constraints)

    solution.reference = reference

    return solution


def lMOMA(model, reference=None, constraints=None, reactions=None, solver=None):
    """ Run a (linear version of) Minimization Of Metabolic Adjustment (lMOMA) simulation:

    Arguments:
        model (CBModel): a constraint-based model
        reference (dict): reference flux distribution (optional)
        constraints (dict): environmental or additional constraints (optional)
        reactions (list): list of reactions to include in the objective (optional, default: all)
        solver (Solver): solver instance instantiated with the model, for speed (optional)

    Returns:
        Solution: solution
    """

    if reference is None:
        wt_solution = pFBA(model)
        reference = wt_solution.values
    else:
        reactions = reference.keys()

    if reactions is None:
        reactions = model.reactions.keys()

    if not solver:
        solver = solver_instance(model)

    if not hasattr(solver, 'lMOMA_flag'):
        solver.lMOMA_flag = True
        for r_id in reactions:
            d_pos, d_neg = r_id + '_d+', r_id + '_d-'
            solver.add_variable(d_pos, 0, inf, update=False)
            solver.add_variable(d_neg, 0, inf, update=False)
        solver.update()
        for r_id in reactions:
            d_pos, d_neg = r_id + '_d+', r_id + '_d-'
            solver.add_constraint('c' + d_pos, {r_id: -1, d_pos: 1}, '>', -reference[r_id], update=False)
            solver.add_constraint('c' + d_neg, {r_id: 1, d_neg: 1}, '>', reference[r_id], update=False)
        solver.update()

    objective = dict()
    for r_id in reactions:
        d_pos, d_neg = r_id + '_d+', r_id + '_d-'
        objective[d_pos] = 1
        objective[d_neg] = 1

    solution = solver.solve(objective, minimize=True, constraints=constraints)

    solution.reference = reference

    return solution


def ROOM(model, reference=None, constraints=None, wt_constraints=None, reactions=None, solver=None,
         delta=0.03, epsilon=0.001, pool_size=0):
    """ Run a Regulatory On/Off Minimization (ROOM) simulation:

    Arguments:
        model (CBModel): a constraint-based model
        reference (dict): reference flux distribution or flux ranges (optional)
        constraints (dict): environmental or additional constraints (optional)
        wt_constraints (dict): constraints to calculate wild-type phenotype (optional)
        reactions (list): list of reactions to include in the objective (optional, default: all)
        solver (Solver): solver instance instantiated with the model, for speed (optional)
        delta (float): relative tolerance (default: 0.03)
        epsilon (float): absolute tolerance (default: 0.001)
        pool_size (int): set to nonzero value to return a pool of solutions (default: 0)

    Returns:
        Solution: solution
    """

    U = 1e6
    L = -1e6

    if not solver:
        solver = solver_instance(model)

    if reference is None:
        wt_solution = pFBA(model, constraints=wt_constraints, solver=solver)
        reference = wt_solution.values

    if reactions is None:
        reactions = reference.keys()

    objective = dict()
    if not hasattr(solver, 'ROOM_flag'):
        solver.ROOM_flag = True

        for r_id in reactions:
            y_i = 'y_' + r_id
            solver.add_variable(y_i, 0, 1, vartype=VarType.BINARY, update=False)
            objective[y_i] = 1
        solver.update()

        for r_id in reactions:
            y_i = 'y_' + r_id
            if isinstance(reference[r_id], Iterable):
                w_i_min = reference[r_id][0] if reference[r_id][0] != -inf else -1000
                w_i_max = reference[r_id][1] if reference[r_id][1] != inf else 1000
            else:
                w_i_min = reference[r_id]
                w_i_max = reference[r_id]
            w_u = w_i_max + delta * abs(w_i_max) + epsilon
            w_l = w_i_min - delta * abs(w_i_min) - epsilon
            solver.add_constraint('c' + r_id + '_u', {r_id: 1, y_i: (w_u - U)}, '<', w_u, update=False)
            solver.add_constraint('c' + r_id + '_l', {r_id: 1, y_i: (w_l - L)}, '>', w_l, update=False)
        solver.update()

    solution = solver.solve(objective, minimize=True, constraints=constraints, pool_size=pool_size)

    if pool_size == 0:
        solution.reference = reference

    return solution
