from ..solvers import solver_instance
from ..solvers.solution import Status
from .simulation import FBA
from .thermodynamics import llFBA
from warnings import warn
from math import inf


def FVA(model, obj_frac=0, reactions=None, constraints=None, loopless=False, internal=None, solver=None):
    """ Run Flux Variability Analysis (FVA).

    Arguments:
        model (CBModel): a constraint-based model
        obj_frac (float): minimum fraction of the maximum growth rate (default 0.0, max: 1.0)
        reactions (list): list of reactions to analyze (default: all)
        constraints (dict): additional constraints (optional)
        loopless (bool): run looplessFBA internally (very slow) (default: false)
        internal (list): list of internal reactions for looplessFBA (optional)
        solver (Solver): pre-instantiated solver instance (optional)

    Returns:
        dict: flux variation ranges
    """

    _constraints = {}
    if constraints:
        _constraints.update(constraints)

    if not solver:
        solver = solver_instance(model)

    if obj_frac > 0:
        target = model.biomass_reaction
        solution = FBA(model, objective=target, constraints=constraints, solver=solver)
        _constraints[target] = (obj_frac * solution.fobj, inf)

    if not reactions:
        reactions = model.reactions.keys()

    variability = {r_id: [None, None] for r_id in reactions}

    for r_id in reactions:
        if loopless:
            solution = llFBA(model, r_id, True, constraints=_constraints, internal=internal,
                             solver=solver, get_values=False)
        else:
            solution = FBA(model, r_id, True, constraints=_constraints, solver=solver, get_values=False)

        if solution.status == Status.OPTIMAL:
            variability[r_id][0] = solution.fobj
        elif solution.status == Status.UNBOUNDED:
            variability[r_id][0] = -inf
        elif solution.status == Status.INF_OR_UNB:
            variability[r_id][0] = -inf
        elif solution.status == Status.INFEASIBLE:
            warn('Infeasible solution status')
        else:
            warn('Unknown solution status')

    for r_id in reactions:
        if loopless:
            solution = llFBA(model, r_id, False, constraints=_constraints, internal=internal,
                             solver=solver, get_values=False)
        else:
            solution = FBA(model, r_id, False, constraints=_constraints, solver=solver, get_values=False)

        if solution.status == Status.OPTIMAL:
            variability[r_id][1] = solution.fobj
        elif solution.status == Status.UNBOUNDED:
            variability[r_id][1] = inf
        elif solution.status == Status.INF_OR_UNB:
            variability[r_id][1] = inf
        elif solution.status == Status.INFEASIBLE:
            warn('Infeasible solution status')
        else:
            warn('Unknown solution status')

    return variability


def blocked_reactions(model, constraints=None, reactions=None, abstol=1e-9):
    """ Find all blocked reactions in a model

    Arguments:
        model (CBModel): a constraint-based model
        constraints (dict): additional constraints (optional)
        reactions (list): List of reactions which will be tested (default: None, test all reactions)
        abstol (float): absolute tolerance (default: 1e-9)

    Returns:
        list: blocked reactions
    """

    variability = FVA(model, reactions=reactions, constraints=constraints)

    return [r_id for r_id, (lb, ub) in variability.items() if (abs(lb) + abs(ub)) < abstol]

