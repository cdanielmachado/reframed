from ..solvers import solver_instance
from ..solvers.solution import Status
from .simulation import FBA
from .thermodynamics import looplessFBA
from warnings import warn
from math import inf
from numpy import linspace


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
            solution = looplessFBA(model, {r_id: 1}, True, constraints=_constraints, internal=internal,
                                   solver=solver, get_values=False)
        else:
            solution = FBA(model, {r_id: 1}, True, constraints=_constraints, solver=solver, get_values=False)

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
            solution = looplessFBA(model, {r_id: 1}, False, constraints=_constraints, internal=internal,
                                   solver=solver, get_values=False)
        else:
            solution = FBA(model, {r_id: 1}, False, constraints=_constraints, solver=solver, get_values=False)

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


def flux_envelope(model, r_x, r_y, steps=10, constraints=None):
    """ Calculate the flux envelope for a pair of reactions.

    Arguments:
        model (CBModel): the model
        r_x (str): reaction on x-axis
        r_y (str): reaction on y-axis
        steps (int): number of steps to compute (default: 10)
        constraints (dict): custom constraints to the FBA problem

    Returns:
        tuple: x values, y_min values, y_max values
    """

    x_range = FVA(model, reactions=[r_x], constraints=constraints)
    xmin, xmax = x_range[r_x]
    xvals = linspace(xmin, xmax, steps).tolist()
    ymins, ymaxs = [None] * steps, [None] * steps

    if constraints is None:
        _constraints = {}
    else:
        _constraints = {}
        _constraints.update(constraints)

    for i, xval in enumerate(xvals):
        _constraints[r_x] = xval
        y_range = FVA(model, reactions=[r_y], constraints=_constraints)
        ymins[i], ymaxs[i] = y_range[r_y]

    return xvals, ymins, ymaxs

