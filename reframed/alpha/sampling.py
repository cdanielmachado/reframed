from ..solvers import solver_instance
from ..solvers.solution import Status
from ..core.environment import Environment
from ..cobra.simulation import FBA

from random import lognormvariate
from warnings import warn
from math import inf


def rpFBA(model, min_growth=0.1, reactions=None, constraints=None, solver=None):

    if not solver:
        solver = solver_instance(model)

    if not reactions:
        reactions = model.reactions.keys()

    if not constraints:
        constraints = {model.biomass_reaction: min_growth}
    else:
        constraints = constraints.copy()
        constraints[model.biomass_reaction] = min_growth

    if not hasattr(solver, 'rpFBA_flag'):
        solver.rpFBA_flag = True

        for r_id in reactions:
            if model.reactions[r_id].reversible:
                pos, neg = r_id + '+', r_id + '-'
                solver.add_variable(pos, 0, inf)
                solver.add_variable(neg, 0, inf)
        solver.update()

        for r_id in reactions:
            if model.reactions[r_id].reversible:
                pos, neg = r_id + '+', r_id + '-'
                solver.add_constraint('c' + pos, {r_id: -1, pos: 1}, '>', 0)
                solver.add_constraint('c' + neg, {r_id: 1, neg: 1}, '>', 0)
        solver.update()

    objective = dict()

    for r_id in reactions:
        weight = lognormvariate(0, 1)
        if model.reactions[r_id].reversible:
            pos, neg = r_id + '+', r_id + '-'
            objective[pos] = weight
            objective[neg] = weight
        else:
            objective[r_id] = weight

    sol = solver.solve(objective, minimize=True, constraints=constraints)

    return sol


def sampling(model, n=1, min_growth=0.1, reactions=None, constraints=None, as_df=False):

    solver = solver_instance(model)
    results = []

    for i in range(n):
        sol = rpFBA(model, min_growth=min_growth, reactions=reactions, constraints=constraints, solver=solver)

        if sol.status == Status.OPTIMAL:
            if as_df:
                results.append(sol.to_dataframe())
            else:
                results.append(sol.values)

    if as_df and len(results) > 0:
        import pandas as pd
        results = pd.concat(results, axis=1).T.reset_index(drop=True)

    return results


def random_medium(model, min_growth=0.1, exchange_reactions=None, validate=False, constraints=None,
                  solver=None, abstol=1e-9):

    if not exchange_reactions:
        exchange_reactions = model.get_exchange_reactions()

    sol = rpFBA(model, min_growth=min_growth, reactions=exchange_reactions,
                   constraints=constraints, solver=solver)

    medium = {r_id for r_id in model.get_exchange_reactions() if sol.values[r_id] < -abstol}

    if validate:
        env = Environment.from_reactions(medium).apply(model, inplace=False, exclusive=True)
        if constraints:
            env.update(constraints)
        sol = FBA(model, constraints=env)
        if sol.status != Status.OPTIMAL:
            warn("Failed to validate medium.")

    return medium
