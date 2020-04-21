from ..core.environment import Environment
from ..cobra.simulation import FBA
from ..solvers import solver_instance
from ..solvers.solution import Status
from warnings import warn


def auxotrophies(model, constraints=None, min_rel_growth=1e-2, min_abs_growth=1e-5):

    solver = solver_instance(model)

    if not constraints:
        constraints = Environment.complete(model, inplace=False)

    sol = FBA(model, constraints=constraints, solver=solver, get_values=False)

    if sol.status != Status.OPTIMAL or sol.fobj < min_abs_growth:
        warn("Organism does not grow in given medium.")
        return
    else:
        ref_growth = sol.fobj

    auxo = []

    for r_id in model.get_exchange_reactions():
        tmp_lb, tmp_ub = constraints[r_id]
        constraints[r_id] = (0, tmp_ub)

        sol = FBA(model, constraints=constraints, solver=solver, get_values=False)
        if sol.fobj < min_abs_growth or sol.fobj < min_rel_growth * ref_growth:
            auxo.append(r_id)

        constraints[r_id] = (tmp_lb, tmp_ub)

    return auxo