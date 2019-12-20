from .solution import CommunitySolution
from ..solvers.solution import Status
from ..core.model import ReactionType
from ..solvers import solver_instance
from warnings import warn
from math import inf, isinf


def SteadyCom(community, constraints=None, solver=None):
    """ Implementation of SteadyCom (Chan et al 2017).

    Args:
        community (CommunityModel): community model
        constraints (dict): environmental or additional constraints (optional)
        solver (Solver): solver instance instantiated with the model, for speed (optional)

    Returns:
        CommunitySolution: solution object
    """

    if solver is None:
        solver = build_problem(community)

    objective = {community.merged_model.biomass_reaction: 1}

    sol = binary_search(solver, objective, minimize=False, constraints=constraints)

    solution = CommunitySolution(community, sol)
    solution.solver = solver

    return solution


def SteadyComVA(community, obj_frac=1.0, constraints=None, solver=None):
    """ Abundance Variability Analysis using SteadyCom (Chan et al 2017).

    Args:
        community (CommunityModel): community model
        obj_frac (float): minimum fraction of the maximum growth rate (default 1.0)
        constraints (dict): environmental or additional constraints (optional)
        solver (Solver): solver instance instantiated with the model, for speed (optional)

    Returns:
        dict: species abundance variability
    """

    if solver is None:
        solver = build_problem(community)

    objective = {community.merged_model.biomass_reaction: 1}

    sol = binary_search(solver, objective, constraints=constraints)
    growth = obj_frac * sol.values[community.merged_model.biomass_reaction]
    solver.update_growth(growth)

    variability = {org_id: [None, None] for org_id in community.organisms}

    for org_id in community.organisms:
        sol2 = solver.solve({f"x_{org_id}": 1}, minimize=True, get_values=False, constraints=constraints)
        variability[org_id][0] = sol2.fobj

    for org_id in community.organisms:
        sol2 = solver.solve({f"x_{org_id}": 1}, minimize=False, get_values=False, constraints=constraints)
        variability[org_id][1] = sol2.fobj

    return variability


def build_problem(community, growth=1, bigM=1000):

    solver = solver_instance()
    model = community.merged_model

    # create biomass variables
    for org_id in community.organisms:
        solver.add_variable(f"x_{org_id}", 0, 1, update=False)

    # create all community reactions
    for r_id, reaction in model.reactions.items():
        if reaction.reaction_type == ReactionType.EXCHANGE:
            solver.add_variable(r_id, reaction.lb, reaction.ub, update=False)
        else:
            lb = -inf if reaction.lb < 0 else 0
            ub = inf if reaction.ub > 0 else 0
            solver.add_variable(r_id, lb, ub, update=False)

    solver.update()

    # sum biomass = 1
    solver.add_constraint("abundance", {f"x_{org_id}": 1 for org_id in community.organisms},
                          rhs=1, update=False)

    # S.v = 0
    table = model.metabolite_reaction_lookup()
    for m_id in model.metabolites:
        solver.add_constraint(m_id, table[m_id], update=False)

    # organism-specific constraints
    for org_id, organism in community.organisms.items():

        for r_id, reaction in organism.reactions.items():
            if (org_id, r_id) not in community.reaction_map:
                continue

            new_id = community.reaction_map[(org_id, r_id)]

            # growth = mu * X
            if r_id == organism.biomass_reaction:
                solver.add_constraint(f"g_{org_id}", {f"x_{org_id}": growth, new_id: -1}, update=False)
            # lb * X < R < ub * X
            else:
                lb = -bigM if isinf(reaction.lb) else reaction.lb
                ub = bigM if isinf(reaction.ub) else reaction.ub

                if lb != 0:
                    solver.add_constraint(f"lb_{new_id}", {f"x_{org_id}": lb, new_id: -1}, '<', 0, update=False)

                if ub != 0:
                    solver.add_constraint(f"ub_{new_id}", {f"x_{org_id}": ub, new_id: -1}, '>', 0, update=False)

    solver.update()

    def update_growth(value):
        # TODO: find a solution that is not CPLEX specific
        coefficients = [(f"g_{x}", f"x_{x}", value) for x in community.organisms]
        solver.problem.linear_constraints.set_coefficients(coefficients)

    solver.update_growth = update_growth

    return solver


def binary_search(solver, objective, obj_frac=1, minimize=False, max_iters=20, abs_tol=1e-3, constraints=None):

    previous_value = 0
    value = 1
    fold = 2
    feasible = False
    last_feasible = 0

    for i in range(max_iters):
        diff = value - previous_value

        if diff < abs_tol:
            break

        if feasible:
            last_feasible = value
            previous_value = value
            value = fold*diff + value
        else:
            if i > 0:
                fold = 0.5
            value = fold*diff + previous_value

        solver.update_growth(value)
        sol = solver.solve(objective, get_values=False, minimize=minimize, constraints=constraints)

        feasible = sol.status == Status.OPTIMAL

    if feasible:
        solver.update_growth(obj_frac * value)
    else:
        solver.update_growth(obj_frac * last_feasible)

    sol = solver.solve(objective, minimize=minimize, constraints=constraints)

    if i == max_iters - 1:
        warn("Max iterations exceeded.")

    return sol
