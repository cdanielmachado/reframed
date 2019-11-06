from .solution import CommunitySolution
from ..core.elements import molecular_weight
from ..solvers.solution import Status
from ..core.model import ReactionType
from ..solvers import solver_instance
from warnings import warn
from math import inf, isinf

SPONTANEOUS = {'G_s0001', 'G_S0001', 'G_s_0001', 'G_S_0001', 'G_spontaneous', 'G_SPONTANEOUS',
               's0001', 'S0001', 's_0001', 'S_0001', 'spontaneous', 'SPONTANEOUS'}


def SteadyCom(community, constraints=None, solver=None):
    return SteadierCom(community, constraints=constraints, w_e=0, w_r=0, solver=solver)


def SteadyComVA(community, obj_frac=1, constraints=None, solver=None):
    return SteadierComVA(community, obj_frac=obj_frac, constraints=constraints, w_e=0, w_r=0, solver=solver)


def SteadierCom(community, fixed_growth=None, w_e=0.001, w_r=0.5, min_uptake=False, parsimony=False,
                constraints=None, solver=None):

    if fixed_growth is None:
        growth = 1
    else:
        growth = fixed_growth

    if solver is None:
        solver = build_problem(community, growth=growth, w_e=w_e, w_r=w_r, min_uptake=min_uptake, parsimony=parsimony)

    objective = {}
    minimize = True

    if min_uptake:
        tmp_obj = min_uptake_objective(community.merged_model)
        objective.update(tmp_obj)

    if parsimony:
        tmp_obj = {r_id: 1 for org_vars in solver.enz_vars.values() for r_id in org_vars}
        objective.update(tmp_obj)

    if not objective:
        objective = {community.merged_model.biomass_reaction: 1}
        minimize = False

    if fixed_growth is None:
        sol = binary_search(solver, objective, minimize=minimize, constraints=constraints)
    else:
        sol = solver.solve(objective, minimize=minimize, constraints=constraints)

    solution = CommunitySolution(community, sol)
    solution.solver = solver

    if min_uptake:
        solution.total_uptake = sum(-sol.values[r_id]*objective["f_" + r_id] for r_id in solver.upt_vars)

    return solution


def SteadierComVA(community, fixed_growth=None, obj_frac=1, w_e=0.001, w_r=0.5, constraints=None, solver=None):

    if solver is None:
        solver = build_problem(community, w_e=w_e, w_r=w_r)

    objective = {community.merged_model.biomass_reaction: 1}

    if fixed_growth is None:
        sol = binary_search(solver, objective, constraints=constraints)
        growth = obj_frac * sol.values[community.merged_model.biomass_reaction]
    else:
        growth = fixed_growth

    solver.update_growth(growth)

    variability = {org_id: [None, None] for org_id in community.organisms}

    for org_id in community.organisms:
        sol2 = solver.solve({f"x_{org_id}": 1}, minimize=True, get_values=False, constraints=constraints)
        variability[org_id][0] = sol2.fobj

    for org_id in community.organisms:
        sol2 = solver.solve({f"x_{org_id}": 1}, minimize=False, get_values=False, constraints=constraints)
        variability[org_id][1] = sol2.fobj

    return variability


def build_problem(community, growth=1, w_e=0.001, w_r=0.5, min_uptake=False, parsimony=False, bigM=1000):

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

    # temporary variables for unidirectional values of uptake fluxes
    upt_vars = []
    solver.upt_vars = upt_vars

    if min_uptake:
        for r_id in model.get_exchange_reactions():
            ub = -model.reactions[r_id].lb
            if ub > 0:
                solver.add_variable('f_' + r_id, 0, ub, update=False)
                upt_vars.append(r_id)

    # temporary variables for computing absolute values of enzymatic reactions
    enz_vars = {}
    solver.enz_vars = enz_vars
    tmp = {}

    if w_e > 0 or parsimony:
        for org_id, organism in community.organisms.items():
            enz_vars[org_id] = []
            tmp[org_id] = []

            for r_id, reaction in organism.reactions.items():
                if (org_id, r_id) not in community.reaction_map:
                    continue

                new_id = community.reaction_map[(org_id, r_id)]

                # test if reaction is enzymatic
                if reaction.gpr is not None and len(set(reaction.get_genes()) & SPONTANEOUS) == 0:
                    if reaction.reversible:
                        pos, neg = new_id + '+', new_id + '-'
                        solver.add_variable(pos, 0, inf, update=False)
                        solver.add_variable(neg, 0, inf, update=False)
                        enz_vars[org_id].append(pos)
                        enz_vars[org_id].append(neg)
                        tmp[org_id].append(new_id)
                    else:
                        enz_vars[org_id].append(new_id)

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

        if w_e > 0 or parsimony:
            # constrain absolute values
            for r_id in tmp[org_id]:
                pos, neg = r_id + '+', r_id + '-'
                solver.add_constraint('c' + pos, {r_id: -1, pos: 1}, '>', 0, update=False)
                solver.add_constraint('c' + neg, {r_id: 1, neg: 1}, '>', 0, update=False)

            # protein allocation constraints
            alloc_constr = {r_id: w_e for r_id in enz_vars[org_id]}
            org_growth = community.reaction_map[(org_id, organism.biomass_reaction)]
            alloc_constr[org_growth] = w_r
            alloc_constr[f"x_{org_id}"] = -1
            solver.add_constraint(f"prot_{org_id}", alloc_constr, '<', 0, update=True)

    # constrain uptake fluxes to negative part of exchange reactions
    if min_uptake:
        for r_id in upt_vars:
            solver.add_constraint('c_' + r_id, {r_id: 1, 'f_' + r_id: 1}, '>', 0, update=False)

    solver.update()

    def update_growth(value):
        # TODO: find a solution that is not CPLEX specific
        coefficients = [(f"g_{x}", f"x_{x}", value) for x in community.organisms]
        solver.problem.linear_constraints.set_coefficients(coefficients)

    solver.update_growth = update_growth

    return solver


def min_uptake_objective(model):

    objective = {}

    for r_id in model.get_exchange_reactions():

        if model.reactions[r_id].lb < 0:
            compounds = model.reactions[r_id].get_substrates()
            metabolite = model.metabolites[compounds[0]]
            formula = metabolite.metadata.get('FORMULA', '')
            weight = molecular_weight(formula)

            if weight is not None:
                objective['f_' + r_id] = weight

    return objective


def binary_search(solver, objective, obj_frac=1, minimize=False, max_iters=100, abs_tol=1e-3, constraints=None):

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
