from .solution import CommunitySolution
from ..core.elements import molecular_weight
from ..solvers.solution import Status
from ..core.model import ReactionType
from ..solvers import solver_instance
from warnings import warn
from math import inf, isinf
from random import lognormvariate

SPONTANEOUS = {'G_s0001', 'G_S0001', 'G_s_0001', 'G_S_0001', 'G_spontaneous', 'G_SPONTANEOUS',
               's0001', 'S0001', 's_0001', 'S_0001', 'spontaneous', 'SPONTANEOUS'}


def SteadierCom(community, objective1=None, objective2=None, growth=None, abundance=None, proteome=False,
                constraints=None, solver=None, w_e=0.001, w_r=0.5, obj1_tol=0.01):
    if objective2 is None:
        objectives = [objective1]
    else:
        objectives = [objective1, objective2]

    if "abundance" in objectives and abundance is None:
        raise RuntimeError("Experimental abundance values must be provided when using the abundance objective.")

    if solver is None:
        solver = build_problem(community, growth=growth, abundance=abundance, proteome=proteome,
                               min_uptake=("uptake" in objectives), parsimony=("parsimony" in objectives), w_e=w_e,
                               w_r=w_r)

    for i, objective in enumerate(objectives):
        if objective is None or objective == "growth":
            obj_func = {community.merged_model.biomass_reaction: 1}
            minimize = False
        elif objective == "abundance":
            obj_func = {x: 1 for x in solver.abd_vars}
            minimize = True
        elif objective == "parsimony":
            if abundance is None:
                obj_func = {r_id: 1 for org_vars in solver.enz_vars.values() for r_id in org_vars}
            else:
                obj_func = {r_id: 1 / abundance[org_id] for org_id, org_vars in solver.enz_vars.items()
                            for r_id in org_vars}
            minimize = True
        elif objective == "uptake":
            obj_func = solver.upt_vars
            minimize = True
        elif isinstance(objective, dict):
            obj_func = objective
            minimize = False
        else:
            raise RuntimeError(f"Invalid objective: {objective}.")

        if growth is None:
            sol = binary_search(solver, obj_func, minimize=minimize, constraints=constraints)
        else:
            sol = solver.solve(obj_func, minimize=minimize, constraints=constraints)

        #        solver.write_to_file(f"obj{i+1}.lp")

        if sol.status == Status.OPTIMAL:
            if i == len(objectives) - 1:
                solution = CommunitySolution(community, sol.values)
                solution.solver = solver
            else:
                if minimize:
                    solver.add_constraint("obj1", obj_func, '<', sol.fobj * (1 + obj1_tol))
                else:
                    solver.add_constraint("obj1", obj_func, '>', sol.fobj * (1 - obj1_tol))
                solver.update()
        else:
            #            warn("Failed to find optimal solution.")
            return

    return solution


def SteadierComVA(community, growth=None, obj_frac=1, proteome=False, w_e=0.001, w_r=0.5, constraints=None):
    solver = build_problem(community, proteome=proteome, w_e=w_e, w_r=w_r)

    objective = {community.merged_model.biomass_reaction: 1}

    if growth is None:
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


def SteadierComSample(community, n=10, growth=None, obj_frac=1, proteome=False, w_e=0.001, w_r=0.5, constraints=None):
    solver = build_problem(community, proteome=proteome, w_e=w_e, w_r=w_r)

    objective = {community.merged_model.biomass_reaction: 1}

    if growth is None:
        sol = binary_search(solver, objective, constraints=constraints)
        growth = obj_frac * sol.values[community.merged_model.biomass_reaction]

    solver.update_growth(growth)

    sols = []

    for _ in range(n):
        objective = {f"x_{org_id}": lognormvariate(0, 1) for org_id in community.organisms}
        sol = solver.solve(objective, minimize=False, constraints=constraints)
        sols.append(CommunitySolution(community, sol.values))

    return sols


def build_problem(community, growth=None, abundance=None, proteome=False, min_uptake=False, parsimony=False,
                  w_e=0.001, w_r=0.5, bigM=1000):
    solver = solver_instance()
    model = community.merged_model

    if growth is None:
        growth = 1

    if abundance is not None:
        norm = sum(abundance.values())
        abundance = {org_id: val / norm for org_id, val in abundance.items()}

    # temporary variables for abundance constraints
    abd_vars = []
    solver.abd_vars = abd_vars

    # create biomass variables
    for org_id in community.organisms:
        solver.add_variable(f"x_{org_id}", 0, 1, update=False)

        # temporary variables for abundance constraints
        if abundance and org_id in abundance:
            d_pos, d_neg = f"d_{org_id}_+", f"d_{org_id}_-"
            solver.add_variable(d_pos, 0, 1, update=False)
            solver.add_variable(d_neg, 0, 1, update=False)
            abd_vars.extend([d_pos, d_neg])

    # create all community reactions
    for r_id, reaction in model.reactions.items():
        if reaction.reaction_type == ReactionType.EXCHANGE:
            solver.add_variable(r_id, reaction.lb, reaction.ub, update=False)
        else:
            lb = -inf if reaction.lb < 0 else 0
            ub = inf if reaction.ub > 0 else 0
            solver.add_variable(r_id, lb, ub, update=False)

    # temporary variables for unidirectional values of uptake fluxes
    upt_vars = {}
    solver.upt_vars = upt_vars

    if min_uptake:
        for r_id in model.get_exchange_reactions():
            if model.reactions[r_id].lb < 0:
                weight = get_mol_weight(model, r_id)
                if weight is not None:
                    ub = -model.reactions[r_id].lb
                    solver.add_variable('f_' + r_id, 0, ub, update=False)
                    upt_vars['f_' + r_id] = weight

    # temporary variables for computing absolute values of enzymatic reactions
    enz_vars = {}
    solver.enz_vars = enz_vars
    tmp = {}

    if proteome or parsimony:
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

        if proteome or parsimony:
            # constrain absolute values
            for r_id in tmp[org_id]:
                pos, neg = r_id + '+', r_id + '-'
                solver.add_constraint('c' + pos, {r_id: -1, pos: 1}, '>', 0, update=False)
                solver.add_constraint('c' + neg, {r_id: 1, neg: 1}, '>', 0, update=False)

        if proteome:
            # protein allocation constraints
            alloc_constr = {r_id: w_e for r_id in enz_vars[org_id]}
            org_growth = community.reaction_map[(org_id, organism.biomass_reaction)]
            alloc_constr[org_growth] = w_r
            alloc_constr[f"x_{org_id}"] = -1
            solver.add_constraint(f"prot_{org_id}", alloc_constr, '<', 0, update=True)

        if abundance and org_id in abundance:
            d_pos, d_neg = f"d_{org_id}_+", f"d_{org_id}_-"
            solver.add_constraint('c' + d_pos, {f"x_{org_id}": -1, d_pos: 1}, '>', -abundance[org_id], update=False)
            solver.add_constraint('c' + d_neg, {f"x_{org_id}": 1, d_neg: 1}, '>', abundance[org_id], update=False)

    # constrain uptake fluxes to negative part of exchange reactions
    if min_uptake:
        for f_id, weight in upt_vars.items():
            r_id = f_id[2:]
            solver.add_constraint('c_' + r_id, {r_id: 1, f_id: 1}, '>', 0, update=False)

    solver.update()

    def update_growth(value):
        # TODO: find a solution that is not CPLEX specific
        coefficients = [(f"g_{x}", f"x_{x}", value) for x in community.organisms]
        solver.problem.linear_constraints.set_coefficients(coefficients)

    solver.update_growth = update_growth

    return solver


def get_mol_weight(model, r_id):
    compounds = model.reactions[r_id].get_substrates()
    metabolite = model.metabolites[compounds[0]]
    formula = metabolite.metadata.get('FORMULA', '')
    return molecular_weight(formula)


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
            value = fold * diff + value
        else:
            if i > 0:
                fold = 0.5
            value = fold * diff + previous_value

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
