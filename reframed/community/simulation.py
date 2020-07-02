from ..core.model import ReactionType
from ..solvers import solver_instance
from math import inf, isinf
from .solution import CommunitySolution
from ..solvers.solution import Status
from numpy.random import lognormal


def build_problem(community, growth=0.1, abundance=None):

    bigM = 1000
    solver = solver_instance()
    model = community.merged_model

    if growth is None and abundance is None:
        raise RuntimeError("If growth is a variable, species abundance must be given.")

    if abundance is None:  # create abundance variables
        for org_id in community.organisms:
            solver.add_variable(f"x_{org_id}", 0, 1, update=False)

    # add all community model reactions
    for r_id, reaction in model.reactions.items():
        if reaction.reaction_type == ReactionType.EXCHANGE:
            solver.add_variable(r_id, reaction.lb, reaction.ub, update=False)
        else:
            lb = -inf if reaction.lb < 0 else 0
            ub = inf if reaction.ub > 0 else 0
            solver.add_variable(r_id, lb, ub, update=False)

    solver.update()

    # sum biomass = 1
    if abundance is None:
        solver.add_constraint("abundance", {f"x_{org_id}": 1 for org_id in community.organisms}, rhs=1, update=False)

    # S.v = 0
    table = model.metabolite_reaction_lookup()
    for m_id in model.metabolites:
        solver.add_constraint(m_id, table[m_id], update=False)

    mu = model.biomass_reaction

    # organism-specific constraints
    for org_id, organism in community.organisms.items():

        for r_id, reaction in organism.reactions.items():
            if (org_id, r_id) not in community.reaction_map:
                continue

            new_id = community.reaction_map[(org_id, r_id)]

            # growth_i = mu * X_i
            if r_id == organism.biomass_reaction:

                if growth is None: # growth is variable, abundance is fixed
                    solver.add_constraint(f"g_{org_id}", {mu: abundance[org_id], new_id: -1}, update=False)
                elif abundance is None: # growth is fixed, abundance is variable
                    solver.add_constraint(f"g_{org_id}", {f"x_{org_id}": growth, new_id: -1}, update=False)
                else: # growth and abundance are fixed
                    solver.add_constraint(f"g_{org_id}", {new_id: 1}, '=', growth * abundance[org_id], update=False)

            # lb * X < R < ub * X
            else:
                lb = -bigM if isinf(reaction.lb) else reaction.lb
                ub = bigM if isinf(reaction.ub) else reaction.ub

                if lb != 0:

                    if abundance is None:
                        solver.add_constraint(f"lb_{new_id}", {f"x_{org_id}": lb, new_id: -1}, '<', 0, update=False)
                    else:
                        solver.add_constraint(f"lb_{new_id}", {new_id: 1}, '>', lb * abundance[org_id], update=False)

                if ub != 0:
                    if abundance is None:
                        solver.add_constraint(f"ub_{new_id}", {f"x_{org_id}": ub, new_id: -1}, '>', 0, update=False)
                    else:
                        solver.add_constraint(f"lb_{new_id}", {new_id: 1}, '<', ub * abundance[org_id], update=False)

    solver.update()

    return solver


def simulate(community, objective=None, growth=0.1, abundance=None, allocation=False, constraints=None,
             w_e=0.001, w_r=0.5, solver=None):

    if abundance:
        norm = sum(abundance.values())
        abundance = {org_id: abundance.get(org_id, 0) / norm for org_id in community.organisms}

    if not solver:
        solver = build_problem(community, growth=growth, abundance=abundance)

    if allocation:
        _ = allocation_constraints(community, solver, w_e=w_e, w_r=w_r, abundance=abundance)

    if not objective:
        objective = community.merged_model.biomass_reaction

    sol = solver.solve(objective, minimize=False, constraints=constraints)

    if sol.status == Status.OPTIMAL:
        sol = CommunitySolution(community, sol.values)

    return sol


def sample(community, n=100, growth=0.1, abundance=None, allocation=False, constraints=None,
             w_e=0.001, w_r=0.5, solver=None):

    if abundance:
        norm = sum(abundance.values())
        abundance = {org_id: abundance.get(org_id, 0) / norm for org_id in community.organisms}

    if not solver:
        solver = build_problem(community, growth=growth, abundance=abundance)

    if not allocation:
        w_e = 0
        w_r = 0

    enz_vars = allocation_constraints(community, solver, w_e=w_e, w_r=w_r, abundance=abundance)

    sols = []
    if abundance:
        w1 = {org_id: 1 / abundance[org_id] if abundance[org_id] > 0 else 0 for org_id in community.organisms}

    for _ in range(n):

        if not abundance:
            w1 = {org_id: lognormal(0, 1) for org_id in community.organisms}

        objective = {vi: w1[org_id] * lognormal(0, 1) for org_id, v_org in enz_vars.items() for vi in v_org}

        sol = solver.solve(objective, minimize=True, constraints=constraints)

        if sol.status == Status.OPTIMAL:
            sol = CommunitySolution(community, sol.values)
            sols.append(sol)

    return sols


def allocation_constraints(community, solver, w_e=0.001, w_r=0.5, abundance=None):
    enz_vars = {}
    split_vars = {}

    for org_id, organism in community.organisms.items():
        enz_vars[org_id] = []
        split_vars[org_id] = {}

        for r_id, reaction in organism.reactions.items():
            if (org_id, r_id) not in community.reaction_map:
                continue

            new_id = community.reaction_map[(org_id, r_id)]

            # test if reaction is enzymatic
            if reaction.gpr is not None and len(reaction.get_genes()) > 0:
                if reaction.reversible:
                    pos, neg = new_id + '+', new_id + '-'
                    solver.add_variable(pos, 0, inf, update=False)
                    solver.add_variable(neg, 0, inf, update=False)
                    enz_vars[org_id].append(pos)
                    enz_vars[org_id].append(neg)
                    split_vars[org_id][new_id] = (pos, neg)
                else:
                    enz_vars[org_id].append(new_id)

    solver.update()

    for org_id, organism in community.organisms.items():

        # constrain absolute values
        for r_id, (pos, neg) in split_vars[org_id].items():
            solver.add_constraint('c' + pos, {r_id: -1, pos: 1}, '>', 0, update=False)
            solver.add_constraint('c' + neg, {r_id: 1, neg: 1}, '>', 0, update=False)

        # protein allocation constraints
        alloc_constr = {r_id: w_e for r_id in enz_vars[org_id]}
        org_growth = community.reaction_map[(org_id, organism.biomass_reaction)]
        alloc_constr[org_growth] = w_r
        if abundance:
            solver.add_constraint(f"prot_{org_id}", alloc_constr, '<', abundance[org_id], update=False)
        else:
            alloc_constr[f"x_{org_id}"] = -1
            solver.add_constraint(f"prot_{org_id}", alloc_constr, '<', 0, update=False)

    solver.update()

    return enz_vars


def fit_abundance(community, abundance, growth=0.1, constraints=None, allocation=False,
             w_e=0.001, w_r=0.5, solver=None):

    if not solver:
        solver = build_problem(community, growth=growth)

    if allocation:
        _ = allocation_constraints(community, solver, w_e=w_e, w_r=w_r, abundance=None)

    abd_vars = []

    for org_id in community.organisms:
        d_pos, d_neg = f"d_{org_id}_+", f"d_{org_id}_-"
        solver.add_variable(d_pos, 0, 1, update=False)
        solver.add_variable(d_neg, 0, 1, update=False)
        abd_vars.extend([d_pos, d_neg])

    solver.update()

    norm = sum(abundance.values())
    for org_id in community.organisms:
        value = abundance.get(org_id, 0) / norm
        d_pos, d_neg = f"d_{org_id}_+", f"d_{org_id}_-"
        solver.add_constraint('c' + d_pos, {f"x_{org_id}": -1, d_pos: 1}, '>', -value, update=False)
        solver.add_constraint('c' + d_neg, {f"x_{org_id}": 1, d_neg: 1}, '>', value, update=False)

    solver.update()

    objective = {var: 1 for var in abd_vars}
    sol = solver.solve(objective, minimize=True, constraints=constraints)

    fitted = None

    if sol.status == Status.OPTIMAL:
        fitted = {org_id: sol.values[f"x_{org_id}"] for org_id in community.organisms}

    return fitted
