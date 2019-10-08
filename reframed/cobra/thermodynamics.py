
from ..solvers import solver_instance
from ..solvers.solver import VarType
from ..solvers.solution import Status
import numpy as np
from scipy.linalg import svd
from scipy import compress
from math import sqrt, inf

from warnings import warn


def nullspace(M, eps=1e-12):
    M = np.array(M)
    u, s, vh = svd(M)
    padding = M.shape[1]-s.shape[0]
    null_mask = np.concatenate(((s <= eps), np.ones((padding,), dtype=bool)), axis=0)
    N = compress(null_mask, vh, axis=0)
    return N


def llFBA(model, objective=None, minimize=False, constraints=None, internal=None, solver=None, get_values=True):
    """ Run a loopless FBA simulation (as defined in Schellenberger et al, 2011):

    Arguments:
        model (CBModel): a constraint-based model
        objective (dict): objective coefficients (optional)
        minimize (bool): sense of optimization (maximize by default)
        constraints (dict: environmental or additional constraints (optional)
        internal (list): list of internal reactions (optional)
        solver (Solver): solver instance instantiated with the model, for speed (optional)
        get_values (bool): set to false for speedup if you only care about the objective value (optional, default: True)

    Returns:
        Solution: solution
    """

    M = 1e4
    abstol = 1e-12

    if not solver:
        solver = solver_instance(model)

    if not objective:
        objective = model.get_objective()

    if not hasattr(solver, 'll_FBA_flag'):
        solver.ll_FBA_flag = True

        if not internal:
            internal = [r_id for r_id, reaction in model.reactions.items()
                        if len(reaction.stoichiometry) > 1]

        Sint = [[model.reactions[r_id].stoichiometry[m_id]
                 if m_id in model.reactions[r_id].stoichiometry else 0
                 for r_id in internal]
                for m_id in model.metabolites]

        Nint = nullspace(Sint)

        for r_id in internal:
            a, g = 'a' + r_id, 'g' + r_id
            solver.add_variable(g, update=False)
            solver.add_variable(a, 0, 1, vartype=VarType.BINARY, update=False)
        solver.update()

        for r_id in internal:
            a, g = 'a' + r_id, 'g' + r_id
            solver.add_constraint('c1_' + r_id, {a: M, r_id: -1}, '<', M, update=False)
            solver.add_constraint('c2_' + r_id, {a: -M, r_id: 1}, '<', 0, update=False)
            solver.add_constraint('c3_' + r_id, {a: M + 1, g: 1}, '>', 1, update=False)
            solver.add_constraint('c4_' + r_id, {a: M + 1, g: 1}, '<', M, update=False)
        solver.update()

        for i, row in enumerate(Nint):
            expr = {'g' + r_id: coeff for r_id, coeff in zip(internal, row) if abs(coeff) > abstol}
            solver.add_constraint(f'n{i}', expr, '=', 0, update=False)

        solver.update()

    if not constraints:
        constraints = dict()

    solution = solver.solve(objective, minimize=minimize, constraints=constraints, get_values=get_values)

    return solution


def TFA(model, deltaG0, sdeltaG0=None, concentrations=None, concentration_max=1e-2, concentration_min=1e-5,
        max_fold_change=10, excluded=None, temperature=298.15, objective=None, minimize=False, constraints=None,
        ignore_model_bounds=True, solver=None, get_values=True):
    """ Run a thermodynamic flux analysis (TFA) simulation.

        Implements a modified version of the TMFA method defined in Henry et al, 2007
        using reversible constraints (as defined by Hoppe et al, 2007) to avoid reaction decomposition.

    Args:
        model (CBModel): a constraint-based model
        deltaG0 (dict): deltaG0 values
        sdeltaG0 (dict): standard deviation of deltaG0 (optional)
        concentrations (dict): measured metabolite concentrations (optional)
        concentration_max (float): max concentration of unmeasured metabolites (default: 1e-2 M)
        concentration_min (float): min concentration of unmeasured metabolites (default: 1e-5 M)
        max_fold_change (float): allowable fold change of measured metabolites (default: 10)
        excluded (list): metabolites to exclude from calculations (ex: waters, protons, etc)
        temperature (float): default 298.15 Kelvin
        objective (dict): objective coefficients (optional)
        minimize (bool): sense of optimization (maximize by default)
        constraints (dict: environmental or additional constraints (optional)
        ignore_model_bounds: ignore model bounds for reactions with measured deltaG0 (default: True)
        solver (Solver): solver instance instantiated with the model, for speed (optional)
        get_values (bool): set to false for speedup if you only care about the objective value (optional, default: True)

    Returns:
        Solution: solution

    """

    bigM = 1000

    ln_min = np.log(concentration_min)
    ln_max = np.log(concentration_max)

    R = 0.00831
    RT = temperature * R

    if not solver:

        if ignore_model_bounds:
            model = model.copy()
            for r_id in model.reactions:
                if r_id in deltaG0:
                    model.set_flux_bounds(r_id, -bigM, bigM)

        solver = solver_instance(model)

    if not objective:
        objective = model.get_objective()

    if not concentrations:
        concentrations = {}

    if not sdeltaG0:
        sdeltaG0 = {r_id: 0 for r_id in deltaG0}

    if not excluded:
        excluded = []

    included = []

    if not hasattr(solver, 'tFBA_flag'):
        solver.tFBA_flag = True

        for r_id in model.reactions:
            if r_id in deltaG0:
                solver.add_variable('y_' + r_id, 0, 1, vartype=VarType.BINARY, update=False)
                solver.add_variable('dG_' + r_id, update=False)
                dG0_min, dG0_max = deltaG0[r_id] - sdeltaG0[r_id], deltaG0[r_id] + sdeltaG0[r_id]
                solver.add_variable('dG0_' + r_id, dG0_min, dG0_max, update=False)

                for m_id in model.reactions[r_id].stoichiometry:
                    if m_id not in excluded and m_id not in included:
                        if m_id in concentrations:
                            ln_min_j = np.log(concentrations[m_id] / sqrt(max_fold_change))
                            ln_max_j = np.log(concentrations[m_id] * sqrt(max_fold_change))
                        else:
                            ln_min_j, ln_max_j = ln_min, ln_max
                        solver.add_variable('ln_' + m_id, ln_min_j, ln_max_j, update=False)
                        included.append(m_id)

        solver.update()

        for r_id in model.reactions:
            if r_id in deltaG0:
                solver.add_constraint('lb_' + r_id, {r_id: 1, 'y_' + r_id: bigM}, '>', 0, update=False)
                solver.add_constraint('ub_' + r_id, {r_id: 1, 'y_' + r_id: bigM}, '<', bigM, update=False)
                solver.add_constraint('lb_dG_' + r_id, {'dG_' + r_id: -1, 'y_' + r_id: bigM}, '>', 0,
                                      update=False)
                solver.add_constraint('ub_dG_' + r_id, {'dG_' + r_id: -1, 'y_' + r_id: bigM}, '<', bigM,
                                      update=False)
                lhs = {'ln_' + m_id: RT * coeff for m_id, coeff in model.reactions[r_id].stoichiometry.items()
                       if m_id in included}
                lhs.update({'dG0_' + r_id: 1, 'dG_' + r_id: -1})
                solver.add_constraint('dGsum_' + r_id, lhs, '=', 0, update=False)

        solver.update()

    solution = solver.solve(objective, minimize=minimize, constraints=constraints, get_values=get_values)

    return solution


def TVA(model, deltaG0, sdeltaG0=None, concentrations=None, concentration_max=1e-2, concentration_min=1e-5,
        max_fold_change=10, excluded=None, temperature=298.15, constraints=None, ignore_model_bounds=True,
        reactions=None):
    """ Thermodynamic flux variability analysis (TVA) based on the thermodynamics flux analysis (TFA) method.

    Args:
        model (CBModel): a constraint-based model
        deltaG0 (dict): deltaG0 values
        sdeltaG0 (dict): standard deviation of deltaG0 (optional)
        concentrations (dict): measured metabolite concentrations (optional)
        concentration_max (float): max concentration of unmeasured metabolites (default: 1e-2 M)
        concentration_min (float): min concentration of unmeasured metabolites (default: 1e-5 M)
        max_fold_change (float): allowable fold change of measured metabolites (default: 10)
        excluded (list): metabolites to exclude from calculations (ex: waters, protons, etc)
        temperature (float): default 298.15 Kelvin
        constraints (dict): environmental or additional constraints (optional)
        ignore_model_bounds (bool): ignore model bounds for reactions with measured deltaG0 (default: True)
        reactions (list): list of reactions to analyse (default: all)

    Returns:
        dict: flux variation ranges

    """

    bigM = 1000

    if ignore_model_bounds:
        model = model.copy()

        for r_id in model.reactions:
            if r_id in deltaG0:
                model.set_flux_bounds(r_id, -bigM, bigM)

    solver = solver_instance(model)

    if not reactions:
        reactions = model.reactions.keys()

    bounds = {}

    for r_id in reactions:

        objective = {r_id: 1}
        sol = TFA(model, deltaG0,
                  sdeltaG0=sdeltaG0,
                  concentrations=concentrations,
                  concentration_max=concentration_max,
                  concentration_min=concentration_min,
                  max_fold_change=max_fold_change,
                  excluded=excluded,
                  temperature=temperature,
                  objective=objective,
                  minimize=True,
                  constraints=constraints,
                  solver=solver,
                  get_values=False)

        if sol.status == Status.INFEASIBLE:
            warn("Problem is infeasible")
            return

        lb = sol.fobj if sol.status == Status.OPTIMAL else -inf

        sol = TFA(model, deltaG0,
                  sdeltaG0=sdeltaG0,
                  concentrations=concentrations,
                  concentration_max=concentration_max,
                  concentration_min=concentration_min,
                  max_fold_change=max_fold_change,
                  excluded=excluded,
                  temperature=temperature,
                  objective=objective,
                  minimize=False,
                  constraints=constraints,
                  solver=solver,
                  get_values=False)

        ub = sol.fobj if sol.status == Status.OPTIMAL else inf

        bounds[r_id] = (lb, ub)

    return bounds


def NET(model, deltaG0, sdeltaG0=None, concentrations=None, concentration_max=1e-2, concentration_min=1e-5,
        fold_change=10, excluded=None, temperature=298.15, reaction_directions=None,
        get_dG_range=True, get_concentration_range=False, verbose=False):
    """ Implementation of Network-embedded thermodynamic (NET) analysis (Kummel et al 2006).

        Differs from the original implementation by (optionally) integrating standard deviation for deltaG0.

    Args:
        model (CBModel): a constraint-based model
        deltaG0 (dict): deltaG0 values
        sdeltaG0 (dict): standard deviation of deltaG0 (optional)
        concentrations (dict): measured metabolite concentrations (optional)
        concentration_max (float): max concentration of unmeasured metabolites (default: 1e-2 M)
        concentration_min (float): min concentration of unmeasured metabolites (default: 1e-5 M)
        fold_change (float): allowable fold change of measured metabolites (default: 10)
        excluded (list): metabolites to exclude from calculations (ex: waters, protons, etc)
        temperature (float): default 298.15 Kelvin
        reaction_directions (dict): pre-defined reaction directions: +1 forward, -1 negative (optional)
        get_dG_range (bool): calculate dG ranges (default: True)
        get_concentration_range (bool): calculate metabolite concentration ranges (default: False)
        verbose (bool): Print progress information if true
    Returns:

    """

    solver = solver_instance(model)

    ln_min = np.log(concentration_min)
    ln_max = np.log(concentration_max)

    R = 0.00831
    RT = temperature * R

    if not concentrations:
        concentrations = {}

    if not reaction_directions:
        reaction_directions = {}

    if not sdeltaG0:
        sdeltaG0 = {r_id: 0 for r_id in deltaG0}

    if not excluded:
        excluded = []

    included = []

    for r_id in model.reactions:
        if r_id in deltaG0:
            dG_min, dG_max = -inf, inf
            if r_id in reaction_directions:
                if reaction_directions[r_id] < 0:
                    dG_min, dG_max = 0, inf
                elif reaction_directions[r_id] > 0:
                    dG_min, dG_max = -inf, 0

            solver.add_variable('dG_' + r_id, dG_min, dG_max, update=False)
            dG0_min, dG0_max = deltaG0[r_id] - sdeltaG0[r_id], deltaG0[r_id] + sdeltaG0[r_id]
            solver.add_variable('dG0_' + r_id, dG0_min, dG0_max, update=False)

            for m_id in model.reactions[r_id].stoichiometry:
                if m_id not in excluded and m_id not in included:
                    if m_id in concentrations:
                        ln_min_j = np.log(concentrations[m_id] / sqrt(fold_change))
                        ln_max_j = np.log(concentrations[m_id] * sqrt(fold_change))
                    else:
                        ln_min_j, ln_max_j = ln_min, ln_max
                    solver.add_variable('ln_' + m_id, ln_min_j, ln_max_j, update=False)
                    included.append(m_id)

    solver.update()

    for r_id in model.reactions:
        if r_id in deltaG0:
            lhs = {'ln_' + m_id: RT * coeff for m_id, coeff in model.reactions[r_id].stoichiometry.items()
                   if m_id in included}
            lhs.update({'dG0_' + r_id: 1, 'dG_' + r_id: -1})
            solver.add_constraint('dGsum_' + r_id, lhs, '=', 0, update=False)

    solver.update()

    if get_dG_range:
        dG_range = {}

        for r_i, r_id in enumerate(model.reactions):

            if r_i % 10 == 0 and verbose:
                print(f"{i}/{len(model.reactions)}")

            if r_id in deltaG0:
                sol_min = solver.solve({'dG_' + r_id: 1}, minimize=True)
                if sol_min.status == Status.INFEASIBLE:
                    warn("Problem is infeasible")
                    return
                sol_max = solver.solve({'dG_' + r_id: 1}, minimize=False)
                dG_min = sol_min.fobj if sol_min.status == Status.OPTIMAL else None
                dG_max = sol_max.fobj if sol_max.status == Status.OPTIMAL else None
                dG_range[r_id] = (dG_min, dG_max)

    if get_concentration_range:
        x_range = {}

        for m_id in model.metabolites:
            if m_id in included:
                sol_min = solver.solve({'ln_' + m_id: 1}, minimize=True)
                if sol_min.status == Status.INFEASIBLE:
                    warn("Problem is infeasible")
                    return
                sol_max = solver.solve({'ln_' + m_id: 1}, minimize=False)
                x_min = np.exp(sol_min.fobj) if sol_min.status == Status.OPTIMAL else None
                x_max = np.exp(sol_max.fobj) if sol_max.status == Status.OPTIMAL else None
                x_range[m_id] = (x_min, x_max)

    if get_dG_range and get_concentration_range:
        return dG_range, x_range
    elif get_dG_range:
        return dG_range
    elif get_concentration_range:
        return x_range