""" This module implements gene and reaction deletion methods.

Author: Daniel Machado

"""
from .simulation import FBA, pFBA, CAFBA, lMOMA, ROOM
from ..solvers import solver_instance
from ..solvers.solution import Status


def gene_knockout(model, genes, method='FBA', reference=None, constraints=None, solver=None, ignore_silent=False):
    """ Simulate the knockout of a set of genes.

    Arguments:
        model (CBModel): model
        genes (list): genes to delete
        method (str): any available simulation method: FBA (default), pFBA, CAFBA, lMOMA, ROOM
        reference (dict): reference flux distribution for lMOMA or ROOM (optional)
        constraints (dict): additional constraints
        solver (Solver): solver instance instantiated with the model, for speed (optional)
        ignore_silent (bool): ignore knockout if no reactions are affected (default: False)

    Returns:
        Solution: solution
    """

    if isinstance(genes, str):
        genes = [genes]

    inactive_reactions = deleted_genes_to_reactions(model, genes)

    if inactive_reactions or not ignore_silent:
        solution = reaction_knockout(model, inactive_reactions, method, reference, constraints, solver)
    else:
        solution = None

    return solution


def deleted_genes_to_reactions(model, genes):
    """ Convert a set of deleted genes to the respective deleted reactions.

    Arguments:
        model (CBModel): model
        genes (list): genes to delete

    Returns:
        list: list of deleted reactions
    """

    if isinstance(genes, str):
        genes = [genes]

    active_genes = set(model.genes) - set(genes)
    active_reactions = model.evaluate_gprs(active_genes)
    inactive_reactions = set(model.reactions) - set(active_reactions)

    return inactive_reactions


def hard_knockout(model, genes, inplace=False, rename=None):
    """ Create a mutant model with permanent gene knockouts.

    Arguments:
        model (CBModel): model
        genes (list): genes to delete
        inplace (bool): do not create a model copy (default: False)
        rename (str): give a new id to the model (optional)

    Returns:
        model: mutant model
    """

    if not inplace:
        model = model.copy()

    for r_id in deleted_genes_to_reactions(model, genes):
        model.set_flux_bounds(r_id, 0, 0)

    if rename is not None:
        model.id = rename

    if not inplace:
        return model


def reaction_knockout(model, reactions, method='FBA', reference=None, constraints=None, solver=None):
    """ Simulate the knockout of a set of reactions.

    Arguments:
        model (CBModel): model
        reactions (list): reactions to delete
        method (str): any available simulation method: FBA (default), pFBA, CAFBA, lMOMA, ROOM
        reference (dict): reference flux distribution for lMOMA or ROOM (optional)
        constraints (dict): additional constraints
        solver (Solver): solver instance instantiated with the model, for speed (optional)

    Returns:
        Solution: solution
    """

    if isinstance(reactions, str):
        reactions = [reactions]

    _constraints = {}

    if constraints:
        _constraints.update(constraints)

    for r_id in reactions:
        _constraints[r_id] = 0

    methods = {
        'FBA': FBA,
        'pFBA': pFBA,
        'CAFBA': CAFBA,
        'lMOMA': lMOMA,
        'ROOM': ROOM
    }

    if method in {'FBA', 'pFBA', 'CAFBA'}:
        function = methods[method]
        solution = function(model, constraints=_constraints, solver=solver)
    elif method in {'lMOMA', 'ROOM'}:
        function = methods[method]
        solution = function(model, reference, constraints=_constraints, solver=solver)
    else:
        raise RuntimeError(f"Method not available: {method}")

    return solution


def essential_genes(model, min_growth=0.01, constraints=None, solver=None):
    """ Compute the set of essential genes for the given experimental conditions.

    Arguments:
        model (CBModel): model
        min_growth (float): minimum fraction of growth rate to consider a deletion non-letal (default: 0.01)
        constraints (dict): environmental or additional constraints (optional)
        solver (Solver): solver instance instantiated with the model, for speed (optional)

    Returns:
        list: essential genes
    """
    return essentiality(model, 'genes', min_growth, constraints, solver)


def essential_reactions(model, min_growth=0.01, constraints=None, solver=None):
    """ Compute the set of essential reactions for the given experimental conditions.

    Arguments:
        model (CBModel): model
        min_growth (float): minimum fraction of growth rate to consider a deletion non-letal (default: 0.01)
        constraints (dict): environmental or additional constraints (optional)
        solver (Solver): solver instance instantiated with the model, for speed (optional)

    Returns:
        list: essential reactions
    """

    return essentiality(model, 'reactions', min_growth, constraints, solver)


def essentiality(model, kind, min_growth=0.01, constraints=None, solver=None):
    """ Generic interface for computing gene or reaction essentiality.

    Arguments:
        model (CBModel): model
        kind (str): 'genes' or 'reactions'
        min_growth (float): minimum fraction of growth rate to consider a deletion non-letal (default: 0.01)
        constraints (dict): environmental or additional constraints (optional)
        solver (Solver): solver instance instantiated with the model, for speed (optional)

    Returns:
        list: essential elements
    """

    if solver is None:
        solver = solver_instance(model)

    wt_solution = FBA(model, constraints=constraints, solver=solver)
    wt_growth = wt_solution.fobj

    if kind == 'genes':
        elements = model.genes
    else:
        kind = 'reactions'
        elements = model.reactions

    essential = []

    for elem in elements:
        if kind == 'genes':
            solution = gene_knockout(model, [elem], constraints=constraints, solver=solver, ignore_silent=True)
        else:
            solution = reaction_knockout(model, [elem], constraints=constraints, solver=solver)

        if (solution is not None
                and ((solution.status == Status.OPTIMAL and solution.fobj < min_growth * wt_growth)
                     or solution.status == Status.INFEASIBLE)):
            essential.append(elem)

    return essential
