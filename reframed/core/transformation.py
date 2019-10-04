from math import inf, isinf

from .model import Reaction
from .cbmodel import CBModel, CBReaction
from ..cobra.variability import blocked_reactions


def fix_reversibility(model):
    """ Make reaction reversibility consistent with the bounds. """

    for reaction in model.reactions.values():
        reaction.reversible = reaction.lb < 0


def clean_bounds(model, threshold=1000):
    """ Replace large bounds with infinity. """

    for reaction in model.reactions.values():
        if reaction.lb <= -threshold:
            reaction.lb = -inf
        if reaction.ub >= threshold:
            reaction.ub = inf


def apply_bounds(model, default_lb=-1000, default_ub=1000):
    """ Apply artificial bounds for unbounded reactions (opposite of `clean_bounds`). """

    for reaction in model.reactions.values():
        if isinf(reaction.lb):
            reaction.lb = default_lb

        if isinf(reaction.ub):
            reaction.ub = default_ub


def make_irreversible(model, inplace=True, reactions=None):
    """ Splits all reversible reactions into forward and backward directions.

    Arguments:
        model : Model (or CBmodel)
        inplace (bool): change model inplace (default), otherwise create a copy first
        reactions (list) : split only reactions in this list (optional)

    Returns:
        dict: mapping of old reaction ids to splitted reaction ids
    """

    if not inplace:
        model = model.copy()

    if reactions is None:
        reactions = model.reactions.keys()

    mapping = dict()

    for r_id, reaction in model.reactions.copy().items():
        if reaction.reversible and r_id in reactions:
            fwd_id = reaction.id + '_f'
            bwd_id = reaction.id + '_b'
            mapping[r_id] = (fwd_id, bwd_id)
            bwd_stoichiometry = [(m_id, -coeff) for m_id, coeff in reaction.stoichiometry.items()]

            if isinstance(model, CBModel):
                lb, ub = reaction.lb, reaction.ub
                lb_fwd = max(0, lb)
                ub_fwd = max(0, ub)
                lb_bwd = max(-ub, 0)
                ub_bwd = max(-lb, 0)
                obj = reaction.objective
                obj_fwd = obj if obj >= 0 else 0
                obj_bwd = -obj if obj < 0 else 0
                r_fwd = CBReaction(fwd_id, reaction.name, False, reaction.stoichiometry, reaction.regulators,
                                   lb_fwd, ub_fwd, obj_fwd, reaction.gpr)
                r_bwd = CBReaction(bwd_id, reaction.name, False, bwd_stoichiometry, reaction.regulators,
                                   lb_bwd, ub_bwd, obj_bwd, reaction.gpr)
            else:
                r_fwd = Reaction(fwd_id, reaction.name, False, reaction.stoichiometry, reaction.regulators)
                r_bwd = Reaction(bwd_id, reaction.name, False, bwd_stoichiometry, reaction.regulators)

            model.add_reaction(r_fwd)
            model.add_reaction(r_bwd)
            model.remove_reaction(r_id)

    if inplace:
        return mapping
    else:
        model.mapping = mapping
        return model


def simplify(model, reactions=None, clean_compartments=True, inplace=True):
    """ Removes all blocked reactions in a constraint based model

    Arguments:
        model (CBModel): model
        reactions (list): List of reactions which will be checked for being blocked (default: None - check all reactions)
        clean_compartments (bool): remove empty compartments (default: True)
        inplace (bool): change model in place (default), otherwise create a copy first

    Returns:
        CBModel: simplified model (if not in place)
    """

    if not inplace:
        model = model.copy()

    model.remove_reactions(blocked_reactions(model, reactions=reactions))
    model.remove_metabolites(disconnected_metabolites(model), safe_delete=False)
    model.remove_genes(disconnected_genes(model))
    if clean_compartments:
        model.remove_compartments(empty_compartments(model))

    if not inplace:
        return model


def disconnected_metabolites(model):
    m_r_table = model.metabolite_reaction_lookup()
    return [m_id for m_id, reactions in m_r_table.items() if not reactions]


def disconnected_genes(model):
    disconnected = set(model.genes)
    for reaction in model.reactions.values():
        disconnected -= set(reaction.get_genes())
    return disconnected


def empty_compartments(model):
    return set(model.compartments) - {met.compartment for met in model.metabolites.values()}
