from math import inf, isinf

from .model import Reaction, Compartment, Metabolite
from .cbmodel import CBModel, CBReaction, GPRAssociation
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
        model : Model (or CBModel)
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


def simplify(model, reactions=None, constraints=None, clean_compartments=True, inplace=True):
    """ Removes all blocked reactions in a constraint based model

    Arguments:
        model (CBModel): model
        reactions (list): List of reactions that can be tested for removal (default: test all reactions)
        constraints (dict): additional constraints (optional)
        clean_compartments (bool): remove empty compartments (default: True)
        inplace (bool): change model in place (default), otherwise create a copy first

    Returns:
        CBModel: simplified model (if not in place)
    """

    if not inplace:
        model = model.copy()

    model.remove_reactions(blocked_reactions(model, reactions=reactions, constraints=constraints))
    model.remove_metabolites(disconnected_metabolites(model), safe_delete=False)
    model.update()
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


def rename(model, inplace=True, fmt_mets=None, fmt_rxns=None, fmt_comps=None, fmt_genes=None):
    """
    Rename model identifiers using formatting functions that transform the old identifiers.

    Arguments:
        model (Model): model
        inplace (bool): change model in place (default), otherwise create a copy first
        fmt_mets (function): function to format metabolite identifiers
        fmt_rxns (function): function to format reaction identifiers
        fmt_comps (function): function to format compartment identifiers
        fmt_genes (function): function to format gene identifiers (only valid for instances of CBModel)
    """

    if not inplace:
        model = model.copy()

    if fmt_genes and not isinstance(model, CBModel):
        raise RuntimeError("fmt_genes can only be used with constraint-based models (CBModel instance)")

    if fmt_comps:
        for c_id in list(model.compartments.keys()):
            comp = model.compartments[c_id]
            comp.id = fmt_comps(c_id)
            model.compartments[comp.id] = comp
            del model.compartments[c_id]

    if fmt_genes:
        for g_id in list(model.genes.keys()):
            gene = model.genes[g_id]
            gene.id = fmt_genes(g_id)
            model.genes[gene.id] = gene
            del model.genes[g_id]

    for m_id in list(model.metabolites.keys()):
        met = model.metabolites[m_id]
        if fmt_comps:
            met.compartment = fmt_comps(met.compartment)
        if fmt_mets:
            met.id = fmt_mets(m_id)
            model.metabolites[met.id] = met
            del model.metabolites[m_id]

    for r_id in list(model.reactions.keys()):
        rxn = model.reactions[r_id]
        if fmt_mets:
            for m_id in list(rxn.stoichiometry.keys()):
                rxn.stoichiometry[fmt_mets(m_id)] = rxn.stoichiometry[m_id]
                del rxn.stoichiometry[m_id]
        if rxn.gpr is not None and fmt_genes:
            for protein in rxn.gpr.proteins:
                protein.genes = [fmt_genes(g_id) for g_id in protein.genes]
        if fmt_rxns:
            rxn.id = fmt_rxns(r_id)
            model.reactions[rxn.id] = rxn
            del model.reactions[r_id]

    if not inplace:
        return model


SPONTANEOUS = {'G_s0001', 'G_S0001', 'G_s_0001', 'G_S_0001', 'G_spontaneous', 'G_SPONTANEOUS', 'G_UNKNOWN',
               's0001', 'S0001', 's_0001', 'S_0001', 'spontaneous', 'SPONTANEOUS', 'UNKNOWN'}


def gpr_transform(model, inplace=True, add_proteome=False, gene_prefix='G_', usage_prefix='u_', pseudo_genes=None):
    """ Transformation method that integrates GPR associations directly into the stoichiometric matrix.

    Notes:
        This method extends the stoichiometric matrix where genes become pseudo-metabolites (rows)
        and enzyme usage variables become pseudo-reactions columns.

        See http://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1005140 for details.

    Args:
        model (CBModel): original model
        inplace (bool): change model in place (default), otherwise create a copy first
        add_proteome (bool): add a global proteome constraints (default: False)
        gene_prefix (str): prefix used for gene ids (default: 'G_')
        usage_prefix (str): prefix to create enzyme usage variables (default: 'u_')
        pseudo_genes (list): ignore pseudo/fake genes in the model (default: 's0001')

    Returns:
        CBModel: only when inplace=False
        list: list of newly created enzyme usage variables
    """

    if not inplace:
        model = model.copy()

    mapping_rev = make_irreversible(model)
    mapping_iso = split_isozymes(model)
    u_reactions = genes_to_species(model, add_proteome=add_proteome, gene_prefix=gene_prefix,
                                   usage_prefix=usage_prefix, pseudo_genes=pseudo_genes)

    model.convert_fluxes = lambda x, net=True: merge_fluxes(x, mapping_rev, mapping_iso, net)
    model.convert_constraints = lambda x: convert_constraints(x, mapping_rev, mapping_iso)
    model.convert_id_to_expr = lambda x, coeff=1: convert_id_to_expr(x, mapping_rev, mapping_iso, coeff=coeff)
    model.mapping_rev = mapping_rev
    model.mapping_iso = mapping_iso
    model.u_reactions = u_reactions

    if not inplace:
        return model


def split_isozymes(model):
    mapping = dict()

    for r_id, reaction in model.reactions.copy().items():

        if reaction.gpr is not None and len(reaction.gpr.proteins) > 1:
            mapping[r_id] = []
            for i, protein in enumerate(reaction.gpr.proteins):
                r_id_new = '{}_iso{}'.format(reaction.id, i + 1)
                mapping[r_id].append(r_id_new)
                gpr_new = GPRAssociation()
                gpr_new.proteins.append(protein)
                reaction_new = CBReaction(r_id_new, reaction.name,
                                          reversible=reaction.reversible,
                                          stoichiometry=reaction.stoichiometry,
                                          regulators=reaction.regulators,
                                          lb=reaction.lb, ub=reaction.ub,
                                          objective=reaction.objective,
                                          gpr_association=gpr_new)
                model.add_reaction(reaction_new)
            model.remove_reaction(r_id)

    return mapping


def genes_to_species(model, add_proteome=False, gene_prefix='G_', usage_prefix='u_', pseudo_genes=None):

    if pseudo_genes is None:
        pseudo_genes = SPONTANEOUS

    new_reactions = []
    compartment = Compartment('genes', 'gene pool')
    model.add_compartment(compartment)

    if add_proteome:
        proteome = Metabolite("proteome", "proteome", "genes")
        model.add_metabolite(proteome)
        r_synthesis = CBReaction("proteome_synth", "proteome synthesis", False, {"proteome": 1})
        model.add_reaction(r_synthesis)

    for gene in model.genes.values():
        if gene.id in pseudo_genes:
            continue
        model.add_metabolite(Metabolite(gene.id, gene.id, 'genes'))
        r_id = usage_prefix + gene.id[len(gene_prefix):]
        stoichiometry = {gene.id: 1}
        if add_proteome:
            stoichiometry["proteome"] = -1
        reaction = CBReaction(r_id, r_id, False, stoichiometry)
        model.add_reaction(reaction)
        new_reactions.append(r_id)

    for r_id, reaction in model.reactions.items():

        if reaction.gpr is not None:
            if len(reaction.gpr.proteins) > 1:
                print('error: isozymes not split:', r_id)
                return
            elif len(reaction.gpr.proteins) == 1:
                for g_id in reaction.gpr.proteins[0].genes:
                    if g_id not in pseudo_genes:
                        reaction.stoichiometry[g_id] = -1

    return new_reactions


def merge_fluxes(fluxes, mapping_rev, mapping_iso, net=True):

    fluxes = fluxes.copy()

    for r_id, r_ids in mapping_iso.items():
        fluxes[r_id] = sum([fluxes[r_id2] for r_id2 in r_ids])
        for r_id2 in r_ids:
            del fluxes[r_id2]

    for r_id, (fwd_id, bwd_id) in mapping_rev.items():
        if net:
            fluxes[r_id] = fluxes[fwd_id] - fluxes[bwd_id]
        else:
            fluxes[r_id] = fluxes[fwd_id] + fluxes[bwd_id]
        del fluxes[fwd_id]
        del fluxes[bwd_id]

    return fluxes


def convert_constraints(constraints, mapping_rev, mapping_iso):
    constraints = constraints.copy()

    for r_id, (fwd_id, bwd_id) in mapping_rev.items():
        if r_id in constraints:
            x = constraints[r_id]
            lb, ub = x if isinstance(x, tuple) else (x, x)
            lb_fwd = max(0, lb) if lb is not None else 0
            ub_fwd = max(0, ub) if ub is not None else None
            lb_bwd = max(-ub, 0) if ub is not None else 0
            ub_bwd = max(-lb, 0) if lb is not None else None
            constraints[fwd_id] = (lb_fwd, ub_fwd)
            constraints[bwd_id] = (lb_bwd, ub_bwd)
            del constraints[r_id]

    for r_id, r_ids in mapping_iso.items():
        if r_id in constraints:
            x = constraints[r_id]
            ub = x[1] if isinstance(x, tuple) else x
            for r_id2 in r_ids:
                constraints[r_id2] = (0, ub)
            del constraints[r_id]

    return constraints


def convert_id_to_expr(r_id, mapping_rev, mapping_iso, coeff=1):

    if r_id in mapping_rev:
        r_id_f, r_id_b = mapping_rev[r_id]
        expr1 = {r_id_f: coeff, r_id_b: -coeff}
    else:
        expr1 = {r_id: coeff}

    expr2 = {}

    for r_id2, val in expr1.items():
        if r_id2 in mapping_iso:
            expr2.update({r_id3: val for r_id3 in mapping_iso[r_id2]})
        else:
            expr2[r_id2] = val

    return expr2
