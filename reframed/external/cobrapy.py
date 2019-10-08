from ..core.model import Compartment, Metabolite
from ..core.cbmodel import CBModel, CBReaction, Gene
from ..io.sbml import parse_gpr_rule


def to_cobrapy(model):
    """ Convert model to a *cobrapy.Model* object

    Arguments:
        model (CBModel): constraint-based model

    Returns:
        cobrapy.Model
    """

    try:
        import cobra as cb
    except ImportError:
        raise RuntimeError("CobraPy is not installed.")

    cb_model = cb.Model(model.id)
    cb_model.compartments = {comp.id: comp.name for comp in model.compartments.values()}

    cb_mets = []
    for met in model.metabolites.values():
        cb_met = cb.Metabolite(id=met.id, name=met.name, compartment=met.compartment,
                               formula=met.metadata.get('FORMULA', ''))
        cb_mets.append(cb_met)
    cb_model.add_metabolites(cb_mets)

    cb_rxns = []
    for rxn in model.reactions.values():
        cb_rxn = cb.Reaction(id=rxn.id, name=rxn.name, lower_bound=rxn.lb, upper_bound=rxn.ub)
        cb_rxns.append(cb_rxn)
    cb_model.add_reactions(cb_rxns)

    for r_id, rxn in model.reactions.items():
        cb_rxn = cb_model.reactions.get_by_id(r_id)

        cb_rxn.add_metabolites(rxn.stoichiometry)
        if rxn.gpr is not None:
            cb_rxn.gene_reaction_rule = str(rxn.gpr)

    cb_model.objective = {cb_model.reactions.get_by_id(r_id): coeff for r_id, coeff in model.get_objective().items()}

    return cb_model


def from_cobrapy(cb_model):
    """ Create a model from a *cobrapy.Model* object

    Arguments:
        cb_model (cobrapy.Model): constraint-based model

    Returns:
        CBModel
    """

    try:
        import cobra as cb
    except ImportError:
        raise RuntimeError("CobraPy is not installed.")

    model = CBModel(cb_model.id)

    for c_id, name in cb_model.compartments.items():
        comp = Compartment(c_id, name)
        model.add_compartment(comp)

    for cb_met in cb_model.metabolites:
        met = Metabolite(cb_met.id, cb_met.name, cb_met.compartment)
        model.add_metabolite(met)

    for cb_gene in cb_model.genes:
        gene = Gene(cb_gene.id, cb_gene.name)
        model.add_gene(gene)

    for cb_rxn in cb_model.reactions:
        stoichiometry = {cb_met.id: val for cb_met, val in cb_rxn.metabolites.items()}
        gpr = parse_gpr_rule(cb_rxn.gene_reaction_rule)

        rxn = CBReaction(cb_rxn.id, name=cb_rxn.name, reversible=(cb_rxn.lower_bound < 0),
                         stoichiometry=stoichiometry, lb=cb_rxn.lower_bound, ub=cb_rxn.upper_bound,
                         objective=cb_rxn.objective_coefficient, gpr_association=gpr)
        model.add_reaction(rxn)

    return model


