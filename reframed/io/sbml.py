import libsbml as sb
from ..core.model import Model, Metabolite, Reaction, Compartment, ReactionType, RegulatorType
from ..core.cbmodel import CBModel, Gene, Protein, GPRAssociation, CBReaction
from ..core.transformation import fix_reversibility, clean_bounds
from enum import Enum
from collections import OrderedDict
from math import inf, isinf, isnan
from sympy.parsing.sympy_parser import parse_expr
from sympy import to_dnf, Or, And
from sympy.logic.boolalg import is_dnf
import os
import re
from html import escape
from warnings import warn


class Flavor(Enum):
    COBRA = 'cobra'  # legacy cobra format
    FBC2 = 'fbc2'  # sbml-fbc2 format
    BIGG = 'bigg'  # fbc2 with BiGG notation


class CobraTags(Enum):
    LB_TAG = 'LOWER_BOUND'
    UB_TAG = 'UPPER_BOUND'
    OBJ_TAG = 'OBJECTIVE_COEFFICIENT'
    GPR_TAG = 'GENE_ASSOCIATION'


class CobraDefaults(Enum):
    LOWER_BOUND = -1000
    UPPER_BOUND = 1000
    LOWER_BOUND_ID = 'cobra_default_lb'
    UPPER_BOUND_ID = 'cobra_default_ub'
    ZERO_BOUND_ID = 'cobra_0_bound'


class SBO(Enum):
    ACTIVATOR_TAG = 'SBO:0000459'
    INHIBITOR_TAG = 'SBO:0000020'


class ExchangeDetection(Enum):
    PATTERN = 'pattern'
    UNBALANCED = 'unbalanced'
    BOUNDARY = 'boundary'


DEFAULT_SBML_LEVEL = 3
DEFAULT_SBML_VERSION = 1

non_alphanum = re.compile(r'\W+')
re_type = type(non_alphanum)


def load_sbml(filename):
    """ Loads an SBML file.

        Arguments:
            filename (str): SBML file path

        Returns:
            SBMLModel
    """

    if not os.path.exists(filename):
        raise IOError("Model file was not found")

    reader = sb.SBMLReader()
    document = reader.readSBML(str(filename))
    sbml_model = document.getModel()

    if sbml_model is None:
        document.printErrors()
        raise IOError(f'Failed to load model {filename}.')

    return sbml_model


def load_model(filename):
    """ Loads an basic metabolic model.

        Arguments:
            filename (str): SBML file path

        Returns:
            Model
    """

    sbml_model = load_sbml(filename)
    model = Model(sbml_model.getId())
    load_compartments(sbml_model, model)
    load_metabolites(sbml_model, model)
    load_reactions(sbml_model, model)

    return model


def load_cbmodel(filename, flavor=Flavor.FBC2.value, exchange_detection=None, external_compartment=True,
                 load_gprs=True, load_metadata=True, reversibility_check=True, use_infinity=True):
    """ Loads a constraint-based model.

    Arguments:
        filename (str): SBML file path
        flavor (str): adapt to different modeling conventions (optional, see Notes)
        exchange_detection (str): detect exchange reactions (optional, see Notes)
        external_compartment (bool or str): identify external compartment (optional, see Notes)
        load_gprs (bool): load GPR associations (default: True)
        load_metadata (bool): load metadata from annotations field (default: True)
        reversibility_check (bool): fix consistency between reversibility attribute and lower bounds (default: True)
        use_infinity (bool): replace large bounds with +/- infinity (default: True)

    Notes:
        Currently supported flavors:
            * 'cobra': legacy format from the cobra toolbox
            * 'fbc2': sbml-fbc2 format (default)
            * 'bigg': fbc2 using BiGG conventions

        Supported exchange detection modes:
            * 'unbalanced': Reactions that either have only substrates or products
            * 'boundary': Detect using the SBML metabolite boundary value
            * <regular expression>: Matches the reaction ID with a regular expression

        Identifying the external compartment:
            * default: do not identify
            * True: identify using exchange reaction detection
            * str: manually provide the compartment id

        Returns:
            CBModel
    """

    if exchange_detection is None:
        if flavor == Flavor.BIGG.value:
            exchange_detection = re.compile(r'^R_EX_')
        else:
            exchange_detection = 'unbalanced'
    elif exchange_detection not in {'unbalanced', 'boundary'}:
        try:
            exchange_detection = re.compile(exchange_detection)
        except:
            raise RuntimeError("Exchange detection must be: 'unbalanced', 'boundary', or a valid regular expression.")

    sbml_model = load_sbml(filename)
    model = CBModel(sbml_model.getId())
    load_compartments(sbml_model, model, load_metadata)
    load_metabolites(sbml_model, model, flavor, load_metadata)
    load_reactions(sbml_model, model, True, load_metadata)

    if flavor == Flavor.COBRA.value:
        load_cobra_bounds(sbml_model, model)
        load_cobra_objective(sbml_model, model)
        if load_gprs:
            load_cobra_gpr(sbml_model, model)

    elif flavor in {Flavor.BIGG.value, Flavor.FBC2.value}:
        load_fbc2_bounds(sbml_model, model)
        load_fbc2_objective(sbml_model, model)
        if load_gprs:
            load_fbc2_gpr(sbml_model, model)
    else:
        warn(f"Invalid flavor {flavor}. Current options are: 'cobra', 'fbc2', and 'bigg'.")

    for r_id, reaction in model.reactions.items():
        reaction.reaction_type = reaction_type_detection(r_id, model, sbml_model, exchange_detection)

    if exchange_detection is not None and len(model.get_exchange_reactions()) == 0:
        warn("Exchange reactions were not detected.")

    if external_compartment is not None:
        detect_external_compartment(model, external_compartment)

    if reversibility_check:
        fix_reversibility(model)

    if use_infinity:
        clean_bounds(model)

    return model


def detect_external_compartment(model, external_compartment):
    if isinstance(external_compartment, str):
        if external_compartment in model.compartments:
            model.compartments[external_compartment].external = True
        else:
            raise RuntimeError(f"Compartment {external_compartment} not in the model.")
    elif isinstance(external_compartment, bool) and external_compartment is True:

        m_r_lookup = model.metabolite_reaction_lookup()
        ext_comps = [model.metabolites[m_id].compartment
                     for m_id, r_ids in m_r_lookup.items() for r_id in r_ids
                     if model.reactions[r_id].reaction_type == ReactionType.EXCHANGE]

        ext_comp = max(set(ext_comps), key=ext_comps.count)

        model.compartments[ext_comp].external = True


def load_compartments(sbml_model, model, load_metadata=True):
    for compartment in sbml_model.getListOfCompartments():
        model.add_compartment(load_compartment(compartment, load_metadata=load_metadata))


def load_compartment(compartment, load_metadata=True):
    size = compartment.getSize()
    if isnan(size) or isinf(size):
        size = 1.0

    comp = Compartment(compartment.getId(), compartment.getName(), False, size)

    if load_metadata:
        extract_metadata(compartment, comp)
    return comp


def load_metabolites(sbml_model, model, flavor=None, load_metadata=True):
    for species in sbml_model.getListOfSpecies():
        model.add_metabolite(load_metabolite(species, flavor, load_metadata=load_metadata))


def load_metabolite(species, flavor=None, load_metadata=True):
    metabolite = Metabolite(species.getId(), species.getName(), species.getCompartment())

    if flavor in {Flavor.BIGG, Flavor.FBC2}:
        fbc_species = species.getPlugin('fbc')
        if fbc_species.isSetChemicalFormula():
            formula = fbc_species.getChemicalFormula()
            metabolite.metadata['FORMULA'] = formula

        if fbc_species.isSetCharge():
            charge = fbc_species.getCharge()
            metabolite.metadata['CHARGE'] = str(charge)

    if load_metadata:
        extract_metadata(species, metabolite)

    return metabolite


def load_reactions(sbml_model, model, cb=False, load_metadata=True):
    for reaction in sbml_model.getListOfReactions():
        rxn = load_reaction(reaction, cb, load_metadata)
        model.add_reaction(rxn)


def load_reaction(reaction, cb=False, load_metadata=True):
    stoichiometry = OrderedDict()
    modifiers = OrderedDict()
    substrates = []
    products = []

    for reactant in reaction.getListOfReactants():
        m_id = reactant.getSpecies()
        substrates.append(m_id)
        coeff = -reactant.getStoichiometry()

        if m_id not in stoichiometry:
            stoichiometry[m_id] = coeff
        else:
            stoichiometry[m_id] += coeff

    for product in reaction.getListOfProducts():
        m_id = product.getSpecies()
        products.append(m_id)
        coeff = product.getStoichiometry()

        if m_id not in stoichiometry:
            stoichiometry[m_id] = coeff
        else:
            stoichiometry[m_id] += coeff
        if stoichiometry[m_id] == 0.0:
            del stoichiometry[m_id]

    for modifier in reaction.getListOfModifiers():
        m_id = modifier.getSpecies()
        sboterm = modifier.getSBOTermID()

        if sboterm == SBO.ACTIVATOR_TAG:
            kind = RegulatorType.ACTIVATOR
        elif sboterm == SBO.INHIBITOR_TAG:
            kind = RegulatorType.INHIBITOR
        else:
            kind = RegulatorType.UNKNOWN
        modifiers[m_id] = kind

    if cb:
        rxn = CBReaction(reaction.getId(), name=reaction.getName(), reversible=reaction.getReversible(),
                         stoichiometry=stoichiometry, regulators=modifiers)
    else:
        rxn = Reaction(reaction.getId(), name=reaction.getName(), reversible=reaction.getReversible(),
                       stoichiometry=stoichiometry, regulators=modifiers)

    if load_metadata:
        extract_metadata(reaction, rxn)

    return rxn


def reaction_type_detection(r_id, model, sbml_model, exchange_detection):

    reaction = model.reactions[r_id]
    substrates = reaction.get_substrates()
    products = reaction.get_products()

    # test exchange reaction
    if exchange_detection == ExchangeDetection.UNBALANCED.value:
        if len(substrates) * len(products) == 0:
            return ReactionType.EXCHANGE

    elif exchange_detection == ExchangeDetection.BOUNDARY.value:
        if all(sbml_model.getSpecies(m_id).getBoundaryCondition() for m_id in substrates):
            return ReactionType.EXCHANGE

        if all(sbml_model.getSpecies(m_id).getBoundaryCondition() for m_id in products):
            return ReactionType.EXCHANGE

    elif isinstance(exchange_detection, re_type):
        if exchange_detection.search(r_id) is not None:
            return ReactionType.EXCHANGE

    # test sink reactions (unbalanced and internal)
    if len(substrates) * len(products) == 0:
        return ReactionType.SINK

    # test transport reactions (multiple compartments)
    if len(model.get_reaction_compartments(r_id)) > 1:
        return ReactionType.TRANSPORT

    # test enzymatic reactions (GPR associated)
    if len(reaction.get_genes()) > 0:
        return ReactionType.ENZYMATIC

    return ReactionType.OTHER


def load_cobra_bounds(sbml_model, model):
    for reaction in sbml_model.getListOfReactions():
        default_lb = -inf if reaction.getReversible() else 0
        lb = get_cb_parameter(reaction, CobraTags.LB_TAG.value, default_lb)
        ub = get_cb_parameter(reaction, CobraTags.UB_TAG.value)
        model.set_flux_bounds(reaction.getId(), lb, ub)


def load_cobra_objective(sbml_model, model):
    objective = OrderedDict()
    for reaction in sbml_model.getListOfReactions():
        coeff = get_cb_parameter(reaction, CobraTags.OBJ_TAG.value, default_value=0)
        if coeff:
            objective[reaction.getId()] = coeff
    model.set_objective(objective)


def get_cb_parameter(reaction, tag, default_value=None):
    param_value = default_value
    kinetic_law = reaction.getKineticLaw()
    if kinetic_law:
        parameter = kinetic_law.getParameter(tag)
        if parameter:
            param_value = parameter.getValue()
    return param_value


def load_cobra_gpr(sbml_model, model):
    genes = set()
    gprs = OrderedDict()

    for reaction in sbml_model.getListOfReactions():
        rxn = model.reactions[reaction.getId()]
        rule = rxn.metadata.pop(CobraTags.GPR_TAG.value, None)
        if rule:
            gpr = parse_gpr_rule(rule, prefix='G_')
            for protein in gpr.proteins:
                genes |= set(protein.genes)
            gprs[reaction.getId()] = gpr
        else:
            gprs[reaction.getId()] = None

    for gene in sorted(genes):
        model.add_gene(Gene(gene, gene[2:]))

    for r_id, gpr in gprs.items():
        model.set_gpr_association(r_id, gpr, add_genes=False)


def sanitize_id(identifier):
    return non_alphanum.sub('_', identifier)


def parse_gpr_rule(rule, prefix=None):
    if not rule:
        return None

    rule = rule.replace('(', '( ').replace(')', ' )')

    def replacement(token):
        if token.lower() == 'and':
            return '&'
        elif token.lower() == 'or':
            return '|'
        elif token == '(' or token == ')':
            return token
        elif prefix is not None and not token.startswith(prefix):
            return prefix + sanitize_id(token)
        else:
            return sanitize_id(token)

    rule = ' '.join(map(replacement, rule.split()))

    expr = parse_expr(rule)

    if not is_dnf(expr):
        expr = to_dnf(expr)

    gpr = GPRAssociation()

    if type(expr) is Or:
        for sub_expr in expr.args:
            protein = Protein()
            if type(sub_expr) is And:
                protein.genes = [str(gene) for gene in sub_expr.args]
            else:
                protein.genes = [str(sub_expr)]
            gpr.proteins.append(protein)
    elif type(expr) is And:
        protein = Protein()
        protein.genes = [str(gene) for gene in expr.args]
        gpr.proteins = [protein]
    else:
        protein = Protein()
        protein.genes = [str(expr)]
        gpr.proteins = [protein]

    return gpr


def load_fbc2_bounds(sbml_model, model):
    params = {param.getId(): param.getValue() for param in sbml_model.getListOfParameters()}

    for reaction in sbml_model.getListOfReactions():
        fbc_rxn = reaction.getPlugin('fbc')
        lb = fbc_rxn.getLowerFluxBound()
        ub = fbc_rxn.getUpperFluxBound()
        model.set_flux_bounds(reaction.getId(), params[lb], params[ub])


def load_fbc2_objective(sbml_model, model):
    fbcmodel = sbml_model.getPlugin('fbc')
    active_obj = fbcmodel.getActiveObjective()
    objective = OrderedDict()
    for rxn_obj in active_obj.getListOfFluxObjectives():
        r_id = rxn_obj.getReaction()
        coeff = rxn_obj.getCoefficient()
        if coeff:
            objective[r_id] = coeff
    model.set_objective(objective)


def load_fbc2_gpr(sbml_model, model):
    fbcmodel = sbml_model.getPlugin('fbc')

    for gene in fbcmodel.getListOfGeneProducts():
        model.add_gene(Gene(gene.getId(), gene.getName()))

    for reaction in sbml_model.getListOfReactions():
        fbcrxn = reaction.getPlugin('fbc')
        gpr_assoc = fbcrxn.getGeneProductAssociation()
        if gpr_assoc:
            gpr = parse_fbc_association(gpr_assoc.getAssociation(), reaction.getId())
            model.set_gpr_association(reaction.getId(), gpr, add_genes=False)
        else:
            model.set_gpr_association(reaction.getId(), None)


def parse_fbc_association(gpr_assoc, reaction_id):
    gpr = GPRAssociation()
    parsing_error = False

    if gpr_assoc.isFbcOr():
        for item in gpr_assoc.getListOfAssociations():
            protein = Protein()
            if item.isFbcAnd():
                for subitem in item.getListOfAssociations():
                    if subitem.isGeneProductRef():
                        protein.genes.append(subitem.getGeneProduct())
                    else:
                        parsing_error = True
            elif item.isGeneProductRef():
                protein.genes.append(item.getGeneProduct())
            else:
                parsing_error = True
            gpr.proteins.append(protein)
    elif gpr_assoc.isFbcAnd():
        protein = Protein()
        for item in gpr_assoc.getListOfAssociations():
            if item.isGeneProductRef():
                protein.genes.append(item.getGeneProduct())
            else:
                parsing_error = True
        gpr.proteins = [protein]
    elif gpr_assoc.isGeneProductRef():
        protein = Protein()
        protein.genes = [gpr_assoc.getGeneProduct()]
        gpr.proteins = [protein]
    else:
        parsing_error = True

    if parsing_error:
        warn(f"Gene association for reaction {reaction_id} is not DNF")
    else:
        return gpr


def extract_metadata(sbml_elem, elem):
    notes = sbml_elem.getNotes()

    if notes:
        recursive_node_parser(notes, elem.metadata)


def recursive_node_parser(node, cache):
    node_data = node.getCharacters()
    if ':' in node_data:
        key, value = node_data.split(':', 1)
        cache[key.strip()] = value.strip()

    for i in range(node.getNumChildren()):
        recursive_node_parser(node.getChild(i), cache)


def save_model(model, filename):
    """ Save a model to an SBML file.

    Arguments:
        model (Model): model
        filename (str): file path
    """

    document = sb.SBMLDocument(DEFAULT_SBML_LEVEL, DEFAULT_SBML_VERSION)
    sbml_model = document.createModel(model.id)
    save_compartments(model, sbml_model)
    save_metabolites(model, sbml_model)
    save_reactions(model, sbml_model)
    writer = sb.SBMLWriter()
    writer.writeSBML(document, filename)


def save_cbmodel(model, filename, flavor=Flavor.FBC2.value):
    """ Save a constraint-based model to an SBML file.

    Arguments:
        model (Model): model
        filename (str): file path
        flavor (str): (optional, currently available: 'cobra', 'fbc2', 'bigg')
    """

    document = sb.SBMLDocument(DEFAULT_SBML_LEVEL, DEFAULT_SBML_VERSION)
    sbml_model = document.createModel(model.id)

    if flavor in {Flavor.BIGG.value, Flavor.FBC2.value}:
        document.enablePackage(sb.FbcExtension.getXmlnsL3V1V2(), 'fbc', True)
        fbc_model = sbml_model.getPlugin('fbc')
        fbc_model.setStrict(True)
        document.setPackageRequired('fbc', False)

    save_compartments(model, sbml_model)
    save_metabolites(model, sbml_model, flavor)
    save_reactions(model, sbml_model)
    save_cb_parameters(model, sbml_model, flavor)
    save_gpr_associations(model, sbml_model, flavor)

    save_metadata(model, sbml_model)
    writer = sb.SBMLWriter()
    writer.writeSBML(document, filename)


def save_compartments(model, sbml_model):
    for compartment in model.compartments.values():
        sbml_compartment = sbml_model.createCompartment()
        sbml_compartment.setId(compartment.id)
        sbml_compartment.setName(compartment.name)
        sbml_compartment.setSize(compartment.size)
        sbml_compartment.setConstant(True)
        save_metadata(compartment, sbml_compartment)


def save_metabolites(model, sbml_model, flavor=None):
    for metabolite in model.metabolites.values():
        species = sbml_model.createSpecies()
        species.setId(metabolite.id)
        species.setName(metabolite.name)
        species.setCompartment(metabolite.compartment)
        species.setHasOnlySubstanceUnits(True)

        if flavor in {Flavor.BIGG.value, Flavor.FBC2.value}:
            fbc_species = species.getPlugin('fbc')

            if 'FORMULA' in metabolite.metadata:
                try:
                    fbc_species.setChemicalFormula(metabolite.metadata['FORMULA'])
                except:
                    pass
            if 'CHARGE' in metabolite.metadata:
                try:
                    charge = int(metabolite.metadata['CHARGE'])
                    fbc_species.setCharge(charge)
                except:
                    pass

        save_metadata(metabolite, species)


def save_reactions(model, sbml_model):
    for reaction in model.reactions.values():
        sbml_reaction = sbml_model.createReaction()
        sbml_reaction.setId(reaction.id)
        sbml_reaction.setName(reaction.name)
        sbml_reaction.setReversible(reaction.reversible)
        sbml_reaction.setFast(False)
        save_metadata(reaction, sbml_reaction)

        for m_id, coeff in reaction.stoichiometry.items():
            if coeff < 0:
                speciesReference = sbml_reaction.createReactant()
                speciesReference.setSpecies(m_id)
                speciesReference.setStoichiometry(-coeff)
                speciesReference.setConstant(True)
            elif coeff > 0:
                speciesReference = sbml_reaction.createProduct()
                speciesReference.setSpecies(m_id)
                speciesReference.setStoichiometry(coeff)
                speciesReference.setConstant(True)
        for m_id, kind in reaction.regulators.items():
            speciesReference = sbml_reaction.createModifier()
            speciesReference.setSpecies(m_id)
            if kind == RegulatorType.ACTIVATOR:
                speciesReference.setSBOTerm(SBO.ACTIVATOR_TAG)
            if kind == RegulatorType.INHIBITOR:
                speciesReference.setSBOTerm(SBO.INHIBITOR_TAG)


def save_cb_parameters(model, sbml_model, flavor):
    if flavor in {Flavor.BIGG.value, Flavor.FBC2.value}:
        save_fbc_fluxbounds(model, sbml_model)
        save_fbc_objective(model, sbml_model)
    else:
        save_cobra_parameters(model, sbml_model)


def save_gpr_associations(model, sbml_model, flavor):
    if flavor in {Flavor.BIGG.value, Flavor.FBC2.value}:
        save_fbc_gprs(model, sbml_model)
    else:
        save_cobra_gprs(model, sbml_model)


def save_cobra_parameters(model, sbml_model):
    for r_id, reaction in model.reactions.items():
        sbml_reaction = sbml_model.getReaction(r_id)
        kineticLaw = sbml_reaction.createKineticLaw()
        kineticLaw.setFormula('0')

        lb = CobraDefaults.LOWER_BOUND.value if isinf(reaction.lb) else reaction.lb
        lbParameter = kineticLaw.createParameter()
        lbParameter.setId(CobraTags.LB_TAG.value)
        lbParameter.setValue(lb)

        ub = CobraDefaults.UPPER_BOUND.value if isinf(reaction.ub) else reaction.ub
        ubParameter = kineticLaw.createParameter()
        ubParameter.setId(CobraTags.UB_TAG.value)
        ubParameter.setValue(ub)

        objParameter = kineticLaw.createParameter()
        objParameter.setId(CobraTags.OBJ_TAG.value)
        objParameter.setValue(reaction.objective)


def save_cobra_gprs(model, sbml_model):
    for r_id, reaction in model.reactions.items():
        if reaction.gpr:
            reaction.metadata[CobraTags.GPR_TAG.value] = str(reaction.gpr)
            sbml_reaction = sbml_model.getReaction(r_id)
            save_metadata(reaction, sbml_reaction)


def save_fbc_fluxbounds(model, sbml_model):
    default_lb = sbml_model.createParameter()
    default_lb.setId(CobraDefaults.LOWER_BOUND_ID.value)
    default_lb.setValue(CobraDefaults.LOWER_BOUND.value)
    default_lb.setConstant(True)

    default_ub = sbml_model.createParameter()
    default_ub.setId(CobraDefaults.UPPER_BOUND_ID.value)
    default_ub.setValue(CobraDefaults.UPPER_BOUND.value)
    default_ub.setConstant(True)

    zero_bound = sbml_model.createParameter()
    zero_bound.setId(CobraDefaults.ZERO_BOUND_ID.value)
    zero_bound.setValue(0)
    zero_bound.setConstant(True)

    for r_id, reaction in model.reactions.items():
        fbcrxn = sbml_model.getReaction(r_id).getPlugin('fbc')

        if reaction.lb <= CobraDefaults.LOWER_BOUND.value:
            fbcrxn.setLowerFluxBound(CobraDefaults.LOWER_BOUND_ID.value)
        elif reaction.lb == 0:
            fbcrxn.setLowerFluxBound(CobraDefaults.ZERO_BOUND_ID.value)
        else:
            lb_id = f"{r_id}_lower_bound"
            lb_param = sbml_model.createParameter()
            lb_param.setId(lb_id)
            lb_param.setValue(reaction.lb)
            lb_param.setConstant(True)
            fbcrxn.setLowerFluxBound(lb_id)

        if reaction.ub >= CobraDefaults.UPPER_BOUND.value:
            fbcrxn.setUpperFluxBound(CobraDefaults.UPPER_BOUND_ID.value)
        elif reaction.ub == 0:
            fbcrxn.setUpperFluxBound(CobraDefaults.ZERO_BOUND_ID.value)
        else:
            ub_id = f"{r_id}_upper_bound"
            ub_param = sbml_model.createParameter()
            ub_param.setId(ub_id)
            ub_param.setValue(reaction.ub)
            ub_param.setConstant(True)
            fbcrxn.setUpperFluxBound(ub_id)


def save_fbc_objective(model, sbml_model):
    fbcmodel = sbml_model.getPlugin('fbc')
    obj = fbcmodel.createObjective()
    obj.setId('objective')
    fbcmodel.setActiveObjectiveId('objective')
    obj.setType('maximize')
    for r_id, reaction in model.reactions.items():
        if reaction.objective:
            r_obj = obj.createFluxObjective()
            r_obj.setReaction(r_id)
            r_obj.setCoefficient(reaction.objective)


def save_fbc_gprs(model, sbml_model):
    fbcmodel = sbml_model.getPlugin('fbc')
    for gene in model.genes.values():
        gene_prod = fbcmodel.createGeneProduct()
        gene_prod.setId(gene.id)
        gene_prod.setName(gene.name)
        gene_prod.setLabel(gene.name)

    for r_id, reaction in model.reactions.items():
        if reaction.gpr:
            fbcrxn = sbml_model.getReaction(r_id).getPlugin('fbc')
            gpr_assoc = fbcrxn.createGeneProductAssociation()

            if len(reaction.gpr.proteins) > 1:
                gpr_assoc = gpr_assoc.createOr()

            for protein in reaction.gpr.proteins:
                if len(protein.genes) > 1:
                    protein_assoc = gpr_assoc.createAnd()
                else:
                    protein_assoc = gpr_assoc

                for gene in protein.genes:
                    gene_ref = protein_assoc.createGeneProductRef()
                    gene_ref.setGeneProduct(gene)


def save_metadata(elem, sbml_elem):
    if elem.metadata:
        try:
            notes = [f'<p>{key}: {escape(value)}</p>'
                     for key, value in elem.metadata.items()]
            note_string = '<html>' + ''.join(notes) + '</html>'
            note_xml = sb.XMLNode.convertStringToXMLNode(note_string)
            note_xml.getNamespaces().add('http://www.w3.org/1999/xhtml')
            sbml_elem.setNotes(note_xml)
        except AttributeError:
            warn(f"Unable to save metadata for object {sbml_elem.getId()}")
