from collections import OrderedDict
from math import inf
from .model import Model, Metabolite, Reaction, Compartment, AttrOrderedDict
import re
from .parser import ReactionParser
from warnings import warn


class Gene(object):
    """ Base class for modeling genes. """

    def __init__(self, gene_id, name=None):
        """
        Arguments:
            gene_id (str): a valid unique identifier
            name (str): common gene name
        """
        self.id = gene_id
        self.name = name if name is not None else gene_id
        self.metadata = OrderedDict()

    def __str__(self):
        return self.name

    def __repr__(self):
        return str(self)


class Protein(object):
    """ Base class for modeling proteins.

        One protein is composed of a list of genes encoding one or more subunits.
    """

    def __init__(self):
        self.genes = []
        self.metadata = OrderedDict()

    def __str__(self):
        protein_str = ' and '.join(self.genes)

        if len(self.genes) > 1:
            protein_str = '(' + protein_str + ')'

        return protein_str

    def __repr__(self):
        return str(self)


class GPRAssociation(object):
    """ Base class for modeling Gene-Protein-Reaction associations.

        Each GPR association is composed by a list of proteins that can catalyze a reaction.
        Each protein is encoded by one or several genes.
    """

    def __init__(self):
        self.proteins = []
        self.metadata = OrderedDict()

    def __str__(self):

        gpr_str = ' or '.join(map(str, self.proteins))

        if len(self.proteins) > 1:
            gpr_str = '(' + gpr_str + ')'

        return gpr_str

    def __repr__(self):
        return str(self)

    def get_genes(self):
        """ Return the set of all associated genes. """

        return {gene for protein in self.proteins for gene in protein.genes}

    def remove_gene(self, gene_id):
        for protein in self.proteins:
            if gene_id in protein.genes:
                protein.genes.remove(gene_id)

        self.proteins = [protein for protein in self.proteins if len(protein.genes) > 0]


class CBReaction(Reaction):

    def __init__(self, reaction_id, name=None, reversible=True, stoichiometry=None, regulators=None,
                 lb=-inf, ub=inf, objective=0, gpr_association=None, reaction_type=None):

        Reaction.__init__(self, reaction_id, name=name, reversible=reversible, stoichiometry=stoichiometry,
                          regulators=regulators, reaction_type=reaction_type)

        self.lb = 0 if reversible == False and lb < 0 else lb
        self.ub = ub
        self.objective = objective
        self.gpr = gpr_association
        self._bool_function = None

    def set_flux_bounds(self, lb, ub):
        self.lb, self.ub = lb, ub
        self.reversible = bool(lb < 0)

    def set_gpr_association(self, gpr_association):
        self.gpr = gpr_association

    def set_objective(self, value):
        self.objective = value

    def get_genes(self):
        if self.gpr is not None:
            return self.gpr.get_genes()
        else:
            return []

    def evaluate_gpr(self, active_genes):
        """ Boolean evaluation of the GPR association for a given set of active genes.

        Arguments:
            active_genes (list): list of active genes

        Returns:
            bool: is the reaction active
        """

        if self._bool_function is None:
            self._gpr_to_function()

        return self._bool_function(active_genes)

    def _gpr_to_function(self):

        if not self.gpr:
            rule = 'True'
        else:
            rule = ' ' + str(self.gpr).replace('(', '( ').replace(')', ' )') + ' '
            for gene in self.get_genes():
                rule = rule.replace(' ' + gene + ' ', ' x[\'' + gene + '\'] ')
        self._bool_function = eval('lambda x: ' + rule)

    def to_string(self, metabolite_names=None):
        """ Print a reaction to a text based representation.

        Arguments:
            metabolite_names (dict): replace metabolite id's with names (optional)

        Returns:
            str: reaction string
        """

        str_rxn = Reaction.to_string(self, metabolite_names)

        if self.lb != -inf and (self.reversible or self.lb != 0.0) or self.ub != inf:
            str_rxn += f" [{self.lb}, {self.ub}]"

        if self.objective:
            str_rxn += f' @{self.objective}'

        return str_rxn

    def __str__(self):
        return self.to_string()

    def __repr__(self):
        return str(self)


class CBModel(Model):
    """ This class implements a constraint-based model."""

    def __init__(self, model_id):
        """
        Arguments:
            model_id (string): a valid unique identifier
        """
        Model.__init__(self, model_id)
        self.genes = AttrOrderedDict()
        self._g_r_lookup = None
        self._biomass_reaction = None

    def update(self):
        Model.update(self)
        self._g_r_lookup = None
        self._biomass_reaction = None

    @property
    def biomass_reaction(self):
        if self._biomass_reaction is None:
            self._detect_biomass_reaction()

        return self._biomass_reaction

    @biomass_reaction.setter
    def biomass_reaction(self, r_id):
        if r_id not in self.reactions:
            raise RuntimeError(f"Reaction {r_id} is not in the model")
        self._biomass_reaction = r_id

    def _detect_biomass_reaction(self):

            matches = [r_id for r_id, rxn in self.reactions.items() if rxn.objective]

            if matches:
                self._biomass_reaction = matches[0]
                if len(matches) > 1:
                    warn("Ambiguous biomass reaction (model has multiple objectives).")

                id_name = self._biomass_reaction + self.reactions[self._biomass_reaction].name
                if not re.search("biomass|growth", id_name, re.IGNORECASE):
                    warn(f"Suspicious biomass identifier: {self._biomass_reaction}")
            else:
                raise RuntimeError(f"No biomass reaction identified from model objective.")

    def add_gene(self, gene, replace=True):
        """ Add a gene metabolite to the model.
        If a gene with the same id exists, it will be replaced.

        Arguments:
            gene (Gene): gene
            replace (bool): replace previous gene with same id (default: True)
       """

        if gene.id in self.genes and not replace:
            warn(f"Gene {gene.id} already exists, ignoring.")
        else:
            self.genes[gene.id] = gene

    def add_reaction_from_str(self, reaction_str, compartment=None):
        """ Parse a reaction from a string and add it to the model.

        Arguments:
            reaction_str (str): string representation a the reaction
            compartment (str): reaction compartment id (optional)

        Notes:
            If the metabolites specified in the reaction are not yet in the model, they will be automatically added.
            If the compartment id is not given, it will use the first available compartment.
        """

        if not self._parser:
            self._parser = ReactionParser()

        if compartment is None:
            compartment = list(self.compartments.keys())[0]

        r_id, reversible, stoichiometry, lb, ub, obj_coeff = \
            self._parser.parse_reaction(reaction_str, kind='cb')

        for m_id in stoichiometry:
            if m_id not in self.metabolites:
                self.add_metabolite(Metabolite(m_id, m_id, compartment=compartment))

        reaction = CBReaction(r_id, r_id, reversible, stoichiometry, None, lb, ub, obj_coeff)
        self.add_reaction(reaction)
        self._needs_update = True

        return r_id

    def add_ratio_constraint(self, r_num, r_den, ratio):
        """ Add a flux ratio constraint to the model.

        Arguments:
            r_num (str): id of the numerator
            r_den (str): id of the denominator
            ratio (float): ratio value

        Returns:
            str : identifier of the pseudo-metabolite
        """

        if r_num not in self.reactions or r_den not in self.reactions:
            raise KeyError(f"Invalid reactions in ratio {r_num}/{r_den}")

        pseudo_c_id = "pseudo"
        pseudo_m_id = f"ratio_{r_num}_{r_den}"

        if pseudo_c_id not in self.compartments:
            self.add_compartment(Compartment(pseudo_c_id))

        self.add_metabolite(Metabolite(pseudo_m_id, compartment=pseudo_c_id))
        self.reactions[r_num].stoichiometry[pseudo_m_id] = 1
        self.reactions[r_den].stoichiometry[pseudo_m_id] = -ratio
        return pseudo_m_id

    def get_reactions_by_gene(self, g_id):
        """ Get a list of reactions associated with a given gene.

        Args:
            g_id (str): gene id

        Returns:
            list: reactions catalyzed by any proteins (or subunits) encoded by this gene
        """
        g_r_lookup = self.gene_to_reaction_lookup()
        return g_r_lookup[g_id]

    def get_objective(self):
        return {r_id: rxn.objective for r_id, rxn in self.reactions.items() if rxn.objective}

    def remove_gene(self, gene_id):
        """ Remove a gene from the model.

        Arguments:
            gene_id (str) : gene id
        """
        self.remove_genes([gene_id])

    def remove_genes(self, genes, safe_delete=True):
        """ Remove a list of genes from the model.
            safe_delete (bool): also remove genes from reaction associations (default: True)

        Arguments:
            genes (list) : gene ids
        """

        if safe_delete:
            g_r_lookup = self.gene_to_reaction_lookup()

        for gene_id in genes:
            if gene_id in self.genes:
                del self.genes[gene_id]
            else:
                warn(f"No such gene '{gene_id}'")

            if safe_delete:
                for r_id in g_r_lookup[gene_id]:
                    self.reactions[r_id].gpr.remove_gene(gene_id)

    def remove_ratio_constraint(self, r_num, r_den):
        """ Remove a flux ratio constraint from the model.

        Arguments:
            r_num (str): id of the numerator
            r_den (str): id of the denominator

        """

        pseudo_m_id = f"ratio_{r_num}_{r_den}"
        if pseudo_m_id in self.metabolites:
            self.remove_metabolite(pseudo_m_id)
        else:
            raise RuntimeError(f"No ratio constraint for {r_num}/{r_den}")

    def set_flux_bounds(self, r_id, lb=None, ub=None):
        """ Define flux bounds for one reaction

        Arguments:
            r_id (str): reaction id
            lb (float): lower bound
            ub (float): upper bound
        """
        if r_id not in self.reactions:
            warn(f"Reaction {r_id} not found")
            return

        if lb is not None:
            self.reactions[r_id].lb = lb
            self.reactions[r_id].reversible = bool(lb < 0)

        if ub is not None:
            self.reactions[r_id].ub = ub

    def set_gpr_association(self, r_id, gpr, add_genes=True):
        """ Set GPR association for a given reaction:

        Arguments:
            r_id (str): reaction id
            gpr (GPRAssociation): GPR association
            add_genes (bool): check if associated genes need to be added to the model
        """

        if r_id not in self.reactions:
            raise KeyError(f"Reaction {r_id} not found")

        self.reactions[r_id].gpr = gpr

        if add_genes and gpr is not None:
            for gene_id in gpr.get_genes():
                if gene_id not in self.genes:
                    self.add_gene(Gene(gene_id))

    def set_objective(self, coefficients):
        """ Define objective coefficients for a list of reactions

        Arguments:
            coefficients (dict): dictionary of reactions and coefficients

        """
        for r_id, coeff, in coefficients.items():
            if r_id not in self.reactions:
                raise KeyError(f"Reaction {r_id} not found")
            self.reactions[r_id].objective = coeff

    def gene_to_reaction_lookup(self):
        """ Build a dictionary from genes to associated reactions.

        Returns:
            dict: gene to reaction mapping

        """
        if not self._g_r_lookup:
            self._g_r_lookup = {g_id: [] for g_id in self.genes}

            for r_id, rxn in self.reactions.items():
                genes = rxn.get_genes()
                for g_id in genes:
                    self._g_r_lookup[g_id].append(r_id)

        return self._g_r_lookup

    def evaluate_gprs(self, active_genes):
        """ Boolean evaluation of the GPR associations for a given set of active genes.

        Arguments:
            active_genes (list): list of active genes

        Returns:
            list: list of active reactions
        """
        genes_state = {gene: gene in active_genes for gene in self.genes}
        return [r_id for r_id, rxn in self.reactions.items() if rxn.evaluate_gpr(genes_state)]

