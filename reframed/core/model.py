from enum import Enum
from collections import OrderedDict
from copy import copy, deepcopy
from .parser import ReactionParser
from warnings import warn


class Compartment(object):
    """ Base class for modeling compartments. """

    def __init__(self, comp_id, name=None, external=False, size=1.0):
        """
        Arguments:
            comp_id (str): a valid unique identifier
            name (str): compartment name (optional)
            external (bool): is external (default: false)
            size (float): compartment size (default: 1.0)
        """
        self.id = comp_id
        self.name = name if name is not None else comp_id
        self.size = size
        self.external = external
        self.metadata = OrderedDict()

    def __str__(self):
        return self.name


class Metabolite(object):
    """ Base class for modeling metabolites. """

    def __init__(self, met_id, name=None, compartment=None):
        """
        Arguments:
            met_id (str): a valid unique identifier
            name (str): common metabolite name
            compartment (str): compartment containing the metabolite
        """
        self.id = met_id
        self.name = name if name is not None else met_id
        self.compartment = compartment
        self.metadata = OrderedDict()

    def __str__(self):
        return self.name


class ReactionType(Enum):
    """ Enumeration of possible reaction types. """
    ENZYMATIC = 0
    TRANSPORT = 1
    EXCHANGE = 2
    SINK = 3
    OTHER = 4


class RegulatorType(Enum):
    """ Enumeration of possible reaction regulator types. """
    ACTIVATOR = 0
    INHIBITOR = 1
    UNKNOWN = 2


class Reaction(object):
    """ Base class for modeling reactions. """

    def __init__(self, reaction_id, name=None, reversible=True, stoichiometry=None, regulators=None,
                 reaction_type=None):
        """
        Arguments:
            reaction_id (str): a valid unique identifier
            name (str): common reaction name
            reversible (bool): reaction reversibility (default: True)
            stoichiometry (dict): stoichiometry
            regulators (dict): reaction regulators
            reaction_type (ReactionType): reaction type
        """
        self.id = reaction_id
        self.name = name if name is not None else reaction_id
        self.reversible = reversible
        self.reaction_type = reaction_type if reaction_type is not None else ReactionType.OTHER
        self.stoichiometry = OrderedDict()
        self.regulators = OrderedDict()
        self.metadata = OrderedDict()

        if stoichiometry:
            self.stoichiometry.update(stoichiometry)
        if regulators:
            self.regulators.update(regulators)

    def __str__(self):
        return self.to_string()

    def get_substrates(self):
        """ Get list of reaction substrates

        Returns:
            list: reaction substrates
        """

        return [m_id for m_id, coeff in self.stoichiometry.items() if coeff < 0]

    def get_products(self):
        """ Get list of reaction products

        Returns:
            list: reaction products
        """

        return [m_id for m_id, coeff in self.stoichiometry.items() if coeff > 0]

    def get_activators(self):
        """ Get list of reaction activators

        Returns:
            list: reaction activators
        """

        return [m_id for m_id, kind in self.regulators.items() if kind == RegulatorType.ACTIVATOR]

    def get_inhibitors(self):
        """ Get list of reaction inhibitors

        Returns:
            list: reaction inhibitors
        """

        return [m_id for m_id, kind in self.regulators.items() if kind == RegulatorType.INHIBITOR]

    def to_equation(self, metabolite_names=None):
        """ Returns reaction equation string

        Arguments:
            metabolite_names (dict): replace metabolite id's with names (optional)

        Returns:
            str: reaction string
        """

        if metabolite_names:
            def met_repr(m_id):
                return metabolite_names[m_id]
        else:
            def met_repr(m_id):
                return m_id

        left = ' + '.join(met_repr(m_id) if coeff == -1.0 else str(-coeff) + ' ' + met_repr(m_id)
                          for m_id, coeff in self.stoichiometry.items() if coeff < 0)
        arrow = '<->' if self.reversible else '-->'
        right = ' + '.join(met_repr(m_id) if coeff == 1.0 else str(coeff) + ' ' + met_repr(m_id)
                           for m_id, coeff in self.stoichiometry.items() if coeff > 0)
        return f"{left} {arrow} {right}"

    def to_string(self, metabolite_names=None):
        """ Returns reaction as a string

        Arguments:
            metabolite_names (dict): replace metabolite id's with names (optional)

        Returns:
            str: reaction string
        """
        return self.id + ': ' + self.to_equation(metabolite_names=metabolite_names)


class AttrOrderedDict(OrderedDict):
    """Helper class to extend ordered dictionaries with indexing"""

    def __init__(self, *args, **nargs):
        super(AttrOrderedDict, self).__init__(*args, **nargs)

    def __getattr__(self, name):
        if not name.startswith('_'):
            return self[name]
        super(AttrOrderedDict, self).__getattr__(name)

    def __setattr__(self, name, value):
        if not name.startswith('_'):
            self[name] = value
        else:
            super(AttrOrderedDict, self).__setattr__(name, value)

    def __dir__(self):
        return dir(OrderedDict) + list(self.keys())

    def __copy__(self):
        my_copy = AttrOrderedDict()
        for key, val in self.items():
            my_copy[key] = copy(val)
        return my_copy

    def __deepcopy__(self, memo):
        my_copy = AttrOrderedDict()
        for key, val in self.items():
            my_copy[key] = deepcopy(val)
        return my_copy


class Model(object):
    """ Base class for all metabolic models."""

    def __init__(self, model_id):
        """
        Arguments:
            model_id (str): a valid unique identifier
        """
        self.id = model_id
        self.metabolites = AttrOrderedDict()
        self.reactions = AttrOrderedDict()
        self.compartments = AttrOrderedDict()
        self.metadata = OrderedDict()
        self._m_r_lookup = None
        self._reg_lookup = None
        self._s_matrix = None
        self._parser = None
        self._needs_update = False

    def copy(self):
        return deepcopy(self)

    def update(self):
        self._m_r_lookup = None
        self._reg_lookup = None
        self._s_matrix = None
        self._needs_update = False

    def add_compartment(self, compartment, replace=True):
        """ Add a compartment to the model.

        Arguments:
            compartment (Compartment): compartment to add
            replace (bool): replace previous compartment with same id (default: True)
        """
        if compartment.id in self.compartments and not replace:
            raise RuntimeError(f"Compartment {compartment.id} already exists.")
        self.compartments[compartment.id] = compartment

    def add_metabolite(self, metabolite, replace=True):
        """ Add a metabolite to the model.

        Arguments:
            metabolite (Metabolite): metabolite to add
            replace (bool): replace previous metabolite with same id (default: True)
        """

        if metabolite.id in self.metabolites and not replace:
            raise RuntimeError(f"Metabolite {metabolite.id} already exists.")

        if metabolite.compartment not in self.compartments:
            raise RuntimeError(f"Metabolite {metabolite.id} has invalid compartment {metabolite.compartment}.")

        self.metabolites[metabolite.id] = metabolite
        self._needs_update = True

    def add_reaction(self, reaction, replace=True):
        """ Add a reaction to the model.

        Arguments:
            reaction (Reaction): reaction to add
            replace (bool): replace previous reaction with same id (default: True)
        """
        if reaction.id in self.reactions and not replace:
            raise RuntimeError(f"Reaction {reaction.id} already exists.")
        self.reactions[reaction.id] = reaction
        self._needs_update = True

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

        r_id, reversible, stoichiometry = self._parser.parse_reaction(reaction_str)

        for m_id in stoichiometry:
            if m_id not in self.metabolites:
                self.add_metabolite(Metabolite(m_id, m_id, compartment=compartment))

        reaction = Reaction(r_id, r_id, reversible, stoichiometry)
        self.add_reaction(reaction)
        self._needs_update = True

        return r_id

    def get_reactions_by_type(self, reaction_type):
        return [rxn.id for rxn in self.reactions.values() if rxn.reaction_type == reaction_type]

    def get_exchange_reactions(self):
        return self.get_reactions_by_type(ReactionType.EXCHANGE)

    def get_compartment_metabolites(self, c_id):
        if c_id not in self.compartments.keys():
            raise RuntimeError(f"No such compartment: {c_id}")

        return [m_id for m_id, met in self.metabolites.items() if met.compartment == c_id]

    def get_external_metabolites(self, from_reactions=False):
        # TODO: a unit test should assert that result is the same from reactions and from compartments

        if from_reactions:
            external = [m_id for r_id in self.get_exchange_reactions()
                        for m_id in self.reactions[r_id].stoichiometry]
        else:
            external = [m_id for m_id, met in self.metabolites.items()
                        if self.compartments[met.compartment].external]
        return external

    def get_reaction_compartments(self, r_id):
        return {self.metabolites[m_id].compartment for m_id in self.reactions[r_id].stoichiometry}

    def get_metabolite_producers(self, m_id, reversible=False):
        """ Return the list of reactions producing a given metabolite

        Arguments:
            m_id (str): metabolite id
            reversible (bool): also include reversible consumers

        Returns:
            list: producing reactions
        """
        table = self.metabolite_reaction_lookup()

        producers = []
        for r_id, coeff in table[m_id].items():
            if coeff > 0 or reversible and self.reactions[r_id].reversible:
                producers.append(r_id)

        return producers

    def get_metabolite_consumers(self, m_id, reversible=False):
        """ Return the list of reactions consuming a given metabolite

        Arguments:
            m_id (str): metabolite id
            reversible (bool): also include reversible producers

        Returns:
            list: consuming reactions
        """
        table = self.metabolite_reaction_lookup()

        consumers = []
        for r_id, coeff in table[m_id].items():
            if coeff < 0 or reversible and self.reactions[r_id].reversible:
                consumers.append(r_id)

        return consumers

    def get_metabolite_reactions(self, m_id):
        """ Return the list of reactions associated with a given metabolite

        Arguments:
            m_id (str): metabolite id

        Returns:
            list: associated reactions
        """
        table = self.metabolite_reaction_lookup()

        return list(table[m_id].keys())

    def get_activation_targets(self, m_id):
        table = self.regulatory_lookup()
        return [r_id for r_id, kind in table[m_id].items() if kind == RegulatorType.ACTIVATOR]

    def get_inhibition_targets(self, m_id):
        table = self.regulatory_lookup()
        return [r_id for r_id, kind in table[m_id].items() if kind == RegulatorType.INHIBITOR]

    def remove_compartment(self, c_id):
        """ Remove a compartment from the model.

        Arguments:
            c_id (str): compartment id
        """
        self.remove_compartments([c_id])

    def remove_compartments(self, c_ids):
        """ Remove a compartment from the model.

        Arguments:
            c_ids (list): compartment ids
        """

        for c_id in c_ids:
            if c_id in self.compartments:
                del self.compartments[c_id]
            else:
                warn(f"No such compartment {c_id}")

        metabolites = [m_id for m_id, met in self.metabolites.items() if met.compartment in c_ids]
        self.remove_metabolites(metabolites)

    def remove_metabolite(self, m_id):
        """ Remove a metabolite from the model.

        Arguments:
            m_id (str): metabolite id
        """
        self.remove_metabolites([m_id])

    def remove_metabolites(self, id_list, safe_delete=True):
        """ Remove a list of metabolites from the model.

        Arguments:
            id_list (list): metabolite ids
            safe_delete (bool): also remove metabolites from reactions (default: True)
        """

        if safe_delete:
            m_r_lookup = self.metabolite_reaction_lookup()
            reactions = set()

        for m_id in list(id_list):
            if m_id in self.metabolites:
                del self.metabolites[m_id]
            else:
                warn(f"No such metabolite {m_id}")

            if safe_delete:
                for r_id in m_r_lookup[m_id]:
                    del self.reactions[r_id].stoichiometry[m_id]
                    reactions.add(r_id)

        if safe_delete:
            to_delete = [r_id for r_id in reactions if len(self.reactions[r_id].stoichiometry) == 0]
            self.remove_reactions(to_delete)

        self._needs_update = True

    def remove_reaction(self, r_id):
        """ Remove a reaction from the model.

        Arguments:
            r_id (str): reaction id
        """
        self.remove_reactions([r_id])

    def remove_reactions(self, id_list):
        """ Remove a list of reactions from the model.

        Arguments:
            id_list (list of str): reaction ids
        """
        for r_id in id_list:
            if r_id in self.reactions:
                del self.reactions[r_id]
            else:
                warn(f"No such reaction {r_id}")
        self._needs_update = True

    def metabolite_reaction_lookup(self):
        if not self._m_r_lookup or self._needs_update:
            self._m_r_lookup = {m_id: {} for m_id in self.metabolites}

            for r_id, reaction in self.reactions.items():
                for m_id, coeff in reaction.stoichiometry.items():
                    self._m_r_lookup[m_id][r_id] = coeff

        return self._m_r_lookup

    def regulatory_lookup(self):
        if not self._reg_lookup or self._needs_update:
            self._reg_lookup = {m_id: {} for m_id in self.metabolites}

            for r_id, reaction in self.reactions.items():
                for m_id, kind in reaction.regulators.items():
                    self._reg_lookup[m_id][r_id] = kind

        return self._reg_lookup

    def stoichiometric_matrix(self):
        """ Return a stoichiometric matrix (as a nested list)

        Returns:
            list: stoichiometric matrix
        """

        if not self._s_matrix or self._needs_update:
            self._s_matrix = [[reaction.stoichiometry[m_id] if m_id in reaction.stoichiometry else 0
                               for reaction in self.reactions.values()]
                              for m_id in self.metabolites]

        return self._s_matrix

    def print_reaction(self, r_id, use_names=False):
        """ Print a reaction to a text based representation.

        Arguments:
            r_id (str): reaction id
            use_names (bool): print metabolite names instead of ids (default: False)

        Returns:
            str: reaction string
        """

        if use_names:
            metabolite_names = {m_id: met.name for m_id, met in self.metabolites.items()}
        else:
            metabolite_names = None

        return self.reactions[r_id].to_string(metabolite_names)

    def to_string(self, use_names=False):
        """ Print the model to a text based representation.

        Arguments:
            use_names (bool): print metabolite names instead of ids (default: False)

        Returns:
            str: model as a string
        """

        return '\n'.join(self.print_reaction(r_id, use_names) for r_id in self.reactions)

    def __str__(self):
        return self.to_string()