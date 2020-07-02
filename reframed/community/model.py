from ..core.model import AttrOrderedDict, ReactionType, Compartment, Metabolite
from ..core.cbmodel import CBModel, CBReaction, Gene, Protein, GPRAssociation
from math import inf
from warnings import warn


class Community(object):
    def __init__(self, community_id, models, copy_models=False):
        self.id = community_id
        self.organisms = AttrOrderedDict()
        self._merged_model = None
        self.reaction_map = None
        self.metabolite_map = None

        model_ids = {model.id for model in models}

        if len(model_ids) < len(models):
            warn("Model ids are not unique, repeated models will be discarded.")

        for model in models:
            self.organisms[model.id] = model.copy() if copy_models else model

    def size(self):
        return len(self.organisms)

    @property
    def merged_model(self):
        if self._merged_model is None:
            self._merged_model = self.merge_models()

        return self._merged_model

    def merge_models(self):
        comm_model = CBModel(self.id)
        old_ext_comps = []
        ext_mets = []
        self.reaction_map = {}
        self.metabolite_map = {}

        # default IDs
        ext_comp_id = "ext"
        biomass_id = "community_biomass"
        comm_growth = "community_growth"

        # create external compartment

        comp = Compartment(ext_comp_id, "extracellular environment", external=True)
        comm_model.add_compartment(comp)

        # community biomass

        met = Metabolite(biomass_id, "Total community biomass", ext_comp_id)
        comm_model.add_metabolite(met)

        rxn = CBReaction(comm_growth, name="Community growth rate",
                         reversible=False, stoichiometry={biomass_id: -1},
                         lb=0, ub=inf, objective=1)

        comm_model.add_reaction(rxn)

        # add each organism

        for org_id, model in self.organisms.items():

            def rename(old_id):
                return f"{old_id}_{org_id}"

            # add internal compartments

            for c_id, comp in model.compartments.items():
                if comp.external:
                    old_ext_comps.append(c_id)
                else:
                    new_comp = Compartment(rename(c_id), comp.name)
                    comm_model.add_compartment(new_comp)

            # add metabolites

            for m_id, met in model.metabolites.items():
                if met.compartment not in old_ext_comps:  # if is internal
                    new_id = rename(m_id)
                    new_met = Metabolite(new_id, met.name, rename(met.compartment))
                    new_met.metadata = met.metadata.copy()
                    comm_model.add_metabolite(new_met)
                    self.metabolite_map[(org_id, m_id)] = new_id

                elif m_id not in comm_model.metabolites:  # if is external but was not added yet
                    new_met = Metabolite(m_id, met.name, ext_comp_id)
                    new_met.metadata = met.metadata.copy()
                    comm_model.add_metabolite(new_met)
                    ext_mets.append(new_met.id)

            # add genes

            for g_id, gene in model.genes.items():
                new_id = rename(g_id)
                new_gene = Gene(new_id, gene.name)
                new_gene.metadata = gene.metadata.copy()
                comm_model.add_gene(new_gene)

            # add internal reactions

            for r_id, rxn in model.reactions.items():

                if rxn.reaction_type == ReactionType.EXCHANGE:
                    continue

                new_id = rename(r_id)
                new_stoichiometry = {
                    m_id if m_id in ext_mets else rename(m_id): coeff
                    for m_id, coeff in rxn.stoichiometry.items()
                }

                if r_id == model.biomass_reaction:
                    new_stoichiometry[biomass_id] = 1

                if rxn.gpr is None:
                    new_gpr = None
                else:
                    new_gpr = GPRAssociation()
                    new_gpr.metadata = rxn.gpr.metadata.copy()

                    for protein in rxn.gpr.proteins:
                        new_protein = Protein()
                        new_protein.genes = [rename(g_id) for g_id in protein.genes]
                        new_protein.metadata = protein.metadata.copy()
                        new_gpr.proteins.append(new_protein)

                new_rxn = CBReaction(
                    new_id,
                    name=rxn.name,
                    reversible=rxn.reversible,
                    stoichiometry=new_stoichiometry,
                    reaction_type=rxn.reaction_type,
                    lb=rxn.lb,
                    ub=rxn.ub,
                    gpr_association=new_gpr
                )

                comm_model.add_reaction(new_rxn)
                new_rxn.metadata = rxn.metadata.copy()
                self.reaction_map[(org_id, r_id)] = new_id

        # Add exchange reactions

        for m_id in ext_mets:
            rxn = CBReaction(f"R_EX_{m_id}", reversible=True, stoichiometry={m_id: -1},
                             reaction_type=ReactionType.EXCHANGE)
            comm_model.add_reaction(rxn)

        return comm_model

