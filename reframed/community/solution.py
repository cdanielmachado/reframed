from ..solvers.solution import print_values, print_balance
from ..core.elements import molecular_weight
import pandas as pd


class CommunitySolution(object):

    def __init__(self, community, solution):
        self.community = community
        self.solution = solution
        self.growth = None
        self.abundance = None
        self.exchange = None
        self.internal = None
        self.normalized = None
        self._exchange_map = None
        self.parse_solution()

    def __str__(self):
        growth = f"Community growth: {self.growth}\n"
        abundance = "\n".join(f"{org_id}\t{val}" for org_id, val in self.abundance.items())
        return growth + abundance

    def __repr__(self):
        return str(self)

    @property
    def exchange_map(self):
        if self._exchange_map is None:
            self._exchange_map = self.compute_exchanges()

        return self._exchange_map

    def parse_solution(self):
        model = self.community.merged_model
        solution = self.solution
        reaction_map = self.community.reaction_map

        self.growth = solution.values[model.biomass_reaction]

        self.abundance = {}
        for org_id, organism in self.community.organisms.items():
            growth_i = self.community.reaction_map[(org_id, organism.biomass_reaction)]
            self.abundance[org_id] = solution.values[growth_i] / self.growth

        self.exchange = {r_id: solution.values[r_id]
                         for r_id in model.get_exchange_reactions()}

        self.internal = {}
        self.normalized = {}

        for org_id, organism in self.community.organisms.items():

            abundance = self.abundance[org_id]

            fluxes = {r_id: solution.values[reaction_map[(org_id, r_id)]]
                      for r_id in organism.reactions
                      if (org_id, r_id) in reaction_map}

            rates = {r_id: fluxes[r_id] / abundance if abundance > 0 else 0
                     for r_id in organism.reactions if r_id in fluxes}

            self.internal[org_id] = fluxes
            self.normalized[org_id] = rates

    # calculate overall exchanges (organism x metabolite) -> rate

    def compute_exchanges(self):
        model = self.community.merged_model
        reaction_map = self.community.reaction_map
        exchanges = {}

        for m_id in model.get_external_metabolites():

            for org_id, organism in self.community.organisms.items():
                rate = 0
                if m_id not in organism.metabolites:
                    continue

                for r_id in organism.get_metabolite_reactions(m_id):
                    if (org_id, r_id) not in reaction_map:
                        continue

                    flux = self.solution.values[reaction_map[(org_id, r_id)]]

                    if flux != 0:
                        coeff = organism.reactions[r_id].stoichiometry[m_id]
                        rate += coeff*flux

                if rate != 0:
                    exchanges[(org_id, m_id)] = rate

        return exchanges

    def cross_feeding(self, as_df=False, abstol=1e-6):
        exchanges = self.compute_exchanges()
        cross_all = []

        for m_id in self.community.merged_model.get_external_metabolites():
            r_out = {x: r for (x, m), r in exchanges.items() if m == m_id and r > abstol}
            r_in = {x: -r for (x, m), r in exchanges.items() if m == m_id and -r > abstol}

            total_in = sum(r_in.values())
            total_out = sum(r_out.values())
            total = max(total_in, total_out)

            if total_in > total_out:
                r_out[None] = total_in - total_out
            if total_out > total_in:
                r_in[None] = total_out - total_in

            cross = [(o1, o2, m_id, r1 * r2 / total) for o1, r1 in r_out.items() for o2, r2 in r_in.items()]
            cross_all.extend(cross)

        if as_df:
            cross_all = pd.DataFrame(cross_all, columns=["donor", "receiver", "compound", "rate"])

        return cross_all

    def mass_flow(self, as_df=False, abstol=1e-6):

        def get_mass(x):
            met = self.community.merged_model.metabolites[x]
            formula = met.metadata.get('FORMULA', '')
            mw = molecular_weight(formula)
            return 0.001 * mw

        entities = list(self.community.organisms) + [None]
        flow = {(o1, o2): 0 for o1 in entities for o2 in entities}

        for o1, o2, m_id, rate in self.cross_feeding():
            flow[(o1, o2)] += get_mass(m_id) * rate

        flow = {key: val for key, val in flow.items() if val > abstol}

        if as_df:
            flow = [(o1, o2, val) for (o1, o2), val in flow.items()]
            flow = pd.DataFrame(flow, columns=["donor", "receiver", "flow"])

        return flow

    def print_external_fluxes(self, pattern=None, sort=False, abstol=1e-9):
        print_values(self.exchange, pattern=pattern, sort=sort, abstol=abstol)

    def print_internal_fluxes(self, org_id, normalized=False, pattern=None, sort=False, abstol=1e-9):
        if normalized:
            print_values(self.normalized[org_id], pattern=pattern, sort=sort, abstol=abstol)
        else:
            print_values(self.internal[org_id], pattern=pattern, sort=sort, abstol=abstol)

    def print_external_balance(self, m_id, sort=False, percentage=False, equations=False, abstol=1e-9):

        print_balance(self.solution.values, m_id, self.community.merged_model, sort=sort, percentage=percentage, equations=equations,
                      abstol=abstol)

    def print_exchanges(self, m_id=None, abstol=1e-9):

        model = self.community.merged_model
        exchange_rxns = set(model.get_exchange_reactions())

        if m_id:
            mets = [m_id]
        else:
            mets = model.get_external_metabolites()

        for m_id in mets:

            entries = []

            ext_rxns = set(model.get_metabolite_reactions(m_id)) & exchange_rxns

            for r_id in ext_rxns:

                flux = self.solution.values[r_id]
                coeff = model.reactions[r_id].stoichiometry[m_id]
                rate = coeff*flux

                if rate > abstol:
                    entries.append(('=> *   ', "in", rate))
                elif rate < -abstol:
                    entries.append(('   * =>', "out", rate))

            for org_id in self.community.organisms:

                if (org_id, m_id) not in self.exchange_map:
                    continue

                rate = self.exchange_map[(org_id, m_id)]
                if rate > abstol:
                    entries.append(('O --> *', org_id, rate))
                elif rate < -abstol:
                    entries.append(('* --> O', org_id, rate))

            if entries:
                print(m_id)
                entries.sort(key=lambda x: x[2])

                for (sense, org_id, rate) in entries:
                    print(f'[ {sense} ] {org_id:<12} {rate:< 10.6g}')