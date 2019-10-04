from ..solvers.solution import print_values, print_balance


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

        self.abundance = {org_id: solution.values[f"x_{org_id}"]
                          for org_id in self.community.organisms}

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
                    entries.append(('--> o', r_id, rate))
                elif rate < -abstol:
                    entries.append(('o -->', r_id, rate))

            for org_id in self.community.organisms:

                if (org_id, m_id) not in self.exchange_map:
                    continue

                rate = self.exchange_map[(org_id, m_id)]
                if rate > abstol:
                    entries.append(('--> o', org_id, rate))
                elif rate < -abstol:
                    entries.append(('o -->', org_id, rate))

            if entries:
                print(m_id)
                entries.sort(key=lambda x: -abs(x[2]))

                for (sense, org_id, rate) in entries:
                    print(f'[ {sense} ] {org_id:<12} {rate:< 10.6g}')