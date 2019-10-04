from enum import Enum
import re


class Status(Enum):
    """ Enumeration of possible solution status. """
    OPTIMAL = 'Optimal'
    UNKNOWN = 'Unknown'
    SUBOPTIMAL = 'Suboptimal'
    UNBOUNDED = 'Unbounded'
    INFEASIBLE = 'Infeasible'
    INF_OR_UNB = 'Infeasible or Unbounded'


def print_values(value_dict, pattern=None, sort=False, abstol=1e-9):

    values = [(key, value) for key, value in value_dict.items() if abs(value) > abstol]

    if pattern:
        re_expr = re.compile(pattern)
        values = [x for x in values if re_expr.match(x[0]) is not None]

    if sort:
        values.sort(key=lambda x: x[1])

    entries = (f'{r_id:<12} {val: .6g}'for (r_id, val) in values)

    print('\n'.join(entries))


def print_balance(values, m_id, model, sort=False, percentage=False, equations=False, abstol=1e-9):
    inputs = model.get_metabolite_producers(m_id)
    outputs = model.get_metabolite_consumers(m_id)

    fwd_in = [(r_id, model.reactions[r_id].stoichiometry[m_id] * values[r_id], '--> o')
              for r_id in inputs if values[r_id] > 0]
    rev_in = [(r_id, model.reactions[r_id].stoichiometry[m_id] * values[r_id], 'o <--')
              for r_id in outputs if values[r_id] < 0]
    fwd_out = [(r_id, model.reactions[r_id].stoichiometry[m_id] * values[r_id], 'o -->')
               for r_id in outputs if values[r_id] > 0]
    rev_out = [(r_id, model.reactions[r_id].stoichiometry[m_id] * values[r_id], '<-- o')
                for r_id in inputs if values[r_id] < 0]

    flux_in = [x for x in fwd_in + rev_in if x[1] > abstol]
    flux_out = [x for x in fwd_out + rev_out if -x[1] > abstol]

    if sort:
        flux_in.sort(key=lambda x: x[1], reverse=True)
        flux_out.sort(key=lambda x: x[1], reverse=False)

    if percentage:
        turnover = sum([x[1] for x in flux_in])
        flux_in = [(x[0], x[1] / turnover, x[2]) for x in flux_in]
        flux_out = [(x[0], x[1] / turnover, x[2]) for x in flux_out]
        print_format = '[ {} ] {:<12} {:< 10.2%}'
    else:
        print_format = '[ {} ] {:<12} {:< 10.6g}'

    if equations:
        print_format += '\t{}'
        lines = (print_format.format(x[2], x[0], x[1], model.reactions[x[0]].to_equation())
                 for x in flux_in + flux_out)
    else:
        lines = (print_format.format(x[2], x[0], x[1]) for x in flux_in + flux_out)

    print('\n'.join(lines))


class Solution(object):
    """ Stores the results of an optimization.

    Instantiate without arguments to create an empty Solution representing a failed optimization.
    """

    def __init__(self, status=Status.UNKNOWN, message=None, fobj=None, values=None, shadow_prices=None, reduced_costs=None):
        self.status = status
        self.message = message
        self.fobj = fobj
        self.values = values
        self.shadow_prices = shadow_prices
        self.reduced_costs = reduced_costs

    def __str__(self):
        return f"Objective: {self.fobj}\nStatus: {self.status.value}\n"

    def show_values(self, pattern=None, sort=False, abstol=1e-9):
        """ Show solution results.

        Arguments:
            pattern (str): show only reactions that contain pattern (optional)
            sort (bool): sort values by magnitude (default: False)
            abstol (float): abstolute tolerance to hide null values (default: 1e-9)

        Returns:
            str: printed table with variable values
        """

        if self.values:
            print_values(self.values, pattern=pattern, sort=sort, abstol=abstol)

    def show_shadow_prices(self, pattern=None, sort=False, abstol=1e-9):
        """ Show shadow prices.

        Arguments:
            pattern (str): show only metabolites that contain pattern (optional)
            sort (bool): sort values by magnitude (default: False)
            abstol (float): abstolute tolerance to hide null values (default: 1e-9)

        Returns:
            str: printed table with shadow prices
        """

        if self.shadow_prices:
            print_values(self.shadow_prices, pattern=pattern, sort=sort, abstol=abstol)

    def show_reduced_costs(self, pattern=None, sort=False, abstol=1e-9):
        """ Show reduced costs.

        Arguments:
            pattern (str): show only reactions that contain pattern (optional)
            sort (bool): sort values by magnitude (default: False)
            abstol (float): abstolute tolerance to hide null values (default: 1e-9)

        Returns:
            str: printed table with shadow prices
        """

        if self.reduced_costs:
            print_values(self.reduced_costs, pattern=pattern, sort=sort, abstol=abstol)

    def show_metabolite_balance(self, m_id, model, sort=False, percentage=False, equations=False, abstol=1e-9):
        """ Show metabolite balance details.

        Arguments:
            m_id (str): metabolite id
            model (CBModel): model that generated the solution
            sort (bool): sort values by magnitude (default: False)
            percentage (bool): show percentage of total turnover instead of flux (default: False)
            equations (bool): show reaction equations (default: False)
            abstol (float): abstolute tolerance to hide null values (default: 1e-9)

        Returns:
            str: formatted output
        """

        if self.values:
            print_balance(self.values, m_id, model, sort=sort, percentage=percentage, equations=equations, abstol=abstol)

    def get_metabolites_turnover(self, model):
        """ Calculate metabolite turnover.

        Arguments:
            model (CBModel): model that generated the solution

        Returns:
            dict: metabolite turnover rates
        """

        if not self.values:
            return None

        m_r_table = model.metabolite_reaction_lookup()
        t = {m_id: 0.5*sum([abs(coeff * self.values[r_id]) for r_id, coeff in neighbours.items()])
             for m_id, neighbours in m_r_table.items()}
        return t

    def show_metabolite_turnover(self, model, pattern=None, sort=False, abstol=1e-9):
        """ Show solution results.

        Arguments:
            model (CBModel): model that generated the solution
            pattern (str): show only reactions that contain pattern (optional)
            sort (bool): sort values by magnitude (default: False)
            abstol (float): abstolute tolerance to hide null values (default: 1e-9)

        Returns:
            str: printed table
        """

        if self.values:
            turnover = self.get_metabolites_turnover(model)
            print_values(turnover, pattern=pattern, sort=sort, abstol=abstol)

