from re import compile
from collections import OrderedDict
from math import inf


class ReactionParser(object):

    def __init__(self):
        id_re = '[a-zA-Z]\w*'
        pos_float_re = '\d+(?:\.\d+)?(?:e[+-]?\d+)?'
        float_re = '-?\d+(?:\.\d+)?(?:e[+-]?\d+)?'

        compound = '(?:' + pos_float_re + '\s+)?' + id_re
        expression = compound + '(?:\s*\+\s*' + compound + ')*'
        bounds = '\[\s*(?P<lb>' + float_re + ')?\s*,\s*(?P<ub>' + float_re + ')?\s*\]'
        objective = '@' + float_re
        reaction = '^(?P<reaction_id>' + id_re + ')\s*:' + \
                   '\s*(?P<substrates>' + expression + ')?' + \
                   '\s*(?P<direction>-->|<->)' + \
                   '\s*(?P<products>' + expression + ')?' + \
                   '\s*(?P<bounds>' + bounds + ')?' + \
                   '\s*(?P<objective>' + objective + ')?$'

        self.regex_compound = compile('(?P<coeff>' + pos_float_re + '\s+)?(?P<met_id>' + id_re + ')')
        self.regex_bounds = compile(bounds)
        self.regex_reaction = compile(reaction)

    def parse_reaction(self, reaction_str, kind=None):
        match = self.regex_reaction.match(reaction_str)

        if not match:
            raise SyntaxError('Unable to parse: ' + reaction_str)

        r_id = match.group('reaction_id')
        reversible = match.group('direction') == '<->'
        substrates = match.group('substrates')
        products = match.group('products')

        stoichiometry = OrderedDict()

        if substrates:
            left_coeffs = self.parse_coefficients(substrates, sense=-1)
            stoichiometry.update(left_coeffs)

        if products:
            right_coeffs = self.parse_coefficients(products, sense=1)
            for m_id, val in right_coeffs:
                if m_id in stoichiometry:
                    new_val = val + stoichiometry[m_id]
                    stoichiometry[m_id] = new_val
                else:
                    stoichiometry[m_id] = val

        if kind is None:
            return r_id, reversible, stoichiometry

        if kind == 'cb':
            bounds = match.group('bounds')
            lb, ub = self.parse_bounds(bounds, reversible)
            objective = match.group('objective')
            obj_coeff = float(objective[1:]) if objective else 0
            return r_id, reversible, stoichiometry, lb, ub, obj_coeff

    def parse_coefficients(self, expression, sense):
        coefficients = []
        terms = expression.split('+')

        for term in terms:
            match = self.regex_compound.match(term.strip())
            coeff = sense * float(match.group('coeff')) if match.group('coeff') else sense
            m_id = match.group('met_id')
            coefficients.append((m_id, coeff))

        return coefficients

    def parse_bounds(self, expression, reversible):
        lb = -inf if reversible else 0.0
        ub = inf

        if expression:
            match = self.regex_bounds.match(expression)
            if match.group('lb'):
                lb = float(match.group('lb'))
            if match.group('ub'):
                ub = float(match.group('ub'))

        return lb, ub
