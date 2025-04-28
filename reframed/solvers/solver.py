from enum import Enum
from math import inf
from .solution import Solution, Status


class VarType(Enum):
    """ Enumeration of possible variable types. """
    BINARY = 'binary'
    INTEGER = 'integer'
    CONTINUOUS = 'continuous'


class Parameter(Enum):
    """ Enumeration of parameters common to all solvers. """
    TIME_LIMIT = 0
    FEASIBILITY_TOL = 1
    INT_FEASIBILITY_TOL = 2
    OPTIMALITY_TOL = 3
    MIP_REL_GAP = 4
    MIP_ABS_GAP = 5
    POOL_SIZE = 6
    POOL_GAP = 7


class Solver(object):
    """ Abstract class representing a generic solver.

    All solver interfaces should implement the methods defined in this class.
    """

    def __init__(self, model=None):
        self.problem = None
        self.model = model
        self.variables = []
        self.constraints = []
        self.objective = {}
        self.minimize = True
        self._cached_vars = {}
        self._cached_constrs = {}

    def add_variable(self, var_id, lb=-inf, ub=inf, vartype=VarType.CONTINUOUS):
        """ Add a variable to the current problem.

        Arguments:
            var_id (str): variable identifier
            lb (float): lower bound
            ub (float): upper bound
            vartype (VarType): variable type (default: CONTINUOUS)
        """

        self._cached_vars[var_id] = (lb, ub, vartype)

    def add_variables(self, var_dict):
        """ Solver specific implementation """
        pass

    def add_constraint(self, constr_id, lhs, sense='=', rhs=0):
        """ Add a constraint to the current problem.

        Arguments:
            constr_id (str): constraint identifier
            lhs (dict): variables and respective coefficients
            sense (str): constraint sense (any of: '<', '=', '>'; default '=')
            rhs (float): right-hand side of equation (default: 0)
        """
        
        self._cached_constrs[constr_id] = (lhs, sense, rhs)

    def add_constraints(self, constr_dict):
        """ Solver specific implementation """
        pass

    def remove_constraint(self, constr_id):
        """ Solver specific implementation """
        pass

    def update(self):
        """ Update internal structure. Used for efficient lazy updating. """
        
        if len(self._cached_vars) > 0:
            self.add_variables(self._cached_vars)
            self._cached_vars = {}

        if len(self._cached_constrs) > 0: 
            self.add_constraints(self._cached_constrs)
            self._cached_constrs = {}


    def build_problem(self, model):
        """ Create problem structure for a given model.

        Arguments:
            model : CBModel
        """

        for r_id, reaction in model.reactions.items():
            self.add_variable(r_id, reaction.lb, reaction.ub)

        table = model.metabolite_reaction_lookup()

        for m_id in model.metabolites:
            self.add_constraint(m_id, table[m_id])

        self.update()

    def solve(self, objective=None, minimize=True, model=None, constraints=None, get_values=True, shadow_prices=False, reduced_costs=False,
              allow_suboptimal=False):
        """ Solve the optimization problem.

        Arguments:
            objective (dict): linear objective (optional)
            minimize (bool): solve a minimization problem (default: True)
            model (CBModel): model (optional, leave blank to reuse previous model structure)
            constraints (dict): additional constraints (optional)
            get_values (bool or list): set to false for speedup if you only care about the objective value (default: True)
            shadow_prices (bool): return shadow prices if available (default: False)
            reduced_costs (bool): return reduced costs if available (default: False)
            allow_suboptimal (bool): return suboptimal solution if found (default: False)

        Returns:
            Solution: solution
        """

        if model:
            self.build_problem(model)
        else:
            self.update()

        if objective is not None:
            self.set_objective(objective, minimize)

        if constraints:
            old_bounds = self.set_temporary_bounds(constraints)

        status = self.internal_solve()

        if status == Status.OPTIMAL or status == Status.SUBOPTIMAL and allow_suboptimal:
            solution = self.get_solution(status, get_values, shadow_prices, reduced_costs)
        else:
            solution = Solution(status)

        if constraints:
            self.reset_bounds(old_bounds)

        return solution

    def set_objective(self, objective, minimize=True):
        """ Set a predefined objective for this problem.

        Args:
            linear (str or dict): linear coefficients (or a single variable to optimize)
            minimize (bool): solve a minimization problem (default: True)

        """
        pass


    def set_temporary_bounds(self, bounds):
        pass

    def reset_bounds(self, bounds):
        pass

    def internal_solve(self):
        pass

    def get_solution(self, status, get_values=True, shadow_prices=False, reduced_costs=False):
        pass

    def set_parameter(self, parameter, value):
        """ Set a parameter value for this optimization problem

        Arguments:
            parameter (Parameter): parameter type
            value (float): parameter value
        """

        raise Exception('Not implemented for this solver.')


    def set_logging(self, enabled=False):
        """ Enable or disable log output:

        Arguments:
            enabled (bool): turn logging on (default: False)
        """

        raise Exception('Not implemented for this solver.')

    def write_to_file(self, filename):
        """ Write problem to file:

        Arguments:
            filename (str): file path
        """

        raise Exception('Not implemented for this solver.')



