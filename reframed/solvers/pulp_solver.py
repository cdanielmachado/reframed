from math import inf
from .solver import Solver, VarType#, Parameter, default_parameters
from .solution import Solution, Status
from pulp import LpProblem, LpVariable, LpConstraint, LpAffineExpression, lpSum, value, getSolver, LpSolverDefault
from pulp.constants import *

status_mapping = {
    LpSolutionOptimal: Status.OPTIMAL,
    LpSolutionInfeasible: Status.INFEASIBLE,
    LpSolutionNoSolutionFound: Status.INF_OR_UNB,
    LpSolutionUnbounded: Status.UNBOUNDED,
    LpSolutionIntegerFeasible: Status.SUBOPTIMAL
}


vartype_mapping = {
    VarType.BINARY: LpBinary,
    VarType.INTEGER: LpInteger,
    VarType.CONTINUOUS: LpContinuous
}


class PuLPSolver(Solver):
    """ Implements the solver interface using PuLP. """

    def __init__(self, model=None, interface=None):
        Solver.__init__(self)

        self.problem = LpProblem()
        self.variables = {}
        self.constraints = {}

        if interface is None:
            self.interface = LpSolverDefault
        else:
            self.interface = getSolver(interface)

#        self.set_parameters(default_parameters)
#        self.set_logging(False)

        if model:
            self.build_problem(model)

    def add_variable(self, var_id, lb=-inf, ub=inf, vartype=VarType.CONTINUOUS, update=True):
        """ Add a variable to the current problem.

        Arguments:
            var_id (str): variable identifier
            lb (float): lower bound
            ub (float): upper bound
            vartype (VarType): variable type (default: CONTINUOUS)
            update (bool): update problem immediately (default: True)
        """

        # fix infinity
        lb = None if lb == -inf else lb
        ub = None if ub == inf else ub


        var = LpVariable(var_id, lowBound=lb, upBound=ub, cat=vartype_mapping[vartype])
        self.var_ids.append(var_id)
        self.variables[var_id] = var
        self.problem.addVariable(var)

        #TODO: implement lazy caching 

    def add_constraint(self, constr_id, lhs, sense='=', rhs=0, update=True):
        """ Add a constraint to the current problem.

        Arguments:
            constr_id (str): constraint identifier
            lhs (dict): variables and respective coefficients
            sense (str): constraint sense (any of: '<', '=', '>'; default '=')
            rhs (float): right-hand side of equation (default: 0)
            update (bool): update problem immediately (default: True)
        """

        sense_map = {
            '=': LpConstraintEQ,
            '<': LpConstraintLE,
            '>': LpConstraintGE,
        }

        expr = LpAffineExpression({self.variables[var_id]: val for var_id, val in lhs.items()})
        
        constr = LpConstraint(e=expr, sense=sense_map[sense], name=constr_id, rhs=rhs)
        
        self.constr_ids.append(constr_id)
        self.constraints[constr_id] = constr
        self.problem.addConstraint(constr)

        #TODO: implement lazy caching 

    def remove_variable(self, var_id):
        """ Remove a variable from the current problem.

        Arguments:
            var_id (str): variable identifier
        """
        self.var_ids.remove(var_id)
        del self.variables[var_id]

    def remove_variables(self, var_ids):
        """ Remove variables from the current problem.

        Arguments:
            var_ids (list): variable identifiers
        """
        for var_id in var_ids:
            self.remove_variable(var_id)

    def remove_constraint(self, constr_id):
        """ Remove a constraint from the current problem.

        Arguments:
            constr_id (str): constraint identifier
        """
        self.constr_ids.remove(constr_id)
        del self.constraints[constr_id]

    def remove_constraints(self, constr_ids):
        """ Remove constraints from the current problem.

        Arguments:
            constr_ids (list): constraint identifiers
        """
        for constr_id in constr_ids:
            self.remove_constraint(constr_id)

    def list_variables(self):
        """ Get a list of the variable ids defined for the current problem.

        Returns:
            list: variable ids
        """
        return self.var_ids

    def list_constraints(self):
        """ Get a list of the constraint ids defined for the current problem.

        Returns:
            list: constraint ids
        """
        return self.constr_ids

    def set_bounds(self, bounds_dict):
        """ Set lower and upper bounds from tuple dictionary
        Arguments:
            bounds_dict (dict): lower and upper bounds
        """

        for var_id, (lb, ub) in bounds_dict.items():
            self.variables[var_id].lowBound = lb
            self.variables[var_id].upBound = ub
 
    def update(self):
        """ Update internal structure. Used for efficient lazy updating. """
        pass

    def set_objective(self, linear=None, quadratic=None, minimize=True):
        """ Set a predefined objective for this problem.

        Args:
            linear (str or dict): linear coefficients (or a single variable to optimize)
            quadratic (dict): quadratic coefficients (optional)
            minimize (bool): solve a minimization problem (default: True)

        Notes:
            Setting the objective is optional. It can also be passed directly when calling **solve**.

        """

        if quadratic is not None:
            raise Exception('PuLP wrapper does not support quadratic objectives.')
        
        objective = lpSum([coeff * self.variables[var_id] for var_id, coeff in linear.items() if coeff != 0])
        self.problem.setObjective(objective)
        self.problem.sense = LpMinimize if minimize else LpMaximize

    def build_problem(self, model):
        """ Create problem structure for a given model.

        Arguments:
            model : CBModel
        """

        for r_id, reaction in model.reactions.items():
            self.add_variable(r_id, reaction.lb, reaction.ub, update=False)
        self.update()

        table = model.metabolite_reaction_lookup()
        for m_id in model.metabolites:
            self.add_constraint(m_id, table[m_id], update=False)
        self.update()

    def solve(self, linear=None, quadratic=None, minimize=True, model=None, constraints=None, get_values=True,
              shadow_prices=False, reduced_costs=False, pool_size=0, pool_gap=None):
        """ Solve the optimization problem.

        Arguments:
            linear (dict): linear objective (optional)
            quadratic (dict): quadratic objective (optional)
            minimize (bool): solve a minimization problem (default: True)
            model (CBModel): model (optional, leave blank to reuse previous model structure)
            constraints (dict): additional constraints (optional)
            get_values (bool or list): set to false for speedup if you only care about the objective value (default: True)
            shadow_prices (bool): return shadow prices if available (default: False)
            reduced_costs (bool): return reduced costs if available (default: False)
            pool_size (int): calculate solution pool of given size (only for MILP problems)
            pool_gap (float): maximum relative gap for solutions in pool (optional)

        Returns:
            Solution: solution
        """

        if model:
            self.build_problem(model)

        if constraints:
            changed_lb, changed_ub = self.temporary_bounds(constraints)

        self.set_objective(linear, quadratic, minimize)

        if pool_size > 1: 
            raise Exception('Solution pool not implemented for PuLP wrapper.')
        else:
            status = self.problem.solve(self.interface)

            status = status_mapping.get(status, Status.UNKNOWN)
            message = str(status)

            if status == Status.OPTIMAL:
                fobj = value(self.problem.objective)
                values, s_prices, r_costs = None, None, None

            if get_values:
                if isinstance(get_values, list):
                    values = {r_id: value(self.variables[r_id]) for r_id in get_values}
                else:
                    values = {var_id: value(var) for var_id, var in self.variables.items()}

                if shadow_prices:
                    pass #TODO: implement

                if reduced_costs:
                    pass #TODO: implement

                solution = Solution(status, message, fobj, values, s_prices, r_costs)
            else:
                solution = Solution(status, message)

        if constraints:
            self.reset_bounds(changed_lb, changed_ub)

        return solution

    def get_solution_pool(self, get_values=True):
        """ Return a solution pool for MILP problems.
        Must be called after using solve with pool_size argument > 0.

        Arguments:

            get_values (bool or list): set to false for speedup if you only care about the objective value (default: True)

        Returns:
            list: list of Solution objects

        """
        raise Exception('Not implemented for this solver.')

    def set_parameter(self, parameter, value):
        """ Set a parameter value for this optimization problem

        Arguments:
            parameter (Parameter): parameter type
            value (float): parameter value
        """

        raise Exception('Not implemented for this solver.')

    def set_parameters(self, parameters):
        """ Set values for multiple parameters

        Arguments:
            parameters (dict of Parameter to value): parameter values
        """

        for parameter, value in parameters.items():
            self.set_parameter(parameter, value)

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

        self.problem.writeLP(filename)



