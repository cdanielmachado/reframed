from optlang import Model, Variable, Constraint, Objective
from optlang.symbolics import Zero, add
from .solver import Solver, VarType, Parameter, default_parameters
from .solution import Solution, Status
from math import inf
from warnings import warn
from collections.abc import Iterable


status_mapping = {
    'optimal': Status.OPTIMAL,
    'unbounded': Status.UNBOUNDED,
    'infeasible': Status.INFEASIBLE,
    'infeasible_or_unbounded': Status.INF_OR_UNB,
    'suboptimal': Status.SUBOPTIMAL,
}


class OptLangSolver(Solver):
    """ Implements the gurobi solver interface. """

    def __init__(self, model=None):
        Solver.__init__(self)
        self.problem = Model()

        self.parameter_mapping = {
            Parameter.TIME_LIMIT: self.problem.configuration.timeout,
            Parameter.FEASIBILITY_TOL: self.problem.configuration.tolerances.feasibility,
            Parameter.OPTIMALITY_TOL: 1e-9,
            Parameter.INT_FEASIBILITY_TOL: self.problem.configuration.tolerances.integrality,
        }

        self.set_parameters(default_parameters)
        self.set_logging(False)

        if model:
            self.build_problem(model)

    def add_variable(self, var_id, lb=-inf, ub=inf, vartype='continuous', update=True):
        """ Add a variable to the current problem.

        Arguments:
            var_id (str): variable identifier
            lb (float): lower bound
            ub (float): upper bound
            vartype (VarType): variable type (default: CONTINUOUS)
            update (bool): update problem immediately (default: True)
        """

        if var_id in self.var_ids:
            var = self.problem.variables[var_id]
            var.lb = lb
            var.ub = ub
            var.type = vartype
        else:
            var = Variable(var_id, lb=lb, ub=ub, type=vartype)
            self.problem.add(var)
            self.var_ids.append(var_id)

        if update:
            self.problem.update()

    def add_constraint(self, constr_id, lhs, sense='=', rhs=0, update=True):
        """ Add a constraint to the current problem.

        Arguments:
            constr_id (str): constraint identifier
            lhs (dict): variables and respective coefficients
            sense (str): constraint sense (any of: '<', '=', '>'; default '=')
            rhs (float): right-hand side of equation (default: 0)
            update (bool): update problem immediately (default: True)
        """

        if constr_id in self.constr_ids:
            self.problem.remove(constr_id)

        if sense == '=':
            constr = Constraint(Zero, lb=rhs, ub=rhs, name=constr_id)
        elif sense == '>':
            constr = Constraint(Zero, lb=rhs, name=constr_id)
        elif sense == '<':
            constr = Constraint(Zero, ub=rhs, name=constr_id)
        else:
            raise RuntimeError(f"Invalid constraint direction: {sense}")

        self.problem.add(constr)
        self.constr_ids.append(constr_id)

        expr = {self.problem.variables[r_id]: coeff for r_id, coeff in lhs.items() if coeff}
        self.problem.constraints[constr_id].set_linear_coefficients(expr)

        if update:
            self.problem.update()

    def remove_variable(self, var_id):
        """ Remove a variable from the current problem.

        Arguments:
            var_id (str): variable identifier
        """
        self.remove_variables([var_id])

    def remove_variables(self, var_ids):
        """ Remove variables from the current problem.

        Arguments:
            var_ids (list): variable identifiers
        """

        for var_id in var_ids:
            if var_id in self.var_ids:
                self.problem.remove(var_id)
                self.var_ids.remove(var_id)

    def remove_constraint(self, constr_id):
        """ Remove a constraint from the current problem.

        Arguments:
            constr_id (str): constraint identifier
        """
        self.remove_constraints([constr_id])

    def remove_constraints(self, constr_ids):
        """ Remove constraints from the current problem.

        Arguments:
            constr_ids (list): constraint identifiers
        """

        for constr_id in constr_ids:
            if constr_id in self.constr_ids:
                self.problem.remove(constr_id)
                self.constr_ids.remove(constr_id)

    def set_objective(self, linear=None, quadratic=None, minimize=True):
        """ Set a predefined objective for this problem.

        Args:
            linear (dict): linear coefficients (optional)
            quadratic (dict): quadratic coefficients (optional)
            minimize (bool): solve a minimization problem (default: True)

        Notes:
            Setting the objective is optional. It can also be passed directly when calling **solve**.

        """

        if linear is None:
            linear = {}

        if quadratic is None:
            quadratic = {}

        if linear and not quadratic:
            objective = {}

            if isinstance(linear, str):
                objective = {self.problem.variables[linear]: 1}
                if linear not in self.var_ids:
                    warn(f"Objective variable not previously declared: {linear}")
            else:
                for r_id, val in linear.items():
                    if r_id not in self.var_ids:
                        warn(f"Objective variable not previously declared: {r_id}")
                    elif val != 0:
                        objective[self.problem.variables[r_id]] = val

            self.problem.objective = Objective(Zero, direction=('min' if minimize else 'max'), sloppy=True)
            self.problem.objective.set_linear_coefficients(objective)
        else:
            objective = []

            for r_id, val in linear.items():
                if r_id not in self.var_ids:
                    warn(f"Objective variable not previously declared: {r_id}")
                elif val != 0:
                    objective.append(val * self.problem.variables[r_id])

            for (r_id1, r_id2), val in quadratic.items():
                if r_id1 not in self.var_ids:
                    warn(f"Objective variable not previously declared: {r_id1}")
                elif r_id2 not in self.var_ids:
                    warn(f"Objective variable not previously declared: {r_id2}")
                elif val != 0:
                    objective.append(val * self.problem.variables[r_id1] * self.problem.variables[r_id2])

            objective_expr = add(objective)
            self.problem.objective = Objective(objective_expr, direction=('min' if minimize else 'max'), sloppy=True)

    def solve(self, linear=None, quadratic=None, minimize=None, model=None, constraints=None, get_values=True,
              shadow_prices=False, reduced_costs=False, pool_size=0, pool_gap=None):
        """ Solve the optimization problem.

        Arguments:
            linear (str or dict): linear coefficients (or a single variable to optimize)
            quadratic (dict): quadratic objective (optional)
            minimize (bool): solve a minimization problem (default: True)
            model (CBModel): model (optional, leave blank to reuse previous model structure)
            constraints (dict): additional constraints (optional)
            get_values (bool or list): set to false for speedup if you only care about the objective (default: True)
            shadow_prices (bool): return shadow prices if available (default: False)
            reduced_costs (bool): return reduced costs if available (default: False)
            pool_size (int): calculate solution pool of given size (only for MILP problems)
            pool_gap (float): maximum relative gap for solutions in pool (optional)

        Returns:
            Solution: solution
        """

        if model:
            self.build_problem(model)

        problem = self.problem

        if constraints:
            old_constraints = {}
            for r_id, x in constraints.items():
                lb, ub = x if isinstance(x, tuple) else (x, x)
                if r_id in self.var_ids:
                    lpvar = problem.variables[r_id]
                    old_constraints[r_id] = (lpvar.lb, lpvar.ub)
                    lpvar.lb, lpvar.ub = lb, ub
                else:
                    warn(f"Constrained variable '{r_id}' not previously declared")
            problem.update()

        self.set_objective(linear, quadratic, minimize)

        # run the optimization
        if pool_size > 1:
            raise RuntimeError("OptLang interface does not support solution pools.")

        problem.optimize()

        status = status_mapping.get(problem.status, Status.UNKNOWN)
        message = str(problem.status)

        if status == Status.OPTIMAL:
            fobj = problem.objective.value
            values, s_prices, r_costs = None, None, None

            if get_values:
                values = dict(problem.primal_values)

                if isinstance(get_values, Iterable):
                    values = {x: values[x] for x in get_values}

            if shadow_prices:
                s_prices = dict(problem.shadow_prices)

            if reduced_costs:
                r_costs = dict(problem.reduced_costs)

            solution = Solution(status, message, fobj, values, s_prices, r_costs)
        else:
            solution = Solution(status, message)

        # restore values of temporary constraints
        if constraints:
            for r_id, (lb, ub) in old_constraints.items():
                lpvar = problem.variables[r_id]
                lpvar.lb, lpvar.ub = lb, ub
            problem.update()

        return solution

    def set_parameter(self, parameter, value):
        """ Set a parameter value for this optimization problem

        Arguments:
            parameter (Parameter): parameter type
            value (float): parameter value
        """

        if parameter in self.parameter_mapping:
            self.parameter_mapping[parameter] = value
        else:
            raise RuntimeError('Parameter unknown (or not yet supported).')

    def set_logging(self, enabled=False):
        """ Enable or disable log output:

        Arguments:
            enabled (bool): turn logging on (default: False)
        """

        self.problem.configuration.verbosity = 3 if enabled else 0

    def write_to_file(self, filename):
        """ Write problem to file:

        Arguments:
            filename (str): file path
        """

        with open(filename, "w") as f:
            f.write(self.problem.to_lp())

