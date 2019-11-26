from collections import Iterable
from .solver import Solver, VarType, Parameter, default_parameters
from .solution import Solution, Status
from gurobipy import Model as GurobiModel, GRB, quicksum
from math import inf
from warnings import warn


def infinity_fix(val):
    if val == inf:
        return GRB.INFINITY
    elif val == -inf:
        return -GRB.INFINITY
    else:
        return val


status_mapping = {
    GRB.OPTIMAL: Status.OPTIMAL,
    GRB.UNBOUNDED: Status.UNBOUNDED,
    GRB.INFEASIBLE: Status.INFEASIBLE,
    GRB.INF_OR_UNBD: Status.INF_OR_UNB
}

vartype_mapping = {
    VarType.BINARY: GRB.BINARY,
    VarType.INTEGER: GRB.INTEGER,
    VarType.CONTINUOUS: GRB.CONTINUOUS
}

parameter_mapping = {
    Parameter.TIME_LIMIT: GRB.Param.TimeLimit,
    Parameter.FEASIBILITY_TOL: GRB.Param.FeasibilityTol,
    Parameter.INT_FEASIBILITY_TOL: GRB.Param.IntFeasTol,
    Parameter.OPTIMALITY_TOL: GRB.Param.OptimalityTol,
    Parameter.MIP_ABS_GAP: GRB.Param.MIPGapAbs,
    Parameter.MIP_REL_GAP: GRB.Param.MIPGap,
    Parameter.POOL_SIZE: GRB.Param.PoolSolutions,
    Parameter.POOL_GAP: GRB.Param.PoolGap
}


class GurobiSolver(Solver):
    """ Implements the gurobi solver interface. """

    def __init__(self, model=None):
        Solver.__init__(self)
        self.problem = GurobiModel()
        self.set_logging(False)
        self.set_parameters(default_parameters)
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

        lb = infinity_fix(lb)
        ub = infinity_fix(ub)

        if var_id in self.var_ids:
            var = self.problem.getVarByName(var_id)
            var.setAttr('lb', lb)
            var.setAttr('ub', ub)
            var.setAttr('vtype', vartype_mapping[vartype])
        else:
            self.problem.addVar(name=var_id, lb=lb, ub=ub, vtype=vartype_mapping[vartype])
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

        grb_sense = {'=': GRB.EQUAL,
                     '<': GRB.LESS_EQUAL,
                     '>': GRB.GREATER_EQUAL}

        if constr_id in self.constr_ids:
            constr = self.problem.getConstrByName(constr_id)
            self.problem.remove(constr)

        expr = quicksum(coeff * self.problem.getVarByName(r_id) for r_id, coeff in lhs.items() if coeff)

        self.problem.addConstr(expr, grb_sense[sense], rhs, constr_id)
        self.constr_ids.append(constr_id)

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
                self.problem.remove(self.problem.getVarByName(var_id))
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
                self.problem.remove(self.problem.getConstrByName(constr_id))
                self.constr_ids.remove(constr_id)

    def update(self):
        """ Update internal structure. Used for efficient lazy updating. """
        self.problem.update()

    def set_objective(self, linear=None, quadratic=None, minimize=True):
        """ Set a predefined objective for this problem.

        Args:
            linear (dict): linear coefficients (optional)
            quadratic (dict): quadratic coefficients (optional)
            minimize (bool): solve a minimization problem (default: True)

        Notes:
            Setting the objective is optional. It can also be passed directly when calling **solve**.

        """

        lin_obj = []
        quad_obj = []

        if linear:

            if isinstance(linear, str):
                lin_obj = [1.0 * self.problem.getVarByName(linear)]
                if linear not in self.var_ids:
                    warn(f"Objective variable not previously declared: {linear}")
            else:
                lin_obj = []
                for r_id, val in linear.items():
                    if r_id not in self.var_ids:
                        warn(f"Objective variable not previously declared: {r_id}")
                    elif val != 0:
                        lin_obj.append(val * self.problem.getVarByName(r_id))

        if quadratic:
            quad_obj = []
            for (r_id1, r_id2), val in quadratic.items():
                if r_id1 not in self.var_ids:
                    warn(f"Objective variable not previously declared: {r_id1}")
                elif r_id2 not in self.var_ids:
                    warn(f"Objective variable not previously declared: {r_id2}")
                elif val != 0:
                    quad_obj.append(val * self.problem.getVarByName(r_id1) * self.problem.getVarByName(r_id2))

        obj_expr = quicksum(quad_obj + lin_obj)
        sense = GRB.MINIMIZE if minimize else GRB.MAXIMIZE

        self.problem.setObjective(obj_expr, sense)

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
                    lpvar = problem.getVarByName(r_id)
                    old_constraints[r_id] = (lpvar.lb, lpvar.ub)
                    lpvar.lb = infinity_fix(lb)
                    lpvar.ub = infinity_fix(ub)
                else:
                    warn(f"Constrained variable '{r_id}' not previously declared")
            problem.update()

        self.set_objective(linear, quadratic, minimize)

        # run the optimization
        if pool_size <= 1:

            problem.optimize()

            status = status_mapping.get(problem.status, Status.UNKNOWN)
            message = str(problem.status)

            if status == Status.OPTIMAL:
                fobj = problem.ObjVal
                values, s_prices, r_costs = None, None, None

                if get_values:
                    if isinstance(get_values, Iterable):
                        get_values = list(get_values)
                        values = {r_id: problem.getVarByName(r_id).X for r_id in get_values}
                    else:
                        values = {r_id: problem.getVarByName(r_id).X for r_id in self.var_ids}

                if shadow_prices:
                    s_prices = {m_id: problem.getConstrByName(m_id).Pi for m_id in self.constr_ids}

                if reduced_costs:
                    r_costs = {r_id: problem.getVarByName(r_id).RC for r_id in self.var_ids}

                solution = Solution(status, message, fobj, values, s_prices, r_costs)
            else:
                solution = Solution(status, message)

        else:

            problem.setParam(GRB.Param.PoolSearchMode, 2)
            self.set_parameter(Parameter.POOL_SIZE, pool_size)

            if pool_gap:
                self.set_parameter(Parameter.POOL_GAP, pool_gap)

            problem.optimize()

            status = status_mapping.get(problem.status, Status.UNKNOWN)

            if status == Status.OPTIMAL or status == Status.UNKNOWN:
                solution = self.get_solution_pool()
            else:
                solution = []

        # restore values of temporary constraints
        if constraints:
            for r_id, (lb, ub) in old_constraints.items():
                lpvar = problem.getVarByName(r_id)
                lpvar.lb, lpvar.ub = lb, ub
            problem.update()

        return solution

    def get_solution_pool(self, get_values=True):
        """ Return a solution pool for MILP problems.
        Must be called after using solve with pool_size argument > 0.

        Arguments:
            get_values (bool or list): set to false for speedup if you only care about the objective (default: True)

        Returns:
            list: list of Solution objects

        """
        solutions = []

        for i in range(self.problem.SolCount):
            self.problem.setParam(GRB.param.SolutionNumber, i)
            obj = self.problem.PoolObjVal
            if get_values:
                if isinstance(get_values, Iterable):
                    get_values = list(get_values)
                    values = {r_id: self.problem.getVarByName(r_id).Xn for r_id in get_values}
                else:
                    values = {r_id: self.problem.getVarByName(r_id).Xn for r_id in self.var_ids}
            else:
                values = None
            sol = Solution(fobj=obj, values=values)
            solutions.append(sol)

        return solutions

    def set_parameter(self, parameter, value):
        """ Set a parameter value for this optimization problem

        Arguments:
            parameter (Parameter): parameter type
            value (float): parameter value
        """

        if parameter in parameter_mapping:
            grb_param = parameter_mapping[parameter]
            self.problem.setParam(grb_param, value)
        else:
            raise Exception('Parameter unknown (or not yet supported).')

    def set_logging(self, enabled=False):
        """ Enable or disable log output:

        Arguments:
            enabled (bool): turn logging on (default: False)
        """

        self.problem.setParam('OutputFlag', 1 if enabled else 0)

    def write_to_file(self, filename):
        """ Write problem to file:

        Arguments:
            filename (str): file path
        """

        self.problem.write(filename)
