
from .solver import Solver, VarType, Parameter
from .solution import Solution, Status
from cplex import Cplex, infinity, SparsePair
import sys
from math import inf
from warnings import warn


def infinity_fix(val):
    if val == inf:
        return infinity
    elif val == -inf:
        return -infinity
    else:
        return float(val)


class CplexSolver(Solver):
    """ Implements the solver interface using CPLEX. """

    def __init__(self, model=None):
        Solver.__init__(self)
        self.problem = Cplex()
        self._cached_lower_bounds = {}
        self._cached_upper_bounds = {}

        self.status_mapping = {
            self.problem.solution.status.optimal: Status.OPTIMAL,
            self.problem.solution.status.optimal_tolerance: Status.OPTIMAL,
            self.problem.solution.status.unbounded: Status.UNBOUNDED,
            self.problem.solution.status.infeasible: Status.INFEASIBLE,
            self.problem.solution.status.infeasible_or_unbounded: Status.INF_OR_UNB,
            self.problem.solution.status.MIP_optimal: Status.OPTIMAL,
            self.problem.solution.status.MIP_unbounded: Status.UNBOUNDED,
            self.problem.solution.status.MIP_infeasible: Status.INFEASIBLE,
            self.problem.solution.status.MIP_infeasible_or_unbounded: Status.INF_OR_UNB
        }

        self.vartype_mapping = {
            VarType.BINARY: self.problem.variables.type.binary,
            VarType.INTEGER: self.problem.variables.type.integer,
            VarType.CONTINUOUS: self.problem.variables.type.continuous
        }

        self.parameter_mapping = {
            Parameter.TIME_LIMIT: self.problem.parameters.timelimit,
            Parameter.FEASIBILITY_TOL: self.problem.parameters.simplex.tolerances.feasibility,
            Parameter.OPTIMALITY_TOL: self.problem.parameters.simplex.tolerances.optimality,
            Parameter.INT_FEASIBILITY_TOL: self.problem.parameters.mip.tolerances.integrality,
            Parameter.MIP_ABS_GAP: self.problem.parameters.mip.tolerances.mipgap,
            Parameter.MIP_REL_GAP: self.problem.parameters.mip.tolerances.absmipgap,
            Parameter.POOL_SIZE: self.problem.parameters.mip.limits.populate,
            Parameter.POOL_GAP: self.problem.parameters.mip.pool.relgap
        }

        self.set_logging(False)

        if model:
            self.build_problem(model)


    def add_variables(self, var_dict):

        var_ids = list(var_dict.keys())

        lbs = [infinity_fix(lb) for (lb, _, _) in var_dict.values()]
        ubs = [infinity_fix(ub) for (_, ub, _) in var_dict.values()]
        vartypes = [self.vartype_mapping[vartype] for (_, _, vartype) in var_dict.values()]

        self.problem.variables.add(names=var_ids, lb=lbs, ub=ubs, types=vartypes)

        self.variables.extend(var_ids)
        self._cached_lower_bounds.update(dict(zip(var_ids, lbs)))
        self._cached_upper_bounds.update(dict(zip(var_ids, ubs)))


    def add_constraints(self, constr_dict):


        constr_ids = list(constr_dict.keys())
        lhs_all = [SparsePair(ind=list(lhs.keys()), val=list(lhs.values())) for (lhs, _, _) in constr_dict.values()]
        map_sense = {'=': 'E', '<': 'L', '>': 'G'}
        sense_all = [map_sense[sense] for (_, sense, _) in constr_dict.values()]
        rhs_all = [rhs for (_, _, rhs) in constr_dict.values()]

        self.problem.linear_constraints.add(lin_expr=lhs_all, senses=sense_all, rhs=rhs_all, names=constr_ids)
        self.constraints.extend(constr_ids)


    def set_objective(self, objective, minimize=True):

        self.update()

        if isinstance(objective, str):
            objective = {objective: 1.0}

        updated_coeffs = {}

        for var_id in self.variables:
            if objective.get(var_id, 0) != self.objective.get(var_id, 0):
                updated_coeffs[var_id] = objective.get(var_id, 0)

        if updated_coeffs:
            self.problem.objective.set_linear(list(updated_coeffs.items()))

        self.objective = objective

        if minimize != self.minimize:
            sense = self.problem.objective.sense.minimize if minimize else self.problem.objective.sense.maximize
            self.problem.objective.set_sense(sense)
            self.sense = minimize


    def internal_solve(self):
        self.problem.solve()
        status = self.status_mapping.get(self.problem.solution.get_status(), Status.UNKNOWN)
        return status


    def get_solution(self, status, get_values=True, shadow_prices=False, reduced_costs=False):

        fobj = self.problem.solution.get_objective_value()

        if get_values:
            if isinstance(get_values, list):
                values = dict(zip(get_values, self.problem.solution.get_values(get_values)))
            else:
                values = dict(zip(self.variables, self.problem.solution.get_values()))
        else:
            values = None

        if shadow_prices:
            s_prices = dict(zip(self.constraints, self.problem.solution.get_dual_values(self.self.constraints)))
        else:
            s_prices = None

        if reduced_costs:
            r_costs = dict(zip(self.variables, self.problem.solution.get_reduced_costs(self.variables)))
        else:
            r_costs = None

        solution = Solution(status, fobj=fobj, values=values, shadow_prices=s_prices, reduced_costs=r_costs)

        return solution

    def set_temporary_bounds(self, constraints):

        lb_new, ub_new = {}, {}
        lb_reset, ub_reset = [], []

        for r_id, x in constraints.items():
            if r_id in self.variables:
                if isinstance(x, tuple):
                    lb, ub = infinity_fix(x[0]), infinity_fix(x[1])
                else:
                    lb, ub = infinity_fix(x), infinity_fix(x)
                
                if lb != self._cached_lower_bounds[r_id]:
                    lb_new[r_id] = lb
                    lb_reset.append(r_id)
                if ub != self._cached_upper_bounds[r_id]:
                    ub_new[r_id] = ub
                    ub_reset.append(r_id)          
            else:
                warn(f"Constrained variable not previously declared: {r_id}")

        if len(lb_new) > 0:
            self.problem.variables.set_lower_bounds(lb_new.items())

        if len(ub_new) > 0:
            self.problem.variables.set_upper_bounds(ub_new.items())

        return (lb_reset, ub_reset)


    def reset_bounds(self, bounds):

        lb_reset, ub_reset = bounds

        lb_old = [(r_id, self._cached_lower_bounds[r_id]) for r_id in lb_reset]
        if len(lb_old) > 0:
            self.problem.variables.set_lower_bounds(lb_old)

        ub_old = [(r_id, self._cached_upper_bounds[r_id]) for r_id in ub_reset]
        if len(ub_old) > 0:
            self.problem.variables.set_upper_bounds(ub_old)


    def set_parameter(self, parameter, value):
        """ Set a parameter value for this optimization problem

        Arguments:
            parameter (Parameter): parameter type
            value (float): parameter value
        """

        if parameter in self.parameter_mapping:
            self.parameter_mapping[parameter].set(value)
        else:
            raise RuntimeError('Parameter unknown (or not yet supported).')

    def set_logging(self, enabled=False):
        """ Enable or disable log output:

        Arguments:
            enabled (bool): turn logging on (default: False)
        """

        if enabled:
            self.problem.set_log_stream(sys.stdout)
            self.problem.set_error_stream(sys.stderr)
            self.problem.set_warning_stream(sys.stderr)
            self.problem.set_results_stream(sys.stdout)
        else:
            self.problem.set_log_stream(None)
            self.problem.set_error_stream(None)
            self.problem.set_warning_stream(None)
            self.problem.set_results_stream(None)

    def write_to_file(self, filename):
        """ Write problem to file:

        Arguments:
            filename (str): file path
        """

        self.problem.write(filename)

