from .solver import Solver, VarType, Parameter
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

sense_mapping = {
    '=': GRB.EQUAL, 
    '<': GRB.LESS_EQUAL, 
    '>': GRB.GREATER_EQUAL
}



class GurobiSolver(Solver):
    """ Implements the gurobi solver interface. """

    def __init__(self, model=None):
        Solver.__init__(self)
        self.problem = GurobiModel()
        self.set_logging(False)
        if model:
            self.build_problem(model)

    def add_variables(self, var_dict):

        for var_id, (lb, ub, vartype) in var_dict.items():
            self.problem.addVar(name=var_id, lb=infinity_fix(lb), ub=infinity_fix(ub), vtype=vartype_mapping[vartype])
        
        self.variables.extend(var_dict.keys())
        self.problem.update()


    def add_constraints(self, constr_dict):

        for constr_id, (lhs, sense, rhs) in constr_dict.items():
            expr = quicksum(coeff * self.problem.getVarByName(r_id) for r_id, coeff in lhs.items() if coeff)
            self.problem.addConstr(expr, sense_mapping[sense], rhs, constr_id)
        
        self.constraints.extend(constr_dict.keys())
        self.problem.update()

    def remove_constraint(self, constr_id):
        """ Remove a constraint from the current problem.

        Arguments:
            constr_id (str): constraint identifier
        """

        if constr_id in self.constraints:
            self.problem.remove(self.problem.getConstrByName(constr_id))
            self.constraints.remove(constr_id)


    def set_objective(self, objective, minimize=True):

        if isinstance(objective, str):
            objective = {objective: 1.0}

        obj_expr = quicksum([coeff * self.problem.getVarByName(r_id) for r_id, coeff in objective.items() if coeff != 0])
        sense = GRB.MINIMIZE if minimize else GRB.MAXIMIZE

        self.problem.setObjective(obj_expr, sense)

        self.objective = objective
        self.minimize = minimize

    def internal_solve(self):
        self.problem.optimize()
        status = status_mapping.get(self.problem.status, Status.UNKNOWN)
        return status
    
    def get_solution(self, status, get_values=True, shadow_prices=False, reduced_costs=False):

        fobj = self.problem.ObjVal

        if get_values:
            if isinstance(get_values, list):
                values = {r_id: self.problem.getVarByName(r_id).X for r_id in get_values}
            else:
                values = {r_id: self.problem.getVarByName(r_id).X for r_id in self.variables}
        else:
            values = None

        if shadow_prices:
            s_prices = {m_id: self.problem.getConstrByName(m_id).Pi for m_id in self.constraints}
        else:
            s_prices = None

        if reduced_costs:
            r_costs = {r_id: self.problem.getVarByName(r_id).RC for r_id in self.variables}
        else:
            r_costs = None

        return Solution(status, fobj=fobj, values=values, shadow_prices=s_prices, reduced_costs=r_costs)
    


    def set_temporary_bounds(self, constraints):

        old_constraints = {}
        for r_id, x in constraints.items():
            lb, ub = x if isinstance(x, tuple) else (x, x)
            if r_id in self.variables:
                lpvar = self.problem.getVarByName(r_id)
                old_constraints[r_id] = (lpvar.lb, lpvar.ub)
                lpvar.lb = infinity_fix(lb)
                lpvar.ub = infinity_fix(ub)
            else:
                warn(f"Constrained variable '{r_id}' not previously declared")
        
        self.problem.update()
        
        return old_constraints


    def reset_bounds(self, bounds):

        for r_id, (lb, ub) in bounds.items():
            lpvar = self.problem.getVarByName(r_id)
            lpvar.lb, lpvar.ub = lb, ub
        self.problem.update()


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
