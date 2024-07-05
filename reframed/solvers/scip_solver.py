from math import inf
import operator
from .solver import Solver, VarType, Parameter
from .solution import Solution, Status
from pyscipopt import Model, quicksum
from warnings import warn


status_mapping = {
   "optimal": Status.OPTIMAL,
   "infeasible": Status.INFEASIBLE,
   "inforunbd": Status.INF_OR_UNB,
   "unbounded": Status.UNBOUNDED,
   "timelimit": Status.SUBOPTIMAL,
   "nodelimit": Status.SUBOPTIMAL,
   "timelimit": Status.SUBOPTIMAL,
   "totalnodelimit": Status.SUBOPTIMAL,
   "stallnodelimit": Status.SUBOPTIMAL,
   "gaplimit": Status.SUBOPTIMAL,
   "memlimit": Status.SUBOPTIMAL,
   "sollimit": Status.SUBOPTIMAL,
   "bestsollimit": Status.SUBOPTIMAL,
   "restartlimit": Status.SUBOPTIMAL,
}


parameter_mapping = {
    Parameter.TIME_LIMIT: 'limits/time',
    Parameter.FEASIBILITY_TOL: 'numerics/feastol',
    Parameter.OPTIMALITY_TOL: 'numerics/barrierconvtol',
    Parameter.MIP_ABS_GAP: 'limits/absgap',
    Parameter.MIP_REL_GAP: 'limits/gap',
    Parameter.POOL_SIZE: 'limits/maxsol',
}


vartype_mapping = {
    VarType.BINARY: 'B',
    VarType.INTEGER: 'I',
    VarType.CONTINUOUS: 'C',
}


sense_mapping = {
    '=': operator.eq,
    '<': operator.le,
    '>': operator.ge,
}

class SCIPSolver(Solver):
    """ Implements the solver interface using SCIP. """

    def __init__(self, model=None):
        Solver.__init__(self)
        self.problem = Model()
        self.problem.hideOutput()
        self.problem.enableReoptimization()
        self.updatable = True
        self._vars_dict = {}
        self._cons_dict = {}

        if model:
            self.build_problem(model)

    def add_variables(self, var_dict):

        for var_id, (lb, ub, vartype) in var_dict.items():

            lb = None if lb == -inf else lb
            ub = None if ub == inf else ub

            self._vars_dict[var_id] = self.problem.addVar(name=var_id, lb=lb, ub=ub, vtype=vartype_mapping[vartype])
        
        self.variables.extend(var_dict.keys())

    def add_constraints(self, constr_dict):

        for constr_id, (lhs, sense, rhs) in constr_dict.items():
            expr = quicksum(self._vars_dict[var_id] * coeff for var_id, coeff in lhs.items())
            constr = sense_mapping[sense](expr,rhs)
            self._cons_dict[constr_id] = self.problem.addCons(constr, name=constr_id)

        self.constraints.extend(constr_dict.keys())
 

    def set_objective(self, objective, minimize=True):
        
        self.check_stage()

        if isinstance(objective, str):
            objective = {objective: 1.0}
        
        objective = quicksum(coeff * self._vars_dict[var_id] for var_id, coeff in objective.items() if coeff != 0)

        #TODO: this will fail to update the coefficient of variables that were added after the problem was solved the first time!
        self.problem.chgReoptObjective(objective, sense='minimize' if minimize else 'maximize')

        self.objective = objective
        self.minimize = minimize


    def internal_solve(self):
        
        self.problem.optimize()
        status = self.problem.getStatus()
        self.updatable = False
        return status_mapping.get(status, Status.UNKNOWN)


    def get_solution(self, status, get_values=True, shadow_prices=False, reduced_costs=False):

        _solution = self.problem.getBestSol()

        fobj = self.problem.getObjVal()

        if get_values:
            if isinstance(get_values, list):
                values = {r_id: _solution[self._vars_dict[r_id]] for r_id in get_values}
            else:
                values = {var_id: _solution[var] for var_id, var in self._vars_dict.items()}
        else:
            values = None

        if shadow_prices:
            s_prices = {constr_id: self.problem.getDualSolVal(constr) for constr_id, constr in self._cons_dict.items()}
        else:
            s_prices = None

        if reduced_costs:
            r_costs = {var_id: self.problem.getVarRedCost(var) for var_id, var in self._vars_dict.items()}
        else:
            r_costs = None

        return Solution(status, fobj=fobj, values=values, shadow_prices=s_prices, reduced_costs=r_costs)

    
    def set_temporary_bounds(self, constraints):

        self.check_stage()

        old_bounds = {}

        for r_id, x in constraints.items():
            lb, ub = x if isinstance(x, tuple) else (x, x)
            old_bounds[r_id] = (self._vars_dict[r_id].getLbOriginal(), self._vars_dict[r_id].getUbOriginal())
            self.problem.chgVarLb(self._vars_dict[r_id], None if lb == -inf else lb)
            self.problem.chgVarUb(self._vars_dict[r_id], None if ub == inf else ub)

        return old_bounds
    

    def reset_bounds(self, old_bounds):

        self.check_stage()

        for r_id, (lb, ub) in old_bounds.items():
             self.problem.chgVarLb(self._vars_dict[r_id], lb)
             self.problem.chgVarUb(self._vars_dict[r_id], ub)


    def set_parameter(self, parameter, value):
        """ Set a parameter value for this optimization problem

        Arguments:
            parameter (Parameter): parameter type
            value (float): parameter value
        """

        if parameter in parameter_mapping:
            scip_param = parameter_mapping[parameter]
            self.problem.setParam(scip_param, value)
        else:
            raise Exception('Parameter unknown (or not yet supported).')
        

    def write_to_file(self, filename):
        """ Write problem to file:

        Arguments:
            filename (str): file path
        """

        self.problem.writeProblem(filename)

    def set_logging(self, enabled=False):
        """ Enable or disable log output:

        Arguments:
            enabled (bool): turn logging on (default: False)
        """

        self.problem.hideOutput(quiet=(not enabled))

    def check_stage(self):
        if not self.updatable:
            self.problem.freeReoptSolve()
        self.updatable = True

    def print_stage(self):
        stages = [
            'INIT', 
            'PROBLEM', 
            'TRANSFORMING', 
            'TRANSFORMED', 
            'INITPRESOLVE', 
            'PRESOLVING', 
            'EXITPRESOLVE', 
            'PRESOLVED', 
            'INITSOLVE', 
            'SOLVING', 
            'SOLVED', 
            'EXITSOLVE', 
            'FREETRANS', 
            'FREE', 
        ]
        print('Current stage:', stages[self.problem.getStage()])