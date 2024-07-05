from enum import Enum
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

class Stage(Enum):
    INIT         =  0        # SCIP data structures are initialized, no problem exists */
    PROBLEM      =  1        # the problem is being created and modified */
    TRANSFORMING =  2        # the problem is being transformed into solving data space */
    TRANSFORMED  =  3        # the problem was transformed into solving data space */
    INITPRESOLVE =  4        # presolving is initialized */
    PRESOLVING   =  5        # the problem is being presolved */
    EXITPRESOLVE =  6        # presolving is exited */
    PRESOLVED    =  7        # the problem was presolved */
    INITSOLVE    =  8        # the solving process data is being initialized */
    SOLVING      =  9        # the problem is being solved */
    SOLVED       = 10        # the problem was solved */
    EXITSOLVE    = 11        # the solving process data is being freed */
    FREETRANS    = 12        # the transformed problem is being freed */
    FREE         = 13 


class SCIPSolver(Solver):
    """ Implements the solver interface using SCIP. """

    def __init__(self, model=None):
        Solver.__init__(self)
        self.problem = Model()
        self.problem.hideOutput()
        self.problem.enableReoptimization()
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
        
        if self.problem.getStage() > Stage.PROBLEM.value:
            self.problem.freeReoptSolve()

        if isinstance(objective, str):
            objective = {objective: 1.0}
        
        objective = quicksum(coeff * self._vars_dict[var_id] for var_id, coeff in objective.items() if coeff != 0)

        #TODO: this will fail to update the coefficient of variables that were added after the problem was solved the first time!
        # this is an issue for methods like pFBA
        self.problem.chgReoptObjective(objective, sense='minimize' if minimize else 'maximize')

        self.objective = objective
        self.minimize = minimize


    def internal_solve(self):
        
        self.problem.optimize()
        status = self.problem.getStatus()
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

        if self.problem.getStage() > Stage.PROBLEM.value:
            raise Exception('This wont work with SCIP.')

        old_bounds = {}

        for r_id, x in constraints.items():
            lb, ub = x if isinstance(x, tuple) else (x, x)
            old_bounds[r_id] = (self._vars_dict[r_id].getLbOriginal(), self._vars_dict[r_id].getUbOriginal())
            self.problem.chgVarLb(self._vars_dict[r_id], None if lb == -inf else lb)
            self.problem.chgVarUb(self._vars_dict[r_id], None if ub == inf else ub)

        return old_bounds
    

    def reset_bounds(self, old_bounds):

        if self.problem.getStage() > Stage.PRESOLVED:
            print('reset_bounds')
            self.problem.freeReoptSolve()

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


    def print_stage(self):
        print('Current stage:', Stage(self.problem.getStage()).name)