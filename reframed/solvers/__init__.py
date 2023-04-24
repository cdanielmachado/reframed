solvers = dict()

try:
    from .gurobi_solver import GurobiSolver
    solvers['gurobi'] = GurobiSolver
except ImportError:
    pass


try:
    from .cplex_solver import CplexSolver
    solvers['cplex'] = CplexSolver
except ImportError:
    pass

try:
    from .optlang_solver import OptLangSolver
    solvers['optlang'] = OptLangSolver
except ImportError:
    pass

try:
    from .pulp_solver import SCIP_Solver
    solvers['scip'] = SCIP_Solver
except ImportError:
    pass

try:
    from .pulp_solver import HiGHS_Solver
    solvers['highs'] = HiGHS_Solver
except ImportError:
    pass

try:
    from .pulp_solver import PulpCplex
    solvers['pulp_cplex'] = PulpCplex
except ImportError:
    pass

try:
    from .pulp_solver import PulpGurobi
    solvers['pulp_gurobi'] = PulpGurobi
except ImportError:
    pass


default_solver = None


def get_default_solver():

    global default_solver

    if default_solver:
        return default_solver

    solver_order = ['gurobi', 'cplex', 'highs', 'scip', 'optlang']

    for solver in solver_order:
        if solver in list(solvers.keys()):
            default_solver = solver
            break

    if not default_solver:
        raise RuntimeError("No solver available.")

    return default_solver


def set_default_solver(solvername):
    """ Sets default solver.

    Arguments:
        solvername : (str) solver name (currently available: 'gurobi', 'cplex')
    """

    global default_solver

    if solvername.lower() in list(solvers.keys()):
        default_solver = solvername.lower()
    else:
        raise RuntimeError(f"Solver {solvername} not available.")


def solver_instance(model=None):
    """ Returns a new instance of the currently selected solver.

    Arguments:
        model : CBModel (optional) -- immediatly instantiate problem with given model

    Returns:
        Solver
    """

    solver = get_default_solver()

    if solver:
        return solvers[solver](model)
