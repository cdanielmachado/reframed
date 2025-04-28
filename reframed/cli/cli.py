import argparse
import textwrap
from math import inf
import pandas as pd

from reframed import set_default_solver, load_cbmodel, Environment, FBA, pFBA, FVA
from reframed.cobra.knockout import gene_knockout
from reframed.solvers import solver_instance
from reframed.solvers.solution import Status


def load_media_db(filename):
    """ Load media library file. """

    data = pd.read_csv(filename, sep='\t')

    has_bounds = 'bound' in data.columns

    media_db = {}
    for medium, data_i in data.groupby('medium'):
        if has_bounds:
            media_db[medium] = dict(data_i[['compound', 'bound']].values)
        else:
            media_db[medium] = list(data_i['compound'])

    return media_db, has_bounds

def create_env(model, medium, media_db, media_has_bounds, max_uptake):
    env = Environment.from_compounds(media_db[medium], max_uptake=max_uptake)
    env = env.apply(model, inplace=False, exclusive=True, warning=False)

    if media_has_bounds:
        for cpd, bound in media_db[medium].items():
            r_id = f'R_EX_{cpd}_e'
            if r_id in env:
                env[r_id] = (-bound, inf)

    return env

def load_constraints(filename):
    data = pd.read_csv(filename, sep='\t', header=None)
    return {x[0]: (x[1], x[2]) for x in data.values}


def save_solutions(solutions, output=None):

    if not output:
        output = 'result.tsv'

    if not solutions:
        print('No solutions found.')
    else:
        df = pd.DataFrame(solutions)
        df.sort_index(inplace=True)
        df.to_csv(output, header=(len(df.columns) > 1), sep='\t')

def run(modelfile, objective=None, method='FBA', growth_frac=0, knockout=None, constraints=None, output=None, 
        media=None, max_uptake=10, mediadb=None):
    
    model = load_cbmodel(modelfile)

    solver = solver_instance(model)

    if media is None:
        media = [None]
    else:
        media = media.split(',')

        if mediadb is None:
            raise RuntimeError('Media database file must be provided.')
        else:   
            media_db, media_has_bounds = load_media_db(mediadb)

    if constraints:
        add_constraints = load_constraints(constraints)

    if knockout:
        knockout = knockout.split(',')

    solutions = {}


    for medium in media:
        
        if medium:
            env = create_env(model, medium, media_db, media_has_bounds, max_uptake)
        else:
            medium = 'default'
            env = {}

        if constraints is not None:
            env.update(add_constraints)

        if knockout:
            for genes in knockout:
                to_delete = [f'G_{x}' for x in genes.split('+')]
                sol = gene_knockout(model, to_delete, method=method, constraints=env)#, solver=solver) can't reuse solver until remove_constraint is implemented
                if sol and (sol.status == Status.OPTIMAL or sol.status == Status.SUBOPTIMAL):
                    solutions[f'{medium}_{genes}'] = sol.values
        else:
            sol = None

            if method == 'FBA':
                sol = FBA(model, objective=objective, constraints=env, solver=solver)
            elif method == 'pFBA':
                sol = pFBA(model, objective=objective, constraints=env, solver=solver)
            elif method == 'FVA':
                var = FVA(model, constraints=env, obj_frac=growth_frac, solver=solver)
                solutions[f'{medium}_lb'] = {x: y[0] for x, y in var.items()}
                solutions[f'{medium}_ub'] = {x: y[1] for x, y in var.items()}

            if sol and (sol.status == Status.OPTIMAL or sol.status == Status.SUBOPTIMAL):
                solutions[f'{medium}'] = sol.values

    save_solutions(solutions, output)


def main():

    parser = argparse.ArgumentParser(description="Run flux balance analysis (FBA) for a given COBRA model",
                                     formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument('model', metavar='MODEL', help="SBML model")
    parser.add_argument('--obj', '--objective', metavar='REACTION ID', dest='objective', 
                        help="Select a different objective function (than the one specified in the model)")

    method = parser.add_mutually_exclusive_group()

    method.add_argument('-p', '--parsimonious', dest='parsimonious', action='store_true', help="Use parsimonious FBA (pFBA)")
    method.add_argument('--fva', metavar='FRACTION', type=float, nargs='?', const=0, default=None,
                        help="Run flux variability analysis at a given fraction of the optimal objective (default: 0.0, max: 1.0)")
    method.add_argument('--lmoma', action='store_true', help="Use linear MOMA to simulate gene knockouts")
    
    parser.add_argument('--ko', '--knockout', metavar='GENES', dest='knockout', help=textwrap.dedent(
        """
        List of genes to knockout
        Separate with ',' for individual deletions or '+' for simultaneous deletions
        
        Multiple combinations are possible:
            --knockout A,B,A+B (delete A, delete B, delete A and B)
        """
    ))

    parser.add_argument('-c', '--constraints', dest='constraints', metavar='FILENAME', help=textwrap.dedent(
        """
        TSV file with additional constraints to the solution space.
        Can also be used to override the uptake rates from the growth medium. 

        Example:
            R_PGI\t0\t1000
            R_PFK\t0\t0
            R_EX_o2_e\t-20\t0
        """
    ))

    parser.add_argument('-o', '--output', metavar='FILENAME',dest='output', help="Output file")
    parser.add_argument('-m', '--medium', dest='medium', help=textwrap.dedent(
        """
        Specify a growth medium for simulation (otherwise uses default conditions in the SBML file)
        Multiple media can be specified (comma separated) to run several simulations. Example: -m M9,LB 
        """
    ))
    parser.add_argument('-u', '--uptake', type=float, default=10.0, metavar='RATE', dest='uptake', 
                        help="Maximum uptake rate for compounds in the medium (default: 10 mmol/gDW/h).")
    parser.add_argument('--mediadb', metavar='FILENAME', help=textwrap.dedent(
        """
        Media database (TSV file).
        
        Two columns are mandatory, 'medium' and 'compound' (but other columns can be present).
        Each compound will have the maximum uptake rate defined by the -u argument.
        
        Example:

        medium\tcompound\tnotes
        M1\tglc__D\tmain carbon source
        M1\to2\tadded for aerobic growth
        M1\tnh4\tnitrogen source
        M2\tglyc\tcarbon source for medium 2
        M2\tnh4\tnitrogen source

        Alternatively, individual flux bounds can be specified for each compound in a new column
        called 'bound'. Column order is not important.

        Example:

        medium\tcompound\tbound
        M1\tglc__D\t10
        M1\to2\t20
        M1\tnh4\t1000
        M2\tglyc\t10
        M2\tnh4\t1000

        """
    ))
    parser.add_argument('--solver', help="Select LP solver (options: gurobi, cplex, scip)")

    args = parser.parse_args()

    if args.solver is not None:
        set_default_solver(args.solver)

    method = 'FBA'

    if args.parsimonious:
        method = 'pFBA'
    elif args.lmoma:
        method = 'lMOMA'
    elif args.fva is not None:
        method = 'FVA'

    if method == 'FVA' and args.knockout:
        raise RuntimeError('FVA with gene knockouts not currently implemented.')
    
    if args.objective and args.knockout:
        raise RuntimeError('Changing objective for gene knockouts not currently implemented.')

    if args.objective and method == 'FVA':
        raise RuntimeError('Changing objective for FVA not currently implemented.')
    
    if method == 'lMOMA' and not args.knockout: 
        raise RuntimeError('lMOMA requires at least one gene knockout.')

    run(
        modelfile=args.model,
        objective=args.objective,
        method=method,
        growth_frac=args.fva,
        knockout=args.knockout,
        constraints=args.constraints,
        output=args.output,
        media=args.medium,
        max_uptake=args.uptake,
        mediadb=args.mediadb,
    )


if __name__ == '__main__':
    main()