from ..solvers import solver_instance
from ..solvers.solution import Status
from ..cobra.knockout import essentiality, reaction_knockout
from ..core.environment import Environment
from ..cobra.simulation import FBA
from warnings import warn


def gapMake2(model, growth_envs, no_growth_envs, kind='reactions'):

    g_essential = []
    ng_essential = []
    solver = solver_instance(model)

    for env in growth_envs:
        essential = essentiality(model, kind,  constraints=env, solver=solver)
        g_essential.append(essential)

    for env in no_growth_envs:
        essential = essentiality(model, kind,  constraints=env, solver=solver)
        ng_essential.append(essential)

    return g_essential, ng_essential


def gapMake(model, growth_envs, no_growth_envs, reactions=None, min_growth=0.001):

    if reactions is None:
        reactions = model.reactions

    solver = solver_instance(model)

    def is_feasible(solution):
        return solution.status == Status.OPTIMAL and solution.fobj > min_growth

    def is_essential(x, y):
        solution = reaction_knockout(model, [x], constraints=y, solver=solver)
        return not is_feasible(solution)

    essential = set()
    non_essential = set(reactions)

    for env_id, env in growth_envs.items():
        sol = FBA(model, constraints=env, solver=solver)
        if not is_feasible(sol):
            warn(f"Infeasible growth in medium {env_id}")
        else:
            for r_id in non_essential:
                if is_essential(r_id, env):
                    essential.add(r_id)
        non_essential -= essential

    cut_sets = {}

    for env_id, env in no_growth_envs.items():
        sol = FBA(model, constraints=env, solver=solver)
        if not is_feasible(sol):
            continue

        for r_id in non_essential:
            if is_essential(r_id, env):
                if r_id in cut_sets:
                    cut_sets[r_id].append(env_id)
                else:
                    cut_sets[r_id] = [env_id]
                break

    return cut_sets


def env_generator(model, medium, media_db):
    cpds = media_db.query(f"medium == '{medium}'")["compound"]
    env = Environment.from_compounds(cpds)
    env = env.apply(model, exclusive=True, inplace=False, warning=False)
    return env


def gapMakeWrapper(model, g_media, ng_media, media_db):

    if len(ng_media) == 0:
        warn("Nothing to do.")
        return dict()

    g_envs = {medium: env_generator(model, medium, media_db) for medium in g_media}
    ng_envs = {medium: env_generator(model, medium, media_db) for medium in ng_media}
    reactions = set(model.reactions) - set(model.get_exchange_reactions())
    return gapMake(model, g_envs, ng_envs, reactions)

