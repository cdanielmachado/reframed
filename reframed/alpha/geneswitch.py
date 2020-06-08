from ..cobra.simulation import FBA, ROOM, pFBA
from ..cobra.variability import FVA
from .GPRtransform import gpr_transform
from math import inf


def GeneSwitch(model, target, constraints=None, production=0.5, delta=0.1, epsilon=0.1, solutions=1, use_pool=False):
    """**GeneSwitch** gives minimal number of genes for which the respective enzymes need to undergo flux change
        to make phenotype switch form growth optimised to producing optimised.
        This is an adaptation of the ROOM algorithm using the GPR transformation approach (Machado et al, 2016).

    Author: Daniel Machado, Astrid Stubbusch (Original method: Paula Jouhten)

    Arguments:
        model (CBModel): a constraint-based model
        target (str): target reaction id
        constraints (dict): environmental or additional constraints (optional)
        production (float): minimum fraction of maximum target production rate (default: 0.5)
        delta (float): relative tolerance (default: 0.1)
        epsilon (float): absolute tolerance (default: 0.1)
        solutions (int): number of solutions to compute (default: 1)
        use_pool (bool) use solver solution pool (default: False)

    Returns:
        dict: significantly up / down regulated genes and respective ratios
        dict: respective Solution objects
    """

    if constraints is None:
        constraints = {}

    mut_sol = FBA(model, constraints=constraints, objective={target: 1})
    v_max = mut_sol.values[target]
    switch = {target: (production * v_max, v_max)}

    model_ext = gpr_transform(model, add_proteome=True, inplace=False)
    constraints = model_ext.convert_constraints(constraints)
    switch = model_ext.convert_constraints(switch)

    wt_sol = pFBA(model_ext, constraints=constraints)
    model_ext.reactions.proteome_synth.ub = wt_sol.values["proteome_synth"]

    reference = FVA(model_ext, obj_frac=0.9, reactions=model_ext.u_reactions, constraints=constraints)

    constraints.update(switch)

    sols = ROOM(model_ext, reference=reference, constraints=constraints, reactions=model_ext.u_reactions,
                delta=delta, epsilon=epsilon, solutions=solutions, use_pool=use_pool)

    if solutions == 1:
        sols = [sols]

    switches = []

    for sol in sols:

        up = {}
        down = {}

        for g_id in model.genes:
            if 'y_u_' + g_id[2:] not in sol.values:
                continue

            y_val = sol.values['y_u_' + g_id[2:]]
            u_mut = sol.values['u_' + g_id[2:]]
            u_wt_lb = reference['u_' + g_id[2:]][0]
            u_wt_ub = reference['u_' + g_id[2:]][1]

            if y_val > 0.5:
                if u_mut > u_wt_ub + epsilon:
                    up[g_id] = u_mut / u_wt_ub if u_wt_ub > 0 else inf
                elif u_mut < u_wt_lb - epsilon:
                    try:
                        down[g_id] = u_mut / u_wt_lb
                    except:
                        print(u_mut, u_wt_lb)

        sol.values_ext = sol.values
        sol.values = model_ext.convert_fluxes(sol.values)

        switches.append((up, down))

    return switches, sols


def allostericSwitch(model, target, constraints=None, production=0.5, delta=0.1, epsilon=0.1, solutions=1,
                     use_pool=False, self_inhibition=True):
    """**allostericSwitch** is an extended version of **GeneSwitch** that finds allosteric regulation targets that
    could require de-regulation in order to reach the desired target production.

    Author: Daniel Machado

    Arguments:
        model (CBModel): a constraint-based model
        target (str): target reaction id
        constraints (dict): environmental or additional constraints (optional)
        production (float): minimum fraction of maximum target production rate (default: 0.5)
        delta (float): relative tolerance (default: 0.1)
        epsilon (float): absolute tolerance (default: 0.1)
        solutions (int): number of solutions to compute (default: 1)
        use_pool (bool): use solver solution pool (default: False)
        self_inhibition (bool): include reaction self-inhibition by its substrates or products (default: True)

    Returns:
        dict: up/down-regulation targets and inhibitory interactions
        dict: respective Solution objects
    """

    ubiquity_threshold = 20

    wt_sol = pFBA(model, constraints=constraints)
    t0 = wt_sol.get_metabolites_turnover(model)
    v0 = wt_sol.values

    switches, sols = GeneSwitch(model, target, constraints=constraints, production=production,
                                delta=delta, epsilon=epsilon, solutions=solutions, use_pool=use_pool)

    reg_switches = []

    def ubiquitous(m_id):
        return len(model.get_metabolite_reactions(m_id)) > ubiquity_threshold

    for (up, down), sol in zip(switches, sols):
        t = sol.get_metabolites_turnover(model)
        v = sol.values
        inhibition = {}

        for r_id, rxn in model.reactions.items():
            if v[r_id] < v0[r_id] * (1 + delta) + epsilon:
                continue

            for m_id in rxn.get_inhibitors():
                if ubiquitous(m_id):
                    continue

                if not self_inhibition and m_id in rxn.stoichiometry:
                    continue

                if t[m_id] > t0[m_id] * (1 + delta) + epsilon:
                    if r_id not in inhibition:
                        inhibition[r_id] = []
                    inhibition[r_id].append(m_id)

        reg_switches.append((up, down, inhibition))

    return reg_switches, sols