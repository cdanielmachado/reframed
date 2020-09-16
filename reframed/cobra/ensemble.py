from collections import OrderedDict

from ..solvers import solver_instance
from ..core.transformation import disconnected_metabolites
from ..solvers.solution import Status
from ..io.sbml import load_cbmodel, save_cbmodel
from .simulation import FBA, pFBA


class EnsembleModel(object):
    """ An Ensemble Model represents a collection of models that differ by a few reactions.
    They can be used to account for uncertainty in the structure of the metabolic network.

    """

    def __init__(self, model, size, reaction_states=None):
        self.model = model.copy()
        self.size = size
        self.reaction_states = {}

        if reaction_states:
            for r_id, states in reaction_states.items():
                assert r_id in model.reactions
                assert len(states) == size
                self.reaction_states[r_id] = states[:]
        else:
            self.reaction_states = {r_id: [True]*size for r_id in model.reactions}

    def get_reaction_states(self, i):
        return {r_id: self.reaction_states[r_id][i] if r_id in self.reaction_states else True
                for r_id in self.model.reactions}

    def get_constraints(self, i):
        return {r_id: (0, 0) for r_id in self.model.reactions if r_id in self.reaction_states and
                not self.reaction_states[r_id][i]}

    def simplify(self):
        inactive = [r_id for r_id, states in self.reaction_states.items() if not any(states)]
        self.model.remove_reactions(inactive)
        del_metabolites = disconnected_metabolites(self.model)
        self.model.remove_metabolites(del_metabolites)

        for r_id in inactive:
            del self.reaction_states[r_id]


def simulate_ensemble(ensemble, method='FBA', constraints=None, solver=None, get_fluxes=True):
    """ Simulate an Ensemble Model

    Args:
        ensemble (EnsembleModel): ensemble model
        method (str): simulation method (default: 'FBA')
        constraints (dict): additional constraints (optional)
        solver (Solver): solver instance (optional)
        get_fluxes (bool): if True returns flux distributions for all models,
            otherwise only the objective function values are returned (default: True)

    Returns:
        ensemble flux distributions or ensemble objective function values

    """

    if method not in ['FBA', 'pFBA']:
        print('Method not available:', method)
        return

    if not solver:
        solver = solver_instance(ensemble.model)

    if get_fluxes:
        flux_sample = OrderedDict([(r_id, [None] * ensemble.size) for r_id in ensemble.model.reactions])
    else:
        objective = [None] * ensemble.size

    for i in range(ensemble.size):
        current = ensemble.get_constraints(i)

        if constraints:
            current.update(constraints)

        if method == 'FBA':
            sol = FBA(ensemble.model, constraints=current, solver=solver, get_values=get_fluxes)

        if method == 'pFBA':
            sol = pFBA(ensemble.model, constraints=current, solver=solver)

        if sol.status == Status.OPTIMAL:
            if get_fluxes:
                for r_id in ensemble.model.reactions:
                    flux_sample[r_id][i] = sol.values[r_id]
            else:
                objective[i] = sol.fobj

    if get_fluxes:
        return flux_sample
    else:
        return objective


def save_ensemble(ensemble, outputfile, **kwargs):
    """ Save ensemble model as an SBML file.

    Args:
        ensemble (EnsembleModel): model ensemble
        outputfile (str): output file
        **kwargs (dict): additional arguments to *save_cbmodel* method

    """

    for r_id, states in ensemble.reaction_states.items():
        state_as_str = ' '.join([str(int(x)) for x in states])
        ensemble.model.reactions[r_id].metadata['ENSEMBLE_STATE'] = state_as_str

    save_cbmodel(ensemble.model, outputfile, **kwargs)


def load_ensemble(inputfile, **kwargs):
    """ Load ensemble model from SBML file.

    Args:
        inputfile (str): input file
        **kwargs (dict): additional arguments to *load_cbmodel* method

    Returns:
        EnsembleModel: ensemble model

    """

    model = load_cbmodel(inputfile, **kwargs)
    reaction_states = {}

    for r_id, rxn in model.reactions.items():
        if 'ENSEMBLE_STATE' in rxn.metadata:
            state_as_str = rxn.metadata['ENSEMBLE_STATE']
            states = [bool(int(x)) for x in state_as_str.split()]
            reaction_states[r_id] = states

    sizes = list(map(len, reaction_states.values()))

    if len(set(sizes)) > 1:
        print('Error: reactions have different ensemble size')
        return

    return EnsembleModel(model, sizes[0], reaction_states)
