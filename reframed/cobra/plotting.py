from .variability import FVA
from .simulation import FBA
import numpy as np


def flux_envelope(model, r_x, r_y, steps=10, constraints=None):
    """ Calculate the flux envelope for a pair of reactions.

    Arguments:
        model (CBModel): the model
        r_x (str): reaction on x-axis
        r_y (str): reaction on y-axis
        steps (int): number of steps to compute (default: 10)
        constraints (dict): custom constraints to the FBA problem

    Returns:
        tuple: x values, y_min values, y_max values
    """

    x_range = FVA(model, reactions=[r_x], constraints=constraints)
    xmin, xmax = x_range[r_x]
    xvals = np.linspace(xmin, xmax, steps)
    ymins, ymaxs = np.zeros(steps), np.zeros(steps)

    if constraints is None:
        _constraints = {}
    else:
        _constraints = {}
        _constraints.update(constraints)

    for i, xval in enumerate(xvals):
        _constraints[r_x] = xval
        y_range = FVA(model, reactions=[r_y], constraints=_constraints)
        ymins[i], ymaxs[i] = y_range[r_y]

    return xvals, ymins, ymaxs


def plot_flux_envelope(model, r_x, r_y, steps=10, substrate=None, constraints=None,
                       label_x=None, label_y=None, flip_x=False, flip_y=False,
                       plot_kwargs=None, fill_kwargs=None, ax=None):
    """ Plots the flux envelope for a pair of reactions.

    Arguments:
        model (CBModel): the model
        r_x (str): reaction on x-axis
        r_y (str): reaction on y-axis
        steps (int): number of steps to compute (default: 20)
        substrate (str): compute yields for given substrate instead of rates (optional)
        constraints (dict): additional simulation constraints
        label_x (str): x label (optional, uses reaction name by default)
        label_y (str): y label (optional, uses reaction name by default)
        flip_x (bool): flip direction of r_x (default: False)
        flip_y (bool): flip direction of r_y (default: False)
        plot_kwargs (dict): additional parameters to *pyplot.plot* (optional)
        fill_kwargs (dict): additional parameters to *pyplot.fill_between* (optional)
        ax (matplotlib.Axes): plot over existing axes (optional)

    Returns:
        matplotlib.Axes: axes object
    """

    try:
        import matplotlib.pyplot as plt
    except ImportError:
        raise RuntimeError("Matplotlib is not installed.")

    offset = 0.03

    if ax is None:
        _, ax = plt.subplots()

    if not plot_kwargs:
        plot_kwargs = {'color': 'k'}

    if not fill_kwargs:
        fill_kwargs = {'color': 'k', 'alpha': 0.1}

    xvals, ymins, ymaxs = flux_envelope(model, r_x, r_y, steps, constraints)

    if flip_x:
        xvals, ymins, ymaxs = -xvals, ymins[::-1], ymaxs[::-1]

    if flip_y:
        ymins, ymaxs = -ymaxs, -ymins

    if substrate:
        sol = FBA(model)
        uptk = abs(sol.values[substrate])
        xvals, ymins, ymaxs = xvals / uptk, ymins / uptk, ymaxs / uptk

    ax.plot(xvals, ymins, **plot_kwargs)
    ax.plot(xvals, ymaxs, **plot_kwargs)
    ax.plot([xvals[0], xvals[0]], [ymins[0], ymaxs[0]], **plot_kwargs)
    ax.plot([xvals[-1], xvals[-1]], [ymins[-1], ymaxs[-1]], **plot_kwargs)

    ax.fill_between(xvals, ymins, ymaxs, **fill_kwargs)

    ax.set_xlabel(label_x) if label_x else ax.set_xlabel(model.reactions[r_x].name)
    ax.set_ylabel(label_y) if label_y else ax.set_ylabel(model.reactions[r_y].name)

    xmin, xmax = min(xvals), max(xvals)
    dx = offset * (xmax - xmin)
    ax.set_xlim((xmin - dx, xmax + dx))

    ymin, ymax = min(ymins), max(ymaxs)
    dy = offset * (ymax - ymin)
    ax.set_ylim((ymin - dy, ymax + dy))

    return ax

