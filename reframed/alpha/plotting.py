from math import isinf
import numpy as np


def plot_flux_bounds(range1, range2=None, values1=None, values2=None, reactions=None,
                     threshold=100, ax=None, figsize=None):
    """ Plot and compare flux ranges (although other types of ranges also supported).

    Args:
        range1 (dict): flux ranges
        range2 (dict): alternative ranges (optional)
        values1 (dict): unique values to plot inside first range (optional)
        values2 (dict): unique values to plot inside second range (optional)
        reactions (list): only display reactions in this list (optional)
        threshold (float): limit unbounded ranges (inf) to this value (default: 100)
        ax (matplotlib.Axes): plot over existing axes (optional)
        figsize (float, float): figure size if new axes is created (optional)

    """

    try:
        import matplotlib.pyplot as plt
    except ImportError:
        raise RuntimeError("Matplotlib is not installed.")

    if reactions is None or len(reactions) == 0:
        reactions = list(range1.keys())

    if ax is None:
        _, ax = plt.subplots(1, 1, figsize=figsize)

    def bounded_left(x):
        return -threshold if isinf(x) else x

    def bounded_right(x):
        return threshold if isinf(x) else x

    if range2 or values2:
        idx1 = np.arange(0.85, 2 * len(reactions), 2)
        idx2 = np.arange(0, 2 * len(reactions), 2)
    else:
        idx1 = np.arange(len(reactions))
        idx2 = None

    lb1 = np.array([bounded_left(range1[key][0]) for key in reactions])
    ub1 = np.array([bounded_right(range1[key][1]) for key in reactions])
    ax.barh(idx1, left=lb1, width=(ub1 - lb1), color='#7ba6ed', alpha=0.5)

    xmin = min(lb1)
    xmax = max(ub1)

    if range2:
        lb2 = np.array([bounded_left(range2[key][0]) for key in reactions])
        ub2 = np.array([bounded_right(range2[key][1]) for key in reactions])
        ax.barh(idx2, left=lb2, width=(ub2 - lb2), color='#43c6c2', alpha=0.5)

        xmin = min(xmin, min(lb2))
        xmax = max(xmax, max(ub2))

    if values1:
        v1 = np.array([values1[key] for key in reactions])
        ax.scatter(v1, idx1, color='#1527cf')

        xmin = min(xmin, min(v1))
        xmax = max(xmax, max(v1))

    if values2:
        v2 = np.array([values2[key] for key in reactions])
        ax.scatter(v2, idx2, color='#059c97')

        xmin = min(xmin, min(v2))
        xmax = max(xmax, max(v2))

    if idx2 is None:
        ax.set_yticks(idx1)
        ax.set_yticklabels(reactions)
        ax.set_ylim(-0.5, len(reactions)-0.5)
    else:
        ax.set_yticks(idx1-0.5)
        ax.set_yticklabels(reactions)
        ax.set_ylim(-1, 2 * len(reactions))

    ax.set_xlim(xmin * 1.05, xmax * 1.05)

    return ax
