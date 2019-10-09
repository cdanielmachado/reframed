def escher_maps():
    """ List of all maps available in **escher**

    Returns:
        list: map names
    """

    try:
        import escher
    except ImportError:
        raise RuntimeError("Escher is not installed.")

    maps = escher.list_available_maps()
    return [entry['map_name'] for entry in maps]


def fluxes2escher(fluxes, map_name=None, fmt_func=None, **kwargs):
    """ Build escher map for a given flux distribution

    Args:
        fluxes (dict): flux distribution
        map_name (str): name of **escher** map (for a list of available maps run *escher_maps*)
        fmt_func (str or function): python format string (see Notes)
        kwargs: additional arguments passed to *escher.Builder* (see Escher's documentation for details).

    Notes:
        The format function parameter is used to convert reaction ids to BiGG ids.
        By default it removes the 'R_' prefix ('R_EX_h2o_e' is transformed to 'EX_h2o_e').

    Returns:
        escher.plots.Builder: escher map object
    """

    try:
        import escher
    except ImportError:
        raise RuntimeError("Escher is not installed.")

    if map_name is None:
        map_name = 'e_coli_core.Core metabolism'

    if fmt_func is None:
        def fmt_func(x):
            return x[2:]

    fluxes = {fmt_func(r_id): val for r_id, val in fluxes.items()}

    return escher.Builder(map_name=map_name, reaction_data=fluxes, **kwargs)

