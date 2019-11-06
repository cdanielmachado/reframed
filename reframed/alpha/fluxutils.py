import re


def compare_fluxes(original, other, tolerance=1e-6, abstol=1e-9, sort=False, intersection=True, pattern=None):

    common = sorted(set(original.keys()) & set(other.keys()))

    if intersection:
        only_left = []
        only_right = []
    else:
        only_left = sorted(set(original.keys()) - set(other.keys()))
        only_right = sorted(set(other.keys()) - set(original.keys()))

    difference = [(r_id, abs(original[r_id] - other[r_id])) for r_id in common]
    flux_left = [(r_id, original[r_id]) for r_id in only_left]
    flux_right = [(r_id, other[r_id]) for r_id in only_right]

    if pattern is not None:
        re_expr = re.compile(pattern)
        difference = [x for x in difference if re_expr.search(x[0]) is not None]
        flux_left = [x for x in flux_left if re_expr.search(x[0]) is not None]
        flux_right = [x for x in flux_right if re_expr.search(x[0]) is not None]

    if sort:
        difference.sort(key=lambda x: x[1], reverse=True)
        flux_left.sort(key=lambda x: x[1], reverse=True)
        flux_right.sort(key=lambda x: x[1], reverse=True)

    for r_id, val in difference:
        if val > tolerance:
            x1 = original[r_id] if abs(original[r_id]) > abstol else 0
            x2 = other[r_id] if abs(other[r_id]) > abstol else 0
            print(f'{r_id: <16} {x1: < 10.3g} {x2: < 10.3g}')

    for r_id, val in flux_left:
        if abs(val) > tolerance:
            x = original[r_id] if abs(original[r_id]) > abstol else 0
            print(f'{r_id: <16} {x: < 10.3g}   --')

    for r_id, val in flux_right:
        if abs(val) > tolerance:
            x = other[r_id] if abs(other[r_id]) > abstol else 0
            print(f'{r_id: <16}   --       {x: < 10.3g}')