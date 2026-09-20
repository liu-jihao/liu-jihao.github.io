"""Finite parameter domain for Section 16; no producer or private state imports."""
from fractions import Fraction as F

SCHEMA = "double-twentynine-v1"
PARAMETERS = {
    "lambda": "1731/2000", "mu": "23/25", "theta": "19/20",
    "terminal": "49/50", "square": "87/100", "old_price": "25/29",
    "generic_error": "1/10000", "terminal_error": "1/1000",
    "fixed_order_epsilon": "1/500", "scalar_margin": "1/2000",
    "terminal_margin": "1/100", "intrinsic_reserve": "1/100000",
}

GROUPS = [
    "d<=2", "d3e3", "d3e4", "d3e5", "d3e6m4", "d3e6m5to9",
    "d4e4", "d4e5", "d4e6m3", "d4e6m4", "d4e6m5", "d4e6m6", "d4e6m7",
]
SQUARE_GROUPS = GROUPS[:4] + ["d3e6"] + GROUPS[6:]


def group_for(v, square=False):
    d, e, m = v
    if d <= 2:
        return "d<=2"
    if d == 3 and e == 6:
        return "d3e6" if square else ("d3e6m4" if m == 4 else "d3e6m5to9")
    return f"d{d}e{e}" + (f"m{m}" if d == 4 and e == 6 else "")


def state_id(v):
    return "-".join(map(str, v))


def query_domain(states):
    """All finite comparisons used in the article, including zero caps.

    Each entry is (id, state, baseline, cap, source order, old-deficit cap).
    The data contain a value and an attaining flag for each such entry.
    This is the specification of the requested comparisons, not an evaluator.
    """
    lam, mu, theta = (F(PARAMETERS[k]) for k in ("lambda", "mu", "theta"))
    scenarios = [
        ("terminal", F(49, 50), F(4), {}),
        ("generic-lambda", lam, F(4), {}), ("generic-mu", mu, F(4), {}),
        ("square", F(87, 100), F(7, 2), {(4, 6, 3): F(7, 3), (3, 6, 4): F(5, 4)}),
        ("low-b", lam, F(3), {}),
    ]
    for family, overrides in [
        ("factorial", {(4, 6, 3): F(0), (4, 6, 5): F(6, 5)}),
        ("distinct", {(4, 6, 3): F(2), (4, 6, 5): F(1), (4, 6, 7): F(0)}),
        ("singular", {(4, 6, 3): F(2)}),
    ]:
        for name, gamma in [("lambda", lam), ("mu", mu)]:
            scenarios.append((family + "-" + name, gamma, F(4), overrides))
    for exponent in range(3, 8):
        for name, gamma in [("mu", mu), ("theta", theta)]:
            scenarios.append((f"smooth-{exponent}-{name}", gamma, F(4),
                              {(4, 6, 3): F(2) + F(1, exponent)}))
    result = []
    for name, gamma, bmax, overrides in scenarios:
        for v, row in states.items():
            cap = min(row["cap"], bmax, overrides.get(v, row["cap"]))
            result.append((name + "/" + state_id(v), v, 1/gamma, cap, gamma, bmax))
    for m in (3, 4, 5):
        v = (4, 6, m)
        result.append((f"affine/{m}", v, states[v]["q"], states[v]["cap"], None, None))
    result += [
        ("quartic-low", (3, 6, 4), 1/lam, F(4, 3), lam, F(4)),
        ("noncartier-quartic", (4, 6, 4), 1/lam, F(3, 2), lam, F(4)),
        ("degenerate-empty", (1, 1, 1), F(1), F(-1), None, None),
        ("degenerate-point", (1, 1, 1), F(1), F(0), None, None),
        ("degenerate-segment", (1, 1, 1), F(1), F(1), None, None),
    ]
    return result
