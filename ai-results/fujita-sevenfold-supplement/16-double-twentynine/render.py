"""Deterministic readable views of the canonical proof data; no optimization."""
from pathlib import Path
import argparse
import json


def tex(value):
    if "/" in value:
        a, b = value.split("/")
        return r"\frac{" + a + "}{" + b + "}"
    return value


def views(tables, arithmetic):
    out = {}
    lines = [r"\begin{center}", r"\small", r"\begin{tabular}{rrrrccc}",
             r"\(d_\nu\)&\(e_\nu\)&\(m_\nu\)&\(T_\nu\)&\(C_\nu\)&\(q_\nu\)&",
             r"\(J_\nu(50/49,C_\nu)\)\\ \hline"]
    for row in tables["states"]:
        cells = [str(row[k]) for k in ("d", "e", "m", "moment")]
        cells += [r"\(" + tex(row[k]) + r"\)" for k in ("cap", "q", "terminal_j")]
        lines.append("&".join(cells) + r"\\")
    lines += [r"\end{tabular}", r"\end{center}"]
    out["states.tex"] = "\n".join(lines) + "\n"
    for key, heading, columns in [
        ("generic", r"First state or group&\(\gamma=\lambda\)&\(\gamma=\mu\)\\ \hline", "lrr"),
        ("square", r"First state or group&Ceiling numerator\\ \hline", "lr"),
    ]:
        lines = [r"\begin{center}", r"\small", r"\begin{tabular}{" + columns + "}", heading]
        lines += [row["label"] + "&" + "&".join(map(str, row["numerators"])) + r"\\"
                  for row in tables[key]]
        lines += [r"\end{tabular}", r"\end{center}"]
        out[key + ".tex"] = "\n".join(lines) + "\n"
    lines = [r"\[", r"\begin{array}{c|rr}", r"\text{first state and cap}&\lambda&\mu\\ \hline"]
    for row in tables["replacements"]:
        lines.append(f"(4,6,{row['m']}),\\ {row['cap']}&" +
                     "&".join(map(str, row["numerators"])) + r"\\")
    lines += [r"\end{array}", r"\]"]
    out["replacements.tex"] = "\n".join(lines) + "\n"
    claims = arithmetic["claims"]
    states = {row["id"]: row for row in tables["states"]}
    lines = [r"\begin{array}{c|ccc}", r"m&3&4&5\\ \hline",
             "q_m&" + "&".join(tex(states[f"4-6-{m}"]["q"]) for m in (3, 4, 5)) + r"\\",
             "p_m&" + "&".join(tex(claims[f"affine_{m}"]) for m in (3, 4, 5)), r"\end{array}"]
    out["affine.tex"] = "\n".join(lines) + "\n"
    lines = [r"\begin{array}{c|ccc}", r"&m=3&m=4&m=5\\ \hline",
             r"\gamma=\lambda&" + "&".join(tex(claims[f"generic_order_{m}"]) for m in (3, 4, 5)),
             r"\end{array}", r"\quad>\frac{33}{100};\qquad",
             r"\gamma=\mu,\ m=3:\quad " + tex(claims["mu_order_3"]) + r">\frac{49}{100}."]
    out["generic-orders.tex"] = "\n".join(lines) + "\n"
    lines = [r"\begin{array}{c|ccc}", r"m&3&4&5\\ \hline",
             r"\gamma=\lambda&" + "&".join(tex(claims[f"power_order_{m}"]) for m in (3, 4, 5)) + r"\\",
             r"\text{strict lower bound}&33/100&7/20&17/50 .", r"\end{array}"]
    out["power-orders.tex"] = "\n".join(lines) + "\n"
    lines = [r"\begin{array}{c|c|c}", r"m_P&\eta\text{ used}&C_0-\Phi_\eta(1)-\Phi_\eta(m_P-1)\\ \hline"]
    for m, eta, key in [(3, "33/100", "cubic_reserve"), (4, "7/20", "singular_quartic_reserve"),
                        (5, "17/50", "singular_quintic_reserve")]:
        lines.append(str(m)+"&"+tex(eta)+"&"+tex(claims[key])+r"\\")
    lines.append(r"\end{array}")
    out["singular-reserves.tex"] = "\n".join(lines) + "\n"
    return out


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--data", type=Path, required=True)
    p.add_argument("--output", type=Path, required=True)
    args = p.parse_args()
    tables = json.loads((args.data / "tables.json").read_text())
    arithmetic = json.loads((args.data / "arithmetic.json").read_text())
    args.output.mkdir(parents=True, exist_ok=True)
    for name, text in views(tables, arithmetic).items():
        (args.output / name).write_text(text)
    print("Rendered eight views; mathematical acceptance requires check.py")


if __name__ == "__main__":
    main()
