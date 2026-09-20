#!/usr/bin/env python3
"""Rebuild the collapsible "Research" topic lists on index.html from publications.html.

    python3 scripts/build_research_topics.py

Each topic is a <details> block (no JavaScript): the keyword is the <summary>, and opening it
shows Jihao's papers on that topic, newest first, taken verbatim from publications.html, together with
the notes and the AI-generated papers listed in EXTRAS below.
To file a new paper under topics, add its publication number to TOPICS_OF below and rerun.
Algebraic geometry that fits no other topic goes under `generalag`; work outside algebraic geometry is
filed under its own field (`groups`, `comb`, ...; add a new field to TOPICS when needed); the Danus system
report has its own topic `danus`.
A paper with no entry in TOPICS_OF is not listed here; it stays on the Publications page.
"""
import html
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]

TOPICS = [
    ("mmp", "Minimal model program"),
    ("gpairs", "Generalized pairs"),
    ("fol", "Foliations"),
    ("sing", "Singularities: minimal log discrepancies, thresholds, ACC"),
    ("fano", "Fano varieties and K-stability"),
    ("cy", "Calabi–Yau varieties and log Calabi–Yau pairs"),
    ("bdd", "Boundedness"),
    ("explicit", "Explicit birational geometry"),
    ("surf", "Surfaces"),
    ("generalag", "General algebraic geometry"),
    ("groups", "Group theory and operator algebras"),
    ("comb", "Combinatorics and number theory"),
    ("danus", "Danus"),
]

TOPICS_OF = {
    "72": "fano", "71": "surf", "70": "gpairs mmp", "69": "fano", "64": "fano sing",
    "62": "fano", "61": "fano", "60a": "sing", "60b": "sing", "58": "surf",
    "57": "sing bdd", "56": "fol bdd", "55": "sing", "54": "cy", "53": "fol explicit",
    "51": "sing", "50": "explicit surf", "49": "bdd fol", "48": "mmp gpairs", "47": "fano sing cy",
    "46": "fol mmp", "45": "fol", "44a": "mmp bdd explicit", "44b": "mmp explicit", "43": "fol mmp",
    "42": "mmp", "41": "fol mmp bdd", "40": "sing explicit", "39": "gpairs mmp surf", "38": "fol mmp",
    "37": "sing fano", "36": "fol mmp", "35": "fol bdd", "34": "fano sing explicit",
    "33": "fol mmp", "32": "mmp", "31": "fano explicit", "30": "fol gpairs mmp", "29": "explicit surf",
    "28": "fol sing", "27": "fol sing", "26": "explicit sing surf", "25": "gpairs", "24": "fol sing surf",
    "23": "fol sing", "22": "bdd mmp", "21": "gpairs mmp", "20": "mmp sing", "19": "sing",
    "18": "sing bdd surf", "17": "mmp bdd", "16": "gpairs mmp", "15": "sing explicit", "14": "explicit surf",
    "13": "gpairs mmp", "12": "sing", "11": "gpairs mmp", "10": "sing explicit surf", "9": "sing surf",
    "8": "bdd", "7": "sing explicit", "6": "sing bdd", "5": "sing", "4": "sing gpairs",
    "3": "sing", "2": "gpairs mmp", "1": "fol surf",
    "66": "danus",
    "52": "surf", "68": "generalag", "67": "generalag", "65": "cy", "63": "generalag", "59": "comb",
}

# Items that are not on publications.html: notes (notes.html) and AI-generated papers that no human
# has verified yet (ai-results.html). The Ehrhart note is publication 69 and is not repeated here.
# (title, link, label shown after the title, arXiv-style yymm used only for sorting, topics)
EXTRAS = [
    ("K-polystability of nilpotent orbit and Slodowy slice bases",
     "ai-results/k-polystability-of-nilpotent-orbit-and-slodowy-slice-bases-2026-09-20.pdf",
     "AI-generated, not verified by a human", 2609, "fano"),
    ("Singular points of rank-one del Pezzo surfaces in characteristic three",
     "notes/del-pezzo-characteristic-three-2026-09-20.pdf",
     "Note", 2609, "surf sing"),
    ("Nonhyperlinear groups exist", "ai-results/nonhyperlinear-groups-exist-2026-09-20.pdf",
     "AI-generated, not verified by a human", 2609, "groups"),
    ("Fujita freeness on complex projective sevenfolds", "ai-results/fujita-freeness-on-sevenfolds-2026-09-20.pdf",
     "AI-generated, not verified by a human", 2609, "explicit"),
    ("Fujita freeness in dimension six", "ai-results/fujita-freeness-in-dimension-six-2026-09-20.pdf",
     "AI-generated, not verified by a human", 2609, "explicit"),
    ("Syzygy bundles on Picard rank one varieties need not be stable", "notes/SyzygyBundleCounterexample.pdf",
     "Note; AI-generated, not verified by a human", 2607, "generalag"),
    ("Factorial asymptotics of the Matryoshka numbers", "notes/Matryoshka.pdf", "Note", 2606, "comb"),
    ("On a conjecture of Esser, Totaro, and Wang", "notes/ET74.pdf", "Note", 2605, "explicit cy"),
]


def parse_publications():
    s = (ROOT / "publications.html").read_text(encoding="utf-8")
    out = []
    for b in re.findall(r'      <div class="entry">\n.*?\n      </div>\n', s, re.S):
        m = re.search(r'entry-number">([^<]*)<', b)
        if not m:
            continue
        num = m.group(1).strip(". ")
        title = re.search(r'entry-title">(.*?)</div>', b, re.S).group(1)
        title = re.sub(r'<span class="entry-number">.*?</span>\s*', "", title).strip()
        am = re.search(r'entry-authors">(.*?)</div>', b, re.S)
        authors = [a.strip() for a in re.sub("<[^>]+>", "", am.group(1)).split(",")] if am else []
        co = [a for a in authors if html.unescape(a) != "Jihao Liu" and a]
        vm = re.search(r'entry-venue">(.*?)</div>', b, re.S)
        venue = re.sub("<[^>]+>", "", vm.group(1)).strip() if vm else ""
        venue = re.sub(r"\s*arXiv:\S+.*$", "", venue).strip()
        if venue.startswith("Note") or venue.lower().startswith("arxiv"):
            venue = ""
        links = re.findall(r'href="([^"]+)"', b)
        link = next((l for l in links if "arxiv.org" in l), links[0] if links else "")
        ym = re.search(r"arxiv\.org/abs/(\d{4})\.", link)
        out.append({"num": num, "title": title, "co": co, "venue": venue.rstrip(". "), "link": link,
                    "yymm": int(ym.group(1)) if ym else None})
    last = 9999
    for q in out:                       # publications are listed newest first
        q["yymm"] = last = q["yymm"] if q["yymm"] is not None else last
    return out


def render(pubs):
    known = {p["num"] for p in pubs}
    missing = sorted(set(TOPICS_OF) - known)
    if missing:
        raise SystemExit("TOPICS_OF refers to publication numbers not on publications.html: %s" % missing)
    lines = ['      <p class="research-interests">',
             "        Open a topic to see my papers on it.",
             "      </p>"]
    for key, label in TOPICS:
        items = [(p["yymm"], 0, -i, p) for i, p in enumerate(pubs) if key in TOPICS_OF.get(p["num"], "").split()]
        for j, (title, link, tag, yymm, topics) in enumerate(EXTRAS):
            if key in topics.split():
                if not (ROOT / link).exists():
                    raise SystemExit("EXTRAS link does not exist: %s" % link)
                items.append((yymm, 1, -j, {"title": title, "link": link, "co": [], "venue": tag}))
        items = [it[3] for it in sorted(items, key=lambda t: t[:3], reverse=True)]
        lines.append('      <details class="topic">')
        count = "%d paper%s" % (len(items), "" if len(items) == 1 else "s")
        lines.append('        <summary>%s <span class="topic-count">(%s)</span></summary>' % (label, count))
        lines.append('        <ul class="topic-list">')
        for p in items:
            t = '<a href="%s">%s</a>' % (p["link"], p["title"]) if p["link"] else p["title"]
            co = " (with %s)" % ", ".join(p["co"]) if p["co"] else ""
            ven = ' <span class="topic-venue">%s.</span>' % p["venue"] if p["venue"] else ""
            lines.append("          <li>%s%s.%s</li>" % (t, co, ven))
        lines.append("        </ul>")
        lines.append("      </details>")
    return "\n".join(lines) + "\n"


def main():
    pubs = parse_publications()
    idx = ROOT / "index.html"
    s = idx.read_text(encoding="utf-8")
    a, b = "      <!-- RESEARCH:START -->\n", "      <!-- RESEARCH:END -->\n"
    if a not in s or b not in s:
        raise SystemExit("index.html lacks the RESEARCH:START / RESEARCH:END markers")
    s = s[: s.index(a) + len(a)] + render(pubs) + s[s.index(b):]
    idx.write_text(s, encoding="utf-8")
    unfiled = [p["num"] for p in pubs if p["num"] not in TOPICS_OF]
    print("topics rebuilt; papers not filed under any topic:", ", ".join(unfiled))


if __name__ == "__main__":
    main()
