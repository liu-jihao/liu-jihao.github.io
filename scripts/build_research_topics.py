#!/usr/bin/env python3
"""Rebuild the collapsible "Research" topic lists on index.html from publications.html.

    python3 scripts/build_research_topics.py

Each topic is a <details> block (no JavaScript): the keyword is the <summary>, and opening it
shows Jihao's papers on that topic, newest first, taken verbatim from publications.html.
To file a new paper under topics, add its publication number to TOPICS_OF below and rerun.
Papers with no entry in TOPICS_OF (work outside birational geometry) are simply not listed here;
they stay on the Publications page.
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
    ("explicit", "Explicit birational geometry and boundedness"),
]

TOPICS_OF = {
    "72": "fano", "71": "mmp", "70": "gpairs mmp", "69": "fano", "64": "fano sing",
    "62": "fano", "61": "fano", "60a": "sing", "60b": "sing", "58": "mmp",
    "57": "sing explicit", "56": "fol explicit", "55": "sing", "54": "fano", "53": "fol explicit",
    "51": "sing", "50": "explicit", "49": "explicit fol", "48": "mmp gpairs", "47": "fano sing",
    "46": "fol mmp", "45": "fol", "44a": "mmp explicit", "44b": "mmp explicit", "43": "fol mmp",
    "42": "mmp", "41": "fol mmp", "40": "sing explicit", "39": "gpairs mmp", "38": "fol mmp",
    "37": "sing fano", "36": "fol mmp", "35": "fol explicit", "34": "fano sing explicit",
    "33": "fol mmp", "32": "mmp", "31": "fano explicit", "30": "fol gpairs mmp", "29": "explicit",
    "28": "fol sing", "27": "fol sing", "26": "explicit sing", "25": "gpairs", "24": "fol sing",
    "23": "fol sing", "22": "explicit mmp", "21": "gpairs mmp", "20": "mmp sing", "19": "sing",
    "18": "sing", "17": "mmp explicit", "16": "gpairs mmp", "15": "sing explicit", "14": "explicit",
    "13": "gpairs mmp", "12": "sing", "11": "gpairs mmp", "10": "sing explicit", "9": "sing",
    "8": "explicit", "7": "sing explicit", "6": "sing", "5": "sing", "4": "sing gpairs",
    "3": "sing", "2": "gpairs mmp", "1": "fol",
}


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
        out.append({"num": num, "title": title, "co": co, "venue": venue.rstrip(". "), "link": link})
    return out


def render(pubs):
    known = {p["num"] for p in pubs}
    missing = sorted(set(TOPICS_OF) - known)
    if missing:
        raise SystemExit("TOPICS_OF refers to publication numbers not on publications.html: %s" % missing)
    lines = ['      <p class="research-interests">',
             "        Birational geometry. Open a topic to see my papers on it.",
             "      </p>"]
    for key, label in TOPICS:
        items = [p for p in pubs if key in TOPICS_OF.get(p["num"], "").split()]
        lines.append('      <details class="topic">')
        lines.append('        <summary>%s <span class="topic-count">(%d papers)</span></summary>' % (label, len(items)))
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
