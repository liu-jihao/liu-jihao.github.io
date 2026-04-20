#!/usr/bin/env python3
"""Inject favicon links and Open Graph / Twitter card meta tags into every page.

Idempotent: if the marker block is already present, the file's existing block
is replaced in place. If not, it is inserted right after the description meta
tag.

Run:  .venv/bin/python3 scripts/inject_head_tags.py
"""

from __future__ import annotations

import re
from html import escape
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
SITE_URL = "https://jihaoliu.org"
SITE_NAME = "Jihao Liu"

START = "<!-- META:START -->"
END = "<!-- META:END -->"

# (filename, og_title, og_description)
PAGES: list[tuple[str, str, str]] = [
    ("index.html",
     "Jihao Liu — Peking University",
     "Algebraic geometer at Peking University (BICMR). Birational geometry, singularities, foliations, and explicit birational geometry."),
    ("publications.html",
     "Publications — Jihao Liu",
     "Publications and preprints of Jihao Liu in algebraic geometry, birational geometry, and foliations."),
    ("teaching.html",
     "Teaching — Jihao Liu",
     "Teaching at Peking University, Northwestern University, and the University of Utah."),
    ("talks.html",
     "Invited Talks — Jihao Liu",
     "Conference, workshop, and seminar talks given by Jihao Liu."),
    ("conferences.html",
     "Conferences and Seminars — Jihao Liu",
     "Conferences and seminars co-organized by Jihao Liu."),
    ("collaborators.html",
     "Collaborators — Jihao Liu",
     "Mathematical collaborators of Jihao Liu, aggregated from his publications and preprints."),
    ("pku-ag-seminar.html",
     "PKU Algebraic Geometry Seminar — Spring 2026",
     "Schedule of the Peking University Algebraic Geometry Seminar, Spring 2026."),
    ("pku-ag-seminar-2024-fall.html",
     "PKU Algebraic Geometry Seminar — Fall 2024",
     "Archive: Peking University Algebraic Geometry Seminar, Fall 2024."),
    ("pku-ag-seminar-2025-spring.html",
     "PKU Algebraic Geometry Seminar — Spring 2025",
     "Archive: Peking University Algebraic Geometry Seminar, Spring 2025."),
    ("ffm-conference.html",
     "Conference on Foliation, Fibration, and Moduli — Westlake, August 2026",
     "International conference at Westlake University on the birational geometry of foliations, fibrations, and their moduli. August 24–29, 2026, Hangzhou."),
    ("404.html",
     "Page not found — Jihao Liu",
     "The page you're looking for doesn't exist on jihaoliu.org."),
]


def block_for(page: str, title: str, description: str) -> str:
    page_url = f"{SITE_URL}/{page}" if page != "index.html" else f"{SITE_URL}/"
    img_url = f"{SITE_URL}/assets/og-image.png"
    t = escape(title, quote=True)
    d = escape(description, quote=True)
    return "\n".join([
        f"  {START}",
        '  <link rel="icon" type="image/svg+xml" href="/assets/favicon.svg" />',
        '  <link rel="icon" type="image/png" sizes="32x32" href="/assets/favicon-32.png" />',
        '  <link rel="icon" type="image/x-icon" href="/assets/favicon.ico" />',
        '  <link rel="apple-touch-icon" sizes="180x180" href="/assets/favicon-180.png" />',
        f'  <link rel="canonical" href="{page_url}" />',
        '',
        '  <meta property="og:type" content="website" />',
        f'  <meta property="og:site_name" content="{SITE_NAME}" />',
        f'  <meta property="og:title" content="{t}" />',
        f'  <meta property="og:description" content="{d}" />',
        f'  <meta property="og:url" content="{page_url}" />',
        f'  <meta property="og:image" content="{img_url}" />',
        '  <meta property="og:image:width" content="1200" />',
        '  <meta property="og:image:height" content="630" />',
        '',
        '  <meta name="twitter:card" content="summary_large_image" />',
        f'  <meta name="twitter:title" content="{t}" />',
        f'  <meta name="twitter:description" content="{d}" />',
        f'  <meta name="twitter:image" content="{img_url}" />',
        f"  {END}",
    ])


_existing_block = re.compile(
    rf"^[ \t]*{re.escape(START)}.*?{re.escape(END)}\n",
    re.DOTALL | re.MULTILINE,
)
_description_line = re.compile(
    r'^([ \t]*)<meta name="description"[^>]*/>\n',
    re.MULTILINE,
)


def update_file(path: Path, title: str, description: str) -> str:
    text = path.read_text(encoding="utf-8")
    block = block_for(path.name, title, description) + "\n"

    if _existing_block.search(text):
        new_text, n = _existing_block.subn(block, text, count=1)
        action = "updated"
    else:
        m = _description_line.search(text)
        if not m:
            return f"SKIPPED (no <meta description> anchor): {path.name}"
        insert_at = m.end()
        new_text = text[:insert_at] + block + text[insert_at:]
        n = 1
        action = "inserted"

    if n != 1:
        return f"SKIPPED (regex matched {n} times): {path.name}"

    if new_text != text:
        path.write_text(new_text, encoding="utf-8")
        return f"{action}: {path.name}"
    return f"unchanged: {path.name}"


def main() -> None:
    for fname, title, desc in PAGES:
        p = ROOT / fname
        if not p.exists():
            print(f"missing: {fname}")
            continue
        print(update_file(p, title, desc))


if __name__ == "__main__":
    main()
