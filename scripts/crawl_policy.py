#!/usr/bin/env python3
"""Maintain voluntary crawler exclusions; this is not access control."""
from pathlib import Path
import argparse
import os
import re

ROOT = Path(__file__).resolve().parents[1]
ROBOTS_TXT = "User-agent: *\nDisallow: /\n"
ROBOTS_META = '  <meta name="robots" content="noindex, nofollow, nosnippet, noimageindex" />'


def apply_policy(source):
    """Normalize the head's robots directive without changing page content."""
    match = re.search(r'(<head\b[^>]*>)(.*?)(</head\s*>)', source, re.I | re.S)
    if not match:
        raise ValueError('HTML page has no head element')
    head = re.sub(
        r'^[ \t]*<meta\b(?=[^>]*\bname\s*=\s*[\"\']robots[\"\'])[^>]*>[ \t]*\n?',
        '', match.group(2), flags=re.I | re.M,
    )
    charset = re.search(r'<meta\b[^>]*\bcharset\s*=[^>]*>[^\S\n]*\n?', head, re.I)
    leading_newline = re.match(r'[ \t]*\n', head)
    at = charset.end() if charset else (leading_newline.end() if leading_newline else 0)
    separator = '' if head[:at].endswith('\n') else '\n'
    head = head[:at] + separator + ROBOTS_META + '\n' + head[at:]
    return source[:match.start(2)] + head + source[match.end(2):]


def html_pages():
    # Include source templates as well as root and nested published pages.
    # Skip hidden build/virtualenv directories and non-HTML source directories.
    for directory, names, files in os.walk(ROOT):
        names[:] = sorted(name for name in names if not name.startswith('.'))
        if Path(directory) == ROOT:
            names[:] = [name for name in names if name not in ('scripts', 'data')]
        for name in sorted(files):
            if name.endswith('.html'):
                yield Path(directory) / name


def build(check=False):
    outputs = {ROOT / 'robots.txt': ROBOTS_TXT}
    outputs.update((path, apply_policy(path.read_text())) for path in html_pages())
    stale = [path for path, text in outputs.items() if not path.is_file() or path.read_text() != text]
    if check and stale:
        raise SystemExit('Crawler policy is stale: ' + ', '.join(str(p.relative_to(ROOT)) for p in stale))
    if not check:
        for path in stale:
            path.write_text(outputs[path])
    print(f'Crawler policy checked: robots.txt and {len(outputs) - 1} HTML files' +
          ('.' if check else f'; {len(stale)} files updated.'))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--check', action='store_true')
    build(parser.parse_args().check)
