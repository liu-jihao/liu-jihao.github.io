#!/usr/bin/env python3
"""Build the liuxiv pages from data/liuxiv.json (stdlib only; local; never publishes).

Every AI result and Note has a permanent identifier liuxiv:YYMM.NNNN (data/liuxiv.json). This script
writes, from that file and data/navigation.json only:

  liuxiv.html                  the liuxiv listing page, arXiv style: one heading per day of first posting, entries
                               numbered from the top (identifier, title, authors, first line of the Comments,
                               Subjects, machine-verified, human-verified, links);
  abs/<YYMM.NNNN>/index.html   one arXiv-style abstract page per record (abstract, Comments, Subjects,
                               Machine-verified, Human-verified, Links, Cite as, Listed on, history);
  papers/<slug>.html           a static redirect (meta refresh, rel=canonical, visible link) from every
                               former permanent page to its abstract page.

ai-results.html and notes.html are not written: they keep their entries, numbers and anchors, and the
abstract pages link back to them. PDFs are never changed; pdf/<id>.pdf and pdf/<id>v<k>.pdf are byte-identical copies. liuxiv.html is passed
through build_chrome.render_page, so the shared navigation leaves it unchanged. Idempotent.

  python3 scripts/build_liuxiv.py           write the pages
  python3 scripts/build_liuxiv.py --check   report stale pages, write nothing
"""
from pathlib import Path
from html import escape, unescape
import argparse, datetime, json, re, shutil, sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'scripts'))
from build_chrome import render_page  # noqa: E402

DATA = ROOT / 'data/liuxiv.json'
NAV = ROOT / 'data/navigation.json'
SITE = 'https://jihaoliu.org'
HOST = 'jihaoliu.org'
PREFIX = 'liuxiv:'
PAGE = 'liuxiv.html'
ID_RE = re.compile(r'^\d{4}\.\d{4}$')
PDF_RE = re.compile(r'^\d{4}\.\d{4}(?:v\d+)?\.pdf$')  # the short PDF addresses pdf/<id>.pdf and pdf/<id>v<k>.pdf
MARK = {'yes': ('✓', 'lx-yes', 'yes'), 'no': ('✗', 'lx-no', 'no'), None: ('–', 'lx-none', 'not stated')}


def plain(html):
    return re.sub(r'\s+', ' ', unescape(re.sub(r'<[^>]+>', '', html or ''))).strip()


def root_relative(html):
    """Links of the site entries are relative to the site root; under /abs/ they need a leading '/'."""
    def fix(m):
        url = m.group(2)
        if not url or re.match(r'(?:[a-zA-Z][a-zA-Z0-9+.-]*:|/|#)', url):
            return m.group(0)
        return f'{m.group(1)}="/{url}"'
    return re.sub(r'\b(href|src)="([^"]*)"', fix, html)


def retarget(html, permanent, prefix):
    """Point links at former permanent pages papers/<slug>.html to abstract pages; text unchanged."""
    def fix(m):
        path = m.group(2).split('#', 1)[0].lstrip('/')
        if path in permanent:
            return f'{m.group(1)}="{prefix}abs/{permanent[path]}/"'
        return m.group(0)
    return re.sub(r'\b(href)="([^"]*)"', fix, html)


def link_systems(html, links):
    """Jihao 11-11: every mention of Danus, Rethlas or Pharos links to the system's paper; text inside links is left alone."""
    if not links:
        return html
    out, depth = [], 0
    for part in re.split(r'(<[^>]+>)', html):
        if part.startswith('<'):
            if re.match(r'<a\b', part):
                depth += 1
            elif part.startswith('</a'):
                depth = max(0, depth - 1)
            out.append(part)
        elif depth == 0:
            # Jihao 11-12: Danus_v1, Danus_v3 (and any _vN) link to the system's page as a whole token
            out.append(re.sub(r'(?<![\w])(' + '|'.join(map(re.escape, links)) + r')(_v\d+)?(?![\w])',
                              lambda m: f'<a href="{links[m.group(1)]}">{m.group(0)}</a>', part))
        else:
            out.append(part)
    return ''.join(out)


def mark(value, hidden=True):
    """The check mark, cross or dash; with hidden=True a visually hidden word for screen readers."""
    sym, c, word = MARK[value]
    extra = f'<span class="visually-hidden"> {word}</span>' if hidden else ''
    return f'<span class="lx-mark {c}" title="{word}">{sym}{extra}</span>'


# ---------------------------------------------------------------- shared chrome (pages under /abs/ and /papers/)

def header(nav, active):
    def link(pair):
        href, label = pair
        cls = ' class="active"' if href == active else ''
        return f'          <a href="/{escape(href, quote=True)}"{cls}>{escape(label)}</a>'
    return '\n'.join([
        '  <a class="skip-link" href="#main-content">Skip to content</a>',
        '  <div class="container">',
        '    <header class="site-header">',
        '      <div class="site-title"><a href="/index.html">Jihao Liu&nbsp;<span style="font-weight: 400; color: var(--color-text-muted);">(刘济豪)</span></a></div>',
        '      <div class="site-subtitle">Peking University · School of Mathematical Sciences</div>',
        '      <nav class="main-nav" aria-label="Main navigation">',
        '        <div class="nav-primary">',
        *map(link, nav['primary']),
        '        </div>',
        '        <details class="nav-more">',
        '          <summary>Research &amp; community</summary>',
        '          <div class="nav-secondary">',
        *map(link, nav['secondary']),
        '          </div>',
        '        </details>',
        '      </nav>',
        '    </header>',
    ])


FOOTER = '''    <footer class="site-footer">
      <p><a href="#main-content">Back to top</a> &middot; <a href="{root}cv/cv.pdf">Curriculum vitae (PDF)</a></p>
      <p>&copy; 2026 Jihao Liu.</p>
      <p class="site-notice" lang="zh-Hans">本主页全部内容未经本人书面许可，禁止任何形式转载、复制、截图转发、摘录二次发布。任何公共平台转载、摘抄均属侵权。</p>
      <p class="site-notice">All content on this website may not be reproduced, copied, forwarded as screenshots, or excerpted and republished in any form without my written permission. Any reposting or excerpting on a public platform constitutes an infringement.</p>
    </footer>
  </div>
</body>
</html>
'''

ICONS = '''  <link rel="icon" type="image/svg+xml" href="/assets/favicon.svg" />
  <link rel="icon" type="image/png" sizes="32x32" href="/assets/favicon-32.png" />
  <link rel="icon" type="image/x-icon" href="/assets/favicon.ico" />
  <link rel="apple-touch-icon" sizes="180x180" href="/assets/favicon-180.png" />'''


def head(title, desc, url, og_type, root):
    t, d = escape(title, quote=True), escape(desc, quote=True)
    return f'''<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="UTF-8" />
  <meta name="viewport" content="width=device-width, initial-scale=1.0" />
  <title>{escape(title, quote=False)}</title>
  <meta name="description" content="{d}" />
  <!-- META:START -->
{ICONS}
  <link rel="canonical" href="{url}" />

  <meta property="og:type" content="{og_type}" />
  <meta property="og:site_name" content="Jihao Liu" />
  <meta property="og:title" content="{t}" />
  <meta property="og:description" content="{d}" />
  <meta property="og:url" content="{url}" />
  <meta property="og:image" content="{SITE}/assets/og-image.png" />
  <meta property="og:image:width" content="1200" />
  <meta property="og:image:height" content="630" />

  <meta name="twitter:card" content="summary_large_image" />
  <meta name="twitter:title" content="{t}" />
  <meta name="twitter:description" content="{d}" />
  <meta name="twitter:image" content="{SITE}/assets/og-image.png" />
  <!-- META:END -->
  <link rel="stylesheet" href="{root}assets/fonts.css" />
  <link rel="stylesheet" href="{root}assets/style.css" />
</head>
<body>
'''


# ---------------------------------------------------------------- data

def load():
    data = json.loads(DATA.read_text())
    nav = json.loads(NAV.read_text())
    if data.get('schema') != 'liuxiv-test-2':
        raise SystemExit('data/liuxiv.json: expected schema liuxiv-test-2')
    if not any(href == PAGE for href, _ in nav['primary'] + nav['secondary']):
        raise SystemExit('data/navigation.json has no link to ' + PAGE)
    records = data['records']
    ids = [r['id'] for r in records]
    if len(ids) != len(set(ids)):
        raise SystemExit('duplicate ids in data/liuxiv.json')
    for r in records:
        if not ID_RE.match(r['id']):
            raise SystemExit('malformed id ' + r['id'])
        for k in ('abstract_display', 'comments_html', 'machine_verified', 'human_verified', 'links', 'history_display', 'listing'):
            if k not in r:
                raise SystemExit(f'{r["id"]}: field {k} missing')
        h = r['history_display']
        if h and h[-1]['pdf_path'] != r['old_urls']['pdf_current']:
            raise SystemExit('current PDF is not the last version shown for ' + r['id'])
    permanent = {r['old_urls']['permanent_page']: r['id'] for r in records if r['old_urls']['permanent_page']}
    return data, nav, records, permanent


def subjects(r):
    if r['old']['section']:
        return [r['old']['section']]
    t = r.get('topics') or {}
    if t.get('labels'):
        return t['labels']
    pub = r.get('publication') or {}
    return pub.get('topic_labels') or []


def subjects_html(r, cls):
    """Jihao 11-11: the main subject in bold, the short sub-topics in parentheses, e.g. Birational geometry (Explicit birational geometry, Calabi–Yau varieties)."""
    sd = r.get('subjects_display')
    if not sd:
        subj = subjects(r)
        return f'<span class="{cls}">{escape(subj[0])}</span>' + ''.join('; ' + escape(x) for x in subj[1:]) if subj else ''
    out = f'<span class="{cls}">{escape(sd["primary"])}</span>'
    if sd['sub']:
        out += ' (' + ', '.join(escape(x) for x in sd['sub']) + ')'
    for a in sd.get('also') or []:          # Jihao 11-12: further subjects after a semicolon, e.g. Dynamical systems; Foliations
        out += '; ' + escape(a['name'])
        if a.get('sub'):
            out += ' (' + ', '.join(escape(x) for x in a['sub']) + ')'
    return out


def page_title(r):
    pp = r.get('permanent_page') or {}
    if pp.get('title_tag'):
        return unescape(pp['title_tag'])            # e.g. "Shokurov's ... — Jihao Liu"
    return plain(r['title_html']) + ' — Jihao Liu'


def description(r):
    pp = r.get('permanent_page') or {}
    if pp.get('meta_description'):
        return unescape(pp['meta_description'])     # the attribute value of the former permanent page
    text = plain(r['abstract_display']['html'][0]) if r['abstract_display']['html'] else plain(r['title_html'])
    if len(text) > 300:
        text = text[:300].rsplit(' ', 1)[0].rstrip(',;:') + '…'
    return text


def links_html(r, with_notes):
    out = []
    for x in r['links']:
        a = f'<a href="{escape(x["href"], quote=True)}">{escape(x["label"])}</a>'
        if with_notes and x.get('note'):
            a += f' ({escape(x["note"])})'
        out.append(a)
    return out


def access_line(r):
    pp = r.get('permanent_page') or {}
    if r.get('access_links_html'):   # 11-12: same form as the other papers with ancillary files
        line = r['access_links_html']
    elif pp.get('links_html'):
        line = pp['links_html']  # root-relative already: PDF (N pages), ancillary files and README
    else:
        line = f'<a href="/{escape(r["old_urls"]["pdf_current"], quote=True)}">PDF</a> ({r["pages_current"]} pages)'
    for x in r['extra_links']:
        if x['href'].startswith('http') or x['href'] in line or '/' + x['href'] in line:
            continue
        line += f' &middot; <a href="{escape(x["href"], quote=True)}">{x["label"]}</a>'
    zips = [a for a in r['ancillary'] if a['path'].endswith('.zip') and a.get('sha256') and a['path'] in line]
    for a in zips:  # the site prints the SHA-256 of the zip archive of the ancillary files
        line, n = re.subn(rf'(href="/{re.escape(a["path"])}">Ancillary files</a> \(zip, [0-9.]+&nbsp;MB)',
                          lambda m: m.group(1) + f', SHA-256 <span class="abs-hash">{a["sha256"]}</span>', line, count=1)
        if n != 1:
            raise SystemExit(f'{r["id"]}: cannot place the SHA-256 of {a["path"]}')
    line = root_relative(line)
    old = f'href="/{escape(r["old_urls"]["pdf_current"], quote=True)}"'
    if old not in line:
        raise SystemExit(f'{r["id"]}: the access line has no link to the current PDF')
    return line.replace(old, f'href="/pdf/{r["id"]}.pdf"')


def dateline(r):
    h = r['history_display']
    if not h:
        return None
    if len(h) == 1:
        return f'[Posted on {escape(h[0]["time"])}]'
    return (f'[Posted on {escape(h[0]["time"])} (<a href="/abs/{r["id"]}v1/">v1</a>), '
            f'last revised {escape(h[-1]["time"])} (this version, v{len(h)})]')


def history_items(r, this=None):
    """Submission history as on arXiv: [vk] opens the page of version k."""
    items, n = [], len(r['history_display'])
    this = this or n
    for v in r['history_display']:
        bits = [f'<strong><a href="/abs/{r["id"]}v{v["v"]}/">[v{v["v"]}]</a></strong> {escape(v["time"])}']
        if v['change']:
            bits.append(escape(v['change']))
        bits.append(f'{v["pages"]} pages')
        if v['sha256']:
            bits.append(f'SHA-256 <span class="abs-hash">{v["sha256"]}</span>')
        if v.get('timestamp_html'):          # Jihao 11-12 (Comments routing): a public timestamp of this version sits by its hash
            bits.append(v['timestamp_html'])
        cls = ' class="abs-current"' if v['v'] == this else ''
        items.append(f'        <li{cls}>' + ' &middot; '.join(bits) + '</li>')
    return '\n'.join(items)


# ---------------------------------------------------------------- abstract pages

def abs_page(data, nav, r, permanent):
    rid = r['id']
    url = f'{SITE}/abs/{rid}/'
    tx = lambda h: root_relative(retarget(h, permanent, '/'))
    title = f'[{PREFIX}{rid}] {page_title(r)}'
    lst = r['listing']

    body = [
        '    <main id="main-content" class="lx-page lx-abs" tabindex="-1">',
        '      <div class="abs-header">',
        f'        <span class="abs-id">{PREFIX}{rid}</span>',
        f'        <span class="abs-browse"><a href="/{PAGE}">liuxiv</a></span>',
        '      </div>',
    ]
    if r.get('first_version_ai'):
        body.append(f'      <div class="lx-ainote">{data["liuxiv_page"]["ai_note_page_html"]}</div>')
    dl = dateline(r)
    if dl:
        body.append(f'      <p class="abs-dateline">{dl}</p>')
    body.append(f'      <h1 class="abs-title">{r["title_html"]}</h1>')
    if r['authors_html']:
        body.append(f'      <div class="entry-authors abs-authors">{tx(r["authors_html"])}</div>')
        if r.get('authorship_note_html'):   # Jihao 11-12: the authorship note (collaboration welcome)
            body.append(f'      <div class="abs-authorship">{tx(r["authorship_note_html"])}</div>')
    body.append(f'      <p class="abs-access">{access_line(r)}</p>')
    ab = r['abstract_display']['html']
    if ab:
        body.append('      <div class="abs-abstract">')
        for i, h in enumerate(ab):
            label = '<span class="abs-descriptor">Abstract:</span> ' if i == 0 else ''
            body.append(f'        <p>{label}{tx(h)}</p>')
        body.append('      </div>')

    sl = data.get('system_links') or {}
    comments = ''.join(f'<p>{link_systems(tx(c), sl)}</p>' for c in r['comments_html'])
    for m in r['comments_more']:
        comments += f'<details class="abs-more"><summary>{escape(m["summary"])} (full text)</summary><div>{tx(m["html"])}</div></details>'
    mv, hv = r['machine_verified'], r['human_verified']
    explain = f'<a href="/{PAGE}#verification">'
    rows = [('Comments:', comments)]
    subj_html = subjects_html(r, 'primary-subject')
    if subj_html:
        rows.append(('Subjects:', subj_html))
    for label, key in (('Agent system:', 'agent_system_html'), ('Base model:', 'base_model_html'), ('Writing:', 'writing_html'), ('Writing model:', 'writing_model_html')):
        # Jihao 11-12: the four rows always shown, "-" when not known
        rows.append((label, link_systems(r[key], sl) if r.get(key) else '-'))
    for name, f in (('Machine-verified', mv), ('Human-verified', hv)):
        cell = mark(f['value']) if f['value'] else mark(None, hidden=False) + ' <span class="lx-note">not stated</span>'
        if f.get('display_text'):           # Jihao 11-12 (2609.0001): the text "(see comments above)" instead of the mark
            cell = f'<span class="lx-note">{escape(f["display_text"])}</span>'
        if f.get('display_note'):           # Jihao 11-12 (Comments routing): who verified, or "in progress"
            cell += f' <span class="lx-note">({tx(f["display_note"])})</span>'
        rows.append((f'{explain}{name}</a>:', cell))
    links = links_html(r, True)
    rows.append(('Links:', ' &middot; '.join(links) if links else '<span class="lx-note">none</span>'))
    nv = len(r['history_display'])
    cite = f'<span class="abs-cite"><a href="/abs/{rid}/">{PREFIX}{rid}</a></span>'
    if nv:
        cite += f'<span class="abs-url">(or <a href="/abs/{rid}v{nv}/">{PREFIX}{rid}v{nv}</a> for this version)</span>'
    rows.append(('Cite as:', cite + f'<span class="abs-url"><a href="/abs/{rid}/">{url}</a></span>'))
    body.append('      <table class="abs-meta">')
    for k, v in rows:
        body.append(f'        <tr><th scope="row">{k}</th><td>{v}</td></tr>')
    body.append('      </table>')
    if r['history_display']:
        body += ['      <div class="abs-subhist">', '      <h2 class="abs-history-title">Submission history</h2>',
                 '      <ul class="abs-history">', history_items(r), '      </ul>', '      </div>']
    body += ['    </main>', '']
    return head(title, description(r), url, 'article', '/') + header(nav, PAGE) + '\n\n' + '\n'.join(body) + FOOTER.format(root='/')


def version_page(data, nav, r, k, permanent):
    """The page of version k, as arXiv's abs/<id>v<k>. The last version shows the full abstract page; an earlier
    version shows its own date, change, pages, SHA-256 and PDF only, since the site keeps the abstract and the
    comments of the current version only."""
    rid, h = r['id'], r['history_display']
    n, v = len(h), h[k - 1]
    if k == n:
        page = abs_page(data, nav, r, permanent)
        for a, b in ((f'<span class="abs-id">{PREFIX}{rid}</span>', f'<span class="abs-id">{PREFIX}{rid}v{k}</span>'),
                     (f'<title>[{PREFIX}{rid}] ', f'<title>[{PREFIX}{rid}v{k}] ')):
            if page.count(a) != 1:
                raise SystemExit(f'{rid}v{k}: cannot mark the version on the abstract page')
            page = page.replace(a, b)
        return page
    url = f'{SITE}/abs/{rid}v{k}/'
    t1, tn = escape(h[0]['time']), escape(h[-1]['time'])
    if k == 1:
        dl = f'[Posted on {t1} (this version), latest version {tn} (<a href="/abs/{rid}/">v{n}</a>)]'
    else:
        dl = (f'[Posted on {t1} (<a href="/abs/{rid}v1/">v1</a>), revised {escape(v["time"])} (this version, v{k}), '
              f'latest version {tn} (<a href="/abs/{rid}/">v{n}</a>)]')
    title_html = v.get('title_html') or r['title_html']
    authors_html = v['authors_html'] if 'authors_html' in v else r['authors_html']
    access = f'<a href="/pdf/{rid}v{k}.pdf">PDF</a> ({v["pages"]} pages)' if v['pdf_path'] else f'{v["pages"]} pages (the PDF of this version is not kept)'
    body = [
        '    <main id="main-content" class="lx-page lx-abs" tabindex="-1">',
        '      <div class="abs-header">',
        f'        <span class="abs-id">{PREFIX}{rid}v{k}</span>',
        f'        <span class="abs-browse"><a href="/{PAGE}">liuxiv</a></span>',
        '      </div>',
    ] + ([f'      <div class="lx-ainote">{data["liuxiv_page"]["ai_note_page_html"]}</div>'] if r.get('first_version_ai') else []) + [
        f'      <p class="abs-dateline">{dl}</p>',
        f'      <h1 class="abs-title">{title_html}</h1>',
    ]
    if authors_html:
        body.append(f'      <div class="entry-authors abs-authors">{root_relative(retarget(authors_html, permanent, "/"))}</div>')
        if r.get('authorship_note_html'):   # Jihao 11-12: the authorship note holds for every version
            body.append(f'      <div class="abs-authorship">{root_relative(retarget(r["authorship_note_html"], permanent, "/"))}</div>')
    body.append(f'      <p class="abs-access">{access}</p>')
    ab = v.get('abstract_html') or []
    if ab:
        body.append('      <div class="abs-abstract">')
        for i2, h2 in enumerate(ab):
            label = '<span class="abs-descriptor">Abstract:</span> ' if i2 == 0 else ''
            body.append(f'        <p>{label}{root_relative(retarget(h2, permanent, "/"))}</p>')
        body.append('      </div>')
    vnote = f' {link_systems(root_relative(retarget(v["comments_html"], permanent, "/")), data.get("system_links") or {})}' if v.get('comments_html') else ''   # Jihao 11-12: what changed in this version
    solves = f' {link_systems(root_relative(retarget(v["solves_html"], permanent, "/")), data.get("system_links") or {})}' if v.get('solves_html') else ''   # Jihao 11-12: what this version solves
    rows = [('Comments:', f'{v["pages"]} pages.{solves}{vnote}')]
    subj_html = subjects_html(r, 'primary-subject')
    if subj_html:
        rows.append(('Subjects:', subj_html))
    sl = data.get('system_links') or {}
    for label, key in (('Agent system:', 'agent_system_html'), ('Base model:', 'base_model_html'), ('Writing:', 'writing_html'), ('Writing model:', 'writing_model_html')):
        # Jihao 11-12: the four rows always shown, "-" when not known; they describe the record (Writing = who wrote the first draft), so every version shows them
        val = v[key] if key in v else r.get(key)   # Jihao 11-12: the writer may differ between versions
        rows.append((label, link_systems(val, sl) if val else '-'))
    # Jihao 11-12: every version shows the verification rows and the links
    explain = f'<a href="/{PAGE}#verification">'
    hvv = v.get('human_verified') or r['human_verified']
    for name, f in (('Machine-verified', r['machine_verified']), ('Human-verified', hvv)):
        cell = mark(f['value']) if f['value'] else mark(None, hidden=False) + ' <span class="lx-note">not stated</span>'
        if f.get('display_text'):   # 2609.0001: no mark; the explanation is in the comments of the latest version
            cell = f'<span class="lx-note">(see the comments of <a href="/abs/{rid}/">the latest version</a>)</span>'
        elif f.get('display_note'):
            cell += f' <span class="lx-note">({root_relative(retarget(f["display_note"], permanent, "/"))})</span>'
        rows.append((f'{explain}{name}</a>:', cell))
    links = links_html(r, True)
    rows.append(('Links:', ' &middot; '.join(links) if links else '<span class="lx-note">none</span>'))
    rows.append(('Cite as:', f'<span class="abs-cite"><a href="/abs/{rid}v{k}/">{PREFIX}{rid}v{k}</a></span><span class="abs-url"><a href="/abs/{rid}v{k}/">{url}</a></span>'))
    body.append('      <table class="abs-meta">')
    for a2, b2 in rows:
        body.append(f'        <tr><th scope="row">{a2}</th><td>{b2}</td></tr>')
    body.append('      </table>')
    body += ['      <div class="abs-subhist">', '      <h2 class="abs-history-title">Submission history</h2>',
             '      <ul class="abs-history">', history_items(r, k), '      </ul>', '      </div>', '    </main>', '']
    title = f'[{PREFIX}{rid}v{k}] {page_title(r)}'
    return head(title, description(r), url, 'article', '/') + header(nav, PAGE) + '\n\n' + '\n'.join(body) + FOOTER.format(root='/')


def redirect_page(nav, r):
    rid = r['id']
    url = f'{SITE}/abs/{rid}/'
    return f'''<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="UTF-8" />
  <meta name="viewport" content="width=device-width, initial-scale=1.0" />
  <title>{escape(page_title(r), quote=False)}</title>
  <meta http-equiv="refresh" content="0; url=/abs/{rid}/" />
  <link rel="canonical" href="{url}" />
{ICONS}
  <link rel="stylesheet" href="/assets/fonts.css" />
  <link rel="stylesheet" href="/assets/style.css" />
</head>
<body>
{header(nav, PAGE)}

    <main id="main-content" tabindex="-1">
      <h1>{r["title_html"]}</h1>
      <p>This page has moved to <a href="/abs/{rid}/">{HOST}/abs/{rid}/</a> ({PREFIX}{rid}).</p>
    </main>

''' + FOOTER.format(root='/')


# ---------------------------------------------------------------- the liuxiv page

def day_heading(t):
    """'2026-11-10 13:53:46 UTC' -> 'Tue, 10 Nov 2026' (the day headings of arXiv listings)."""
    d = datetime.date.fromisoformat(t[:10])
    return f'{d:%a}, {d.day} {d:%b %Y}'


def liuxiv_page(data, nav, records, permanent):
    """An arXiv-style listing: one heading per day of first posting (UTC), entries numbered [1], [2], ... from the top."""
    p = data['liuxiv_page']
    local = lambda h: re.sub(r'\b(href|src)="/(?!/)', r'\1="', retarget(h, permanent, ''))  # root page: relative links
    groups = []
    for r in sorted(records, key=lambda r: r['id'], reverse=True):
        h = r['history_display']
        if h:
            key, heading, anchor = h[0]['time'][:10], day_heading(h[0]['time']), 'd' + h[0]['time'][:10]
        else:  # no posting time is shown for this entry (Notes 6): a month heading from its identifier
            d = datetime.date(2000 + int(r['id'][:2]), int(r['id'][2:4]), 1)
            key, heading, anchor = f'{d:%Y-%m}-00', f'{d:%B %Y} (posting date not listed)', f'm{d:%Y-%m}'
        if groups and groups[-1][0] == key:
            groups[-1][3].append(r)
        else:
            if any(g[0] == key for g in groups):
                raise SystemExit(f'liuxiv.html: entries of {key} are not consecutive in identifier order')
            groups.append((key, heading, anchor, [r]))
    months = []
    for key, heading, anchor, rs in groups:
        m = key[:7]
        if not months or months[-1][0] != m:
            d = datetime.date(int(m[:4]), int(m[5:]), 1)
            months.append((m, f'{d:%b %Y}', anchor))
    parts, n = [], 0
    for key, heading, anchor, rs in groups:
        count = '1 entry' if len(rs) == 1 else f'{len(rs)} entries'
        parts.append(f'      <h3 class="lx-day" id="{anchor}">{heading} <span class="lx-count">({count})</span></h3>')
        parts.append('      <dl class="lx-list">')
        for r in rs:
            n += 1
            rid = r['id']
            pdf = f'pdf/{rid}.pdf'
            parts.append(f'        <dt><span class="lx-num">[{n}]</span> <a href="abs/{rid}/" title="Abstract">{PREFIX}{rid}</a>'
                         f' <span class="lx-formats">[<a href="abs/{rid}/">abs</a>, <a href="{pdf}">pdf</a>]</span></dt>')
            parts.append('        <dd>')
            parts.append(f'          <div class="lx-ttl"><span class="visually-hidden">Title: </span>{r["title_html"]}</div>')
            if r['authors_html']:
                am = r.get('authorship_listing_html')   # Jihao 11-12: the listing marks the authorship
                star = '' if not am else ('<a class="lx-star" href="#authorship" title="Authorship is temporary">*</a>' if am == '*'
                                          else f' <span class="lx-authnote">{local(am)}</span>')
                parts.append(f'          <div class="lx-auth">{local(r["authors_html"])}{star}</div>')
            if r['comments_html']:
                comm = link_systems(local(r['comments_html'][0]), data.get('system_links') or {})
                parts.append(f'          <div class="lx-comm"><span class="lx-desc">Comments:</span> {comm}</div>')
            subj_html = subjects_html(r, 'lx-primary')
            if subj_html:
                parts.append(f'          <div class="lx-subj"><span class="lx-desc">Subjects:</span> {subj_html}</div>')
            hv = r['human_verified']   # Jihao 11-12 (2609.0001): a display_text replaces the mark
            hv_cell = '<span class="lx-note">' + escape(hv['display_text']) + '</span>' if hv.get('display_text') else mark(hv['value'])
            parts.append('          <div class="lx-ver">'
                         f'<span class="lx-desc"><a href="#verification">Machine-verified</a>:</span> {mark(r["machine_verified"]["value"])}'
                         f' <span class="lx-desc"><a href="#verification">Human-verified</a>:</span> {hv_cell}</div>')
            links = links_html(r, False)
            if links:
                parts.append(f'          <div class="lx-lnk"><span class="lx-desc">Links:</span> {" &middot; ".join(local(x) for x in links)}</div>')
            parts.append('        </dd>')
        parts.append('      </dl>')
    intro = '\n'.join(f'      <p class="lx-intro">{x}</p>' for x in p['intro_html'])
    index = ' | '.join(f'<a href="#{a}">{label}</a>' for _, label, a in months)
    main = f'''    <main id="main-content" class="lx-page" tabindex="-1">
      <h1>{escape(p["title"])}</h1>
{intro}
      <div class="lx-ainote" id="ai-generated">{p["ai_note_html"]}</div>
      <p class="lx-intro" id="verification">{p["verification_html"]}</p>
      <p class="lx-intro" id="authorship">{p["authorship_legend_html"]}</p>
      <p class="lx-index">Total of {n} entries: {index}</p>

{chr(10).join(parts)}
    </main>

'''
    nav_html = header(nav, PAGE).replace('href="/', 'href="')  # root page: relative links, as on the other root pages
    page = head(f'{p["title"]} — Jihao Liu', p['description'], f'{SITE}/{PAGE}', 'website', '') + nav_html + '\n\n' + main + FOOTER.format(root='')
    return render_page(page, PAGE, nav)


# ---------------------------------------------------------------- driver

def pdf_copies(records):
    """Short arXiv-style PDF addresses: pdf/<id>.pdf is the current version, pdf/<id>v<k>.pdf version k.
    They are byte-identical copies; the old addresses under papers/, notes/ and ai-results/ stay as they are."""
    want = {}
    for r in records:
        want[ROOT / 'pdf' / f'{r["id"]}.pdf'] = ROOT / r['old_urls']['pdf_current']
        for v in r['history_display']:
            if v['pdf_path']:
                want[ROOT / 'pdf' / f'{r["id"]}v{v["v"]}.pdf'] = ROOT / v['pdf_path']
    for dst, src in want.items():
        if not src.is_file():
            raise SystemExit(f'missing PDF {src.relative_to(ROOT)} for {dst.name}')
    return want


def outputs():
    data, nav, records, permanent = load()
    out = {ROOT / PAGE: liuxiv_page(data, nav, records, permanent)}
    for r in records:
        out[ROOT / 'abs' / r['id'] / 'index.html'] = abs_page(data, nav, r, permanent)
        for v in r['history_display']:
            out[ROOT / 'abs' / f'{r["id"]}v{v["v"]}' / 'index.html'] = version_page(data, nav, r, v['v'], permanent)
        if r['old_urls']['permanent_page']:
            out[ROOT / r['old_urls']['permanent_page']] = redirect_page(nav, r)
    wanted = {p.parent for p in out if p.parent.parent == ROOT / 'abs'}
    orphans = sorted(p for p in (ROOT / 'abs').glob('*/index.html') if p.parent not in wanted) if (ROOT / 'abs').is_dir() else []
    copies = pdf_copies(records)
    managed = sorted(p for p in (ROOT / 'pdf').glob('*.pdf') if PDF_RE.match(p.name)) if (ROOT / 'pdf').is_dir() else []
    orphans += [p for p in managed if p not in copies]
    return out, orphans, copies


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--check', action='store_true')
    args = ap.parse_args()
    out, orphans, copies = outputs()
    stale = [p for p, text in out.items() if not p.is_file() or p.read_text() != text]
    stale_pdf = [d for d, s in copies.items() if not d.is_file() or d.read_bytes() != s.read_bytes()]
    rel = lambda p: str(p.relative_to(ROOT))
    if args.check:
        if stale or orphans or stale_pdf:
            raise SystemExit('liuxiv pages are stale: ' + ', '.join(map(rel, stale + orphans + stale_pdf)))
        print(f'liuxiv pages match data/liuxiv.json ({len(out)} files, {len(copies)} PDF copies).')
        return
    for p in stale:
        p.parent.mkdir(parents=True, exist_ok=True)
        p.write_text(out[p])
    for d in stale_pdf:
        d.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(copies[d], d)
    for p in orphans:
        p.unlink()
        if not any(p.parent.iterdir()):
            p.parent.rmdir()
    print(f'liuxiv: {len(out)} files checked, {len(stale)} written, {len(stale_pdf)} PDF copies written, {len(orphans)} orphans removed.')


if __name__ == '__main__':
    main()
