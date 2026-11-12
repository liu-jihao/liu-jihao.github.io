#!/usr/bin/env python3
"""Read-only checks of the liuxiv build (test build 2) against a base ref (default: main).

  1  old URLs: every file the base site serves still exists; every former permanent page
     papers/<slug>.html redirects to the abstract page of the same paper; ai-results.html and
     notes.html equal the base except for the one shared navigation line, so every anchor is kept;
  2  every PDF is byte-identical to the base; only the intended files changed (each other root page
     differs from the base by the navigation line only; the stylesheet only gained lines at its end);
     the served PDFs have the SHA-256 recorded in the data;
  3  every record is on liuxiv.html exactly once and has exactly one abstract page;
  4  every internal link (href, src, meta refresh) of all root pages, abstract pages and redirects
     resolves, fragments included;
  5  every posting time and SHA-256 that the base prints for a record is on its abstract page, and the
     abstract page shows no time or SHA-256 that the base does not print for that record, except the
     exact posting times taken from git where the site gives only a date or a month;
  6  no Update paragraph, Disclaimer, Status line or Timestamp paragraph of the base is left on an
     abstract page; Notes no. 6 shows nothing that was removed from the site; every check mark and
     cross has a source sentence that the base prints for that record, and the pages agree with the
     data; old numbers of AI results and Notes are gone from the abstract pages; the abstracts are the
     site's own paragraphs or the PDF's abstract (checked with pdftotext when it is installed);
  7  the generator is idempotent (build_liuxiv.py --check) and left no placeholders.

Prints a JSON report; exits 1 if any check fails.   python3 scripts/check_liuxiv.py [--base main]
"""
from pathlib import Path
from html.parser import HTMLParser
from html import unescape
from urllib.parse import urlsplit, unquote
import argparse, hashlib, json, re, shutil, subprocess, sys, traceback, unicodedata

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'scripts'))
import build_liuxiv as B  # noqa: E402

TIME = re.compile(r'\d{4}-\d{2}-\d{2} \d{2}:\d{2}:\d{2} UTC')
HEX = re.compile(r'(?<![0-9a-f])[0-9a-f]{64}(?![0-9a-f])')
NAV_LINE = '          <a href="liuxiv.html">liuxiv</a>'
OLD_PAGES = ('ai-results.html', 'notes.html', 'ai-results-about.html')  # no longer updated; still served, no longer linked
FROZEN_NOTE = re.compile(r'^ *<p class="lx-frozen">.*?</p>\n', re.M)
LIST_LINK = re.compile(r'<a href="(?:liuxiv|notes|ai-results)\.html"[^>]*>[^<]*</a>\n?(?: *)')


def norm_root(name, text):
    """A root page with the intended edits undone: navigation links to liuxiv, notes and ai-results removed,
    the no-longer-updated notice removed, and the one retargeted link of the Fujita supplement put back."""
    text = FROZEN_NOTE.sub('', text)
    text = LIST_LINK.sub('', text)
    return text.replace('href="abs/2609.0003/"', 'href="ai-results.html#paper-3"')


def git(*args, binary=False):
    r = subprocess.run(['git', '-C', str(ROOT), *args], capture_output=True)
    if r.returncode:
        raise SystemExit('git ' + ' '.join(args) + ': ' + r.stderr.decode())
    return r.stdout if binary else r.stdout.decode()


class Page(HTMLParser):
    def __init__(self, text):
        super().__init__(convert_charrefs=True)
        self.ids, self.urls, self.tags, self.refresh, self.canonical = [], [], [], [], []
        self.feed(text)

    def handle_starttag(self, tag, attrs):
        a = dict(attrs)
        self.tags.append(tag)
        if 'id' in a:
            self.ids.append(a['id'])
        for k in ('href', 'src'):
            if a.get(k) and not (tag == 'link' and a.get('rel') == 'canonical'):
                self.urls.append(a[k])
        if tag == 'link' and a.get('rel') == 'canonical':
            self.canonical.append(a.get('href'))
        if tag == 'meta' and (a.get('http-equiv') or '').lower() == 'refresh':
            m = re.match(r'\s*\d+\s*;\s*url=(.+)$', a.get('content', ''), re.I)
            if m:
                self.refresh.append(m.group(1).strip())
                self.urls.append(m.group(1).strip())


def plain(html):
    return re.sub(r'\s+', ' ', unescape(re.sub(r'<[^>]+>', ' ', html or ''))).strip()


def tight(html):
    return re.sub(r'\s+', '', unescape(re.sub(r'<[^>]+>', '', html or '')))


def entries(html):
    out = {m.group(2): m.group(1) for m in re.finditer(r'(      <div class="entry" id="([^"]+)">\n.*?\n      </div>)(?=\n)', html, re.S)}
    # 11-13: a record first listed under Publications (the Popa–Schnell counterexample, liuxiv:2610.0007)
    out.update({m.group(2): m.group(1) for m in re.finditer(r'(      <article class="entry" id="([^"]+)"[^>]*>\n.*?\n      </article>)(?=\n)', html, re.S)})
    return out


def main_of(html):
    m = re.search(r'<main\b.*?</main>', html, re.S)
    return m.group(0) if m else ''


def resolve(page_path, url):
    u = urlsplit(url)
    if u.scheme or u.netloc:
        return None, u.fragment
    if not u.path:
        return page_path, u.fragment
    target = (ROOT / unquote(u.path).lstrip('/')) if u.path.startswith('/') else (page_path.parent / unquote(u.path))
    target = target.resolve()
    if target.is_dir() or u.path.endswith('/'):
        target = target / 'index.html'
    return target, u.fragment


def served(path, excluded):
    return not any(path == e.rstrip('/') or path.startswith(e.rstrip('/') + '/') for e in excluded)


def norm_letters(s):
    """Letters and digits only, after compatibility decomposition (ligatures, ℝ, accents)."""
    return re.sub(r'[^a-z0-9]', '', unicodedata.normalize('NFKD', unescape(s)).lower())


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--base', default='main')
    base = ap.parse_args().base
    data, nav, records, permanent = B.load()
    rep, fail = {}, []

    def check(name, ok, **detail):
        rep[name] = {'ok': bool(ok), **detail}
        if not ok:
            fail.append(name)

    base_files = [f for f in git('ls-tree', '-r', '--name-only', base).split('\n') if f]
    excluded = re.findall(r'^\s*-\s*(\S+)\s*$', (ROOT / '_config.yml').read_text(), re.M)
    base_served = [f for f in base_files if served(f, excluded)]
    show = lambda p: git('show', f'{base}:{p}')
    abs_html = {r['id']: (ROOT / 'abs' / r['id'] / 'index.html').read_text() if (ROOT / 'abs' / r['id'] / 'index.html').is_file() else '' for r in records}
    lx_html = (ROOT / B.PAGE).read_text() if (ROOT / B.PAGE).is_file() else ''

    def base_texts(r):
        """The base site's list entry and permanent page <main> for a record."""
        lst = r['listing']
        e = entries(show(lst['page']))[lst['anchor']]
        for al in r.get('listing_aliases', []):  # a paper listed twice on the base (Notes no. 8 = first version of AI results no. 30)
            e += '\n' + entries(show(al['page']))[al['anchor']]
        pp = main_of(show(r['old_urls']['permanent_page'])) if r['old_urls']['permanent_page'] else ''
        return e, pp

    try:
        # ---- 1: old URLs
        missing = [f for f in base_served if not (ROOT / f).exists()]
        check('1a_base_served_files_exist', not missing, base_served_files=len(base_served), missing=missing)

        old_pp = sorted(f for f in base_files if re.fullmatch(r'papers/[^/]+\.html', f))
        bad = []
        for f in old_pp:
            rid = permanent.get(f)
            page = Page((ROOT / f).read_text())
            target = f'/abs/{rid}/'
            old_h1 = re.search(r'<h1>(.*?)</h1>', show(f), re.S)
            new_h1 = re.search(r'<h1 class="abs-title">(.*?)</h1>', abs_html.get(rid, ''), re.S)
            old_title = re.search(r'<title>(.*?)</title>', show(f), re.S).group(1)
            new_title = re.search(r'<title>(.*?)</title>', (ROOT / f).read_text(), re.S).group(1)
            problems = []
            if not rid:
                problems.append('no record')
            if page.refresh != [target]:
                problems.append(f'refresh {page.refresh}')
            if page.canonical != [f'{B.SITE}{target}']:
                problems.append(f'canonical {page.canonical}')
            if page.urls.count(target) < 2:
                problems.append('no visible link')
            if not (old_h1 and new_h1 and plain(old_h1.group(1)) == plain(new_h1.group(1))):
                problems.append('title differs')
            if unescape(old_title) != unescape(new_title):
                problems.append('<title> differs')
            if problems:
                bad.append((f, problems))
        check('1b_former_permanent_pages_redirect', not bad and len(old_pp) == len(permanent),
              former_permanent_pages=len(old_pp), problems=bad)

        bad = []
        for page in ('ai-results.html', 'notes.html'):
            old, new = show(page), (ROOT / page).read_text()
            if new.count('class="lx-frozen"') != 1:
                bad.append((page, 'no single notice that the page is no longer updated'))
            if norm_root(page, new) != norm_root(page, old):
                bad.append((page, 'differs from the base beyond the navigation and the notice'))
            if main_of(norm_root(page, new)) != main_of(norm_root(page, old)):
                bad.append((page, '<main> differs beyond the notice'))
            if Page(old).ids != Page(new).ids:
                bad.append((page, 'anchors differ'))
        check('1c_listing_pages_as_on_base', not bad, problems=bad,
              anchors=sum(len(Page(show(p)).ids) for p in ('ai-results.html', 'notes.html')))

        pdfs_base = [f for f in base_served if f.endswith('.pdf')]
        check('1d_old_pdf_paths_exist', all((ROOT / f).is_file() for f in pdfs_base), pdf_paths=len(pdfs_base))

        # ---- 2: PDFs unchanged, only intended files changed
        stat = git('diff', '--stat', base, '--', '*.pdf', ':!pdf/')
        tree = {}
        for line in git('ls-tree', '-r', base).split('\n'):
            if line:
                meta, path = line.split('\t', 1)
                tree[path] = meta.split()[2]
        pdf_all = sorted(str(p.relative_to(ROOT)) for p in ROOT.rglob('*.pdf') if '.git' not in p.parts)
        copies = [f for f in pdf_all if re.fullmatch(r'pdf/\d{4}\.\d{4}(?:v\d+)?\.pdf', f)]
        pdf_now = [f for f in pdf_all if f not in copies]
        blob_bad = [f for f in pdf_now if tree.get(f) != git('hash-object', f).strip()]
        want = {}
        for r in records:  # the short addresses: pdf/<id>.pdf = current version, pdf/<id>v<k>.pdf = version k
            want[f'pdf/{r["id"]}.pdf'] = r['old_urls']['pdf_current']
            for v in r['history_display']:
                if v['pdf_path']:
                    want[f'pdf/{r["id"]}v{v["v"]}.pdf'] = v['pdf_path']
        copy_bad = [f for f in copies if f not in want or (ROOT / f).read_bytes() != (ROOT / want[f]).read_bytes()] + [f for f in want if f not in copies]
        check('2a_pdfs_identical', stat.strip() == '' and not blob_bad and set(pdf_now) == {f for f in base_files if f.endswith('.pdf')} and not copy_bad,
              git_diff_stat_pdf=stat.strip(), pdfs=len(pdf_now), differing=blob_bad, short_copies=len(copies), short_copy_problems=copy_bad)

        changed = [l.split('\t') for l in git('diff', '--name-status', base).split('\n') if l]
        untracked = [l[3:] for l in git('status', '--porcelain', '--untracked-files=all').split('\n') if l.startswith('??')]
        intended = lambda p: (p in ('liuxiv.html', 'assets/style.css', 'assets/search-index.js', 'data/navigation.json', 'data/liuxiv.json',
                                    'data/liuxiv_map.tsv', 'scripts/build_liuxiv.py', 'scripts/check_liuxiv.py', 'scripts/build_discovery.py', 'scripts/check_site.py',
                                    'templates/search.html')   # 11-13: the search page's quick links go to liuxiv instead of the two frozen pages
                              or re.fullmatch(r'pdf/\d{4}\.\d{4}(?:v\d+)?\.pdf', p)
                              or re.fullmatch(r'[^/]+\.html', p) or re.fullmatch(r'papers/[^/]+\.html', p) or re.fullmatch(r'abs/\d{4}\.\d{4}(?:v\d+)?/index\.html', p))
        odd = [c for c in changed if c[0] not in ('M', 'A') or not intended(c[-1])] + [['??', p] for p in untracked if not intended(p)]
        nav_only, nav_bad = 0, []
        for c in changed:
            p = c[-1]
            if c[0] == 'M' and re.fullmatch(r'[^/]+\.html', p):
                if norm_root(p, (ROOT / p).read_text()) == norm_root(p, show(p)):
                    nav_only += 1
                else:
                    nav_bad.append(p)
        css_old, css_new = show('assets/style.css'), (ROOT / 'assets/style.css').read_text()
        sha_bad = [(r['id'], v['v']) for r in records for v in r['versions']
                   if v['pdf_path'] and hashlib.sha256((ROOT / v['pdf_path']).read_bytes()).hexdigest() != v['sha256']]
        check('2b_only_intended_files_changed', not odd and not nav_bad and css_new.startswith(css_old),
              changed_vs_base={k: sum(1 for c in changed if c[0] == k) for k in sorted({c[0] for c in changed})},
              untracked=len(untracked), root_pages_differing_by_the_navigation_line_only=nav_only,
              root_pages_differing_otherwise=nav_bad, unexpected=odd, stylesheet_append_only=css_new.startswith(css_old))
        check('2c_served_pdf_sha256_match_data', not sha_bad, mismatches=sha_bad)

        # ---- 3: each record once on liuxiv.html and one abstract page
        bad = []
        rows = re.findall(r'<dt><span class="lx-num">\[\d+\]</span> ?<a href="abs/(\d{4}\.\d{4})/" title="Abstract">liuxiv:(\d{4}\.\d{4})</a>.*?</dt>\s*<dd>(.*?)</dd>', lx_html, re.S)
        row_ids = [a for a, b, _ in rows]
        if [a for a, b, _ in rows] != [b for a, b, _ in rows]:
            bad.append('row link and identifier differ')
        if sorted(row_ids) != sorted(r['id'] for r in records) or len(row_ids) != len(set(row_ids)):
            bad.append('rows of liuxiv.html are not the records, each once')
        abs_all = sorted(p.name for p in (ROOT / 'abs').iterdir() if p.is_dir())
        abs_dirs = [x for x in abs_all if 'v' not in x]
        ver_want = sorted(f'{r["id"]}v{v["v"]}' for r in records for v in r['history_display'])
        if [x for x in abs_all if 'v' in x] != ver_want:
            bad.append('version pages abs/<id>v<k>/ are not exactly the versions of the records')
        for r in records:
            t = abs_html[r['id']]
            p = Page(t)
            if t.count(f'<span class="abs-id">liuxiv:{r["id"]}</span>') != 1:
                bad.append((r['id'], 'abstract page id'))
            if [p.tags.count(x) for x in ('h1', 'main', 'title')] != [1, 1, 1] or len(p.ids) != len(set(p.ids)) or 'main-content' not in p.ids:
                bad.append((r['id'], 'abstract page structure'))
            if p.canonical != [f'{B.SITE}/abs/{r["id"]}/']:
                bad.append((r['id'], 'canonical'))
        check('3_each_record_once_on_liuxiv_and_one_abstract_page', not bad and abs_dirs == sorted(r['id'] for r in records),
              records=len(records), rows_on_liuxiv_html=len(rows), abstract_pages=len(abs_dirs), problems=bad)

        # ---- 4: internal links
        pages = sorted({ROOT / 'abs' / r['id'] / 'index.html' for r in records} | {ROOT / f for f in old_pp} | set(ROOT.glob('*.html')))
        cache, bad, n = {}, [], 0
        for page_path in pages:
            for url in Page(page_path.read_text()).urls:
                target, frag = resolve(page_path, url)
                if target is None:
                    continue
                n += 1
                if not target.is_file():
                    bad.append((str(page_path.relative_to(ROOT)), url, 'missing'))
                    continue
                if frag and target.suffix == '.html':
                    cache.setdefault(target, set(Page(target.read_text()).ids))
                    if unquote(frag) not in cache[target]:
                        bad.append((str(page_path.relative_to(ROOT)), url, 'missing fragment'))
        check('4_internal_links_resolve', not bad, pages=len(pages), local_links_checked=n, problems=bad)

        # ---- 5: times and SHA-256: kept, and nothing new
        bad, counts, git_times = [], {'base_times': 0, 'base_hashes': 0}, []
        for r in records:
            e, pp = base_texts(r)
            old_t, old_h = set(TIME.findall(e + pp)), set(HEX.findall(e + pp))
            counts['base_times'] += len(old_t)
            counts['base_hashes'] += len(old_h)
            page = main_of(abs_html[r['id']])
            # Jihao 11-11: the Danus verifier acceptance sentences are not shown; their times are exempt
            dropped = set(TIME.findall(' '.join(re.findall(r'accepted by the Danus verifier[^.]*', plain(e + ' ' + pp)))))
            for x in sorted((old_t | old_h) - dropped):
                if x not in page:
                    bad.append((r['id'], 'abstract page lacks', x))
            allowed = {h['time'] for h in r['history_display'] if h['time_source'] != 'stated on the site'}
            for x in sorted(set(TIME.findall(page)) - old_t):
                if x in allowed:
                    git_times.append((r['id'], r['old']['label'], x))
                else:
                    bad.append((r['id'], 'time not printed on the base', x))
            # Jihao 11-11 (every PDF has a SHA-256): a hash may also be the SHA-256 of one of the record's version PDFs, computed from the file
            served_h = {hashlib.sha256((ROOT / v['pdf_path']).read_bytes()).hexdigest() for v in r['history_display'] if v.get('pdf_path')}
            git_h = {v['sha256'] for v in r['history_display'] if not v.get('pdf_path') and v.get('sha256') and v.get('sha256_source', '').startswith('git ')}
            # 11-12: the SHA-256 of an ancillary zip may also be shown when it is computed from the file (2609.0003, as for 2611.0006 and 2611.0007)
            anc_h = {hashlib.sha256((ROOT / a['path']).read_bytes()).hexdigest() for a in r['ancillary'] if a['path'].endswith('.zip') and (ROOT / a['path']).is_file()}
            for x in sorted(set(HEX.findall(page)) - old_h - served_h - git_h - anc_h):
                bad.append((r['id'], 'SHA-256 neither printed on the base nor of a version PDF or ancillary zip', x))
            for v in r['history_display']:
                if v.get('pdf_path') and v.get('sha256') != hashlib.sha256((ROOT / v['pdf_path']).read_bytes()).hexdigest():
                    bad.append((r['id'], f'v{v["v"]} SHA-256 differs from its file'))
            hm = re.search(r'<ul class="abs-history">(.*?)</ul>', page, re.S)
            items = re.findall(r'<li[^>]*>(.*?)</li>', hm.group(1)) if hm else []
            if len(items) != len(r['history_display']):
                bad.append((r['id'], 'history length'))
            for v, it in zip(r['history_display'], items):
                if (f'[v{v["v"]}]' not in it or v['time'] not in it or f'{v["pages"]} pages' not in it
                        or (v['sha256'] or '') not in it or bool(v['sha256']) != ('SHA-256' in it)
                        or f'href="/abs/{r["id"]}v{v["v"]}/"' not in it):
                    bad.append((r['id'], f'history v{v["v"]}'))
                vp = (ROOT / 'abs' / f'{r["id"]}v{v["v"]}' / 'index.html')
                vt = vp.read_text() if vp.is_file() else ''
                if v['pdf_path'] and f'href="/pdf/{r["id"]}v{v["v"]}.pdf"' not in vt and v['v'] != len(r['history_display']):
                    bad.append((r['id'], f'version page v{v["v"]} has no PDF link'))
        check('5_times_and_hashes', not bad, problems=bad, git_times_shown=git_times, **counts)

        # ---- 6: what the abstract pages no longer carry, Notes 6, check marks, old numbers, abstracts
        bad = []
        for r in records:
            e, pp = base_texts(r)
            page = main_of(abs_html[r['id']])
            pt = tight(page)
            blocks = re.findall(r'<p\b[^>]*>(.*?)</p>', e + pp, re.S)
            for b in blocks:
                lab = plain(b)[:12]
                if lab.startswith(('Update (', 'Disclaimer:', 'Status:', 'Timestamp.')) and tight(b) in pt:
                    bad.append((r['id'], 'base paragraph left on the abstract page', plain(b)[:60]))
            for word in ('Disclaimer:', 'Update (', 'Status:', 'Timestamp.', 'Readers should independently check', 'This is an AI-generated manuscript'):
                if word in plain(page):
                    bad.append((r['id'], 'abstract page contains', word))
        check('6a_no_update_disclaimer_status_timestamp_paragraphs', not bad, problems=bad)

        n6 = [r for r in records if r['old']['label'] == 'Notes no. 6'][0]
        page = main_of(abs_html[n6['id']])
        row = [x for a, b, x in rows if a == n6['id']]
        # Jihao 11-11 asked for the posting date of Notes 6 (first posting, from git); the rest stays removed
        # Jihao 11-12: the authorship note of this paper ("The author is collaborating with Yoshinori Gongyo and Christopher D. Hacon ...")
        # speaks of the author, so the page shows the author of the current PDF and of arXiv:2610.35913, Jihao Liu; nothing else
        removed = {'Authors withheld': 'Authors withheld' in page,
                   'an authors line other than Jihao Liu': len(re.findall(r'class="entry-authors abs-authors">', page)) != len(re.findall(r'class="entry-authors abs-authors">Jihao Liu</div>', page)),
                   'posted here': 'posted here' in plain(page),
                   'AI results no. 6 or its old PDF path': 'ai-results/nonvanishing' in page or 'paper-6' in page,
                   'a SHA-256 other than the draft hash and the current PDF': set(HEX.findall(page)) - {'cf3ab57a51040c7939b56cd1dd3acbfdc6dd472cee2131b80e7cdf87bdb14bb8', hashlib.sha256((ROOT / n6['old_urls']['pdf_current']).read_bytes()).hexdigest()} != set(),
                   'a time other than the first posting': set(TIME.findall(page)) - {'2026-10-06 22:43:54 UTC'} != set()}
        check('6b_notes6_shows_nothing_removed', not any(removed.values()), shows={k: v for k, v in removed.items()})

        bad, table = [], []
        for r in records:
            e, pp = base_texts(r)
            base_plain = plain(e + ' ' + pp) + ' || ' + re.sub(r'\s+', ' ', unescape(re.sub(r'<[^>]+>', '', e + ' ' + pp)))
            page = abs_html[r['id']]
            row = [x for a, b, x in rows if a == r['id']]
            cells = re.findall(r'-verified</a>:</span> (<span class="lx-mark [^"]*" title="[^"]*">.<|<span class="lx-note">[^<]*<)', row[0]) if row else []
            for k, name in (('machine_verified', 'Machine-verified'), ('human_verified', 'Human-verified')):
                f = r[k]
                if f.get('display_text'):   # Jihao 11-12 (2609.0001): the text "(see comments above)" replaces the mark on both pages
                    want = f'<span class="lx-note">{unescape(f["display_text"])}'
                    if f'{name}</a>:</th><td>{want}</span>' not in page:
                        bad.append((r['id'], name, 'abstract page text differs from the data'))
                    if len(cells) != 2 or want not in cells[0 if k == 'machine_verified' else 1]:
                        bad.append((r['id'], name, 'liuxiv.html text differs from the data'))
                    continue
                sym = B.MARK[f['value']][0]
                in_pdf = False   # 11-13: a source sentence may instead be printed in the current PDF (liuxiv:2610.0007, its acknowledgements)
                if (f.get('source_where') or '').endswith('of the PDF') and f['source_text']:
                    pdf_txt = subprocess.run(['pdftotext', str(ROOT / r['old_urls']['pdf_current']), '-'], capture_output=True, text=True).stdout
                    in_pdf = re.sub(r'\s+', ' ', f['source_text']) in re.sub(r'\s+', ' ', pdf_txt)
                if f['value'] is not None and not (f.get('source_where') or '').startswith('Jihao 11-11') and not in_pdf and (not f['source_text'] or f['source_text'] not in base_plain):
                    bad.append((r['id'], name, 'source sentence not printed on the base'))
                m = re.search(rf'{name}</a>:</th><td><span class="lx-mark [^"]*" title="[^"]*">(.)<', page)
                if not m or m.group(1) != sym:
                    bad.append((r['id'], name, 'abstract page mark differs from the data'))
                idx = 0 if k == 'machine_verified' else 1
                if len(cells) != 2 or f'>{sym}<' not in cells[idx]:
                    bad.append((r['id'], name, 'liuxiv.html mark differs from the data'))
            table.append((r['id'], r['machine_verified']['value'], r['human_verified']['value']))
        not_stated = [(r['id'], r['old']['label'], [n for k, n in (('machine_verified', 'machine'), ('human_verified', 'human')) if r[k]['value'] is None]) for r in records
                      if r['machine_verified']['value'] is None or r['human_verified']['value'] is None]
        check('6c_marks_have_sources_and_match', not bad, problems=bad, not_stated=not_stated,
              counts={k: {v: sum(1 for t in table if t[i] == v) for v in ('yes', 'no', None)} for i, k in ((1, 'machine'), (2, 'human'))})

        bad = []
        # allowed: numbers of the Publications page, OpenAI's paper number, and the "Listed on" row
        allowed_ctx = re.compile(r'(Publications no\. \d+|no\. 75 in my list of publications|no\. 77 under Publications|Paper No\. 038|(?:Formerly listed|Listed) on: [^,]*, no\. \d+)')
        for r in records:
            text = re.sub(r' ([,.;:])', r'\1', plain(main_of(abs_html[r['id']])).replace(' ', ' '))
            spans = [(a.start(), a.end()) for a in allowed_ctx.finditer(text)]
            for m in re.finditer(r'\b[Nn]o\. ?\d+', text):
                if not any(a <= m.start() < b for a, b in spans):
                    bad.append((r['id'], text[max(0, m.start() - 40):m.end() + 40]))
            for phrase in ('proof below', 'update below', 'PDF below'):
                if phrase in text:
                    bad.append((r['id'], f'"{phrase}"'))
        check('6d_old_numbers_replaced', not bad, problems=bad)

        # Jihao 11-11: every abstract is the one printed in the current PDF, verbatim, including its AI-declaration sentences.
        # Compared letter by letter (tags, entities, accents and everything but A-Z removed): the shown abstract equals the
        # extracted PDF text, and that text occurs in pdftotext's output of the first pages of the PDF.
        def letters(s):
            s = unescape(re.sub(r'<[^>]+>', '', s))
            s = ''.join(ch for ch in unicodedata.normalize('NFKD', s) if not unicodedata.combining(ch))
            return re.sub(r'[^A-Za-z]', '', s)
        bad, pdf_checked = [], 0
        have_pdftotext = shutil.which('pdftotext') is not None
        for r in records:
            a = r['abstract_display']
            page = abs_html[r['id']]
            m = re.search(r'<div class="abs-abstract">(.*?)</div>', page, re.S)
            shown = re.sub(r'<span class="abs-descriptor">Abstract:</span>', '', m.group(1)) if m else ''
            if a.get('source') != 'pdf' or not a.get('pdf_text'):
                bad.append((r['id'], 'abstract not taken from the PDF'))
                continue
            if letters(shown) != letters(a['pdf_text']):
                bad.append((r['id'], 'shown abstract differs from the PDF abstract'))
            if have_pdftotext:
                out = subprocess.run(['pdftotext', '-f', '1', '-l', '3', str(ROOT / r['old_urls']['pdf_current']), '-'], capture_output=True, text=True).stdout
                pdf_checked += 1
                if letters(a['pdf_text']) not in letters(out):
                    bad.append((r['id'], 'recorded PDF abstract not found in the PDF'))
        # earlier versions: the version page shows that version's abstract, equal letter by letter to its recorded PDF text,
        # and for a served PDF the text occurs in that PDF
        for r in records:
            h = r['history_display']
            for v in h[:-1]:
                vp = ROOT / 'abs' / f'{r["id"]}v{v["v"]}' / 'index.html'
                m = re.search(r'<div class="abs-abstract">(.*?)</div>', vp.read_text(), re.S) if vp.is_file() else None
                if not v.get('abstract_html') or not m:
                    bad.append((r['id'], f'v{v["v"]}: no abstract on the version page'))
                    continue
                shown = re.sub(r'<span class="abs-descriptor">Abstract:</span>', '', m.group(1))
                if letters(shown) != letters(v.get('pdf_text', '')):
                    bad.append((r['id'], f'v{v["v"]}: abstract differs from its PDF text'))
                if have_pdftotext and v.get('pdf_path'):
                    out = subprocess.run(['pdftotext', '-f', '1', '-l', '3', str(ROOT / v['pdf_path']), '-'], capture_output=True, text=True).stdout
                    pdf_checked += 1
                    if letters(v['pdf_text']) not in letters(out):
                        bad.append((r['id'], f'v{v["v"]}: recorded abstract not found in its PDF'))
        check('6e_abstracts', not bad, problems=bad, pdf_abstracts_checked=pdf_checked if have_pdftotext else 'pdftotext not installed',
              with_ai_sentences=[r['id'] for r in records if r['abstract_display'].get('ai_sentences')])

        # ---- 7: idempotent, no placeholders
        run = subprocess.run([sys.executable, '-B', str(ROOT / 'scripts/build_liuxiv.py'), '--check'], capture_output=True, text=True)
        holders = []
        for p in [ROOT / 'abs' / r['id'] / 'index.html' for r in records] + [ROOT / f for f in old_pp] + [ROOT / B.PAGE]:
            t = p.read_text()
            for pat in (r'@@', r'\{\{', r'\}\}', r'\{r\[', r'\bNone\b', r'\bnull\b', r'\bundefined\b', r'\bNaN\b'):
                for m in re.finditer(pat, t):
                    holders.append((str(p.relative_to(ROOT)), pat, t[max(0, m.start() - 40):m.end() + 40]))
        check('7_idempotent_no_placeholders', run.returncode == 0 and not holders, generator_check=run.stdout.strip() or run.stderr.strip(), placeholders=holders)

    except Exception:  # a crash is a failed check, never a silent pass
        rep['exception'] = traceback.format_exc()
        fail.append('exception')

    rep['summary'] = {'failed': fail, 'all_ok': not fail}
    print(json.dumps(rep, indent=1, ensure_ascii=False))
    sys.exit(1 if fail else 0)


if __name__ == '__main__':
    main()
