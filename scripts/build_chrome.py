#!/usr/bin/env python3
"""Regenerate shared navigation and local progressive-enhancement scripts."""
from pathlib import Path
import argparse,html,json,re
ROOT=Path(__file__).resolve().parents[1]
def render_page(source,name,data):
    def link(pair):
        href,label=pair
        active=href==name or (href=='pku-ag-seminar.html' and name.startswith('pku-ag-seminar-')) or (href=='ai-results.html' and name=='ai-results-about.html')
        attrs=' class="active" aria-current="page"' if active and href==name else ' class="active"' if active else ''
        return f'          <a href="{html.escape(href,quote=True)}"{attrs}>{html.escape(label)}</a>'
    nav=['      <nav class="main-nav" aria-label="Main navigation">','        <div class="nav-primary">']
    nav.extend(map(link,data['primary']))
    nav+=['        </div>','        <details class="nav-more">','          <summary>Research &amp; community</summary>','          <div class="nav-secondary">']
    nav.extend(map(link,data['secondary']))
    nav+=['          </div>','        </details>','      </nav>']
    source=re.sub(r'      <nav class="main-nav".*?</nav>','\n'.join(nav),source,count=1,flags=re.S)
    source=re.sub(r'  <script[^>]+src="assets/(?:search-index|discovery)\.js"[^>]*></script>\n?', '',source)
    scripts=('  <script src="assets/search-index.js" defer></script>\n' if name in ('publications.html','search.html') else '')+'  <script src="assets/discovery.js" defer></script>\n'
    return source.replace('</head>',scripts+'</head>').replace('<main>','<main id="main-content" tabindex="-1">')
def build(check=False):
    data=json.loads((ROOT/'data/navigation.json').read_text());stale=[]
    for page in ROOT.glob('*.html'):
        source=render_page(page.read_text(),page.name,data)
        if source!=page.read_text():stale.append(page.name)
        if not check:page.write_text(source)
    if check and stale:raise SystemExit('Shared navigation is stale: '+', '.join(stale))
    if check:print('Shared navigation matches data/navigation.json.')
if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--check',action='store_true');args=parser.parse_args();build(args.check)
