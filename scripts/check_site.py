#!/usr/bin/env python3
"""Local structural and link checks for all root website pages (stdlib only)."""
from pathlib import Path
from html.parser import HTMLParser
from urllib.parse import urlsplit,unquote
import json,sys,re
ROOT=Path(__file__).resolve().parents[1]
class Page(HTMLParser):
 def __init__(self,text):
  super().__init__();self.ids=[];self.urls=[];self.tags=[];self.feed(text)
 def handle_starttag(self,tag,attrs):
  a=dict(attrs);self.tags.append(tag)
  if 'id' in a:self.ids.append(a['id'])
  for k in ('href','src'):
   if a.get(k):self.urls.append(a[k])
pages={p:Page(p.read_text()) for p in ROOT.glob('*.html')}; errors=[]; links=0
for path,page in pages.items():
 for tag in ('h1','main','title'):
  if page.tags.count(tag)!=1:errors.append(f'{path.name}: expected one {tag}')
 if len(page.ids)!=len(set(page.ids)):errors.append(f'{path.name}: duplicate IDs')
 if 'main-content' not in page.ids:errors.append(f'{path.name}: missing main landmark ID')
 for url in page.urls:
  u=urlsplit(url)
  if u.scheme or u.netloc:continue
  target=(ROOT/u.path.lstrip('/') if u.path.startswith('/') else path.parent/unquote(u.path)).resolve() if u.path else path
  if target.is_dir():target=target/'index.html'
  links+=1
  if not target.is_file():errors.append(f'{path.name}: missing {url}');continue
  if u.fragment and target.suffix=='.html':
   other=pages.get(target) or Page(target.read_text())
   if unquote(u.fragment) not in other.ids:errors.append(f'{path.name}: missing fragment {url}')
# Embedded font URLs need the same check as HTML references.
for css in (ROOT/'assets').glob('*.css'):
 for raw in re.findall(r'url\([\"\']?([^\)\"\']+)',css.read_text()):
  if urlsplit(raw).scheme:continue
  if not (css.parent/raw).is_file():errors.append(f'{css.name}: missing {raw}')
data=json.loads((ROOT/'data/content.json').read_text());publications=data['publications'];ids=[p['id'] for p in publications]
if len(ids)!=len(set(ids)):errors.append('Duplicate publication data IDs')
known={p['number'] for p in publications}
if not {'60a','60b','44a','44b'}<=known:errors.append('Paired publications lost')
for p in publications:
 if p['id'] not in pages[ROOT/'publications.html'].ids:errors.append('Unrendered publication '+p['id'])
index=json.loads((ROOT/'assets/search-index.js').read_text().split('=',1)[1].strip().rstrip(';'))
records=index['records']; search_ids=[r['id'] for r in records]
if len(search_ids)!=len(set(search_ids)):errors.append('Duplicate search record IDs')
for record in records:
 u=urlsplit(record['url']);target=ROOT/u.path
 if u.scheme or u.netloc or target not in pages:errors.append('Nonlocal or missing search result '+record['url']);continue
 if u.fragment and unquote(u.fragment) not in pages[target].ids:errors.append('Missing search result anchor '+record['url'])
result={'pages':len(pages),'local_links_checked':links,'publication_records':len(publications),'search_records':len(records),'errors':errors}
print(json.dumps(result,indent=2))
if errors:sys.exit(1)
