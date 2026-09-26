#!/usr/bin/env python3
"""Rebuild homepage topic lists from the canonical data/content.json.

Run after scripts/build_content.py. No publication metadata lives in this script.
Topic membership and existing public notes/AI disclaimers also live in content.json.
Only the RESEARCH:START / RESEARCH:END block in index.html is replaced.
"""
import argparse, html, json, re
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]

def render(data):
    lines=['      <p class="research-interests">Open a topic to see my papers on it.</p>']
    for topic in data['topics']:
        key=topic['id']; items=[]
        for i,p in enumerate(data['publications']):
            if key not in p['topics']:continue
            link=next((l['url'] for l in p['links'] if 'arxiv.org' in l['url']),p['links'][0]['url'] if p['links'] else '')
            m=re.search(r'arxiv\.org/abs/(\d{4})\.',link); yymm=int(m.group(1)) if m else p['year']%100*100
            venue=p['venue_html']
            if venue.startswith('Note.'): venue=''
            items.append((yymm,0,-i,{'title':p['title_html'],'link':link,'co':[a for a in p['authors'] if a!=data['profile']['name']],'venue':venue}))
        for i,p in enumerate(data['research_extras']):
            if key not in p['topics']:continue
            if not (ROOT/p['url']).exists():raise SystemExit('Missing extra: '+p['url'])
            items.append((p['yymm'],1,-i,{'title':p['title_html'],'link':p['url'],'co':[],'venue':p['disclaimer']}))
        items=[x[3] for x in sorted(items,key=lambda x:x[:3],reverse=True)]
        lines.extend(['      <details class="topic">',f'        <summary>{html.escape(topic["label"])} <span class="topic-count">({len(items)} papers)</span></summary>','        <ul class="topic-list">'])
        for p in items:
            title=f'<a href="{html.escape(p["link"],quote=True)}">{p["title"]}</a>'
            co=' (with '+html.escape(', '.join(p['co']))+')' if p['co'] else ''
            venue=' <span class="topic-venue">'+p['venue'].rstrip('. ')+'.</span>' if p['venue'] else ''
            lines.append(f'          <li>{title}{co}.{venue}</li>')
        lines.extend(['        </ul>','      </details>'])
    return '\n'.join(lines)+'\n'

def main():
    parser=argparse.ArgumentParser();parser.add_argument('--check',action='store_true');args=parser.parse_args()
    data=json.loads((ROOT/'data/content.json').read_text());idx=ROOT/'index.html';source=idx.read_text()
    a,b='      <!-- RESEARCH:START -->\n','      <!-- RESEARCH:END -->\n'
    if a not in source or b not in source:raise SystemExit('Missing homepage Research markers')
    result=source[:source.index(a)+len(a)]+render(data)+source[source.index(b):]
    if args.check:
        if result!=source:raise SystemExit('Homepage topics are stale')
        print('Homepage topics match data/content.json.')
    else:idx.write_text(result);print('Homepage topics rebuilt from data/content.json.')
if __name__=='__main__':main()
