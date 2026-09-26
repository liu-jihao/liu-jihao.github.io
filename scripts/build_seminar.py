#!/usr/bin/env python3
"""Render current seminar dates against the explicit shared content review date."""
from pathlib import Path
from datetime import date
from html import escape
import argparse,json,re
ROOT=Path(__file__).resolve().parents[1]
def build(check=False):
    data=json.loads((ROOT/'data/seminar-current.json').read_text());asof=json.loads((ROOT/'data/content.json').read_text())['updated']
    groups=[('upcoming','Upcoming',sorted((r for r in data['records'] if r['date']>=asof),key=lambda r:r['date'])),('past-talks','Past talks this semester',sorted((r for r in data['records'] if r['date']<asof),key=lambda r:r['date'],reverse=True))]
    lines=[]
    for key,label,records in groups:
        lines.append(f'      <h2 id="{key}">{label}</h2>')
        if not records:lines.append('      <p>For further dates, see the <a href="https://bicmr.pku.edu.cn/content/lists/35.html">BICMR seminar announcements</a>.</p>')
        for r in records:
            d=date.fromisoformat(r['date']);display=f'{d:%B} {d.day}, {d:%Y}'
            lines.extend([f'      <div class="seminar-entry" id="{r["id"]}">',f'        <div class="seminar-date"><time datetime="{r["date"]}">{display}</time> · {r["time"]} (Beijing time)</div>',r['body_html'],f'        <p class="course-meta">{escape(r["venue"])} · <a href="{escape(r["source_url"],quote=True)}">Official announcement</a></p>' if 'bicmr.pku.edu.cn/content/show/' in r['source_url'] else f'        <p class="course-meta">{escape(r["venue"])}</p>'])
            if r.get('note'):lines.append(f'        <p class="schedule-note">{escape(r["note"])}</p>')
            lines.append('      </div>')
    page=ROOT/'pku-ag-seminar.html';s=page.read_text();a=s.index('      <h2 id="upcoming">');b=s.index('    </main>',a)
    result=s[:a]+'\n\n'.join(lines)+'\n'+s[b:]
    if check:
        if result!=s:raise SystemExit('Current seminar page is stale')
        print('Current seminar matches its source and explicit review date.')
    else:page.write_text(result)
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--check',action='store_true');args=p.parse_args();build(args.check)
