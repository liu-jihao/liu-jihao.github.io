#!/usr/bin/env python3
"""Build a local, public-text-only search index and static search page (stdlib)."""
from pathlib import Path
from html.parser import HTMLParser
from html import escape
import argparse, hashlib, json, re
ROOT=Path(__file__).resolve().parents[1]
VOID={'area','base','br','col','embed','hr','img','input','link','meta','param','source','track','wbr'}
class Node:
    def __init__(self,tag='',attrs=None,start=0,raw=''):
        self.tag=tag;self.attrs=dict(attrs or []);self.children=[];self.start=start;self.raw=raw
    def all(self,predicate):
        result=[]
        for child in self.children:
            if isinstance(child,Node):
                if predicate(child):result.append(child)
                result.extend(child.all(predicate))
        return result
    def text(self):
        return re.sub(r'\s+',' ',''.join(c.text()+' ' if isinstance(c,Node) else c for c in self.children)).strip()
    def cls(self,name):return name in self.attrs.get('class','').split()
    def first(self,name):
        return next(iter(self.all(lambda n:n.cls(name))),None)
class Tree(HTMLParser):
    def __init__(self,text):
        super().__init__(convert_charrefs=True);self.root=Node();self.stack=[self.root]
        self.starts=[0]
        for line in text.splitlines(keepends=True):self.starts.append(self.starts[-1]+len(line))
        self.feed(text)
    def handle_starttag(self,tag,attrs):
        line,col=self.getpos();node=Node(tag,attrs,self.starts[line-1]+col,self.get_starttag_text());self.stack[-1].children.append(node)
        if tag not in VOID:self.stack.append(node)
    def handle_startendtag(self,tag,attrs):self.handle_starttag(tag,attrs);self.handle_endtag(tag)
    def handle_endtag(self,tag):
        for i in range(len(self.stack)-1,0,-1):
            if self.stack[i].tag==tag:self.stack=self.stack[:i];break
    def handle_data(self,data):self.stack[-1].children.append(data)
def plain(text):return Tree(text).root.text()
def year(text):
    found=re.findall(r'\b(?:19|20)\d{2}\b',text)
    return found[0] if found else ''
def text_of(node,cls):
    target=node.first(cls);return target.text() if target else ''
def build(check=False):
    data=json.loads((ROOT/'data/content.json').read_text());records=[];results={}
    def add(id,title,text,url,category,year_value='',**extra):
        records.append(dict(id=id,title=plain(title),text=plain(text),url=url,category=category,year=str(year_value),**extra))
    for p in data['publications']:
        links=' '.join(l['label']+' '+l['url'] for l in p['links'])
        add(p['id'],p['title_html'],p['number']+' '+', '.join(p['authors'])+' '+p['venue_html']+' '+links,'publications.html#'+p['id'],'publications',p['year'],status=p['cv']['section'],topics=p['topics'],number=p['number'])
    for group in data['teaching']:
        for row in group['items']:
            add(row['id'],row['title_html'],group['heading']+' '+row['meta_html']+' '+' '.join(r['label'] for r in row.get('resources',[])),'teaching.html#'+row['id'],'teaching',year(row['meta_html']))
    for group in data['talks']:
        for row in group['items']:
            add(row['id'],row['title_html'],group['heading']+' '+row['meta_html']+' '+' '.join(r['label'] for r in row.get('resources',[])),'talks.html#'+row['id'],'talks',year(row['meta_html']))
    for row in data['organized']:
        if row.get('cv_only'):continue
        # Index only entries rendered on the public activity page.
        if f'id="{row["id"]}"' not in (ROOT/'conferences.html').read_text():continue
        add('organized-'+row['id'],row['title_html'],row['date_html']+' '+' '.join(row['details_html']),'conferences.html#'+row['id'],'activities',year(row['date_html']))
    structured={'publications.html','teaching.html','talks.html','conferences.html'}
    for page in sorted(ROOT.glob('*.html')):
        if page.name in structured|{'search.html','404.html'}:continue
        source=page.read_text();tree=Tree(source).root;main=next(iter(tree.all(lambda n:n.tag=='main')),None)
        if not main:continue
        replacements=[]
        if page.name.startswith('pku-ag-seminar'):
            pageyear=year(page.name) or '2026'
            for entry in main.all(lambda n:n.cls('seminar-entry')):
                date=text_of(entry,'seminar-date');speaker=text_of(entry,'seminar-speaker');title=text_of(entry,'seminar-title') or speaker or date
                stable=entry.attrs.get('id') or 'seminar-'+hashlib.sha256((date+'|'+speaker+'|'+title).encode()).hexdigest()[:12]
                if not entry.attrs.get('id'):replacements.append((entry.start,entry.raw,entry.raw[:-1]+f' id="{stable}">'))
                if len(entry.all(lambda n:n.cls('seminar-title')))>1: title=date+' '+pageyear+' · '+(speaker or 'Seminar program')
                add(page.name+'-'+stable,title,entry.text()+' '+pageyear,page.name+'#'+stable,'seminars',pageyear)
        elif page.name in ('notes.html','ai-results.html'):
            for index,entry in enumerate(main.all(lambda n:n.cls('entry')),1):
                title=text_of(entry,'entry-title');number=text_of(entry,'entry-number').strip('. ')
                stable=entry.attrs.get('id') or ('note-' if page.name=='notes.html' else 'paper-')+(number or str(index))
                if not entry.attrs.get('id'):replacements.append((entry.start,entry.raw,entry.raw[:-1]+f' id="{stable}">'))
                add(page.name+'-'+stable,title,entry.text(),page.name+'#'+stable,'notes' if page.name=='notes.html' else 'ai',year(text_of(entry,'entry-venue')))
        else:
            heading=next(iter(main.all(lambda n:n.tag=='h1')),None)
            text=main.text()
            if page.name=='index.html':text=' '.join(n.text() for n in main.all(lambda n:n.tag=='p')[:2])+' '+text_of(main,'contact-block')
            add('page-'+page.stem,heading.text() if heading else page.stem,text,page.name,'pages',year(text) if page.name=='ffm-conference.html' else '')
        for start,old,new in sorted(replacements,reverse=True):source=source[:start]+new+source[start+len(old):]
        if replacements:results[page]=source
    # The published index contains only rendered text and public navigation metadata, never source notes/confirmation fields.
    payload={'version':1,'updated':data['updated'],'topics':data['topics'],'records':records}
    results[ROOT/'assets/search-index.js']='/* Generated from public website text. No analytics or remote requests. */\nwindow.SITE_SEARCH = '+json.dumps(payload,ensure_ascii=False,separators=(',',':')).replace('<','\\u003c')+';\n'
    template=(ROOT/'templates/search.html').read_text()
    from build_chrome import render_page
    results[ROOT/'search.html']=render_page(template,'search.html',json.loads((ROOT/'data/navigation.json').read_text()))
    stale=[str(p.relative_to(ROOT)) for p,s in results.items() if not p.exists() or p.read_text()!=s]
    if check:
        if stale:raise SystemExit('Discovery outputs are stale: '+', '.join(stale))
        print(f'Search index matches {len(records)} public records.')
    else:
        for p,s in results.items():p.write_text(s)
        print(f'Built search index: {len(records)} public records.')
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--check',action='store_true');args=p.parse_args();build(args.check)
