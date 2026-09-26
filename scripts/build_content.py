#!/usr/bin/env python3
"""Generate the existing static pages and self-contained CV from data/content.json.

Run from any directory: python3 scripts/build_content.py [--check]
No network, packages, JavaScript runtime, or remote writes are required.
Only page <main> interiors are replaced; shared navigation/footer attributes are preserved.
"""
from pathlib import Path
from html.parser import HTMLParser
from html import escape, unescape
from datetime import date
import argparse, hashlib, json, re

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / 'data/content.json'


def plain(s):
    return unescape(re.sub('<[^>]+>', '', s)).strip()


def escape_tex(s):
    specials = {'\\':r'\textbackslash{}','&':r'\&','%':r'\%','$':r'\$','#':r'\#','_':r'\_','{':r'\{','}':r'\}','~':r'\textasciitilde{}','^':r'\textasciicircum{}',
                '–':'--','—':'---','−':'-','\u00a0':'~','ℝ':r'$\mathbb{R}$','ε':r'$\varepsilon$','δ':r'$\delta$','ν':r'$\nu$','≤':r'$\le$'}
    s = ''.join(specials.get(c,c) for c in s)
    return re.sub(r'([\u3400-\u9fff]+)',r'{\\cnfont \1}',s)


class TexHTML(HTMLParser):
    def __init__(self):
        super().__init__(convert_charrefs=True)
        self.parts=[]; self.stack=[]
    def handle_data(self,s): self.parts.append(escape_tex(s))
    def handle_starttag(self,tag,attrs):
        attrs=dict(attrs)
        if tag=='a':
            url=attrs['href']
            if not re.match(r'https?://|mailto:',url): url='https://jihaoliu.org/'+url
            self.parts.append(r'\href{'+url.replace('%',r'\%').replace('#',r'\#')+'}{'); self.stack.append((tag,'}'))
        elif tag in ('em','i','b','strong'):
            self.parts.append(r'\emph{' if tag in ('em','i') else r'\textbf{'); self.stack.append((tag,'}'))
        elif tag in ('sup','sub'):
            self.parts.append(r'\textsuperscript{' if tag=='sup' else r'\textsubscript{'); self.stack.append((tag,'}'))
        elif tag=='br': self.parts.append(r'\\ ')
    def handle_endtag(self,tag):
        if self.stack and self.stack[-1][0]==tag: self.parts.append(self.stack.pop()[1])


def latex(s):
    parser=TexHTML(); parser.feed(s)
    return re.sub(r'\s+', ' ', ''.join(parser.parts)).strip()


def jump(items,label):
    links=''.join(f'<a href="#{escape(i)}">{escape(t)}</a>' for i,t in items)
    return f'      <nav class="section-jump" aria-label="{escape(label)}">{links}</nav>\n'


def render_publications(data):
    pubs=data['publications']; years=list(dict.fromkeys(p['year'] for p in pubs))
    lines=['      <h1>Publications and Preprints</h1>',
           '      <p class="page-intro">In reverse chronological order by arXiv posting date. Related versions retain their shared publication number.</p>',
           jump([(f'publications-{year}',str(year)) for year in years],'Publication years').rstrip()]
    current=None
    for p in pubs:
        if p['year']!=current:
            if current is not None: lines.append('      </section>')
            current=p['year']; lines.extend([f'      <section aria-labelledby="publications-{current}">',f'      <h2 class="year-heading" id="publications-{current}">{current}</h2>'])
        titleid=p['id']+'-title'
        lines.extend([f'      <article class="entry" id="{p["id"]}" aria-labelledby="{titleid}">',
                      f'        <h3 class="entry-title" id="{titleid}"><span class="entry-number">{p["number"]}.</span> {p["title_html"]}</h3>',
                      f'        <div class="entry-authors">{escape(", ".join(p["authors"]))}</div>'])
        if p['venue_html']: lines.append(f'        <div class="entry-venue">{p["venue_html"]}</div>')
        anchors=' '.join(f'<a href="{escape(l["url"],quote=True)}">{escape(l["label"])}</a>' for l in p['links'])
        anchors+=f' <a class="entry-permalink" href="#{p["id"]}" aria-label="Permanent link to publication {p["number"]}">Permalink</a>'
        lines.extend([f'        <div class="entry-links">{anchors}</div>','      </article>'])
    lines.append('      </section>')
    return '\n'.join(lines)+'\n'


def render_activities(groups,title,numbered=False):
    lines=[f'      <h1>{title}</h1>','      <p class="page-intro">In reverse chronological order within each section.</p>',jump([(g['id'],g['heading']) for g in groups],title+' sections').rstrip()]
    for g in groups:
        lines.extend([f'      <section aria-labelledby="{g["id"]}">',f'      <h2 class="year-heading" id="{g["id"]}">{escape(g["heading"])}</h2>',f'      <ul class="{"talks-list" if numbered else "teaching-list"}">'])
        for item in g['items']:
            prefix=f'<span class="entry-number">{item["number"]}.</span> ' if numbered else ''
            meta=item['meta_html']
            context=escape('Event information: '+plain(item['title_html']),quote=True)
            meta=re.sub(r'<a (href="[^"]+")>link</a>', lambda m:f'<a {m.group(1)} aria-label="{context}">link</a>',meta)
            lines.extend([f'        <li id="{item["id"]}">',f'          <div class="course-name">{prefix}{item["title_html"]}</div>',f'          <div class="course-meta">{meta}</div>',
                          *(['          <div class="course-meta">Historical link unavailable.</div>'] if item.get('unavailable_links') and not item.get('link_preservation') else []),
                          *(["          <div class=\"entry-links\">"+' · '.join(f'<a href="{escape(r["url"],quote=True)}">{escape(r["label"])}</a>' for r in item['resources'])+'</div>'] if item.get('resources') else []),'        </li>'])
        lines.extend(['      </ul>','      </section>'])
    return '\n'.join(lines)+'\n'


def organized_details(event,data):
    parts=list(event['details_html'])
    if event.get('coorganizers_profile_key'):
        parts.append('With '+', '.join(data['profile'][event['coorganizers_profile_key']]))
    return parts


def render_organized(data):
    lines=['      <h1>Conferences and Seminars</h1>',
           '      <p class="page-intro">Conferences and seminar series I have co-organized, including ongoing series and past events.</p>',
           '      <ul class="conferences-list">']
    for q in data['organized']:
        if not q.get('show_on_web',True):continue
        lines.extend([f'        <li id="organized-{q["id"]}">',f'          <div class="course-name">{q["title_html"]}</div>',f'          <div class="course-meta">'+ ' · '.join([q['date_html']]+organized_details(q,data))+'</div>'])
        for extra in q['website_extra_html']:lines.append(f'          <div class="course-meta">{extra}</div>')
        if q.get('unavailable_links'):lines.append('          <div class="course-meta">Historical link unavailable.</div>')
        lines.append('        </li>')
    lines.append('      </ul>')
    return '\n'.join(lines)+'\n'

def render_contact(data):
    p=data['profile']; d=date.fromisoformat(data['updated']); updated=f'{d.day} {d:%B %Y}'
    return f'''      <div class="contact-block contact-grid">
        <div class="contact-label">Email:</div>
        <div class="contact-value">{escape(p['email_display'])}</div>

        <div class="contact-label">Office:</div>
        <div class="contact-value">
          <span class="office-en">{escape(p['office_en']).replace(' 78402','&nbsp;78402')}</span><br>
          <span class="office-cn">{escape(p['office_zh']).replace(' 78402','&nbsp;78402')}</span>
        </div>

        <div class="contact-label">CV:</div>
        <div class="contact-value">
          <a href="cv/cv.pdf">cv.pdf</a>
          &nbsp;<span style="color: var(--color-text-muted); font-size: 0.9em;">(updated {updated})</span>
        </div>

        <div class="contact-label">Google Scholar:</div>
        <div class="contact-value">
          <a href="{escape(p['scholar_url'],quote=True)}">profile</a>
          &nbsp;<span style="color: var(--color-text-muted); font-size: 0.9em;">{escape(p.get('scholar_citation_label',''))}</span>
        </div>
      </div>
'''

def pub_latex(p,data):
    cv=p['cv']; title=latex(p['title_html']).rstrip('.')+'.'
    authors=', '.join(r'\me' if a=='Jihao Liu' else escape_tex(a) for a in p['authors'])+'.'
    venue=p['venue_html']
    # The two numbers are section-local in the CV, but website numbers stay stable.
    if p['number'] in ('60a','60b'):
        other_number='60b' if p['number']=='60a' else '60a'
        other=next(q['cv']['label'] for q in data['publications'] if q['number']==other_number)
        venue=re.sub(r'no\.&nbsp;60[ab]', 'no. '+other,venue)
    parts=[r'\paper{'+title+'}',authors]
    if venue: parts.append(latex(venue).rstrip('.')+'.')
    anchors=[]
    for l in p['links']:
        label=l['label']
        if label.lower()=='doi':
            label='doi:'+re.sub(r'https?://(?:dx\.)?doi\.org/','',l['url'])
            url=l['url'].replace('%',r'\%').replace('#',r'\#')
            anchors.append(r'\idlink{'+url+'}{'+label+'}')
        else:
            anchors.append(latex(f'<a href="{escape(l["url"],quote=True)}">{escape(label)}</a>'))
    if anchors: parts.append(', '.join(anchors)+'.')
    return r'\pubitem{'+cv['label']+'}{'+str(cv['year'])+'}{%\n  '+'\n  '.join(parts)+'}\n'


def render_cv(data):
    p=data['profile']; d=date.fromisoformat(data['updated']); updated=f'{d.day} {d:%B %Y}'
    lines=[r'''% ---------- Heading ----------
\noindent
\begin{minipage}[t]{0.66\textwidth}
  \vspace{0pt}
  {\fontsize{30pt}{34pt}\selectfont\color{textcol}'''+escape_tex(p['name'])+r'''}\;%
  {\fontsize{20pt}{24pt}\selectfont\color{muted}\cnfont '''+p['name_zh']+r'''}\\[8pt]
  {\normalsize\color{accent}\scshape Curriculum Vitae}\\[3pt]
  {\footnotesize\color{muted}Last updated '''+updated+r'''}
\end{minipage}%
\hfill%
\begin{minipage}[t]{0.32\textwidth}
  \vspace{6pt}\raggedleft\small
  '''+'\\\\\n  '.join(escape_tex(s) for s in p['affiliation_lines'])+r'''\\[4pt]
  \href{mailto:'''+p['email']+'}{'+escape_tex(p['email'])+r'''}\\
  \href{'''+p['website']+'}{'+escape_tex(p['website'].replace('https://',''))+r'''}
\end{minipage}\par
\vspace{10pt}{\color{rule}\hrule height 0.4pt}\par\vspace{14pt}

\cvsection{Research Interests}{research-interests}
'''+escape_tex(p['research_interests'])]
    for key,title in [('employment','Employment'),('education','Education')]:
        lines.append(r'\cvsection{'+title+'}{'+key+'}')
        for e in p[key]:
            lines.append(r'\entry{'+escape_tex(e['dates'])+r'}{\textbf{'+escape_tex(e['institution'])+'}, '+escape_tex(e['unit'])+r'\\ '+escape_tex(e['details'])+'}')
    for key,title,note in [('published','Published Papers','Listed by year of journal publication, most recent first.'),('accepted','Papers to Appear','Accepted but not yet published. Listed by arXiv posting date, most recent first.'),('preprint','Preprints','Listed by arXiv posting date, most recent first; related versions are grouped together.')]:
        entries=sorted((p for p in data['publications'] if p['cv']['section']==key),key=lambda p:(p['cv']['order'],p['number']))
        if not entries: continue
        labels=[p['cv']['label'] for p in entries]
        assert len(labels)==len(set(labels)), f'Duplicate CV labels in {key}'
        lines.append(r'\cvsection{'+title+'}{'+key+'}\n'+r'{\small\color{muted}\noindent '+note+r'}\par\vspace{4pt}')
        lines.extend(pub_latex(p,data) for p in entries)
    lines.append(r'\cvsection{Conference and Workshop Talks}{talks}'+'\n'+r'{\small\color{muted}\noindent Selected; seminar talks omitted.}\par\vspace{4pt}')
    for row in data['talks'][0]['items']:
        title=row['title_html']
        if row.get('cv_title_language')=='english_parentheses': title=plain(title).split('(',1)[1].rsplit(')',1)[0]
        segments=re.split(r'\s*·\s*',row['meta_html']); when=plain(segments[0]); year=re.search(r'\d{4}',when).group(0); month=when.split(' ')[0][:3]
        # Put the existing event link around its name rather than an isolated "link" label.
        rest=segments[1:]; am=re.search(r'<a href="([^"]+)">(.*?)</a>',rest[-1]) if rest else None
        if am:
            rest=rest[:-1]
            if rest: rest[0]=f'<a href="{am.group(1)}">{rest[0]}</a>'
        lines.append(r'\entry{'+month+' '+year+r'}{\emph{'+latex(title).rstrip('.')+r'.}\\ '+latex(', '.join(rest)).rstrip('.')+'.}')
    lines.append(r'\cvsection{Conferences and Seminars Organized}{organized}')
    for q in sorted(data['organized'],key=lambda q:q['cv_order']):
        details='. '.join(q.get('cv_details_html',organized_details(q,data)))
        lines.append(r'\entry{'+escape_tex(q['cv_date'])+r'}{'+latex(q['title_html']).rstrip('.')+r'.\\ '+latex(q['date_html']).rstrip('.')+'. '+latex(details).rstrip('.')+'.}')
    if data['teaching']:
        first_group_space=min(24,2*len(data['teaching'][0]['items'])+5)
        lines.append(r'\Needspace{'+str(first_group_space)+r'\baselineskip}')
    lines.append(r'\cvsection{Teaching}{teaching}')
    for group_number,group in enumerate(data['teaching'],1):
        # Keep the current short institutional lists together when possible,
        # without making a future list longer than a page unbreakable.
        needed=min(24,2*len(group['items'])+2)
        lines.append(r'\Needspace{'+str(needed)+r'\baselineskip}')
        lines.append(r'\cvsubsection{'+escape_tex(group['heading'])+'}{teaching-'+str(group_number)+'}')
        # Merge repeated Utah courses in the CV while keeping all terms in the source/site.
        rows=group['items']
        if group['heading']=='University of Utah':
            merged={}
            for row in rows:
                meta=plain(row['meta_html']).split(' · '); key=(row['title_html'],meta[1]); merged.setdefault(key,[]).append(meta[0])
            for (title,role),terms in merged.items():
                lines.append(r'\entry{'+escape_tex(', '.join(terms))+r'}{'+latex(title)+r'~---~'+escape_tex(role.lower())+'.}')
        else:
            for row in rows:
                meta=re.split(r'\s*·\s*',row['meta_html'])
                lines.append(r'\entry{'+latex(meta[0])+r'}{'+latex(row['title_html'])+r'~---~'+latex('; '.join(meta[1:])).rstrip('.')+'.}')
    lines.append(data['cv_only_sections']['refereeing_latex'].replace(r'\section*{Refereeing}',r'\cvsection{Refereeing}{refereeing}'))
    template=(ROOT/'cv/template.tex').read_text()
    rendered=template.replace('@@CONTENT@@','\n\n'.join(lines))
    return rendered.replace('@@PDF_ID@@',hashlib.sha256(rendered.encode()).hexdigest()[:32])


def main():
    parser=argparse.ArgumentParser(); parser.add_argument('--check',action='store_true')
    parser.add_argument('--cv-only',action='store_true',help='Regenerate/check only cv/cv.tex, without touching website pages')
    args=parser.parse_args()
    data=json.loads(DATA.read_text())
    ids=[p['id'] for p in data['publications']]
    assert len(ids)==len(set(ids)), 'Duplicate publication IDs'
    assert {p['cv']['section'] for p in data['publications']} <= {'published','accepted','preprint'}
    results={}
    if not args.cv_only:
        for file,content in [('publications.html',render_publications(data)),('teaching.html',render_activities(data['teaching'],'Teaching')),('talks.html',render_activities(data['talks'],'Invited Talks',True)),('conferences.html',render_organized(data))]:
            source=(ROOT/file).read_text()
            results[ROOT/file]=re.sub(r'(<main\b[^>]*>).*?(</main>)',lambda m:m.group(1)+'\n'+content+'    '+m.group(2),source,count=1,flags=re.S)
    results[ROOT/'cv/cv.tex']=render_cv(data)
    if not args.cv_only:
        index=(ROOT/'index.html').read_text()
        a,b='      <!-- CONTACT:START -->\n','      <!-- CONTACT:END -->\n'
        if a not in index or b not in index: raise SystemExit('Missing homepage Contact markers')
        results[ROOT/'index.html']=index[:index.index(a)+len(a)]+render_contact(data)+index[index.index(b):]
    stale=[str(path.relative_to(ROOT)) for path,text in results.items() if not path.exists() or path.read_text()!=text]
    if args.check:
        if stale: raise SystemExit('Generated files are stale: '+', '.join(stale))
        print(('CV matches' if args.cv_only else 'Generated pages and CV match')+' data/content.json.')
    else:
        for path,text in results.items(): path.write_text(text)
        print('Built cv/cv.tex from shared content.' if args.cv_only else f'Built {len(data["publications"])} publication records, teaching, talks, and cv/cv.tex.')

if __name__=='__main__': main()
