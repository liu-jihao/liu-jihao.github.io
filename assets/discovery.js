/* Local search and progressive enhancement. No requests, tracking or storage. */
(() => {
  'use strict';
  const catalogue = window.SITE_SEARCH;
  const normal = value => String(value || '').normalize('NFKD').replace(/[\u0300-\u036f]/g, '').toLowerCase().replace(/ł/g, 'l').replace(/[^\p{L}\p{N}]+/gu, ' ').trim();
  const tokens = value => normal(value).split(/\s+/).filter(Boolean);
  const matches = (record, query) => tokens(query).every(t => normal(record.title + ' ' + record.text + ' ' + (record.number || '')).includes(t));
  const categories = {publications:'Publications & preprints',teaching:'Teaching',talks:'Invited talks',activities:'Organized events',seminars:'PKU seminar archive',notes:'Notes',ai:'AI results — not human-verified',pages:'Other pages'};
  const statuses = {published:'Published',accepted:'To appear',preprint:'Preprint'};
  const element = (tag, attrs = {}, text = '') => {
    const el = document.createElement(tag);
    Object.entries(attrs).forEach(([k,v]) => el.setAttribute(k, v));
    if (text) el.textContent = text;
    return el;
  };
  function inputField(form, label, id, type='search') {
    const field=element('div',{class:'filter-field'}), input=element('input',{type,id,name:id,autocomplete:'off'});
    field.append(element('label',{for:id},label),input);form.append(field);return input;
  }
  function selectField(form,label,id,options,empty) {
    const field=element('div',{class:'filter-field'}), select=element('select',{id,name:id});
    select.append(element('option',{value:''},empty));
    for (const [value,name] of options) select.append(element('option',{value},name));
    field.append(element('label',{for:id},label),select);form.append(field);return select;
  }
  function state(form,controls,keys,render) {
    let timer;
    const fromURL=()=>{
      const params=new URLSearchParams(location.search);
      for (const [key,el] of Object.entries(controls)) {
        const value=params.get(key)||'';
        el.value=el.tagName==='SELECT' && ![...el.options].some(o=>o.value===value) ? '' : value;
      }
      render();
    };
    const update=()=>{
      const url=new URL(location.href);
      keys.forEach(key=>controls[key].value ? url.searchParams.set(key,controls[key].value) : url.searchParams.delete(key));
      history.replaceState(null,'',url.pathname+url.search+url.hash);render();
    };
    form.addEventListener('submit',e=>{e.preventDefault();clearTimeout(timer);update();});
    form.addEventListener('input',()=>{clearTimeout(timer);timer=setTimeout(update,100);});
    form.addEventListener('change',()=>{clearTimeout(timer);update();});
    form.addEventListener('reset',()=>{clearTimeout(timer);setTimeout(update,0);});
    window.addEventListener('popstate',fromURL);fromURL();
    document.addEventListener('keydown',e=>{
      const editing=/INPUT|TEXTAREA|SELECT/.test(document.activeElement.tagName)||document.activeElement.isContentEditable;
      if(e.key==='/'&&!editing&&!e.ctrlKey&&!e.metaKey&&!e.altKey&&!e.isComposing){e.preventDefault();controls.q.focus();}
    });
  }
  if(catalogue && location.pathname.endsWith('/publications.html')) {
    const records=catalogue.records.filter(r=>r.category==='publications');
    const numberedWorks=new Set(records.map(r=>String(r.number).replace(/[a-z]+$/i,''))).size;
    const form=element('form',{class:'discovery-form publications-filter',role:'search','aria-label':'Filter publications'});
    const q=inputField(form,'Title, author or arXiv / DOI','publication-query');q.setAttribute('aria-keyshortcuts','/');
    const year=selectField(form,'arXiv year','publication-year',[...new Set(records.map(r=>r.year))].sort().reverse().map(y=>[y,y]),'All years');
    const status=selectField(form,'Status','publication-status',Object.entries(statuses),'All statuses');
    const topic=selectField(form,'Research topic','publication-topic',catalogue.topics.map(t=>[t.id,t.label]),'All topics');
    form.append(element('button',{type:'reset',class:'reset-filters'},'Clear filters'));
    const result=element('p',{class:'filter-status',role:'status','aria-live':'polite','aria-atomic':'true'});
    const empty=element('p',{class:'empty-results',hidden:''},'No matching papers. Try fewer terms or clear the filters.');
    const container=element('div',{class:'discovery-panel'});container.append(form,result,empty);
    const yearNav=document.querySelector('nav.section-jump');yearNav.before(container);
    const entries=new Map(records.map(r=>[r.id,document.getElementById(r.id)]));
    state(form,{q,year,status,topic},['q','year','status','topic'],()=>{
      let count=0;
      for(const r of records){const visible=matches(r,q.value)&&(!year.value||r.year===year.value)&&(!status.value||r.status===status.value)&&(!topic.value||r.topics.includes(topic.value));entries.get(r.id).hidden=!visible;if(visible)count++;}
      for(const section of document.querySelectorAll('main > section'))section.hidden=![...section.querySelectorAll('article.entry')].some(e=>!e.hidden);
      for(const a of yearNav.querySelectorAll('a'))a.hidden=document.querySelector(a.getAttribute('href')).closest('section').hidden;
      yearNav.hidden=count===0;empty.hidden=count!==0;
      const chosen=[q.value&&`“${q.value}”`,year.value,statuses[status.value],topic.value&&catalogue.topics.find(t=>t.id===topic.value)?.label].filter(Boolean);
      result.textContent=`Showing ${count} of ${records.length} records${chosen.length?' · '+chosen.join(' · '):' ('+numberedWorks+' numbered works)'}.`;
    });
  }
  const searchHost=document.getElementById('site-search');
  if(catalogue && searchHost) {
    const form=element('form',{class:'discovery-form site-search-form',role:'search','aria-label':'Search this site'});
    const q=inputField(form,'Search terms','site-query');q.setAttribute('aria-keyshortcuts','/');
    const category=selectField(form,'Section','site-category',Object.entries(categories),'All sections');
    const year=selectField(form,'Year','site-year',[...new Set(catalogue.records.map(r=>r.year).filter(Boolean))].sort().reverse().map(y=>[y,y]),'All years');
    form.append(element('button',{type:'reset',class:'reset-filters'},'Clear search'));
    const status=element('p',{class:'filter-status',role:'status','aria-live':'polite','aria-atomic':'true'});
    const list=element('ol',{class:'search-results','aria-label':'Search results'}),empty=element('p',{class:'empty-results',hidden:''},'No matches. Try another name, title, course or year.');
    const more=element('button',{type:'button',class:'load-results',hidden:''},'Show more results');
    let results=[],limit=20;
    function draw(append=false){
      const start=append?list.children.length:0;if(!append)list.replaceChildren();
      for(const r of results.slice(start,limit)){
        const li=element('li'),h=element('h2',{class:'search-result-title'}),a=element('a',{href:r.url},r.title);
        h.append(a);li.append(h,element('p',{class:'result-meta'},categories[r.category]+(r.year?' · '+r.year:'')+(r.status?' · '+statuses[r.status]:'')));
        const excerpt=r.text.length>240?r.text.slice(0,237)+'…':r.text;
        li.append(element('p',{class:'result-excerpt'},excerpt));list.append(li);
      }
      empty.hidden=results.length>0;more.hidden=limit>=results.length;
      status.textContent=`${results.length} result${results.length===1?'':'s'}${q.value?' for “'+q.value+'”':''}${category.value?' · '+categories[category.value]:''}${year.value?' · '+year.value:''}. Showing ${Math.min(limit,results.length)}.`;
    }
    more.addEventListener('click',()=>{const firstNew=list.children.length;limit+=20;draw(true);const target=list.children[firstNew]?.querySelector('a');if(target)target.focus();});
    searchHost.append(form,status,list,empty,more);
    state(form,{q,category,year},['q','category','year'],()=>{
      results=catalogue.records.filter(r=>matches(r,q.value)&&(!category.value||r.category===category.value)&&(!year.value||r.year===year.value));
      if(q.value)results.sort((a,b)=>tokens(q.value).filter(t=>normal(b.title).includes(t)).length-tokens(q.value).filter(t=>normal(a.title).includes(t)).length);
      limit=20;draw();
    });
  }
  // Deep links open their containing native disclosure before moving to the entry.
  function revealHash(){
    let target;try{target=document.getElementById(decodeURIComponent(location.hash.slice(1)));}catch{return;}
    if(!target)return;
    for(let n=target;n;n=n.parentElement)if(n.tagName==='DETAILS')n.open=true;
    if(target.querySelector?.('details.entry-details'))target.querySelector('details.entry-details').open=true;
  }
  revealHash();window.addEventListener('hashchange',revealHash);
  let beforePrint=[];
  window.addEventListener('beforeprint',()=>{beforePrint=[...document.querySelectorAll('main details')].map(n=>[n,n.open]);beforePrint.forEach(([n])=>n.open=true);});
  window.addEventListener('afterprint',()=>{beforePrint.forEach(([n,open])=>n.open=open);beforePrint=[];});
})();
