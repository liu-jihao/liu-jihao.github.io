#!/usr/bin/env python3
"""Portable exact validator for the five original-source certificate families.

Usage: python3 check.py --data data --report report.json
Python >=3.10, standard library only. No producer, private state or hash oracle.
Every input field is type/range checked; every proof record is consumed.
"""
import argparse
from collections import Counter
from fractions import Fraction as F
from hashlib import sha256
import json
from pathlib import Path
import re
import sys

SCHEMA_CHECKS = Counter()
FILES = {'scalar-tail.json', 'volume-two.json', 'two-profile.json', 'incidence.json', 'retained-divisor.json'}


def need(ok, message):
    SCHEMA_CHECKS[message] += 1
    if not ok:
        raise ValueError(message)


def keys(obj, expected):
    need(type(obj) is dict and set(obj) == set(expected.split()), 'exact object keys: ' + expected)


def seq(value, length=None):
    need(type(value) is list and (length is None or len(value) == length), 'list length ' + str(length))
    return value


def integer(value, low=0, high=None):
    need(type(value) is int and value >= low and (high is None or value <= high), 'integer range')
    return value


def rational(value, low=None, high=None):
    need(type(value) is str and re.fullmatch(r'-?(0|[1-9][0-9]*)(/[1-9][0-9]*)?', value) is not None,
         'rational syntax')
    q = F(value)
    need(str(q) == value, 'reduced canonical rational')
    need((low is None or q >= low) and (high is None or q <= high), 'rational range')
    return q


def no_duplicate_keys(pairs):
    result = {}
    for key, value in pairs:
        if key in result:
            raise ValueError('duplicate JSON key: ' + key)
        result[key] = value
    return result


def forbidden_number(value):
    raise ValueError('non-integer JSON number: ' + value)


def read_json(path):
    return json.loads(path.read_text(encoding='utf-8'), object_pairs_hook=no_duplicate_keys,
                      parse_float=forbidden_number, parse_constant=forbidden_number)


def family(obj, name, extra):
    keys(obj, 'schema family ' + extra)
    need(type(obj['schema']) is int and obj['schema'] == 1, 'schema version')
    need(obj['family'] == '09-original-sources/' + name, 'stable family identifier')


def unique(values, message):
    need(len(values) == len(set(values)), message)


def load_scalar(obj):
    family(obj, 'scalar-tail', 'states')
    rows = []
    for i, row in enumerate(seq(obj['states'], 71), 1):
        keys(row, 'id d e m moment cap radius lines')
        need(integer(row['id'], 1, 71) == i, 'scalar ids 1 through 71 in order')
        d = integer(row['d'], 1, 6)
        e = integer(row['e'], d, 7)
        m = integer(row['m'], 1)
        pairs = [tuple(rational(v, 0) for v in seq(pair, 2)) for pair in seq(row['lines'])]
        need(len(pairs) > 0, 'nonempty scalar affine list')
        unique(pairs, 'unique scalar affine lines')
        rows.append(dict(id=i, d=d, e=e, m=m, moment=integer(row['moment']),
                         cap=rational(row['cap'], 0, d), radius=rational(row['radius'], 1), lines=pairs))
    unique([(r['d'], r['e'], r['m']) for r in rows], 'unique scalar states')
    need(sum(len(r['lines']) for r in rows) == 253, 'all 253 scalar affine witnesses')
    return rows


def load_volume(obj):
    family(obj, 'volume-two', 'bounded unbounded endpoint_differences rounding_differences')
    expected = [('D', '0', '1/2', 12), ('D', '1/2', '3/5', 12), ('D', '3/5', '7/10', 12),
                ('H', '0', '1/2', 13), ('H', '1/2', '3/5', 13), ('H', '3/5', '17/25', 13)]
    bounded = []
    for row, (kind, lo, hi, degree) in zip(seq(obj['bounded'], 6), expected):
        keys(row, 'kind interval degree scale coefficients')
        interval = [rational(v, 0) for v in seq(row['interval'], 2)]
        need((row['kind'], interval, integer(row['degree'])) == (kind, [F(lo), F(hi)], degree),
             'entire bounded polynomial interval cover')
        bounded.append(dict(kind=kind, interval=interval, degree=degree,
                            scale=rational(row['scale'], 1),
                            coefficients=[integer(v, 1) for v in seq(row['coefficients'], degree+1)]))
    row = obj['unbounded']
    keys(row, 'shift degree scale coefficients')
    need(rational(row['shift']) == F(17,25) and integer(row['degree']) == 8, 'entire polynomial half-line')
    unbounded = dict(shift=F(17,25), degree=8, scale=rational(row['scale'], 1),
                     coefficients=[integer(v, 1) for v in seq(row['coefficients'], 9)])
    keys(obj['endpoint_differences'], 'divisor_a_le_1 double_a_le_1 divisor_large_at_7_10 divisor_large_limit')
    return dict(bounded=bounded, unbounded=unbounded,
                endpoint_differences={k:rational(v, 0) for k,v in obj['endpoint_differences'].items()},
                rounding_differences=[rational(v, 0) for v in seq(obj['rounding_differences'], 4)])


def tag(value):
    need(type(value) is str and re.fullmatch(r'[1-5]\.[1-5][01]*', value) is not None, 'tree tag syntax')
    return value


def qmd_cells(raw):
    cells = []
    for i, cell in enumerate(seq(raw), 1):
        l, r, name, n = seq(cell, 4)
        l, r = rational(l, 0, F(27,25)), rational(r, 0, F(27,25))
        need(l < r, 'positive Q/M/D cell length')
        need(type(name) is str and re.fullmatch(r'B(1|61/100|11/25|17/50)|[01]:(N[35](lo|hi)|C|E1|Qbal|Mspan|Mref|H3|H3lip|Dlo|Dhi)', name) is not None,
             'paper polynomial dictionary')
        cells.append((l, r, name, integer(n, 1, 968749999), i))
    need(len(cells) > 0, 'nonempty Q/M/D cell list')
    unique([(c[0], c[1]) for c in cells], 'unique Q/M/D cell intervals')
    return cells


def load_qmd(obj):
    family(obj, 'two-profile', 'families')
    expected = {'Q': (F(37,200), 31, 56, 477, 9), 'M': (F(41,200), 13, 38, 321, 11),
                'D': (F(189,1000), 0, 25, 201, 12)}
    output = {}
    for raw, kind in zip(seq(obj['families'], 3), ['Q','M','D']):
        keys(raw, 'kind grid0 grid1 splits rectangles tails')
        need(raw['kind'] == kind, 'exact ordered Q/M/D families')
        h1, ns, nr, nc, nt = expected[kind]
        g0 = [rational(v) for v in seq(raw['grid0'], 6)]
        g1 = [rational(v) for v in seq(raw['grid1'], 6)]
        need(g0 == [F(131,900),F(1),F(2),F(4),F(8),F(16)], 'Q/M/D fixed old slope grid')
        need(g1 == [h1/F(26,25),F(1),F(2),F(4),F(8),F(16)], 'Q/M/D fixed new slope grid')
        splits = []
        for i, s in enumerate(seq(raw['splits'], ns), 1):
            name, coordinate, midpoint = seq(s, 3)
            splits.append((tag(name), integer(coordinate, 0, 1), rational(midpoint, 0, 16), i))
        unique([s[0] for s in splits], 'unique internal tree nodes')
        rectangles, cells = {}, {}
        for i, row in enumerate(seq(raw['rectangles'], nr), 1):
            keys(row, 'tag bounds total cells')
            name = tag(row['tag'])
            need(name not in rectangles, 'unique terminal tree nodes')
            bounds = [rational(v, 0, 16) for v in seq(row['bounds'], 4)]
            need(bounds[0] < bounds[1] and bounds[2] < bounds[3], 'positive slope rectangle')
            rectangles[name] = (bounds, integer(row['total'], 1, 968749999), i)
            cells[name] = qmd_cells(row['cells'])
        tails = {}
        for raw_tail, name in zip(seq(raw['tails'], 2), ['a0','a1']):
            keys(raw_tail, 'profile total cells')
            need(raw_tail['profile'] == name, 'both ordered unbounded slope tails')
            tails[name] = dict(total=integer(raw_tail['total'], 1, 968749999), cells=qmd_cells(raw_tail['cells']))
        need(sum(map(len, cells.values())) == nc, 'exact Q/M/D rectangle cell count')
        need(sum(len(t['cells']) for t in tails.values()) == nt, 'exact Q/M/D unbounded cell count')
        output[kind] = dict(grid0=g0, grid1=g1, splits=splits, rectangles=rectangles, cells=cells, tails=tails)
    return output


def load_incidence(obj):
    family(obj, 'incidence', 'incoming outgoing')
    out = {}
    for table, count, width in [('incoming',9,3), ('outgoing',21,5)]:
        rows = []
        for row in seq(obj[table], count):
            seq(row, width)
            a, b = [rational(v, F(263,1800), 3) for v in row[:2]]
            need(a < b, 'positive incidence slope interval')
            if table == 'outgoing':
                d, e = [rational(v, F(241,100), 5) for v in row[2:4]]
                need(d < e, 'positive incidence prime-slope interval')
                rows.append((a,b,d,e,integer(row[-1], 1, 999999)))
            else:
                rows.append((a,b,integer(row[-1], 1, 999999)))
        unique([r[:-1] for r in rows], 'unique incidence domains')
        out[table] = rows
    return out


def load_retained(obj):
    family(obj, 'retained-divisor', 'beta_intervals rows')
    betas = [[rational(v, F(3,5), 1) for v in seq(b,2)] for b in seq(obj['beta_intervals'],16)]
    need(betas == [[F(3,5)+F(j,40),F(3,5)+F(j+1,40)] for j in range(16)], 'full closed beta interval cover')
    counts = [97,96,97,90,84,79,80,76,75,70,58,62,65,70]
    totals = [78813928,78316263,78478991,78479182,78294000,78461598,77466558,76864869,75668322,65463863,
              83874182,84492074,77012365,65467130]
    rows = []
    for i, (raw, count, total) in enumerate(zip(seq(obj['rows'],14), counts, totals)):
        keys(raw, 'case number a0 a1 cells printed_sum')
        case, number = ('prime3',i+1) if i < 10 else ('no_prime',i-9)
        need(raw['case'] == case and integer(raw['number'],1,10) == number, 'all ordered retained row identifiers')
        a0, a1 = rational(raw['a0'],F(263,1800),3), rational(raw['a1'],F(263,1800),3)
        need(a0 < a1, 'positive retained slope interval')
        need(integer(raw['printed_sum'],1,84499999) == total, 'recorded retained row sum')
        cells = []
        for j, cell in enumerate(seq(raw['cells'],count), 1):
            need(type(cell) is dict and cell.get('kind') in ('B','S'), 'retained cell kind')
            keys(cell, 'l v D kind choices M' if cell['kind'] == 'B' else 'l v D kind rate')
            converted = dict(l=rational(cell['l'],0,F(26,25)), v=rational(cell['v'],0,F(26,25)),
                             D=integer(cell['D'],0,84499999), kind=cell['kind'], line=j)
            need(converted['l'] < converted['v'], 'positive retained cell length')
            if cell['kind'] == 'B':
                need(type(cell['choices']) is str and re.fullmatch(r'[-X0-3]{16}',cell['choices']) is not None,
                     'all sixteen retained branch witnesses')
                converted.update(choices=cell['choices'], M=[integer(v,0) for v in seq(cell['M'],7)])
            else:
                converted['rate'] = rational(cell['rate'],0,1)
            cells.append(converted)
        unique([(c['l'],c['v']) for c in cells], 'unique retained cell domains')
        rows.append(dict(case=case,number=number,a0=a0,a1=a1,cells=cells,printed_sum=total))
    for kind in ['prime3','no_prime']:
        selected = [r for r in rows if r['case'] == kind]
        need(selected[0]['a0'] == F(263,1800) and selected[-1]['a1'] == 3, 'full retained slope domain endpoints')
        need(all(a['a1'] == b['a0'] for a,b in zip(selected,selected[1:])), 'retained slope cover without gaps or overlaps')
    return rows


def serial(value):
    if isinstance(value,F):
        return str(value)
    if isinstance(value,dict):
        return {k:serial(v) for k,v in value.items()}
    if isinstance(value,(list,tuple)):
        return [serial(v) for v in value]
    return value


def write_json(path, value):
    path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(serial(value),indent=2)+'\n')


def table_tex(rows, columns, header):
    return ('% Derived view of canonical JSON, verified by check.py.\n\\[\n\\begin{array}{'+columns+'}\n'+header+
            '\\\\ \\hline\n'+'\\\\\n'.join('&'.join(str(v) for v in row) for row in rows)+'\n\\end{array}\n\\]\n')


def article_views(incidence, retained):
    totals = [[i+1,retained[i]['printed_sum'],retained[i+10]['printed_sum'] if i < 4 else ''] for i in range(10)]
    return {'incidence-incoming.tex':table_tex(incidence['incoming'],'c|c|r','a_0&a_1&N'),
            'incidence-outgoing.tex':table_tex(incidence['outgoing'],'c|c|c|c|r','a_0&a_1&d_0&d_1&N'),
            'retained-totals.tex':table_tex(totals,'c|r|r','\\text{row}&f=F_3&f=0')}


def audit_qmd(families, evidence):
    import two_profile as q
    result = {}
    settings = {'Q':(F(37,200),56,968716327,946259716), 'M':(F(41,200),38,968546109,944079086),
                'D':(F(189,1000),25,937247717,932247877)}
    need((1-q.QMIN)**6 < (1-2*q.QMIN)**5*(1+10*q.QMIN), 'Q crossing lower rational endpoint')
    need((1-q.QMAX)**6 > (1-2*q.QMAX)**5*(1+10*q.QMAX), 'Q crossing upper rational endpoint')
    for kind, f in families.items():
        h, count, maximum, tailmax = settings[kind]
        q.CONTEXT = kind
        q.coverage(f,h,F(26,25),F(16),count)
        records, tails = {}, {}
        for name,(bounds,total,_) in f['rectangles'].items():
            records[name] = q.audit_cells(f['cells'][name],bounds,kind,h,F(26,25),F(27,25),total)
        need(max(v['upper_numerator'] for v in records.values()) == maximum, 'Q/M/D exact family maximum')
        for name,t in f['tails'].items():
            tails[name] = q.audit_cells(t['cells'],None,kind,h,F(26,25),F(27,25),t['total'],name)
        need([tails[s]['upper_numerator'] for s in ['a0','a1']] == [502032677,tailmax], 'both unbounded exact totals')
        for value in list(records.values())+list(tails.values()):
            need(F(value['upper_numerator'],10**9) < F(31,32), 'strict Q/M/D target for every bounded and unbounded domain')
        write_json(evidence/(kind+'.json'),dict(rectangles=records,tails=tails))
        result[kind] = dict(splits=len(f['splits']),rectangles=len(records),rectangle_cells=sum(len(v['cells']) for v in records.values()),
                            tail_cells=sum(len(v['cells']) for v in tails.values()),maximum_numerator=maximum,
                            target='31/32',margin=str(F(31,32)-F(maximum,10**9)),
                            tail_numerators=[502032677,tailmax],evidence=kind+'.json')
        print('PASS two-profile/'+kind,flush=True)
    need(sum(r['rectangle_cells']+r['tail_cells'] for r in result.values()) == 1031, 'all 1031 final cells')
    return dict(families=result,checks=dict(q.CHECKS),interpolation_weights=serial(q.WEIGHTS))


def audit_incidence(data,evidence):
    import incidence as inc
    inc.coverage(data['incoming'],data['outgoing'],F(263,1800))
    result = {}
    for table, rows in data.items():
        details, all_cells = [], {}
        target = 1-F(4,5)**7 if table == 'incoming' else F(1)
        for i,row in enumerate(rows,1):
            a0,a1 = row[:2]
            d0,d1 = (F(241,100),F(241,100)) if table == 'incoming' else row[2:4]
            end = F(26,25) if table == 'incoming' else (d0*F(26,25)-F(1,5))/(d0-1)
            if table == 'outgoing':
                need(end > F(26,25) and F(1,5)+d0*(end-F(26,25)) == end, 'outgoing exact vanishing endpoint')
            exact,cells = inc.certificate(a0,a1,d0,d1,end,F(263,2000),F(3,5))
            need(inc.ceiling_million(exact) == row[-1], 'incidence exact upward rounding')
            need(F(row[-1],10**6) < target, 'incidence strict row target')
            details.append(dict(row=i,domain=serial(row[:-1]),endpoint=str(end),exact_integral=str(exact),
                                upper_millionths=row[-1],half_cells=len(cells)))
            all_cells[str(i)] = cells
        maximum = max(r[-1] for r in rows)
        need(maximum == (790166 if table == 'incoming' else 999409), 'incidence exact maximum')
        write_json(evidence/('incidence-'+table+'.json'),all_cells)
        result[table] = dict(rows=details,count=len(rows),maximum=maximum,target=str(target),
                             margin=str(target-F(maximum,10**6)),evidence='incidence-'+table+'.json')
    return dict(tables=result,checks=dict(inc.checks),interpolation_weights=serial(inc.WEIGHTS))


def audit_retained(rows,evidence):
    import retained as rt
    totals = Counter()
    summaries = []
    with (evidence/'retained-cells.jsonl').open('w') as out:
        for row in rows:
            rt.CONTEXT = f"{row['case']} row {row['number']}"
            rebuilt,branch = rt.partition(row)
            cells = row['cells']
            need([cells[0]['l']]+[c['v'] for c in cells] == rebuilt, 'exact full retained branch/refinement partition')
            need(rebuilt[0] == 0 and rebuilt[-1] == F(26,25), 'retained full order endpoints')
            need(all(c['l'] == rebuilt[i] for i,c in enumerate(cells)), 'retained adjacent cells no gaps or overlaps')
            out.write(json.dumps(serial(dict(record='row',case=row['case'],number=row['number'],a0=row['a0'],a1=row['a1'],
                                            branch_knots=branch,refined_knots=rebuilt)))+'\n')
            counts = Counter()
            for cell in cells:
                l,v = cell['l'],cell['v']
                rt.CONTEXT = f"{row['case']} row {row['number']} cell {cell['line']} [{l},{v}]"
                need(0 < v-l <= F(1,40), 'retained mesh width')
                need(not any(l < point < v for point in [F(9,10),F(99,100)]), 'retained shell breakpoint coverage')
                h,f = rt.profiles(row,l,v)
                rate = F(1) if v <= F(9,10) else F(61,100) if v <= F(99,100) else F(11,25)
                record = dict(record='cell',case=row['case'],number=row['number'],**cell)
                if cell['kind'] == 'S':
                    need(cell['rate'] == rate, 'retained exact shell rate')
                    exact = rate*(v**7-l**7)
                else:
                    coefficients, codes, bins = [], '', []
                    for j,(lowbeta,B) in enumerate(rt.BETAS):
                        poly,code,metadata = rt.interval_polynomial(h,f,B,1-lowbeta,l,v)
                        codes += code
                        bern = rt.bernstein(poly,l,v)
                        need(rt.at(poly,l) >= 0 and rt.at(poly,v) >= 0, 'retained polynomial endpoint nonnegativity')
                        gaps = [M*b.denominator-10**6*b.numerator for M,b in zip(cell['M'],bern)]
                        need(all(g >= 0 for g in gaps), 'all seven Bernstein inequalities in every beta interval')
                        coefficients.append(bern)
                        bins.append(dict(beta_index=j,polynomial=poly,bernstein=bern,positive_cross_products=gaps,**metadata))
                    need(codes == cell['choices'], 'all sixteen retained branch witnesses recomputed')
                    maxima = [max(c[i] for c in coefficients) for i in range(7)]
                    need([rt.ceilq(m*10**6) for m in maxima] == cell['M'], 'exact upward Bernstein coefficient rounding')
                    exact = (v-l)*F(sum(cell['M']),10**6)
                    record['beta_records'] = bins
                    totals['coefficient_inequalities'] += 112
                need(exact <= F(cell['D'],10**8), 'retained integrated enclosure')
                need(rt.ceilq(exact*10**8) == cell['D'], 'retained exact upward cell rounding')
                record.update(exact_enclosure=exact,D_cross_product=cell['D']*exact.denominator-10**8*exact.numerator)
                out.write(json.dumps(serial(record),separators=(',',':'))+'\n')
                counts[cell['kind']] += 1
            total = sum(c['D'] for c in cells)
            need(total == row['printed_sum'], 'retained sum consumes every cell')
            target = F(789,1000) if row['case'] == 'prime3' else F(169,200)
            need(F(total,10**8) < target, 'retained strict row target')
            summaries.append(dict(case=row['case'],number=row['number'],slope_interval=serial([row['a0'],row['a1']]),
                                  cells=len(cells),B=counts['B'],S=counts['S'],total=total,target=str(target),
                                  margin=str(target-F(total,10**8))))
            totals.update(counts)
            print(f"PASS retained/{row['case']}/{row['number']}: {len(cells)} cells",flush=True)
    need(totals['B'] == 622 and totals['S'] == 477 and totals['coefficient_inequalities'] == 69664,
         'all 1099 retained cells and 69664 Bernstein inequalities consumed')
    return dict(rows=summaries,cells=1099,beta_intervals=16,counts=dict(totals),
                polynomial_checks=rt.CHECKS,evidence='retained-cells.jsonl')


def audit_scalar_extras(rows):
    byid = {row['id']:row for row in rows}
    need((F(5,4),F(187,250)) in byid[66]['lines'], 'smooth-root line 5/4,187/250')
    for m,cap,radius in [(8,F(1,4),F(2)),(9,F(1,9),F(9,4))]:
        selected = [r for r in rows if r['d'] == 3 and r['e'] <= 6 and r['m'] == m]
        need(len(selected) > 0 and all(r['cap'] <= cap for r in selected), 'm8/m9 deficit caps')
        need(radius**3 >= m, 'm8/m9 rational radius enclosure')
        for row in selected:
            for nxt in rows:
                if nxt['d'] == 2 and nxt['e'] <= row['e'] and row['cap']-F(nxt['m'],m) > 0:
                    need(m == 8 and nxt['m'] < 2, 'm8/m9 complete adjacent-successor restriction')


def scalar_gaps():
    gaps = {'initial_height':F(5,38)-F(263,2000)-F(1,40000),
            'old_large_slope':F(5,7)-F(9,10)**7*F(10,7),
            'prime_large_slope':F(3,40)-F(21,25)**7/4,
            'prime3_contradiction':F(211,1000)-F(4,5)**7,
            'd_ge3_contradiction':1-F(169,200)-F(21,25)**7/2}
    need(all(v > 0 for v in gaps.values()), 'all early unbounded/contradiction margins positive')
    need(gaps['prime3_contradiction'] == F(803,625000), 'prime3 exact contradiction margin')
    need(gaps['d_ge3_contradiction'] == F(364005211,48828125000), 'd>=3 exact contradiction margin')
    return serial(gaps)


def run(data_dir,report_path):
    need(data_dir.is_dir(), 'data directory exists')
    need({p.name for p in data_dir.iterdir()} == FILES and all(p.is_file() for p in data_dir.iterdir()),
         'exact five data files, no extras or omissions')
    loaded = {p:read_json(data_dir/p) for p in sorted(FILES)}
    rows = load_scalar(loaded['scalar-tail.json'])
    volume = load_volume(loaded['volume-two.json'])
    qmd = load_qmd(loaded['two-profile.json'])
    inc = load_incidence(loaded['incidence.json'])
    retained = load_retained(loaded['retained-divisor.json'])
    evidence = report_path.with_name(report_path.stem+'.evidence')
    evidence.mkdir(parents=True,exist_ok=True)
    from volume_scalar import audit
    vol = audit(volume,rows)
    audit_scalar_extras(rows)
    write_json(evidence/'volume-scalar.json',vol)
    result = {'volume-two':dict(identities=7,coefficients=vol['details']['polynomial_coefficients_independently_matched'],
                               endpoint_differences=vol['details']['endpoint_positive_differences']),
              'scalar-tail':dict(states=71,affine_lines=253,current_row32_slack='1599/6250',m8_m9_completion=True)}
    print('PASS volume-two and scalar-tail',flush=True)
    result['two-profile'] = audit_qmd(qmd,evidence)
    result['incidence'] = audit_incidence(inc,evidence)
    print('PASS incidence: all 30 rows',flush=True)
    result['retained-divisor'] = audit_retained(retained,evidence)
    views = Path(__file__).resolve().parent/'views'
    expected = article_views(inc,retained)
    need({p.name for p in views.iterdir()} == set(expected), 'exact derived article view set')
    for name,content in expected.items():
        need((views/name).read_bytes() == content.encode(), 'canonical data/rendered table identity: '+name)
    return dict(status='PASS',arithmetic='exact fractions.Fraction; explicit checks active under python -O',
                families={'09-original-sources/'+k:v for k,v in result.items()},
                early_strict_gaps=scalar_gaps(),schema_checks=sum(SCHEMA_CHECKS.values()),
                data_sha256={p:sha256((data_dir/p).read_bytes()).hexdigest() for p in sorted(FILES)},
                evidence_directory=evidence.name,
                scope='Finite original-source certificates only; geometric applicability and uniform errors are proved in the article.')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--data',required=True,type=Path)
    parser.add_argument('--report',required=True,type=Path)
    args = parser.parse_args()
    try:
        result = run(args.data,args.report)
    except Exception as exc:
        result = dict(status='FAIL',error_type=type(exc).__name__,error=str(exc),schema_checks=sum(SCHEMA_CHECKS.values()))
        write_json(args.report,result)
        print(json.dumps(result),file=sys.stderr,flush=True)
        return 1
    write_json(args.report,result)
    print('PASS: all five certificate families; report '+str(args.report),flush=True)
    return 0


if __name__ == '__main__':
    sys.exit(main())
