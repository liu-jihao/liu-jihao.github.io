# Quartic-threefold certificates

This package supplies the four finite certificates used in the affine
quartic-threefold endpoint. The article proves the geometric full-image
bounds, the complete source construction, and the implication from the
certificate predicates to the required integral bound. This package contains
the complete numerical records and an exact validator for those predicates.

The article interfaces are `def:quartic-threefold-certificate`,
`lem:quartic-threefold-certificate-correctness`, and
`subsec:quartic-threefold-data`, in `section.tex`. The resulting theorem is
`thm:quartic-threefold-affine`. The validator does not establish the geometric
applicability assertions or certify the historical fact graph.

## Run the validator

Use Python 3.9 or later. Only the Python standard library is required; all
mathematical arithmetic uses `fractions.Fraction` and integers. From this
directory, on an ordinary machine:

```sh
python3 check.py --data data --report ../fresh-11-quartic-threefold.json
```

Exit status zero means that all four families passed. A malformed record,
missing family/bin/piece/description/candidate/cover, duplicate identifier,
out-of-range index, or failed exact predicate gives nonzero status and a
`FAIL` report identifying the family, bin, record line, and predicate when
available. `manifest.json` documents the format and identities; it is not
an input to mathematical acceptance. The checker reads the four delivered
data files and imports no generator, private path, project state, or fact store.

## Dataset inventory and domains

Each stable identifier starts with `11-quartic-threefold/` and ends with the
family key below. Each file is UTF-8 text with LF line endings.

| Family key / file under `data/` | Deficit coverage | Bins | Pieces | Maximum mass ceiling | Gap to 999/1000 |
|---|---|---:|---:|---|---|
| `two-quadrics.txt` | (1, 3/2] | 25 | 693 | 49900103/50000000 | 49897/50000000 |
| `cubic-plane.txt` | (1, 10/7] | 14 | 436 | 19979391/20000000 | 609/20000000 |
| `simple-chains.txt` | PPC (1, 18/13]; PQP (1, 4/3] | 46 | 1618 | 24966091/25000000 | 8909/25000000 |
| `remaining.txt` | LOW, COMMON, STAR, PATH0/1/2 as below | 40 | 2481 | 19975087/20000000 | 4913/20000000 |

The two-component data also bound the closed lower endpoint 1; the theorem
uses deficits strictly greater than 1. Every successive bin endpoint is
retained. PPC changes its index restrictions at 4/3 and 26/19; PQP changes
them at 6/5 and 22/17. Equality belongs to the preceding regime. LOW covers
(1,8/7]; COMMON and STAR each cover (8/7,6/5]; PATH0 covers (8/7,6/5], PATH1
covers (6/5,9/7], and PATH2 covers (9/7,4/3]. Bins cannot cross those boundaries.

For every bin the initial source level is the smallest multiple of 1/10000
strictly greater than `B_+/(K_*+a_* B_+)`. The terminal level is 18/25. Constants
are `a_*=5000/3997`, `c_*=804700/107919`, `T_*=799999/100000`, and `K_*=T_*-c_*`.
All multiplicities k=3,4,5 and the full corresponding component allocations
are checked. The article pays the other successor states before this finite
comparison and retains both cutting errors in its 1/1000000 reserve.

## Record formats and acceptance

The two-component files have Markdown tables. Aggregate fields are bin
number, lower deficit, upper deficit, initial source level, mass ceiling,
and piece count; the cubic/plane format also places `A_C` before the ceiling.
Each piece has eight fields: bin, left source level, right source level,
left upper bound, right upper bound, left maximizing label, right maximizing
label, and number of active row pairs. A label is `(k, region, (i,j))`, using
the displayed list/tuple notation, or `None` for an empty list.

The validator reconstructs all seven rows in each of the two quadric regions
and all nine rows in each of the four cubic/plane regions. It considers every
row pair and intersects its affine coordinates with every row inequality.
The source partition must equal the complete list of feasible-interval
endpoints together with the initial level plus j/200. Singular pairs,
infeasible pairs, and isolated source levels are counted separately. At both
ends of every open piece it checks every active function, the active count,
the maximizing label, and exact upward rounding to multiples of 10^-8.
It then checks the complete trapezoid sum and aggregate ceiling.

The conductor files contain the following records in each bin, in this order.
Bin numbers start at one; all other indices start at zero. Rationals are
signed integers or signed integer numerators over positive denominators.
Blank separator lines carry no record.

| Record | Meaning |
|---|---|
| `BIN n: FAMILY; (lo,hi]; gamma0=g; mass ceiling=R.` | Deficit bin, family, initial level, and ceiling. |
| `GAMMA i=value;...` | Strictly increasing, contiguous dictionary of all source endpoints. |
| `PIECES` | Entries `left-index,right-index,10^8 V_l,10^8 V_u`, separated by semicolons. |
| `EMPTY` | `description;row-indices;positive-weights;negative-numerator`. |
| `POLYGONS` | `description;cyclic-edge-rows;positive-edge-length-numerators`. |
| `SEGMENTS` | `description;zero-Farkas-rows;weights;line-row;endpoint-pair:endpoint-pair`. |
| `POINTS` (optional header) | `description;zero-Farkas-rows;weights;independent-row-pair`. |
| `FUNCTIONS` | `function-id;representative-description;V-or-E;row-pair`. |
| `CANDIDATES` | `description:V-or-E/row-pair/function-id/left-index/right-index;...`. |
| `UPPER COVERS` | `function-id;C-or-L;left-index,right-index;bound(s);left-witness;right-witness`. |
| `MASS SLACK NUMERATOR n` | Exact positive reduced numerator of `R - mass`. |

The article fixes the row order, all allocation-description indices, and
the formulas for the sorted weighted affine terms. The checker reconstructs
these formulas; rows are not trusted merely because they appear in a file.
Every description has exactly one master-domain record. Farkas contradictions
exclude empty domains. The polygon checks verify a full turn, positive edge
lengths, and containment in each selected halfplane. Zero-Farkas certificates
force equality on their supported rows. For a nontrivial segment, direction
`q-p` and the two strict endpoint signs prove containment. The simple-chain
file has 92 such segments and three separate point records. All ten
segment-format records in the remaining file have coincident endpoints; the
nonzero line tangent and opposite endpoint signs restrict each to that point.

All old vertices and all intersections of a master edge with the tightening
row four are reconstructed with their exact source intervals. Their complete
multiset must match `CANDIDATES`. Grouping identical functions never discards
their intervals. Covers must have precisely the same interval union as the
associated candidates. A `C` cover is constant and lies below every overlapping
piece interpolant. An `L` cover is the exact interpolant chord on one piece.

Every endpoint witness is `e,v,u_0,...,u_(N-1),S`. With D=10^e, it must satisfy
`v=floor(D B)`, `u_i=ceil(D max(xi_i+zeta_i gamma,0))`, and
`S=2 D v-sum(m_i u_i^2)>=0`. Term weights are sorted as (1,1,2) or (1,1,1,1).
The floor, ceilings, slack identity, and resulting exact capacity inequality
are all checked. The files contain 14504 and 15062 individual endpoint
witnesses, respectively; repeated numerical witnesses remain at every distinct
record. The historical generator's dormant `Z` fallback is not part of this
format and is rejected.

The common mass is `4 gamma0^3 + 3 sum((u-l)(V_l+V_u))`. All 5228 pieces are
retained, including zero-capacity pieces. Both the ceiling and the positive
mass-slack numerator are checked where supplied. The finite exceptional
source levels are kept in the article's ordinary O(p^2) error; the uniform
symbolic errors, original positive integral degree, ordinary lifting cutoff,
and complete spent costs are also article obligations.



## Public release

The mathematical data, executable checkers, supporting modules and required TeX views in this directory are byte-identical to the frozen source package. `manifest.json` records their identities and the complete family interfaces. Historical source archives, operator discussions, discovery programs and old execution receipts are retained separately and are not runtime inputs. Run the root dispatcher for a fresh complete report. Keep generated reports outside this section directory, whose inventory is fixed by the root manifest.
