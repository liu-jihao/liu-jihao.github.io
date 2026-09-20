# Section 16: double-fivefold endpoint at order 29/25

This package checks the finite arithmetic used in `thm:double-twentynine`.
The article proves the geometric reductions, the complete-source arguments,
the local algebra, and the analytic reductions for unbounded parameters.
Arithmetic acceptance is conditional on those proofs and the earlier cited
interfaces. It is not a formal verification or a verification of the entire
historical fact graph.

From this directory, using Python 3.9 or later and only its standard library:

```sh
python3 check.py --data data --report ../fresh-16-double-twentynine.json
```

Exit status zero means every required record, exact identity, inequality,
coverage condition, and rendered view passed. Malformed data or any failed
condition gives status one and a JSON failure report. Explicit exceptions,
rather than Python assertions, implement acceptance; `python -O` cannot
disable checks. The checker reads only the four files in `data/`, the eight
files in `views/`, and its local `contract.py` and `render.py` modules. It does
not read the provenance archive, producer, historical receipts, source prose,
private paths, live graph, or network.

## Canonical data and readable views

There is one canonical mathematical record stream per family. The JSON
files are UTF-8, schema `double-twentynine-v1`. A rational is an integer
string or a signed integer numerator followed by `/` and a strictly positive
integer denominator. Decimal and nonfinite numbers, zero or negative
denominators, duplicate JSON keys, unknown fields, and missing/duplicate/extra
record IDs are rejected. Unreduced fractions are allowed to preserve a
source's displayed denominator; equality is exact rational equality.

`16-double-twentynine/states` is `data/tables.json` and `data/flags.json`.
The 28 state records have fields `id,d,e,m,moment,beta,cap,q,terminal_j`.
Their IDs are `d-e-m`. The checker independently reconstructs moments by
the cumulative-binomial formula, starts at minimal multiplicity, and stops
at the first nonpositive intrinsic or first-parent cap. It reconstructs the
least upward thousandth root, every one of the 170 legal edges, all skips,
and all 739 flags. `flags.json` lists each entire flag, including one-vertex
flags. It must equal that complete census. The terminal `J(50/49,C)` entry
of every state is checked against independent optimization. In particular,
`(4,6,8)` is excluded by its zero parent cap, although its intrinsic cap is
positive; the `(3,6,8)` and `(3,6,9)` states remain. The separate
embedding-seven caps are checked without adding those states below this
embedding-six parent.

`tables.json` also contains the exact source parameters, 13 generic group
records with two numerator columns, 12 square-tangent group records, and
four two-column replacement records. Group IDs, their state membership,
labels, column counts, and union coverage are fixed by the article's
formulas. The generic and replacement denominator is 10000; the square
denominator is 100000. A displayed strict ceiling is `floor(D*v)+1`,
including when `D*v` is an integer. No competing state is omitted.

`16-double-twentynine/prices` is `data/queries.json`. Each of the 596 records
has fields `id,state,baseline,cap,value,flag,point,price`. `contract.py`
specifies the complete finite query domain, without evaluating its answers:
21 scenarios at all 28 states; the three affine-tail corrections; the
quartic-threefold cap 4/3; the non-Cartier quartic-fourfold cap 3/2; and
explicit empty, point and segment cases. The scenarios are terminal,
generic lambda/mu, square tangent, b<=3, factorial lambda/mu, distinct-factor
lambda/mu, singular-power lambda/mu, and smooth exponents 3 through 7 at
mu and theta. Impossible positive first-center types have cap zero in
their scenario; this encodes the zero-deficit comparison, not geometric
realizability of that center.

The query value is the maximum of the entire flag functional. The checker
enumerates all feasible nonincreasing tuples from the flag's cumulative
caps and zero. Lemma `lem:double-twentynine-certificate` proves that this
contains every vertex of every closed flag polytope. All ties and endpoints
are retained; cap zero is a point, a one-state positive cap is a segment,
and cap below zero is empty. The empty record must have null `value`,
`flag`, `point`, and `price`. Every other record has an attaining feasible
flag and tuple. Its value is compared with complete vertex enumeration and
an independent Bellman optimization on all cap levels. In the successful
run the complete domain has one empty, five zero-cap, and 590 positive-cap
queries; 1,830 distinct flag/cap grids contain 23,365 feasible grid tuples.
Some tuples need not be vertices; retaining them is harmless because they
are feasible and every vertex is included.

`16-double-twentynine/arithmetic` is `data/arithmetic.json`, containing 45
named asserted rational identities. They are literal transcriptions of the
article's source masses, scalar prices, affine corrections, actual fixed
orders, interval losses and reserves. The checker independently evaluates
their defining expressions and then checks the inequalities that use them.
It also checks all 20 nonzero characters for exponents 3 through 7 and the
eight finite small-exponent interval comparisons. These last finite checks
are diagnostics of the stated analytic endpoint formulas, not a proof for
all residual orders or exponents.

`views/` contains eight deterministic TeX renderings of the canonical data:
`states.tex`, `generic.tex`, `square.tex`, `replacements.tex`, `affine.tex`,
`generic-orders.tex`, `power-orders.tex`, and `singular-reserves.tex`.
The article inputs them from `paper/src/` using
`../supplement/16-double-twentynine/views/<name>.tex`. They remain small
readable tables in the article. The checker requires exact equality with
the rendering of the checked JSON; the views are not independent evidence.
The source identities and other formulas remain readable prose/displays.

## Margins and proof boundary

The terminal maximum J is `384791/294000 < 131/100`; the terminal scalar
gap at n=7 is `3149/142100`, exceeding the 1/100 point margin and 1/1000
total weighted cutting error. Other scalar payments reserve 1/10000 total
weighted cutting error and more than 1/2000 additional margin. The smallest
square-tangent gap is `373/261000`. Raising n increases the available gap
by exactly `4/29` per dimension. The direct quartic threefold keeps its
actual augmented coefficient: the explicit c+y comparison converts to the
earlier affine quartic-threefold theorem's own hypothesis and allowance.

The integral-quadric final normalized mass is
`4343146775241/320000000000000 > 1/100`; the uniform cubic reserve is
`126108747621/320000000000000 > 3/10000`; and the small smooth-power final
mass is `8665174751849/2880000000000000 > 1/500`. The distinct Cartier and
singular quartic/quintic reserves are all separately checked. The common
intrinsic reserve 1/100000 is below every applicable reserve. There is no
floating arithmetic in the delivered checker. It does not estimate the
geometric O(p^4) constants: the article fixes their uniform cutoff before
the degree, prime, order or place and keeps ordinary lifting separate.

The following are substantive article proofs, not algorithmic conclusions:
normality and Cohen--Macaulayness; the actual convergent hypersurface and
splitting; the square tangent's nilpotent and repeated-component argument;
the four exhaustive residual factorizations; the complete-intersection and
Gorenstein restrictions; all corrected integral-quadric, principal and
symbolic local interfaces; old discriminant/moduli compatibility and the
coefficient-one identity; the full-source kernel and threshold construction;
and the all-D/all-character reductions by convexity, the derivative sign
`a1<=4*eta`, and monotonicity. The unbounded smooth range uses the proved
cap `2+1/ell<=15/7` for ell>=7. The finite endpoint loop never substitutes
for these analytic proofs. Positive-part endpoints, integer source orders,
the original degree, all ambient dimension increments, and both cutting
errors stay in the article.



## Public release

The mathematical data, executable checkers, supporting modules and required TeX views in this directory are byte-identical to the frozen source package. `manifest.json` records their identities and the complete family interfaces. Historical source archives, operator discussions, discovery programs and old execution receipts are retained separately and are not runtime inputs. Run the root dispatcher for a fresh complete report. Keep generated reports outside this section directory, whose inventory is fixed by the root manifest.
