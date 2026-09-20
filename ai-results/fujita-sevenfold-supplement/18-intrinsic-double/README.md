# Intrinsic double-fivefold certificates

This directory accompanies the article section `sec:18-intrinsic-double`,
“Intrinsic double-fivefold sources.” It contains the exact arithmetic used in
`lem:intrinsic-double-prices`, `thm:intrinsic-double-sources`,
`thm:factorial-double-endpoint`, `lem:intrinsic-double-ambient-certificate`,
and `prop:initial-fourfold-history`. It also preserves the distinct rank-two
companion calculation at deficit at most three and terminal order `19/20`.

## Run the independent checker

From this directory, with Python 3.9 or later and its standard library:

```sh
python3 check.py --data data --report ../fresh-18-intrinsic-double.json
```

No installation, network, live fact graph, absolute project path, optimizer,
candidate producer, or archived file is needed. The checker reads the seven
canonical files in `data/`, reconstructs their asserted quantities using exact
integer/rational arithmetic, and writes a JSON report. Exit status zero means
that every declared family passed. A malformed file, missing or duplicate key,
extra file or record, incomplete cover, changed witness, or failed inequality
causes a nonzero exit. Checks remain active under Python's `-O` option.

The arithmetic report is scoped to these formulas and data. The article proves
why the capacities, scalar recurrences, and full local intervals apply to the
actual complete original section spaces. It retains the entire incoming
boundary, source-specific thresholds, all local modules and components, degree
rounding, and the uniform errors. Arithmetic acceptance is not a formal proof
of the geometric premises or an acceptance of the whole historical fact graph.

## Canonical data and notation

Each file is ASCII with an exact tab-separated header and LF-terminated rows.
Only printable ASCII bytes, TAB field delimiters and LF record delimiters are
allowed. The final record must also end in LF. Carriage returns (including
CRLF), vertical tabs, form feeds, the other ASCII control separators and
non-ASCII newline forms are rejected before parsing. Rows are split only at
LF; a different separator is never normalized to LF.
Integers have no leading zeros; rational numbers are reduced fractions with a
positive denominator, or integers when the denominator is one. Empty lines,
comments, unconsumed columns, and conflicting duplicate records are rejected.
The manifest gives the expected schema and hashes; mathematical checking does
not rely on those hashes.

| Family ID | Files | Records and interpretation |
| --- | --- | --- |
| `18-intrinsic-double/scalar` | `states.tsv` | All 28 positive-cap states `(d,e,m)`; moment, first-child cap, rational upper radius, `J(50/49,C)`, and both main prices. |
| `18-intrinsic-double/transitions` | `transitions.tsv` | All eight consecutive intrinsic intervals; endpoints, full shift, loss and cumulative remaining mass. |
| `18-intrinsic-double/rank-two` | `rank-two-prices.tsv`, `rank-two-powers.tsv` | Four quantities for each of the same 28 states, plus every integer exponent `M=3,...,9`. |
| `18-intrinsic-double/ambient` | `ambient-prefix.tsv`, `ambient-boxes.tsv` | All 19 old-profile intervals and 237 closed rectangles, with strict integral numerators at denominator `10^9`. Both unbounded ranges are checked by formulas. |
| `18-intrinsic-double/identities` | `identities.tsv` | All 40 named rational endpoint, shift, mass, price and integral witnesses. Each name and its independently calculated value is checked. |

The main constants are `r0=25/27`, `lambda=1731/2000`, `theta=21/20`,
`B*=18/5`, and `eta=999/1000`. In `states.tsv`, the cubic `(4,6,3)`
`P_theta` is the **whole affine price**, not a bound on the complete scalar
tail. The corresponding article branch applies the cubic affine theorem to
the entire child boundary. The direct `(3,6,4)` `P_lambda` entry uses cap one;
the article treats larger deficits with its separate whole affine price.

The historical column name `theta_b3` in `rank-two-prices.tsv` means the
companion terminal order `rho=19/20`, with deficit cap three. Its cubic row is
the raw scalar recurrence, and the paid cubic affine price is a separate named
identity. `lambda_b18_5` uses deficit cap `18/5`; `J_eta` is the scalar tail at
`1/eta`, before adding the old coefficient. These are deliberately separate
columns and parameter regimes. `M` in the power file is an analytic exponent
(written as `mu` in the article), not the global multiplicity symbol.

The two small article tables are rendered views of `states.tsv`. The unabridged
ambient records occur once as the canonical TSV payload. Original expanded
representations are retained separately in private archival material.

## Complete domains and checking rule

The scalar checker independently enumerates dimensions one through four,
embedding dimensions at most six, and every positive-cap multiplicity, with
the curve and codimension-zero cases treated separately. If `L=floor(d/2)+1`
and `H=binomial(c+L-1,L-1)`, the monomial moment obeys
`T(c,m)>=L(m-H)`. A positive cap therefore forces
`m<L*H/(L-d/2)`. This proves the finite multiplicity cutoff. The recurrence
includes every dimension skip and the adjacent-dimension moment bound. It
retains the embedding-six threefold multiplicities eight and nine, including
their different caps. Every radius is checked by its integer power, and all
six numerical columns of all 28 delivered rows are compared with the
independent recurrence, including the `J(50/49,C)` column.

The intrinsic mesh is exactly
`1731/2000,22/25,9/10,23/25,47/50,24/25,49/50,1,21/20`.
The checker reconstructs all shifts and losses and subtracts each loss in
order, requiring every partial reserve to exceed `1/1000`. It also checks the
degree-two transition and degree-at-least-three kernel. A successful later
source is tested at its own threshold from the original boundary; failed
trials are alternatives and are not added to its price.

For the ambient family, the bounded cover is
`[263/1800,3] x [5399983/18360000,20]`. The checker requires containment,
positive rectangle side lengths, disjoint interiors, exact total area, and
complete vertical fibers on every endpoint and every open strip. In the
delivered data these are 29 fibers. The 19 prefix intervals cover the whole
first-coordinate interval without gaps or repeated interiors. Feasible
slopes further satisfy `a1>=eta(a0)/(27/25)`, with the linked height from the
article. Lower profiles use the upper slope before an anchor and the lower
slope after it. An enlarged feasible interval of zero length is retained;
only a strictly empty interval can have no feasible point. Finite closed
coverage includes boundary segments and vertices.

The exact integration routine reconstructs the full old-capacity polynomials,
the new `E1+E2-E3` and `E1` polynomials at their stated validity ranges, and
every baseline strip. It splits at the anchors, every interior multiple of
`1/100`, and all roots for ratios `0,1/7,1/6,4/21,1/3,1/2,1`.
On each cell it integrates every valid polynomial exactly and takes the
minimum of those integrals, an upper bound for the integral of the pointwise
minimum. Coincident roots occur once; constant ratios, zero slopes, the zero
line, zero-length intervals and weak endpoint branches are retained without
division by zero. The delivered prefix and rectangle records reconstruct
31,627 integration cells. This is a reconstructed-cell count, not a number
of additional stored witnesses.

For each row, the supplied numerator must equal `floor(10^9*U)+1`, including
the case when `10^9*U` is integral. The maximum is `899892923`, leaving margin
`107077/10^9` below `9/10`. The two unbounded ranges `a0>=3` and
`a0<=3,a1>=20` give strict numerators `627860905` and `898876676`, respectively.
The latter uses the complete incoming bound `891/1000` and the exact outgoing
integral `(u1-h1)^7/19`, with its positivity and branch hypotheses checked.

The terminal identities also check the small-slope incoming bound `<2/3`,
outgoing bound `<7/50`, large-slope loss `<1/100`, divisorial bound `<3/50`,
and linked new height `>341/1000`. Their geometric use is conditional on the
actual nongeneration and own-`27/25` large-deficit history in
`prop:initial-fourfold-history`; no unconditional volume `>1/10` is asserted.
The article earns that volume from the two actual original observations.

## Analytic scope of the rank-two companion and factorial endpoint

The factorial endpoint uses deficit cap `18/5` and terminal order `999/1000`.
An unpaid first outcome is an actual Cartier quartic or quintic fourfold prime.
Its full fixed order `N` gives both principal shifts at least `2N`. The exact
minimum normalized shift is `2426571/4206250`; the other candidate is
`41786199/70212500`. The principal theorem permits zero-divisor initial forms.
Its two positive weights sum to two for every fixed analytic order `D>=3`,
so the displayed full interval bound applies to all such orders. The exact
remaining mass is
`188420987150292742295489880001/32053856100156250000000000000000>1/250`.

The additional rank-two data preserve Sections 8--12 of the full endpoint
source retained separately in private archival material. They use cap three and terminal order `rho=19/20`.
The finite exponent list does not establish an assertion about all germs
without the following local reduction, supplied by `thm:rank-two-symbolic`:

* Factor the fixed nonzero convergent `psi` in `uv-psi`. An irreducible first
  power is factorial. Multiple distinct factors give `m+2y<=7` for an actual
  non-Cartier successor: hence cubic `y<=2` and quartic `y<=3/2`.
* A singular irreducible power gives `y<=2` for a non-Cartier successor.
  A smooth power `z^M`, `M>=3`, gives `y<=2+1/M`. A Cartier embedding-six
  prime has multiplicity at least four, even if the parent is not factorial.
  Thus every cubic is non-Cartier and has `y<=7/3`.
* At `rho`, every noncubic scalar price is at most `4086217/513000<7.998`;
  the cubic whole affine price is at most `3836449/513000<7.51`.
  At `lambda`, all states except the cubic and quartic fourfold are paid
  by the full recurrence or the direct quartic-threefold affine alternative.
  The increasing cubic affine price at `y=21/10` is
  `1169868403/155790000<7.51`. Hence all `M>=10`, all singular-power cubics,
  and all distinct-factor cubics are paid at this first trial.
* For `3<=M<=9`, the old cyclic coefficient inequality is
  `b+(2-1/M)c_P<=4`, with the same original discriminant and
  `c_P+t*xi=1`. Failed affine payment, its positive denominator, and the
  decreasing ratio at `b=3` give
  `xi>=(M-1)/(2M-1)*R`, where `R=7552629/7814000`.
  All seven delivered rows check this value, `xi>3/8`, and both full shifts
  at least `xi*(1+1/M)>1/2`. On the entire real interval `[3,9]`, the lower
  comparison follows from
  `-7M^2+80M-153=(9-M)(7M-17)>=0`; denominators are positive.
  It gives the uniform shift `R*80/153=839181/1660475>1/2`.
* The distinct-factor non-Cartier quartic is already paid using `y<=3/2`.
  Every remaining quartic is Cartier or lies in a single irreducible power.
  Its fixed coefficient obeys `xi>=2426571/7121600>1/3`.
  A Cartier prime has both shifts at least `2N`. For a single-power prime,
  the full endpoint weights satisfy `alpha+beta=4`, `y<=beta`, and
  `zeta>=beta`. Absence of scalar payment forces
  `y>=6996401/4044285`, so both normalized shifts are at least
  `20989203/35608000>1/2`. These identities are checked explicitly.

After dilution, these lower bounds exceed every nonfixed coefficient, so
`xi=N/p` is an actual integer fixed order divided by the chosen original
degree. In a full character interval, the positive weights
`D/(D-1)` and `(D-2)/(D-1)` sum to two. All components and opposite characters
remain present. Lowering both shifts to `p/2` bounds the entire interval
loss by `2*((rho-1/2)^5-(lambda-1/2)^5)*p^5/120+O_D(p^4)`.
The exact remaining mass is `15402539679/3200000000000>1/250`.
The quotient has the complete `rho` source as its exact kernel at that same
degree. Its own threshold is paid by the terminal prices above.

The fixed analytic germ determines `D` and its uniform error constant before
any prime or computing place is chosen. There is no assertion of a uniform
Cartier index across different germs. The intrinsic cutoff absorbs the
uniform interval error with reserve `1/1000`; the ordinary lifting cutoff is
fixed separately. Both published cutting errors remain reserved in the
scalar margin, and the cubic affine theorem already includes its own errors.
The checker verifies the corrected factorial grouped price `7982791/999000`
and its post-error gap `821/99900>1/400000`.



## Public release

The mathematical data, executable checkers, supporting modules and required TeX views in this directory are byte-identical to the frozen source package. `manifest.json` records their identities and the complete family interfaces. Historical source archives, operator discussions, discovery programs and old execution receipts are retained separately and are not runtime inputs. Run the root dispatcher for a fresh complete report. Keep generated reports outside this section directory, whose inventory is fixed by the root manifest.
