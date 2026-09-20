# Complete original-source certificates

This package supplies the finite rational certificates used in the article's
section **Complete original sources**. Its five canonical inputs are the JSON
files in `data/`. The geometric hypotheses, original-degree quantifiers, full
ideal and layer arguments, uniform errors, and reasons the finite comparisons
imply generation remain in the article.

From this directory, with Python 3.10 or later and no third-party packages:

```sh
python3 check.py --data data --report ../fresh-09-original-sources.json
```

Exit zero means that every family passed. Malformed data, an incomplete cover,
an inapplicable polynomial branch, or a failed rational inequality produces
exit one and a `FAIL` report. Detailed exact arithmetic is written to
`report.evidence/`. The checker uses `fractions.Fraction` and integer arithmetic;
it has no floating-point decisions, numerical optimizer, candidate search,
private fact graph, network dependency, or host-specific import. Every check is
an explicit conditional exception, including the adapted Bernstein engine;
Python's `-O` option does not disable the checks. Hashes record identity and
provenance; they are not used to accept mathematical data.

## Canonical inputs and their proof interfaces

| Dataset identifier | Article correctness argument | Complete input and checks |
| --- | --- | --- |
| `09-original-sources/scalar-tail` | `lem:original-sources-smooth-tail` | `scalar-tail.json`: 71 distinct states, 253 affine lines; exhaustive positive-deficit state enumeration, all moments and rational radii, every permitted adjacent successor and dimension skip, and every affine recurrence inequality. The smooth root, threefold multiplicities eight and nine, current sharper row-32 slack `1599/6250`, and the `39/20` gates are also checked. |
| `09-original-sources/volume-two` | `lem:original-sources-larger-loss` | `volume-two.json`: six bounded polynomial identities, one translated half-line identity, 90 coefficients, four endpoint differences, four rounding differences. The displayed polynomials are rebuilt independently and every coefficient is matched exactly and checked positive. |
| `09-original-sources/two-profile` | `lem:original-sources-two-profile-loss` | `two-profile.json`: all Q/M/D grids, 44 ordered bisections, 119 terminal rectangles, 999 bounded and 32 unbounded integration cells, totaling 1,031. Each cell's domain, applicable degree-six formula, positive-part signs, exact integral, strict upward integer rounding, and total are checked. Both unbounded slope ranges are mandatory for each family. |
| `09-original-sources/incidence` | `lem:original-sources-prime-incidence` | `incidence.json`: nine incoming fixed-prime-slope segments and 21 outgoing rectangles. The validator reconstructs all affine roots and all applicable competitors, halves the resulting intervals, integrates each competitor exactly, and reproduces every upward-rounded row total. Coverage is checked on every atomic rectangle of the full arrangement, not just by area. |
| `09-original-sources/retained-divisor` | `lem:original-sources-twentysix-nondivisorial`, using `lem:original-sources-retained-capacity` | `retained-divisor.json`: 14 slope partitions, all 16 closed beta intervals, 1,099 cells (622 Bernstein and 477 shell). The checker reconstructs all branch roots and mesh refinements, all clipped conditional integrals and branch strings, every one of 69,664 Bernstein inequalities, all coefficient and cell roundings, and every row sum. |

The early 1,099-cell retained-divisor certificate and the final 1,031-cell
two-profile certificate have different files, identifiers, parameters and
checkers. They are not interchangeable.

All five top-level objects have `schema: 1` and their full dataset identifier in
`family`. Rational numbers are reduced strings such as `"263/1800"` or `"0"`;
integer witnesses and indices are JSON integers, never booleans or floats.
Unknown or missing object keys, duplicate JSON keys, extra or missing data
files, unknown family or branch names, noncanonical fractions, incorrect list
lengths, duplicate domains, and indices outside their prescribed ranges are
rejected. Every declared field is parsed and used. Data integrity does not
depend on the manifest's expected hashes.

`scalar-tail.json` stores each row as `(id,d,e,m,moment,cap,radius,lines)` with
explicit field names; each line is a pair `(slope,intercept)`. The historical
table's `mu`, `B`, and `r` become `moment`, `cap`, and `radius`; the article
uses `T`, beta, and `R`. All three historical 71-row tables were compared in
full. Multiplicities eight and nine at embedding dimension at most six are
kept distinct from the separate embedding-seven states.

`volume-two.json` stores the D intervals `[0,1/2]`, `[1/2,3/5]`, `[3/5,7/10]`
with degree 12, the H intervals `[0,1/2]`, `[1/2,3/5]`, `[3/5,17/25]` with
degree 13, and the L translate `17/25+y` with degree 8. Coefficients are in
ascending powers. The old scale names `M` and `scale` both become `scale`;
the old root-rounding names `k,C` mean the article's kappa and `C_0`.

`two-profile.json` stores ordered splits `[tag,coordinate,midpoint]`, each leaf's
four bounds and total, and cells `[left,right,polynomial,upper_numerator]`.
The old line numbers are documentary locations and are not certificate data;
the checker reports canonical within-table row numbers. Original slope is
anchored at `9/10` with height `131/1000`. The new anchor is `26/25`, with
Q/M/D heights `37/200`, `41/200`, `189/1000`. The common order endpoint is
`27/25`. Closed bounded grids start at height/anchor and end at 16; separate
`a0` and `a1` tables cover the two half-lines beyond 16, including their
intersection. A strict cell numerator is `floor(10^9 * integral) + 1`, even
when the integral is already an integer multiple of `10^-9`.

`incidence.json` stores `[a0,a1,N]` incoming rows and `[a0,a1,d0,d1,N]`
outgoing rows. Here the old height is `263/2000`, the prime height is `1/5`,
and the valuation-debit upper bound is `3/5`. Incoming rows fix `d=241/100`
and integrate through `26/25`. Outgoing rows cover
`[263/1800,3] x [241/100,5]` and end at
`(d0*(26/25)-1/5)/(d0-1)`. The three legal partial divisions and every
capacity competitor are exactly those displayed in the article. The row
integer is `ceil(10^6 * integral)`. The geometrical proof handles actual
debit zero, smaller initial codimension by zero padding, slopes outside the
bounded rectangles, and eventual zero sources.

`retained-divisor.json` uses case names `prime3` and `no_prime` for the
article's floors `F_3` and zero. All 14 original slope intervals and their
ordered cells are present. `B` cells supply a 16-character choice string and
seven integers `M`; `S` cells supply a rational shell rate. `D` is the cell
numerator with denominator `10^8`. Middle-choice codes are `0=BB`, `1=BC`,
`2=CB`, `3=CC` (centroid/convex in the two halves), `-` for no middle slice,
and `X` for an empty enlarged region. All beta endpoints are supplied:
`b_j=3/5+j/40`, `0<=j<=16`. The enlarged coefficients are
`B=b_(j+1)` and `C=1-b_j`; their sum `41/40` is never renormalized.
Coefficient enclosures use denominator `10^6`. The checker verifies the
enclosures and their exact upward roundings, then the exact upward rounding
of every integrated cell bound. Zero polynomials, empty slices, coincident
affine boundaries, and endpoints are explicitly handled by the branch rules;
only intervals of positive length contribute integrals.

## Bounds actually preserved

| Family | Rounded maximum / target | Strict margin |
| --- | --- | --- |
| Q | `968716327/10^9 < 31/32` | `33673/10^9` |
| M | `968546109/10^9 < 31/32` | `203891/10^9` |
| D | `937247717/10^9 < 31/32` | `31502283/10^9` |
| Incidence incoming | `790166/10^6 < 1-(4/5)^7` | `297/2500000` |
| Incidence outgoing | `999409/10^6 < 1` | `591/10^6` |
| Retained floor `F_3` | `78813928/10^8 < 789/1000` | `10759/12500000` |
| Retained zero floor | `84492074/10^8 < 169/200` | `3963/50000000` |

The common old unbounded Q/M/D numerator is `502032677`; the new ones are
`946259716`, `944079086`, `932247877`. The early contradiction margins are
`803/625000` and `364005211/48828125000`. All four volume-two endpoint
differences and the polynomial positivity tests are in the acceptance report.
Strict geometric error margins are spent only as explained in the article;
the finite checker does not assert new geometric facts.



## Public release

The mathematical data, executable checkers, supporting modules and required TeX views in this directory are byte-identical to the frozen source package. `manifest.json` records their identities and the complete family interfaces. Historical source archives, operator discussions, discovery programs and old execution receipts are retained separately and are not runtime inputs. Run the root dispatcher for a fresh complete report. Keep generated reports outside this section directory, whose inventory is fixed by the root manifest.
