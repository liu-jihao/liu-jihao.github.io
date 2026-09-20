# Section 14: reduced tangents and the three-profile source comparison

This package contains the complete rational witnesses used by `lem:double-geometry-two-tails` and `lem:double-geometry-rational-loss`, and an independent exact checker. The article proves the geometric and asymptotic reductions leading to these finite inequalities. In particular, its original observing degree is divisible by 25, its auxiliary degrees are k_0 divisible by 50 and k_1,k_2 divisible by 25, and its output counts hold for every sufficiently large integer q. Executing the checker does not establish those geometric hypotheses or replace native proof review.

From this directory, with Python 3.10 or later and no third-party packages, run:

```sh
python3 check.py --data data --report ../fresh-14-double-geometry.json
```

The checker reads exactly `data/three-profile.txt` and `data/tail-table.txt`. It imports no generator, optimizer, private project state or network resource. All mathematical operations use integers and `fractions.Fraction`; checks are explicit conditions and remain enabled with `python3 -O`. Exit 0 means all declared data, coverage and arithmetic checks passed. Malformed data or a failed check returns nonzero and writes a failure report when its destination is writable. The JSON report identifies every family, actual counts, exact maxima, strict gaps, unbounded formulas and the arithmetic scope. `receipts/check-final.json` is the saved successful run; it is evidence of that execution, not a substitute for reproducing it.

The canonical dataset identifiers, relative to the supplement root, are:

| Identifier | Canonical file | Content |
|---|---|---|
| `14-double-geometry/tail-table` | `data/tail-table.txt` | All 71 states and 253 affine witnesses |
| `14-double-geometry/three-profile` | `data/three-profile.txt` | All three grids and split trees, 510 parameter boxes, 4,729 integration records, and the unbounded formulas |

There is one authoritative record stream per dataset. The historical source and old article views are retained separately as private archival evidence; they are not additional independent certificates.

## Exact text formats

Both files use UTF-8/ASCII and LF newlines. Every rational is a reduced signed integer or a reduced `numerator/denominator` with positive denominator. Decimals, exponents, leading zeroes and nonreduced fractions are rejected. Integer fields are nonnegative decimal integers without leading zeroes. Every proof-bearing nonblank line must be consumed; unexpected files in `data/` are rejected.

`tail-table.txt` has the literal eight-column Markdown header and separator shown in its first two lines, followed by exactly 71 rows in increasing ID order. Fields are `id | d | E | m | mu | B | r | upper lines (slope,intercept)`. Here E is embedding dimension, B is the numerical deficit cap, and r is an upward radius bound. The last field is a nonempty semicolon-and-space separated list of exact `(slope,intercept)` pairs. The checker enforces the complete state set, each row's expected number of distinct affine witnesses, moments, caps and minimal upward thousandth radii. It checks every affine line against every compatible smaller-dimensional successor, using both cap endpoints and every in-range affine intersection. Adjacent caps and dimension skips remain distinct. Empty successor intervals contribute nothing; zero deficit has zero actual continuation cost, already covered by the direct point alternative. Parallel lines require no division; endpoint ties are retained.

`three-profile.txt` contains exactly these ordered blocks for Q, then M, then D:

1. `FAMILY <family>.hmin=<rational>.`, then one `Grid0=...;grid1=...;grid2=....` line containing three comma-separated rational grids.
2. The literal `Splits: parent tag | coordinate | midpoint` header and all split records. Coordinates are 0,1,2. A root tag `i.j.k` uses one-based consecutive grid intervals. Appending `0` or `1` selects the lower or upper child. Each parent must be an existing leaf and each split lies strictly inside it.
3. The literal ten-column `Rectangles:` header and every leaf record: tag; l0,u0,l1,u1,l2,u2; linked h1,h2; total strict numerator. These are closed parameter boxes, not integration cells.
4. For every leaf, `Rectangle <tag>. l | r | polynomial | upper numerator`, followed by all its four-field integration records. Every list partitions [0,29/25], with positive interval lengths and consecutive equal endpoints.
5. `Unbounded a0>=16: l | r | polynomial | upper numerator`, its five integration records and its total. The last two records give the exact a1>=8 post26 loss and reserve, then the a2>=8 post27 loss and reserve. Their complete spelling is fixed by the stream and checked by the parser.

Blank lines in this stream have no meaning. Duplicate sections, grids, tags, missing intervals, mismatched coordinates, unknown labels, bad rational tokens and unconsumed records are errors. Parsing never skips an unrecognized line as prose. Per-family counts are fixed:

| Family | Root boxes | Splits | Leaf boxes | Bounded integration records | First-slope tail records | Maximum numerator |
|---|---:|---:|---:|---:|---:|---:|
| Q | 150 | 22 | 172 | 1445 | 5 | 9986634132 |
| M | 150 | 38 | 188 | 1827 | 5 | 9989635195 |
| D | 150 | 0 | 150 | 1442 | 5 | 9796511857 |

Thus 510 counts parameter boxes; 4714 counts bounded integration records; 4729 includes the 15 first-slope tail records. The other two unbounded cases per family are checked as rational formulas, not counted as integration records. The common strict-numerator denominator is 10^10, and the target is 9990000000/10^10. The smallest bounded margin is 364805/10^10. The full earlier reserve is 760031169920421/10^16 > 3/40. The checker also verifies all three unbounded loss formulas and their margins above 1/1000, the double-payment error allowance, cubic coefficient/slack bounds and profile-height gaps.

## Mathematical interface

The constants, grids and maps are those of `def:double-geometry-comparison-data`: anchors 9/10,26/25,27/25; terminal order 29/25; h0=131/1000, h*=1609/5000 and epsilon=1/100000. The middle families have hmin=37/200,41/200,189/1000, respectively. On each box evaluate the linked maxima at the lower earlier slopes, then use the upper profile slope before its anchor and the lower slope afterwards. This decreases each supporting line over the entire closed box, including boundary faces, points and segments. No feasible region is discarded: even a box with empty intersection with the linked slope constraints is checked on its full rectangular relaxation.

Each label has the exact polynomial and admissibility range in the article's rational-loss proof. The allowed baseline labels are B1, B61/100, B11/25 and B17/50; profile labels use a prefix 0,1,2 followed by N5lo/N5hi, N3lo/N3hi, C, E1, Qbal, Mspan, Mref, H3, H3lip, Dlo or Dhi, as allowed by that profile's family. Signed capacity polynomials stay intact. The Q crossing brackets are 16020581569/90194313216 and 8010290789/45097156608. The checker verifies their signs and uses no Q branch in the gap.

The checker splits each proposed merged interval at every applicable anchor, baseline jump, positive-part zero and branch boundary. Rational midpoint evaluations determine branches on these exact subintervals, where c(t)=max(0,a-d/t) is monotone; this is not numerical sampling. All reconstructed polynomial coefficient vectors must agree within a merged interval. The article explains boundary continuity, t=0 and the zero-denominator case when a equals a branch value. Integration is exact, and each integer must equal `floor(10^10 * integral) + 1`, including when the integral is zero or integral after scaling. This is strict upward rounding, not ordinary ceiling. Every interval upper bound is summed; no maximizing-box sample replaces full coverage. An accepted polynomial need not minimize the integral among candidates: admissibility and the checked strict total suffice.



## Public release

The mathematical data, executable checkers, supporting modules and required TeX views in this directory are byte-identical to the frozen source package. `manifest.json` records their identities and the complete family interfaces. Historical source archives, operator discussions, discovery programs and old execution receipts are retained separately and are not runtime inputs. Run the root dispatcher for a fresh complete report. Keep generated reports outside this section directory, whose inventory is fixed by the root manifest.
