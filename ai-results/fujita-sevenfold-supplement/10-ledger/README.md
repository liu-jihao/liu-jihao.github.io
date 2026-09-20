# Complete original-source ledger certificates

This directory supplies the exact certificates used in the section **The
exhaustive original-source ledger**. The article proves their geometric
applicability and why the checking conditions imply its scalar and integral
bounds. This is an arithmetic certificate checker, not a formalization of the
geometric proof or an audit of every historical source fact.

With Python 3.9 or later and no third-party packages, run from this directory:

```sh
python3 check.py --data data --report ../fresh-10-ledger.json
```

Exit status zero means every family passed. Malformed data, missing coverage,
unknown records or a failed inequality cause a nonzero exit and a FAIL report.
The report contains all scalar recurrence checks and all exact cell integrals,
as well as domain counts, totals, tail fractions and strict margins. The
checker reads exactly `data/scalar.json` and `data/profiles.json`; it does not
read the manifest, provenance, reproduction code, old reports, a fact graph or
the network. Hashes establish file identity separately from arithmetic validity.

## Families and article interface

| Stable identifier | Complete domain and records | Article interface |
|---|---|---|
| `10-ledger/scalar` | 71 states, all 253 affine witnesses, every allowed lower-dimensional successor | `lem:ledger-complete-scalar`, `eq:ledger-line-test` |
| `10-ledger/P2` | `[131/900,16] x [8/45,5]`; 20 initial rectangles, 24 splits, 44 leaves, 340 order cells | `lem:ledger-profile-exclusion`, `lem:ledger-certificate-checking` |
| `10-ledger/T` | `[131/900,16] x [103/540,4]`; 20 initial rectangles, 62 splits, 82 leaves, 717 order cells | Same |
| `10-ledger/old-tail` | `a0 >= 16`, every `a1`; five cells partition `[0,13/5]` | Old-slope part of `lem:ledger-profile-exclusion` |
| `10-ledger/P2-new-tail` | `a1 >= 5`, every `a0`; all orders `t >= 27/25` | New-slope part of that lemma |
| `10-ledger/T-new-tail` | `a1 >= 4`, every `a0`; all orders `t >= 27/25` | New-slope part of that lemma |
| `10-ledger/source-margins` | Three source-specific height comparisons with error `1/1000000` | `eq:ledger-finite-profile-comparison` and the actual-price argument |

The original point orders are `r=9/10`, `u=27/25`, and `Q=13/5`.
The old height is `131/1000`. P2 is the singular divisor at the actual
`27/25` source, with height `24/125`. T is the full reduced Cohen--Macaulay
embedding-seven fivefold at that same source, with height `103/500`.
The historical paid-types producer also investigated a T family at height
`24/125`; that is **not** the T certificate used here. The CM producer supplies
the T family in this delivery. Earlier source histories are unrestricted;
the article retains both actual computing places, full comparison ideals,
full ordinary tangents, all original degrees, and the uniform cutoffs and
remainders required to pass from layer dimensions to volume.

## Canonical formats

Both files are UTF-8 JSON with no duplicate object keys. Rationals are
canonical reduced strings, such as `"0"`, `"5"` and `"131/900"`, with positive
denominators. Integer fields use JSON integers, never booleans or floats.
Every object has exactly its specified fields. Extra input files, fields,
duplicate identifiers, conflicting records and out-of-range values fail.

`scalar.json` has schema `10-ledger/scalar/v1` and an ordered `rows` array.
Each row has integer fields `id,d,e,m,moment`, rational `cap,radius`, and all
`lines`, each a two-element array `[slope,intercept]`. IDs are exactly 1--71
in lexicographic `(d,e,m)` order. Expected witness counts are checked row by
row. The article's scalar table is a complete derived view of this file;
conversion checks equality against all of its entries, including wrapped
affine lists. The data, not that archival view, are the checker's input.

`profiles.json` has schema `10-ledger/profiles/v1`, the fixed `constants`,
two ordered `families` (P2 then T), and `old_tail`. Each family has exactly
`id,height,grid0,grid1,splits,rectangles,new_tail`. A split is
`{tag,coordinate,midpoint}`, with coordinate 0 or 1. The initial tag `i.j`
uses one-based consecutive grid intervals. Append `0` for a lower half and
`1` for an upper half. A rectangle has `tag,bounds,upper,cells`, with bounds
`[l0,u0,l1,u1]`; each cell has `l,r,name,upper`. The rectangle tags are ordered
lexicographically. All cell numerators and totals have denominator `10^9`.
Every cell numerator is **floor of the scaled exact integral plus one**,
including when the integral is zero or already an exact multiple.

The `old_tail` object has `slope_min,upper,cells`. Each `new_tail` object has
`slope_min,denominator,loss,target,volume_floor`, all rational strings. The
checker binds every constant to the source-specific article formulas; these
fields cannot move the expected domain or weaken its target.

## What is checked, and why it suffices

For scalars, the checker fills the Hilbert capacities to compute the minimum
moment, enumerates the entire positive-cap state set, and verifies the least
upward thousandth radius. For each affine witness it considers **every**
embedding-compatible lower-dimensional state. Adjacent successors include
the multiplicity drop; skipped dimensions do not inherit that restriction.
On a positive cap it maximizes the minimum of all successor lines minus the
specified slope at both endpoints and all pairwise intersections. Parallel
lines create no breakpoint. A nonpositive cap admits no positive-deficit
successor; the direct point successor is included separately. The complete
recurrence in `eq:ledger-line-test` proves the bound by induction on dimension.
The state enumeration terminates because the mean filled degree is
nondecreasing in multiplicity and tends to infinity.

For profiles, exact binary subdivision reconstructs all terminal rectangles
from the entire root grids. Every split must consume an active leaf, create
both exact halves and match the reported terminal bounds. This proves
coverage of the closed roots; shared edges and vertices are retained.
There are no empty leaves or zero-length order intervals in these datasets.
Each leaf's ordered cells must start at 0, end at Q and meet without gaps
or overlaps. All cells, including zero polynomials, are checked and summed.

On a rectangle the endpoint slope is the upper slope below its own anchor
and the lower slope above it. It minimizes the supporting height; all
capacities are nonincreasing. Every polynomial name denotes an individual
upper bound for the **minimum of all competing capacities**, so it need not
be the optimizing branch. Before integrating, the checker verifies its
entire legal branch using monotonicity of `a-(av-h)/t`, and checks the signs
of every affine positive-part base at both cell endpoints. Interior anchor
changes of the selected slope and positive-part crossings are rejected.
Order endpoints do not change the integral; parameter edges are covered
pointwise by the same bounds.

The allowed names are `B1`, `B61/100`, `B11/25`, `B17/50`, `Bhigh`, `Bzero`,
`0:N5lo`, `0:N5hi`, `0:C`, `0:E1`, `1:E2`, and `1:Tlow` (T only).
Their formulas, source inequalities and branches are all printed in the
article. A constant B rate is checked against every baseline interval it
meets; `17/50` also dominates the decreasing high baseline. Each legal
expression is a polynomial of degree at most six. The checker uses exact
rational seven-node interpolation and checks its moments 0--6, independently
of the producer's coefficient-expansion integration. It verifies each strict
cell inequality, exact floor-plus-one numerator, sum and worst rectangle.
There is no floating-point error; all upward rounding is charged in totals.
The article's `lem:ledger-certificate-checking` proves the implication for
every parameter, including edges and vertices.

The maximum upper numerators are `999872274` for P2 and `999778839` for T.
They lie below the target `9999/10000` by `27726/10^9` and `121161/10^9`.
The resulting volume margins exceed `127726/10^9` and `221161/10^9`.
The old tail has upper numerator `502032677`, leaving more than
`497967323/10^9`. For the new tails the checker verifies the E2 branch,
positive denominator, exact zero of the affine positive part and integral
identity `(u-2h)^7/(2a1-1)`. At slope thresholds 5 and 4 this gives respectively
`4191719943087/476837158203125` and
`3622557586593623/427246093750000000`, each strictly below `1/100`.
The denominator increases thereafter, so the whole ray is covered. With the
article's available volume `>1/32`, both tails leave more than `17/800`.
Exact margins of the fractions themselves are included in the report.



## Public release

The mathematical data, executable checkers, supporting modules and required TeX views in this directory are byte-identical to the frozen source package. `manifest.json` records their identities and the complete family interfaces. Historical source archives, operator discussions, discovery programs and old execution receipts are retained separately and are not runtime inputs. Run the root dispatcher for a fresh complete report. Keep generated reports outside this section directory, whose inventory is fixed by the root manifest.
