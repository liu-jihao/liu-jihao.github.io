# Exact paired-profile certificate

Family ID: **19-large-deficit/paired-profile-v1**. Format version: **1**.

This package supplies the exact data for the article's
`lem:large-deficit-integral-certificate`. Its conditions and correctness
proof are `def:large-deficit-certificate` and
`lem:large-deficit-certificate-correctness`. Its geometric applications are
`prop:linked-volume-reserve` and `thm:large-deficit-original`.
The article proves capacity admissibility, geometric applicability,
complete-source convergence and uniform finite-degree error estimates.

The single canonical dataset, `data/paired-profile.txt`, retains all 301
boxes, 4,258 finite integral rows, 11 unbounded-first-slope rows, both
unbounded estimates, both coefficient dictionaries and the complete tree.
There is no encoded duplicate. Saved checking reports are outputs, never
inputs to acceptance.

## Run

Python 3.9 or later and its standard library suffice. No network, installed
dependency, private path, Danus state, fact graph, generator or optimizer
is required. From this package directory:

    python3 check.py --data data --report ../fresh-19-large-deficit.json

Both arguments also accept absolute paths. Exit zero means all checks
succeeded. Malformed data, missing/extra files or records, a failed
inequality or incomplete coverage gives a nonzero exit and a JSON failure
report. Explicit exceptions keep every predicate active under Python
optimization; the checker does not use disableable assert statements.

## Exact record format

The stream is UTF-8 with LF newlines and a terminal newline. Empty lines
are separators; every nonempty line has a prescribed grammar and position.
No prose record is silently ignored. Rational fields are signed integers
or `integer/positive-integer`, never decimal or floating-point literals.
The delivered bytes preserve the original spelling of every rational.
Integer row/reserve numerators have denominator `10^12`.

The records occur in this order:

1. Heading, denominator declaration, two crossing cutoffs and their exact
   signed differences, six chord pairs and nine call pairs.
2. Boxes numbered consecutively 0 through 300. Each specifies closed
   `a` and `d` intervals, `eta` and `d_lower`. The historical field
   **d means the article's supporting slope delta**, never dimension.
   The two lower fields mean `eta_-` and `delta_-`.
3. Each box has an incoming row list and integer sum, or the explicit
   `Use incoming 31/32.` record. Every box has an outgoing row list,
   outgoing sum, reserve upper numerator and total. A row is
   `l,u,density,N`. The lists partition all of `[0,27/25]` and
   `[27/25,29/25]`, respectively.
4. Largest total, all eleven `a>=3` rows and their total, and the exact
   `d>=20` rational bound.
5. Coverage heading, root rectangle, node conventions, leaf count, area,
   coarse rounding records and one JSON tree. `[0,left,right]` bisects
   the first coordinate; `[1,left,right]` bisects the second. Both
   closed children contain the common midpoint. A string `"j"` names
   Box j. No empty-leaf marker occurs or is accepted in this dataset.

The article prints all polynomial definitions, tangent points and valid
ranges. `B0`, if present, is the full density `7t^6`, an alias of `B1`;
zero is labelled `Z`. The checker conservatively accepts `Oh` only for
ratio at most `1/2`, which contains all delivered Oh rows. Other ranges
are exactly those in the article. Every coefficient pair is reconstructed
from its rational tangency point before use, including unused entries.

## Why checking proves the numerical conclusion

The fixed constants are `r=9/10`, `v=27/25`, `T=29/25`,
`epsilon=1/1000000`, `h=263/2000` and `H=27/85-epsilon`.
The domain is `a>=h/r` and
`delta>=max(H,h+(v-r)*a-epsilon)/v`. The finite root is
`[263/1800,3] x [5399983/18360000,20]`. The remaining ranges are
`a>=3` and `delta>=20`.

Each box's lower height, effective slope and linked reserve height are
reconstructed from its endpoints. The upper slope is used before each
profile's anchor and the lower slope after it. The checker verifies every
density's ratio interval and affine positive-part signs on the whole row.
Baseline jumps are checked on separate subintervals. A positive part
crossing zero inside a row is rejected; endpoint zeros and identically
zero expressions are treated exactly. Rows starting at zero use the full
jet polynomial, without division by zero.

For `p(t)=sum p_j*t^j` the checker computes
`sum p_j*(u^(j+1)-l^(j+1))/(j+1)` over the rationals and requires
`0<=integral<N/10^12`. An upper numerator need not be the smallest
possible rounding. Every reserve has a weak upper bound using the same
full linked height. Complete partitions, all sums, the seven incoming
caps and every total are checked independently of the printed totals.
No floating-point tolerance enters.

Recursive bisection reconstructs every leaf rectangle. It requires 300
internal nodes, all 301 indices exactly once, and equality of the derived
rectangle with its box. Area is an additional check, not the coverage
proof. Feasibility is computed by exact intersection with
`v*delta>=H` and `v*delta>=h+(v-r)*a-epsilon`. All 301 supplied
intersections have positive area, with closed edges and vertices
retained. The same bounds remain valid on a feasible segment or point;
the clipping routine also checks explicit empty/segment/point examples.
No set is omitted for having zero area. Time-order rows have positive
length and retain common endpoints by continuity.

Both unbounded ranges use the article's analytic reductions. The eleven
old-slope-three rows may cross the unused anchor v. The delta tail uses
the endpoint capacity and a zero reserve. The checker also evaluates the
initial finite-comparison slack, linked-height lower slack and both
small-slope terminal integrals, with their crossing and positivity tests.

The largest finite total is `999894402477/10^12` at Box 215.
Its gap below `9999/10000` is `5597523/10^12`. The resulting source
reserve is at least `105597523/10^12`, strictly greater than `1/10000`.
The full linked height is unchanged between this reserve and the later
profile comparison. The exact sum of the unbounded-a row upper bounds
with `r^7+1/16` is `667049244495/10^12`, strictly below the printed
upper rounding `667049244496/10^12`.



## Public release

The mathematical data, executable checkers, supporting modules and required TeX views in this directory are byte-identical to the frozen source package. `manifest.json` records their identities and the complete family interfaces. Historical source archives, operator discussions, discovery programs and old execution receipts are retained separately and are not runtime inputs. Run the root dispatcher for a fresh complete report. Keep generated reports outside this section directory, whose inventory is fixed by the root manifest.
