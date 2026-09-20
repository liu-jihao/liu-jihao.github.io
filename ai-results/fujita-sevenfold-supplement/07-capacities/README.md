# Section 07 capacity arithmetic

This package checks the finite arithmetic supporting **Full local images and
ambient capacity**. It requires Python 3.10 or later and only the standard
library. From this directory run:

```sh
python3 check.py --data data --report ../fresh-07-capacities.json
```

Exit status is zero only when every required record passes. A malformed input,
missing/duplicate/unknown record or field, invalid rational, incorrect value,
coefficient, domain or sign produces a failure report and nonzero exit status.
Checks use explicit exceptions, so Python's `-O` option cannot disable them.
The checker reads only `data/certificate.json` and `data/bindings.json`.
It does not import the extractor, historical checker, optimizer, live fact
graph, article files, or any private state, and makes no network requests.

## Data and scope

`data/certificate.json` is the canonical arithmetic record stream, schema
`fujita-capacities-v1`. Integers are JSON integers; rationals are strings
`integer` or `integer/positive-integer`. Printed unreduced fractions are
allowed and evaluated exactly. Decimal numbers, NaN, infinity, denominator
zero and JSON booleans used as integers are rejected. Coefficient vectors
are in ascending degree/index order. Each record has one unique `id` and
its `family`, followed by the exact fields implemented in `check.py`.
The checker fixes the complete expected ID set independently of the data.

Family IDs have the prefix `07-capacities/`:

| Family | Records | Tested interface |
| --- | ---: | --- |
| spanning | 4 | All coefficients, including zero coefficients, for `(n,t,s)=(7,1,2),(7,1,3),(7,0,4),(7,0,3)` |
| pairs | 4 | Both coordinate and balanced branches at `(n,c)=(6,1/5),(7,1/6),(7,3/20),(7,7/40)` |
| constants | 13 | The displayed extrema, strip, correction, rounding and full signed capacity values listed in the article's arithmetic interface |
| polynomial | 3 | The full quartic derivative identity, its four derivative values and two endpoint values |
| source-tables | 17 | Eight complete rows at `51/50`, five at `26/25`, four at `29/25`; all source/cap parameters and displayed slack/deficit witnesses; the own-29 full cap costs and all four Hilbert numerators and first moments |
| fourfold-majorant | 22 | Every row 28–49, including all easy and cap-dependent comparisons, the equality row and the closed deficit intervals |
| baseline | 9 | Every divisorial multiplicity 1–6 on the entire unbounded range `u>=1`; all three embedding-seven fivefold branches on `[29/25,34/25]` |
| integrals | 5 | All eight coefficients of `(t-j)^7`, `j=0,...,4`, checked by differentiating and evaluating at the knot; signed polynomial identities on all 20 combinations of four capacities with `[0,1],...,[3,4],[4,infinity)` |
| crossings | 2 | Exact opposite signs at `16020581569/90194313216` and `8010290789/45097156608`; bracket order is also checked |
| envelopes | 12 | Six joins of `B(7,r)`; eight competing expressions for `E(7,5),Q5,T5` at each of three endpoints; all three `K5` candidates and their domain flag at `0,4/25,1/2` |

Total: **91 records in 10 families**. Every list is complete for its declared
finite domain. `check-report.json` records the exact values, margins and
counts. Extra checks include all denominator signs before division, the
exceptional early `(4,7,4)` slack split, both level-strip margins, the
double-fivefold equality endpoint and the factor-two comparison.

The mathematical correctness and scope argument is
`lem:capacities-finite-certificate` in `finite-certificate.tex`. The original
proofs retain all continuum optimizations, zero/coincident weights and their
limits, compact and unbounded parameter ranges, all competing capacities,
full primary/ordinary/tangent images, the own-source geometric hypotheses,
the whole scalar census, both cutting errors and uniform degree/lattice
errors. A finite list of values does **not** check those proofs. The global
scalar census and unrelated historical height/surface computations are not
claimed as checked by this package. Acceptance is not formal verification
or whole-paper acceptance.

The primitives use the zero branch for `t<=j` and the recorded polynomial
for `t>=j`, with matching value zero at the knot. Splitting an interval at
all contained knots gives complete coverage; empty active intervals and
point intervals contribute zero, and endpoints are included. Each cell,
including the unbounded last cell, is checked by polynomial coefficient
equality, not by sampling its endpoints. The union-volume proof in the
article supplies nonnegativity and monotonicity. Scaling
`t` by a positive `eta` gives the article's integral. There is no quadrature
or floating-point error. The continuous and degree-rounding errors in the
article remain `O((ell+1)^5)` per layer and `O(q^6)` after summation.



## Public release

The mathematical data, executable checkers, supporting modules and required TeX views in this directory are byte-identical to the frozen source package. `manifest.json` records their identities and the complete family interfaces. Historical source archives, operator discussions, discovery programs and old execution receipts are retained separately and are not runtime inputs. Run the root dispatcher for a fresh complete report. Keep generated reports outside this section directory, whose inventory is fixed by the root manifest.
