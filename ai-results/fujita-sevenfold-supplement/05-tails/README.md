This supplement checks the two finite scalar certificates in “Complete descendant costs.” Run from this directory with Python 3.9 or later (tested with CPython 3.10.19); only the standard library is required:

```sh
python3 check.py --data data --report ../fresh-05-tails.json
```

Exit 0 means both complete families passed. Malformed data, missing coverage, a failed inequality, or a wrong claimed support value produces a nonzero exit; certificate failures produce a JSON report with `status: FAIL`. An unwritable output path also fails. The checker neither imports a generator nor reads the provenance archive, private project files, other supplements, or the network. Its only mathematical inputs are the two files in `data/`. Hashes identify the delivered version; acceptance depends on the mathematical checks, not matching a stored hash.

The canonical family identifiers and article interfaces are:

| Family | Canonical file | Article correctness argument | Coverage |
| --- | --- | --- | --- |
| `05-tails/ambient71` | `data/ambient71.json` | `def:tails-sevenfold-certificate`, `lem:tails-sevenfold-certificate`, `eq:tails-line-test` | All 71 states; all 253 affine pairs; 1,377 allowed state edges |
| `05-tails/double-fivefold28` | `data/double-fivefold28.json` | `lem:tails-support-recursion`, `lem:tails-double-fivefold-certificate` | All 28 first successors of the double fivefold; 170 subsequent edges; 739 finite chains |

Both files are UTF-8 JSON objects with exactly `schema_version`, `family`, `parameters`, and `rows`. Version is the integer 1. Parameters are fixed by the checker: degree floor 1, grid 1000, smooth curves only; the ambient family has dimension 7 and intrinsic caps, and the restricted family has parent state `(5,6,2)`, parent moment 1, baseline `50/49`, and a cap applied only to the first successor. The parser rejects extra or missing fields, duplicate JSON keys, floating-point or nonfinite numbers, wrong parameters, booleans in integer fields, and malformed rationals. The data directory must contain exactly the two named ordinary files.

Rows are in lexicographic `(dimension, embedding dimension, multiplicity)` order. All integer fields are JSON integers. Rational fields are strings in reduced form with a positive denominator: `0`, `1`, `283/200`; `2/2`, `1/1`, leading zeros, decimal notation, and signed denominators are rejected. Bounds and radii are reconstructed, so a numerically incorrect but syntactically valid entry is also rejected.

- Ambient rows contain `id` (1 through 71), `state` (three integers), `moment`, `cap`, `radius`, and `lines` (all slope/intercept pairs, in increasing lexicographic order). The fixed per-row pair counts are given in `manifest.json` and `LINE_COUNTS` in the checker. They prevent deletion of even a redundant pair. Duplicate pairs and missing/extra states are rejected.
- Restricted rows contain `id` (1 through 28), `ambient_id`, `state`, `moment`, `intrinsic_cap`, `first_cap`, `radius`, and `support`. The ambient index is checked against the independently enumerated state. The two cap fields are distinct: for a fourfold `first_cap = min(intrinsic_cap,(8-m)/2)`; for lower dimensions it is the intrinsic cap. Subsequent transitions always use the current parent's intrinsic data.

For each embedding codimension `c>0`, the checker sets `ell=floor(d/2)+1` and enumerates `c+1 <= m <= floor((2*ell*S_c(ell-1)-1)/(2*ell-d))`, retaining positive `d-2*T/m`. This bound follows from `T >= ell*(m-S_c(ell-1))` in the article. The smooth state is handled separately; dimension one retains only `(1,1,1)`. Moments are independently reconstructed by filling graded monomial capacities. Integer binary search reconstructs the smallest `j/1000` satisfying `j**d >= m*1000**d`, so both the direction and minimality of every rounded radius are checked without floating point.

For ambient rows, every dimension drop and every allowed embedding is considered. Adjacent caps subtract the child's multiplicity divided by its parent's multiplicity; skips have no adjacent penalty. The checker minimizes all successor lines before maximizing on the closed deficit interval, using both endpoints and every in-range intersection of unequal-slope lines. Parallel lines require no division, repeated intersections are deduplicated, a zero-length interval has one endpoint, and a nonexistent positive successor leaves the point alternative. It checks every instance of the article's affine inequality: 6,889 line/successor comparisons and 34,839 endpoint/intersection evaluations. The minimum affine slack is zero; strictness is required in the separate payment inequalities, not in every certificate line.

For the restricted family, the same independent enumeration produces 30 intrinsic states; the first-successor cap excludes exactly `(4,6,8)` and `(4,6,9)`, leaving counts `1+5+13+9=28`. Recursion on decreasing dimension evaluates the complete support function at `a=50/49` and each supplied `first_cap`. A second evaluator expands every finite chain and maximizes the resulting linear form over the monotone deficit polytope. Its finite coordinate grid consists of zero and all cumulative caps; the article proves that it contains every vertex. Both methods must equal every supplied support entry, not just the largest. Zero intervals are checked for all 28 states, and the subunit-baseline curve value at `a=20/21,h=1` must be `1/21`.

The checker also verifies the two embedding-six repairs, all printed scalar census maxima at orders `26/25`, `29/25`, `34/25`, every divisorial fixed-mass comparison, and all six dimensional maxima at order `13/5`. Important margins in the fresh checker report are `27/2600`, `47/5200`, `49/4875` for the first gates, the common `37/11600` and `91/17000` for the later exclusions, and `53/2600` below `8-1/100` after charging `1/1000` for both weighted cutting errors. The restricted maximum is `384791/294000`, with gaps `11/58800` to `1309/1000` and `349/294000` to `131/100`.

These finite checks certify the numerical premises of the cited local lemmas. The article supplies their induction and geometric applicability: actual Cohen--Macaulay centers, degree floors, persistent parents, named old-prime caps, complete original sources, and both cutting errors. Arbitrary retained-root data and those geometric assertions are not numerical inputs of these two scalar families. This is neither a new formal proof nor a certification of the historical fact graph. 



## Public release

The mathematical data, executable checkers, supporting modules and required TeX views in this directory are byte-identical to the frozen source package. `manifest.json` records their identities and the complete family interfaces. Historical source archives, operator discussions, discovery programs and old execution receipts are retained separately and are not runtime inputs. Run the root dispatcher for a fresh complete report. Keep generated reports outside this section directory, whose inventory is fixed by the root manifest.
