# Cubic-fourfold certificate supplement

This package checks the finite arithmetic in the two cubic-fourfold endpoints. It preserves the distinct affine input and actual whole-boundary multiplicity input. Geometric proofs, complete-source constructions, uniformity and cutting errors remain in the article. Arithmetic success does not replace native mathematical verification.

From this directory, with Python 3.10 or later:

    python3 check.py --data data --report ../fresh-12-cubic-fourfold.json

Only the standard library is required. The checker reads exactly the eight files in data/. It imports no producer, archive, Danus module, live fact graph, private path, optimizer or network resource. It exits 0 only after all checks succeed; failures exit 1 with a diagnostic report. Explicit exceptions keep validation active under Python -O.

## Canonical data and article

The JSON files in data/ are canonical records. The three TeX files are derived views; every checker run compares their text with the rendering of those records after newline normalization. The article inputs them from paper/src/ using ../supplement/12-cubic-fourfold/data/<name>.tex. The small readable tables remain printed; they are not independent duplicate certificates.

| Stable family | Records and article interface |
| --- | --- |
| 12-cubic-fourfold/states | 19 states, 7 first-excluded multiplicities, 2 embedding-seven comparisons; def:cubic-fourfold-costs |
| 12-cubic-fourfold/affine | 16 exact nonquartic corrections, the low-deficit quartic correction and high-deficit affine entry; eq:cubic-fourfold-affine-table |
| 12-cubic-fourfold/integral-prices | 17 rows with 3 exact costs/ceilings; eq:cubic-fourfold-integral-prices |
| 12-cubic-fourfold/reduced-prices | 17 rows with 3 exact costs/ceilings; eq:cubic-fourfold-reduced-prices |
| 12-cubic-fourfold/claims | 26 exact identities and their supporting strict inequalities |

Schema version is integer 1. Each family has exactly the fields implemented by check.py; unexpected fields are errors. Rationals are reduced strings p/q with positive denominator, or integer strings when the denominator is 1. Decimal JSON numbers, NaN, infinity, noncanonical fractions, duplicate JSON fields, missing/duplicate/extra/reordered rows and unconsumed fields are rejected. State coordinates and ceiling numerators are integers, not booleans.

The state record contains (d,e,m), moment, intrinsic cap, signed cubic parent cap and radius. The checker independently fills degree capacities binom(c+j-1,j), c=e-d, until length m. It enumerates positive states for d=2,3 and e<=6, together with the smooth curve. The first nonpositive cap for each positive-codimension pair excludes all later lengths because the minimum average degree is nondecreasing. Every radius is the unique upward multiple of 1/1000 whose d-th power reaches m.

Intrinsic and parent caps remain distinct. The m=8,9 intrinsic caps are 1/4,1/9; their cubic parent caps are 0,-1/3. Their full tails are separately checked against 2z,9z/4. The embedding-seven comparisons are 1/2,1/3 and are not admitted to this local state set.

| Price family | Columns (gamma,b,sigma,n) |
| --- | --- |
| integral | (3/4,2,0,7), (3/4,8/3,0,7), (41/50,8/3,0,7) |
| reduced | (3/4,5/2,0,7), (381/500,5/2,0,7), (39/50,5/2,0,7) |

The cost is (25/27)(n-b-sigma)+b/gamma+J_v(1/gamma,min(b,C_v)). Every ceiling N must satisfy N-1 < 10000*cost <= N; an integral scaled cost is not increased by one. The checker consumes both exact costs and all 102 printed integers, recomputes maxima and checks every paid/unpaid state at 7999/1000.

## Completeness and proof boundary

The article's lem:cubic-fourfold-certificate proves the rule used by check.py. The checker expands every path with decreasing dimension and nonincreasing embedding dimension, letting each prefix end at a point. It constructs all inequalities 0<=y_l<=...<=y_0<=a and y_i<=c_(v_(i-1),v_i). It solves every full-rank subset of active constraints by exact rational Gaussian elimination and tests all original inequalities. The full linear cost is evaluated at every feasible vertex. No recurrence or producer is imported.

A bounded nonempty polyhedron has a maximizing vertex. Every vertex has a full-rank subset of active normals, also on a proper affine subspace. Singular subsets lose no vertices. Negative caps produce empty domains; zero caps, zero deficit, segments, singletons, endpoints and ties are retained. Zero-child paths give a shorter path's value. Every dimension skip remains represented.

The final run covers 145 distinct state paths, 332 path/cap domains, 755 vertices and 114 surplus requests. Domains comprise 64 empty, 119 singleton, 40 segment and 109 higher-dimensional cases. The fresh checker report records every path domain, vertex, maximizing path and exact maximum. These are counts of arithmetic requests, not geometric centers.

The article retains the reductions over every allowed b,sigma,n: monotonicity in b, stabilized surplus above b=2, the negative sigma slope after b<=5/2+sigma/2, and the favorable dimension slope 25/27. The checker checks the rational endpoint/stabilization conditions, source-mass identities and both forced-prime orders. The geometry justifying their use is not inferred from one endpoint.

The affine correction maximum is 451/1000, the low-deficit quartic correction 477/1000, and the high-deficit entry margin 1357633/107919000>1/2000 on 1<z<=4/3, where the relevant affine slope decreases. The existing weighted cutting allowance 1/2000 and final margin 1/400000 are retained.

The integral source keeps initial mass 13/256 and final mass 51227/1250000>1/25. The two reduced-source masses are 581177533/62500000000>9/1000 and 19512331/3906250000>1/250. The latter stays positive after intrinsic normalized error 1/500. The final complete price 93433/11700, and every paid case, still pay after total weighted cutting error 1/10000 with the theorem's 1/2000 margin.

Original-degree ceilings, uniform cubic errors, separate ambient lifting, actual own thresholds, full incoming boundary, whole tangent/conductor schemes and branch-degree applicability are proved in the article. Exact fractions introduce no floating-point error.



## Public release

The mathematical data, executable checkers, supporting modules and required TeX views in this directory are byte-identical to the frozen source package. `manifest.json` records their identities and the complete family interfaces. Historical source archives, operator discussions, discovery programs and old execution receipts are retained separately and are not runtime inputs. Run the root dispatcher for a fresh complete report. Keep generated reports outside this section directory, whose inventory is fixed by the root manifest.
