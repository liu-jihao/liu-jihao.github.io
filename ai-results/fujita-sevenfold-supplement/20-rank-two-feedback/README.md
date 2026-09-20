# Rank-two finite source feedback

This supplement carries two canonical certificate families for `thm:rank-two-twentyseven` and `prop:rank-two-source-induction` in Section 20. The article retains the complete geometric proof, full symbolic and principal modules, continuous reductions for all factor orders and characters, original-source induction, both degree cutoffs, all old coefficients and both cutting errors. Its local correctness lemma is `lem:rank-two-feedback-finite-certificate`.

From this directory, using Python 3.9 or later with only its standard library:

```sh
python3 check.py --data data --report ../fresh-20-rank-two-feedback.json
```

## Data and exact conventions

`data/` must contain exactly `state_rows.json` and `mesh_rows.json`. Each is a JSON array in the prescribed increasing order. Integers use JSON integer syntax; exact rational values use reduced strings such as `"283/200"`, with positive denominator, or `"2"` for integers. No floating-point, nonfinite, duplicate-key, unknown-field, missing-record or extra-record input is accepted. Every field is compared with an independently reconstructed exact value, including secondary witnesses and all maximizing names.

`20-rank-two-feedback/states` consists of 28 objects in lexicographic state order. Each has exactly these fields:

| Field | Meaning |
|---|---|
| `state` | Integer triple `[d,e,m]`. |
| `T` | Least Artinian first moment. |
| `intrinsic_cap` | `d-2*T/m`; this remains the cap used in descendant recursion. |
| `root_cap` | Minimum with `(8-m)/2` for a first fourfold, otherwise the intrinsic cap. |
| `radius` | The article's rational upper radius, with `radius**d >= m`. |
| `lambda_price` | `kappa(lambda)+J(1/lambda,min(18/5,root_cap))`, with the separate first-deficit restriction to one only for `(3,6,4)`. |
| `theta_price` | Complete scalar price at `theta`, except `(4,6,3)`, which carries `kappa(theta)+(1317/1000-1/theta)*7/3`. |
| `theta_kind` | Exactly `complete scalar tail` or, for that cubic, `cubic affine endpoint`. |
| `descendants` | Every allowed lower-dimensional state, in increasing order, as objects with exactly `state` and `cap`; 170 records in total. |

The checker enumerates states by the article's finite first-moment bound, then applies every positive-cap condition. It reconstructs all dimension skips and adjacent caps and evaluates the complete support recurrence. It checks every moment, radius, cap, successor and both price columns. It does not transplant the distinct `J(50/49,C)` column or the `theta=21/20` certificate from other sections. The article's small state table is a derived view of these records.

`20-rank-two-feedback/mesh` consists of exactly 235 objects with `j=1731,...,1965`. The implicit right endpoint is `(j+1)/2000`, giving 236 endpoints in total. Each record has exactly `j`, `gamma`, `Y`, `expressions`, `L`, `E`, `maximizers` and `rounded_prefix_reserve`. Here `gamma=j/2000`; `Y` has exactly the keys `"3"`, `"4"`, `"5"`; `L` is the maximum of zero and all included expressions; `E` is the upper integer; and the prefix reserve is `1-2*(1731/2000)^5-sum(E_1731,...,E_j)/10^9`.

All eight expression keys map to the article's seven family rows as follows. The two smooth cubic expressions are separate candidates.

| Key | Inclusion condition | Article expression |
|---|---|---|
| `smooth_cubic_extreme` | `Y_3 < 7/3` | `C_E`. |
| `smooth_cubic_balanced` | `Y_3 < 7/3` | `C_B`. |
| `other_cubic` | `Y_3 < 2` | `Phi(4 R_3/5)+Phi(8 R_3/5)`. |
| `distinct_quartic` | `Y_4 < 3/2` | `Phi(21 R_4/25)+Phi(63 R_4/25)`. |
| `cartier_quartic` | `Y_4 < 2` | `2 Phi(8 R_4/5)`. |
| `power_quartic` | `Y_4 < 2` | `2 Phi(gamma Y_4/(4-Y_4))`. |
| `cartier_quintic` | `Y_5 < 3/2` | `2 Phi(42 R_5/25)`. |
| `power_quintic` | `Y_5 < 3/2` | `2 Phi(gamma(1+Y_5)/(4-Y_5))`. |

The constants are `B_*=18/5`, `r_0=25/27`, `lambda=1731/2000`, `theta=983/1000`, and interval width `1/2000`. The article defines `kappa`, `A_mu,C_mu,D_mu,R_mu,Y_mu`, `Phi`, `chi`, `xi_*`, `C_E` and `C_B`; `check.py` reconstructs precisely these formulas. Each record must supply the exact included key set, all rational values, the maximum and every tie. There are 1266 included expression records. A supplied maximizer alone cannot pass.

The rounding rule is ordinary ceiling: `E=ceil(10^9 L)`. Every delivered row satisfies the stronger strict bracket `E-1 < 10^9 L < E`, so floor-plus-one agrees for these nonintegral scaled maxima. The checker tests both sides, every competing value, every prefix reserve, sum `28148081`, and final reserve `8402435542349/16000000000000000`, whose excess over `1/2000` is `402435542349/16000000000000000`. Positive-part signs, ties and equality at interval endpoints are exact. Empty candidate sets have maximum zero; the actual delivered intervals all have positive maxima. No sampled power loop occurs.

The checks also retain the initial and terminal whole-price inequalities, scalar intercepts, continuous-reduction denominator signs and endpoint constants, direct quartic affine entry and cubic descendant margins. Exact arithmetic has no floating-point error. Rounding only increases the loss. The article separately fixes `120*(C_0+235*C_1')/p < 1/4000` and the ordinary restriction cutoff before the same 2000-divisible degree; its geometric injection then meets the correctness lemma's hypotheses. The total weighted cutting allowance remains below `1/1000` and the final point margin exceeds `1/400000`. Arithmetic checking does not prove normality, matrix factorization, geometric realization, the full-source injection, or the uniform analytic error estimate; those arguments and imported hypotheses remain in the article.



## Public release

The mathematical data, executable checkers, supporting modules and required TeX views in this directory are byte-identical to the frozen source package. `manifest.json` records their identities and the complete family interfaces. Historical source archives, operator discussions, discovery programs and old execution receipts are retained separately and are not runtime inputs. Run the root dispatcher for a fresh complete report. Keep generated reports outside this section directory, whose inventory is fixed by the root manifest.
