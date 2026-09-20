# Ambient profile endpoints: certificate supplement

This directory contains the 25 complete data objects used in the ambient profile endpoint proofs. Each stable identifier is `17-ambient-endpoints/<key>`, with its payload in `data/<key>.txt`. The article retains the geometric arguments, domains, competing functions, analytic slope tails, source errors, and conclusions. Definition `def:ambient-endpoints-finite-witness` and Lemma `lem:ambient-endpoints-certificate-correctness` give the common acceptance argument; each data reference follows its specific formulas and margin.

Run from this directory with Python 3.10 or later:

```sh
python3 check.py --data data --report ../fresh-17-ambient-endpoints.json
```

Only the Python standard library is required. Mathematical decisions use `fractions.Fraction`, integer arithmetic, exact polynomial integration and exact Bernstein coefficients. Do not use Python's `-O` option: assertions are part of the checks, and the driver rejects disabled assertions. The fresh checker report records its results. The complete check takes several minutes; progress is printed after each family and during the larger trees. Arithmetic caches are cleared between families. Exit zero from `check.py` means that all 25 objects passed; a missing object, duplicate record, malformed tree, failed domain test, false integral or insufficient margin causes a nonzero exit and a failure report. There is no successful partial mode in the delivery driver.

## Data and decoding

All bytes inside every original verbatim block were preserved. Blocks belonging to one certificate were joined in their original order with one newline between blocks. TeX wrappers, font instructions, and surrounding proof prose are not payload. The manifest records the payload identities and family definitions; the original-to-canonical transcription record is retained separately. These identities establish transcription, not mathematical correctness.

Physical wrapping has no mathematical meaning. A continued record beginning with indentation is concatenated without inserting a character: a wrap may divide a word or number. Readers then recognize their exact record grammar and separators. Deterministic trees use balanced parentheses and ignore whitespace inside the tree; the reserve's indented replacement table has a separate parser. The initial-surface object begins directly with rational index zero because the first dictionary heading was outside its original verbatim block; its reader enters rational-dictionary mode explicitly. All indices, roots, children and leaf traces are checked, not inferred from historical acceptance.

The formats remain distinct:

| Format | Objects | Exact acceptance interface |
| --- | --- | --- |
| Linked profiles | Five `reserve-*`, three `point-*` | Prescribed closed slope grids; every binary split; linked heights with error `1/100000`; local rational/polynomial dictionaries; complete traces. Polynomials include the factor 7. Strict complete-row upper denominators are `10^10`. |
| Retained parameter | `retained-parameter` | Fourteen slope intervals and all sixteen enlarged incidence intervals; 1,099 cells; Bernstein coefficient bounds over `10^6`, contribution bounds over `10^8`; both incoming tests. This is the **early 1099** certificate. |
| Deterministic trees | Five `initial-fivefold-*` | Reconstruct each prescribed bisection, integrate every competing bound on the prescribed time subdivision, and test its complete leaf criterion. The reserve includes all three original trees and all fifteen replacement subtrees. |
| Affine source traces | Four easy-source objects and `source-initial-threefold` | Physical root parameters, explicit splits, strict geometric empty tests where present, exact branch validity, and complete adjacent time cells. Polynomials omit the factor 7; strict cell denominator is `10^9`. |
| Indexed/extended source traces | `source-initial-surface`, `source-quintuple-surface` | Full four- or six-coordinate domains; dictionary or literal records; every empty and nonempty leaf; the source-specific retained terms and full mass criterion. The quintuple case checks all three horizons and the residual root degree. |
| Distinct divisors | `distinct-divisors` | All 31 cells, exact capacity domains, polynomial positivity, four histories, and their corresponding full prime-mass ceilings. The later distinct divisor is allowed to be singular. |
| Persistent prime | Three `persistent-*` objects, seven families | Complete written trees and literal witnesses; physical version 1 and normalized version 5 are kept separate. Reconstruct every effective box, including the interior kink of the linked slope cap and fixed coordinates. Check each cell and its recorded mass at its recorded horizon. |

Every family has its own entry in `manifest.json`, with its canonical data, source identifiers, mathematical specification and exact data hash. Original sources are retained separately as private archival material; the executable modules do not read them.

## Source histories and margins

The five embedding-seven reserve objects prove loss strictly below `63/64`, earning `V(29/25)>1/64`. Their later nonsmooth tests have target `999/1000` through `13/5`. The separate initial-fivefold source history earns `V(27/25)>6/25` and then `V(29/25)>1/500`; the latter uses target `499/500`, including fifteen replacements in the generic tree. Its nonsmooth, incidence and payment tests use `9999/10000`, with their own incoming cap `499/500`. The eight retained-incidence slabs contain 10,756 leaves; their largest written upper bound is `9998671/10^7`, leaving `329/10^7` to their stated threshold.

The early retained-parameter incoming maxima are `78813928/10^8<789/1000` and `84492074/10^8<169/200`. These are separate tests with different outgoing/mass arguments. The **final 1031** source calculation belongs to the earlier original-source theorem and is not substituted here.

The regular-source trace has actual target `19/20`, stronger than its old generic footer `99/100`. The initial-surface and quintuple-surface closing tests depend on the method: time-29 integral rows require `<999/1000`; certain earlier-horizon rows add their full mass and require `<1`. They must not be assigned the stronger target. The persistent NE maximum is `9998993/10^7`, leaving `7/10^7` to `9999/10000`; every other persistent family also has its own positive gap. Exact maxima, method counts, cell counts and minimum rounding slack are recorded in the full report.

These numerical gaps are not the complete-series approximation errors. The article independently fixes the original observing degrees and geometric error budgets before centers or computing places are chosen. In particular the `25` divisibility requirement belongs to the eligible observations, while the earned source dimensions remain valid for every sufficiently large output integer.

## Independence and scope

The portable readers consume the delivered text directly. Arithmetic modules were adapted from independent exact audits after inspecting their formulas. Their old source-file, submission-log, producer-JSON and receipt dependencies were removed. The persistent arithmetic audit formerly read saved boxes; the new reader instead reconstructs the prescribed tree and physical/normalized map and checks every *written* witness. A generator's cached assertion is never an input. The manifest lists each adaptation, including the separate former wrappers, so this distinction can be reviewed.

The full check is an executable audit of these 25 finite certificates and the scalar/degenerate identities explicitly implemented in their modules. The payment module also checks the two degree-12 Bernstein bounds printed in the article. It does not establish the geometric hypotheses, replace the article's proofs or verify the whole historical fact graph. In particular earlier-section source certificates and other unchanged short inline calculations remain with their stated proofs. Native section verification must assess the entire article-to-data-to-checker interface.



## Public release

The mathematical data, executable checkers, supporting modules and required TeX views in this directory are byte-identical to the frozen source package. `manifest.json` records their identities and the complete family interfaces. Historical source archives, operator discussions, discovery programs and old execution receipts are retained separately and are not runtime inputs. Run the root dispatcher for a fresh complete report. Keep generated reports outside this section directory, whose inventory is fixed by the root manifest.
