# Triple-fivefold certificates

This package checks the finite inequalities used in Section 13, thm:original-triple-fivefold. It includes every distinct numerical record used from the first-source, normal-moment, integral-transition, and supporting-slope certificates. The geometric reductions, complete original sources, both cutting errors, and separate ambient/intrinsic cutoffs remain in the article.

Use Python 3.9 or later. Only the standard library is required; there is no installation step or network access.

    python3 check.py --data data --report ../fresh-13-triple-fivefold.json

Run from this package directory, or use absolute data/report paths from any working directory. Exit 0 means every mathematical and format check passed; exit 1 means malformed data, incomplete coverage, or a failed inequality. The report identifies each family, exact counts, bounds, margins and input identities. The checker imports no conversion script, historical checker, optimizer or private project state. Predicates remain active with Python optimization enabled.

## Canonical data and strict schema

All JSON objects have exactly the fields described below and implemented in check.py. Extra keys and duplicate JSON keys are rejected. Indices, counts and dimensions are JSON integers, not booleans or floating-point numbers. Rational quantities are strings matching -?(0|[1-9][0-9]*)(/[1-9][0-9]*)?. An omitted denominator means one. Denominators are positive; signs belong to numerators. Arithmetic normalizes values, so equivalent rational spellings cannot conceal duplicate witnesses.

The data directory contains exactly four JSON files and the derived scalar TeX view. Indexed record lists contain every expected index exactly once, in order. Polygon point sets need no particular order, but duplicates and conflicting values are rejected.

| Dataset ID | File and schema | Complete content |
| --- | --- | --- |
| 13-triple-fivefold/tails | data/tails.json, triple-fivefold-tails-v1 | 30 states and 77 affine pairs; independent moment census, radii, all eligible descendants and skips |
| 13-triple-fivefold/first-source | data/first-source.json, triple-fivefold-first-source-v1 | 29 unrestricted and 10 constrained polygons, all 466 subdivided witnesses, 39 summary maximizers, 2 zero-cap records, 4 normal roots with 19 descendants each, 8 quartic low-deficit witnesses, and both source masses |
| 13-triple-fivefold/integral | data/integral.json, triple-fivefold-integral-v1 | 19 initial comparisons, 17 intervals, all 170 constrained comparisons, 17 losses and prefix masses, terminal level/bound and retained historical jump cap |
| 13-triple-fivefold/slopes | data/slopes.json, triple-fivefold-slopes-v1 | 4 height identities; 39 bins in families of 4,11,11,13; all 47 local integrals and every global/combined value |

data/tail-view.tex is a derived presentation, not a second canonical certificate. The checker parses all 30 rendered rows and compares every state, moment, cap, radius and affine pair with tails.json. This file is byte-identical to the article's retained tail-data.tex.

tails.json has fields schema, rows. Each row has row, state, moment, cap, radius, lines. The state is [d,e,mu]; lines is the complete list of [a,c] pairs for H(y)=min(a*y+c). Row IDs are 1–6, 8–20, 28–38 in increasing state order. The checker independently enumerates all states with d<=4, e<=6 and 2*T(e-d,mu)<d*mu. The finite multiplicity cutoff and complete recurrence are proved in the article.

first-source.json has fields schema, gamma, theta, constrained_order, initial_mass, reached_mass, summaries, restricted_summaries, roots, polygons, empty_domains, quartic_low.

- gamma=3997/5000, theta=803/1000, constrained_order=27/25 are independently fixed in the checker.
- A polygon has mode, row, state, cap, vertices, witnesses, upper. Mode is unrestricted or restricted. Vertices are [b,z]; witnesses are [b,z,Phi(b,z)]. Both entire sets are independently reconstructed and compared.
- An unrestricted summary has row, state, cap, upper, point; a restricted summary omits cap. Each point=[b,z] must attain the independently recomputed entire upper value.
- empty_domains contains exactly the unrestricted and restricted row-38 records with mode, row, cap. Positive deficit is impossible; the respective closed segment and point are still enumerated.
- A root has mu, E, radius, cap, descendants, surplus, upper, gap. Its 19 descendants have row, state, cap, empty, points, upper. A positive cap carries every [z,H(z)-radius*z] record; a nonpositive cap has empty=true, no points, upper="0". Every inherited cap is recomputed. The root surplus includes zero.
- quartic_low contains all [z,Phi(3,z)] witnesses on [0,1].

integral.json has fields schema, first, intervals, comparisons, losses, terminal_level, terminal_upper, jump_whole_cap.

- first uses the unrestricted summary shape for the 19 states of dimension at most three, with row 15 restricted to deficit at most one.
- An interval has interval, left, right, order, upper, point. Endpoints and proposed multiplicity-weighted fixed orders are the 17 triples displayed in the article.
- A comparison has interval, row, state, cap, upper, point. The complete expected index set is {1,...,17} x {28,...,37}; the zero-cap row 38 is separately checked at each interval. All polygon candidates are independently reconstructed. A proposed maximizer alone is never an upper-bound proof.
- A loss has interval, loss, remaining. The exact positive-part fifth-power formula, every prefix reserve, and the final positive reserve are checked.
- terminal_level=23/25 and terminal_upper=2481721/310500 are checked over all 30 states. jump_whole_cap=777619/107919 preserves a distinct original scalar datum, although the article uses the broader preceding quartic affine endpoint for each actual jump.

slopes.json has fields schema, heights, families.

- A height has mu, D, derivative, lower, buffer, in the article's notation.
- A family has mu, height, domain, count, ceiling, bins. Multiplicities are exactly 3,4,5,6. The article's heights, entire domains, counts and ceilings are independently fixed by the checker.
- A bin has bin, left, right, test, local, global_value, upper. Consecutive closed endpoints meet and the first/last equal the proved entire domain. The local list has three entries for multiplicity 3 and one for each other multiplicity, in the article's pattern order.
- Integrals are recomputed by splitting at the profile onset, anchor and every affine zero, then expanding and integrating the fourth-degree polynomial. This is independent of the source certificate's primitive. Zero shifts, zero slopes, negative slopes, zero-length intervals and endpoint ties are retained.

## Checking argument and scope

The article's lem:triple-fivefold-polygon-check proves exhaustive candidate coverage on bounded polygons, segments, points and empty domains. Independent boundary normals enumerate all vertices; horizontal equalities of every competing tail function subdivide into affine pieces. Singular line pairs are skipped as intersections without removing either inequality. Coincident boundary/equality lines retain their endpoints. A zero-area domain is not discarded.

Scalar inequalities are checked by induction in dimension with complete inherited adjacent caps and every permitted direct skip. Normal roots use the same interval endpoint rule. The integral checker reconstructs 2,886 polygon candidates for 170 comparisons and checks all 30 terminal domains. All prefix masses remain positive at one original degree. The slope checker verifies allocation vertex sets, every exact local/global value, and the whole closed slope partition.

The smallest constrained integral price gap below eight is 5471/53001000 > 1/10000; the article reserves less than 1/100000 for approximate cutting errors. The final integral mass is 51447479205198179/45000000000000000000 > 0. The reached 803/1000 mass is 18145162291861493/3125000000000000000 > 1/200. The last slope family's ceiling is below the strict target by 301/20000000. Other exact margins and counts appear in evidence/check-report.json.

These are finite arithmetic checks. They do not establish geometric hypotheses, infer five-dimensional lower mass from a three-dimensional theorem, change a source statement or certify the whole manuscript. A fresh native section verdict must judge the new article/data/checker interface.



## Public release

The mathematical data, executable checkers, supporting modules and required TeX views in this directory are byte-identical to the frozen source package. `manifest.json` records their identities and the complete family interfaces. Historical source archives, operator discussions, discovery programs and old execution receipts are retained separately and are not runtime inputs. Run the root dispatcher for a fresh complete report. Keep generated reports outside this section directory, whose inventory is fixed by the root manifest.
