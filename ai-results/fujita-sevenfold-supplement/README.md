# Computational supplement for the seven-dimensional Fujita manuscript

This public release contains the complete active finite-arithmetic data and independent checkers for the thirteen packages described in Appendix A.4. It has 75 certificate families. The section definitions, geometric applicability arguments, uniform error estimates and analytic reductions are in the article; a successful computation does not replace those proofs.

## Run

Use Python 3.10 or later with its standard library, from this directory:

```sh
python3 -B check.py --manifest manifest.json --report runs/full-001/aggregate.json
```

Use a new report path for each run. The dispatcher verifies the exact delivered inventories, runs every section, requires fresh complete reports with all original deterministic assertions, and checks that inputs remain unchanged. Reports are written under `runs/`, outside the fixed section inventories. No network, private source archive, external package, optimizer or live research system is needed. Do not use Python's `-O` option. A zero exit and `status: pass` with `full_package_run: true` mean every declared arithmetic package passed; they are not a formal proof of the paper.

For an explicitly partial run, add `--section 05-tails` (repeat for more sections). A subset is always recorded as partial. The full original run took about eleven minutes on its original machine, mostly in Section 17; current runtime depends on hardware. That historical timing is not a validation of this release. This release completed a fresh full run in 341.327 seconds on 2026-09-20 UTC. `validation-summary.json` records all thirteen successful sections and unchanged inputs; their original fresh reports and Section 09 evidence are under `validation/`.

## Contents and identity

All 122 active section files (canonical data, checker code, modules and required mathematical TeX views), plus the root dispatcher, are copied without changing any byte. The owner manifests, root manifest and README files were repackaged for public distribution. Mathematical report predicates and error allowances are unchanged. The root manifest binds every delivered section file by byte count and SHA-256; each owner manifest describes its family domains and formats.

| Section package | Families | Unchanged active files |
| --- | ---: | ---: |
| 05-tails | 2 | 3 |
| 07-capacities | 10 | 3 |
| 09-original-sources | 5 | 13 |
| 10-ledger | 7 | 3 |
| 11-quartic-threefold | 4 | 5 |
| 12-cubic-fourfold | 5 | 9 |
| 13-triple-fivefold | 4 | 6 |
| 14-double-geometry | 2 | 3 |
| 16-double-twentynine | 3 | 15 |
| 17-ambient-endpoints | 25 | 49 |
| 18-intrinsic-double | 5 | 8 |
| 19-large-deficit | 1 | 2 |
| 20-rank-two-feedback | 2 | 3 |

The original historical archives, candidate-discovery programs, duplicate expanded TeX, operational records, operator discussions and old execution receipts are retained privately. They are not imported or read by the delivered checkers. Their omission does not remove any active arithmetic input. Documentary source identifiers and hashes remain in the manifests to record origin; they are not runtime paths or certificates of mathematical truth.

Sections 09 and 16 require their delivered `views/`; Section 17 requires all `modules/*.py`. The TeX files in Sections 12 and 13 `data/` are active mathematical inputs and must also be retained. Section 09 generates a new `.evidence/` directory beside its report. Generated reports should never be inserted into a fixed section directory.

Section 09 produces a large derived arithmetic trace. The fresh trace is compressed as `validation/09-original-sources/report.evidence.tar.gz`; extract it in that directory to restore the report-relative `report.evidence/` paths. It contains only newly computed numerical output, no original historical archive. Re-running the checker also regenerates it. Historical module docstrings describe the origins of the unchanged code; use the portable commands documented here.
