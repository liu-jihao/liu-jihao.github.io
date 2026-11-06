# Large-mld certificate — 2026-10-31 12:15 UTC

## Mathematical status

The chamber exclusion is **accepted as fact 1f3b8f21353c6e92**:

> For a complex projective integral klt surface with ample K_X and
> mld(X)>=5/46, one has K_X^2>1/48983.

The submission is anchored to shared finding **a392e535fe5b4b8e**. Verification
accepted it in `verifier/runs/20261031T121402Z_fc3452a6216e/log.md`.
Use **1f3b8f21353c6e92** as the predecessor for subsequent project arguments.
The full universal volume/equality target remains unsolved because the small-mld
branch remains separate.

The earlier continuous-relaxation obstruction is accepted as
**1844fb29d81dd125**. It gives a fractional point at v=1/49000 and explains why
integer multiplicities/dimensions matter. It is not in conflict with the integer
exclusion certificate below.

## Decisive artifacts

- Submitted theorem: `xhigh2_large_mld_exclusion_statement.md`.
- Written reduction and certificate rules: `xhigh2_large_mld_exclusion_proof.md`.
- Explicit rational proof certificate:
  **`xhigh2_scip_exact_45.complete.vipr.opt`**.
- Certificate size: **1327627770 bytes**.
- Certificate SHA256:
  **42a2669f4670ad31f1ae1869db820462ae847312b2ec7a530261bd614e743df5**.
- Losslessly compressed archive:
  `xhigh2_scip_exact_45.complete.vipr.opt.gz`, **47184015 bytes**.
  Decompression was checked against the same SHA256. The original remains.
- Independent exact checker result: `xhigh2_vipr_final_check.json` and `.log`.
  Return code zero; “Successfully verified infeasibility.”; 89.9096 CPU seconds.
- Independently regenerated input audit: `xhigh2_vipr_model_audit.py` and `.json`.
  Exact coefficient match: 1368 local rows, 1396 integer variables,
  111 non-bound constraints, 2792 bounds, 2903 total original constraints.
  The audit was also rerun directly on the final `.complete.vipr.opt` artifact
  and passed, so this identification does not rely on header preservation by
  the optimizer.

The proof has **136707 retained derivations**. It contains no feasible-solution
or objective-cutoff premise (`sol`). All derivations use rational linear
combinations, integer rounding, integral branch assumptions and exhaustive
binary disjunction discharge. The final absurdity has no undischarged assumptions.

Supplementary structure audit (`xhigh2_vipr_structure_audit.json`) also passed:
454 integral branch assumptions; 90960 linear combinations; 45066 rounding
steps; 227 binary disjunction discharges. All references are backward and
distinct. The compiled checker source matches the pinned upstream byte for
byte, SHA256 1ac2c9f174f2eaf3420f4fc9a6ffebcc2583196be497e9ae5a9a52ff8dfcde56.
An intentionally invalid rounding certificate is rejected by this checker,
while the valid rational-input toy certificate is accepted.

The descriptive parse `xhigh2_vipr_explain.json` identifies 177 distinct branch
variables, including multiplicities and plurigenera. It is only a navigation
aid to the certificate; no shorter independent symbolic exclusion has been
extracted from it.

## Model and mathematical scope

The independent regeneration enumerates every parameter tuple of accepted
36a63bb3547e6223, checks Qb=k and the exact local corrections/group orders/cost,
then applies exactly the mld and local-cost filters. It matches every local row
and every rational constraint. **Even individual RR-zero rows are retained.**
No packet deletion, support cap, rank-one assumption or guessed coefficient
bound is used.

The model has P2,...,P17=0 and binary P18,...,P45. It retains the free-rank
ledger, actual-group BMY, length, strict-derived local upper bounds and global
cost, both volume bounds including the accepted positive floor, power section
injections, gcd-saturation cuts involving caps through degree156, and all
product cuts a+b<=45.

Geometric predecessors in submission:
36a63bb3547e6223, 7808ce58c8bab2e0, f81c3dfc71b97399,
4e9502305388be64, a49bfb45b25bbe79, 3bafca4b69c87d88,
e883a2cebaee22c6, dadf467fc1e53459.

## Reproducibility and toolchain

The earlier PySCIPOpt wheel lacked exact support. The official Linux SCIP10.1.0
bundle **does** include exact solving and proof export. It was downloaded into
`xhigh2_scip_exact_bundle/` without any system installation. The official
archive SHA256 is 330a884428af4d644370cdf36ba405b879a9e17341c5ad2646917ba8972767eb.

- Model export: `[job runner] computation/xhigh2_exact_gmi_probe.py 45 export`.
- Exact solving: `[job runner] computation/xhigh2_scip_exact_run.py`.
  Exact mode is enabled before input; presolve disabled; one thread; 180s limit.
  It solved infeasible after 455 nodes and 160.17 solver seconds.
- Raw proof: `xhigh2_scip_exact_45.vipr`, 3915388670 bytes.
- `xhigh2_vipr_trim.py`: retains ancestors and necessary global bounds.
- `xhigh2_vipr_complete.py`: completes the four remaining weak derivations
  with explicit exact bound multipliers. It does not validate the whole proof.
- `xhigh2_vipr_final_check.py`: runs official trimming then official checking.
  `checkonly` reuses the final trimmed file and reruns the independent checker.
- `xhigh2_vipr_model_audit.py`: regenerates the full mathematical model and
  compares the certificate's base constraints exactly.

VIPR source is pinned to scipopt/vipr commit
99eec197a49b2bd9ad33c1f4c7afabe93caa6511 under `xhigh2_vipr/`.
The unchanged upstream `viprchk.cpp` and `viprttn.cpp` were compiled with an
isolated Zig0.13 compiler because the host has no C++ compiler. GMP text I/O
requires a two-function libc++ compatibility shim in
`xhigh2_vipr/build/gmp_stream_shim.cpp`; its exact contents were included in
the verifier submission, and it changes no proof logic. The upstream checker
also passed a simple rational-input toy certificate.

At this pinned commit, viprttn returns -1 even after successful output because
its local `rs` is initialized to -1 and never updated. Its exit status was not
used as proof. The independent viprchk returned zero on the resulting file.

All computational runs used `[job runner]` DEFAULT, never heavy. No system
package, credential, backend, or service was changed. All compute processes
have finished. Several large intermediate proof files remain for traceability;
the definitive certificate and hash above are the ones submitted.
