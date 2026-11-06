Ancillary files of "The minimal volume of a stable surface" (Jihao Liu and Wenfei Liu).

large_mld_proof.vipr.xz    xz compression of large_mld_proof.vipr, a proof that the integer system of
                           the paper has no integer solution, in the VIPR format, version 1.0
                           (Cheung, Gleixner, Steffy, Verifying Integer Programming Results, IPCO 2017,
                           arXiv:1611.08832). 1396 integer variables, 2901 constraints (109 constraints
                           and 2792 variable bounds), 136707 derivations: 90960 linear combinations,
                           45066 integer roundings, 454 branching assumptions, 227 discharges of two
                           complementary assumptions. The last derivation is 0 >= 1 with no
                           undischarged assumption. No sol derivation, no objective cutoff, no
                           floating-point data.
local_rows.json            The 1368 local types with their invariants.
regenerate_and_compare.py  Regenerates from the formulas of the paper the 1368 local types, their
                           invariants and all 2901 constraints, with exact rational arithmetic, and
                           compares them with local_rows.json and with the constraint section of the
                           proof file, as multisets of normalized exact coefficient vectors.
viprchk/                   The checker viprchk, file code/viprchk.cpp at commit
                           99eec197a49b2bd9ad33c1f4c7afabe93caa6511 of https://github.com/scipopt/vipr,
                           unmodified, a two-function text input/output adapter for GMP rationals, and
                           a build script.
check_tables.py            Lists, from the conditions in the proofs, the parameter tuples of the three
                           tables of the appendix, compares them with tables.json, and recomputes the
                           rows from the definitions with exact rational arithmetic (the rows with an
                           unbounded interval at their first nine points; on the whole interval they
                           are polynomial identities by the paper).
tables.json                The three tables of the appendix as printed.
provenance/                How the proof file was made: see below. Also the script structure_audit.py of
                           the third check below, with its log structure_audit.log.
SHA256SUMS                 Hashes of all files here.

To check (needs a C++ compiler, GMP, xz, Python 3.9 or later):

    xz -dk large_mld_proof.vipr.xz
    sh viprchk/build.sh ./viprchk.out        # GMP_PREFIX=/opt/homebrew on macOS with Homebrew
    ./viprchk.out large_mld_proof.vipr       # prints: Successfully verified infeasibility.
    python3 regenerate_and_compare.py local_rows.json large_mld_proof.vipr
                                             # prints exact_model_match: true
    python3 check_tables.py tables.json      # prints all_match: true
    python3 provenance/structure_audit.py large_mld_proof.vipr
                                             # prints all_asm_integral_and_references_backward: true

SHA256 of large_mld_proof.vipr (1267596770 bytes):
    c3ea6f374e31789e5d559fd3929605e339536e7a93941559fc27b8d26fd50650

On 2026-11-02, on an Apple Silicon Mac with GMP 6.3.0, viprchk returned 0 and printed the line above
in about 42 seconds using about 60 MB of memory (its warning about non-ascending indices in one
derivation is harmless), regenerate_and_compare.py reported an exact match in about 3 seconds, and
structure_audit.py, which confirms that all branching assumptions have integral coefficients and
right-hand sides and that each line refers only to earlier lines, printed the line above in about
8 seconds.
viprchk rejects provenance/xhigh2_scip_toy_bad.vipr and accepts provenance/xhigh2_scip_toy.vipr.

How the proof file was made. The solver is used only to find the proof; correctness rests on the
three checks above. SCIP 10.1.0 with SoPlex 8.1.0, in exact mode enabled before reading the rational
MPS input (provenance/large_mld_model.mps), presolving disabled, one thread, explored 455
branch-and-bound nodes in 160 seconds. The raw proof of 416688 derivations was reduced to the
derivations needed for the final contradiction, and four derivations that SCIP left in a weak form
were completed by explicit multipliers of variable bounds. Finally the two constraints that no
derivation refers to and that the paper does not need (a Picard number inequality and a bound on the
sum of the numbers of exceptional curves) were removed by provenance/drop_unused_constraints.py,
which renumbers the references. The file before that step is
provenance/large_mld_proof_before_dropping.vipr.xz, with
provenance/regenerate_and_compare_before_dropping.py. The scripts and logs with the prefix xhigh2_
are records of these steps, unchanged except that local directory paths and the name of the local job
runner are replaced by [project directory] and [job runner].
