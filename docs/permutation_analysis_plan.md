# RMC Permutation Null Analysis

> **Branch:** `rp/permutation_null_analysis`, cut from `kc/mane_plus_clinical_freeze`.
> First action on approval: create the branch and commit this plan to it as
> `docs/permutation_analysis_plan.md`. (Chosen over `temp/`, which is currently
> untracked and reads as scratch — say if you'd rather it sit next to the working
> document.)

## Context

A reviewer of the RMC paper raised two statistical objections to the interpretation of MCR-level o/e ratios (recorded in `temp/RMC permutation analysis working document.md`):

1. **Variance artifact.** MCRs are necessarily no larger than transcripts, and smaller regions have higher variance in realized o/e. So MCR o/e should have heavier tails than transcript o/e *even under a null of no regional variation in constraint* — which would partly explain Figure 1d and the finding that depleted MCRs tend to be small.
2. **Winner's curse.** The break search explicitly partitions transcripts to *maximize* o/e differences between segments, biasing discovered MCRs toward o/e values more divergent from the transcript mean than the truth.

Both objections are mechanically real in this codebase. The break call is a max over every candidate breakpoint compared against a single-test threshold (`rmc/utils/constraint.py:1032-1034`, `P_VALUE = 0.001`), and there is **no multiple-testing correction anywhere in the repo** — not across breakpoints within a section, not across transcripts.

The reviewer's proposed remedy: permute position labels within each transcript, rerun the pipeline, and see what MCRs and o/e ratios the null produces.

**Intended outcome:** a quantified null distribution that lets the "MCR size" supplement section state how much of the size–constraint relationship is statistical. Published RMC calls, MPC, and missense badness are **unchanged**. This is a diagnostic, not a recalibration.

## Settled decisions

| Decision | Choice |
|---|---|
| Deliverable | Diagnostic only; no re-filtering of published calls |
| Null construction | Joint row permutation — shuffle the `(observed, expected)` pair together to a random locus within each section |
| Search fidelity | Full procedure: all iterative single-break rounds **and** the simultaneous two-break search |
| Allocation | 100 genes × 100 permutations, sampled from freeze 2's searched set |
| Stratification | Proportional: 36 genes from the RMC-positive pool, 64 from the RMC-negative pool |
| Sampling frame | All searched transcripts (not restricted to those above `MIN_EXP_MIS`) |
| Basis freeze | Freeze 2; permutation outputs written under **freeze 99** |
| Rounds | Run to natural exhaustion |
| Packing | All permuted sections in one packed run per round |
| Randomness | `hl.rand_unif` with a fixed recorded seed |
| Region size metric | Coding base pairs |
| Real comparison arm | Matched 100 genes (primary) + genome-wide freeze 2 (secondary) + an unpermuted control replicate |
| Headline statistics | Size-stratified excess **and** a tail-mass table |
| Interpretation threshold | Not pre-committed; interpret after seeing results |

Two properties worth noting because they simplify the design:

- **Transcript-level o/e is exactly invariant** under this permutation — the joint shuffle preserves each section's `observed` and `expected` totals. Only the regional decomposition changes, so the null isolates precisely what the reviewer asked about.
- **Chi-square ties are not a concern.** The statistic depends on *cumulative expected*, a continuous float floored at `1e-09` (`rmc/utils/constraint.py:698`), so it strictly increases per locus. Exact ties are measure-zero.

## Key structural facts this plan relies on

- The search between `prep-constraint` and `finalize` is **self-contained**: it joins only on `section` and never reads `transcript_ref`, `transcript_cds`, or the VEP context. The only external joins are in `create_constraint_prep_ht` (before) and `format_rmc_browser_ht` (after).
- `section` is the **sole carrier of transcript identity** during the search — `create_constraint_prep_ht` drops the `transcript` field at `rmc/utils/constraint.py:720`. `merge_rmc_hts` re-derives it as `section.split("_")[0]` (`:1484`).
- The scan is **dict-valued, keyed by section** (`rmc/utils/constraint.py:239`), so many sections sharing loci is already normal — that is how overlapping canonical transcripts work today. 10,000 packed sections sits below production freeze 2's ~19,000.
- **All six section parse sites are positional `split("_")[0|1|2]`.** The permutation label must therefore use a non-underscore separator.
- Round N+1 receives **only break-found sections**, split at their breakpoints (`rmc/pipeline/regional_constraint.py:379-399`). Two-break cost is `Σ nₛ²`, so splitting a section into halves gives `n²/4`. Round 1 is the peak; no extra batching is needed for later rounds.

## Phase 0 — Gene sampling

New code in `rmc/utils/permutation.py`.

Partition freeze 2's searched transcripts into RMC-positive and RMC-negative pools using **`rmc_browser.versions[2]`** (`rmc/resources/rmc.py:553`) rather than `rmc_results`. `rmc_browser` is keyed by `transcript` with a `regions` array, so the positive pool is a direct read, and its globals — written by `add_globals_rmc_browser` (`rmc/utils/constraint.py:2226`) — carry `transcripts_no_rmc` inside the `transcripts` / `all_transcripts` structs, giving both pools and the full searched denominator from one resource.

Note the access pattern: `transcripts_no_rmc` is a field *of those global structs*, not a top-level global — reading it as one raises `AttributeError`, which is the documented bug in `create_rmc_release_downloads` (`rmc/utils/constraint.py:2869-2871`).

Sample 36 positive and 64 negative with a recorded seed. **Persist the sample as JSON** — the transcript lists, the seed, the pool sizes, and the realized RMC-positive proportion — rather than as a `HailExpression`, so it is diffable, reviewable, and readable without a Hail session.

## Phase 1 — Permuted constraint prep (freeze 99)

Write directly to `constraint_prep.versions[99].path`, reusing freeze 2's prep as input so all upstream work (context, mu, expected) is inherited unchanged.

For each sampled transcript and each `p` in `0..99`:

- Section string becomes `{transcript}-perm{p:03d}_{start}_{stop}` — **`-perm`, never `_perm`**. A `_` here silently drops the label at `merge-single-simul` and collapses all 100 replicates onto one section key.
- Collect the section's `(observed, expected)` structs into an array, sort by a random key (`hl.rand_unif` with the fixed seed, or `hl.shuffle`), and zip against the section's loci in ascending order. Explode back to one row per locus.
- Add one unpermuted control replicate per gene, labelled `{transcript}-ctrl_{start}_{stop}`, with values left in place.

Key by `("locus", "section")` to match the production schema. Output is roughly 15M rows (100 genes × ~1,500 loci × 101 replicates), against production's ~28M.

## Phase 2 — Break search

### Prerequisite: three blockers from rmc_production#215

[Issue 215](https://github.com/broadinstitute/rmc_production/issues/215) documents that the two-breaks steps assume every round has a non-empty result of each kind, which is why the break search has always needed a human watching each round. **All three failure modes are guaranteed to hit a run-to-exhaustion permutation job** — under the null, later rounds will routinely have no over-threshold sections and no simultaneous breaks at all. For calibration, the issue records that **freeze 2 needed 16 rounds** and freeze 3 needed 10+; at six invocations per round that is ~96 manual steps, which is not viable here.

These must be fixed before Phase 2 can run unattended. All three are pre-existing production bugs, so the fixes are useful beyond this analysis.

1. **`split-sections` writes only the size classes it has sections for.** Once a round has zero over-threshold sections, `run_batches_dataproc.py --run-all-sections` dies reading a missing `sections_to_simul_over_threshold.he`. Fix in `rmc/utils/simultaneous_breaks.py:133-149`: always write both files, writing an empty set for an empty class.
2. **`merge_hts.py` raises when no section had two simultaneous breaks.** `merge_simul_break_temp_hts` raises `DataException: All temp tables had 0 rows` (`rmc/utils/constraint.py:1221-1224`), but this is a legitimate outcome in later rounds — and `merge-single-simul` already handles the resulting missing `final_results/merged.ht` through its `simul_exists` branch (`rmc/pipeline/regional_constraint.py:264-273`). Fix: add `--allow-no-results` so the check stays useful interactively but a driver can pass it, log, and exit 0 without writing.
3. **A round cannot be retried after a partial failure.** `run_batches_dataproc.py:136-140` raises if raw results already exist, and the issue's author confirmed by testing that `--read-if-exists` does *not* bypass it — the only workaround is deleting `raw_results/` by hand. Fix: give the raw-results write an escape hatch, either treating a complete existing output as a no-op or adding a real `--overwrite` that clears and recomputes.

Fix 3 pairs with the deterministic-grouping change below: sorting sections before grouping is what makes an idempotent re-run meaningful, since otherwise `simul_break_dataproc_7.ht` may hold a different section set on the retry than it did on the original attempt.

### Round loop

Add a `run-rounds` subcommand to the new driver that executes the cycle below until `merge-single-simul` writes no merged break-found table, emitting per-round counts (sections broken / not broken, split sizes, two-break sections) as a JSON artifact. This is exactly the encapsulation issue 215 lists as a follow-up, and it is reusable for production runs.

Each round, against `--freeze 99`:

1. `regional_constraint.py search-for-single-break --search-num N`
2. `two_breaks/prepare_transcripts.py create-grouped-ht --search-num N`
3. `two_breaks/prepare_transcripts.py split-sections --search-num N`
4. `two_breaks/run_batches_dataproc.py --search-num N --group-size 100`
5. `two_breaks/merge_hts.py --search-num N`
6. `regional_constraint.py merge-single-simul --search-num N`

Then `regional_constraint.py finalize --freeze 99` **without** `--filter-outliers` and **without** `--filter-to-canonical` — both filter against real transcript sets that labelled IDs won't match.

`--group-size` must be set explicitly. The group loop at `rmc/pipeline/two_breaks/run_batches_dataproc.py:129-156` is serial and each call is ~5-6 Hail passes; omitting it means one group per section, i.e. thousands of sequential jobs.

## Phase 3 — Region assembly with an AA-free exon snapper

New code in `rmc/utils/permutation.py`. This replaces `format_rmc_browser_ht`, which cannot be used: its exon-snapping step is an amino-acid *repair* whose trigger is `hl.is_missing(ht.start_aa)` (`rmc/utils/constraint.py:2115-2131`), it reads coverage stats unconditionally (`:2779`), it raises on unknown transcripts (`:1735`), and its `group_by("transcript")` (`:2845`) would collapse all replicates into one row.

The need is real regardless: region starts are set to `breakpoint + 1` (`rmc/pipeline/regional_constraint.py:361`), which lands in an intron whenever a breakpoint falls on an exon's last coding base.

Steps:

1. Read `rmc_results.versions[99]`. Split `transcript` into `real_transcript` and `perm_label` on `-`.
2. **Lift `_create_exon_position_ht` to module scope.** It currently sits nested inside `fix_region_start_stop_aas` at `rmc/utils/constraint.py:2005`, but it takes `cds_ht` as a parameter and contains no AA logic — it is reusable verbatim once importable.
3. Snap boundaries geometrically against `transcript_cds` joined on `real_transcript`: any region start not on a CDS position moves to the next coding position, any stop to the previous. Do **not** carry over the `-1`/`+1` offsets at `:2122`/`:2132` — those encode a freeze-7-specific off-by-one assumption (`:1953-1958`), not a general snap.
4. Compute `coding_bp` by interval intersection against `transcript_cds`, not by genomic span — a region crossing two exons otherwise counts the intron. Use interval arithmetic rather than `explode_intervals_to_loci` (`rmc/utils/constraint.py:124`) at this row count.
5. Annotate `section_oe = section_obs / section_exp`.

**Gene-replicates with no break at all are not in `rmc_results`.** `merge_rmc_hts` skips round 1 (`rmc/utils/constraint.py:1454`) and requires ≥2 rounds (`:1443-1447`); replicates that never broke are recorded by `create_no_breaks_he` (`:1118`). Both sources must be unioned to compute the break-discovery rate correctly.

## Phase 4 — Comparisons (design now, build after generation)

All four use the Phase 3 table. The matched real arm is freeze 2's regions for the same 100 genes; the secondary arm is all genome-wide freeze 2 regions.

1. **o/e marginal distribution.** Real vs null MCR o/e, overlaid. Matched primary, genome-wide secondary. This is the direct Figure 1d rebuttal.
2. **o/e vs region size — headline.** Bin by `coding_bp` (deciles of the pooled null size distribution, reported with per-bin counts). Within each bin report **size-stratified excess**: real tail mass minus null tail mass.
3. **Null break-discovery rate.** Fraction of gene-replicates with ≥1 region, against the real 36%. Report separately for the RMC-positive and RMC-negative strata, which differ in power.
4. **Most-extreme o/e per gene.** Distribution of each gene-replicate's minimum region o/e, real vs null. Isolates winner's-curse magnitude.

Plus a **tail-mass table**: fraction of regions with o/e below 0.2 / 0.4 / 0.6, real vs null.

## Code changes

| File | Change |
|---|---|
| `rmc/resources/rmc.py:29` | Add `99` to `FREEZES` with a prominent comment marking it as the permutation null. Leave `CURRENT_FREEZE = 2`. Without this, the whole search runs and then `constraint_prep.versions[99]` dies on a bare `KeyError`. |
| `rmc/pipeline/two_breaks/run_batches_dataproc.py:111-122` | Sort `sections_to_run` before grouping. It is currently `list(frozenset(...))`, whose order is not stable across runs, while output files are named by index — so a partially failed round cannot be safely resumed and the `file_exists` guard at `:137-141` can raise on sections that were never processed. |
| `rmc/utils/simultaneous_breaks.py:133-149` | **Issue 215 case 1.** Always write both size-class `.he` files, empty set included, so `--run-all-sections` works on every round. |
| `rmc/utils/constraint.py:1221-1224` + `rmc/pipeline/two_breaks/merge_hts.py` | **Issue 215 case 2.** Add `--allow-no-results` so an empty simultaneous-break round logs and exits 0 instead of raising. |
| `rmc/pipeline/two_breaks/run_batches_dataproc.py:136-140` | **Issue 215 case 3.** Give the raw-results write an escape hatch so a failed round can be retried without hand-deleting `raw_results/`. |
| `rmc/utils/constraint.py:2005` | Lift `_create_exon_position_ht` from a nested function to module scope. No behaviour change. |
| `rmc/utils/permutation.py` | **New.** Sampling, permuted prep construction, geometric snapper, comparison functions. |
| `rmc/pipeline/run_permutation.py` | **New.** Argparse driver with subcommands for `sample-genes`, `prep-permuted`, `run-rounds`, `assemble-regions`, `compare`. |

The four issue-215 / determinism fixes touch shared production code. They are behaviour-preserving on the happy path (new opt-in flags, an empty-file write, a sort) and each should be a separate reviewable commit, since they are independently useful to the production break search.

## Landmines

- **Round-1 cache short-circuit.** `rmc/pipeline/regional_constraint.py:140-152` reads `{CONSTRAINT_PREFIX}/4.1/{freeze}/constraint_prep.ht` and, if present, skips `process_sections` entirely, reusing cached chi-squares. It is dead today only because of a path mismatch. **Verify that path is absent for freeze 99 before round 1.**
- **Unscoped temp path.** `rmc/utils/simultaneous_breaks.py:455` names `_tmp_repart.ht` by `section_group[0]` alone — no freeze, no round, `overwrite=True`. Do not run freeze 99 concurrently with any other pipeline run.
- **Merge glob.** `merge_simul_break_temp_hts` selects temp tables by the substring `"dataproc"` (`rmc/utils/constraint.py:1196-1218`). Do not mix `--run-all-sections` with `--run-sections-over/under-threshold` within a round, or sections get double-counted.
- Labels must contain no `_` and no `/` — section strings become GCS object names in the two-break temp paths.

## Verification

**Smoke test first: 5 genes × 3 permutations, plus the same 5 genes unpermuted.** Run the full Phase 1–3 chain end to end before committing to the 10,000-section run.

Assertions:

1. **Totals preserved.** For every permuted section, `section_obs` and `section_exp` must equal the corresponding real transcript's values exactly. A permutation that changes either is a bug — this is the single strongest correctness check available.
2. **Labels survive the search.** After `merge-single-simul` in round 2, confirm every section string still carries its `-perm###` label and that `split("_")` yields exactly 3 parts with `[1]` and `[2]` numeric.
3. **Unpermuted control reproduces production.** The `-ctrl` replicate must yield the same regions as published freeze 2 for those genes. If it does not, the harness is not faithful to the pipeline and nothing downstream is trustworthy.
4. **Snapper correctness.** Every snapped region start and stop must fall on a CDS position of its real transcript; `coding_bp` must be ≤ genomic span and equal to it for single-exon regions.
5. **No silent NA.** Count missing `start_coordinate` / `stop_coordinate` / `coding_bp` after snapping — the production path fails this way silently.
6. **Round-1 cache absent.** `gsutil ls gs://regional_missense_constraint/constraint/4.1/99/constraint_prep.ht` must return nothing before the first round.
7. **Empty-round handling.** The smoke test must run to exhaustion, which means its final rounds will exercise all three issue-215 paths — a round with no over-threshold sections, a round with no simultaneous breaks, and a deliberate mid-round kill followed by a retry. If the loop cannot survive those unattended on 5 genes, it will not survive 100.

Record the master seed, the 100-gene sample, and the per-round section counts alongside the outputs.
