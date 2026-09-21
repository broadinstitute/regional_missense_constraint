# RMC Permutation Null Analysis

> **Branch:** `rp/permutation_null_analysis`, rebased onto `main` at `e0d8fb6` (PR #324, freezes 3-5 defined).
> **Status:** plan revised 2026-09-17 against `docs/permutation_analysis_plan_review.md`. Line references point at the current tree.
> **Review disposition** (decided by the RMC team, not by the review):
>
> | Item | Decision |
> |---|---|
> | C1 stale rounds shared across smoke/full run | Agreed. Smoke test gets freeze 98, full run freeze 99; pre-flight asserts both empty; termination read from job logs. |
> | C2 two more empty-round failure modes | Finding accepted, remedy rejected. **No fixes to the two-break flow.** The driver routes around every empty-round failure mode externally (Phase 2). |
> | C3 typed empty set | Moot: nothing in this plan writes an empty `.he`. |
> | C4 `-ctrl` cannot match published freeze 2 | Superseded. The basis is **freeze 3**, the rerun of freeze 2 under the PR #325 fix. Freeze 3 does not exist yet; the full-run `-ctrl` assertion waits for it. |
> | C5 browser globals are not the searched set | Disagreed. `no_breaks.he` includes transcripts that should not have been searched (the v4 freeze 2 search was messy), so the **browser globals define the pool** by decision. The sample is reviewed by the team before the smoke test. |
> | C6 seeding | Agreed. |
> | C7 snapper offsets | Agreed. Starts only, via the lifted helper. |
> | C8 text inaccuracies | Corrected. |
> | E1-E4 | Agreed. |

## Context

A reviewer of the RMC paper raised two statistical objections to the interpretation of MCR-level o/e ratios (recorded in `temp/RMC permutation analysis working document.md`):

1. **Variance artifact.** MCRs are necessarily no larger than transcripts, and smaller regions have higher variance in realized o/e. So MCR o/e should have heavier tails than transcript o/e *even under a null of no regional variation in constraint*, which would partly explain Figure 1d and the finding that depleted MCRs tend to be small.
2. **Winner's curse.** The break search explicitly partitions transcripts to *maximize* o/e differences between segments, biasing discovered MCRs toward o/e values more divergent from the transcript mean than the truth.

Both objections are mechanically real in this codebase. The break call is a max over every candidate breakpoint compared against a single-test threshold (`rmc/utils/constraint.py:1033`, `P_VALUE = 0.001` at `rmc/resources/rmc.py:64`), and there is **no multiple-testing correction anywhere in the repo**, neither across breakpoints within a section nor across transcripts.

The reviewer's proposed remedy: permute position labels within each transcript, rerun the pipeline, and see what MCRs and o/e ratios the null produces.

**Intended outcome:** a quantified null distribution that lets the "MCR size" supplement section state how much of the size-constraint relationship is statistical. Published RMC calls, MPC, and missense badness are **unchanged**. This is a diagnostic, not a recalibration.

## Settled decisions

| Decision | Choice |
|---|---|
| Deliverable | Diagnostic only; no re-filtering of published calls |
| Null construction | Joint row permutation: shuffle the `(observed, expected)` pair together to a random locus within each section |
| Search fidelity | Full procedure: all iterative single-break rounds **and** the simultaneous two-break search, run to natural exhaustion |
| Allocation | 100 genes × 100 permutations, plus one unpermuted control replicate per gene |
| Sampling frame | **QC-pass browser globals** of the basis freeze (`transcripts.all_transcripts`); positive pool `transcripts.rmc_transcripts`, negative pool `transcripts.transcripts_no_rmc` |
| Stratification | Proportional to the frame: 36 positive / 64 negative (freeze 2 frame is 6,361 / 17,841 = 35.7%) |
| Sample review | The 100 transcripts are written as JSON and **reviewed by the team before the smoke test starts** |
| Basis freeze | **Freeze 3** (rerun of freeze 2 under the PR #325 two-break fix). Not yet generated. Until it exists, freeze 2 stands in wherever noted below |
| Output freezes | **98** for the smoke test, **99** for the full run. Never reused, never shared |
| Production code touched | Two behaviour-preserving changes only: `FREEZES += [98, 99]` and lifting `_create_exon_position_ht` to module scope. **The two-break flow is not modified**; the driver works around it |
| Packing | All labelled sections in one prep table; one break-search round loop per output freeze |
| Randomness | `hl.init(global_seed=SEED)` plus a single `hl.rand_unif(seed=SEED)` call site nested inside `hl.range(100).map(...)` |
| Region size metrics | Coding base pairs, **plus** `section_exp` and locus count (expected count is what drives the variance objection) |
| Real comparison arm | Matched: the `-ctrl` replicates (primary). Genome-wide: freeze 3 `rmc_results` pushed through the same Phase 3 code (secondary) |
| Headline statistics | Size-stratified excess **and** a tail-mass table |
| Interpretation threshold | Not pre-committed; interpret after seeing results |

Two properties worth noting because they simplify the design:

- **Transcript-level o/e is exactly invariant** under this permutation: the joint shuffle preserves each section's `observed` and `expected` totals (up to float summation order; see verification 1). Only the regional decomposition changes, so the null isolates precisely what the reviewer asked about.
- **Chi-square ties are not a concern.** The statistic depends on *cumulative expected*, a continuous float floored at `1e-09` (`rmc/utils/constraint.py:698`), so it strictly increases per locus. Exact ties are measure-zero.

## Key structural facts this plan relies on

- The search between `prep-constraint` and `finalize` is **self-contained**: it joins only on `section` and never reads `transcript_ref`, `transcript_cds`, or the VEP context. The only external joins are in `create_constraint_prep_ht` (before) and `format_rmc_browser_ht` (after).
- `section` is the **sole carrier of transcript identity** during the search. `create_constraint_prep_ht` drops the `transcript` field at `rmc/utils/constraint.py:720`; `merge_rmc_hts` re-derives it as `section.split("_")[0]` (`:1484`).
- The scan is **dict-valued, keyed by section** (`rmc/utils/constraint.py:239`), so many sections sharing loci is already normal. That is how overlapping canonical transcripts work today. 10,100 packed sections sits below production freeze 2's ~19,000.
- **All six section parse sites are positional `split("_")[0|1|2]`** (`merge_hts.py:56`; `regional_constraint.py:303-320`, `:355-363`; `constraint.py:1151`, `:1484-1486`, `:1535`). The permutation label must therefore use a non-underscore separator.
- Round N+1 receives **only break-found sections**, split at their breakpoints (`rmc/pipeline/regional_constraint.py:379-399`). Two-break cost is `Σ nₛ²`, so splitting a section into halves gives `n²/4`. Round 1 is the peak; no extra batching is needed for later rounds.
- Region starts are `breakpoint + 1` (`regional_constraint.py:356`; `:304`/`:313` for simultaneous breaks) where the breakpoint is a prep locus, so a coding site. Region stops are breakpoints or CDS ends. **Only starts can ever be non-coding, and a non-coding start is exactly one past an exon stop by construction.**
- Every pipeline step calls `hl.init` and `hl.copy_log` in its own `finally` (`regional_constraint.py:125`, `prepare_transcripts.py:43`/`:65`, `run_batches_dataproc.py:38`, `merge_hts.py:26`), and two steps shell out to `gsutil` (`constraint.py:1186`, `:1252`). **Steps cannot be called as functions from one process**; each is its own Dataproc job.

## Phase 0: Gene sampling and team review

New code in `rmc/utils/permutation.py`; `sample-genes` subcommand in `rmc/pipeline/permutation.py`.

**Frame.** Read the globals of `rmc_browser` for the basis freeze (`rmc/resources/rmc.py:555`). Use the **QC-pass struct** `transcripts` (`rmc/utils/constraint.py:2316-2321`): `all_transcripts` is the frame, `rmc_transcripts` the positive pool, `transcripts_no_rmc` the negative pool. This struct has the same three fields in the current code and in the on-disk freeze 2 table (which predates the June 2026 globals restructure of the `all_transcripts` struct), so read it by name and do not assume the nine-field `all_transcripts` layout. `transcripts_no_rmc` is a field of the struct, not a top-level global; reading it as one raises `AttributeError`, which is the documented bug in `create_rmc_release_downloads`.

Until freeze 3 exists, sample from **freeze 2's** browser globals. The frame (QC-pass canonical transcripts) is identical between freezes 2 and 3; only the positive/negative split can move for the handful of transcripts whose two-break results change under PR #325. Handling below.

**Draw.** Sort each pool, then `random.Random(SEED).sample(sorted(pool), k)` with `k = 36` and `64`. Set iteration order is not stable across processes, so the sort is required for reproducibility.

**Review JSON**, committed to the repo as `docs/permutation/sample_freeze99.json` so it is diffable and reviewable without a Hail session. Per transcript:

- `transcript`, `stratum` (`rmc` / `no_rmc`), `hgnc_symbol` and `chrom` from `transcript_ref` (`rmc/resources/reference_data.py:76`)
- `cds_bp` from `transcript_cds` (`:141`), and `n_loci`, `total_obs`, `total_exp` from the freeze 2 `constraint_prep` section
- `n_regions_freeze2` for the positive stratum
- flags: `total_exp < 32` (cannot break in any arm: both sides of a break need `MIN_EXP_MIS = 16`, `rmc/resources/rmc.py:97`; each such gene costs 101 wasted replicates, E4); `n_loci >= 5000` (over-threshold in the two-break search; its cost is paid 101 times, so the review should see it); `in_pr325_set` (one of the ten transcripts whose freeze 2 sections PR #325 changes, listed in the review's C4)

Plus header fields: `seed`, `global_seed`, pool sizes, frame size, realized positive fraction, basis freeze used for the strata.

**Review gate.** The team approves the list before anything runs. Any transcript the team rejects is replaced by the next seeded draw from the same pool, and the replacement is recorded in the JSON.

**Freeze 3 reconciliation.** When `rmc_results.versions[3]` exists, recompute each sampled transcript's stratum against freeze 3. A transcript that changed stratum is swapped for the next seeded draw from the correct pool, recorded in the JSON. Because the full run's strata are frozen only at this point, **the full run should not start before freeze 3's `rmc_results` exists**; the smoke test can.

**Smoke-test subset.** Five genes from the approved list: three positive, two negative, spanning short and long CDS, and **none from the PR #325 set**, so that freeze 2 is a valid provisional reference for the `-ctrl` check while freeze 3 is pending.

## Phase 1: Permuted constraint prep (freeze 98 / 99)

`prep-permuted` subcommand in `rmc/pipeline/permutation.py`, run as one Dataproc job.

**Input.** `constraint_prep.versions[2]`. Freeze 3 is a rerun of freeze 2 from this same prep (the PR #325 fix is in two-break search code only), so freeze 2's prep is the right input whether or not freeze 3 writes its own copy. All upstream work (context, mu, expected) is inherited unchanged. Filter to the sampled transcripts' sections.

**Construction** (C6). Group by section, collecting loci in ascending order and the `(observed, expected)` structs in the same order, then:

```python
hl.init(..., global_seed=SEED)

def seeded_shuffle(a, seed):
    return hl.sorted(a, key=lambda _: hl.rand_unif(0.0, 1.0, seed=seed))

g = g.annotate(
    perms=hl.range(100).map(
        lambda p: hl.struct(p=p, rows=hl.zip(g.loci, seeded_shuffle(g.vals, SEED)))
    )
)
```

Explode `perms`, then `rows`. This is the one construction that gives distinct permutations per gene and per replicate from a fixed seed:

- `hl.shuffle(a, seed=...)` **ignores** `seed` (since the Hail 0.2.106 RNG rewrite; unchanged on `main`).
- Passing the same explicit seed at separate Python call sites yields the *same* permutation at every site, so a Python `for p in range(100)` collapses all 100 replicates into identical copies with no error. One call site nested in `hl.range(100).map(...)` is what makes them differ.
- `global_seed` is never set by any `hl.init` in this repo; set it explicitly so the value is recorded rather than defaulted.

**Labels.** Section becomes `{transcript}-perm{p:03d}_{start}_{stop}`, **`-perm`, never `_perm`**. A `_` here is silently parsed as the transcript boundary at every `split("_")` site and collapses all replicates onto one section key. Add one unpermuted control replicate per gene, `{transcript}-ctrl_{start}_{stop}`, with values left in place. Labels contain no `_` and no `/` (section strings become GCS object names in the two-break temp paths, `simultaneous_breaks.py:455`).

**Partitioning** (E2). A grouped table of ~100 rows explodes to ~15M rows spread over only the partitions that held those rows, and `key_by` keeps the partition count. Write with an explicit `repartition` at production density (freeze 2: 4,382 partitions for ~28M rows, about 6,400 rows per partition): about 2,400 partitions for the full run, about 10 for the smoke test. Later rounds inherit their partitioning from round 1.

Key by `("locus", "section")` to match the production schema and write to `constraint_prep.versions[F].path`. Record `SEED`, `global_seed`, Hail version, row count, and partition count in `docs/permutation/prep_freeze{F}.json`. The written prep table is the real record of the permutation.

## Phase 2: Break search, driven externally

### Why a driver, and what it must not do

[Issue 215](https://github.com/broadinstitute/rmc_production/issues/215) documents that the two-break steps assume every round has a non-empty result of each kind, which is why the break search has always needed a human watching each round (freeze 2 needed 16 rounds). Under the null, later rounds will routinely have no over-threshold sections, no two-break results, or no sections entering the two-break step at all. **Every one of the failure modes below is guaranteed to hit a run-to-exhaustion job.**

**Decision: the two-break flow is not fixed.** The driver (`rmc/pipeline/run_permutation_rounds.py`, pure Python, no Hail import) runs on a laptop or small VM, submits each step as a blocking `hailctl dataproc submit` job with the `rmc` package as `--pyfiles`, captures stdout/stderr and the job ID, inspects GCS state between steps with `gcloud storage ls`, and decides what to run next. It refuses any `--freeze` outside `{98, 99}`.

The five failure modes and how the driver routes around each:

| # | Failure mode | Where | Driver workaround |
|---|---|---|---|
| 1 | `split-sections` writes only the size classes it has sections for, so `--run-all-sections` dies reading a missing `.he` | `simultaneous_breaks.py:134-153`; `run_batches_dataproc.py:55-76` | Never use `--run-all-sections`. After step 3, check which of `sections_to_simul_over_threshold.he` / `..._under_threshold.he` exist and submit one `run_batches_dataproc.py` job per existing class with the matching flag (`:219-234`). This mirrors production, which ran the classes separately. |
| 2 | `merge_hts.py` raises `DataException: All temp tables had 0 rows` when no section had two breaks | `constraint.py:1221-1224` | Run it; treat that exact message in the captured output as the legitimate empty outcome. Any other non-zero exit halts. |
| 3 | A round cannot be retried: `run_batches_dataproc.py` raises if a raw-results table already exists, and the script has no `--read-if-exists` flag (`:163-251`) | `run_batches_dataproc.py:137-141` | On retry, delete the round's `raw_results/` with `gcloud storage rm -r` before resubmitting. Because grouping is `list()` of a frozenset (`:55-76`, unstable across processes), a partial resume is never safe anyway; the two-break step of a round always restarts from scratch. |
| 4 | `merge_simul_break_temp_hts` shells out to `gsutil ls raw_results/` before the zero-rows check, so when zero sections enter the two-break step the directory does not exist and the job dies with `CalledProcessError` | `constraint.py:1185-1192` | Skip steps 4 and 5 entirely when neither `.he` exists (case 1 already established this). `raw_results/` is then never created and `merge_hts.py` is never called. |
| 5 | `create_no_breaks_he`, called last inside `finalize`, reads `round1/final_results/merged.ht` unconditionally | `constraint.py:1132-1134`; `regional_constraint.py:459` | `rmc_results` is written before this call (`regional_constraint.py:452-456`). Tolerate a `finalize` failure only when `all_rmc.ht/_SUCCESS` exists and the output names `round1/final_results/merged.ht`. **Phase 3 never reads `no_breaks.he`**; it derives no-break replicates itself. |

Everything the driver decides on is either file existence in a freeze namespace that pre-flight guarantees is clean, or a log line the production scripts already emit.

### Round loop

For each round N, against `--freeze F --overwrite`:

1. `regional_constraint.py search-for-single-break --search-num N`. Always. Parse nothing; it writes `break_found.ht` and `no_break_found.ht` even when empty (`regional_constraint.py:225-248`).
2. `two_breaks/prepare_transcripts.py create-grouped-ht --search-num N`. Always; it writes the round's `prep/` directory even for zero rows, which `check_break_search_round_nums` needs at `finalize` (`constraint.py:1305-1319`).
3. `two_breaks/prepare_transcripts.py split-sections --search-num N`. Always. Parse the log line `Found %i sections under threshold and %i sections over threshold`.
4. Two-break search, only for classes with a `.he` (mode 1). Over-threshold: `run_batches_dataproc.py --run-sections-over-threshold --group-size 100`. Under-threshold: `--run-sections-under-threshold --group-size 2000` (E1: freeze 2 round 1 ran 15,812 under-threshold sections in 8 groups, about 2,000 per group; 100 sections per group would pay the per-group overhead at `simultaneous_breaks.py:395-481` roughly twenty times more often). Parse `Found %i transcripts or transcript sections to search`, compute the expected group count, and confirm that many `_SUCCESS` markers exist under `raw_results/` before continuing, so a crashed job is not mistaken for an empty result.
5. `two_breaks/merge_hts.py --search-num N`, only if step 4 ran (modes 2 and 4).
6. `regional_constraint.py merge-single-simul --search-num N`. Always. **Termination** is the log line `No sections in round N had breakpoints (neither in single nor in simultaneous search)` (`regional_constraint.py:394-399`), cross-checked against the absence of `round{N}_merged_break_found.ht`. Disagreement halts.

After termination: `regional_constraint.py finalize --freeze F --overwrite --overwrite-temp`, **without** `--filter-outliers` and **without** `--filter-to-canonical` (they filter against real transcript sets that labelled IDs never match, `:429-438`). `--overwrite-temp` forces the `freeze{F}_search_union.ht` checkpoint on the 4-day bucket to be rewritten (`constraint.py:1471-1475`). Handle the `create_no_breaks_he` outcome per mode 5.

### Records, failure, resume

- `docs/permutation/rounds_freeze{F}.json`: per round, per step, the job ID, wall time, exit status, the parsed counts (sections broken / not broken, under / over, groups, two-break sections), and the decision taken (skipped / ran / expected-empty). Job IDs are what let a failed step be diagnosed after the fact.
- Any unrecognised failure halts the driver, writes `state_freeze{F}.json` (round, step, job ID, output tail), and exits non-zero.
- `--resume` re-runs from the recorded step with `--overwrite`; for step 4 it first deletes `raw_results/` (mode 3). A round is never partially resumed inside step 4.
- Run the driver in `tmux` under `caffeinate -i`, or on a small GCE VM; a 16-round run spans days and an interruption costs at most one step.

### Cluster

The driver takes `--cluster`, checks it exists before every submission, and if given `--cluster-start-cmd` recreates it when missing. The two-break step is meant for an autoscaling cluster (`run_batches_dataproc.py:37`); the cluster's `--max-age` must cover the longest single step (round 1 two-break search, hours), and `--max-idle` only needs to cover the gap between back-to-back jobs.

### Pre-flight (C1)

Run before Phase 1 for freeze F, and again before round 1:

- Assert absent (or delete under an explicit `--purge`): `model/4.1/F/` (before Phase 1 only), `constraint/4.1/F/`, `temp/single_breaks/F/`, `temp/simul_breaks/F/`, `temp/freezeF/`, and `gs://gnomad-tmp-4day/rmc/freezeF_*`. `get_break_search_round_nums` globs every `round(\d+)/` directory (`constraint.py:1249-1264`) and `merge_rmc_hts` unions every round from 2 onward (`:1454-1470`), so a stale higher-numbered round from an earlier run is merged into the results silently. Separate freeze numbers plus this check make that impossible.
- Assert `constraint/4.1/F/constraint_prep.ht` is absent: it is the dead round-1 cache path (`regional_constraint.py:140-152`) that would skip `process_sections` and reuse cached chi-squares if it ever existed.
- Assert `model/4.1/F/constraint_prep.ht` exists before round 1.

## Phase 3: Region assembly with an AA-free exon snapper

`assemble-regions` subcommand in `rmc/pipeline/permutation.py`; helpers in `rmc/utils/permutation.py`. This replaces `format_rmc_browser_ht`, which cannot be used on labelled transcripts: `get_transcript_cds_intervals` raises on transcripts missing from `transcript_cds` (`constraint.py:1735`), and `rmc_coverage_stats_ht` is read unconditionally (`:2795`). (Its `group_by("transcript")` at `:2861` would be fine: labelled IDs are distinct.)

Steps:

1. Read `rmc_results.versions[F]`. Split `transcript` on `-` into `real_transcript` and `replicate` (`ctrl` or `perm###`).
2. **Lift `_create_exon_position_ht` to module scope** (`constraint.py:2005-2061`). Two small changes are needed to make it importable: `freeze` becomes an explicit parameter instead of defaulting from the enclosing scope (`:2009`), and the `exon_start`/`exon_stop` annotation the caller performs at `:1999-2002` moves inside or is documented as a precondition. `fix_region_start_stop_aas` keeps calling it; behaviour unchanged.
3. **Snap starts only** (C7). For each region, look up `start - 1` in the exon-stop-keyed table (`fix_missing_start=True`, keyed at `:2041`) joined on `real_transcript`; if it hits, the start was one past an exon stop, so replace it with `next_exon_start`; otherwise keep it. This is exactly the browser's lookup at `:2115-2126`, so snapped coordinates equal the browser's amino-acid-fixed coordinates, and `rmc_browser` region sizes are directly comparable. Stops are never snapped.
4. `coding_bp` by interval intersection against the real transcript's CDS intervals (`transcript_cds`, collected per transcript; ~10 intervals each). A region crossing two exons must not count the intron.
5. `n_loci` per region: read `merge_round_no_break_ht(N, F)` for rounds 2 onward (the same inputs `merge_rmc_hts` uses), aggregate `hl.agg.count()` per section, join on the rebuilt section string. Carry `section_obs`, `section_exp`, `section_oe = section_obs / section_exp`, `section_chisq` from `rmc_results` (E4).
6. **No-break replicates.** Gene-replicates that never broke are absent from `rmc_results` (`merge_rmc_hts` starts at round 2, `constraint.py:1454`, and requires two rounds, `:1443-1447`). Derive them as the labelled section set from the freeze-F prep minus the transcripts present in `rmc_results`. Do not read `no_breaks.he` (Phase 2, mode 5).

The genome-wide real arm (freeze 3 `rmc_results`) goes through steps 1 and 3-5 of this same code, so both arms get identical coordinate handling.

## Phase 4: Comparisons (design now, build after generation)

All four use the Phase 3 table. The matched real arm is the `-ctrl` replicates; the secondary arm is genome-wide freeze 3.

1. **o/e marginal distribution.** Real vs null MCR o/e, overlaid. Matched primary, genome-wide secondary. This is the direct Figure 1d rebuttal.
2. **o/e vs region size, headline.** Bin by `coding_bp` (deciles of the pooled null size distribution, reported with per-bin counts) and repeat binned by `section_exp`. Within each bin report **size-stratified excess**: real tail mass minus null tail mass.
3. **Null break-discovery rate.** Fraction of gene-replicates with at least one region, against the real fraction in the QC-pass frame (freeze 3 browser globals once they exist; 35.7% in freeze 2). Report separately for the positive and negative strata, which differ in power, and report how many of the 100 genes fall below `total_exp = 32` (E4).
4. **Most-extreme o/e per gene.** Distribution of each gene-replicate's minimum region o/e, real vs null. Isolates winner's-curse magnitude.

Plus a **tail-mass table**: fraction of regions with o/e below 0.2 / 0.4 / 0.6, real vs null.

## Code changes

| File | Change |
|---|---|
| `rmc/resources/rmc.py:29` | Add `98` and `99` to `FREEZES` with a prominent comment: permutation null smoke test and full run; not production freezes. Leave `CURRENT_FREEZE` alone. Without this every `versions[F]` lookup dies on a bare `KeyError`; `--freeze` has no `choices` and no path builder validates it. |
| `rmc/utils/constraint.py:2005` | Lift `_create_exon_position_ht` to module scope with an explicit `freeze` parameter. No behaviour change for `fix_region_start_stop_aas`. |
| `rmc/utils/permutation.py` | **New.** Pool derivation from browser globals, seeded sampling and review JSON, permuted-prep construction, start snapper, `coding_bp`, `n_loci`, comparison functions. |
| `rmc/pipeline/permutation.py` | **New.** Hail entry points: `sample-genes`, `prep-permuted`, `assemble-regions`, `compare`. Each runs as one Dataproc job. |
| `rmc/pipeline/run_permutation_rounds.py` | **New.** The local driver: pre-flight, round loop, `finalize`, records, resume. No Hail import. Refuses freezes outside `{98, 99}`. |
| `docs/permutation/` | **New.** `sample_freeze99.json` (reviewed and committed), `prep_freeze{F}.json`, `rounds_freeze{F}.json`. |

Nothing in `rmc/pipeline/two_breaks/`, `rmc/utils/simultaneous_breaks.py`, or the round/merge logic in `rmc/utils/constraint.py` is modified. The issue-215 fixes remain a separate, later decision.

## Landmines

- **Round-1 cache short-circuit.** `regional_constraint.py:140-152` reads `constraint/4.1/{freeze}/constraint_prep.ht` (not the resource path under `model/`) and, if present, skips `process_sections`. Dead today only because of that path mismatch. Pre-flight asserts it is absent.
- **Merge glob.** `merge_simul_break_temp_hts` selects temp tables by the substrings `"under"` or `"dataproc"` (`constraint.py:1196-1200`), which all three output naming schemes match (`run_batches_dataproc.py:130-135`). Running the over and under classes separately is fine; running `--run-all-sections` in the same round as either class flag would double count. The driver never uses `--run-all-sections`.
- **Labels** must contain no `_` and no `/`.
- **`_tmp_repart.ht`** (`simultaneous_breaks.py:455`) is named by the group's first section, which carries the `-perm`/`-ctrl` label, so it cannot collide with a production run. Harmless.
- **Smoke and full run must never share a freeze number**, and a freeze number must never be reused after a failed or abandoned run without `--purge`.

## Verification

**Smoke test first: freeze 98, the five approved smoke-test genes × 3 permutations, plus their `-ctrl` replicates.** Run Phases 1-3 end to end before the 10,100-section run on freeze 99.

Assertions:

1. **Totals preserved.** For every replicate, `section_obs` summed over regions equals the real transcript's total exactly (integers), and `section_exp` matches within a relative tolerance of `1e-9` (`hl.agg.sum` over floats in a different row order differs in the last bits). A permutation that changes either is a bug; this is the single strongest correctness check available.
2. **Replicates differ.** For every gene, no two `-perm` replicates have the same `(locus, observed, expected)` row sequence, and none equals `-ctrl`. This directly catches the C6 collapse, which produces correct totals and no error.
3. **Labels survive the search.** After `merge-single-simul` in round 2, every section string still carries its `-perm###`/`-ctrl` label, and `split("_")` yields exactly three parts with `[1]` and `[2]` numeric.
4. **Unpermuted control reproduces the basis freeze.** Smoke test: `-ctrl` regions must equal published freeze 2 regions for the five genes (valid because the PR #325 set is excluded). Full run: `-ctrl` regions must equal freeze 3 regions for all 100 genes. `section_exp` compared with the same `1e-9` tolerance. If this fails, the harness is not faithful to the pipeline and nothing downstream is trustworthy.
5. **Snapper correctness.** Every snapped region start falls on a CDS position of its real transcript; every stop already does; `coding_bp <= stop - start + 1`, with equality for single-exon regions; snapped `-ctrl` coordinates equal the basis freeze's `rmc_browser` coordinates for the same regions.
6. **No silent NA.** Zero missing `start`, `stop`, `coding_bp`, `n_loci` after Phase 3.
7. **Pre-flight clean.** Both round-directory sets and the round-1 cache path absent for freeze 98 before round 1 (and again for 99).
8. **Empty-round handling.** The smoke test must run to exhaustion unattended, which will exercise the driver's decision paths: a round with only one size class, a round with no two-break results (mode 2), a round with no sections entering the two-break step (modes 1 and 4), the terminal round, and `finalize` with or without a round-1 two-break result (mode 5). Add one deliberate kill during step 4 followed by `--resume` (mode 3). If the loop cannot survive these on 5 genes, it will not survive 100.

Record the master seed, `global_seed`, the 100-gene sample, and the per-round records alongside the outputs.

## Open dependency

**Freeze 3 does not exist yet.** It blocks: the strata reconciliation that freezes the sample (so the full run start), assertion 4 for the full run, and Phase 4's real arms. It does not block: Phase 0 sampling and review (from freeze 2 globals), the smoke test on freeze 98, or the driver and Phase 1-3 code.
