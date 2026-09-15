# Review of `docs/permutation_analysis_plan.md`

- **Reviewed:** the plan at `rp/permutation_null_analysis` commit `67f64667` (cut from `05de7d4f`, which already includes the PR #325 two-breaks fix).
- **Date:** 2026-09-15.
- **Line references** point at that commit's tree: `https://github.com/broadinstitute/regional_missense_constraint/blob/67f64667/<path>#L<line>`. The working tree on `kc/mane_plus_clinical_freeze` has newer commits on `rmc/utils/constraint.py`, so local line numbers may differ by a few.
- **Hail version:** every version-dependent claim below was run on Hail 0.2.134, 0.2.137, and 0.2.139 with identical results. The behaviour in question dates from the random-number-generator rewrite in 0.2.106, so any version from 0.2.106 onward behaves the same. Everything else is repo code, GCS state, or sampling logic and does not depend on Hail.

## Summary

The plan's structural claims about the search hold up (see "What checks out"). The review found seven correctness problems, two of which would corrupt results without any error, and four efficiency problems.

**Fix before any run:**

1. The smoke test and the full run share freeze 99, and `finalize` merges every round directory it finds, so leftovers from the smoke test end up in the full run's results (C1).
2. Two empty-round failure modes are missing from the fix list: `create_no_breaks_he` reads round 1's `merged.ht` unconditionally, and `merge_hts` shells out to `gsutil ls` on a `raw_results/` directory that may not exist (C2).

**Fix during implementation:** C3 to C7 and E1 to E4.

## Correctness

### C1. Smoke test and full run share freeze 99, and stale rounds get merged

**What.** Round outputs live under `temp/single_breaks/99/round<N>/`, `temp/simul_breaks/99/round<N>/`, and `temp/freeze99/round<N>_merged_break_found.ht`, none of which are on the 4-day bucket. A full run that needs fewer rounds than the smoke test inherits the smoke test's higher-numbered rounds.

**Evidence.**
- `get_break_search_round_nums` globs every `round(\d+)/` directory under the freeze's bucket (`rmc/utils/constraint.py:1249-1264`); `check_break_search_round_nums` only requires the numbers to be consecutive from 1 and to match between single and simultaneous search (`:1305-1319`), which stale rounds satisfy.
- `merge_rmc_hts` unions the no-break table of every round from 2 onward (`:1454-1470`).
- `search-for-single-break` overwrites only the rounds the new run reaches (`rmc/pipeline/regional_constraint.py:225-248`), and `merge-single-simul` writes nothing on a round with no breaks (`:394-399`). A stale `round<N>_merged_break_found.ht` from a longer earlier run therefore survives, and any driver that ends the loop on file existence reads it as "breaks found" and keeps going on stale input.
- `freeze99_search_union.ht` on the 4-day bucket is `_read_if_exists` unless `--overwrite-temp` is passed (`:1471-1475`).

**Fix.** Give the smoke test its own freeze number (add both to `FREEZES`), or delete `temp/single_breaks/99`, `temp/simul_breaks/99`, `temp/freeze99`, `model/4.1/99`, and `constraint/4.1/99` before the full run and pass `--overwrite-temp` to `finalize`. Have `run-rounds` decide termination from the merge step's own result (return value or log line), not from file existence.

### C2. Two empty-round failure modes beyond the three in the plan

**(a) `create_no_breaks_he`.** Called by `finalize` (`regional_constraint.py:459`), it reads `round1/final_results/merged.ht` unconditionally (`constraint.py:1127-1135`) and then builds `hl.literal(simul_sections)` (`:1146`). With fix 2 in place, that file does not exist when round 1 had no two-break result, which is realistic for a 5-gene smoke test. Guard it with `file_exists` the way `merge_round_no_break_ht` does (`:1384-1393`), and skip the filter when the set is empty, because `hl.literal` of an empty set raises (C3).

**(b) `merge_simul_break_temp_hts`.** It runs `gsutil ls` on `raw_results/` through `subprocess.check_output` (`:1185-1192`) before the zero-rows check the plan targets (`:1221-1224`). When a round sends zero sections to the two-break step (every remaining section found a single break), `run_batches_dataproc` loops over zero groups (`rmc/pipeline/two_breaks/run_batches_dataproc.py:129-156`), the directory is never created, `gsutil ls` exits non-zero, and the step dies with `CalledProcessError`, not the `DataException` that `--allow-no-results` is designed around. The flag must treat a missing directory as "no results", and `merge_hts.py` must then skip the read at line 52 and the `sections.he` write at 63-67.

### C3. Fix 1's "empty set" must be typed

`hl.literal(set())` and `hl.literal(frozenset())` raise `ExpressionException: cannot impute array elements`, and so does a filter built on either. `hl.experimental.write_expression(hl.empty_set(hl.tstr), path)` reads back as an empty set, and `list()` of it gives zero groups in `run_batches_dataproc` (`:111-122`). The repo already notes this limitation at `constraint.py:1986-1989`. Verified on 0.2.134, 0.2.137, and 0.2.139 (Appendix B).

### C4. Verification assertion 3 will fail for legitimate reasons

**What.** "The `-ctrl` replicate must yield the same regions as published freeze 2" cannot hold for every gene, and a failure would be misread as a harness bug.

**Evidence.**
- Freeze 2's search ran in May 2025: `model/4.1/2/constraint_prep.ht` and `temp/single_breaks/2/round1/break_found.ht` are dated 2025-05-13, `temp/simul_breaks/2/round1/final_results/merged.ht` 2025-05-14, and `constraint/4.1/2/all_rmc.ht` 2025-06-25 (Appendix A).
- PR #325 (merged 2026-09-03, in the plan's base commit) changes ten freeze 2 sections under the fixed code, and its description calls that a lower bound because a moved round-1 breakpoint changes which sections enter round 2:

  | Round | Transcripts |
  |---|---|
  | 1 | ENST00000246186, ENST00000253693, ENST00000360835, ENST00000381668, ENST00000559447 |
  | 2 | ENST00000286835, ENST00000325324, ENST00000396946 |
  | 3 | ENST00000303924 |
  | 4 | ENST00000280772 |

**Fix.** Drop those ten transcripts from the sampling frame, or record them as expected mismatches. The `-ctrl` replicate is the correct matched real arm regardless (the plan already uses it that way); the comparison with published freeze 2 is only a fidelity check and needs the known exceptions.

Assertion 1 ("must equal exactly") also needs a tolerance: `hl.agg.sum` over floats in a different row order differs in the last bits. Compare `section_exp` with a relative tolerance such as 1e-9.

### C5. The browser globals are not the searched set

**What.** Phase 0 derives both pools and the denominator from `rmc_browser.versions[2]` globals, but `add_globals_rmc_browser` scopes every set to the constraint HT's canonical transcripts plus `extra_transcripts` (`constraint.py:2266-2285`), not to what the search actually ran on.

**Evidence** (read from GCS, Appendix A):

| Frame | RMC-positive | Denominator | Fraction |
|---|---|---|---|
| QC-pass browser globals (`transcripts.*`) | 6,361 | 17,841 | 35.7% |
| All browser globals (`all_transcripts.*`) | 6,739 | 19,375 | 34.8% |
| Searched (`no_breaks.he` + browser transcripts) | 6,739 | 20,007 | 33.7% |

`no_breaks.he` for freeze 2 holds 13,268 transcripts against 12,636 in `all_transcripts.transcripts_no_rmc`; the 632 difference is the same gap `freeze3_run_log.md` records under "Freeze 4 verification". The plan's "36%" is the QC-pass frame, while its stated frame ("all searched transcripts") gives 33.7%.

**Fix.** Build the negative pool from `no_breaks.he` and the positive pool from the browser HT's transcripts (6,739 including outliers; use `transcripts.rmc_transcripts` if outliers should be excluded), or derive both from the `constraint_prep.versions[2]` section prefixes. Use one frame for the 36/64 split and for the "real rate" comparison, and state whether outlier transcripts are in or out. Sort each pool before seeded sampling (`random.Random(seed).sample(sorted(pool), k)`); set iteration order changes from process to process.

### C6. The randomness line does not give a reproducible permutation as written

**What.** "`hl.rand_unif` with the fixed seed, or `hl.shuffle`" leaves two ways to get it wrong.

**Evidence** (Appendix B, identical on 0.2.134, 0.2.137, 0.2.139):
- `hl.shuffle(a, seed=...)` ignores `seed`; its body is `sorted(a, key=lambda _: hl.rand_unif(0.0, 1.0))`. It has been that way since the 0.2.106 RNG rewrite (0.2.105 still passed the seed) and is unchanged on Hail `main`.
- A random value is determined by the per-call seed, the position (row, array element, loop index), and the session's `global_seed`. Passing the same explicit seed at separate Python call sites therefore yields the same permutation at every site, so a Python `for p in range(100)` collapses all replicates into copies. One call site nested in `hl.range(100).map(...)` gives distinct permutations per gene and per replicate.
- `global_seed` is never set by any `hl.init` in this repo. An unset value currently behaves as 0 only because the `rng_nonce` flag defaults to `"0x0"`, which makes the random-nonce branch in `hail/context.py` dead code; Hail warns about exactly this whenever a seed is passed without `global_seed`.
- Values do not depend on partitioning: a table read with 4, 2, or 1 partitions gives the same values.

**Fix.**

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

Record both `SEED` and `global_seed` in the JSON artifact. The written freeze-99 prep table remains the real record of the permutation.

### C7. The snapper section misreads the offsets

**What.** The plan says the `-1`/`+1` at `constraint.py:2122`/`:2132` "encode a freeze-7-specific off-by-one assumption" and should not be carried over.

**Evidence.**
- Region starts are always `breakpoint + 1`, and the breakpoint is a prep locus, so a coding site (`regional_constraint.py:356`, and `:304`/`:313` for simultaneous breaks). A non-CDS start is therefore exactly one past an exon stop by construction, not by assumption.
- Region stops are breakpoints or CDS ends, so always coding. Only starts can ever need snapping.
- `_create_exon_position_ht` is keyed by exon stop (`constraint.py:2041`), so a "geometric" snap must look up `start - 1` anyway, which is the `-1` the plan says to drop.
- The helper is not reusable verbatim: it defaults `freeze` from the enclosing scope (`:2009`) and expects `exon_start`/`exon_stop` to have been annotated by the caller (`:1999-2002`).

**Fix.** Lift the helper, reuse it for starts, and drop the stop path. Because the browser's amino-acid fix and this snap produce the same coordinates, `rmc_browser` region sizes are directly comparable for the genome-wide arm.

### C8. Small inaccuracies in the plan text

- `format_rmc_browser_ht`'s `group_by("transcript")` (`:2845`) would not collapse replicates: the labelled transcript IDs are distinct. The CDS raise at `:1735` and the unconditional coverage read at `:2779` are the real blockers.
- `run_batches_dataproc.py` has no `--read-if-exists` flag at all (argparse at `:172-247`), so fix 3 is the only route to a retry.
- The `_tmp_repart.ht` landmine cannot collide with a production run: the path is named by the group's first section (`simultaneous_breaks.py:455`), which carries the `-perm`/`-ctrl` label. Harmless to keep, but it does not require serialising freeze 99 against other runs.
- "Region starts are set to breakpoint + 1 (`regional_constraint.py:361`)": the `+ 1` is at `:356`.

## Efficiency

### E1. `--group-size 100` is far more fragmented than production

Freeze 2 round 1 sent 16,037 sections to the two-break step (225 over the 5,000-site threshold, 15,812 under) in 2 + 8 groups, about 2,000 sections per under-threshold group (Appendix A). Each group pays a full read and `hl.literal` filter of the grouped table, a `count`, a `repartition(n_rows)`, a checkpoint, the search checkpoint, another `count`, a `group_by` plus checkpoint in `annotate_max_chisq_per_section`, and a final repartition and write (`simultaneous_breaks.py:395-481`). One hundred groups for the null's 10,100 sections pays that overhead roughly twenty times more often than production did. Use a group size near 2,000; fix 3 makes retrying a large group safe.

### E2. The exploded prep will be badly partitioned

A grouped table of ~100 rows explodes to ~15M rows spread over only the partitions that held those rows, and `key_by` keeps the partition count. Freeze 2's prep has 4,382 partitions. Write the freeze-99 prep with an explicit `repartition`, or pass `--n-partitions` to round 1, which the search already supports (`regional_constraint.py:157-161`). Later rounds read the previous round's merged table, whose partitioning derives from round 1.

### E3. The driver design is unspecified in a way that matters

Each of the six steps calls `hl.init` (`regional_constraint.py:125`, `prepare_transcripts.py:43` and `:65`, `run_batches_dataproc.py:38`, `merge_hts.py:26`) and `hl.copy_log` in a `finally`, and `merge_hts` shells out to `gsutil` (`constraint.py:1186`). `run-rounds` therefore cannot call the scripts' `main()` functions in one process. It has to run locally and submit each step as a blocking `hailctl dataproc submit` job, as the freeze 3 loop did, record the job IDs (the run log's gotcha), and the cluster's max-age has to cover a multi-hour round 1.

### E4. Part of Phase 3 may be unnecessary

The single-break statistic is missing unless both sides have at least `MIN_EXP_MIS` (16) expected (`constraint.py:977-980`), and the two-break statistic requires it of all three windows (`simultaneous_breaks.py:188-190`). Expected count, not base pairs, is what drives the reviewer's variance objection, and `section_exp` plus the locus count come free from the search output (`FINAL_ANNOTATIONS`, `rmc/resources/rmc.py:239-243`). Report them alongside coding base pairs. Sampled transcripts with total expected below 32 can never break in any arm and contribute 100 wasted replicates each; at least report how many of the 100 genes fall there.

## What checks out

Verified against the plan branch's code or against GCS on 2026-09-15:

- `P_VALUE = 0.001` (`rmc.py:62`), the single-test break call (`constraint.py:1032-1034`), and no multiple-testing correction anywhere in the repo.
- The `1e-09` expected floor at `constraint.py:698`; per-locus rows keyed by `locus, section` with fields `observed`, `expected` only, and `transcript` dropped at `:720` (freeze 2's prep schema confirms this).
- The dict-valued scan keyed by section (`:239`).
- All six section parse sites are positional `split("_")` (`merge_hts.py:56`; `regional_constraint.py:303-320` and `:355-363`; `constraint.py:1151`, `:1484-1486`, `:1535`), so a `-perm` label in `[0]` survives every round.
- Round N+1 reads only the previous round's merged break-found table (`regional_constraint.py:169-174`, `:379-393`).
- The round-1 cache path `constraint/4.1/<freeze>/constraint_prep.ht` (`:140-152`) differs from the resource path `model/4.1/<freeze>/constraint_prep.ht`, so it is dead, and it is absent for freeze 99.
- `split-sections` writes only non-empty size classes (`simultaneous_breaks.py:134-153`); the zero-rows `DataException` (`constraint.py:1221-1224`); the raw-results guard (`run_batches_dataproc.py:137-141`); `merge-single-simul`'s `simul_exists` branch (`regional_constraint.py:277-278`, `:328-333`); the `"under"`/`"dataproc"` merge glob (`constraint.py:1196-1200`).
- `FREEZES` at `rmc.py:29`; `constraint_prep.versions[99]` raises a bare `KeyError` today (checked by importing the module).
- `sections_to_run` is `list()` of a frozenset (`run_batches_dataproc.py:56-76`), so group membership is not stable across processes; `--group-size` defaults to one section per group (`:120-122`).
- `finalize` makes no transcript-table joins (`regional_constraint.py:403-459`), and `--filter-outliers`/`--filter-to-canonical` are opt-in (`:429-438`).
- `merge_rmc_hts` requires at least two rounds (`constraint.py:1443-1447`) and starts at round 2 (`:1454`); `rmc_browser` is keyed by transcript with a `regions` array (`:2845`), `transcripts_no_rmc` sits inside the two global structs (`:2316-2335`), and `create_rmc_release_downloads` documents the top-level read bug (`:2869-2871`).
- Section strings become GCS object names (`simultaneous_breaks.py:339`, `:455`; `constraint.py:914`), so the no-`/` rule is right.
- PR #325 (`kc/fix_two_breaks_index_gap`) is an ancestor of the plan branch.

## Appendix A. GCS state and freeze 2 counts

Listed with `gsutil` (run with `CLOUDSDK_PYTHON` pointing at a Python 3.11 env; the SDK's bundled 3.9 fails on import) and read with local Hail (`gcs_requester_pays_configuration="broad-mpg-gnomad"`).

| Object | Written |
|---|---|
| `gs://regional_missense_constraint/model/4.1/2/constraint_prep.ht` | 2025-05-13 |
| `gs://regional_missense_constraint/temp/single_breaks/2/round1/break_found.ht` | 2025-05-13 |
| `gs://regional_missense_constraint/temp/simul_breaks/2/round1/final_results/merged.ht` | 2025-05-14 |
| `gs://regional_missense_constraint/temp/freeze2/round1_merged_break_found.ht` | 2025-05-14 |
| `gs://regional_missense_constraint/constraint/4.1/2/all_rmc.ht` | 2025-06-25 |
| `gs://regional_missense_constraint/constraint/4.1/2/no_breaks.he` | 2025-06-25 |
| `gs://regional_missense_constraint/constraint/4.1/2/rmc_browser.ht` | 2026-04-22 |
| `gs://regional_missense_constraint/constraint/4.1/99/constraint_prep.ht` | absent |

Freeze 2 has round directories 1 to 16 under `temp/simul_breaks/2/`; round 16 has `prep/` and `raw_results/` but no `final_results/`. Round 1's `raw_results/` holds 10 tables: `simul_break_dataproc_0..1.ht` (over threshold) and `simul_break_dataproc_under_0..7.ht`.

| Freeze 2 quantity | Value |
|---|---|
| `constraint_prep` partitions | 4,382 |
| `constraint_prep` row fields / key | `locus, observed, expected, section` / `locus, section` |
| browser `transcripts.{all, rmc, no_rmc}` | 17,841 / 6,361 / 11,480 |
| browser `all_transcripts.{all, rmc, no_rmc, outlier}` | 19,375 / 6,739 / 12,636 / 1,534 |
| browser rows (transcripts with regions) | 6,739 |
| `no_breaks.he` | 13,268 |
| round 1 sections to two-break search (over / under) | 225 / 15,812 |

## Appendix B. Hail checks

Run on 0.2.134 (`gnomad-qc` env), 0.2.137 (`hail-0.2.137`), and 0.2.139 (`hail`), all with `hl.init(global_seed=0)` unless stated. Results were identical across versions.

```python
SEED = 42; vals = hl.range(6)
def seeded_shuffle(a, seed):
    return hl.sorted(a, key=lambda _: hl.rand_unif(0.0, 1.0, seed=seed))

# hl.shuffle ignores seed: three call sites, same seed, three different orders
[hl.eval(hl.shuffle(vals, seed=SEED)) for _ in range(3)]
# -> [[5,3,4,0,1,2], [2,0,3,4,5,1], [4,0,5,1,2,3]]

# Seed really passed: separate call sites collapse to one permutation
[hl.eval(seeded_shuffle(vals, SEED)) for _ in range(3)]
# -> [[1,3,5,2,4,0], [1,3,5,2,4,0], [1,3,5,2,4,0]]

# One call site nested in hl.range: distinct per outer index
hl.eval(hl.range(3).map(lambda p: seeded_shuffle(vals, SEED)))
# -> [[1,4,0,5,2,3], [0,2,3,4,1,5], [5,4,2,1,0,3]]

# Table of 3 rows x nested 3: nine distinct permutations from one seed
t = hl.utils.range_table(3, 1).annotate(perms=hl.range(3).map(lambda p: seeded_shuffle(vals, SEED)))
# vs. dict of separate call sites, same seed: p0 == p1 == p2 within every row

# global_seed: 0 and unset give the same values; 1 differs
hl.eval(hl.rand_unif(0, 1, seed=SEED))   # 0.4943 with global_seed=0 or unset; 0.8100 with global_seed=1

# Partitioning: read with 4 / _n_partitions=2 / repartition(1) -> identical values

# Empty sets
hl.literal(set())                         # ExpressionException: cannot impute array elements
hl.experimental.write_expression(hl.empty_set(hl.tstr), p); hl.eval(hl.experimental.read_expression(p))  # set()
```

Plan-style construction (group by section, collect `(obs, exp)`, shuffle, zip to ascending loci, explode) was also run on a 2-section toy table: each replicate had a different order and every replicate's `section_obs`/`section_exp` totals matched the original.

## Appendix C. Hail versions

- `hl.shuffle` passed its `seed` through 0.2.105 and stopped at 0.2.106, when the RNG was rewritten (changelog 0.2.106; `_seeded_func` gained `static_rng_uid`). Unchanged through `main`.
- `hail/context.py` in 0.2.134, 0.2.137, and 0.2.139 all carry the same `global_seed is None` branch and the same `"rng_nonce": ("HAIL_RNG_NONCE", "0x0")` default in `hail/backend/backend.py`.
