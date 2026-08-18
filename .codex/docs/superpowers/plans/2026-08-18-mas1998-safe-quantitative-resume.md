# MAS1998 Safe Quantitative Resume Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Make interrupted MAS1998 quantitative runs safely resumable at the job level, report abandoned partial jobs accurately, remove unrelated scaffold tests, and preserve all completed numerical evidence.

**Architecture:** Keep each solver invocation restart-from-zero because the custom executioner has no checkpoint state. Extend the existing monitor to distinguish a live solver from an abandoned partial output, then let the shell runner use that single state classifier to skip verified completed jobs, refuse duplicate live jobs, archive stale attempts, and run the remaining queue under an exclusive directory lock. Keep numerical solvers unchanged apart from committing the existing formatting-only edit.

**Tech Stack:** Python 3 standard library, Bash, GNU `flock` and `time`, MOOSE C++17 application, MOOSE TestHarness, GoogleTest.

---

### Task 1: Detect Live and Stale Quantitative Jobs

**Files:**
- Modify: `/home/fangxiaozhong/git_repo/mas1998_benchmark/scripts/monitor_quantitative_reproduction.py`
- Modify: `/home/fangxiaozhong/git_repo/mas1998_benchmark/scripts/tests/test_monitor_quantitative_reproduction.py`

- [x] **Step 1: Add failing process-aware state tests**

Add tests that pass an explicit set of active output paths to `inspect_job()`. Require a partial job with an active output path to be `RUNNING`, the same files without an active process to be `STALE`, and a complete CSV without a successful GNU time exit status to be `INCOMPLETE`.

- [x] **Step 2: Run the focused monitor test and observe failure**

Run:

```bash
conda run -n moose --no-capture-output python -m unittest \
  scripts.tests.test_monitor_quantitative_reproduction -v
```

Expected: failure because `inspect_job()` does not accept active output paths and has no `STALE` state.

- [x] **Step 3: Implement process-aware classification**

Read `/proc/[0-9]*/cmdline`, split NUL-delimited arguments, and collect exact paths from arguments beginning with `Executioner/output_file=`. Classify jobs as:

```text
nonzero exit status                         -> FAILED
zero exit status and complete progress      -> DONE
zero exit status and incomplete progress    -> INCOMPLETE
no exit status, active output path           -> RUNNING
no exit status, artifacts but no process     -> STALE or INCOMPLETE
no artifacts                                 -> PENDING
```

Expose `--job-state NAME` so the runner can consume the same classifier without parsing the human-readable table. Include `STALE` in snapshot totals.

- [x] **Step 4: Pass monitor tests and inspect the abandoned run**

Run the focused test command and:

```bash
conda run -n moose --no-capture-output python \
  scripts/monitor_quantitative_reproduction.py \
  runs/mas1998_quantitative_2026-08-17 --with-convergence
```

Expected: the two completed `64x32` cases are `DONE`, interrupted `128x64` cases are `STALE`, and unstarted cases are `PENDING`.

### Task 2: Add Safe Job-Level Resume to the Runner

**Files:**
- Modify: `/home/fangxiaozhong/git_repo/mas1998_benchmark/scripts/run_quantitative_reproduction.sh`
- Create: `/home/fangxiaozhong/git_repo/mas1998_benchmark/scripts/tests/test_run_quantitative_reproduction.py`

- [x] **Step 1: Add failing resume dry-run tests**

Create a temporary output directory containing one complete successful production summary and one stale partial summary. Run the script with `--dry-run --resume` and assert that it prints `SKIP` for the complete job, schedules the stale job and all pending jobs, emits exactly eight commands in default nine-job mode, and does not modify existing files. Also assert that a non-resume dry-run remains nine commands.

- [x] **Step 2: Run the focused runner test and observe failure**

Run:

```bash
conda run -n moose --no-capture-output python -m unittest \
  scripts.tests.test_run_quantitative_reproduction -v
```

Expected: failure because `--resume` is not recognized.

- [x] **Step 3: Implement preflight and resume selection**

Add `--resume`. For a normal run, continue rejecting a non-empty output directory. For a resume run:

```text
DONE        -> preserve and skip
RUNNING     -> abort before changing files
PENDING     -> schedule
STALE       -> archive and schedule from t=0
FAILED      -> archive and schedule from t=0
INCOMPLETE  -> archive and schedule from t=0
```

Acquire `flock -n` on `$output_dir/.runner.lock` before archiving or launching. Refuse a non-empty resume directory that contains none of the selected job artifacts. Make `--dry-run --resume` perform selection without creating locks, directories, or archives.

- [x] **Step 4: Preserve stale attempts and invalidate derived outputs**

Move each rerun job's existing CSV, summary, log, and time report beneath:

```text
previous_attempts/YYYYmmddTHHMMSS[_N]/JOB/
```

If any job is rerun, move existing `errors_*.csv`, `reference_convergence_64x32.csv`, and `mas1998_quantitative_diagonal.png` into the same attempt root so stale derived evidence cannot be mistaken for the resumed result.

- [x] **Step 5: Pass runner tests and shell checks**

Run:

```bash
conda run -n moose --no-capture-output python -m unittest \
  scripts.tests.test_run_quantitative_reproduction -v
bash -n scripts/run_quantitative_reproduction.sh
scripts/run_quantitative_reproduction.sh --dry-run --jobs 4 --with-convergence
scripts/run_quantitative_reproduction.sh --dry-run --resume --jobs 4 \
  --with-convergence --output-dir runs/mas1998_quantitative_2026-08-17
```

Expected: tests pass; shell syntax passes; fresh mode prints eleven commands; resume mode skips two jobs and prints nine commands without modifying the run directory.

### Task 3: Remove Scaffold Files and Update Documentation

**Files:**
- Delete: `/home/fangxiaozhong/git_repo/mas1998_benchmark/unit/src/SampleTest.C`
- Delete: `/home/fangxiaozhong/git_repo/mas1998_benchmark/unit/include/place_holder`
- Delete: `/home/fangxiaozhong/git_repo/mas1998_benchmark/test/tests/kernels/simple_diffusion/`
- Modify: `/home/fangxiaozhong/git_repo/mas1998_benchmark/README.md`
- Preserve and stage: `/home/fangxiaozhong/git_repo/mas1998_benchmark/src/utils/MAS1998BenchmarkSolver.C`
- Modify: `/home/fangxiaozhong/git_repo/moose/.codex/docs/architecture.md`
- Modify: `/home/fangxiaozhong/git_repo/moose/.codex/docs/numerical_mapping.md`

- [x] **Step 1: Delete only the confirmed scaffold artifacts**

Remove the generic arithmetic GoogleTests, empty unit include placeholder, and generic simple-diffusion TestHarness directory. Keep the MAS1998 test application, format-hook installer, characteristic inputs, and ignored generated KPP library.

- [x] **Step 2: Document restart semantics precisely**

Document this recovery command:

```bash
nohup scripts/run_quantitative_reproduction.sh \
  --resume --jobs 4 --with-convergence \
  --output-dir runs/mas1998_quantitative_2026-08-17 \
  >> runs/mas1998_quantitative_2026-08-17.launch.log 2>&1 &
```

State that resume is job-level, not checkpoint-level: successful jobs are skipped, while stale/failed/incomplete jobs are archived and restart from `t=0`. Document `STALE`, locking, active-job refusal, and derived-output invalidation in both architecture and numerical mapping documents.

### Task 4: Verify and Commit

**Files:** all files above.

- [x] **Step 1: Run script verification in the established conda environment**

Run all script tests, Ruff checks, Ruff formatting checks, shell syntax, fresh/convergence dry runs, and resume dry run. Expected: all pass without changing the active output directory.

- [x] **Step 2: Rebuild and run focused MOOSE tests**

Run from the app:

```bash
conda run -n moose --no-capture-output make -j4
conda run -n moose --no-capture-output make -C unit -j4
conda run -n moose --no-capture-output ./unit/mas1998_benchmark-unit-opt \
  --gtest_filter='MAS1998*'
conda run -n moose --no-capture-output ./run_tests \
  --re='mas1998_physics|mas1998_spherical_fv_advection' -j1
```

Expected: builds succeed; all MAS1998 unit tests pass; the focused TestHarness checks pass after scaffold removal.

- [x] **Step 3: Check diffs and commit the app**

Run `git diff --check`, inspect staged names, and commit all authorized app changes, including the pre-existing `MAS1998BenchmarkSolver.C` formatting edit, with:

```bash
git commit -m 'Add safe resume for MAS1998 quantitative runs'
```

- [x] **Step 4: Commit MOOSE documentation separately**

Stage only `.codex/docs/architecture.md`, `.codex/docs/numerical_mapping.md`, and this plan. Preserve `HETerogeneous-vectorized-or-Parallel/`, then commit with:

```bash
git commit -m 'Document MAS1998 quantitative run recovery'
```

## Self-Review

- Spec coverage: job-level recovery, duplicate-run protection, stale-attempt preservation, scaffold deletion, the authorized C++ edit, verification, and both required MOOSE documents are covered.
- Scope: no checkpoint implementation or numerical solver behavior change is introduced.
- Placeholder scan: every task has explicit files, behavior, commands, and expected outcomes; no deferred implementation placeholder remains.
- Type consistency: monitor state names are identical in Python tests, shell selection, snapshots, and documentation.
