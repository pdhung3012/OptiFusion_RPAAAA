# Open frontier LLM teachers — A/B contribution

This directory extends OptiFusion's teacher coverage with **6 open frontier
LLMs** accessed via OpenRouter, run on **all 30 PolyBench/C kernels** (the 25
in your closed-teacher set plus the 5 you'd scoped out — `correlation`,
`covariance`, `deriche`, `floyd-warshall`, `nussinov`).

Each teacher was run under **two prompt conditions**, giving a clean A/B on
the expert-knowledge document from `expert_knowledge_documentation_genral.txt`:

- `open-<slug>/` — prompt *without* the expert-knowledge doc
- `open-withkd-<slug>/` — prompt *with* the expert-knowledge doc (same template
  position as your closed `kd-context-*` runs)

All generated code was verified end-to-end with the same thorough-mode harness
used for Phase 0:
compile → forbidden-region hash check → bit-identical output check on
SMALL+LARGE datasets × threads {1, 2, 4}. Only files that pass all gates get
a headline speedup recorded in that directory's `summary.csv`.

## Teachers

| directory | OpenRouter slug | reasoning config |
|---|---|---|
| `open-deepseek-r1-0528/` | `deepseek/deepseek-r1-0528` | default |
| `open-qwen3.8-max-0902/` | `qwen/qwen3.8-max-0902` | `{enabled, effort: medium}` |
| `open-glm-5.3/` | `z-ai/glm-5.3` | `{enabled, effort: medium}` |
| `open-nemotron-3-super-120b-a12b/` | `nvidia/nemotron-3-super-120b-a12b` | `{effort: medium}` |
| `open-kimi-k2.7-code/` | `moonshotai/kimi-k2.7-code` | default |
| `open-hermes-4-405b/` | `nousresearch/hermes-4-405b` | `{enabled}` |

Add `open-withkd-` prefix for the KD-doc arm. Same model slugs, same
reasoning config, same temperature (0.2), same `max_tokens` (32000).

## Headline findings

### 1. Teacher-dependent effect of the expert-knowledge doc

The doc **hurts** strong teachers (closed frontier, Qwen-tier) and **rescues**
weaker open teachers. Not a universal win or loss.

| teacher | no-doc median | with-doc median | Δmedian |
|---|---:|---:|---:|
| **Hermes-4-405B** | 0.92× | **1.51×** | **+65%** |
| **GLM-5.3** | 0.84× | **1.36×** | **+61%** |
| DeepSeek-R1-0528 | 1.00× | 1.14× | +14% |
| Qwen3.8-Max-0902 | 1.23× | 1.23× | +0% |
| Kimi-K2.7-Code | 1.06× | 1.03× | −3% |
| Nemotron-3-Super | 1.00× | 0.93× | −7% |

For comparison, OptiFusion's own closed-teacher data in this repo:

| teacher | no-doc median | with-doc median | Δmedian |
|---|---:|---:|---:|
| gpt-5.2 | 1.64× | 1.12× | −32% |
| gemini | 1.48× | 1.07× | −28% |

### 2. Oracle ensemble matches gpt-5.2 without a frontier closed model

Pick the best correct teacher per kernel ("oracle best-of-6"):

| source | n | median | geomean | max |
|---|---:|---:|---:|---:|
| closed gpt-5.2 alone (no-doc) | 25 | 1.64× | 0.95× | 5.1× |
| open no-doc oracle (best-of-6) | 26 | 1.44× | 1.15× | 15.1× |
| open with-doc oracle (best-of-6) | 26 | **2.00×** | 1.44× | 15.3× |
| open union oracle (best across both conditions) | 26 | **2.08×** | 1.53× | 15.3× |

The with-doc open oracle (2.00×) **beats gpt-5.2's median (1.64×) by 22%**, and
the union oracle climbs another 4% on top. The big single-kernel wins are also
~3× larger: open teachers hit 15× on `2mm` and `covariance`, closed caps at
5.1×.

### 3. Per-teacher win spread — ensembles matter more than any single teacher

In the open no-doc sweep, 6 teachers took the per-kernel win on 26 kernels. No
teacher won more than 7. This is why the oracle gains so much over any
individual:

| teacher | wins (no-doc) |
|---|---:|
| Qwen3.8-Max-0902 | 7 |
| Nemotron-3-Super | 6 |
| DeepSeek-R1-0528 | 4 |
| Kimi-K2.7-Code | 4 |
| GLM-5.3 | 3 |
| Hermes-4-405B | 2 |

## Notes on scope

- **5 kernels you'd scoped out** (`datamining/correlation`,
  `datamining/covariance`, `medley/deriche`, `medley/floyd-warshall`,
  `medley/nussinov`) are now covered. On LARGE dataset, `correlation` and
  `covariance` are where pluto hits its headline 19-20× speedups, so including
  them matters for any benchmark that wants to compare LLMs against classical
  polyhedral tools.
- **Teacher-output format normalization.** Some open teachers followed the
  prompt's "Optimize only code in regions delimited by `#pragma scop` /
  `#pragma endscop`" literally and returned only the scope body. We normalize
  with a scope-patcher (`polyBench/scripts/patch_scope.py`) that splices the
  returned body into the original PolyBench kernel file so verification runs on
  a complete .c file. The field `output_format` in each `summary.csv` records
  whether the teacher returned `full_file` or `body_only`.
- **Max tokens.** Reasoning teachers (Kimi, DeepSeek-R1, Hermes) can produce
  >10K reasoning tokens before emitting any code. We used `max_tokens=32000`
  per call; both no-doc and with-doc sweeps finished at 180/180 with 0 errors
  after a 10-minute wall-clock retry cutoff.

## Reproduction

All scripts live in the `polyBench/scripts/` tree next to this repo:

- `patch_scope.py` — body-only → full-file normalizer
- `teacher_sweep.py` — the main 6 × 30 harness driver
  (pass `--with-kd` for the doc condition)
- `retry_kimi_failures.py` — hard wall-clock cutoff wrapper for tail-latency
  failures
- `stage_for_optifusion.py` / `stage_withkd_for_optifusion.py` — reorganize
  flat teacher outputs into this repo's nested PolyBench layout

Raw sweep CSVs: `polyBench/runs/teacher_sweep.csv` (no-doc) and
`polyBench/runs/teacher_sweep_withkd.csv` (with-doc).

**Total API spend: $10.86** across both sweeps (360 calls).
