# Teacher sweep v2 — KD-ready data

This directory holds the **training-data-ready** output of the teacher sweep.
Each cell is **one JSON file** at
`<slug>/[withkd/]<kernel>.json` with a self-contained record of the request,
full response, extracted code, and verification. The patched `.c` file sits
next to each `.json` so the PolyBench harness can rerun verification standalone.

## Why v2 exists

v1 (`teacher_outputs/` + `runs/teacher_sweep[_withkd].csv`) captured only the
final `.c` code and token counts. For the KD plan we need the **reasoning
trace** (for SFT on CoT + student-rewrite-teacher methods) and **per-token
top-K logprobs** on the output stream (for white-box KL distillation).

The audit on Nov 23 established:
- `msg.reasoning` and `msg.reasoning_details` were silently returned by
  OpenRouter and discarded by v1's extractor
- `logprobs=True, top_logprobs=K` was never requested by v1
- 12 of 360 v1 cells also hit `finish_reason=length` at the 32K cap

v2 fixes all three.

## Data model (one JSON per cell)

```jsonc
{
  "meta": {
    "teacher": "DeepSeek-R1-0528",
    "slug": "deepseek/deepseek-r1-0528",
    "kernel": "3mm",
    "kernel_path": "code/PolyBenchC/linear-algebra/kernels/3mm/3mm.c",
    "with_kd_doc": false,
    "timestamp_utc": "2026-11-23T...",
    "sweep_version": "v2-kd-ready"
  },
  "request": {
    "system_message": "You are an expert HPC compiler optimization assistant.",
    "user_prompt": "...FULL prompt as sent to the model...",
    "kd_doc_chars": 0,
    "temperature": 0.2,
    "max_tokens": 64000,
    "reasoning_config": {},
    "top_logprobs": null         // null, 5, or 20 per teacher
  },
  "response": {
    "content_raw": "...",        // msg.content pre-regex (preserves surrounding text)
    "reasoning_text": "...",     // msg.reasoning — full CoT as single string
    "reasoning_details": [...],  // msg.reasoning_details structured list
    "refusal": null,
    "finish_reason": "stop",     // or "length" if truncated
    "logprobs": [                // null if teacher doesn't support logprobs
      {
        "token": "void",
        "logprob": -0.0014,
        "bytes": [118, 111, 105, 100],
        "top_logprobs": [
          {"token": "void", "logprob": -0.001, "bytes": [...]},
          ...                    // up to K alternatives per output token
        ]
      },
      ...                        // one entry per output token
    ],
    "usage": {
      "prompt_tokens": 1732,
      "completion_tokens": 4822,
      "reasoning_tokens": 3697,
      "total_tokens": 6554
    }
  },
  "extracted": {
    "code_c_path": "teacher_outputs_v2/.../3mm.c",  // the patched .c on disk
    "output_format": "body_only",                    // or "full_file"
    "code_chars": 1456
  },
  "verification": {
    "correct": true,
    "forbidden_edited": false,
    "all_threads_ok": "yes",
    "headline_speedup": 1.15,
    "scale_speedup": 2.18,
    "final_line": "FINAL: correct=yes headline_speedup=1.15x ...",
    "rc": 0
  },
  "latency_s": 57.44,
  "cost_usd": 0.0057,
  "error": null
}
```

## Field → KD method mapping

| field | SeqKD | SFT on reasoning | Student rewrites teacher trace | White-box top-K KL |
|---|:-:|:-:|:-:|:-:|
| `response.content_raw` / `extracted.code_c_path` | ✓ | ✓ | ✓ | ✓ |
| `response.reasoning_text` | — | **required** | **required** | — |
| `response.logprobs` | — | — | — | **required** |
| `request.user_prompt` | ✓ | ✓ | ✓ | ✓ |
| `verification.correct` | **filter** | **filter** | **filter** | **filter** |
| `verification.headline_speedup` | optional (RFT) | optional (RFT) | optional (RFT) | optional (RFT) |

Reject any cell with `correct=false` or `finish_reason=length` before training.

## Per-teacher logprobs support (verified 2026-11-23)

| teacher | top_logprobs cap | notes |
|---|---|---|
| DeepSeek-R1-0528 | — | provider returns null |
| Nemotron-3-Super | — | provider returns null |
| Qwen3.8-Max-0902 | 5 | 400 error above 5 |
| GLM-5.3 | 20 | silent drop on ~half of OR providers — see Provider pinning below |
| Kimi-K2.7-Code | 20 | silent drop on some OR providers — see Provider pinning below |
| Hermes-4-405B | 20 | Nebius-only endpoint |

## White-box KL applicability (committee audit 2026-11-23)

**Only Qwen3.8-Max-0902 shares tokenizer family with the Qwen2.5-Coder-1.5B
student.** Top-K logprobs from the other 5 teachers are over FOREIGN vocabularies
(Llama-3.1-based for Hermes; separate vocabs for GLM, Kimi, DeepSeek, Nemotron),
so standard token-level KL divergence is undefined against the Qwen student.

Usable KD methods per teacher:

| teacher | SeqKD (code) | SFT-on-reasoning | Student-rewrites-trace | White-box top-K KL |
|---|:-:|:-:|:-:|:-:|
| DeepSeek-R1-0528 | ✓ | ✓ | ✓ | ✗ (no logprobs) |
| Qwen3.8-Max-0902 | ✓ | ✓ | ✓ | ✓ **only teacher where valid** |
| GLM-5.3 | ✓ | ✓ | ✓ | ✗ (tokenizer mismatch) |
| Nemotron-3-Super | ✓ | ✓ | ✓ | ✗ (no logprobs) |
| Kimi-K2.7-Code | ✓ | ✓ | ✓ | ✗ (tokenizer mismatch) |
| Hermes-4-405B | ✓ | ⚠️ see below | ⚠️ see below | ✗ (tokenizer mismatch) |

To use the GLM/Kimi/Hermes/DeepSeek/Nemotron logprobs as KL targets you need
cross-tokenizer alignment (ULD, DSKD-style methods) or self-hosted teachers
with re-tokenization — those are Phase-2 work. We continue capturing top-K
logprobs for all logprob-capable teachers for future ablation.

## Provider pinning (retry pass, 2026-11-23)

OpenRouter routes each call across many backend providers. For GLM-5.3 (48
endpoints) and Kimi-K2.7-Code (8 endpoints), many providers silently drop
`logprobs` and/or `reasoning` fields even though they claim to support them.
`require_parameters=True` alone does not catch this.

The **retry pass** uses explicit `provider.order:[...] + allow_fallbacks:False`
with a per-teacher allowlist of providers verified (2026-11-23) to return
non-null logprobs AND non-null reasoning on a probe request:

| teacher | validated provider allowlist |
|---|---|
| GLM-5.3 | Together, Modal, DigitalOcean, Novita (verified) |
| Kimi-K2.7-Code | Moonshot AI, CoreWeave, Inceptron (verified) |
| Hermes-4-405B | Nebius (sole endpoint) |

Every retried cell records `response.provider_used` so downstream training
can stratify or filter on provider identity. Original (non-retried) cells
from the initial sweep do NOT have this field and may be mixed-provider.

## Hermes-4-405B reasoning: dual provenance

Hermes-4's Nebius endpoint silently drops `message.reasoning` when given a
generic system prompt — the model instead reasons INLINE in `content`.
This was confirmed by `scripts/probe_hermes_think_mode.py`, which found that
a **Nous-style `<think>` system prompt** wakes up the native reasoning
channel (1584 chars reasoning + 402 reasoning_tokens vs. 0 with the generic
prompt).

Two categories of Hermes cells therefore coexist in the dataset. Downstream
training can filter by `response.reasoning_source`:

- `reasoning_source: "inline_heuristic_last_marker_split"` — reasoning was
  post-hoc extracted from `content_raw` via regex split at the last
  `#include` / fenced block / `kernel_` function definition. Use with
  caution for SFT-on-reasoning; may contain formatting artifacts and quoted
  code snippets.
- `reasoning_source: "native_channel_nous_prompt"` — real reasoning channel
  surfaced by the Nous-style system prompt. Use normally.
- `reasoning_source: null` — on non-Hermes teachers, reasoning came from
  `message.reasoning` as the model produced it (no post-processing).

## Hermes-4-405B post-processing (2026-11-24)

Hermes was a disk-size outlier (1.9 GB vs 2.5-145 MB for the five other
teachers) driven by two factors: (a) full-file rewrites 10K+ output tokens
long vs body-only patches ~1-2 KB, and (b) 20 top-K logprob alternatives per
output token serialized with `indent=2`. A committee consultation and
planning fork unanimously recommended keeping Hermes but normalizing it; see
`scripts/consult_committee_hermes.py` for the memo and responses.

Post-processing applied (`scripts/normalize_hermes.py`):

1. **Full-file → body-only patches.** 40 cells emitted full-file rewrites.
   Each was re-splicer-extracted (the `#pragma scop...endscop` region only)
   into the original kernel, then re-verified with the thorough-mode harness
   (SMALL+LARGE × threads {1,2,4}). Conversions were committed only when
   correctness was preserved or improved; otherwise the original full-file
   `.c` was restored. Net: 40 converted, **0 regressions**, 1 improvement
   (`jacobi-1d` with-doc: incorrect→correct). The `extracted.output_format`
   now reads `body_only_converted_from_full_file` on these cells.
2. **`response.logprobs` dropped.** Hermes top-K are over the Llama-3.1 vocab
   and unusable for white-box KL against the Qwen student (see tokenizer
   matrix above). `response._logprobs_dropped` records the reason.
3. **Compact JSON serialization.** Dropped `indent=2`, saves ~55% file size.
4. **`reasoning_text` and `content_raw` are preserved verbatim** — no
   truncation. The SFT-on-reasoning and SeqKD signals are unaffected.
5. **Idempotency marker:** `meta.hermes_normalized: true` and
   `meta.hermes_normalized_at_utc`. Re-running the script is a no-op.
6. **Backups:** `<kernel>.c.fullfile_original` next to each converted `.c`.

Net disk: 1.9 GB → 4.6 MB (99.8% reduction). Correct cells: 50 → 51/60.

## Per-method filter policy (committee recommendation)

Different KD methods have different correctness requirements:

| method | filter |
|---|---|
| SFT on code (SeqKD) | `verification.correct == True` ONLY |
| SFT on reasoning | `verification.correct == True` ONLY (wrong reasoning teaches wrong reasoning) |
| White-box top-K KL | KEEP failures (teacher distribution includes errors) |
| Rejection fine-tuning | separate pass needed (N=4 @ T=0.7 sampling) |

Reject any cell with `finish_reason=length` or `error != null` before training,
across all methods.

## Directory layout

```
teacher_outputs_v2/
  <slug_safe>/
    <kernel>.json                       # no-doc arm: full record
    <kernel>.c                          # no-doc arm: patched compilable file
    withkd/
      <kernel>.json                     # with-doc arm
      <kernel>.c
```

Thin CSV views for leaderboard queries:
- `runs/teacher_sweep_v2.csv` (no-doc arm)
- `runs/teacher_sweep_v2_withkd.csv` (with-doc arm)

Each CSV row has a `json_path` column pointing to the full per-cell JSON.

## Reproduction

```bash
uv run python scripts/teacher_sweep_v2.py              # no-doc arm
uv run python scripts/teacher_sweep_v2.py --with-kd    # with-doc arm
```

## Accepted gaps (committee consultation 2026-11-24)

Per-teacher correctness after retry + Hermes normalization:

| teacher | no-doc | with-doc |
|---|:-:|:-:|
| DeepSeek-R1-0528 | 26/30 | 26/30 |
| Nemotron-3-Super | 26/30 | 25/30 |
| Qwen3.8-Max-0902 | 26/30 | 26/30 |
| GLM-5.3 | 19/30 | 23/30 |
| Kimi-K2.7-Code | 26/30 | 25/30 |
| Hermes-4-405B | 26/30 | 25/30 |

Union across all 6 teachers × both arms covers **26/30 kernels**. The 4
universal-fail kernels (`bicg`, `cholesky`, `gramschmidt`, `seidel-2d`) are
genuine structural blockers — wavefront / triangular / loop-carried
dependencies the PolyBench harness refuses to parallelize naïvely. No open
teacher produces correct OpenMP for any of them under our verification
protocol, with or without the KD doc.

**Decision:** we accept the 26/30 ceiling. Pushing closed teachers or custom
prompts to fill those 4 holes risks reward-hacking the dataset for the
student without meaningful correctness gains. If the student pilot later
needs those kernels specifically, we re-run those cells at T=0.7 with N=4
samples (rejection fine-tuning) as a separate pass. See
`scripts/consult_committee_coverage.py` for the memo and responses.

## Known gaps

- Reasoning-token-level logprobs: never exposed by OpenRouter for any model.
  White-box KL on the reasoning stream (not just output) requires a
  self-hosted teacher.
- Multiple-samples-per-cell (for rejection fine-tuning): v2 is still N=1 per
  cell. If RFT becomes the target method, we re-run with N=4 at temperature
  0.7 — separate pass, ~4× cost.
- The 4 universal-fail kernels (see Accepted gaps above).
