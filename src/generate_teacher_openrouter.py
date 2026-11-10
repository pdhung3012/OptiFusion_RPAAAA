"""Generate OpenMP-optimized teacher examples using a priority-ordered prompt."""

import argparse
import json
import os
from pathlib import Path
from openai import OpenAI


def build_prompt(sample, knowledge_content):
    return f"""You are an expert HPC compiler optimizing C code for CPU execution
with OpenMP. Prioritize loop transformations that improve locality and
parallelism, then apply suitable OpenMP directives and clauses.

Task and scope:
- Optimize only code in regions delimited by '#pragma scope' and
  '#pragma endscop'. Preserve these marker lines unchanged.
- Preserve function signatures, program inputs/outputs, and original semantics.
- You MAY restructure loops and adjust associated C statements within the
  marked regions when needed for the priority loop transformations below.
- Do not modify code outside the marked regions. Do not add dependencies,
  data races, invalid OpenMP directives, or unsupported transformations.

HIGHEST PRIORITY — LOOP TRANSFORMATIONS (examine FIRST, in this order):
1. Loop tiling / blocking:
   - First investigate tiling of suitable nested loops to improve cache
     locality, data reuse, and parallel granularity.
   - Choose legal tile sizes based on loop bounds and available knowledge;
     handle boundary tiles and partial iterations correctly.
   - Use an OpenMP tile construct only if the specified target supports it;
     otherwise implement legal C loop blocking with OpenMP parallelization.
2. Loop fusion:
   - Next investigate merging adjacent loops with compatible iteration
     spaces to improve locality and reduce overhead.
   - Fuse ONLY if there is no dependence, ordering, side effect, or
     synchronization requirement that would make the fused loop incorrect.
3. Loop skewing:
   - Next investigate skewing loops with carried dependencies where an
     affine change of iteration space can expose legal parallelism or improve
     locality, potentially in combination with tiling.
   - Derive correct transformed bounds and ensure all original iterations
     execute exactly as required, respecting dependency order.

IMPORTANT: These three are loop transformations, not generally standalone
OpenMP pragmas. Loop fusion and skewing may require rewriting C loops. Do
NOT falsely claim that adding '#pragma omp' alone performs them. Prioritize
them in analysis and implementation when legal and beneficial; do not force
a transformation that is unsafe or counterproductive.

SECOND PRIORITY — OPENMP DIRECTIVE-LEVEL OPTIMIZATIONS:
4. Parallelism: 'parallel for', 'for', 'collapse(n)' on independent loops;
   avoid nested-region overhead and invalid parallel dependencies.
5. Scheduling: 'schedule(static)', 'schedule(dynamic[, chunk])',
   'schedule(guided[, chunk])', or 'schedule(runtime)' as warranted.
6. Synchronization: minimize barriers with 'nowait' only when safe;
   use 'atomic'/'critical' where correctness requires them.
7. Data sharing: use 'private', 'firstprivate', 'lastprivate', 'shared',
   and 'default(none)' correctly; check initialization and lifetime.
8. Reductions: use 'reduction(+:x)', 'reduction(*:x)',
   'reduction(max:x)', or 'reduction(min:x)' for legal reductions;
   consider floating-point reassociation and numerical requirements.
9. SIMD: use 'simd', 'parallel for simd', 'simdlen(n)', 'safelen(n)',
   and 'aligned(...)' only when dependence and alignment permit.
10. Thread management: consider 'num_threads(n)' and
    'proc_bind(close/spread/primary)' when justified; avoid arbitrary
    hardware-specific constants that harm portability.

Optimization decision process:
- FIRST examine tiling, fusion, and skewing, in that priority order. Apply
  all that are legally compatible and likely to help, before considering
  additional OpenMP pragma-level optimizations.
- THEN examine all seven OpenMP categories, selecting only beneficial,
  legal, compatible directives and clauses.
- Maintain correct dependence order, variable scope, data sharing, and
  synchronization. Never assume independent iterations without evidence.
- Do not change observable behavior. If no safe optimization is available,
  return the original code unchanged.
- Keep all transformations strictly within the marked scope regions.

Original code:
{sample["original_code"]}

Expert Knowledge (use when consistent with correctness and scope):
{knowledge_content}

Output requirements:
Return ONLY the complete optimized C code. No explanations, Markdown code
fences, or surrounding text. Preserve the '#pragma scope' and
'#pragma endscop' markers exactly as supplied.
"""


def read_samples(input_file):
    """Accept either a JSON array or newline-delimited JSON objects (JSONL)."""
    content = Path(input_file).read_text(encoding="utf-8")
    try:
        data = json.loads(content)
    except json.JSONDecodeError:
        data = [json.loads(line) for line in content.splitlines() if line.strip()]
    if isinstance(data, dict):
        data = [data]
    if not isinstance(data, list):
        raise ValueError("Input must be a JSON array or JSONL objects")
    return data


def generate_teacher_data(input_file, output_file, model, kd_documentation, api_key_file=None):
    knowledge_content = Path(kd_documentation).read_text(encoding="utf-8")

    # Prefer the supplied key file; otherwise use the OPENROUTER_API_KEY environment variable.
    if api_key_file:
        api_key = Path(api_key_file).expanduser().read_text(encoding="utf-8").strip()
    else:
        api_key = os.getenv("OPENROUTER_API_KEY", "").strip()

    if not api_key:
        raise ValueError("Provide --api_key_file or set OPENROUTER_API_KEY")

    # OpenRouter supports the OpenAI-compatible Chat Completions API.
    client = OpenAI(
        base_url="https://openrouter.ai/api/v1",
        api_key=api_key,
    )
    samples = read_samples(input_file)
    Path(output_file).parent.mkdir(parents=True, exist_ok=True)

    with open(output_file, "w", encoding="utf-8") as dst:
        for index, sample in enumerate(samples):
            prompt = build_prompt(sample, knowledge_content)
            response = client.chat.completions.create(
                model=model,
                messages=[
                    {
                        "role": "system",
                        "content": "You are an expert HPC compiler optimization assistant.",
                    },
                    {"role": "user", "content": prompt},
                ],
            )
            optimized_code = response.choices[0].message.content
            if not optimized_code or not optimized_code.strip():
                raise ValueError(f"Empty model response for sample {index}")

            record = {
                "id": index,
                "code_location": sample["code_location"],
                "completion": optimized_code.strip(),
                "teacher": model,
                "verified": False,
            }
            dst.write(json.dumps(record, ensure_ascii=False) + "\n")
            dst.flush()
            print("Generated:", index)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", default="input_original.json")
    parser.add_argument("--kd_documentation", default="../expert_knowledge_documentation_genral.txt")
    parser.add_argument("--output", default="output_teacher_gemini-4-argon.jsonl")
    parser.add_argument("--model", default="google/gemini-4-argon",
                        help="OpenRouter model ID, e.g. anthropic/claude-opus-5.5")
    parser.add_argument("--api_key_file", default='/home/hungphd/openrouter_key.txt',
                        help="Path to a plain-text OpenRouter API key (otherwise use OPENROUTER_API_KEY)")
    args = parser.parse_args()
    generate_teacher_data(args.input, args.output, args.model,
                          args.kd_documentation, args.api_key_file)
