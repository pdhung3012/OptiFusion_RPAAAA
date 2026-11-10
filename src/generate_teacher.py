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
    # api_key = os.getenv("OPENAI_API_KEY")
    if api_key_file:
        # api_key = Path(api_key_file).read_text(encoding="utf-8").strip()
        f1=open(api_key_file, "r", encoding="utf-8")
        api_key=f1.read().strip()
        f1.close()
        # print('go here')
    if not api_key:
        raise ValueError("Set OPENAI_API_KEY or provide --api_key_file")
    # print(api_key)
    # input('aaaa')
    client = OpenAI(api_key=api_key)
    samples = read_samples(input_file)
    Path(output_file).parent.mkdir(parents=True, exist_ok=True)

    with open(output_file, "w", encoding="utf-8") as dst:
        for index, sample in enumerate(samples):
            prompt = build_prompt(sample, knowledge_content)
            response = client.responses.create(
                model=model,
                input=[
                    {
                        "role": "system",
                        "content": "You are an expert HPC compiler optimization assistant."
                    },
                    {"role": "user", "content": prompt},
                ],
            )
            optimized_code = response.output_text.strip()
            record = {
                "id": index,
                "code_location": sample["code_location"],
                "completion": optimized_code,
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
    parser.add_argument("--output", default="output_teacher_gpt-6-astra.json")
    parser.add_argument("--model", default="gpt-5.2")
    parser.add_argument("--api_key_file", default="/home/hungphd/openai_key.txt",
                        help="Optional path to API key file, instead of OPENAI_API_KEY")
    args = parser.parse_args()
    generate_teacher_data(args.input, args.output, args.model,
                          args.kd_documentation, args.api_key_file)
