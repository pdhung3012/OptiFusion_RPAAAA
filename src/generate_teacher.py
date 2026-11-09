
import os
import json
import argparse
from pathlib import Path
from openai import OpenAI

f1=open("/home/hungphd/openai_key.txt", "r", encoding="utf-8")
api_key=f1.read().strip()
f1.close()

# input("aaaa")

client = OpenAI(api_key=api_key)


def build_prompt(sample, knowledge_content):
    return f"""You are optimizing C code using OpenMP pragma directives ONLY.

Task:
Insert or adjust suitable '#pragma omp' directives within regions delimited by
'#pragma scope' and '#pragma endscop' (preserve these marker lines unchanged).
Do not modify, reorder, remove, or introduce any C statements, declarations,
expressions, loops, variable names, or function signatures. Do not apply
non-pragma source transformations such as tiling, fusion, interchange,
skewing, or unrolling.

IMPORTANT OPENMP OPTIMIZATION STRATEGIES (examine ALL categories):
1. Parallelism optimization:
   - Consider '#pragma omp parallel for', '#pragma omp for', and
     'collapse(n)' for independent nested loops with valid canonical forms.
   - Avoid oversubscription and inappropriate nested parallel regions.
2. Scheduling optimization:
   - Consider 'schedule(static)', 'schedule(dynamic[, chunk])',
     'schedule(guided[, chunk])', or 'schedule(runtime)' based on
     iteration cost, locality, and workload imbalance.
3. Synchronization optimization:
   - Consider 'nowait' when no cross-loop dependency or required barrier exists.
   - Use 'atomic' or 'critical' only when needed for correctness; do not
     remove required synchronization or introduce races.
4. Data-sharing optimization:
   - Use 'private', 'firstprivate', 'lastprivate', 'shared', and
     'default(none)' where valid to express correct variable sharing.
   - Check variable lifetime, initialization, and per-thread state.
5. Reduction optimization:
   - Consider 'reduction(+:x)', 'reduction(*:x)', 'reduction(max:x)',
     'reduction(min:x)', or other appropriate reductions for recognized
     independent accumulation patterns.
   - Do not add a reduction when it is incompatible with existing C code
     or would require changing statements; consider numerical effects.
6. SIMD vectorization:
   - Consider '#pragma omp simd' or '#pragma omp parallel for simd',
     along with 'simdlen(n)', 'safelen(n)', or 'aligned(...)' ONLY when
     loop dependencies and alignment assumptions justify them.
7. Thread management optimization:
   - Consider 'num_threads(n)' and 'proc_bind(close/spread/primary)'
     only when workload characteristics and hardware assumptions justify
     explicit settings; prefer portability over arbitrary constants.

Selection and correctness rules:
- Evaluate every category, but apply ONLY those that are demonstrably safe,
  applicable, and likely beneficial; do NOT force all seven categories.
- Prefer appropriate combined constructs rather than redundant nested pragmas.
- Never parallelize loops with loop-carried dependencies unless the directive
  preserves the original behavior. Avoid data races, deadlocks, incorrect
  barriers, illegal data-sharing attributes, and unsupported directive placement.
- Use syntactically valid OpenMP directives and clauses compatible with the
  target compiler. If a proposed optimization requires C-code changes,
  do not perform it; choose another pragma-only optimization or leave it as is.
- Preserve the original program's semantics and all C code verbatim.
- Keep all edits inside the marked scope regions. If no safe pragma-level
  optimization is applicable, return the original code unchanged.

Original code:
{sample["original_code"]}

Expert Knowledge (use when consistent with the constraints above):
{knowledge_content}

Output requirements:
Return ONLY the complete resulting C code, with no explanations, comments
about your choices, Markdown fences, or surrounding text.
"""

import os
from openai import OpenAI

# Set your OpenAI API key

# os.environ["OPENAI_API_KEY"] = ""


def generate_teacher_data(input_file, output_file, model,kd_documentation):
    Path(output_file).parent.mkdir(
        parents=True, exist_ok=True
    )

    f1=open(kd_documentation, "r", encoding="utf-8")
    knowledge_content=f1.read()
    f1.close()


    with open(input_file, encoding="utf-8") as src, \
         open(output_file, "w", encoding="utf-8") as dst:

        data = json.load(src)
        for index, sample in enumerate(data):
            # sample = json.loads(line)
            prompt = build_prompt(sample,knowledge_content)

            response = client.responses.create(
                model=model,
                input=[
                    {
                        "role": "system",
                        "content": (
                            "You are an expert HPC compiler "
                            "optimization assistant."
                        ),
                    },
                    {
                        "role": "user",
                        "content": prompt,
                    },
                ],
            )

            optimized_code = response.output_text.strip()

            record = {
                "id": index,
                "code_location": sample['code_location'],
                "completion": optimized_code,
                "teacher": model,
                "verified": False,
            }

            dst.write(json.dumps(record) + "\n")
            print("Generated:", index)

knowledge_content=''

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--input", default="input_original.json"
    )
    parser.add_argument(
        "--kd_documentation", default="../expert_knowledge_documentation_genral.txt"
    )

    # parser.add_argument(
    #     "--api_key_location", default="/home/hungphd/openai_key.txt"
    # )

    parser.add_argument(
        "--output", default="output_teacher_gpt-5.2.json"
    )
    parser.add_argument(
        "--model", default="gpt-5.2"
    )

    args = parser.parse_args()

    generate_teacher_data(
        args.input, args.output, args.model,args.kd_documentation
    )
