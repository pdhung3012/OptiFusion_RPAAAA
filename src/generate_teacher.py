
import os
import json
import argparse
from pathlib import Path
from openai import OpenAI

client = OpenAI()


def build_prompt(sample):
    return f"""Optimize the following C code for CPU execution.

Original code:
{sample["input"]}

Expert Knowledge:
{sample["knowledge"]}

Requirements:
1. Preserve the semantics of the original program.
2. Apply suitable OpenMP optimizations.
3. Avoid data races.
4. Keep function signatures unchanged.
5. Return only the optimized C code.
Do not include Markdown code fences.
"""


def generate_teacher_data(input_file, output_file, model):
    Path(output_file).parent.mkdir(
        parents=True, exist_ok=True
    )

    with open(input_file, encoding="utf-8") as src, \
         open(output_file, "w", encoding="utf-8") as dst:

        for line in src:
            sample = json.loads(line)
            prompt = build_prompt(sample)

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
                "id": sample["id"],
                "prompt": prompt,
                "completion": optimized_code,
                "teacher": model,
                "verified": False,
            }

            dst.write(json.dumps(record) + "\n")
            print("Generated:", sample["id"])


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--input", default="data/input.jsonl"
    )
    parser.add_argument(
        "--output", default="data/teacher.jsonl"
    )
    parser.add_argument(
        "--model", default="gpt-5.2"
    )

    args = parser.parse_args()

    generate_teacher_data(
        args.input, args.output, args.model
    )
