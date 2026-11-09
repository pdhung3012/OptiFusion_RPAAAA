
import os
import json
import argparse
from pathlib import Path
from openai import OpenAI

f1=open("/home/hungphd/openai_key.txt", "r", encoding="utf-8")
api_key=f1.read().strip()
f1.close()
print("API Key:", api_key)
# input("aaaa")

client = OpenAI(api_key=api_key)


def build_prompt(sample,knowledge_contenrt):
    return f"""Create the OpenMP-optimized version of the following C code. \
        Please add openmp pragmas to the regions start with '#pragma scope' and end with '#pragma endscop'. \
        The optimized code should be semantically equivalent to the original code, and should avoid data races. \
        Please do not change the content of C code. Return only the optimized code without any additional \
        text or explanations.

Original code:
{sample["original_code"]}

Expert Knowledge:
{knowledge_content}

Requirements:
1. Preserve the semantics of the original program.
2. Apply suitable OpenMP optimizations.
3. Avoid data races.
4. Keep function signatures unchanged.
5. Return only the optimized C code.
Do not include Markdown code fences.
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
