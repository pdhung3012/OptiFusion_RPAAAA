
import os
import json
import argparse

def extract_c_files(input_folder, output_file):
    dataset = []

    for root, dirs, files in os.walk(input_folder):
        for filename in sorted(files):
            if not filename.endswith(".c"):
                continue

            file_path = os.path.join(root, filename)
            relative_path = os.path.relpath(file_path, input_folder)

            with open(file_path, "r", encoding="utf-8", errors="replace") as f:
                code = f.read()

            dataset.append({
                "original_code": code,
                "code_location": relative_path
            })

    dataset.sort(key=lambda x: x["code_location"])

    with open(output_file, "w", encoding="utf-8") as f:
        json.dump(dataset, f, indent=2, ensure_ascii=False)

    print(f"Extracted {len(dataset)} C files.")
    print(f"Output saved to: {output_file}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--input_folder", default="../original-code/")
    parser.add_argument("--output_file", default="input_original.json")
    args = parser.parse_args()

    extract_c_files(args.input_folder, args.output_file)
