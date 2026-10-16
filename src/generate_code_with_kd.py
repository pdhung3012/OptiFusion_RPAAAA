
import json
import argparse
from pathlib import Path


def extract_code(input_json, folder_output):
    output_root = Path(folder_output).resolve()
    output_root.mkdir(parents=True, exist_ok=True)

    with open(input_json, "r", encoding="utf-8") as f:
        # data = json.load(f)
        data = [
        json.loads(line)
        for line in f
        if line.strip()
        ]

    # Support a JSON array or a single JSON object
    if isinstance(data, dict):
        data = [data]

    count = 0

    for item in data:
        code_location = item.get("code_location")
        code = item.get("completion")

        if not code_location or not isinstance(code, str):
            print(f"Skipping invalid item: {item.get('id')}")
            continue

        # Preserve relative directory structure
        output_path = (output_root / code_location).resolve()

        # Prevent paths from escaping output directory
        if not output_path.is_relative_to(output_root):
            print(f"Skipping unsafe path: {code_location}")
            continue

        # Create directories if necessary
        output_path.parent.mkdir(parents=True, exist_ok=True)

        # Write extracted C code
        with open(output_path, "w", encoding="utf-8") as f:
            f.write(code)

        print(f"Saved: {output_path}")
        count += 1

    print(f"\nSuccessfully extracted {count} code files.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Extract C code from JSON records"
    )

    parser.add_argument(
        "--input_json",
        default="output_teacher_gpt-5.2.json",
        help="Path to input JSON file"
    )

    parser.add_argument(
        "--folder_output",
        default="../generated-code/kd-context/",
        help="Output directory"
    )

    args = parser.parse_args()

    extract_code(args.input_json, args.folder_output)
