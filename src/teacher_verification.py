
import json
import argparse
import subprocess
import tempfile
from pathlib import Path


def run_command(command, timeout=60):
    try:
        result = subprocess.run(
            command,
            capture_output=True,
            text=True,
            timeout=timeout,
        )
        return result.returncode == 0, result.stderr
    except (subprocess.TimeoutExpired, OSError) as error:
        return False, str(error)


def main(args):
    Path(args.output).parent.mkdir(
        parents=True, exist_ok=True
    )

    accepted = 0

    with open(args.input, encoding="utf-8") as src, \
         open(args.output, "w", encoding="utf-8") as dst:

        for line in src:
            sample = json.loads(line)

            with tempfile.TemporaryDirectory() as folder:
                source = Path(folder) / "optimized.c"
                source.write_text(
                    sample["completion"], encoding="utf-8"
                )

                compiled, error = run_command([
                    "gcc", "-std=c11", "-fopenmp",
                    "-fsyntax-only", str(source)
                ])

                correct = False
                if compiled and args.checker:
                    correct, error = run_command([
                        args.checker,
                        sample["id"],
                        str(source)
                    ])

                if compiled and correct:
                    sample["verified"] = True
                    dst.write(json.dumps(sample) + "\n")
                    accepted += 1
                else:
                    print(
                        f"Rejected {sample['id']}: "
                        f"{error[:200]}"
                    )

    print("Verified samples:", accepted)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--input", default="data/teacher.jsonl"
    )
    parser.add_argument(
        "--output", default="data/verified.jsonl"
    )
    parser.add_argument(
        "--checker",
        help=(
            "Executable correctness checker accepting "
            "sample ID and generated C file path"
        ),
    )
    main(parser.parse_args())
