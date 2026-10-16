from pathlib import Path
from docx import Document
from docxcompose.composer import Composer
import argparse


def merge_docx(input_folder, output_file):
    input_folder = Path(input_folder)
    output_file = Path(output_file).resolve()

    # Get all DOCX files, sorted alphabetically
    docx_files = sorted(
        [
            f for f in input_folder.glob("*.docx")
            if f.resolve() != output_file
            and not f.name.startswith("~$")
        ]
    )

    if not docx_files:
        print("No DOCX files found.")
        return

    # Use the first document as the master document
    master = Document(docx_files[0])
    composer = Composer(master)

    print(f"Adding: {docx_files[0].name}")

    # Append remaining documents
    for file in docx_files[1:]:
        print(f"Adding: {file.name}")

        # Add page break between documents
        master.add_page_break()

        doc = Document(file)
        composer.append(doc)

    # Save merged document
    output_file.parent.mkdir(parents=True, exist_ok=True)
    composer.save(str(output_file))

    print(f"\nMerged {len(docx_files)} documents.")
    print(f"Output: {output_file}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Merge DOCX files into one document."
    )

    parser.add_argument(
        "--input_folder",
        default="../expert_knowledge_documentation/",
        help="Folder containing DOCX files"
    )

    parser.add_argument(
        "--output_file",
        default="../expert_knowledge_documentation_genral.docx",
        help="Output merged DOCX file"
    )

    args = parser.parse_args()

    merge_docx(args.input_folder, args.output_file)