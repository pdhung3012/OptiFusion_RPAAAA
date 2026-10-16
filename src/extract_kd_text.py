
from docx import Document
from pathlib import Path
import argparse


def extract_docx_text(input_file, output_file):
    doc = Document(input_file)
    extracted_text = []

    # Extract paragraphs and tables in document order
    for block in doc.iter_inner_content():
        if hasattr(block, "text"):
            # Paragraph
            if isinstance(block.text, str):
                if block.text.strip():
                    extracted_text.append(block.text)

        if hasattr(block, "rows"):
            # Table
            for row in block.rows:
                row_text = [
                    cell.text.strip()
                    for cell in row.cells
                ]
                extracted_text.append("\t".join(row_text))

    # Save extracted text
    text = "\n".join(extracted_text)

    Path(output_file).parent.mkdir(
        parents=True, exist_ok=True
    )

    with open(output_file, "w", encoding="utf-8") as f:
        f.write(text)

    print(f"Extracted {len(text)} characters.")
    print(f"Saved to: {output_file}")

    return text


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Extract text from DOCX."
    )

    parser.add_argument(
        "--input_file",
        default="../expert_knowledge_documentation_genral.docx",
        help="Input DOCX file"
    )

    parser.add_argument(
        "--output_file",
        default="../expert_knowledge_documentation_genral.txt",
        help="Output TXT file"
    )

    args = parser.parse_args()

    extract_docx_text(
        args.input_file,
        args.output_file
    )
