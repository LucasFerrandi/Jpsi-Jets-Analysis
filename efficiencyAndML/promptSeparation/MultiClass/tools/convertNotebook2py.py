
import nbformat
from pathlib import Path
import argparse


parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--input")
args = parser.parse_args()
input = args.input
print(input)

src = Path(input)
dst = src.with_suffix(".py")

nb = nbformat.read(src, as_version=4)

with dst.open("w") as f:
    for cell in nb.cells:
        if cell.cell_type == "code":
            f.write(cell.source + "\n\n")

print(f"Created {dst}")