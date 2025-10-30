#!/usr/bin/env python
import argparse
import subprocess
from cellbender.remove_background.downstream import anndata_from_h5

parser = argparse.ArgumentParser(
    description="wrapper for DoubletDetection for doublet detection from transcriptomic data."
)
parser.add_argument("raw_h5", help="Path to raw 10x h5 format.")
parser.add_argument("output_h5", help="Path to output h5 format.")
parser.add_argument(
    "--total_droplets_included",
    nargs="?",
    default=None,
    type=int,
    help="Estimate of total number of droplets (optional).",
)
parser.add_argument(
    "--empty_drop_training_fraction",
    nargs="?",
    default=None,
    type=float,
    help="Fraction of empty droplets to use for training (optional).",
)
parser.add_argument("--filtered", nargs="*", help="Path to filtered cellranger output.")


args = parser.parse_args()

# estimate of number of cells from cellranger
expected_cells = 0
for idx, filtered_file_path in enumerate(args.filtered, start=1):
    adata_filtered = anndata_from_h5(filtered_file_path)
    expected_cells += adata_filtered.shape[0]
del adata_filtered

# Build command with conditional parameters
command_parts = [
    "cellbender remove-background",
    "--cuda",
    f"--input {args.raw_h5}",
    f"--output {args.output_h5}",
    f"--expected-cells {expected_cells}",
]

if args.total_droplets_included is not None:
    command_parts.extend(
        [
            f"--total-droplets-included {args.total_droplets_included}",
        ]
    )
if args.empty_drop_training_fraction is not None:
    command_parts.extend(
        [
            f"--empty-drop-training-fraction {args.empty_drop_training_fraction}",
        ]
    )

command = " \\\n      ".join(command_parts)

print(command)
try:
    result = subprocess.run(
        command, shell=True, check=True, text=True, capture_output=True
    )
    print("Output:", result.stdout)
except subprocess.CalledProcessError as e:
    # Print the error output
    print("Error:", e.stderr)
    # Print the command that caused the error
    print("Failed Command:", e.cmd)
    # Print the return code
    print("Return Code:", e.returncode)
