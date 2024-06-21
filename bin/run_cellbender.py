#!/usr/bin/env python
import os
import argparse
import subprocess
from cellbender.remove_background.downstream import load_anndata_from_input

parser = argparse.ArgumentParser(
    description="wrapper for DoubletDetection for doublet detection from transcriptomic data.")
parser.add_argument('--raw_h5', help='Path to raw 10x h5 format.')
parser.add_argument("--filtered_h5", help = "Path to filtered 10x h5 format.")
parser.add_argument("--output_h5", help = "Path to output h5 format.")

args = parser.parse_args()

# estimate of number of cells from cellranger
adata_filtered = load_anndata_from_input(args.filtered_h5)
expected_cells = adata_filtered.shape[0]

# number of surely empty droplet
total_droplets_included = 2000

command = f"""
cellbender remove-background \
      --cuda \
      --input {args.raw_h5} \
      --output {args.output_h5} \
      --expected-cells {expected_cells} \
      --total-droplets-included {total_droplets_included}
"""

print(command)
try:
    result = subprocess.run(command, shell=True, check=True, text=True, capture_output=True)
    print("Output:", result.stdout)
except subprocess.CalledProcessError as e:
    # Print the error output
    print("Error:", e.stderr)
    # Print the command that caused the error
    print("Failed Command:", e.cmd)
    # Print the return code
    print("Return Code:", e.returncode)
