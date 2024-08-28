import argparse
import muon as mu

parser = argparse.ArgumentParser(
    description="Extract RNA and ATAC modalities from 10x h5 file."
)
parser.add_argument("input_h5", help="Path to the input 10x h5 format file.")
parser.add_argument("output_rna_h5ad", help="Path to the output RNA h5ad file.")
parser.add_argument("output_atac_h5ad", help="Path to the output ATAC h5ad file.")
args = parser.parse_args()

mdata = mu.read_10x_h5(args.input_h5)
mdata.var_names_make_unique()

rna = mdata.mod["rna"]
atac = mdata.mod["atac"]

rna.write_h5ad(args.output_rna_h5ad)
atac.write_h5ad(args.output_atac_h5ad)

print(f"RNA data saved to: {args.output_rna_h5ad}")
print(f"ATAC data saved to: {args.output_atac_h5ad}")
