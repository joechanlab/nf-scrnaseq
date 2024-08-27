import argparse
import muon as mu

parser = argparse.ArgumentParser(description="Preprocessing ATAC data by muon.")
parser.add_argument("input_h5", help="Path to the input 10x h5ad format file.")
parser.add_argument("output_h5ad", help="The path to the output 10x h5ad file")
args = parser.parse_args()

mdata = mu.read_10x_h5(args.input_h5)
mdata.var_names_make_unique()

rna = mdata.mod["rna"]
rna.write_h5ad(args.output_h5ad)
