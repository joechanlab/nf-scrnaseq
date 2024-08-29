import argparse
import os
import scanpy as sc
import SEACells


def parse_arguments():
    parser = argparse.ArgumentParser(description="Run SEACells on a sample")
    parser.add_argument("input_h5ad", help="Path to the input h5ad file")
    parser.add_argument("output_h5ad", help="Path to the output h5ad file")
    parser.add_argument(
        "output_dir",
        help="Path to the output directory to save metacells and assignments",
    )
    parser.add_argument("--n_SEACells", default=25, type=int, help="Number of SEACells")
    parser.add_argument(
        "--build_kernel_on", default="X_pca", help="Key to use for computing metacells"
    )
    parser.add_argument(
        "--n_waypoint_eigs",
        type=int,
        default=10,
        help="Number of eigenvalues to consider when initializing metacells",
    )
    return parser.parse_args()


def main():
    args = parse_arguments()

    ad = sc.read_h5ad(args.input_h5ad)

    model = SEACells.core.SEACells(
        ad,
        build_kernel_on=args.build_kernel_on,
        n_SEACells=args.n_SEACells,
        n_waypoint_eigs=args.n_waypoint_eigs,
        convergence_epsilon=1e-5,
    )

    model.construct_kernel_matrix()
    model.initialize_archetypes()
    model.fit(min_iter=10, max_iter=200)

    SEACell_ad = SEACells.core.summarize_by_SEACell(
        ad, SEACells_label="SEACell", summarize_layer="raw"
    )

    SEACell_ad.write_h5ad(args.output_h5ad)

    model.save_assignments(args.output_dir)

    model.get_hard_assignments().to_csv(
        os.path.join(args.output_dir, "mc_assignments.csv")
    )


if __name__ == "__main__":
    main()
