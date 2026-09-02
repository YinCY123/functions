import argparse
import os
import loompy
import pandas as pd
import scipy.sparse as sparse


def convert_csv_to_loom(csv_path, loom_path, transpose=False):
    print(f"🔄 Step 1: Reading CSV file from {csv_path}...")

    try:
        # Read CSV file
        df = pd.read_csv(csv_path, index_col=0)
    except Exception as e:
        print(f"❌ Error reading CSV: {e}")
        return

    # Handle transposition if cells are rows
    if transpose:
        print("↪️ Transposing matrix (Cells x Genes -> Genes x Cells)...")
        df = df.T

    print(
        f"📊 Matrix Dimensions: {df.shape[0]} Genes and {df.shape[1]} Cells."
    )

    print("⚡ Step 2: Extracting attributes and compressing matrix...")
    gene_names = df.index.astype(str).tolist()
    cell_ids = df.columns.astype(str).tolist()

    # Convert to sparse matrix for memory efficiency
    sparse_matrix = sparse.csr_matrix(df.values)

    # Standard naming conventions for pySCENIC
    row_attrs = {"Gene": gene_names}
    col_attrs = {"CellID": cell_ids}

    print(f"💾 Step 3: Creating Loom file at {loom_path}...")
    try:
        os.makedirs(os.path.dirname(loom_path) or ".", exist_ok=True)
        loompy.create(loom_path, sparse_matrix, row_attrs, col_attrs)
        print("✅ Success! Loom file created successfully.")
    except Exception as e:
        print(f"❌ Error creating Loom file: {e}")


if __name__ == "__main__":
    # Set up command line argument parsing
    parser = argparse.ArgumentParser(
        description="Convert single-cell CSV expression matrix to pySCENIC-compatible Loom file."
    )

    # Positional arguments (required)
    parser.add_argument(
        "-i",
        "--input",
        required=True,
        help="Path to the input CSV expression file.",
    )
    parser.add_argument(
        "-o",
        "--output",
        required=True,
        help="Path where the output Loom file will be saved.",
    )

    # Optional flag
    parser.add_argument(
        "--transpose",
        action="store_true",
        help="Add this flag if your CSV has cells as rows and genes as columns.",
    )

    args = parser.parse_args()

    # Execute conversion
    convert_csv_to_loom(args.input, args.output, transpose=args.transpose)
