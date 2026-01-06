#!/usr/bin/env python3
"""
Cell2location Pipeline Script for Spatial Transcriptomics Analysis
==================================================================

This script performs a complete Cell2location analysis pipeline for spatial
transcriptomics data deconvolution. The workflow consists of two main steps:

1. Reference Model Training:
   - Trains a RegressionModel on single-cell RNA-seq reference data
   - Extracts cell type-specific expression signatures
   - Exports reference signatures for spatial deconvolution

2. Spatial Deconvolution:
   - Trains a Cell2location model on spatial transcriptomics data
   - Uses reference signatures to estimate cell type abundances
   - Exports posterior distributions of cell type locations

Prerequisites:
  - Python packages: scvi-tools, cell2location, scanpy, pandas
  - Input files: h5ad format (AnnData objects)
  - Reference data must contain cell type labels in obs
  - Spatial data should contain batch information in obs (if applicable)

Output:
  - Annotated spatial AnnData object with cell type abundance estimates
  - Optional: Reference signatures CSV file (if --save-reference-signatures is used)

For more information, see: https://cell2location.readthedocs.io/
"""

import argparse
import sys
import os
from pathlib import Path

import scvi
import pandas as pd
import cell2location
from cell2location.models import RegressionModel, Cell2location
from cell2location.utils.filtering import filter_genes
import torch
torch.set_default_tensor_type(torch.FloatTensor) 

def main():
    parser = argparse.ArgumentParser(
        description="""
Run Cell2location analysis pipeline for spatial transcriptomics deconvolution.

This script performs a two-stage analysis:
  1. Trains a RegressionModel on scRNA-seq reference data to extract
     cell type-specific expression signatures
  2. Uses these signatures to deconvolve spatial transcriptomics data
     and estimate cell type abundances at each spatial location

The output is an annotated AnnData object containing posterior distributions
of cell type abundances, which can be used for downstream spatial analysis.
        """,
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  Basic usage with default parameters:
    %(prog)s -r scRNA_reference.h5ad -s spatial_data.h5ad -o results.h5ad

  Custom training epochs and save reference signatures:
    %(prog)s -r ref.h5ad -s spatial.h5ad -o out.h5ad \\
        --ref-epochs 500 --spatial-epochs 50000 \\
        --save-reference-signatures signatures.csv

  Custom annotation keys:
    %(prog)s -r ref.h5ad -s spatial.h5ad -o out.h5ad \\
        --labels-key celltype --batch-key sample_id

  High-precision analysis with more samples:
    %(prog)s -r ref.h5ad -s spatial.h5ad -o out.h5ad \\
        --num-samples 2000 --spatial-epochs 50000

Notes:
  - Training time increases with the number of epochs
  - More samples (--num-samples) provide better posterior estimates but take longer
  - Ensure reference and spatial data share common genes for best results
  - The output h5ad file contains all original data plus new deconvolution results
        """
    )
    
    # Required arguments
    parser.add_argument(
        '-r', '--reference',
        required=True,
        type=str,
        metavar='FILE',
        help='Path to scRNA-seq reference h5ad file. Must contain cell type '
             'labels in obs (specified by --labels-key, default: cell_type)'
    )
    parser.add_argument(
        '-s', '--spatial',
        required=True,
        type=str,
        metavar='FILE',
        help='Path to spatial transcriptomics h5ad file. Should contain '
             'batch information in obs if using batch correction (specified '
             'by --batch-key, default: batch)'
    )
    parser.add_argument(
        '-o', '--output',
        required=True,
        type=str,
        metavar='FILE',
        help='Path to output h5ad file. Will contain original spatial data '
             'plus deconvolution results (cell type abundance estimates)'
    )
    
    # Optional arguments for data keys
    parser.add_argument(
        '--labels-key',
        type=str,
        default='cell_type',
        metavar='KEY',
        help='Column name in adata_ref.obs containing cell type labels. '
             'This should be a categorical column with cell type annotations '
             'for each cell in the reference dataset (default: cell_type)'
    )
    parser.add_argument(
        '--batch-key',
        type=str,
        default='batch',
        metavar='KEY',
        help='Column name in adata_vis.obs containing batch/sample information. '
             'Used for batch correction during spatial deconvolution. '
             'Set to None if no batch correction is needed (default: batch)'
    )
    
    # Optional arguments for training
    parser.add_argument(
        '--ref-epochs',
        type=int,
        default=250,
        metavar='N',
        help='Maximum number of training epochs for the RegressionModel. '
             'More epochs may improve signature quality but increase runtime. '
             'Typical range: 200-500 (default: 250)'
    )
    parser.add_argument(
        '--spatial-epochs',
        type=int,
        default=30000,
        metavar='N',
        help='Maximum number of training epochs for the Cell2location model. '
             'This is the main training step and typically requires many epochs. '
             'Typical range: 20000-50000 (default: 30000)'
    )
    parser.add_argument(
        '--num-samples',
        type=int,
        default=1000,
        metavar='N',
        help='Number of posterior samples to draw when exporting results. '
             'More samples provide better estimates but increase computation time. '
             'Typical range: 500-2000 (default: 1000)'
    )
    
    # Optional: save intermediate results
    parser.add_argument(
        '--save-reference-signatures',
        type=str,
        default=None,
        metavar='FILE',
        help='Optional path to save reference signatures as a CSV file. '
             'This allows reusing signatures without retraining the reference model. '
             'Useful for analyzing multiple spatial datasets with the same reference'
    )
    
    args = parser.parse_args()
    
    # Validate input files exist
    if not os.path.exists(args.reference):
        print(f"Error: Reference file not found: {args.reference}", file=sys.stderr)
        sys.exit(1)
    
    if not os.path.exists(args.spatial):
        print(f"Error: Spatial file not found: {args.spatial}", file=sys.stderr)
        sys.exit(1)
    
    # Create output directory if needed
    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    
    try:
        # Step 1: Load scRNA reference
        print(f"Loading scRNA reference from: {args.reference}")
        adata_ref = scvi.data.read_h5ad(args.reference)
        print(f"  Reference data shape: {adata_ref.shape}")
        
        # Step 2: Train reference model
        print(f"\nTraining RegressionModel (max_epochs={args.ref_epochs})...")
        RegressionModel.setup_anndata(adata=adata_ref, labels_key=args.labels_key)
        mod = RegressionModel(adata_ref)
        mod.train(max_epochs=args.ref_epochs)
        
        # Step 3: Export reference signatures
        print(f"Exporting reference signatures (num_samples={args.num_samples})...")
        adata_ref = mod.export_posterior(adata_ref, sample_kwargs={'num_samples': args.num_samples})
        inf_aver = pd.DataFrame(
            adata_ref.varm['means_per_cluster_mu_fg'],
            index=adata_ref.var_names,
            columns=adata_ref.uns['mod']['factor_names']
        )
        print(f"  Reference signatures shape: {inf_aver.shape}")
        
        # Optionally save reference signatures
        if args.save_reference_signatures:
            print(f"Saving reference signatures to: {args.save_reference_signatures}")
            inf_aver.to_csv(args.save_reference_signatures)
        
        # Step 4: Load spatial data
        print(f"\nLoading spatial data from: {args.spatial}")
        adata_vis = scvi.data.read_h5ad(args.spatial)
        print(f"  Spatial data shape: {adata_vis.shape}")
        
        # Step 5: Train cell2location model
        print(f"\nTraining Cell2location model (max_epochs={args.spatial_epochs})...")
        Cell2location.setup_anndata(adata=adata_vis, batch_key=args.batch_key)
        mod2 = Cell2location(adata_vis, cell_state_df=inf_aver)
        mod2.train(max_epochs=args.spatial_epochs)
        
        # Step 6: Export and save results
        print(f"Exporting posterior (num_samples={args.num_samples})...")
        adata_vis = mod2.export_posterior(adata_vis, sample_kwargs={'num_samples': args.num_samples})
        
        print(f"\nSaving results to: {args.output}")
        adata_vis.write(args.output)
        print("Cell2location analysis completed successfully!")
        
    except Exception as e:
        print(f"Error during analysis: {e}", file=sys.stderr)
        import traceback
        traceback.print_exc()
        sys.exit(1)


if __name__ == "__main__":
    main()
