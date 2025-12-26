#!/usr/bin/env bash
set -euo pipefail

# Optional: create/activate env first
# mamba create -n cell2loc_env python=3.9 -y
# conda activate cell2loc_env

# Install pinned deps (close to tutorial versions)
install_cell2loc_deps() {
  python - <<'PY'
import sys, subprocess
pkgs = [
    "scanpy>=1.9.1",
    "anndata>=0.8.0",
    "cell2location==0.1.3",
    "scvi-tools==0.19.0",
    "torch==1.13.1",
    "pytorch-lightning==1.7.7",
    "pyro-ppl==1.8.4",
]
subprocess.check_call([sys.executable, "-m", "pip", "install", "--quiet"] + pkgs)
PY
}

# Run cell2location using existing AnnData files for scRNA reference and spatial data
# Args:
#   1: SCRNA_H5AD (required) - scRNA reference with obs cell-type labels
#   2: SPATIAL_H5AD (required) - spatial counts (Visium) as AnnData
#   3: OUT_DIR (default: ./cell2loc_run)
#   4: CELL_LABEL_KEY (default: cell_type)
#   5: BATCH_KEY (default: batch; use "None" to disable)
#   6: EXPECTED_CELLS_PER_SPOT (default: 8)
#   7: N_TOP_GENES (default: 3000)
#   8: DEVICE (default: cpu; use cuda if available)
run_cell2loc_anndata() {
  # Defaults
  local SCRNA_H5AD=""
  local SPATIAL_H5AD=""
  local OUT_DIR="./cell2loc_run"
  local CELL_LABEL_KEY="cell_type"
  local BATCH_KEY="batch"
  local EXPECTED_CELLS_PER_SPOT=8
  local N_TOP_GENES=3000
  local DEVICE="cpu"

  # Parse key=value args
  while [[ $# -gt 0 ]]; do
    case "$1" in
      SCRNA_H5AD=*|scrna_h5ad=*) SCRNA_H5AD="${1#*=}" ;;
      SPATIAL_H5AD=*|spatial_h5ad=*) SPATIAL_H5AD="${1#*=}" ;;
      OUT_DIR=*|out_dir=*) OUT_DIR="${1#*=}" ;;
      CELL_LABEL_KEY=*|cell_label_key=*) CELL_LABEL_KEY="${1#*=}" ;;
      BATCH_KEY=*|batch_key=*) BATCH_KEY="${1#*=}" ;;
      EXPECTED_CELLS_PER_SPOT=*|expected_cells_per_spot=*) EXPECTED_CELLS_PER_SPOT="${1#*=}" ;;
      N_TOP_GENES=*|n_top_genes=*) N_TOP_GENES="${1#*=}" ;;
      DEVICE=*|device=*) DEVICE="${1#*=}" ;;
      *)
        echo "Unknown argument: $1"
        return 1
        ;;
    esac
    shift
  done

  if [[ -z "$SCRNA_H5AD" || -z "$SPATIAL_H5AD" ]]; then
    echo "Usage: run_cell2loc_anndata SCRNA_H5AD=/path/to/ref.h5ad SPATIAL_H5AD=/path/to/spatial.h5ad [OUT_DIR=...] [CELL_LABEL_KEY=...] [BATCH_KEY=...] [EXPECTED_CELLS_PER_SPOT=...] [N_TOP_GENES=...] [DEVICE=cpu|cuda]"
    return 1
  fi

  local REF_MODEL_DIR="${OUT_DIR}/ref_nb_model"
  local MAP_MODEL_DIR="${OUT_DIR}/c2l_map"
  mkdir -p "${OUT_DIR}" "${REF_MODEL_DIR}" "${MAP_MODEL_DIR}"

  SCRNA_H5AD="${SCRNA_H5AD}" \
  SPATIAL_H5AD="${SPATIAL_H5AD}" \
  OUT_DIR="${OUT_DIR}" \
  REF_MODEL_DIR="${REF_MODEL_DIR}" \
  MAP_MODEL_DIR="${MAP_MODEL_DIR}" \
  N_TOP_GENES="${N_TOP_GENES}" \
  EXPECTED_CELLS_PER_SPOT="${EXPECTED_CELLS_PER_SPOT}" \
  CELL_LABEL_KEY="${CELL_LABEL_KEY}" \
  BATCH_KEY="${BATCH_KEY}" \
  DEVICE="${DEVICE}" \
  python - <<'PY'
import scanpy as sc
import scvi
import cell2location
from cell2location.models import RegressionModel, Cell2location
import os

SCRNA_H5AD = os.environ["SCRNA_H5AD"]
SPATIAL_H5AD = os.environ["SPATIAL_H5AD"]
OUT_DIR = os.environ["OUT_DIR"]
REF_MODEL_DIR = os.environ["REF_MODEL_DIR"]
MAP_MODEL_DIR = os.environ["MAP_MODEL_DIR"]
N_TOP_GENES = int(os.environ["N_TOP_GENES"])
EXPECTED_CELLS_PER_SPOT = float(os.environ["EXPECTED_CELLS_PER_SPOT"])
CELL_LABEL_KEY = os.environ["CELL_LABEL_KEY"]
BATCH_KEY = os.environ.get("BATCH_KEY")
DEVICE = os.environ.get("DEVICE", "cuda")

print("Loading scRNA reference AnnData...")
adata_ref = sc.read_h5ad(SCRNA_H5AD)
print("Loading spatial AnnData...")
adata_vis = sc.read_h5ad(SPATIAL_H5AD)

common_genes = adata_ref.var_names.intersection(adata_vis.var_names)
adata_ref = adata_ref[:, common_genes].copy()
adata_vis = adata_vis[:, common_genes].copy()

sc.pp.highly_variable_genes(
    adata_ref,
    n_top_genes=N_TOP_GENES,
    flavor="seurat_v3",
    subset=True,
    batch_key=BATCH_KEY if BATCH_KEY != "None" else None,
)

scvi.model.setup_anndata(
    adata_ref,
    labels_key=CELL_LABEL_KEY,
    batch_key=BATCH_KEY if BATCH_KEY != "None" else None,
)
mod = RegressionModel(adata_ref)
mod.train(max_epochs=250, accelerator="gpu" if DEVICE == "cuda" else "cpu")
mod.save(REF_MODEL_DIR, overwrite=True)

adata_ref = mod.export_posterior(adata_ref, sample_kwargs={"num_samples": 1000})
sig_mat = adata_ref.varm["means_per_cluster_mu_fg"]
print("Signature matrix shape:", sig_mat.shape)

Cell2location.setup_anndata(
    adata_vis,
    batch_key=None,
    layer=None,
)

c2l = Cell2location(
    adata_vis,
    cell_state_df=sig_mat,
    N_cells_per_location=EXPECTED_CELLS_PER_SPOT,
    detection_alpha=200,
)
c2l.train(
    max_epochs=30000,
    batch_size=None,
    train_size=1.0,
    accelerator="gpu" if DEVICE == "cuda" else "cpu",
    check_val_every_n_epoch=200,
)
c2l.save(MAP_MODEL_DIR, overwrite=True)

adata_vis = c2l.export_posterior(adata_vis, sample_kwargs={"num_samples": 1000})

summaries = c2l.export_posterior(
    adata_vis,
    sample_kwargs={"num_samples": 1000, "batch_size": adata_vis.n_obs},
)
adata_vis.obs[summaries["means"].keys()] = summaries["means"]["cell_abundance_w_sf"]

adata_ref.write_h5ad(os.path.join(OUT_DIR, "ref_with_signatures.h5ad"))
adata_vis.write_h5ad(os.path.join(OUT_DIR, "spatial_cell2loc_results.h5ad"))
print("Done. Outputs in", OUT_DIR)
PY
}

