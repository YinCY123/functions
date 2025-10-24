library(magrittr)
library(fs)
library(gluedocking)

# my ubuntu
# prepare_for_gluedocking(
#   python_path = "/usr/bin/python",
#   prepare_receptor_script = "/home/yincy/git/bior/functions/molecular_docking/prepare_receptor4.py",
#   prepare_ligand_script = "/home/yincy/git/bior/functions/molecular_docking/prepare_ligand4.py",
#   prepare_split_alt_script = "/home/yincy/git/bior/functions/molecular_docking/prepare_pdb_split_alt_confs.py",
#   obabel_path = "/usr/bin/obabel",
#   vina_path = "/usr/bin/vina"
# )

# CNS server
prepare_for_gluedocking(
  python_path = "/data1/yincy/miniconda3/envs/biopy/bin/python",
  prepare_receptor_script = "/home/yincy/BioHome/github/functions/molecular_docking/prepare_receptor4.py",
  prepare_ligand_script = "/home/yincy/BioHome/github/functions/molecular_docking/prepare_ligand4.py",
  prepare_split_alt_script = "/home/yincy/BioHome/github/functions/molecular_docking/prepare_pdb_split_alt_confs.py",
  obabel_path = "/data1/yincy/tools/openbabel/bin/obabel",
  vina_path = "/data1/yincy/github/AutoDock-Vina-1.2.7/build/linux/release/vina", 
  force = TRUE
)


pdb_ids <- c("1iep", "4hg7")
cids <- c("2244", "5090")

paths <- c("receptors/raw", 
    "receptors/trimmed",
    "receptors/split_alt",
    "receptors/prepared",

    "ligands/raw", 
    "ligands/mol2", 
    "ligands/prepared", 
    
    "configs", 
    
    "docked")

dir_create(paths)

pdb_files <- download_receptor(pdb_ids = pdb_ids, 
  output_dir = "receptors/raw")


ligands_files <- download_ligand(
  cids = cids, 
  output_dir = "ligands/raw", 
  verify_ssl = TRUE
)

# step 1: remove non-protien atoms (water, ligands, etc)
trimmed_files <- trim_receptor(input_paths = "receptors/raw", 
  output_dir = "receptors/trimmed")

# step 2: split alternative conformations (if present)
pdb_files <- split_alt(
  inputs = trimmed_files, 
  output_dir = "receptors/split_alt", 
  keep_label = "A" # only keep the "A" conformation
)

# step 3: convert PDB files to PDBQT format for docking
receptor_pdbqt <- prepare_receptor(
  pdb_files = pdb_files, 
  output_dir = "receptors/prepared"
)


# step 1: convert SDF files to MOL2 format
converted_files <- convert_molecule(
  input_file = "ligands/raw", 
  output_dir = "ligands/mol2", 
  output_format = "mol2"
)

# step 2: convert ligand files to PDBQT format for docking
ligand_pdbqt <- prepare_ligand(
  mol_files = converted_files, 
  output_dir = "ligands/prepared"
)

box_params <- calculate_box(
  pdb_files = receptor_pdbqt, 
  padding = 5
)

config_files <- write_configs(
  receptor_paths = receptor_pdbqt, 
  ligand_paths = ligand_pdbqt, 
  box_df = box_params, 
  output_dir = "configs", 
  exhaustiveness = 8
)

# run docking
source("../../../bior/functions/molecular_docking/run_vina.R")
args(run_vina)

results <- run_vina(
  config_paths = "configs", 
  output_dir = "docked", 
  cpu = 8
)
