#!/usr/bin/env bash
# GROMACS MD workflow: protein-only, protein+ligand (merge), or pre-built complex.
#
# Source or run:
#   source /path/to/gromacs_md_workflow.sh
#   gromacs_md_workflow action=init mdp_dir=./mdp
#
#   # Protein only (default mode=protein)
#   gromacs_md_workflow structure=protein.pdb ff=charmm36-mar2019 workdir=./md_run
#
#   # Protein + ligand: CGenFF .itp + .gro (same coordinate frame as protein; align ligand first)
#   gromacs_md_workflow mode=complex structure=protein.pdb ligand_itp=UNL.itp ligand_gro=UNL.gro ff=charmm36-mar2019
#
#   # Optional: PubChem SDF -> MOL2 (needs obabel), then CGenFF str + cgenff_charmm2gmx.py
#   gromacs_md_workflow mode=complex structure=protein.pdb ligand_sdf=lig.sdf ligand_str=lig.str \
#     cgenff_charmm2gmx=/path/to/cgenff_charmm2gmx.py
#
#   # Pre-merged complex from CHARMM-GUI / manual merge
#   gromacs_md_workflow mode=premade topology=topol.top complex_gro=complex.gro workdir=./md_run
#
# Keys: see _gromacs_md_workflow_help

gromacs_md_workflow() {
  local _action="${GROMACS_MD_ACTION:-run}"
  local _mode="${GROMACS_MD_MODE:-protein}"
  local _gmx="${GMX:-gmx}"
  local _structure="${STRUCTURE:-protein.pdb}"
  local _force_field="${FORCE_FIELD:-}"
  local _water="${WATER:-tip3p}"
  local _box_distance="${BOX_DISTANCE:-1.0}"
  local _box_type="${BOX_TYPE:-dodecahedron}"
  local _ion_pos="${ION_POS:-NA}"
  local _ion_neg="${ION_NEG:-CL}"
  local _mdp_dir="${MDP_DIR:-./mdp}"
  local _workdir="${WORKDIR:-./md_run}"
  local _ligand_itp="${LIGAND_ITP:-}"
  local _ligand_gro="${LIGAND_GRO:-}"
  local _ligand_str="${LIGAND_STR:-}"
  local _ligand_mol2="${LIGAND_MOL2:-}"
  local _ligand_sdf="${LIGAND_SDF:-}"
  local _cgenff_charmm2gmx="${CGENFF_CHARMM2GMX:-cgenff_charmm2gmx.py}"
  local _topology="${TOPOLOGY:-}"
  local _complex_gro="${COMPLEX_GRO:-}"

  local _arg _k _v
  for _arg in "$@"; do
    if [[ "$_arg" == *=* ]]; then
      _k="${_arg%%=*}"
      _v="${_arg#*=}"
      case "$_k" in
        action) _action="$_v" ;;
        mode) _mode="$_v" ;;
        gmx) _gmx="$_v" ;;
        structure) _structure="$_v" ;;
        force_field|ff) _force_field="$_v" ;;
        water) _water="$_v" ;;
        box_distance) _box_distance="$_v" ;;
        box_type) _box_type="$_v" ;;
        ion_pos) _ion_pos="$_v" ;;
        ion_neg) _ion_neg="$_v" ;;
        mdp_dir) _mdp_dir="$_v" ;;
        workdir) _workdir="$_v" ;;
        ligand_itp) _ligand_itp="$_v" ;;
        ligand_gro) _ligand_gro="$_v" ;;
        ligand_str) _ligand_str="$_v" ;;
        ligand_mol2) _ligand_mol2="$_v" ;;
        ligand_sdf) _ligand_sdf="$_v" ;;
        cgenff_charmm2gmx) _cgenff_charmm2gmx="$_v" ;;
        topology) _topology="$_v" ;;
        complex_gro) _complex_gro="$_v" ;;
        *)
          echo "gromacs_md_workflow: unknown key: ${_k}" >&2
          return 2
          ;;
      esac
    else
      case "$_arg" in
        init|--init|-i) _action=init ;;
        help|--help|-h) _action=help ;;
        run) _action=run ;;
        *)
          echo "gromacs_md_workflow: unknown argument: ${_arg} (use key=value or init|help|run)" >&2
          return 2
          ;;
      esac
    fi
  done

  case "$_action" in
    init) _gromacs_md_init_mdps "$_mdp_dir" ;;
    help|-h) _gromacs_md_workflow_help ;;
    run)
      _gromacs_md_run_pipeline \
        "$_mode" "$_gmx" "$_structure" "$_force_field" "$_water" \
        "$_box_distance" "$_box_type" "$_ion_pos" "$_ion_neg" \
        "$_mdp_dir" "$_workdir" \
        "$_ligand_itp" "$_ligand_gro" "$_ligand_str" "$_ligand_mol2" "$_ligand_sdf" \
        "$_cgenff_charmm2gmx" "$_topology" "$_complex_gro"
      ;;
    *)
      echo "gromacs_md_workflow: unknown action: ${_action}" >&2
      return 2
      ;;
  esac
}

_gromacs_md_workflow_help() {
  cat << 'EOF'
gromacs_md_workflow — GROMACS solvate → ions → EM → NVT → NPT → MD.

Parameters are key=value. Env defaults: GMX, STRUCTURE, WORKDIR, MDP_DIR, FORCE_FIELD, WATER,
BOX_*, ION_*, GROMACS_MD_ACTION, GROMACS_MD_MODE, LIGAND_*, TOPOLOGY, COMPLEX_GRO, CGENFF_CHARMM2GMX.

  action=init|run|help
  mode=protein|complex|premade
       protein — pdb2gmx on structure= (PDB/GRO), then pipeline
       complex — pdb2gmx on protein PDB, merge ligand.itp + ligand.gro (CHARMM/CGenFF), then pipeline
       premade — use topology= + complex_gro= (merged top + coords); skip pdb2gmx

  gmx=…   structure=… (protein PDB for protein/complex)   ff= / force_field=   water=tip3p
  box_distance=…  box_type=…  ion_pos=…  ion_neg=…  mdp_dir=…  workdir=…

Complex / ligand (mode=complex):
  ligand_itp=  ligand_gro=     CGenFF GROMACS .itp and matching single-molecule .gro (aligned to protein)
  ligand_sdf=  optional; needs obabel → internal ligprep.mol2 in workdir
  ligand_mol2= use with ligand_str= for automatic conversion
  ligand_str=  CGenFF .str from paramchem.yale.edu (with ligand_mol2) → runs cgenff_charmm2gmx.py
  cgenff_charmm2gmx= path to cgenff_charmm2gmx.py (default: cgenff_charmm2gmx.py on PATH)

Premade (mode=premade):
  topology=    merged topol.top
  complex_gro= merged solute coordinates (protein + ligand in one .gro)

Notes:
  • Ligand must share the same coordinate frame as the protein (dock/superpose before merge).
  • CGenFF conversion: python3 cgenff_charmm2gmx.py expects mol2 + str; outputs basename.itp/.gro.
  • Review MDP tc-grps and position restraints for your topology.

Examples:
  gromacs_md_workflow action=init mdp_dir=./mdp
  gromacs_md_workflow ff=charmm36-mar2019 structure=protein.pdb workdir=./run1
  gromacs_md_workflow mode=complex ff=charmm36-mar2019 structure=protein.pdb ligand_itp=UNL.itp ligand_gro=UNL.gro workdir=./pl
  gromacs_md_workflow mode=premade topology=topol.top complex_gro=complex.gro workdir=./md
EOF
}

_gromacs_md_die() {
  echo "ERROR: $*" >&2
  return 1
}

# Resolve path: absolute, then cwd, then prev (for files created in workdir).
_gromacs_md_resolve_path() {
  local p="$1" prev="${2:-.}"
  if [[ "$p" == /* ]]; then
    printf '%s' "$p"
    return
  fi
  if [[ -f "$PWD/$p" ]]; then
    printf '%s/%s' "$PWD" "${p#./}"
    return
  fi
  printf '%s/%s' "$prev" "$p"
}

# Merge protein+ligand topologies and .gro files (requires python3).
_gromacs_md_merge_complex_files() {
  python3 - "$@" << 'PY'
import os, re, shutil, sys

def read_gro(path):
    with open(path) as f:
        raw = f.read().splitlines()
    if len(raw) < 3:
        sys.exit("invalid .gro: " + path)
    title = raw[0]
    n = int(raw[1].strip())
    atoms = raw[2 : 2 + n]
    box = raw[2 + n] if 2 + n < len(raw) else raw[-1]
    return title, n, atoms, box

def parse_atom_line(line):
    if len(line) < 44:
        sys.exit("unexpected .gro atom line length: " + repr(line[:80]))
    resnr = int(line[0:5])
    resname = line[5:10]
    aname = line[10:15]
    atnr = int(line[15:20])
    rest = line[20:]
    return resnr, resname, aname, atnr, rest

def format_atom(resnr, resname, aname, atnr, rest):
    return f"{(resnr % 100000):5d}{resname}{aname}{(atnr % 1000000):5d}{rest}"

def merge_gro(p_path, l_path, out_path):
    _, np, pat, pbox = read_gro(p_path)
    _, nl, lat, _ = read_gro(l_path)
    max_res = max(int(x[0:5]) for x in pat)
    max_at = max(int(x[15:20]) for x in pat)
    out_lines = ["protein + ligand (merged)", str(np + nl)]
    out_lines.extend(pat)
    for L in lat:
        r, rn, an, _old, rest = parse_atom_line(L)
        max_at += 1
        new_r = max_res + r
        out_lines.append(format_atom(new_r, rn, an, max_at, rest))
    out_lines.append(pbox)
    with open(out_path, "w") as f:
        f.write("\n".join(out_lines) + "\n")

def molname_from_itp(path):
    with open(path) as f:
        for line in f:
            line = line.strip()
            if line.startswith(";") or not line:
                continue
            if line.startswith("[ moleculetype ]"):
                for line in f:
                    line = line.strip()
                    if line.startswith(";"):
                        continue
                    parts = line.split()
                    if parts:
                        return parts[0]
                break
    sys.exit("could not read [ moleculetype ] name from " + path)

def merge_top(top_path, lig_itp_src, workdir, molname):
    dest = os.path.join(workdir, "ligand.itp")
    shutil.copy2(lig_itp_src, dest)
    with open(top_path) as f:
        lines = f.read().splitlines()
    out = []
    inc_done = False
    for line in lines:
        if line.strip() == "[ system ]" and not inc_done:
            out.append('#include "ligand.itp"')
            inc_done = True
        out.append(line)
    lines2 = out
    out = []
    i = 0
    while i < len(lines2):
        line = lines2[i]
        out.append(line)
        if line.strip() == "[ molecules ]":
            i += 1
            while i < len(lines2):
                ln = lines2[i]
                if ln.strip().startswith("["):
                    out.append(f"{molname:<20} 1")
                    out.append(ln)
                    i += 1
                    break
                out.append(ln)
                i += 1
            else:
                out.append(f"{molname:<20} 1")
            continue
        i += 1
    with open(top_path, "w") as f:
        f.write("\n".join(out) + "\n")

def main():
    workdir, top_path, lig_itp, lig_gro, prot_gro = sys.argv[1:6]
    molname = molname_from_itp(lig_itp)
    merge_top(top_path, lig_itp, workdir, molname)
    merge_gro(prot_gro, lig_gro, os.path.join(workdir, "processed.gro"))

if __name__ == "__main__":
    main()
PY
}

_gromacs_md_init_mdps() {
  local mdp_dir="${1:?mdp_dir required}"
  mkdir -p "$mdp_dir"

  cat > "$mdp_dir/ions.mdp" << 'EOF'
; Used only to create ions.tpr — not a dynamics run
integrator  = steep
nsteps      = 0
EOF

  cat > "$mdp_dir/em.mdp" << 'EOF'
; Energy minimization
integrator  = steep
emtol       = 1000.0
emstep      = 0.01
nsteps      = 50000
nstlist     = 10
cutoff-scheme = Verlet
ns_type     = grid
coulombtype = PME
rcoulomb    = 1.0
rvdw        = 1.0
pbc         = xyz
EOF

  cat > "$mdp_dir/nvt.mdp" << 'EOF'
; NVT equilibration (heat with position restraints on protein)
define                  = -DPOSRES
integrator              = md
nsteps                  = 50000
dt                      = 0.002
nstxout-compressed      = 5000
continuation            = no
constraint_algorithm    = lincs
constraints             = h-bonds
cutoff-scheme           = Verlet
ns_type                 = grid
coulombtype             = PME
rcoulomb                = 1.0
rvdw                    = 1.0
DispCorr                = EnerPres
tcoupl                  = V-rescale
tc-grps                 = System
tau_t                   = 0.1
ref_t                   = 300
pcoupl                  = no
pbc                     = xyz
gen_vel                 = yes
gen_temp                = 300
gen_seed                = -1
EOF

  cat > "$mdp_dir/npt.mdp" << 'EOF'
; NPT equilibration (Berendsen barostat — common for equilibration only)
define                  = -DPOSRES
integrator              = md
nsteps                  = 50000
dt                      = 0.002
nstxout-compressed      = 5000
continuation            = yes
constraint_algorithm    = lincs
constraints             = h-bonds
cutoff-scheme           = Verlet
ns_type                 = grid
coulombtype             = PME
rcoulomb                = 1.0
rvdw                    = 1.0
DispCorr                = EnerPres
tcoupl                  = V-rescale
tc-grps                 = System
tau_t                   = 0.1
ref_t                   = 300
pcoupl                  = Berendsen
pcoupltype              = isotropic
tau_p                   = 2.0
ref_p                   = 1.0
compressibility         = 4.5e-5
refcoord_scaling        = com
pbc                     = xyz
gen_vel                 = no
EOF

  cat > "$mdp_dir/md.mdp" << 'EOF'
; Production MD (Parrinello–Rahman + V-rescale)
integrator              = md
nsteps                  = 5000000
dt                      = 0.002
nstxout-compressed      = 5000
nstenergy               = 5000
continuation            = yes
constraint_algorithm    = lincs
constraints             = h-bonds
cutoff-scheme           = Verlet
ns_type                 = grid
coulombtype             = PME
rcoulomb                = 1.0
rvdw                    = 1.0
DispCorr                = EnerPres
tcoupl                  = V-rescale
tc-grps                 = System
tau_t                   = 0.1
ref_t                   = 300
pcoupl                  = Parrinello-Rahman
pcoupltype              = isotropic
tau_p                   = 5.0
ref_p                   = 1.0
compressibility         = 4.5e-5
pbc                     = xyz
gen_vel                 = no
EOF

  echo "Wrote default MDP files to $mdp_dir"
  echo "Review and tune: time steps, temperatures, pressure, nsteps, tc-grps for your topology."
}

# _gromacs_md_run_solvation_md: cwd must be WORKDIR; topol.top + boxed.gro path chain handled inside
_gromacs_md_run_solvation_md() {
  local GMX="$1" MDP_DIR="$2" BOX_DISTANCE="$3" BOX_TYPE="$4" ION_POS="$5" ION_NEG="$6"

  echo "=== editconf: simulation box ==="
  "$GMX" editconf -f processed.gro -o boxed.gro -c -d "$BOX_DISTANCE" -bt "$BOX_TYPE"

  echo "=== solvate ==="
  "$GMX" solvate -cp boxed.gro -cs spc216.gro -o solvated.gro -p topol.top

  echo "=== genion: neutralize ==="
  "$GMX" grompp -f "$MDP_DIR/ions.mdp" -c solvated.gro -p topol.top -o ions.tpr -maxwarn 2
  echo SOL | "$GMX" genion -s ions.tpr -o solv_ions.gro -p topol.top -pname "$ION_POS" -nname "$ION_NEG" -neutral

  echo "=== energy minimization ==="
  "$GMX" grompp -f "$MDP_DIR/em.mdp" -c solv_ions.gro -p topol.top -o em.tpr -maxwarn 2
  "$GMX" mdrun -v -deffnm em

  echo "=== NVT equilibration ==="
  "$GMX" grompp -f "$MDP_DIR/nvt.mdp" -c em.gro -r em.gro -p topol.top -o nvt.tpr -maxwarn 2
  "$GMX" mdrun -deffnm nvt

  echo "=== NPT equilibration ==="
  "$GMX" grompp -f "$MDP_DIR/npt.mdp" -c nvt.gro -t nvt.cpt -r nvt.gro -p topol.top -o npt.tpr -maxwarn 2
  "$GMX" mdrun -deffnm npt

  echo "=== production MD ==="
  "$GMX" grompp -f "$MDP_DIR/md.mdp" -c npt.gro -t npt.cpt -p topol.top -o md.tpr -maxwarn 2
  "$GMX" mdrun -deffnm md
}

# mode gmx structure ff water box* ion* mdp workdir ligand* cgenff topology complex_gro
_gromacs_md_run_pipeline() {
  local MODE="$1" GMX="$2" STRUCTURE="$3" FORCE_FIELD="$4" WATER="$5"
  local BOX_DISTANCE="$6" BOX_TYPE="$7" ION_POS="$8" ION_NEG="$9"
  local MDP_DIR="${10}" WORKDIR="${11}"
  local LIGAND_ITP="${12}" LIGAND_GRO="${13}" LIGAND_STR="${14}" LIGAND_MOL2="${15}" LIGAND_SDF="${16}"
  local CGENFF_SCRIPT="${17}" TOPOLOGY="${18}" COMPLEX_GRO="${19}"

  command -v python3 &>/dev/null || { _gromacs_md_die "python3 required for complex mode merge"; return 1; }
  command -v "$GMX" &>/dev/null || { _gromacs_md_die "GROMACS not found: $GMX"; return 1; }

  local prev="$PWD"
  local mdp_base struct_path
  if [[ "$MDP_DIR" == /* ]]; then
    mdp_base="$MDP_DIR"
  else
    mdp_base="$prev/$MDP_DIR"
  fi

  local f
  for f in ions.mdp em.mdp nvt.mdp npt.mdp md.mdp; do
    [[ -f "$mdp_base/$f" ]] || {
      _gromacs_md_die "Missing $mdp_base/$f — run: gromacs_md_workflow action=init mdp_dir=$MDP_DIR"
      return 1
    }
  done

  if [[ "$STRUCTURE" == /* ]]; then
    struct_path="$STRUCTURE"
  else
    struct_path="$prev/$STRUCTURE"
  fi

  mkdir -p "$WORKDIR" || return 1
  cd "$WORKDIR" || return 1

  case "$MODE" in
    protein)
      [[ -f "$struct_path" ]] || { _gromacs_md_die "Structure not found: $struct_path"; return 1; }
      echo "=== pdb2gmx: topology (protein) ==="
      if [[ -n "$FORCE_FIELD" ]]; then
        "$GMX" pdb2gmx -f "$struct_path" -o processed.gro -water "$WATER" -ff "$FORCE_FIELD" -ignh
      else
        "$GMX" pdb2gmx -f "$struct_path" -o processed.gro -water "$WATER" -ignh
      fi
      ;;
    premade)
      local tp cg
      if [[ "$TOPOLOGY" == /* ]]; then tp="$TOPOLOGY"; else tp="$prev/$TOPOLOGY"; fi
      if [[ "$COMPLEX_GRO" == /* ]]; then cg="$COMPLEX_GRO"; else cg="$prev/$COMPLEX_GRO"; fi
      [[ -f "$tp" ]] || { _gromacs_md_die "topology= file not found: $tp"; return 1; }
      [[ -f "$cg" ]] || { _gromacs_md_die "complex_gro= file not found: $cg"; return 1; }
      echo "=== premade: copy topology and coordinates ==="
      cp -f "$tp" topol.top
      cp -f "$cg" processed.gro
      ;;
    complex)
      [[ -f "$struct_path" ]] || { _gromacs_md_die "Protein structure not found: $struct_path"; return 1; }

      if [[ -n "$LIGAND_SDF" ]]; then
        local sdf_path
        if [[ "$LIGAND_SDF" == /* ]]; then sdf_path="$LIGAND_SDF"; else sdf_path="$prev/$LIGAND_SDF"; fi
        [[ -f "$sdf_path" ]] || { _gromacs_md_die "ligand_sdf not found: $sdf_path"; return 1; }
        command -v obabel &>/dev/null || { _gromacs_md_die "ligand_sdf= requires Open Babel (obabel) for SDF→MOL2"; return 1; }
        echo "=== obabel: SDF → MOL2 ==="
        obabel "$sdf_path" -O ligand_from_sdf.mol2 -h --gen3d
        LIGAND_MOL2="ligand_from_sdf.mol2"
      fi

      if [[ -n "$LIGAND_STR" && -n "$LIGAND_MOL2" ]]; then
        local strp m2p cgenp
        strp=$(_gromacs_md_resolve_path "$LIGAND_STR" "$prev")
        m2p=$(_gromacs_md_resolve_path "$LIGAND_MOL2" "$prev")
        [[ -f "$strp" && -f "$m2p" ]] || { _gromacs_md_die "ligand_str / ligand_mol2 not found ($strp) ($m2p)"; return 1; }
        cgenp="$CGENFF_SCRIPT"
        [[ "$cgenp" == /* ]] || { [[ -f "$prev/$cgenp" ]] && cgenp="$prev/$cgenp"; }
        [[ -f "$cgenp" ]] || cgenp=$(command -v "${CGENFF_SCRIPT##*/}" 2>/dev/null || true)
        [[ -f "$cgenp" ]] || cgenp=$(command -v "$CGENFF_SCRIPT" 2>/dev/null || true)
        [[ -n "$cgenp" && -f "$cgenp" ]] || { _gromacs_md_die "cgenff_charmm2gmx not found: $CGENFF_SCRIPT"; return 1; }
        echo "=== cgenff_charmm2gmx: ligand .mol2 + .str → .itp + .gro (order matches Justin Lemkul tutorial) ==="
        python3 "$cgenp" "$m2p" "$strp" || return 1
        local base
        base=$(basename "$m2p" .mol2)
        if [[ -f "${base}.itp" ]]; then LIGAND_ITP="${base}.itp"; else
          _gromacs_md_die "expected ${base}.itp after cgenff_charmm2gmx; check script output names"
          return 1
        fi
        if [[ -f "${base}.gro" ]]; then LIGAND_GRO="${base}.gro"
        elif [[ -f "${base}_fix.gro" ]]; then LIGAND_GRO="${base}_fix.gro"
        else
          _gromacs_md_die "expected ${base}.gro (or ${base}_fix.gro) after cgenff_charmm2gmx"
          return 1
        fi
      fi

      local lip lgp
      lip=$(_gromacs_md_resolve_path "$LIGAND_ITP" "$prev")
      lgp=$(_gromacs_md_resolve_path "$LIGAND_GRO" "$prev")
      [[ -f "$lip" && -f "$lgp" ]] || {
        _gromacs_md_die "mode=complex needs ligand_itp= and ligand_gro= (or ligand_str=+ligand_mol2= for conversion)"
        return 1
      }

      echo "=== pdb2gmx: protein topology ==="
      if [[ -n "$FORCE_FIELD" ]]; then
        "$GMX" pdb2gmx -f "$struct_path" -o protein_only.gro -water "$WATER" -ff "$FORCE_FIELD" -ignh -p topol.top
      else
        "$GMX" pdb2gmx -f "$struct_path" -o protein_only.gro -water "$WATER" -ignh -p topol.top
      fi

      echo "=== merge: ligand.itp into topol.top + combine .gro → processed.gro ==="
      _gromacs_md_merge_complex_files "$PWD" topol.top "$lip" "$lgp" protein_only.gro || return 1
      ;;
    *)
      _gromacs_md_die "unknown mode: $MODE (use protein, complex, or premade)"
      return 1
      ;;
  esac

  _gromacs_md_run_solvation_md "$GMX" "$mdp_base" "$BOX_DISTANCE" "$BOX_TYPE" "$ION_POS" "$ION_NEG" || return 1

  cd "$prev" || true
  echo "Done. Main trajectory: $WORKDIR/md.xtc  Checkpoint: $WORKDIR/md.cpt"
}

if [[ "${BASH_SOURCE[0]}" == "${0}" ]]; then
  set -euo pipefail
  gromacs_md_workflow "$@"
fi
