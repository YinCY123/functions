
# Single-run convenience function: perform MD for one protein/ligand
run_gromacs() {
    # set -euo pipefail

    # Multiple inputs supported
    local -a PROTEINS=()
    local -a LIGAND_GROS=()
    local -a LIGAND_ITPS=()
    local -a PROTEIN_DIRS=()
    local -a LIGAND_DIRS=()
    local -a LIGAND_ITP_DIRS=()
    local OUTBASE=""
    local FORCEFIELD="amber99sb-ildn"
    local WATERMODEL="spce"
    local CPUS=10

    while [[ $# -gt 0 ]]; do
        case "$1" in
            --protein=*) PROTEINS+=("${1#*=}"); shift ;;
            --ligand-gro=*) LIGAND_GROS+=("${1#*=}"); shift ;;
            --ligand-itp=*) LIGAND_ITPS+=("${1#*=}"); shift ;;
            --protein-dir=*) PROTEIN_DIRS+=("${1#*=}"); shift ;;
            --ligand-dir=*) LIGAND_DIRS+=("${1#*=}"); shift ;;
            --ligand-itp-dir=*) LIGAND_ITP_DIRS+=("${1#*=}"); shift ;;
            --outdir=*) OUTBASE="${1#*=}"; shift ;;
            --forcefield=*) FORCEFIELD="${1#*=}"; shift ;;
            --watermodel=*) WATERMODEL="${1#*=}"; shift ;;
            --cpus=*) CPUS="${1#*=}"; shift ;;
            --protein|--ligand-gro|--ligand-itp|--protein-dir|--ligand-dir|--ligand-itp-dir|--outdir|--forcefield|--watermodel|--cpus)
                if [[ -n "${2:-}" && "${2:0:1}" != "-" ]]; then
                    case "$1" in
                        --protein) PROTEINS+=("$2") ;;
                        --ligand-gro) LIGAND_GROS+=("$2") ;;
                        --ligand-itp) LIGAND_ITPS+=("$2") ;;
                        --protein-dir) PROTEIN_DIRS+=("$2") ;;
                        --ligand-dir) LIGAND_DIRS+=("$2") ;;
                        --ligand-itp-dir) LIGAND_ITP_DIRS+=("$2") ;;
                        --outdir) OUTBASE="$2" ;;
                        --forcefield) FORCEFIELD="$2" ;;
                        --watermodel) WATERMODEL="$2" ;;
                        --cpus) CPUS="$2" ;;
                    esac
                    shift 2
                else
                    echo "Error: $1 requires an argument" >&2; return 1
                fi
                ;;
            -h|--help)
                echo "Usage: gmx_md_sim [--protein FILE ...] [--protein-dir DIR ...] [--ligand-gro FILE ...] [--ligand-dir DIR ...] [--ligand-itp FILE ...] [--ligand-itp-dir DIR ...] [--outdir DIR] [--forcefield NAME] [--watermodel NAME] [--cpus N]"; return 0 ;;
            *) echo "Unknown option: $1" >&2; return 1 ;;
        esac
    done

    # Expand directories
    shopt -s nullglob
    local d f
    for d in "${PROTEIN_DIRS[@]}"; do
        if [[ -d "$d" ]]; then
            for f in "$d"/*.{pdb,PDB}; do [[ -f "$f" ]] && PROTEINS+=("$f"); done
        else
            echo "Warning: protein dir not found: $d" >&2
        fi
    done
    for d in "${LIGAND_DIRS[@]}"; do
        if [[ -d "$d" ]]; then
            for f in "$d"/*.{gro,GRO}; do [[ -f "$f" ]] && LIGAND_GROS+=("$f"); done
        else
            echo "Warning: ligand dir not found: $d" >&2
        fi
    done
    for d in "${LIGAND_ITP_DIRS[@]}"; do
        if [[ -d "$d" ]]; then
            for f in "$d"/*.{itp,ITP}; do [[ -f "$f" ]] && LIGAND_ITPS+=("$f"); done
        else
            echo "Warning: ligand itp dir not found: $d" >&2
        fi
    done
    shopt -u nullglob

    if [[ ${#PROTEINS[@]} -eq 0 ]]; then echo "Error: no proteins provided" >&2; return 1; fi
    if [[ ${#LIGAND_GROS[@]} -eq 0 ]]; then echo "Error: no ligand gro files provided" >&2; return 1; fi

    # Helper: find matching itp by basename, or single-itp fallback
    find_itp_for() {
        local ligand_gro="$1"
        local base="${ligand_gro##*/}"; base="${base%.*}"
        local itp
        for itp in "${LIGAND_ITPS[@]}"; do
            local itp_base="${itp##*/}"; itp_base="${itp_base%.*}"
            if [[ "$itp_base" == "$base" ]]; then printf '%s' "$itp"; return 0; fi
        done
        if [[ ${#LIGAND_ITPS[@]} -eq 1 ]]; then printf '%s' "${LIGAND_ITPS[0]}"; return 0; fi
        return 1
    }

    local WORKDIR="$(pwd)"
    local prot ligand prot_base lig_base combo_name comb_dir proc_gro itp_file
    export OMP_NUM_THREADS="$CPUS"

    for prot in "${PROTEINS[@]}"; do
        if [[ ! -f "$prot" ]]; then echo "Warning: missing protein file: $prot" >&2; continue; fi
        prot_base="$(basename "${prot%.*}")"
        for ligand in "${LIGAND_GROS[@]}"; do
            if [[ ! -f "$ligand" ]]; then echo "Warning: missing ligand gro: $ligand" >&2; continue; fi
            lig_base="$(basename "${ligand%.*}")"
            combo_name="${prot_base}__${lig_base}"
            if [[ -n "$OUTBASE" ]]; then
                comb_dir="$OUTBASE/$combo_name"
            else
                comb_dir="$WORKDIR/$combo_name"
            fi
            mkdir -p "$comb_dir"
            echo "=== Processing combo: $combo_name -> $comb_dir ==="
            echo "Using forcefield=$FORCEFIELD water=$WATERMODEL cpus=$CPUS" | tee -a "$comb_dir/run.log"

            # Protein preprocess for this combo: strip HETATM to avoid nonstandard residues (e.g., ligands like LMT)
            proc_gro="$comb_dir/${prot_base}_processed.gro"
            tmp_prot_pdb="$comb_dir/${prot_base}_protein_only.pdb"
            awk 'BEGIN{printed=0} /^ATOM/{print; printed=1} /^HETATM/{next} END{if(printed==0){exit 1}}' "$prot" > "$tmp_prot_pdb" 2>> "$comb_dir/run.log" || true
            echo "[pdb2gmx]" >> "$comb_dir/run.log"
            gmx pdb2gmx -f "$tmp_prot_pdb" -o "$proc_gro" -water "$WATERMODEL" -ff "$FORCEFIELD" -ignh -missing -chainsep id_or_ter >> "$comb_dir/run.log" 2>&1
            if [[ ! -s "$proc_gro" ]]; then echo "Error: failed pdb2gmx (see $comb_dir/run.log)" >&2; continue; fi

            # Topology
            if [[ -f topol.top ]]; then
                cp topol.top "$comb_dir/topol.top"
            else
                printf '%s\n' "; autogenerated minimal topol" > "$comb_dir/topol.top"
            fi
            if itp_file="$(find_itp_for "$ligand" 2>/dev/null || true)"; then
                if ! grep -q "#include \"${itp_file}\"" "$comb_dir/topol.top" 2>/dev/null; then
                    printf '#include "%s"\n' "$itp_file" >> "$comb_dir/topol.top"
                    echo "  included itp: $itp_file"
                fi
            else
                echo "Warning: no matching itp for $ligand; update $comb_dir/topol.top manually" >&2
            fi
            {
                printf "\n# --- auto-generated molecules block (verify!) ---\n"
                printf "[ molecules ]\n"
                printf "Protein    1\n"
                printf "%s    1\n" "$lig_base"
            } >> "$comb_dir/topol.top"

            # Build complex
            echo "[insert-molecules complex]" >> "$comb_dir/run.log"
            gmx insert-molecules -f "$proc_gro" -ci "$ligand" -nmol 1 -o "$comb_dir/complex.gro" >> "$comb_dir/run.log" 2>&1
            if [[ ! -s "$comb_dir/complex.gro" ]]; then echo "Error: failed to create complex.gro (see $comb_dir/run.log)" >&2; continue; fi

            # New box, solvate
            echo "[editconf newbox]" >> "$comb_dir/run.log"
            gmx editconf -f "$comb_dir/complex.gro" -o "$comb_dir/newbox.gro" -c -d 1.0 -bt cubic >> "$comb_dir/run.log" 2>&1
            if [[ ! -s "$comb_dir/newbox.gro" ]]; then echo "Error: failed to create newbox.gro (see $comb_dir/run.log)" >&2; continue; fi

            echo "[solvate]" >> "$comb_dir/run.log"
            gmx solvate -cp "$comb_dir/newbox.gro" -cs spc216.gro -o "$comb_dir/solvated.gro" -p "$comb_dir/topol.top" >> "$comb_dir/run.log" 2>&1
            if [[ ! -s "$comb_dir/solvated.gro" ]]; then echo "Error: failed to create solvated.gro (see $comb_dir/run.log)" >&2; continue; fi

            # Ions
            echo "[grompp ions]" >> "$comb_dir/run.log"
            gmx grompp -f ions.mdp -c "$comb_dir/solvated.gro" -p "$comb_dir/topol.top" -o "$comb_dir/ions.tpr" >> "$comb_dir/run.log" 2>&1
            if [[ ! -f "$comb_dir/ions.tpr" ]]; then echo "Error: failed to create ions.tpr (see $comb_dir/run.log)" >&2; continue; fi
            printf 'SOL\n' | gmx genion -s "$comb_dir/ions.tpr" -o "$comb_dir/ionized.gro" -p "$comb_dir/topol.top" -pname NA -nname CL -neutral >> "$comb_dir/run.log" 2>&1

            # EM
            echo "[grompp em]" >> "$comb_dir/run.log"
            gmx grompp -f em.mdp -c "$comb_dir/ionized.gro" -p "$comb_dir/topol.top" -o "$comb_dir/em.tpr" >> "$comb_dir/run.log" 2>&1
            if [[ ! -f "$comb_dir/em.tpr" ]]; then echo "Error: failed to create em.tpr (see $comb_dir/run.log)" >&2; continue; fi
            gmx mdrun -v -s "$comb_dir/em.tpr" -deffnm "$comb_dir/em" -nt "$CPUS" >> "$comb_dir/run.log" 2>&1

            # NVT
            echo "[grompp nvt]" >> "$comb_dir/run.log"
            gmx grompp -f nvt.mdp -c "$comb_dir/em.gro" -r "$comb_dir/em.gro" -p "$comb_dir/topol.top" -o "$comb_dir/nvt.tpr" >> "$comb_dir/run.log" 2>&1
            if [[ ! -f "$comb_dir/nvt.tpr" ]]; then echo "Error: failed to create nvt.tpr (see $comb_dir/run.log)" >&2; continue; fi
            gmx mdrun -s "$comb_dir/nvt.tpr" -deffnm "$comb_dir/nvt" -nt "$CPUS" >> "$comb_dir/run.log" 2>&1

            # NPT
            echo "[grompp npt]" >> "$comb_dir/run.log"
            gmx grompp -f npt.mdp -c "$comb_dir/nvt.gro" -t "$comb_dir/nvt.cpt" -p "$comb_dir/topol.top" -o "$comb_dir/npt.tpr" >> "$comb_dir/run.log" 2>&1
            if [[ ! -f "$comb_dir/npt.tpr" ]]; then echo "Error: failed to create npt.tpr (see $comb_dir/run.log)" >&2; continue; fi
            gmx mdrun -s "$comb_dir/npt.tpr" -deffnm "$comb_dir/npt" -nt "$CPUS" >> "$comb_dir/run.log" 2>&1

            # MD
            echo "[grompp md]" >> "$comb_dir/run.log"
            gmx grompp -f md.mdp -c "$comb_dir/npt.gro" -t "$comb_dir/npt.cpt" -p "$comb_dir/topol.top" -o "$comb_dir/md.tpr" >> "$comb_dir/run.log" 2>&1
            if [[ ! -f "$comb_dir/md.tpr" ]]; then echo "Error: failed to create md.tpr (see $comb_dir/run.log)" >&2; continue; fi
            gmx mdrun -s "$comb_dir/md.tpr" -deffnm "$comb_dir/md" -nt "$CPUS" >> "$comb_dir/run.log" 2>&1

            echo "=== Finished combo: $combo_name. Results in: $comb_dir ==="
        done
    done
}
