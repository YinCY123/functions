analyze_md_pub() {
    ############################################################
    #  Molecular Dynamics Analysis Function (GNU‑style flags)
    #
    #  DESCRIPTION:
    #    Performs a full suite of downstream MD analyses for a
    #    protein–ligand complex using GROMACS. Outputs include:
    #    RMSD, RMSF, Rg, SASA, H‑bonds, COM distance, interaction
    #    energy, DSSP secondary structure, clustering, and
    #    MM/PBSA frame extraction.
    #
    #  PARAMETERS (GNU‑style):
    #
    #    --traj <file>
    #        Path to the trajectory file (.xtc, .trr).
    #        Default: md.xtc
    #        This is the main MD trajectory used for all analyses.
    #
    #    --tpr <file>
    #        Path to the GROMACS run input file (.tpr).
    #        Default: md.tpr
    #        Required for structural reference, topology, and
    #        energy term extraction.
    #
    #    --ndx <file>
    #        Path to the index file (.ndx).
    #        Default: index.ndx
    #        Must contain groups: Protein, Backbone, and ligand group.
    #
    #    --lig <group>
    #        Name of the ligand group in the index file.
    #        Default: LIG
    #        Used for RMSD, H‑bond, COM distance, and energy analysis.
    #
    #    --help / -h
    #        Display usage information.
    #
    ############################################################

    # --- Default values ---
    traj="md.xtc"
    tpr="md.tpr"
    ndx="index.ndx"
    lig="LIG"

    # --- Parse GNU-style flags ---
    while [[ $# -gt 0 ]]; do
        case "$1" in
            --traj)
                traj="$2"
                shift 2
                ;;
            --tpr)
                tpr="$2"
                shift 2
                ;;
            --ndx)
                ndx="$2"
                shift 2
                ;;
            --lig)
                lig="$2"
                shift 2
                ;;
            --help|-h)
                echo "Usage: analyze_md_pub [--traj file] [--tpr file] [--ndx file] [--lig group]"
                return 0
                ;;
            *)
                echo "Unknown option: $1"
                return 1
                ;;
        esac
    done

    mkdir -p rmsd rmsf rog sasa hbond dist energy dssp cluster mmpbsa

    # RMSD: protein backbone + ligand
    printf "Backbone\nBackbone\n" | gmx rms -s "$tpr" -f "$traj" -o rmsd/protein_backbone_rmsd.xvg
    printf "$lig\n$lig\n" | gmx rms -s "$tpr" -f "$traj" -o rmsd/ligand_rmsd.xvg

    # RMSF
    printf "Backbone\n" | gmx rmsf -s "$tpr" -f "$traj" -o rmsf/protein_rmsf.xvg -res

    # Radius of gyration
    printf "Protein\n" | gmx gyrate -s "$tpr" -f "$traj" -o rog/protein_rog.xvg

    # SASA
    printf "Protein\n" | gmx sasa -s "$tpr" -f "$traj" -o sasa/protein_sasa.xvg -or sasa/residue_sasa.xvg

    # Hydrogen bonds
    printf "Protein\n$lig\n" | gmx hbond -s "$tpr" -f "$traj" \
        -num hbond/protein_lig_hbnum.xvg -hbn hbond/hbmap.ndx

    # COM distance
    gmx distance -s "$tpr" -f "$traj" \
        -select "com of group \"Protein\" plus com of group \"$lig\"" \
        -oall dist/protein_lig_comdist.xvg

    # Interaction energy
    printf "Protein\n$lig\n" | gmx energy -f md.edr -o energy/protein_lig_energy.xvg

    # DSSP
    gmx do_dssp -s "$tpr" -f "$traj" -o dssp/ss.xpm -sc dssp/ss_count.xvg

    # Clustering
    printf "Backbone\n" | gmx cluster -s "$tpr" -f "$traj" \
        -o cluster/cluster.xpm -g cluster/cluster.log \
        -dist cluster/cluster-dist.xvg -method gromos

    # MM/PBSA frames
    printf "0\n" | gmx trjconv -s "$tpr" -f "$traj" -o mmpbsa/frames.pdb -sep

    echo "Publication-grade analysis complete."
}
