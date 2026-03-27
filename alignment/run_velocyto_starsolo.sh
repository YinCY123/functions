#!/bin/bash

# Run velocyto on one or more STARsolo BAMs using `velocyto run` (not run10x).
# Requirements: coordinate-sorted BAM with CB and UB tags; barcodes file listing
# allowed cells (same strings as the CB tag). Typical STARsolo paths:
#   <STAR_out>/Aligned.sortedByCoord.out.bam
#   <STAR_out>/Solo.out/Gene/filtered/barcodes.tsv

run_velocyto_starsolo() {
    if [[ ${1:-} == "-h" || ${1:-} == "--help" ]]; then
        echo ""
        echo "Run velocyto run on STARsolo-produced BAM file(s)."
        echo ""
        echo "  --gtf FILE              genome annotation (required)"
        echo "  -m, --mask FILE         repeat mask .gtf (recommended)"
        echo ""
        echo "  One of:"
        echo "  --manifest FILE         TSV: bam_path<TAB>bcfile_path<TAB>out_dir (out_dir optional)"
        echo "  --star-dirs PARENT      run each subdirectory of PARENT as a STAR output folder:"
        echo "                          <subdir>/Aligned.sortedByCoord.out.bam"
        echo "                          <subdir>/Solo.out/Gene/filtered/barcodes.tsv"
        echo "  --bam FILE -b FILE      single sample (use with --bcfile and --outputfolder)"
        echo ""
        echo "  --bcfile FILE           with --bam only: filtered barcodes.tsv"
        echo "  -o, --outputfolder DIR  with --bam only: output directory for .loom"
        echo "  -e, --sampleid NAME     sample id embedded in loom (default: basename of out dir)"
        echo "  -l, --logic TEXT        filtering logic (default: Default)"
        echo "  --samtools-threads N    samtools sort threads (default: 8)"
        echo "  --samtools-memory MB    MB per samtools thread (default: 2048)"
        echo "  -t, --dtype TEXT        loom layer dtype (default: uint32)"
        echo "  -d, --dump TEXT         debug dump (default: 0)"
        echo "  --bam-subpath REL       under each --star-dirs subdir, BAM path (see default below)"
        echo "  --bc-subpath REL        under each --star-dirs subdir, barcodes path (see default below)"
        echo ""
        echo "Default --bam-subpath:  Aligned.sortedByCoord.out.bam"
        echo "Default --bc-subpath:   Solo.out/Gene/filtered/barcodes.tsv"
        echo ""
        return 0
    fi

    local gtf
    local m
    local l=Default
    local samtools_threads=8
    local samtools_memory=2048
    local t=uint32
    local d=0
    local manifest
    local star_dirs_parent
    local bam
    local bcfile
    local out
    local sampleid
    local bam_subpath="Aligned.sortedByCoord.out.bam"
    local bc_subpath="Solo.out/Gene/filtered/barcodes.tsv"

    while [[ $# -gt 0 ]]; do
        case $1 in
            --gtf)
                gtf="$2"
                shift 2
                ;;
            -m|--mask)
                m="$2"
                shift 2
                ;;
            -l|--logic)
                l="$2"
                shift 2
                ;;
            --samtools-threads)
                samtools_threads="$2"
                shift 2
                ;;
            --samtools_memory|--samtools-memory)
                samtools_memory="$2"
                shift 2
                ;;
            -t|--dtype)
                t="$2"
                shift 2
                ;;
            -d|--dump)
                d="$2"
                shift 2
                ;;
            --manifest)
                manifest="$2"
                shift 2
                ;;
            --star-dirs)
                star_dirs_parent="$2"
                shift 2
                ;;
            --bam)
                bam="$2"
                shift 2
                ;;
            -b|--bcfile)
                bcfile="$2"
                shift 2
                ;;
            -o|--outputfolder)
                out="$2"
                shift 2
                ;;
            -e|--sampleid)
                sampleid="$2"
                shift 2
                ;;
            --bam-subpath)
                bam_subpath="$2"
                shift 2
                ;;
            --bc-subpath)
                bc_subpath="$2"
                shift 2
                ;;
            *)
                echo "unknown parameter: $1"
                return 1
                ;;
        esac
    done

    if [[ -z "$gtf" || ! -f "$gtf" ]]; then
        echo "Error: --gtf must point to an existing annotation file"
        return 1
    fi

    _run_one() {
        local bam_path="$1"
        local bc_path="$2"
        local out_dir="$3"
        local sid="$4"

        if [[ ! -f "$bam_path" ]]; then
            echo "Error: BAM not found: $bam_path"
            return 1
        fi
        if [[ ! -f "$bc_path" ]]; then
            echo "Error: barcodes file not found: $bc_path"
            return 1
        fi
        mkdir -p "$out_dir"

        local mask_args=()
        if [[ -n "$m" ]]; then
            mask_args=( -m "$m" )
        fi

        echo "Processing: $bam_path -> $out_dir"
        local e_args=()
        if [[ -n "$sid" ]]; then
            e_args=( -e "$sid" )
        fi

        velocyto run \
            -b "$bc_path" \
            -o "$out_dir" \
            "${mask_args[@]}" \
            -l "$l" \
            --samtools-threads "$samtools_threads" \
            --samtools-memory "$samtools_memory" \
            -t "$t" \
            -d "$d" \
            "${e_args[@]}" \
            -vvv \
            "$bam_path" \
            "$gtf"
    }

    local mode_count=0
    [[ -n "$manifest" ]] && mode_count=$((mode_count + 1))
    [[ -n "$star_dirs_parent" ]] && mode_count=$((mode_count + 1))
    [[ -n "$bam" ]] && mode_count=$((mode_count + 1))

    if [[ "$mode_count" -ne 1 ]]; then
        echo "Error: specify exactly one of --manifest, --star-dirs, or --bam"
        return 1
    fi

    if [[ -n "$bam" ]]; then
        if [[ -z "$bcfile" || ! -f "$bcfile" ]]; then
            echo "Error: --bam requires --bcfile pointing to an existing barcodes.tsv"
            return 1
        fi
        if [[ -z "$out" ]]; then
            echo "Error: --bam requires -o/--outputfolder"
            return 1
        fi
        local sid_arg="$sampleid"
        if [[ -z "$sid_arg" ]]; then
            sid_arg="$(basename "$out")"
        fi
        _run_one "$bam" "$bcfile" "$out" "$sid_arg"
        return $?
    fi

    if [[ -n "$manifest" ]]; then
        if [[ ! -f "$manifest" ]]; then
            echo "Error: manifest not found: $manifest"
            return 1
        fi
        local n=0
        while IFS=$'\t' read -r f1 f2 f3 _rest || [[ -n "${f1:-}" ]]; do
            [[ -z "${f1// }" ]] && continue
            [[ "$f1" =~ ^# ]] && continue
            n=$((n + 1))
            if [[ -z "$f1" || -z "$f2" ]]; then
                echo "Error: manifest line $n needs at least two columns (bam, bcfile)"
                return 1
            fi
            local odir
            if [[ -n "$f3" ]]; then
                odir="$f3"
            else
                odir="$(dirname "$f1")/velocyto"
            fi
            local sid_line="$sampleid"
            if [[ -z "$sid_line" ]]; then
                sid_line="$(basename "$odir")"
            fi
            _run_one "$f1" "$f2" "$odir" "$sid_line" || return 1
        done < "$manifest"
        return 0
    fi

    if [[ -n "$star_dirs_parent" ]]; then
        if [[ ! -d "$star_dirs_parent" ]]; then
            echo "Error: --star-dirs not a directory: $star_dirs_parent"
            return 1
        fi
        local sub
        shopt -s nullglob
        for sub in "$star_dirs_parent"/*; do
            [[ -d "$sub" ]] || continue
            local bb="$sub/$bam_subpath"
            local bc="$sub/$bc_subpath"
            local odir="$sub/velocyto_out"
            local sid_sub="$sampleid"
            if [[ -z "$sid_sub" ]]; then
                sid_sub="$(basename "$sub")"
            fi
            _run_one "$bb" "$bc" "$odir" "$sid_sub" || return 1
        done
        return 0
    fi

    return 1
}
