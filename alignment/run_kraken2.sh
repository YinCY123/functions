run_kraken2(){
    # help information
    if [[ "$1" == "-h" || "$1" == "--help" ]];then
        cat << EOF
        Kraken2 help information:

        --db                    kraken database (default: ~/ssd/kraken2/database/k2_minusb_20250402).
        --threads               number of threads to use (default: 8).
        --quick                 quick operation (use first hit or hits).
        --unclassified-out      print unclassified sequences to filename.
        --classified-out        print classified sequences to filename.
        --output                print output to filename (default: stdout); "-" will suppress normal output.
        --confidence            confidence score threshold (deafult: 0.0), must between [0,1].
        --minimum-base-quality  minmum base quality used in classification (default: 0).
        --report                print a report with aggregate counts/clade to file.
        --use-mpa-style         with --report, format report output like kraken 1's kraken-mpa-report.
        --report-zero-counts    with --report, report counts for ALL taxa, even if counts are zero.
        --report-mimimizer-data with --report, report minimizer and distinct minimizer count information in addition to normal kraken report.
        --memory-mapping        avoids loading database into RAM.
        --paired                if the input files are from paired sequencing (default: paired).
        --use-names             print scientific names instead of just taxids.
        --compression           which compression method used to compress the file (default gzip-compressed).
        --minimum-hit-groups    minimum number of hit groups needed to make a call (default: 2).

EOF
        return 0
    fi

    # initialize variable
    local db="/home/yincy/ssd/kraken2/database/k2_minusb_20250402"
    local threads=8


    while [ $# -ge 1 ]; do
        case $1 in
            --db)
                db="$2"
                shift 2
                ;;
            --threads)
                threads="$2"
                shift 2
                ;;
            --unclassified-out)
                unclassified-out="$2"
                shift 2
                ;;
            --classified-out)
                classified-out="$2"
                shift 2
                ;;
            --output)
                output="$2"
                shift 2
                ;;
            --confidence)
                confidence="$2"
                shift 2
                ;;
            --minimum-base-quality)
                minimum-base-quality="$2"
                shift 2
                ;;
            --compression)
                case $2 in
                    gzip|bzip2)
                        compression="$2";;
                    *) echo "The compression is not supported."; return 1 ;;
                esac
                shift 2
                ;;
            --minimum_hit_groups)
                minimum_hit_groups="$2"
                shift 2
                ;;
            --memory_mapping)
                memory_mapping="$2"
                shift 2
                ;;
            --paired)
                paired="$2"
                shift 2
                ;;
            --use_names)
                use_names="$2"
                shift 2
                ;;
            *)
                echo "the parameter is not supported!"
        esac
    done
}