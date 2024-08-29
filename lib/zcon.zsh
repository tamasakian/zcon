#!/usr/bin/env zsh

#### Function library for ZCON

function generate_pep_fasta() {
    ## This function generates a fasta file of protein.
    function usage() {
        cat <<EOS
Usage: generate_pep_fasta <arg1> <arg2> <arg3> <arg4> <arg5>

    arg1: num_orgs
    arg2: orgs
    arg3: num_peps
    arg4: pep_key
    arg5: pep_value

EOS
        exit 1
    }

    function parse_args() {
        if [[ $# -lt 5 ]]; then
            usage
        fi

        num_orgs=$1
        shift
        orgs=()
        for ((i=1; i<=num_orgs; i++)); do
            orgs[i]="${1// /_}"
            shift
        done

        num_peps=$1
        shift
        typeset -g -A peps
        for ((i=1; i<=num_peps; i++)); do
            key="$1"
            shift
            value="$1"
            shift
            peps[$key]="$value"
        done
    }

    function merge_fasta() {
        touch "${taskdir}/database.fasta"
        for org in "${orgs[@]}"; do
            genus=${org%%_*}
            cat "${DATA}/${genus}/${org}.pep.all.fasta" >> "${taskdir}/database.fasta"
        done
    }

    function generate_fasta() {
        touch "${taskdir}/pep.fasta"
        makeblastdb \
            -in "${taskdir}/database.fasta" \
            -dbtype prot \
            -hash_index \
            -parse_seqids

        for pep_key in ${(@k)peps}; do
            tmpfile=$(mktemp)
            blastdbcmd \
                -entry "$pep_key" \
                -db "${taskdir}/database.fasta" \
                -out "$tmpfile"
            python3 -m fasp rename_header \
                "$tmpfile" \
                "$tmpfile" \
                "${peps[$pep_key]}" \
                ""
            cat "$tmpfile" >> "${taskdir}/pep.fasta"
            rm "$tmpfile"
        done
    }

    function main() {
        parse_args "$@"
        make_taskdir
        merge_fasta
        generate_fasta
    }

    main "$@"
}


function construct_pep_msa() {
    ## This function constructs a multiple sequence alignment of protein.
    ## This is based on generate_pep_fasta.
    function usage() {
        cat <<EOS
Usage: construct_pep_msa <arg1> <arg2> <arg3> <arg4> <arg5>

    arg1: num_orgs
    arg2: orgs
    arg3: num_peps
    arg4: pep_key
    arg5: pep_value

EOS
        exit 1
    }

    function construct_msa() {
        mafft --auto "${taskdir}/pep.fasta" > "${taskdir}/pep.aln.fasta"
        mafft --auto --clustalout "${taskdir}/pep.fasta" > "${taskdir}/pep.aln.clustal"
    }

    function main() {
        generate_pep_fasta "$@"
        construct_msa
    }

    main "$@"
}


function trim_pep_msa() {
    ## This function trims a multiple sequence alignment of protein.
    ## This is based on generate_pep_fasta, construct_pep_msa.
    function usage() {
        cat <<EOS
Usage: trim_pep_msa <arg1> <arg2> <arg3> <arg4> <arg5>

    arg1: num_orgs
    arg2: orgs
    arg3: num_peps
    arg4: pep_key
    arg5: pep_value

EOS
        exit 1
    }

    function trim_msa() {
        trimal \
            -in "${taskdir}/pep.aln.fasta" \
            -out "${taskdir}/pep.aln.trim.fasta" \
            -automated1 \
            -htmlout "${taskdir}/pep.aln.trim.html"
    }

    function main() {
        construct_pep_msa "$@"
        trim_msa
    }

    main "$@"
}


function estimate_pep_mltree() {
    ## This function estimates a gene tree by maximum likelihood method.
    ## This is based on generate_pep_fasta, construct_pep_msa, trim_pep_msa.
    function usage() {
        cat <<EOS
Usage: estimate_pep_mltree <arg1> <arg2> <arg3> <arg4> <arg5>

    arg1: num_orgs
    arg2: orgs
    arg3: num_peps
    arg4: pep_key
    arg5: pep_value

EOS
        exit 1
    }

    function estimate_mltree() {
        raxml-ng \
            --msa "${taskdir}/pep.aln.trim.fasta" \
            --all \
            --model Blosum62 \
            --bs-trees 1000 \
            --threads 6
        iqtree2 \
            -s "${taskdir}/pep.aln.trim.fasta" \
            --seqtype AA \
            -m MFP \
            -bb 1000 \
            -alrt 1000

    }

    function main() {
        trim_pep_msa "$@"
        estimate_mltree
    }

    main "$@"
}


function find_pep_motif() {
    ## This function finds individual motif occurrences with generated pep fasta.
    ## This is based on generate_pep_fasta.
    function usage_2() {
        cat <<EOS
Usage: find_pep_motif <arg1> <arg2> <arg3> <arg4> <arg5> <arg6> <arg7>

    arg1: num_orgs
    arg2: orgs
    arg3: num_peps
    arg4: pep_key
    arg5: pep_value
    arg6: motif_db
    arg7: qvalue

EOS
        exit 1
    }

    function parse_args_2() {
        ## This process is for parse_args.
        shift; shift
        for ((i=1; i<=num_orgs; i++)); do
            shift
        done
        for ((i=1; i<=num_peps; i++)); do
            shift; shift
        done

        if [[ $# -lt 2 ]]; then
            usage_2
        fi
        motif_db=$1
        shift
        qvalue=$1
        shift
    }

    function find_motif() {
        fimo \
            --o "${taskdir}/fimo_out" \
            --thresh "$qvalue" \
            --qv-thresh \
            --verbosity 1 \
            "${DATA}/motif_databases/${motif_db}" \
            "${taskdir}/pep.fasta"
    }

    function main() {
        generate_pep_fasta "$@"
        parse_args_2 "$@"
        find_motif
    }

    main "$@"
}


function generate_flanking_dna_fasta() {
    ## This function generates a fasta file of flanking region of genes.
    function usage() {
        cat <<EOS
Usage: generate_flanking_dna_fasta <arg1> <arg2> <arg3> <arg4> <arg5> <arg6>

    arg1: num_orgs
    arg2: orgs
    arg3: num_peps
    arg4: pep_key
    arg5: pep_value
    arg6: bp

EOS
        exit 1
    }

    function parse_args() {
        if [[ $# -lt 5 ]]; then
            usage
        fi

        num_orgs=$1
        shift
        orgs=()
        for ((i=1; i<=num_orgs; i++)); do
            orgs[i]="${1// /_}"
            shift
        done

        num_peps=$1
        shift
        typeset -g -A peps
        for ((i=1; i<=num_peps; i++)); do
            key="$1"
            shift
            value="$1"
            shift
            peps[$key]="$value"
        done

        bp=$1
        shift
    }

    function merge_fasta() {
        touch "${taskdir}/database.fasta"
        for org in "${orgs[@]}"; do
            genus=${org%%_*}
            cat "${DATA}/${genus}/${org}.dna.toplevel.fasta" >> "${taskdir}/database.fasta"
        done
    }

    function merge_gff() {
        touch "${taskdir}/database.gff"
        for org in "${orgs[@]}"; do
            genus=${org%%_*}
            cat "${DATA}/${genus}/${org}.genome.gff" >> "${taskdir}/database.gff"
        done
    }

    function generate_fasta() {
        touch "${taskdir}/flanking_dna.fasta"
        makeblastdb \
            -in "${taskdir}/database.fasta" \
            -dbtype nucl \
            -hash_index \
            -parse_seqids

        for pep_key in ${(@k)peps}; do
            tmpfile=$(mktemp)
            seq="" strand="" start="" end=""
            output=$(python3 -m biotp output_seqid_strand_locs_by_pepid "${taskdir}/database.gff" "$pep_key")
            read seq strand start end <<< "$output"
            blastdbcmd \
                -entry "$seq" \
                -db "${taskdir}/database.fasta" \
                -out "$tmpfile"
            python3 -m fasp slice_sequence_by_flanking_region \
                "$tmpfile" \
                "$tmpfile" \
                "${peps[$pep_key]}" \
                "${bp}bp_flanking_region" \
                "$strand" \
                "$start" \
                "$end" \
                "$bp"
            cat "$tmpfile" >> "${taskdir}/flanking_dna.fasta"
            rm "$tmpfile"
        done
    }

    function main() {
        parse_args "$@"
        make_taskdir
        merge_fasta
        merge_gff
        generate_fasta
    }

    main "$@"
}


function construct_flanking_dna_msa() {
    ## This function construct a multiple sequence alighment of flanking region of genes.
    ## This is based on generate_flanking_dna_fasta.
    function usage() {
        cat <<EOS
Usage: construct_flanking_dna_msa <arg1> <arg2> <arg3> <arg4> <arg5> <arg6>

    arg1: num_orgs
    arg2: orgs
    arg3: num_peps
    arg4: pep_key
    arg5: pep_value
    arg6: bp

EOS
        exit 1
    }

    function construct_msa() {
        mafft --auto "${taskdir}/flanking_dna.fasta" > "${taskdir}/flanking_dna.aln.fasta"
        mafft --auto --clustalout "${taskdir}/flanking_dna.fasta" > "${taskdir}/flanking_dna.aln.clustal"
    }

    function main() {
        generate_flanking_dna_fasta "$@"
        construct_msa
    }

    main "$@"
}


function generate_upstream_dna_fasta() {
    ## This function generates a fasta file of upstream region of genes.
    function usage() {
        cat <<EOS
Usage: generate_upstream_dna_fasta <arg1> <arg2> <arg3> <arg4> <arg5> <arg6>

    arg1: num_orgs
    arg2: orgs
    arg3: num_peps
    arg4: pep_key
    arg5: pep_value
    arg6: bp

EOS
        exit 1
    }

    function parse_args() {
        if [[ $# -lt 6 ]]; then
            usage
        fi

        num_orgs=$1
        shift
        orgs=()
        for ((i=1; i<=num_orgs; i++)); do
            orgs[i]="${1// /_}"
            shift
        done

        num_peps=$1
        shift
        typeset -g -A peps
        for ((i=1; i<=num_peps; i++)); do
            key="$1"
            shift
            value="$1"
            shift
            peps[$key]="$value"
        done

        bp=$1
        shift
    }

    function merge_fasta() {
        touch "${taskdir}/database.fasta"
        for org in "${orgs[@]}"; do
            genus=${org%%_*}
            cat "${DATA}/${genus}/${org}.dna.toplevel.fasta" >> "${taskdir}/database.fasta"
        done
    }

    function merge_gff() {
        touch "${taskdir}/database.gff"
        for org in "${orgs[@]}"; do
            genus=${org%%_*}
            cat "${DATA}/${genus}/${org}.genome.gff" >> "${taskdir}/database.gff"
        done
    }

    function generate_fasta() {
        touch "${taskdir}/upstream_dna.fasta"
        makeblastdb \
            -in "${taskdir}/database.fasta" \
            -dbtype nucl \
            -hash_index \
            -parse_seqids

        for pep_key in ${(@k)peps}; do
            tmpfile=$(mktemp)
            seq="" strand="" start="" end=""
            output=$(python3 -m biotp output_seqid_strand_locs_by_pepid "${taskdir}/database.gff" "$pep_key")
            read seq strand start end <<< "$output"
            blastdbcmd \
                -entry "$seq" \
                -db "${taskdir}/database.fasta" \
                -out "$tmpfile"
            python3 -m fasp slice_sequence_by_upstream_region \
                "$tmpfile" \
                "$tmpfile" \
                "${peps[$pep_key]}" \
                "${bp}bp_upstream_region" \
                "$strand" \
                "$start" \
                "$end" \
                "$bp"
            cat "$tmpfile" >> "${taskdir}/upstream_dna.fasta"
            rm "$tmpfile"
        done
    }

    function main() {
        parse_args "$@"
        make_taskdir
        merge_fasta
        merge_gff
        generate_fasta
    }

    main "$@"
}


function generate_downstream_dna_fasta() {
    ## This function generates a fasta file of downstream region of genes.
    function usage() {
        cat <<EOS
Usage: generate_downstream_dna_fasta <arg1> <arg2> <arg3> <arg4> <arg5> <arg6>

    arg1: num_orgs
    arg2: orgs
    arg3: num_peps
    arg4: pep_key
    arg5: pep_value
    arg6: bp

EOS
        exit 1
    }

    function parse_args() {
        if [[ $# -lt 6 ]]; then
            usage
        fi

        num_orgs=$1
        shift
        orgs=()
        for ((i=1; i<=num_orgs; i++)); do
            orgs[i]="${1// /_}"
            shift
        done

        num_peps=$1
        shift
        typeset -g -A peps
        for ((i=1; i<=num_peps; i++)); do
            key="$1"
            shift
            value="$1"
            shift
            peps[$key]="$value"
        done

        bp=$1
        shift
    }

    function merge_fasta() {
        touch "${taskdir}/database.fasta"
        for org in "${orgs[@]}"; do
            genus=${org%%_*}
            cat "${DATA}/${genus}/${org}.dna.toplevel.fasta" >> "${taskdir}/database.fasta"
        done
    }

    function merge_gff() {
        touch "${taskdir}/database.gff"
        for org in "${orgs[@]}"; do
            genus=${org%%_*}
            cat "${DATA}/${genus}/${org}.genome.gff" >> "${taskdir}/database.gff"
        done
    }

    function generate_fasta() {
        touch "${taskdir}/downstream_dna.fasta"
        makeblastdb \
            -in "${taskdir}/database.fasta" \
            -dbtype nucl \
            -hash_index \
            -parse_seqids

        for pep_key in ${(@k)peps}; do
            tmpfile=$(mktemp)
            seq="" strand="" start="" end=""
            output=$(python3 -m biotp output_seqid_strand_locs_by_pepid "${taskdir}/database.gff" "$pep_key")
            read seq strand start end <<< "$output"
            blastdbcmd \
                -entry "$seq" \
                -db "${taskdir}/database.fasta" \
                -out "$tmpfile"
            python3 -m fasp slice_sequence_by_downstream_region \
                "$tmpfile" \
                "$tmpfile" \
                "${peps[$pep_key]}" \
                "${bp}bp_downstream_region" \
                "$strand" \
                "$start" \
                "$end" \
                "$bp"
            cat "$tmpfile" >> "${taskdir}/downstream_dna.fasta"
            rm "$tmpfile"
        done
    }

    function main() {
        parse_args "$@"
        make_taskdir
        merge_fasta
        merge_gff
        generate_fasta
    }

    main "$@"
}

function generate_intron_dna_fasta() {
    ## This function generates a fasta file of intron of genes.
    function usage() {
        cat <<EOS
Usage:  generate_intron_dna_fasta <arg1> <arg2> <arg3> <arg4> <arg5> <arg6>

    arg1: num_orgs
    arg2: orgs
    arg3: num_peps
    arg4: pep_key
    arg5: pep_value
    arg6: num_introns
    
EOS
        exit 1
    }

    function parse_args() {
        if [[ $# -lt 6 ]]; then
            usage
        fi

        num_orgs=$1
        shift
        orgs=()
        for ((i=1; i<=num_orgs; i++)); do
            orgs[i]="${1// /_}"
            shift
        done

        num_peps=$1
        shift
        typeset -g -A peps
        for ((i=1; i<=num_peps; i++)); do
            key="$1"
            shift
            value="$1"
            shift
            peps[$key]="$value"
        done
        num_introns=$1
        shift
    }

    function merge_fasta() {
        touch "${taskdir}/database.fasta"
        for org in "${orgs[@]}"; do
            genus=${org%%_*}
            cat "${DATA}/${genus}/${org}.dna.toplevel.fasta" >> "${taskdir}/database.fasta"
        done
    }

    function merge_gff() {
        touch "${taskdir}/database.gff"
        for org in "${orgs[@]}"; do
            genus=${org%%_*}
            cat "${DATA}/${genus}/${org}.genome.gff" >> "${taskdir}/database.gff"
        done
    }

    function generate_fasta() {
        python3 -m biotp generate_coordinate_all_introns \
            "${taskdir}/database.gff" \
            "${taskdir}/database_introns.gff"
        python3 -m biotp generate_all_introns \
            "${taskdir}/database.fasta" \
            "${taskdir}/database_introns.gff" \
            "${taskdir}/database_introns.fasta"
        makeblastdb \
            -in "${taskdir}/database_introns.fasta" \
            -dbtype nucl \
            -hash_index \
            -parse_seqids

        for ((i=1; i<=num_introns; i++)); do
            touch "${taskdir}/intron_dna_${i}.fasta"
            for pep_key in ${(@k)peps}; do
                tmpfile=$(mktemp)
                blastdbcmd \
                    -entry "${pep_key}_intron_${i}" \
                    -db "${taskdir}/database_introns.fasta" \
                    -out "$tmpfile"
                python3 -m fasp rename_header \
                    "$tmpfile" \
                    "$tmpfile" \
                    "${peps[$pep_key]}" \
                    "intron_${i}" 
                cat "$tmpfile" >> "${taskdir}/intron_dna_${i}.fasta"
                rm "$tmpfile"
            done
        done
    }

    function main() {
        parse_args "$@"
        make_taskdir
        merge_fasta
        merge_gff
        generate_fasta
    }

    main "$@"
}


