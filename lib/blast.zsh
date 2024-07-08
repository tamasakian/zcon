#!/usr/bin/env zsh

#### Function libs with BLAST

function makeblastdb_genome_by_genus() {
    ## This function is depend on fetch_genome_by_genus, datasets.zsh
    function usage() {
        cat <<EOS
Usage:  makeblastdb_genome_by_genus <arg1>

    arg1: genus

EOS
        exit 1
    }
    
    function makeblastdb_genome() {
        for org in ${org_li[*]}; do
            makeblastdb \
                -in "${DATA}/${genus}/${org}.dna.toplevel.fasta" \
                -dbtype nucl -hash_index -parse_seqids
            makeblastdb \
                -in "${DATA}/${genus}/${org}.cds.all.fasta" \
                -dbtype nucl -hash_index -parse_seqids
            makeblastdb \
                -in "${DATA}/${genus}/${org}.pep.all.fasta" \
                -dbtype prot -hash_index -parse_seqids
        done
    }

    function main() {
        fetch_genome_by_genus "$@"
        makeblastdb_genome
    }

    main "$@"
}

function output_blastp_genus_symbol() {
    function usage() {
        cat <<EOS
Usage:  output_blastp_genus_symbol <arg1> <arg2> <arg3> <arg4>

    arg1: genus
    arg2: symbol
    arg3: symbol_org
    arg4: evalue

EOS
        exit 1
    }

    function parse_args() {
        if [[ $# != 4 ]]; then
            usage
        fi
        genus="$1"
        symbol="$2"
        symbol_org="${3/ /_}"
        evalue=$4
    }

    function make_dir() {
        taskdir=$(make_dir_by_date $TASKFILE)
        echo "taskdir: ${taskdir}"
    }

    function blastp_genus_symbol() {
        echo "Protein BLAST"
        touch "${taskdir}/${genus}.${symbol}.blastp"
        for org in ${org_li[*]}; do
            echo "ref: all [${org}], qry: ${symbol} [${symbol_org}]"
            tmpfile=$(mktemp)
            blastp \
                -outfmt 7 -evalue $evalue \
                -db "${DATA}/${genus}/${org}.pep.all.fasta" \
                -query "${DATA}/${symbol_org}/${symbol}.pep.fasta" \
                -out "$tmpfile"
            cat "$tmpfile" >> "${taskdir}/${genus}.${symbol}.blastp"
            rm "$tmpfile"
        done
    }

    function main() {
        parse_args "$@"
        make_dir
        redeclare_genome_by_genus "$genus"
        blastp_genus_symbol
    }

    main "$@"
}

function construct_pep_blastp_genus_symbol() {
    function usage() {
        cat <<EOS
Usage:  construct_pep_blastp_genus_symbol <arg1> <arg2> <arg3> <arg4> (<arg5> <arg6> ...)

    arg1: genus
    arg2: symbol
    arg3: symbol_org
    arg4: evalue
    arg5: pid
    arg6: pnm
    ...

EOS
        exit 1
    } 

    function parse_args() {
        if [[ $# -lt 6 ]]; then
            usage
        fi
        genus="$1"
        symbol="$2"
        symbol_org="${3/ /_}"
        evalue=$4
        local _pli=("${@:5}")
        for ((i=1; i<${#_pli[@]}; i+=2)); do 
            echo "${i}, ${_pli[i]}, ${_pli[i+1]}"
            pid_li+=("${_pli[i]}")
            pnm_li+=("${_pli[i+1]}")
        done

        echo "${pid_li[1]}"
    }

    function make_dir() {
        taskdir=$(make_dir_by_date $TASKFILE)
        echo "taskdir: ${taskdir}"
    }

    function retrieve_blastp() {
        touch "${taskdir}/${genus}.${symbol}.pep.fasta"
        echo "Protein BLAST"
        for org in ${org_li[*]}; do
            echo "ref: all [${org}], qry: ${symbol} [${symbol_org}]"
            local _blastp=$(blastp -outfmt 6 -evalue $evalue -db ${DATA}/${genus}/${org}.pep.all.fasta -query ${DATA}/${symbol_org}/${symbol}.pep.fasta)
            local _ids=$(echo "$_blastp" | cut -f 2 | sort -u)
            if [ -z "$_ids" ]; then
                echo "No hit."
                continue
            fi
            echo "Hit."; echo "$_ids"
            for _id in ${(f)_ids}; do
                for ((i=1; i<=${#pid_li[@]}; ++i)); do
                    if [[ "$_id" != "${pid_li[i]}" ]]; then
                        continue
                    fi
                    tmpfile=$(mktemp)
                    blastdbcmd \
                        -entry "$_id" \
                        -db "${DATA}/${genus}/${org}.pep.all.fasta" \
                        -out "$tmpfile"
                    python3 -m biotp rename_header \
                        "$tmpfile" \
                        "$tmpfile" \
                        "${pnm_li[i]}" \
                        "" \
                        ""
                    cat "$tmpfile" >> "${taskdir}/${genus}.${symbol}.pep.fasta"
                    rm "$tmpfile"
                done
            done
        done
    }

    function main() {
        parse_args "$@"
        redeclare_genome_by_genus "$genus"
        make_dir
        retrieve_blastp
    }

    main "$@"
}

function addn_construct_pep_blastp_genus_symbol() {
    function usage() {
        cat <<EOS
Usage:  addn_construct_pep_blastp_genus_symbol <arg1> <arg2> (<arg3> <arg4> <arg5> ...) <arg1> <arg2> <arg3> (<arg4> <arg5> ...)

    arg1: cons_genus
    arg2: num [number of add_genus]
    arg3: add1_genus
    arg4: add2_genus
    arg5: add3_genus
    ...
    arg1: symbol
    arg2: symbol_org
    arg3: evalue
    arg4: pid
    arg5: pnm
    ...

EOS
        exit 1
    } 

    function parse_args() {
        if [[ $# -lt 9 ]]; then
            usage
        fi
        cons_genus="$1"
        num=$2

        addn_genus=()
        for ((i=0; i<num; i++)); do
            addn_genus[i]="${@:3+i:1}"
        done

        symbol="${@:3+num:1}"
        symbol_org="${@:4+num:1}"
        symbol_org=${symbol_org/ /_}
        evalue="${@:5+num:1}"

        local _pli=("${@:6+num}")
        pid_li=()
        pnm_li=()
        for ((i=0; i<${#_pli[@]}; i+=2)); do
            pid_li+=("${_pli[i]}")
            pnm_li+=("${_pli[i+1]}")
        done
    }

    function make_dir() {
        taskdir=$(make_dir_by_date $TASKFILE)
        echo "taskdir: ${taskdir}"
    }

    function retrieve_blastp_genus_symbol() {
        touch "${taskdir}/${cons_genus}.${symbol}.pep.fasta"
        echo "Protein BLAST"
        for org in ${org_li[*]}; do
            echo "ref: all [${org}], qry: ${symbol} [${symbol_org}]"
            local _blastp=$(blastp -outfmt 6 -evalue $evalue -db "${DATA}/${cons_genus}/${org}.pep.all.fasta" -query "${DATA}/${symbol_org}/${symbol}.pep.fasta")
            local _ids=$(echo "$_blastp" | cut -f 2 | sort -u)
            if [ -z "$_ids" ]; then
                echo "No hit."
                continue
            fi
            echo "Hit."; echo "$_ids"
            for _id in $_ids; do
                for ((i=0; i<${#pid_li[@]}; ++i)); do
                    if [[ "$_id" != "${pid_li[i]}" ]]; then
                        continue
                    fi
                    tmpfile=$(mktemp)
                    blastdbcmd \
                        -entry "$_id" \
                        -db "${DATA}/${cons_genus}/${org}.pep.all.fasta" \
                        -out "$tmpfile"
                    python3 -m biotp rename_header \
                        "$tmpfile" \
                        "$tmpfile" \
                        "${pnm_li[i]}" \
                        "" \
                        ""
                    cat "$tmpfile" >> "${taskdir}/${cons_genus}.${symbol}.pep.fasta"
                    rm "$tmpfile"
                done
            done
        done
    }

    function unset_genome_by_genus() {
        unset -v acc_li org_li asm_li
    }

    function addn_genus_symbol() {
        for add_genus in ${addn_genus[*]}; do
            unset_genome_by_genus
            redeclare_genome_by_genus "$add_genus"
            echo "Protein BLAST with ${add_genus}"
            for org in ${org_li[*]}; do
                echo "ref: all [${org}], qry: ${symbol} [${symbol_org}]"
                local _blastp; _blastp=$(blastp -outfmt 6 -evalue $evalue -db "${DATA}/${add_genus}/${org}.pep.all.fasta" -query "${DATA}/${symbol_org}/${symbol}.pep.fasta")
                local _ids; _ids=$(echo "$_blastp" | cut -f 2 | sort -u)
                if [ -z "$_ids" ]; then
                    echo "No hit."
                    continue
                fi
                echo "Hit."; echo "$_ids"
                for _id in $_ids; do
                    for ((i=0; i<${#pid_li[@]}; ++i)); do
                        if [[ "$_id" != "${pid_li[i]}" ]]; then
                            continue
                        fi
                        tmpfile=$(mktemp)
                        blastdbcmd \
                            -entry "$_id" \
                            -db "${DATA}/${add_genus}/${org}.pep.all.fasta" \
                            -out "$tmpfile"
                        python3 -m biotp rename_header \
                            "$tmpfile" \
                            "$tmpfile" \
                            "${pnm_li[i]}" \
                            "" \
                            ""
                        cat "$tmpfile" >> ""
                        rm "$tmpfile"
                    done
                done
            done
        done
    }

    function main() {
        parse_args "$@"
        redeclare_genome_by_genus "$cons_genus"
        make_dir
        retrieve_blastp_genus_symbol
        addn_genus_symbol
    }

    main "$@"    
}

function addn_genus_to_mltree_genus_symbol() {
    function usage() {
        cat <<EOS
Usage:  addn_genus_to_mltree_genus_symbol <arg1> <arg2> (<arg3> <arg4> <arg5> ...) <arg1> <arg2> <arg3> (<arg4> <arg5> ...)

    arg1: cons_genus
    arg2: num [number of add_genus]
    arg3: add1_genus
    arg4: add2_genus
    arg5: add3_genus
    ...
    arg1: symbol
    arg2: symbol_org
    arg3: evalue
    arg4: pid
    arg5: pnm
    ...

EOS
        exit 1
    } 

    function parse_args() {
        if [[ $# -lt 9 ]]; then
            usage
        fi
        cons_genus="$1"
        num=$2

        addn_genus=()
        for ((i=0; i<num; i++)); do
            addn_genus[i]="${@:3+i:1}"
        done

        symbol="${@:3+num:1}"
        symbol_org="${@:4+num:1}"
        symbol_org=${symbol_org/ /_}
        evalue="${@:5+num:1}"

        local _pli=("${@:6+num}")
        pid_li=()
        pnm_li=()
        for ((i=0; i<${#_pli[@]}; i+=2)); do
            pid_li+=("${_pli[i]}")
            pnm_li+=("${_pli[i+1]}")
        done
    }

    function make_dir() {
        taskdir=$(make_dir_by_date $TASKFILE)
        echo "taskdir: ${taskdir}"
    }

    function retrieve_blastp_genus_symbol() {
        touch "${taskdir}/${cons_genus}.${symbol}.pep.fasta"
        echo "Protein BLAST"
        for org in ${org_li[*]}; do
            echo "ref: all [${org}], qry: ${symbol} [${symbol_org}]"
            local _blastp=$(blastp -outfmt 6 -evalue $evalue -db "${DATA}/${cons_genus}/${org}.pep.all.fasta" -query "${DATA}/${symbol_org}/${symbol}.pep.fasta")
            local _ids=$(echo "$_blastp" | cut -f 2 | sort -u)
            if [ -z "$_ids" ]; then
                echo "No hit."
                continue
            fi
            echo "Hit."; echo "$_ids"
            for _id in $_ids; do
                for ((i=0; i<${#pid_li[@]}; ++i)); do
                    if [[ "$_id" != "${pid_li[i]}" ]]; then
                        continue
                    fi
                    tmpfile=$(mktemp)
                    blastdbcmd \
                        -entry "$_id" \
                        -db "${DATA}/${cons_genus}/${org}.pep.all.fasta" \
                        -out "$tmpfile"
                    python3 -m biotp rename_header \
                        "$tmpfile" \
                        "$tmpfile" \
                        "${pnm_li[i]}" \
                        "" \
                        ""
                    cat "$tmpfile" >> "${taskdir}/${cons_genus}.${symbol}.pep.fasta"
                    rm "$tmpfile"
                done
            done
        done
    }

    function unset_genome_by_genus() {
        unset -v acc_li org_li asm_li
    }

    function addn_genus_symbol() {
        for add_genus in ${addn_genus[*]}; do
            unset_genome_by_genus
            redeclare_genome_by_genus "$add_genus"
            echo "Protein BLAST with ${add_genus}"
            for org in ${org_li[*]}; do
                echo "ref: all [${org}], qry: ${symbol} [${symbol_org}]"
                local _blastp; _blastp=$(blastp -outfmt 6 -evalue $evalue -db "${DATA}/${add_genus}/${org}.pep.all.fasta" -query "${DATA}/${symbol_org}/${symbol}.pep.fasta")
                local _ids; _ids=$(echo "$_blastp" | cut -f 2 | sort -u)
                if [ -z "$_ids" ]; then
                    echo "No hit."
                    continue
                fi
                echo "Hit."; echo "$_ids"
                for _id in $_ids; do
                    for ((i=0; i<${#pid_li[@]}; ++i)); do
                        if [[ "$_id" != "${pid_li[i]}" ]]; then
                            continue
                        fi
                        tmpfile=$(mktemp)
                        blastdbcmd \
                            -entry "$_id" \
                            -db "${DATA}/${add_genus}/${org}.pep.all.fasta" \
                            -out "$tmpfile"
                        python3 -m biotp rename_header \
                            "$tmpfile" \
                            "$tmpfile" \
                            "${pnm_li[i]}" \
                            "" \
                            ""
                        cat "$tmpfile" >> "${taskdir}/${cons_genus}.${symbol}.pep.fasta"
                        rm "$tmpfile"
                    done
                done
            done
        done
    }

    function construct_msa_blastp() {
        mafft "${taskdir}/${cons_genus}.${symbol}.pep.fasta" > "${taskdir}/${cons_genus}.${symbol}.pep.aln"
    }

    function estimate_mltree() {
        raxml-ng \
            --msa "${taskdir}/${cons_genus}.${symbol}.pep.aln" \
            --all \
            --model Blosum62 \
            --bs-trees 1000 \
            --threads 8 \
            --redo
    }

    function make_tree() {
        Rscript ${SCRIPT}/make_tree.R  \
            "${taskdir}/${cons_genus}.${symbol}.pep.aln.raxml.support" \
            "${taskdir}/${cons_genus}.${symbol}.pep.aln.raxml.support.png"
    }

    function main() {
        parse_args "$@"
        redeclare_genome_by_genus "$cons_genus"
        make_dir
        retrieve_blastp_genus_symbol
        addn_genus_symbol
        construct_msa_blastp
        estimate_mltree
        make_tree
    }

    main "$@"
}

function construct_dna_flanking_region_blastp_genus_symbol() {
    function usage() {
        cat <<EOS
Usage:  construct_dna_flanking_region_blastp_genus_symbol <arg1> <arg2> <arg3> <arg4> <arg5> (<arg6> <arg7> ...) 

    arg1: genus
    arg2: symbol
    arg3: symbol_org
    arg4: evalue
    arg5: bp
    arg6: pid
    arg7: pnm
    ...
EOS
        exit 1
    }

    function parse_args() {
        if [[ $# -lt 7 ]]; then
            usage
        fi
        genus="$1"
        symbol="$2"
        symbol_org="${3/ /_}"
        evalue=$4
        bp=$5
        local _pli=("${@:6}")
        for ((i=0; i<${#_pli[@]}; i+=2)); do
            pid_li+=("${_pli[i]}")
            pnm_li+=("${_pli[i+1]}")
        done
    }

    function make_dir() {
        taskdir=$(make_dir_by_date $TASKFILE)
        echo "taskdir: ${taskdir}"
    }

    function retrieve_blastp_genus_symbol() {
        touch "${taskdir}/${genus}.${symbol}.dna_flanking_region.fasta"
        echo "Protein BLAST"
        for org in ${org_li[*]}; do
            echo "ref: all [${org}], qry: ${symbol} [${symbol_org}]"
            local _blastp=$(blastp -outfmt 6 -evalue $evalue -db ${DATA}/${genus}/${org}.pep.all.fasta -query ${DATA}/${symbol_org}/${symbol}.pep.fasta)
            local _ids=$(echo "$_blastp" | cut -f 2 | sort -u)
            if [ -z "$_ids" ]; then
                echo "No hit."
                continue
            fi
            echo "Hit."; echo "$_ids"
            for _id in $_ids; do
                for ((i=0; i<${#pid_li[@]}; ++i)); do
                    if [[ "$_id" != "${pid_li[i]}" ]]; then
                        continue
                    fi
                    tmpfile=$(mktemp)
                    local _seq; local _strand; local _start; local _end
                    read _seq _strand _start _end <<< $(python3 -m biotp output_seqid_strand_locs_by_proid "${DATA}/${genus}/${org}.genome.gff" "$_id")
                    echo "seqid: ${_seq}, strand: ${_strand}, start: ${_start}, end: ${_end}"
                    blastdbcmd \
                        -entry "$_seq" \
                        -db "${DATA}/${genus}/${org}.dna.toplevel.fasta" \
                        -out "$tmpfile"
                    python3 -m biotp slice_seq_flanking_region \
                        "$tmpfile" \
                        "$tmpfile" \
                        "${pnm_li[i]}" \
                        "${pid_li[i]}" \
                        "${bp}bp flanking region [${org}]" \
                        "$_strand" \
                        "$_start" \
                        "$_end" \
                        "$bp"
                    cat "$tmpfile" >> "${taskdir}/${genus}.${symbol}.dna_flanking_region.fasta"
                    rm "$tmpfile"
                done
            done
        done
    }

    function main() {
        parse_args "$@"
        redeclare_genome_by_genus "$genus"
        make_dir
        retrieve_blastp_genus_symbol
    }

    main "$@"
}

function construct_dna_upstream_region_blastp_genus_symbol() {
    function usage() {
        cat <<EOS
Usage:  construct_dna_upstream_region_blastp_genus_symbol <arg1> <arg2> <arg3> <arg4> <arg5> (<arg6> <arg7> ...)

    arg1: genus
    arg2: symbol
    arg3: symbol_org
    arg4: evalue
    arg5: bp
    arg6: pid
    arg7: pnm
    ...
    
EOS
        exit 1
    }

    function parse_args() {
        if [[ $# -lt 7 ]]; then
            usage
        fi
        genus="$1"
        symbol="$2"
        symbol_org="${3/ /_}"
        evalue=$4
        bp=$5
        local _pli=("${@:6}")
        for ((i=0; i<${#_pli[@]}; i+=2)); do
            pid_li+=("${_pli[i]}")
            pnm_li+=("${_pli[i+1]}")
        done
    }

    function make_dir() {
        taskdir=$(make_dir_by_date $TASKFILE)
        echo "taskdir: ${taskdir}"
    }

    function retrieve_blastp_genus_symbol() {
        touch "${taskdir}/${genus}.${symbol}.dna_upstream_region.fasta"
        echo "Protein BLAST"
        for org in ${org_li[*]}; do
            echo "ref: all [${org}], qry: ${symbol} [${symbol_org}]"
            local _blastp=$(blastp -outfmt 6 -evalue $evalue -db ${DATA}/${genus}/${org}.pep.all.fasta -query ${DATA}/${symbol_org}/${symbol}.pep.fasta)
            local _ids=$(echo "$_blastp" | cut -f 2 | sort -u)
            if [ -z "$_ids" ]; then
                echo "No hit."
                continue
            fi
            echo "Hit."; echo "$_ids"
            for _id in $_ids; do
                for ((i=0; i<${#pid_li[@]}; ++i)); do
                    if [[ "$_id" != "${pid_li[i]}" ]]; then
                        continue
                    fi
                    tmpfile=$(mktemp)
                    local _seq; local _strand; local _start; local _end
                    read _seq _strand _start _end <<< $(python3 -m biotp output_seqid_strand_locs_by_proid "${DATA}/${genus}/${org}.genome.gff" "$_id")
                    echo "seqid: ${_seq}, strand: ${_strand}, start: ${_start}, end: ${_end}"
                    blastdbcmd \
                        -entry "$_seq" \
                        -db "${DATA}/${genus}/${org}.dna.toplevel.fasta" \
                        -out "$tmpfile"
                    python3 -m biotp slice_seq_upstream_region \
                        "$tmpfile" \
                        "$tmpfile" \
                        "${pnm_li[i]}" \
                        "${pid_li[i]}" \
                        "${bp}bp upstream region [${org}]" \
                        "$_strand" \
                        "$_start" \
                        "$_end" \
                        "$bp"
                    cat "$tmpfile" >> "${taskdir}/${genus}.${symbol}.dna_upstream_region.fasta"
                    rm "$tmpfile"
                done
            done
        done
    }

    function main() {
        parse_args "$@"
        redeclare_genome_by_genus "$genus"
        make_dir
        retrieve_blastp_genus_symbol
    }

    main "$@"
}