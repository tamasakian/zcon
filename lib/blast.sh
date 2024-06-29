#!/bin/zsh

#### Function libs with BLAST

function makeblastdb_genome_by_genus() {
    ## This function is depend on fetch_genome_by_genus, datasets.sh
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
        > "${taskdir}/${genus}.${symbol}.blastp"
        echo "Protein BLAST"
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
Usage:  construct_pep_blastp_genus_symbol <arg1> <arg2> <arg3> (<arg4> <arg5> ...)

    arg1: genus
    arg2: symbol
    arg3: symbol_org
    arg4: pid
    arg5: pnm
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
        for ((i=0; i<${#_pli[@]}; i+=2)); do
            pid_li+=("${_pli[i]}")
            pnm_li+=("${_pli[i+1]}")
        done
    }

    function make_dir() {
        taskdir=$(make_dir_by_date $TASKFILE)
        echo "taskdir: ${taskdir}"
    }

    function retrieve_blastp() {
        > "${taskdir}/${genus}.${symbol}.pep.fasta"
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
                    blastdbcmd \
                        -entry "$_id" \
                        -db "${DATA}/${genus}/${org}.pep.all.fasta" \
                        -out "$tmpfile"
                    $PYTHON3 -m biotp rename_header \
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

function make_mltree_genus_symbol() {
    function usage() {
        cat <<EOS
Usage:  make_mltree_genus_symbol <arg1> <arg2> <arg3> (<arg4> <arg5> ...)

    arg1: genus
    arg2: symbol
    arg3: symbol_org
    arg4: pid
    arg5: pnm
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
        > "${taskdir}/${genus}.${symbol}.pep.fasta"
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
                    blastdbcmd \
                        -entry "$_id" \
                        -db "${DATA}/${genus}/${org}.pep.all.fasta" \
                        -out "$tmpfile"
                    $PYTHON3 -m biotp rename_header \
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

    function construct_msa_blastp() {
        mafft "${taskdir}/${genus}.${symbol}.pep.fasta" > "${taskdir}/${genus}.${symbol}.pep.aln"
    }

    function estimate_mltree() {
        raxml-ng \
            --msa "${taskdir}/${genus}.${symbol}.pep.aln" \
            --all \
            --model Blosum62 \
            --bs-trees 500 \
            --threads 8 \
            --redo
    }

    function make_tree() {
        Rscript ${SCRIPT}/make_tree.R  \
            "${taskdir}/${genus}.${symbol}.pep.aln.raxml.support" \
            "${taskdir}/${genus}.${symbol}.pep.aln.raxml.support.png"
    }

    function main() {
        parse_args "$@"
        redeclare_genome_by_genus "$genus"
        make_dir
        retrieve_blastp_genus_symbol
        construct_msa_blastp
        estimate_mltree
        make_tree
    }

    main "$@"
}

function add1_genus_to_mltree_genus_symbol() {
    function usage() {
        cat <<EOS
Usage:  add1_genus_to_mltree_genus_symbol <arg1> <arg2> <arg3> <arg4> <arg5> (<arg6> <arg7> ...)

    arg1: cons_genus
    arg2: add1_genus
    arg3: symbol
    arg4: symbol_org
    arg5: evalue
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
        cons_genus="$1"
        add1_genus="$2"
        symbol="$3"
        symbol_org="${4/ /_}"
        evalue=$5
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
        > "${taskdir}/${cons_genus}.${symbol}.pep.fasta"
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
                    $PYTHON3 -m biotp rename_header \
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

    function add1_genus_symbol() {
        echo "Protein BLAST"
        for org in ${org_li[*]}; do
            echo "ref: all [${org}], qry: ${symbol} [${symbol_org}]"
            local _blastp; _blastp=$(blastp -outfmt 6 -evalue $evalue -db "${DATA}/${add1_genus}/${org}.pep.all.fasta" -query "${DATA}/${symbol_org}/${symbol}.pep.fasta")
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
                        -db "${DATA}/${add1_genus}/${org}.pep.all.fasta" \
                        -out "$tmpfile"
                    $PYTHON3 -m biotp rename_header \
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

    function construct_msa_blastp() {
        mafft "${taskdir}/${cons_genus}.${symbol}.pep.fasta" > "${taskdir}/${cons_genus}.${symbol}.pep.aln"
    }

    function estimate_mltree() {
        raxml-ng \
            --msa "${taskdir}/${cons_genus}.${symbol}.pep.aln" \
            --all \
            --model Blosum62 \
            --bs-trees 500 \
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

        ## add1
        unset_genome_by_genus
        redeclare_genome_by_genus "$add1_genus"
        add1_genus_symbol

        construct_msa_blastp
        estimate_mltree
        make_tree
    }

    main "$@"
}

function add2_genus_to_mltree_genus_symbol() {
    function usage() {
        cat <<EOS
Usage:  add2_genus_to_mltree_genus_symbol <arg1> <arg2> <arg3> <arg4> <arg5> <arg6> (<arg7> <arg8> ...)

    arg1: cons_genus
    arg2: add1_genus
    arg3: add2_genus
    arg4: symbol
    arg5: symbol_org
    arg6: evalue
    arg7: pid
    arg8: pnm
    ...

EOS
        exit 1
    } 

    function parse_args() {
        if [[ $# -lt 8 ]]; then
            usage
        fi
        cons_genus="$1"
        add1_genus="$2"
        add2_genus="$3"
        symbol="$4"
        symbol_org="${5/ /_}"
        evalue=$6
        local _pli=("${@:7}")
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
        > "${taskdir}/${cons_genus}.${symbol}.pep.fasta"
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
                    $PYTHON3 -m biotp rename_header \
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

    function add1_genus_symbol() {
        echo "Protein BLAST"
        for org in ${org_li[*]}; do
            echo "ref: all [${org}], qry: ${symbol} [${symbol_org}]"
            local _blastp; _blastp=$(blastp -outfmt 6 -evalue $evalue -db "${DATA}/${add1_genus}/${org}.pep.all.fasta" -query "${DATA}/${symbol_org}/${symbol}.pep.fasta")
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
                        -db "${DATA}/${add1_genus}/${org}.pep.all.fasta" \
                        -out "$tmpfile"
                    $PYTHON3 -m biotp rename_header \
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

    function add2_genus_symbol() {
        echo "Protein BLAST"
        for org in ${org_li[*]}; do
            echo "ref: all [${org}], qry: ${symbol} [${symbol_org}]"
            local _blastp; _blastp=$(blastp -outfmt 6 -evalue $evalue -db "${DATA}/${add2_genus}/${org}.pep.all.fasta" -query "${DATA}/${symbol_org}/${symbol}.pep.fasta")
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
                        -db "${DATA}/${add2_genus}/${org}.pep.all.fasta" \
                        -out "$tmpfile"
                    $PYTHON3 -m biotp rename_header \
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

    function construct_msa_blastp() {
        mafft "${taskdir}/${cons_genus}.${symbol}.pep.fasta" > "${taskdir}/${cons_genus}.${symbol}.pep.aln"
    }

    function estimate_mltree() {
        raxml-ng \
            --msa "${taskdir}/${cons_genus}.${symbol}.pep.aln" \
            --all \
            --model Blosum62 \
            --bs-trees 500 \
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

        ## add1
        unset_genome_by_genus
        redeclare_genome_by_genus "$add1_genus"
        add1_genus_symbol

        ## add2
        unset_genome_by_genus
        redeclare_genome_by_genus "$add2_genus"
        add2_genus_symbol

        construct_msa_blastp
        estimate_mltree
        make_tree
    }

    main "$@"
}

function add3_genus_to_mltree_genus_symbol() {
    function usage() {
        cat <<EOS
Usage:  add3_genus_to_mltree_genus_symbol <arg1> <arg2> <arg3> <arg4> <arg5> <arg6> <arg7> (<arg8> <arg9> ...)

    arg1: cons_genus
    arg2: add1_genus
    arg3: add2_genus
    arg4: add3_genus
    arg5: symbol
    arg6: symbol_org
    arg7: evalue
    arg8: pid
    arg9: pnm
    ...

EOS
        exit 1
    } 

    function parse_args() {
        if [[ $# -lt 9 ]]; then
            usage
        fi
        cons_genus="$1"
        add1_genus="$2"
        add2_genus="$3"
        add3_genus="$4"
        symbol="$5"
        symbol_org="${6/ /_}"
        evalue=$7
        local _pli=("${@:8}")
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
        > "${taskdir}/${cons_genus}.${symbol}.pep.fasta"
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
                    $PYTHON3 -m biotp rename_header \
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

    function add1_genus_symbol() {
        echo "Protein BLAST"
        for org in ${org_li[*]}; do
            echo "ref: all [${org}], qry: ${symbol} [${symbol_org}]"
            local _blastp; _blastp=$(blastp -outfmt 6 -evalue $evalue -db "${DATA}/${add1_genus}/${org}.pep.all.fasta" -query "${DATA}/${symbol_org}/${symbol}.pep.fasta")
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
                        -db "${DATA}/${add1_genus}/${org}.pep.all.fasta" \
                        -out "$tmpfile"
                    $PYTHON3 -m biotp rename_header \
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

    function add2_genus_symbol() {
        echo "Protein BLAST"
        for org in ${org_li[*]}; do
            echo "ref: all [${org}], qry: ${symbol} [${symbol_org}]"
            local _blastp; _blastp=$(blastp -outfmt 6 -evalue $evalue -db "${DATA}/${add2_genus}/${org}.pep.all.fasta" -query "${DATA}/${symbol_org}/${symbol}.pep.fasta")
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
                        -db "${DATA}/${add2_genus}/${org}.pep.all.fasta" \
                        -out "$tmpfile"
                    $PYTHON3 -m biotp rename_header \
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

    function add3_genus_symbol() {
        echo "Protein BLAST"
        for org in ${org_li[*]}; do
            echo "ref: all [${org}], qry: ${symbol} [${symbol_org}]"
            local _blastp; _blastp=$(blastp -outfmt 6 -evalue $evalue -db "${DATA}/${add3_genus}/${org}.pep.all.fasta" -query "${DATA}/${symbol_org}/${symbol}.pep.fasta")
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
                        -db "${DATA}/${add3_genus}/${org}.pep.all.fasta" \
                        -out "$tmpfile"
                    $PYTHON3 -m biotp rename_header \
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

    function construct_msa_blastp() {
        mafft "${taskdir}/${cons_genus}.${symbol}.pep.fasta" > "${taskdir}/${cons_genus}.${symbol}.pep.aln"
    }

    function estimate_mltree() {
        raxml-ng \
            --msa "${taskdir}/${cons_genus}.${symbol}.pep.aln" \
            --all \
            --model Blosum62 \
            --bs-trees 2000 \
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

        ## add1
        unset_genome_by_genus
        redeclare_genome_by_genus "$add1_genus"
        add1_genus_symbol

        ## add2
        unset_genome_by_genus
        redeclare_genome_by_genus "$add2_genus"
        add2_genus_symbol

        ## add3
        unset_genome_by_genus
        redeclare_genome_by_genus "$add3_genus"
        add3_genus_symbol

        construct_msa_blastp
        estimate_mltree
        make_tree
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
        > "${taskdir}/${cons_genus}.${symbol}.pep.fasta"
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
                    $PYTHON3 -m biotp rename_header \
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
                        $PYTHON3 -m biotp rename_header \
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

function make_njtree_genus_symbol() {
    function usage() {
        cat <<EOS
Usage:  make_njtree_genus_symbol <arg1> <arg2> <arg3> (<arg4> <arg5> ...) 

    arg1: genus
    arg2: symbol
    arg3: symbol_org
    arg4: pid
    arg5: pnm
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
        cp "${DATA}/${symbol_org}/${symbol}.pep.fasta" "${taskdir}/${genus}.${symbol}.pep.fasta"
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
                    blastdbcmd \
                        -entry "$_id" \
                        -db "${DATA}/${genus}/${org}.pep.all.fasta" \
                        -out "$tmpfile"
                    $PYTHON3 -m biotp rename_header \
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

    function construct_msa_blastp() {
        mafft --clustalout "${taskdir}/${genus}.${symbol}.pep.fasta" > "${taskdir}/${genus}.${symbol}.pep.clustal.aln"
    }

    function convert_msa_to_stockholm() {
        msaconverter \
            -i ${taskdir}/${genus}.${symbol}.pep.clustal.aln \
            -o ${taskdir}/${genus}.${symbol}.pep.sth \
            -p clustal \
            -q stockholm \
            -t protein
    }

    function estimate_njtree() {
        rapidnj \
            "${taskdir}/${genus}.${symbol}.pep.sth" \
            -i sth \
            -x "${taskdir}/${genus}.${symbol}.pep.sth.nwk" \
            -t p
    }

    function make_tree() {
        Rscript ${SCRIPT}/make_tree.R  \
            "${taskdir}/${genus}.${symbol}.pep.sth.nwk" \
            "${taskdir}/${genus}.${symbol}.pep.sth.png"
    }

    function main() {
        parse_args "$@"
        redeclare_genome_by_genus "$genus"
        make_dir
        retrieve_blastp_genus_symbol
        construct_msa_blastp
        convert_msa_to_stockholm
        estimate_njtree
        make_tree
    }

    main "$@"
}

function construct_dna_flanking_region_genus_symbol() {
    function usage() {
        cat <<EOS
Usage:  construct_dna_flanking_region <arg1> <arg2> <arg3> <arg4> <arg5> (<arg6> <arg7> ...) 

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
        > "${taskdir}/${genus}.${symbol}.dna_flanking_region.fasta"
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
                    read _seq _strand _start _end <<< $($PYTHON3 -m biotp output_seqid_strand_locs_by_proid "${DATA}/${genus}/${org}.genome.gff" "$_id")
                    echo "seqid: ${_seq}, strand: ${_strand}, start: ${_start}, end: ${_end}"
                    blastdbcmd \
                        -entry "$_seq" \
                        -db "${DATA}/${genus}/${org}.dna.toplevel.fasta" \
                        -out "$tmpfile"
                    $PYTHON3 -m biotp slice_seq_flanking_region \
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

function construct_dna_upstream_region_genus_symbol() {
    function usage() {
        cat <<EOS
Usage:  construct_dna_upstream_region <arg1> <arg2> <arg3> <arg4> <arg5> (<arg6> <arg7> ...)

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
        > "${taskdir}/${genus}.${symbol}.dna_upstream_region.fasta"
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
                    read _seq _strand _start _end <<< $($PYTHON3 -m biotp output_seqid_strand_locs_by_proid "${DATA}/${genus}/${org}.genome.gff" "$_id")
                    echo "seqid: ${_seq}, strand: ${_strand}, start: ${_start}, end: ${_end}"
                    blastdbcmd \
                        -entry "$_seq" \
                        -db "${DATA}/${genus}/${org}.dna.toplevel.fasta" \
                        -out "$tmpfile"
                    $PYTHON3 -m biotp slice_seq_upstream_region \
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

## Work in progress
function counstruct_dna_intron_genus_symbol() {
    function usage() {
        cat <<EOS
Usage:  construct_dna_intron_genus_symbol <arg1> <arg2> <arg3> <arg4> (<arg5> <arg6> ...)

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
        if [[ $# -lt 7 ]]; then
            usage
        fi
        genus="$1"
        symbol="$2"
        symbol_org="${3/ /_}"
        evalue=$4
        local _pli=("${@:5}")
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
        > "${taskdir}/${genus}.${symbol}.dna_intron.fasta"
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
                    read _seq _strand _start _end <<< $($PYTHON3 -m biotp output_seqid_strand_locs_by_proid "${DATA}/${genus}/${org}.genome.gff" "$_id")
                    echo "seqid: ${_seq}, strand: ${_strand}, start: ${_start}, end: ${_end}"
                    blastdbcmd \
                        -entry "$_seq" \
                        -db "${DATA}/${genus}/${org}.dna.toplevel.fasta" \
                        -out "$tmpfile"
                    $PYTHON3 -m biotp slice_seq_upstream_region \
                        "$tmpfile" \
                        "$tmpfile" \
                        "${pnm_li[i]}" \
                        "${pid_li[i]}" \
                        "${bp}bp flanking region [${org}]" \
                        "$_strand" \
                        "$_start" \
                        "$_end" \
                        "$bp"
                    cat "$tmpfile" >> "${taskdir}/${genus}.${symbol}.dna_intron.fasta"
                    rm "$tmpfile"
                done
            done
        done
    }
}


