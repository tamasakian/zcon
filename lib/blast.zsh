#!/usr/bin/env zsh

#### Function libs with BLAST

function makeblastdb_genome_by_genus() {
    ## This function is depend on fetch_genome_by_genus
    ## which is located in datasets.zsh
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

function remakeblastdb_genome_by_genus() {
    function usage() {
        cat <<EOS
Usage:  remakeblastdb_genome_by_genus <arg1>

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
        redeclare_genome_by_genus "$@"
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

    function blastp_genus_symbol() {
        touch "${taskdir}/${genus}.${symbol}.blastp"
        for org in ${org_li[*]}; do
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
        make_taskdir
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

    arg5: protein_id
    arg6: protein_name
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

        typeset -g -A pep_li; local _pli=("${@:5}")
        for ((i=1; i<${#_pli[@]}; i+=2)); do
            pep_li[${_pli[i]}]="${_pli[i+1]}"
                ## key: protein_id
                ## value: protein_name
        done
    }

    function retrieve_blastp() {
        touch "${taskdir}/${genus}.${symbol}.pep.fasta"
        for org in ${org_li[*]}; do
            local _blastp; _blastp=$(blastp -outfmt 6 -evalue $evalue -db ${DATA}/${genus}/${org}.pep.all.fasta -query ${DATA}/${symbol_org}/${symbol}.pep.fasta)
            local _ids; _ids=$(echo "$_blastp" | cut -f 2 | sort -u)
            if [ -z "$_ids" ]; then
                continue
            fi
            for _id in ${(f)_ids}; do
                for pep_id in ${(k)pep_li}; do
                    if [[ "$_id" != "$pep_id" ]]; then
                        continue
                    fi
                    tmpfile=$(mktemp)
                    blastdbcmd \
                        -entry "$pep_id" \
                        -db "${DATA}/${genus}/${org}.pep.all.fasta" \
                        -out "$tmpfile"
                    python3 -m biotp rename_header \
                        "$tmpfile" \
                        "$tmpfile" \
                        "${pep_li[$pep_id]}" \
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
        make_taskdir
        redeclare_genome_by_genus "$genus"
        retrieve_blastp
    }

    main "$@"
}

function addn_construct_pep_blastp_genus_symbol() {
    function usage() {
        cat <<EOS
Usage:  addn_construct_pep_blastp_genus_symbol <arg1> <arg2> (<arg3> <arg4> ... <arg(n+2)>) <arg(n+3)> <arg(n+4)> <arg(n+5)> (<arg(n+6)> <arg(n+7)> ...)

    arg1: cons_genus
    arg2: num [number of add_genus]
    arg3: add1_genus
    arg4: add2_genus
    arg(n+2): addn_genus

    arg(n+3): symbol
    arg(n+4): symbol_org
    arg(n+5): evalue

    arg(n+6): protein_id
    arg(n+7): protein_name
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
        for ((i=1; i<=num; i++)); do
            addn_genus[i]="${@:2+i:1}"
        done

        symbol="${@:3+num:1}"
        symbol_org="${@:4+num:1}"
        symbol_org=${symbol_org/ /_}
        evalue="${@:5+num:1}"

        typeset -g -A pep_li; local _pli=("${@:6+num}")
        for ((i=1; i<${#_pli[@]}; i+=2)); do
            pep_li[${_pli[i]}]="${_pli[i+1]}"
        done
    }

    function retrieve_blastp_genus_symbol() {
        touch "${taskdir}/${cons_genus}.${symbol}.pep.fasta"
        for org in ${org_li[*]}; do
            local _blastp; _blastp=$(blastp -outfmt 6 -evalue $evalue -db "${DATA}/${cons_genus}/${org}.pep.all.fasta" -query "${DATA}/${symbol_org}/${symbol}.pep.fasta")
            local _ids; _ids=$(echo "$_blastp" | cut -f 2 | sort -u)
            if [ -z "$_ids" ]; then
                continue
            fi
            for _id in ${(f)_ids}; do
                for pep_id in ${(k)pep_li}; do
                    if [[ "$_id" != "$pep_id" ]]; then
                        continue
                    fi
                    tmpfile=$(mktemp)
                    blastdbcmd \
                        -entry "$pep_id" \
                        -db "${DATA}/${cons_genus}/${org}.pep.all.fasta" \
                        -out "$tmpfile"
                    python3 -m biotp rename_header \
                        "$tmpfile" \
                        "$tmpfile" \
                        "${pep_li[$pep_id]}" \
                        "" \
                        ""
                    cat "$tmpfile" >> "${taskdir}/${cons_genus}.${symbol}.pep.fasta"
                    rm "$tmpfile"
                done
            done
        done
    }

    function addn_genus_symbol() {
        for add_genus in ${addn_genus[*]}; do
            unset -v acc_li org_li asm_li

            redeclare_genome_by_genus "$add_genus"
            for org in ${org_li[*]}; do
                local _blastp; _blastp=$(blastp -outfmt 6 -evalue $evalue -db "${DATA}/${add_genus}/${org}.pep.all.fasta" -query "${DATA}/${symbol_org}/${symbol}.pep.fasta")
                local _ids; _ids=$(echo "$_blastp" | cut -f 2 | sort -u)
                if [ -z "$_ids" ]; then
                    continue
                fi
                for _id in ${(f)_ids}; do
                    for pep_id in ${(k)pep_li}; do
                        if [[ "$_id" != "$pep_id" ]]; then
                            continue
                        fi
                        tmpfile=$(mktemp)
                        blastdbcmd \
                            -entry "$pep_id" \
                            -db "${DATA}/${add_genus}/${org}.pep.all.fasta" \
                            -out "$tmpfile"
                        python3 -m biotp rename_header \
                            "$tmpfile" \
                            "$tmpfile" \
                            "${pep_li[$pep_id]}" \
                            "" \
                            ""
                        cat "$tmpfile" >> "${taskdir}/${cons_genus}.${symbol}.pep.fasta"
                        rm "$tmpfile"
                    done
                done
            done
        done
    }

    function main() {
        parse_args "$@"
        make_taskdir
        redeclare_genome_by_genus "$cons_genus"
        retrieve_blastp_genus_symbol
        addn_genus_symbol
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

    arg6: protein_id
    arg7: protein_name
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

        typeset -g -A pep_li; local _pli=("${@:6}")
        for ((i=1; i<${#_pli[@]}; i+=2)); do
            pep_li[${_pli[i]}]="${_pli[i+1]}"
        done
    }

    function retrieve_blastp_genus_symbol() {
        touch "${taskdir}/${genus}.${symbol}.dna_flanking_region.fasta"
        for org in ${org_li[*]}; do
            local _blastp; _blastp=$(blastp -outfmt 6 -evalue $evalue -db ${DATA}/${genus}/${org}.pep.all.fasta -query ${DATA}/${symbol_org}/${symbol}.pep.fasta)
            local _ids; _ids=$(echo "$_blastp" | cut -f 2 | sort -u)
            if [ -z "$_ids" ]; then
                continue
            fi
            for _id in ${(f)_ids}; do
                for pep_id in ${(k)pep_li}; do
                    if [[ "$_id" != "$pep_id" ]]; then
                        continue
                    fi
                    tmpfile=$(mktemp)
                    local _seq; local _strand; local _start; local _end
                    read _seq _strand _start _end <<< $(python3 -m biotp output_seqid_strand_locs_by_proid "${DATA}/${genus}/${org}.genome.gff" "$pep_id")
                    blastdbcmd \
                        -entry "$_seq" \
                        -db "${DATA}/${genus}/${org}.dna.toplevel.fasta" \
                        -out "$tmpfile"
                    python3 -m biotp slice_seq_flanking_region \
                        "$tmpfile" \
                        "$tmpfile" \
                        "${pep_li[$pep_id]}" \
                        "$pep_id" \
                        "${bp}bp-flanking-region" \
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
        make_taskdir
        redeclare_genome_by_genus "$genus"
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

    arg6: protein_id
    arg7: protein_name
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

        typeset -g -A pep_li; local _pli=("${@:6}")
        for ((i=1; i<${#_pli[@]}; i+=2)); do
            pep_li[${_pli[i]}]="${_pli[i+1]}"
        done
    }

    function retrieve_blastp_genus_symbol() {
        touch "${taskdir}/${genus}.${symbol}.dna_upstream_region.fasta"
        for org in ${org_li[*]}; do
            local _blastp; _blastp=$(blastp -outfmt 6 -evalue $evalue -db ${DATA}/${genus}/${org}.pep.all.fasta -query ${DATA}/${symbol_org}/${symbol}.pep.fasta)
            local _ids; _ids=$(echo "$_blastp" | cut -f 2 | sort -u)
            if [ -z "$_ids" ]; then
                continue
            fi
            for _id in ${(f)_ids}; do
                for pep_id in ${(k)pep_li}; do
                    if [[ "$_id" != "$pep_id" ]]; then
                        continue
                    fi
                    tmpfile=$(mktemp)
                    local _seq; local _strand; local _start; local _end
                    read _seq _strand _start _end <<< $(python3 -m biotp output_seqid_strand_locs_by_proid "${DATA}/${genus}/${org}.genome.gff" "$pep_id")
                    blastdbcmd \
                        -entry "$_seq" \
                        -db "${DATA}/${genus}/${org}.dna.toplevel.fasta" \
                        -out "$tmpfile"
                    python3 -m biotp slice_seq_upstream_region \
                        "$tmpfile" \
                        "$tmpfile" \
                        "${pep_li[$pep_id]}" \
                        "$pep_id" \
                        "${bp}bp-upstream-region" \
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
        make_taskdir
        redeclare_genome_by_genus "$genus"
        retrieve_blastp_genus_symbol
    }

    main "$@"
}