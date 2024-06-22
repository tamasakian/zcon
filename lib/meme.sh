#!/bin/zsh

#### Function libs with MEME

function download_meme_databases() {

    function download_databases() {
        wget -P "$DOWNLOAD" "https://meme-suite.org/meme/meme-software/Databases/motifs/motif_databases.12.24.tgz"
        wget -P "$DOWNLOAD" "https://meme-suite.org/meme/meme-software/Databases/gomo/gomo_databases.3.2.tgz"
        wget -P "$DOWNLOAD" "https://meme-suite.org/meme/meme-software/Databases/tgene/tgene_databases.1.0.tgz"
    }

    function copy_databases() {
        tar -zxvf "${DOWNLOAD}/motif_databases.12.24.tgz" -C "$DATA"
        tar -zxvf "${DOWNLOAD}/gomo_databases.3.2.tgz"    -C "$DATA"
        tar -zxvf "${DOWNLOAD}/tgene_databases.1.0.tgz"   -C "$DATA"
    }

    function main() {
        download_databases
        copy_databases
    }
    
    main
}

function fimo_dna_upstream_region_genus_symbol() {
    function usage() {
        cat <<EOS
Usage:  fimo_dna_upstream_region <arg1> <arg2> <arg3> <arg4> <arg5> <arg6> <arg7> (<arg8> <arg9> ...)

    arg1: genus
    arg2: symbol
    arg3: symbol_org
    arg4: evalue (blastn)
    arg5: bp
    arg6: qvalue (fimo)
    arg7: database
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
        genus="$1"
        symbol="$2"
        symbol_org="${3/ /_}"
        evalue=$4
        bp=$5
        qvalue=$6
        database="$7"
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
                        "${bp}bp flanking region [${org}]" \
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

    function scan_motif() {
        fimo \
            --o "${taskdir}/fimo_out" \
            --thresh $qvalue \
            --qv-thresh \
            --verbosity 1 \
            "$database" \
            "${taskdir}/${genus}.${symbol}.dna_upstream_region.fasta"
    }
    
    function main() {
        parse_args "$@"
        redeclare_genome_by_genus "$genus"
        make_dir
        retrieve_blastp_genus_symbol
        scan_motif
    }

    main "$@"
}

