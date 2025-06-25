#!/usr/bin/env zsh
# Last updated: 2025-06-26
# Tools: HMMER3 SonicParanoid2 MCScanX bithon

: << 'FUNCTIONS'
run_pfam: Run Pfam domain search on a given FASTA file.
run_sonicparanoid: Run SonicParanoid2 on a given directory of FASTA files.
run_mcscanx: Run MCScanX for collinearity analysis between pairs of species.
FUNCTIONS

function run_pfam() {
    if [[ $# -lt 2 ]]; then
        echo "Usage: run_pfam <sp_name1> <sp_name2> ..." >&2
        return 1
    fi

    local sp_names=("${@}")
    local threads=4
    local taskdir="${TASKFILE}/pfam_$(date +"%Y-%m-%d-%H-%M-%S")"
    mkdir -p "${taskdir}/tmp" "${taskdir}/fasta" "${taskdir}/fasta_prefix" "${taskdir}/input" "${taskdir}/output"

    ## --- FASTA Preparation ---
    for sp in "${sp_names[@]}"; do
        local sp_fs="${sp// /_}"
        local gn="${sp_fs%%_*}"
        local pep="${DATA}/${gn}/${sp_fs}.pep.all.fasta"
        local cds="${DATA}/${gn}/${sp_fs}.cds.all.fasta"

        if [[ -e "${DATA}/Ensembl/${sp_fs}.pep.all.fasta" ]]; then
            bithon ensgls -i "${DATA}/Ensembl/${sp_fs}.pep.all.fasta" -o "${taskdir}/fasta/${sp_fs}.fasta" --header transcript
        else
            mkdir -p "${taskdir}/tmp/${sp_fs}"
            cp "${pep}" "${taskdir}/tmp/${sp_fs}/pep.fasta"
            cp "${cds}" "${taskdir}/tmp/${sp_fs}/cds.fasta"
            bithon gls -i "${taskdir}/tmp/${sp_fs}" -o "${taskdir}/tmp/${sp_fs}"
            cp "${taskdir}/tmp/${sp_fs}/longest.pep.fasta" "${taskdir}/fasta/${sp_fs}.fasta"
            rm -r "${taskdir}/tmp/${sp_fs}"
        fi

        python3 -m fasp prefix_to_sequence_ids \
            "${taskdir}/fasta/${sp_fs}.fasta" \
            "${taskdir}/fasta_prefix/${sp_fs}.fasta" \
            "${sp_fs}"
    done
    cat "${taskdir}/fasta_prefix/"*.fasta > "${taskdir}/input/input.fasta"

    ## --- Pfam Domain Search ---
    hmmscan \
        -o "${taskdir}/output/hmmscan.txt" \
        --domtblout "${taskdir}/output/domtblout.txt" \
        --cpu "$threads" \
        "${DATA}/Pfam/Pfam-A.hmm" \
        "${taskdir}/input/input.fasta"
}

run_sonicparanoid(){
    if [[ $# -lt 2 ]]; then
        echo "Usage: run_sonicparanoid <sp_name1> <sp_name2> ..." >&2
        return 1
    fi

    local sp_names=("${@}")
    local threads=4
    local taskdir="${TASKFILE}/sonicparanoid_$(date +"%Y-%m-%d-%H-%M-%S")"
    mkdir -p "${taskdir}/input" "${taskdir}/output"

    ## --- FASTA Preparation ---
    for sp in "${sp_names[@]}"; do
        local sp_fs="${sp// /_}"
        local gn="${sp_fs%%_*}"
        local pep="${DATA}/${gn}/${sp_fs}.pep.all.fasta"
        local cds="${DATA}/${gn}/${sp_fs}.cds.all.fasta"

        if [[ -e "${DATA}/Ensembl/${sp_fs}.pep.all.fasta" ]]; then
            bithon ensgls -i "${DATA}/Ensembl/${sp_fs}.pep.all.fasta" -o "${taskdir}/input/${sp_fs}.fasta" --header transcript
        else
            mkdir -p "${taskdir}/input/${sp_fs}"
            cp "${pep}" "${taskdir}/input/${sp_fs}/pep.fasta"
            cp "${cds}" "${taskdir}/input/${sp_fs}/cds.fasta"
            bithon gls -i "${taskdir}/input/${sp_fs}" -o "${taskdir}/input/${sp_fs}"
            cp "${taskdir}/input/${sp_fs}/longest.pep.fasta" "${taskdir}/input/${sp_fs}.fasta"
            rm -r "${taskdir}/input/${sp_fs}"
        fi

        python3 -m fasp prefix_to_sequence_ids \
            "${taskdir}/input/${sp_fs}.fasta" \
            "${taskdir}/input/${sp_fs}.fasta" \
            "${sp_fs}"
    done

    ## --- SonicParanoid2 ---
    sonicparanoid -i "${taskdir}/input" -o "${taskdir}/output" -t "$threads"
}

function run_mcscanx() {
    if [[ $# -lt 2 ]]; then
        echo "Usage: run_mcscanx <sp_name1> <sp_name2> ..." >&2
        return 1
    fi

    local sp_names=("${@}")
    local threads=4
    local taskdir="${TASKFILE}/mcscanx_$(date +"%Y-%m-%d-%H-%M-%S")"
    mkdir -p "${taskdir}/tmp" "${taskdir}/fasta" "${taskdir}/fasta_prefix" "${taskdir}/gff" "${taskdir}/bed" "${taskdir}/bed_prefix" "${taskdir}/input" "${taskdir}/output"

    # --- FASTA & GFF Preparation ---
    for sp in "${sp_names[@]}"; do
        local sp_fs="${sp// /_}"
        local gn="${sp_fs%%_*}"

        local pep="${DATA}/${gn}/${sp_fs}.pep.all.fasta"
        local cds="${DATA}/${gn}/${sp_fs}.cds.all.fasta"
        local gff="${DATA}/${gn}/${sp_fs}.genome.gff"

        if [[ -e "${DATA}/Ensembl/${sp_fs}.pep.all.fasta" && -e "${DATA}/Ensembl/${sp_fs}.genome.gff" ]]; then
            bithon ensgls -i "${DATA}/Ensembl/${sp_fs}.pep.all.fasta" -o "${taskdir}/fasta/${sp_fs}.fasta" --header transcript
            cp "${DATA}/Ensembl/${sp_fs}.genome.gff" "${taskdir}/gff/${sp_fs}.gff"
        else
            mkdir -p "${taskdir}/tmp/${sp_fs}"
            cp "${pep}" "${taskdir}/tmp/${sp_fs}/pep.fasta"
            cp "${cds}" "${taskdir}/tmp/${sp_fs}/cds.fasta"
            cp "${gff}" "${taskdir}/gff/${sp_fs}.gff"
            bithon gls -i "${taskdir}/tmp/${sp_fs}" -o "${taskdir}/tmp/${sp_fs}"
            cp "${taskdir}/tmp/${sp_fs}/longest.pep.fasta" "${taskdir}/fasta/${sp_fs}.fasta"
            rm -r "${taskdir}/tmp/${sp_fs}"
        fi

        python3 -m g2bp fasta4mcscanx "${taskdir}/gff/${sp_fs}.gff" "${taskdir}/bed/${sp_fs}.bed" "${taskdir}/fasta/${sp_fs}.fasta"
        python3 -m fasp prefix_to_sequence_ids "${taskdir}/fasta/${sp_fs}.fasta" "${taskdir}/fasta_prefix/${sp_fs}.fasta" "${sp_fs}"
        awk -v sp="${sp_fs}_" 'BEGIN{OFS=FS="\t"} {$4=sp $4; print}' "${taskdir}/bed/${sp_fs}.bed" > "${taskdir}/bed_prefix/${sp_fs}.bed"

    done

    for ((i=1; i<=${#sp_names[@]}-1; i++)); do
        for ((j=i+1; j<=${#sp_names[@]}; j++)); do
            local sp1="${sp_names[i]}"
            local sp2="${sp_names[j]}"
            local sp1_fs="${sp1// /_}"
            local sp2_fs="${sp2// /_}"
            local gn1="${sp1_fs%%_*}"
            local gn2="${sp2_fs%%_*}"
            local pairdir="${taskdir}/${sp1_fs}__${sp2_fs}"
            mkdir -p "$pairdir"

            # --- Blastp and MCScanX ---
            cp "${taskdir}/fasta_prefix/${sp1_fs}.fasta" "${pairdir}/${sp1_fs}.fasta"
            cp "${taskdir}/fasta_prefix/${sp2_fs}.fasta" "${pairdir}/${sp2_fs}.fasta"
            cat "${taskdir}/bed_prefix/${sp1_fs}.bed" "${taskdir}/bed_prefix/${sp2_fs}.bed" > "${pairdir}/${sp1_fs}__${sp2_fs}.gff"
            diamond makedb --in "${pairdir}/${sp2_fs}.fasta" -d "${pairdir}/${sp2_fs}.dmnd"
            diamond blastp \
                -q "${pairdir}/${sp1_fs}.fasta" \
                -d "${pairdir}/${sp2_fs}.dmnd" \
                -o "${pairdir}/${sp1_fs}__${sp2_fs}.blast" \
                --more-sensitive \
                -e 1e-5 \
                --threads "$threads" \
                --outfmt 6
            MCScanX "${pairdir}/${sp1_fs}__${sp2_fs}" \
                -b 2 \
                -e 1e-5
        done
    done
}

