#!/usr/bin/env zsh
# Last updated: 2025-07-01
# Tools: HMMER3 SonicParanoid2 MCScanX bithon

: << 'FUNCTIONS'
run_pfam: Run Pfam domain search on a given FASTA file.
run_sonicparanoid: Run SonicParanoid2 on a given directory of FASTA files.
run_mcscanx: Run MCScanX for collinearity analysis between pairs of species.
run_mcscanx_da: Run MCScanX downstream analyses on collinearity files.
construct_family_tree_mso_spo: Construct a family tree from MSO and SPO files.
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
        --domtblout "${taskdir}/output/domtblout_ga.txt" \
        --cpu "$threads" \
        --cut_ga \
        "${DATA}/Pfam/Pfam-A.hmm" \
        "${taskdir}/input/input.fasta"

    hmmscan \
        -o "${taskdir}/output/hmmscan.txt" \
        --domtblout "${taskdir}/output/domtblout_nc.txt" \
        --cpu "$threads" \
        --cut_nc \
        "${DATA}/Pfam/Pfam-A.hmm" \
        "${taskdir}/input/input.fasta"

    hmmscan \
        -o "${taskdir}/output/hmmscan.txt" \
        --domtblout "${taskdir}/output/domtblout_tc.txt" \
        --cpu "$threads" \
        --cut_tc \
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
    mkdir -p "${taskdir}/fasta" "${taskdir}/input" "${taskdir}/output"

    ## --- FASTA Preparation ---
    for sp in "${sp_names[@]}"; do
        local sp_fs="${sp// /_}"
        local gn="${sp_fs%%_*}"
        local pep="${DATA}/${gn}/${sp_fs}.pep.all.fasta"
        local cds="${DATA}/${gn}/${sp_fs}.cds.all.fasta"

        if [[ -e "${DATA}/Ensembl/${sp_fs}.pep.all.fasta" ]]; then
            bithon ensgls -i "${DATA}/Ensembl/${sp_fs}.pep.all.fasta" -o "${taskdir}/fasta/${sp_fs}.fasta" --header transcript
        else
            mkdir -p "${taskdir}/fasta/${sp_fs}"
            cp "${pep}" "${taskdir}/fasta/${sp_fs}/pep.fasta"
            cp "${cds}" "${taskdir}/fasta/${sp_fs}/cds.fasta"
            bithon gls -i "${taskdir}/fasta/${sp_fs}" -o "${taskdir}/fasta/${sp_fs}"
            cp "${taskdir}/fasta/${sp_fs}/longest.pep.fasta" "${taskdir}/fasta/${sp_fs}.fasta"
            rm -r "${taskdir}/fasta/${sp_fs}"
        fi

        python3 -m fasp prefix_to_sequence_ids \
            "${taskdir}/fasta/${sp_fs}.fasta" \
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
        awk -v sp="${sp_fs}_" 'BEGIN{OFS=FS="\t"} {$2=sp $2; print}' "${taskdir}/bed/${sp_fs}.bed" > "${taskdir}/bed_prefix/${sp_fs}.bed"

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
            cp "${pairdir}/${sp1_fs}__${sp2_fs}.collinearity" "${taskdir}/output/${gn1}__${gn2}.collinearity"
        done
    done
}

function run_mcscanx_da() {
    if [[ $# != 1 ]]; then
        echo "Usage: run_mcscanx_da <taskname>" >&2
        return 1
    fi
    local taskname="$1"
    local taskdir="${TASKFILE}/${taskname}/mcscanx_da_$(date +"%Y-%m-%d-%H-%M-%S")"
    mkdir -p "${taskdir}/output"

    if [[ -d "${TASKFILE}/${taskname}/input" ]]; then
        for input_file in "${TASKFILE}/${taskname}/input/"*.collinearity; do
            if [[ -f "$input_file" ]]; then
                local base_name="${input_file##*/}"
                local output_file="${taskdir}/output/${base_name%.collinearity}.txt"
                perl ${HOME}/bin/MCScanX-1.0.0/downstream_analyses/group_syntenic_genes.pl -i "$input_file" -o "$output_file"
            fi
        done
    fi
}

function construct_family_tree_mso_spo(){
    if [[ $# -lt 3 ]]; then
        echo "Usage: construct_family_tree_mso_spo <taskname> <family> <PF_entry> ..." >&2
        return 1
    fi

    local taskname="$1"
    local family="$2"
    local PF_entry=("${@:3}")
    local taskdir="${TASKFILE}/${taskname}/family_${family}_$(date +"%Y-%m-%d-%H-%M-%S")"
    mkdir -p "${taskdir}/output"

    python3 -m biotp extract_protein_by_entry \
        "${TASKFILE}/${taskname}/input/domtblout.txt" \
        "${taskdir}/output/${family}.txt" \
        "${PF_entry[@]}"

    if [[ ! -s "${taskdir}/output/${family}.txt" ]]; then
        echo "No entries found for family ${family} with PF entries: ${PF_entry[*]}. Please check the input file or PF entries." >&2
        return 1
    fi

    python3 -m fasp seq_extractor \
        "${TASKFILE}/${taskname}/input/input.fasta" \
        "${taskdir}/output/${family}.fasta" \
        "${taskdir}/output/${family}.txt"

    gs2 -e 100 -l "${taskdir}/output/${family}.fasta" > "${taskdir}/output/${family}.tree"

    python3 -m biotp annotate_tree_mso_spo \
        "${taskdir}/output/${family}.tree" \
        "${taskdir}/output/${family}_annotated.tree" \
        "${TASKFILE}/${taskname}/input/syngraph.txt" \
        "${TASKFILE}/${taskname}/input/flat.ortholog_groups.tsv"

}
