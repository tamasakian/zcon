#!/usr/bin/env zsh
# Last updated: 2024-10-30

#### Function libs with GS2 (Matshui et al., 2019)
: << 'FUNCTIONS'
estimate_gstree: Reconstruct gene trees using the Graph Splitting method.
FUNCTIONS


function estimate_gstree() {
    function usage() {
        cat <<EOS
Usage: estimate_gstree <arg1> <arg2> <arg3> <arg4> <arg5>

    arg1: num_orgs  <- Number of organisms.                 (e.g., 1)
    arg2: org       <- Names  of organisms.                 (e.g., "Arabidopsis thaliana")
    arg3: num_peps  <- Number of proteins.                  (e.g., 2)
    arg4: peps      <- Accession IDs and Names of proteins. (e.g., "NP_001326314.1" "AtFLD" "NP_195315.3" "AtFD")
    arg5: threads   <- Number of threads                    (e.g., 2)

EOS
        exit 1
    }

    function parse_args() {
        if [[ $# -lt 8 ]]; then
            usage
        fi

        num_orgs=$1; shift
        orgs=()
        for ((i=1; i<=num_orgs; i++)); do
            orgs[i]="${1// /_}"; shift
            echo "${i}:${orgs[i]}"
        done
        num_peps=$1; shift
        typeset -g -A peps
        for ((i=1; i<=num_peps; i++)); do
            key="$1"; shift
            value="$1"; shift
            peps[$key]="$value"
            echo "${i}: ${key}"
            echo "${i}: ${peps[$key]}"
        done
        threads=$1
    }

    function merge_fasta() {
        touch ${taskdir}/database.fasta
        for org in "${orgs[@]}"; do
            genus=${org%%_*}
            cat ${DATA}/${genus}/${org}.pep.all.fasta >> ${taskdir}/database.fasta
        done
    }

    function generate_fasta() {
        touch ${taskdir}/pep.fasta
        makeblastdb \
            -in ${taskdir}/database.fasta \
            -dbtype prot \
            -hash_index \
            -parse_seqids

        for pep_key in ${(@k)peps}; do
            tmpfile=$(mktemp)
            blastdbcmd \
                -entry $pep_key \
                -db ${taskdir}/database.fasta \
                -out $tmpfile
            python3 -m fasp rename_header \
                $tmpfile \
                $tmpfile \
                ${peps[$pep_key]} \
                ""
            cat $tmpfile >> ${taskdir}/pep.fasta
            rm $tmpfile
        done
    }

    function main() {
        parse_args "$@"
        make_taskdir
        merge_fasta
        generate_fasta
        gs2 \
            -e 100 \
            -t $threads \
            -l \
            ${taskdir}/pep.fasta > ${taskdir}/pep.nwk
    }

    main "$@"
}