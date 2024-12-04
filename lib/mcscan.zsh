#!/usr/bin/env zsh
# Last updated: 2024-12-03
# Tools: JCVI LAST
# Function libs with MCscan in JCVI.
: << 'FUNCTIONS'
calculate_syntenic_depth: Calculate syntenic depth.
search_one_to_one_microsynteny: Search 1:1 microsynteny.
FUNCTIONS

function calculate_syntenic_depth() {
    function usage() {
        cat <<EOS
Usage:  calculate_syntenic_depth <arg1> <arg2> <arg3> <arg4>

    arg1: ref
    arg2: qry
    arg3: ref_feat
    arg4: qry_feat

EOS
        exit 1
    }

    function parse_args() {
        if [[ $# != 4 ]]; then
            usage
        fi
        ref="${1// /_}"
        qry="${2// /_}"
        ref_feat=$3
        qry_feat=$4
        ref_genus=${ref%%_*}
        qry_genus=${qry%%_*}
    }

    function to_bed() {
        cd $taskdir
        python3 -m jcvi.formats.gff bed --type=gene --key=${ref_feat} "${DATA}/${ref_genus}/${ref}.genome.gff" -o "${ref}.bed"
        python3 -m jcvi.formats.gff bed --type=gene --key=${qry_feat} "${DATA}/${qry_genus}/${qry}.genome.gff" -o "${qry}.bed"
        cd $ROOT
    }
    
    function to_cds() {
        cd $taskdir
        python3 -m biotp rename_headers_feature \
            "${DATA}/${ref_genus}/${ref}.cds.all.fasta" \
            "${taskdir}/${ref}.cds.${ref_feat}.fasta" \
            "$ref_feat"
        python3 -m biotp rename_headers_feature \
            "${DATA}/${qry_genus}/${qry}.cds.all.fasta" \
            "${taskdir}/${qry}.cds.${qry_feat}.fasta" \
            "$qry_feat"
        python3 -m jcvi.formats.fasta format "${taskdir}/${ref}.cds.${ref_feat}.fasta" "${ref}.cds"
        python3 -m jcvi.formats.fasta format "${taskdir}/${qry}.cds.${qry_feat}.fasta" "${qry}.cds"
        cd $ROOT
    }

    function calculate_depth() {
        cd $taskdir
        python3 -m jcvi.compara.catalog ortholog "$ref" "$qry" --no_strip_names
        python3 -m jcvi.compara.synteny depth --histogram "${ref}.${qry}.anchors"
        cd $ROOT
    }

    function main() {
        parse_args "$@"
        make_taskdir
        to_bed
        to_cds
        calculate_depth
    }

    main "$@"
}

function search_one_to_one_microsynteny() {
    function usage() {
        cat <<EOS
Usage:  search_one_to_one_macrosynteny <arg1> <arg2> <arg3> <arg4>

    arg1: ref
    arg2: qry
    arg3: ref_feat
    arg4: qry_feat

EOS
        exit 1
    }

    function parse_args() {
        if [[ $# != 4 ]]; then
            usage
        fi
        ref="${1// /_}"
        qry="${2// /_}"
        ref_feat=$3
        qry_feat=$4
        ref_genus=${ref%%_*}
        qry_genus=${qry%%_*}
    }

    function to_bed() {
        cd $taskdir
        python3 -m jcvi.formats.gff bed --type=gene --key=${ref_feat} "${DATA}/${ref_genus}/${ref}.genome.gff" -o "${ref}.bed"
        python3 -m jcvi.formats.gff bed --type=gene --key=${qry_feat} "${DATA}/${qry_genus}/${qry}.genome.gff" -o "${qry}.bed"
        cd $ROOT
    }
    
    function to_cds() {
        cd $taskdir
        python3 -m biotp rename_headers_feature \
            "${DATA}/${ref_genus}/${ref}.cds.all.fasta" \
            "${taskdir}/${ref}.cds.${ref_feat}.fasta" \
            "$ref_feat"
        python3 -m biotp rename_headers_feature \
            "${DATA}/${qry_genus}/${qry}.cds.all.fasta" \
            "${taskdir}/${qry}.cds.${qry_feat}.fasta" \
            "$qry_feat"
        python3 -m jcvi.formats.fasta format "${taskdir}/${ref}.cds.${ref_feat}.fasta" "${ref}.cds"
        python3 -m jcvi.formats.fasta format "${taskdir}/${qry}.cds.${qry_feat}.fasta" "${qry}.cds"
        cd $ROOT
    }

    function search_microsynteny() {
        cd $taskdir
        python3 -m jcvi.compara.catalog ortholog "$ref" "$qry" --no_strip_names
        python3 -m jcvi.compara.synteny mcscan "${ref}.bed" "${ref}.${qry}.lifted.anchors" --iter=1 -o "${ref}.${qry}.i1.blocks"
        cd $ROOT
    }

    function main() {
        parse_args "$@"
        make_taskdir
        to_bed
        to_cds
        search_microsynteny
    }

    main "$@"
}