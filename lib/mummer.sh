#!/usr/bin/env zsh

#### Function libs with MUMmer4

function promer_dna_flanking_region_genus_symbol() {
    function usage() {
        cat <<EOS
Usage:  promer_dna_flanking_region_genus_symbol <arg1> <arg2> <arg3> <arg4> <arg5> (<arg6> <arg7> ...)

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

    function promer_dna_flanking_region() {
        pyvenv -m biotp remove_n_from_seqs \
            "${taskdir}/${genus}.${symbol}.dna_flanking_region.fasta" \
            "${taskdir}/${genus}.${symbol}.dna_flanking_region.fasta"
        pyvenv -m biotp split_multi_into_single \
            "${taskdir}/${genus}.${symbol}.dna_flanking_region.fasta" \
            "${taskdir}"
        for pnm1 in ${pnm_li[*]}; do
            for pnm2 in ${pnm_li[*]}; do
                promer \
                    --prefix "${taskdir}/${pnm1}.${pnm2}" \
                    "${taskdir}/${pnm1}.fasta" \
                    "${taskdir}/${pnm2}.fasta"
                mummerplot \
                    "${taskdir}/${pnm1}.${pnm2}.delta" \
                    --png \
                    --prefix "${taskdir}/${pnm1}.${pnm2}" \
                    --large
            done
        done
    }

    function main() {
        construct_dna_flanking_region_genus_symbol "$@"
        promer_dna_flanking_region
    }

    main "$@"
}