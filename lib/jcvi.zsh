#!/usr/bin/env zsh

#### Function libs with JCVI

function visualize_synteny_depth() {
    function usage() {
        cat <<EOS
Usage:  visualize_synteny_depth <arg1> <arg2> <arg3> <arg4>

    arg1: genus
    arg2: ref
    arg3: qry
    arg4: feat

EOS
        exit 1
    }

    function parse_args() {
        if [[ $# != 4 ]]; then
            usage
        fi
        genus="$1"
        ref="$2"
        qry="$3"    
        feat=$4    
    }

    function to_bed() {
        cd $taskdir
        for org in "$ref" "$qry"; do
            org_us=${org/ /_}
            echo "${org_us}.genome.gff -> ${org_us}.bed"
            python3 -m jcvi.formats.gff bed --type=gene --key=${feat} "${DATA}/${genus}/${org_us}.genome.gff" -o "${org_us}.bed"
        done
        cd $ROOT
    }

    function to_cds() {
        cd $taskdir
        for org in "$ref" "$qry"; do
            org_us=${org/ /_}
            echo "${org_us}.cds.all.fasta -> ${org_us}.cds"
            python3 -m biotp rename_headers_to_features \
                "${DATA}/${genus}/${org_us}.cds.all.fasta" \
                "${taskdir}/${org_us}.cds.all.${feat}.fasta" \
                "$feat" 
            python3 -m jcvi.formats.fasta format ${taskdir}/${org_us}.cds.all.${feat}.fasta ${org_us}.cds
        done
        cd $ROOT

    }

    function search_pairwise_synteny() {
        cd $taskdir
        ref_us=${ref/ /_}
        qry_us=${qry/ /_}
        python3 -m jcvi.compara.catalog ortholog ${ref_us} ${qry_us} --no_strip_names
        python3 -m jcvi.graphics.dotplot ${ref_us}.${qry_us}.anchors
        python3 -m jcvi.compara.synteny depth --histogram ${ref_us}.${qry_us}.anchors
        cd $ROOT
    }

    function main() {
        parse_args "$@"
        make_taskdir
        to_bed
        to_cds
        search_pairwise_synteny
    }

    main "$@"
}

function search_one_to_one_synteny() {
    function usage() {
        cat <<EOS
Usage:  search_one_to_one_synteny <arg1> <arg2> <arg3> <arg4>

    arg1: genus
    arg2: ref
    arg3: qry
    arg4: feature

EOS
        exit 1
    }

    function parse_args() {
        if [[ $# != 4 ]]; then
            usage
        fi
        genus="$1"
        ref="$2"
        qry="$3"
        feat=$4
        ref_us=${ref/ /_}
        qry_us=${qry/ /_}
    }

    function to_bed() {
        cd $taskdir
        for org in "$ref" "$qry"; do
            org_us=${org/ /_}
            echo "${org_us}.genome.gff -> ${org_us}.bed"
            python3 -m jcvi.formats.gff bed --type=gene --key=${feat} "${DATA}/${genus}/${org_us}.genome.gff" -o "${org_us}.bed"
        done
        cd $ROOT
    }
    
    function to_cds() {
        cd $taskdir
        for org in "$ref" "$qry"; do
            org_us=${org/ /_}
            echo "${org_us}.cds.all.fasta -> ${org_us}.cds"
            python3 -m biotp rename_headers_to_features \
                "${DATA}/${genus}/${org_us}.cds.all.fasta" \
                "${taskdir}/${org_us}.cds.${feat}.fasta" \
                "$feat" 
            python3 -m jcvi.formats.fasta format "${taskdir}/${org_us}.cds.${feat}.fasta" "${org_us}.cds"
        done
        cd $ROOT
    }

    function search_microsynteny() {
        cd $taskdir
        python3 -m jcvi.compara.catalog ortholog "${ref_us}" "${qry_us}" --no_strip_names
        python3 -m jcvi.compara.synteny mcscan "${ref_us}.bed" "${ref_us}.${qry_us}.lifted.anchors" --iter=1 -o "${ref_us}.${qry_us}.i1.blocks"
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

function search_one_to_two_synteny() {
    function usage() {
        cat <<EOS
Usage:  search_one_to_two_synteny <arg1> <arg2> <arg3> <arg4>


    arg1: genus
    arg2: ref
    arg3: qry
    arg4: feature

EOS
        exit 1
    }

    function parse_args() {
        if [[ $# != 4 ]]; then
            usage
        fi
        genus="$1"
        ref="$2"
        qry="$3"
        feat=$4
        ref_us=${ref/ /_}
        qry_us=${qry/ /_}
    }

    function to_bed() {
        cd $taskdir
        for org in "$ref" "$qry"; do
            org_us=${org/ /_}
            echo "${org_us}.genome.gff -> ${org_us}.bed"
            python3 -m jcvi.formats.gff bed --type=gene --key=${feat} "${DATA}/${genus}/${org_us}.genome.gff" -o "${org_us}.bed"
        done
        cd $ROOT
    }
    
    function to_cds() {
        cd $taskdir
        for org in "$ref" "$qry"; do
            org_us=${org/ /_}
            echo "${org_us}.cds.all.fasta -> ${org_us}.cds"
            python3 -m biotp rename_headers_to_features \
                "${DATA}/${genus}/${org_us}.cds.all.fasta" \
                "${taskdir}/${org_us}.cds.${feat}.fasta" \
                "$feat" 
            python3 -m jcvi.formats.fasta format "${taskdir}/${org_us}.cds.${feat}.fasta" "${org_us}.cds"
        done
        cd $ROOT
    }

    function search_microsynteny() {
        cd $taskdir
        python3 -m jcvi.compara.catalog ortholog "${ref_us}" "${qry_us}" --no_strip_names
        python3 -m jcvi.compara.synteny mcscan "${ref_us}.bed" "${ref_us}.${qry_us}.lifted.anchors" --iter=2 -o "${ref_us}.${qry_us}.i2.blocks"
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

function search_one_to_three_synteny() {
    function usage() {
        cat <<EOS
Usage:  search_one_to_three_synteny <arg1> <arg2> <arg3> <arg4>

    arg1: genus
    arg2: ref
    arg3: qry
    arg4: feature

EOS
        exit 1
    }

    function parse_args() {
        if [[ $# != 4 ]]; then
            usage
        fi
        genus="$1"
        ref="$2"
        qry="$3"
        feat=$4
        ref_us=${ref/ /_}
        qry_us=${qry/ /_}
    }

    function to_bed() {
        cd $taskdir
        for org in "$ref" "$qry"; do
            org_us=${org/ /_}
            echo "${org_us}.genome.gff -> ${org_us}.bed"
            python3 -m jcvi.formats.gff bed --type=gene --key=${feat} "${DATA}/${genus}/${org_us}.genome.gff" -o "${org_us}.bed"
        done
        cd $ROOT
    }
    
    function to_cds() {
        cd $taskdir
        for org in "$ref" "$qry"; do
            org_us=${org/ /_}
            echo "${org_us}.cds.all.fasta -> ${org_us}.cds"
            python3 -m biotp rename_headers_to_features \
                "${DATA}/${genus}/${org_us}.cds.all.fasta" \
                "${taskdir}/${org_us}.cds.${feat}.fasta" \
                "$feat" 
            python3 -m jcvi.formats.fasta format "${taskdir}/${org_us}.cds.${feat}.fasta" "${org_us}.cds"
        done
        cd $ROOT
    }

    function search_microsynteny() {
        cd $taskdir
        python3 -m jcvi.compara.catalog ortholog "${ref_us}" "${qry_us}" --no_strip_names
        python3 -m jcvi.compara.synteny mcscan "${ref_us}.bed" "${ref_us}.${qry_us}.lifted.anchors" --iter=3 -o "${ref_us}.${qry_us}.i3.blocks"
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

function search_longest_one_to_one_microsynteny() {
    function usage() {
        cat <<EOS
Usage:  search_longest_one_to_one_microsynteny <arg1> <arg2> <arg3> <arg4>

    arg1: genus
    arg2: ref
    arg3: qry
    arg4: feat

EOS
        exit 1
    }

    function output_longest_microsynteny() {
        python3 -m biotp output_longest_one_to_one_microsynteny \
            "${taskdir}/${ref_us}.bed" \
            "${taskdir}/${qry_us}.bed" \
            "${taskdir}/${ref_us}.${qry_us}.i1.blocks" \
            "${taskdir}/longest_microsynteny.csv"
    }

    function main() {
        search_one_to_one_synteny "$@"
        output_longest_microsynteny
    }

    main "$@"
}

function search_longest_one_to_two_microsynteny() {
    function usage() {
        cat <<EOS
Usage:  search_longest_one_to_two_microsynteny <arg1> <arg2> <arg3> <arg4>

    arg1: genus
    arg2: ref
    arg3: qry
    arg4: feat

EOS
        exit 1
    }

    function output_longest_microsynteny() {
        python3 -m biotp output_longest_one_to_two_microsynteny \
            "${taskdir}/${ref_us}.bed" \
            "${taskdir}/${qry_us}.bed" \
            "${taskdir}/${ref_us}.${qry_us}.i2.blocks" \
            "${taskdir}/longest_microsynteny.csv"
    }

    function main() {
        search_one_to_two_synteny "$@"
        output_longest_microsynteny
    }

    main "$@"
}

function search_longest_one_to_three_microsynteny() {
    function usage() {
        cat <<EOS
Usage:  search_longest_one_to_three_microsynteny <arg1> <arg2> <arg3> <arg4>

    arg1: genus
    arg2: ref
    arg3: qry
    arg4: feat

EOS
        exit 1
    }

    function output_longest_microsynteny() {
        python3 -m biotp output_longest_one_to_three_microsynteny \
            "${taskdir}/${ref_us}.bed" \
            "${taskdir}/${qry_us}.bed" \
            "${taskdir}/${ref_us}.${qry_us}.i3.blocks" \
            "${taskdir}/longest_microsynteny.csv"
    }

    function main() {
        search_one_to_three_synteny "$@"
        output_longest_microsynteny
    }

    main "$@"
}

function visualize_one_to_one_microsynteny() {
    function usage() {
        cat <<EOS
Usage:  visualize_one_to_one_microsynteny <arg1> <arg2> <arg3> <arg4> (<arg5> ...)

    arg1: genus
    arg2: ref
    arg3: qry
    arg4: feature
    arg5: seqid
    ...

EOS
        exit 1
    }

    function parse_args() {
        if [[ $# -lt 6 ]]; then
            usage
        fi
        genus="$1"
        ref="$2"
        qry="$3"
        feat=$4
        seq_li=("${@:5}")
    }

    function to_bed() {
        cd $taskdir
        for org in "$ref" "$qry"; do
            org_us=${org/ /_}
            echo "${org_us}.genome.gff -> ${org_us}.bed"
            python3 -m biotp slice_lines_by_seqids \
                "${DATA}/${genus}/${org_us}.genome.gff" \
                "${taskdir}/${org_us}.genome.sliced.gff" \
                "${seq_li[@]}"
            python3 -m jcvi.formats.gff bed --type=gene --key=${feat} "${taskdir}/${org_us}.genome.sliced.gff" -o "${org_us}.bed"
        done
        cd $ROOT
    }
    
    function to_cds() {
        cd $taskdir
        for org in "$ref" "$qry"; do
            org_us=${org/ /_}
            echo "${org_us}.cds.all.fasta -> ${org_us}.cds"
            python3 -m biotp slice_headers_by_ids \
                "${DATA}/${genus}/${org_us}.cds.all.fasta" \
                "${taskdir}/${org_us}.cds.sliced.fasta" \
                "${seq_li[@]}"
            python3 -m biotp rename_headers_to_features \
                "${taskdir}/${org_us}.cds.sliced.fasta" \
                "${taskdir}/${org_us}.cds.sliced.${feat}.fasta" \
                "$feat" 
            python3 -m jcvi.formats.fasta format "${taskdir}/${org_us}.cds.sliced.${feat}.fasta" "${org_us}.cds"
        done
        cd $ROOT
    }

    function search_microsynteny() {
        cd $taskdir
        ref_us=${ref/ /_}
        qry_us=${qry/ /_}
        python3 -m jcvi.compara.catalog ortholog "${ref_us}" "${qry_us}" --no_strip_names
        python3 -m jcvi.compara.synteny mcscan "${ref_us}.bed" "${ref_us}.${qry_us}.lifted.anchors" --iter=1 -o "${ref_us}.${qry_us}.i1.blocks"
        touch "blocks"
        cat "${ref_us}.${qry_us}.i1.blocks" > "blocks"
        touch "blocks.layout"
        echo "# x, y, rotation, ha, va, color, ratio, label" >> "blocks.layout"
        echo "0.5, 0.4, 0, left, center, #009E73, 1, ${seq_li[1]}" >> "blocks.layout"
        echo "0.5, 0.6, 0, left, center, #E69F00, 1, ${seq_li[2]}" >> "blocks.layout"
        echo "# edges" >> "blocks.layout"
        echo "e, 0, 1" >> "blocks.layout"
        python3 -m jcvi.formats.bed merge "${ref_us}.bed" "${qry_us}.bed" -o "${ref_us}.${qry_us}.bed"
        python3 -m jcvi.graphics.synteny blocks "${ref_us}.${qry_us}.bed" "blocks.layout" --glyphstyle=arrow --shadestyle=line
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

function visualize_one_to_two_microsynteny() {
    function usage() {
        cat <<EOS
Usage:  visualize_one_to_two_microsynteny <arg1> <arg2> <arg3> <arg4> (<arg5> ...)


    arg1: genus
    arg2: ref
    arg3: qry
    arg4: feature
    arg5: seqid
    ...

EOS
        exit 1
    }

    function parse_args() {
        if [[ $# -lt 6 ]]; then
            usage
        fi
        genus="$1"
        ref="$2"
        qry="$3"
        feat=$4
        seq_li=("${@:5}")
    }

    function to_bed() {
        cd $taskdir
        for org in "$ref" "$qry"; do
            org_us=${org/ /_}
            echo "${org_us}.genome.gff -> ${org_us}.bed"
            python3 -m biotp slice_lines_by_seqids \
                "${DATA}/${genus}/${org_us}.genome.gff" \
                "${taskdir}/${org_us}.genome.sliced.gff" \
                "${seq_li[@]}"
            python3 -m jcvi.formats.gff bed --type=gene --key=${feat} "${taskdir}/${org_us}.genome.sliced.gff" -o "${org_us}.bed"
        done
        cd $ROOT
    }
    
    function to_cds() {
        cd $taskdir
        for org in "$ref" "$qry"; do
            org_us=${org/ /_}
            echo "${org_us}.cds.all.fasta -> ${org_us}.cds"
            python3 -m biotp slice_headers_by_ids \
                "${DATA}/${genus}/${org_us}.cds.all.fasta" \
                "${taskdir}/${org_us}.cds.sliced.fasta" \
                "${seq_li[@]}"
            python3 -m biotp rename_headers_to_features \
                "${taskdir}/${org_us}.cds.sliced.fasta" \
                "${taskdir}/${org_us}.cds.sliced.${feat}.fasta" \
                "$feat" 
            python3 -m jcvi.formats.fasta format "${taskdir}/${org_us}.cds.sliced.${feat}.fasta" "${org_us}.cds"
        done
        cd $ROOT
    }

    function search_microsynteny() {
        cd $taskdir
        ref_us=${ref/ /_}
        qry_us=${qry/ /_}
        python3 -m jcvi.compara.catalog ortholog "${ref_us}" "${qry_us}" --no_strip_names
        python3 -m jcvi.compara.synteny mcscan "${ref_us}.bed" "${ref_us}.${qry_us}.lifted.anchors" --iter=2 -o "${ref_us}.${qry_us}.i2.blocks"
        touch "blocks"
        cat "${ref_us}.${qry_us}.i2.blocks" > "blocks"
        touch "blocks.layout"
        echo "# x, y, rotation, ha, va, color, ratio, label" >> "blocks.layout"
        echo "0.5, 0.5, 0, left, center, #009E73, 1, ${seq_li[1]}" >> "blocks.layout"
        echo "0.5, 0.3, 0, left, center, #E69F00, 1, ${seq_li[2]}" >> "blocks.layout"
        echo "0.5, 0.7, 0, left, center, #56B4E9, 1, ${seq_li[3]}" >> "blocks.layout"
        echo "# edges" >> "blocks.layout"
        echo "e, 0, 1" >> "blocks.layout"
        echo "e, 0, 2" >> "blocks.layout"
        python3 -m jcvi.formats.bed merge "${ref_us}.bed" "${qry_us}.bed" -o "${ref_us}.${qry_us}.bed"
        python3 -m jcvi.graphics.synteny blocks "${ref_us}.${qry_us}.bed" "blocks.layout" --glyphstyle=arrow --shadestyle=line
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

function visualize_one_to_three_microsynteny() {
    function usage() {
        cat <<EOS
Usage:  visualize_one_to_three_microsynteny <arg1> <arg2> <arg3> <arg4> (<arg5> ...)


    arg1: genus
    arg2: ref
    arg3: qry
    arg4: feature
    arg5: seqid
    ...

EOS
        exit 1
    }

    function parse_args() {
        if [[ $# -lt 6 ]]; then
            usage
        fi
        genus="$1"
        ref="$2"
        qry="$3"
        feat=$4
        seq_li=("${@:5}")
    }

    function to_bed() {
        cd $taskdir
        for org in "$ref" "$qry"; do
            org_us=${org/ /_}
            echo "${org_us}.genome.gff -> ${org_us}.bed"
            python3 -m biotp slice_lines_by_seqids \
                "${DATA}/${genus}/${org_us}.genome.gff" \
                "${taskdir}/${org_us}.genome.sliced.gff" \
                "${seq_li[@]}"
            python3 -m jcvi.formats.gff bed --type=gene --key=${feat} "${taskdir}/${org_us}.genome.sliced.gff" -o "${org_us}.bed"
        done
        cd $ROOT
    }
    
    function to_cds() {
        cd $taskdir
        for org in "$ref" "$qry"; do
            org_us=${org/ /_}
            echo "${org_us}.cds.all.fasta -> ${org_us}.cds"
            python3 -m biotp slice_headers_by_ids \
                "${DATA}/${genus}/${org_us}.cds.all.fasta" \
                "${taskdir}/${org_us}.cds.sliced.fasta" \
                "${seq_li[@]}"
            python3 -m biotp rename_headers_to_features \
                "${taskdir}/${org_us}.cds.sliced.fasta" \
                "${taskdir}/${org_us}.cds.sliced.${feat}.fasta" \
                "$feat" 
            python3 -m jcvi.formats.fasta format "${taskdir}/${org_us}.cds.sliced.${feat}.fasta" "${org_us}.cds"
        done
        cd $ROOT
    }

    function search_microsynteny() {
        cd $taskdir
        ref_us=${ref/ /_}
        qry_us=${qry/ /_}
        python3 -m jcvi.compara.catalog ortholog "${ref_us}" "${qry_us}" --no_strip_names
        python3 -m jcvi.compara.synteny mcscan "${ref_us}.bed" "${ref_us}.${qry_us}.lifted.anchors" --iter=3 -o "${ref_us}.${qry_us}.i3.blocks"
        touch "blocks"
        cat "${ref_us}.${qry_us}.i3.blocks" > "blocks"
        touch "blocks.layout"
        echo "# x, y, rotation, ha, va, color, ratio, label" >> "blocks.layout"
        echo "0.5, 0.6, 0, center, top, #009E73, 1, ${seq_li[1]}" >> "blocks.layout"
        echo "0.1, 0.4, 0, center, bottom, #E69F00, .3, ${seq_li[2]}" >> "blocks.layout"
        echo "0.4, 0.4, 0, center, bottom, #56B4E9, .3, ${seq_li[3]}" >> "blocks.layout"
        echo "0.7, 0.4, 0, center, bottom, #CC79A7, .3, ${seq_li[4]}" >> "blocks.layout"
        echo "# edges" >> "blocks.layout"
        echo "e, 0, 1" >> "blocks.layout"
        echo "e, 0, 2" >> "blocks.layout"
        echo "e, 0, 3" >> "blocks.layout"
        python3 -m jcvi.formats.bed merge "${ref_us}.bed" "${qry_us}.bed" -o "${ref_us}.${qry_us}.bed"
        python3 -m jcvi.graphics.synteny blocks "${ref_us}.${qry_us}.bed" "blocks.layout" --glyphstyle=arrow --shadestyle=line
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

function make_mltree_one_to_one_microsynteny() {
    function usage() {
        cat <<EOS
Usage:  make_mltree_one_to_one_microsynteny <arg1> <arg2> <arg3> <arg4> <arg5> <arg6> <arg7> (<arg8> ...)

    arg1: genus
    arg2: ref
    arg3: qry
    arg4: feat
    arg5: out_genus
    arg6: out
    arg7: out_feat
    arg8: seqid
    ...

EOS
        exit 1
    } 

    function parse_args() {
        if [[ $# -lt 8 ]]; then
            usage
        fi
        genus="$1"
        ref="$2"
        qry="$3"
        feat="$4"
        out_genus="$5"
        out="${6/ /_}"
        out_feat="$7"
        seq_li=("${@:8}")
    }

    function to_bed() {
        cd $taskdir
        for org in "$ref" "$qry"; do
            org_us=${org/ /_}
            echo "${org_us}.genome.gff -> ${org_us}.bed"
            python3 -m biotp slice_lines_by_seqids \
                "${DATA}/${genus}/${org_us}.genome.gff" \
                "${taskdir}/${org_us}.genome.sliced.gff" \
                "${seq_li[@]}"
            python3 -m jcvi.formats.gff bed --type=gene --key=${feat} "${taskdir}/${org_us}.genome.sliced.gff" -o "${org_us}.bed"
        done
        cd $ROOT
    }
    
    function to_cds() {
        cd $taskdir
        for org in "$ref" "$qry"; do
            org_us=${org/ /_}
            echo "${org_us}.cds.all.fasta -> ${org_us}.cds"
            python3 -m biotp slice_headers_by_ids \
                "${DATA}/${genus}/${org_us}.cds.all.fasta" \
                "${taskdir}/${org_us}.cds.sliced.fasta" \
                "${seq_li[@]}"
            python3 -m biotp rename_headers_to_features \
                "${taskdir}/${org_us}.cds.sliced.fasta" \
                "${taskdir}/${org_us}.cds.sliced.${feat}.fasta" \
                "$feat" 
            python3 -m jcvi.formats.fasta format "${taskdir}/${org_us}.cds.sliced.${feat}.fasta" "${org_us}.cds"
        done
        cd $ROOT
    }

    function search_microsynteny() {
        cd $taskdir
        ref_us=${ref/ /_}
        qry_us=${qry/ /_}
        python3 -m jcvi.compara.catalog ortholog ${ref_us} ${qry_us} --no_strip_names
        python3 -m jcvi.compara.synteny mcscan ${ref_us}.bed ${ref_us}.${qry_us}.lifted.anchors --iter=1 -o ${ref_us}.${qry_us}.i1.blocks
        cd $ROOT
    }

    function makeblastdb_microsynteny() {
        cd $taskdir
        python3 -m biotp rename_headers_to_features \
            "${DATA}/${out_genus}/${out}.cds.all.fasta" \
            "${taskdir}/${out}.cds" \
            "${out_feat}"
        ref_us=${ref/ /_}
        qry_us=${qry/ /_}
        for org_us in "$ref_us" "$qry_us" "$out"; do
            makeblastdb \
                -in ${taskdir}/${org_us}.cds \
                -dbtype nucl \
                -hash_index \
                -parse_seqids
        done
        cd $ROOT
    }

    function declare_microsynteny() {
        cd $taskdir
        ref_us=${ref/ /_}
        qry_us=${qry/ /_}
        tmpfile=$(mktemp)
        python3 -m biotp output_blocks \
            "${ref_us}.${qry_us}.i1.blocks" \
                > "$tmpfile"

        ref_li=()
        qry_li=()
        while IFS=" " read -r index ref qry; do
            echo "${index}, ref: ${ref}, qry: ${qry}"
            ref_li[index]="$ref"
            qry_li[index]="$qry"
        done < "$tmpfile"
        rm "$tmpfile"
        cd $ROOT
    }

    function retrieve_microsynteny() {
        outflag_li=()

        cd $taskdir
        echo "Nucleotide BLAST"
        for ((i=1; i<${#ref_li[@]}+1; i++)); do
            echo "${i}, ref: ${ref_li[i]}, qry: ${qry_li[i]}"
            outflag_li[i]=false

            ## add ref genes
            blastdbcmd \
                -entry "${ref_li[i]}" \
                -db ${taskdir}/${ref_us}.cds \
                -out ${ref_li[i]}.cds
            python3 -m biotp rename_header \
                "${ref_li[i]}.cds" \
                "${ref_li[i]}.cds" \
                "${seq_li[0]}" \
                "" \
                ""

            echo "out: all [${out}], qry: ${ref_li[i]} [${ref_us}]"
            local _blastn; _blastn=$(blastn -outfmt 6 -evalue 10 -db "${taskdir}/${out}.cds" -query "${taskdir}/${ref_li[i]}.cds")
            local _id; _ids=$(echo "$_blastn" | cut -f 2 | head -n 2)

            ## add qry genes
            tmpfile=$(mktemp)
            blastdbcmd \
                -entry "${qry_li[i]}" \
                -db "${taskdir}/${qry_us}.cds" \
                -out "$tmpfile"
            python3 -m biotp rename_header \
                "$tmpfile" \
                "$tmpfile" \
                "${seq_li[1]}" \
                "" \
                ""
            cat "$tmpfile" >> "${ref_li[i]}.cds"
            rm "$tmpfile"

            ## flag for outgroup genes
            if [[ ${#_ids[@]} -lt 2 ]]; then
                echo "NOT 2 hits."
                outflag_li[i]=true
                continue
            fi

            ## add outgroup genes
            for _id in $_ids; do
                echo "Hit, out: ${_id}"
                tmpfile=$(mktemp)
                blastdbcmd \
                    -entry "$_id" \
                    -db ${taskdir}/${out_}.cds \
                    -out "$tmpfile"
                python3 -m biotp rename_header \
                    "$tmpfile" \
                    "$tmpfile" \
                    "${out}_${_id}" \
                    "" \
                    ""
                cat "$tmpfile" >> "${ref_li[i]}.cds"
                rm "$tmpfile"
            done
        done
        cd $ROOT
    }

    function construct_msa_microsynteny() {
        for ((i=1; i<${#ref_li[@]}+1; i++)); do
            mafft "${taskdir}/${ref_li[i]}.cds" > "${taskdir}/${ref_li[i]}.cds.aln"
        done
    }

    function estimate_mltree() {
        for ((i=1; i<${#ref_li[@]}+1; i++)); do
            mkdir -p "${taskdir}/trees"
            echo "${i}, outfrag: ${outflag_li[i]}"

            ## no outgroup genes
            if ${outflag_li[i]}; then
                continue
            fi

            raxml-ng \
                --msa "${taskdir}/${ref_li[i]}.cds.aln" \
                --all \
                --model GTR \
                --bs-trees 500 \
                --threads 8 \
                --redo
            Rscript ${SCRIPT}/make_tree.R \
                "${taskdir}/${ref_li[i]}.cds.aln.raxml.support" \
                "${taskdir}/trees/${seq_li[0]}.${ref_li[i]}.cds.aln.raxml.png"
        done
    }

    function main() {
        parse_args "$@"
        make_taskdir
        to_bed
        to_cds
        search_microsynteny
        makeblastdb_microsynteny
        declare_microsynteny
        retrieve_microsynteny
        construct_msa_microsynteny
        estimate_mltree
    }
    main "$@"
}

function make_mltree_one_to_two_microsynteny() {
    function usage() {
        cat <<EOS
Usage:  make_mltree_one_to_two_microsynteny <arg1> <arg2> <arg3> <arg4> <arg5> <arg6> <arg7> (<arg8> ...)

    arg1: genus
    arg2: ref
    arg3: qry
    arg4: feat
    arg5: out_genus
    arg6: out
    arg7: out_feat
    arg8: seqid
    ...

EOS
        exit 1
    } 

    function parse_args() {
        if [[ $# -lt 8 ]]; then
            usage
        fi
        genus="$1"
        ref="$2"
        qry="$3"
        feat="$4"
        out_genus="$5"
        out="${6/ /_}"
        out_feat="$7"
        seq_li=("${@:8}")
        ref_us=${ref/ /_}
        qry_us=${qry/ /_}
    }

    function to_bed() {
        cd $taskdir
        for org in "$ref" "$qry"; do
            org_us=${org/ /_}
            echo "${org_us}.genome.gff -> ${org_us}.bed"
            python3 -m biotp slice_lines_by_seqids \
                "${DATA}/${genus}/${org_us}.genome.gff" \
                "${taskdir}/${org_us}.genome.sliced.gff" \
                "${seq_li[@]}"
            python3 -m jcvi.formats.gff bed --type=gene --key=${feat} "${taskdir}/${org_us}.genome.sliced.gff" -o "${org_us}.bed"
        done
        cd $ROOT
    }
    
    function to_cds() {
        cd $taskdir
        for org in "$ref" "$qry"; do
            org_us=${org/ /_}
            echo "${org_us}.cds.all.fasta -> ${org_us}.cds"
            python3 -m biotp slice_headers_by_ids \
                "${DATA}/${genus}/${org_us}.cds.all.fasta" \
                "${taskdir}/${org_us}.cds.sliced.fasta" \
                "${seq_li[@]}"
            python3 -m biotp rename_headers_to_features \
                "${taskdir}/${org_us}.cds.sliced.fasta" \
                "${taskdir}/${org_us}.cds.sliced.${feat}.fasta" \
                "$feat"
            python3 -m jcvi.formats.fasta format "${taskdir}/${org_us}.cds.sliced.${feat}.fasta" "${org_us}.cds"
        done
        cd $ROOT
    }

    function search_microsynteny() {
        cd $taskdir
        python3 -m jcvi.compara.catalog ortholog ${ref_us} ${qry_us} --no_strip_names
        python3 -m jcvi.compara.synteny mcscan ${ref_us}.bed ${ref_us}.${qry_us}.lifted.anchors --iter=2 -o ${ref_us}.${qry_us}.i2.blocks
        cd $ROOT
    }

    function makeblastdb_microsynteny() {
        touch "${taskdir}/${ref_us}.${qry_us}.cds"
        cat "${taskdir}/${ref_us}.cds" >> "${taskdir}/${ref_us}.${qry_us}.cds"
        cat "${taskdir}/${qry_us}.cds" >> "${taskdir}/${ref_us}.${qry_us}.cds"
        python3 -m biotp rename_headers_to_features \
            "${DATA}/${out_genus}/${out}.cds.all.fasta" \
            "${taskdir}/${out}.cds" \
            "${out_feat}"
        for org in "${ref_us}.${qry_us}" "${out}"; do
            makeblastdb \
                -in "${taskdir}/${org}.cds" \
                -dbtype nucl \
                -hash_index \
                -parse_seqids
        done
    }

    function declare_microsynteny() {
        tmpfile=$(mktemp)
        python3 -m biotp output_blocks \
            "${taskdir}/${ref_us}.${qry_us}.i2.blocks" \
                > "$tmpfile"
        ref_li=()
        qry1_li=()
        qry2_li=()
        while IFS=" " read -r i ref qry1 qry2; do
            echo "${i}, ref: ${ref}, qry1: ${qry1}, qry2: ${qry2}"
            ref_li[i]="$ref"
            qry1_li[i]="$qry1"
            qry2_li[i]="$qry2"
        done < "$tmpfile"
        rm "$tmpfile"
    }

    function retrieve_microsynteny() {
        outflag_li=()

        echo "Nucleotide BLAST"
        for ((i=1; i<${#ref_li[@]}+1; i++)); do
            echo "${i}, ref: ${ref_li[i]}, qry1: ${qry1_li[i]}, qry2: ${qry2_li[i]}"
            outflag_li[i]=false

            ## add ref
            blastdbcmd \
                -entry "${ref_li[i]}" \
                -db "${taskdir}/${ref_us}.${qry_us}.cds" \
                -out "${taskdir}/${ref_li[i]}.cds"
            python3 -m biotp rename_header \
                "${taskdir}/${ref_li[i]}.cds" \
                "${taskdir}/${ref_li[i]}.cds" \
                "${seq_li[1]}" \
                "" \
                ""

            echo "out: all [${out}], qry: ${ref_li[i]} [${ref_us}]"
            local _blastn; _blastn=$(blastn -outfmt 6 -evalue 10 -db "${taskdir}/${out}.cds" -query "${taskdir}/${ref_li[i]}.cds")
            local _id; _id=$(echo "$_blastn" | cut -f 2 | head -n 1)
            
            ## add qry1
            tmpfile=$(mktemp)
            blastdbcmd \
                -entry "${qry1_li[i]}" \
                -db "${taskdir}/${ref_us}.${qry_us}.cds" \
                -out "$tmpfile"
            python3 -m biotp rename_header \
                "$tmpfile" \
                "$tmpfile" \
                "${seq_li[2]}" \
                "" \
                ""
            cat "$tmpfile" >> "${taskdir}/${ref_li[i]}.cds"
            rm "$tmpfile"

            ## add qry2
            tmpfile=$(mktemp)
            blastdbcmd \
                -entry "${qry2_li[i]}" \
                -db "${taskdir}/${ref_us}.${qry_us}.cds" \
                -out "$tmpfile"
            python3 -m biotp rename_header \
                "$tmpfile" \
                "$tmpfile" \
                "${seq_li[3]}" \
                "" \
                ""
            cat "$tmpfile" >> "${taskdir}/${ref_li[i]}.cds"
            rm "$tmpfile"

            ## nothing
            if [ -z "$_id" ]; then
                outflag_li[i]=true
                continue
            fi

            ## add out
            tmpfile=$(mktemp)
            blastdbcmd \
                -entry "$_id" \
                -db "${taskdir}/${out}.cds" \
                -out "$tmpfile"

            ## reverse Nucleotide BLAST
            blastdbcmd \
                -entry "$_id" \
                -db "${taskdir}/${out}.cds" \
                -out "$tmpfile"
            local _blastn; _blastn=$(blastn -outfmt 6 -evalue 10 -db "${taskdir}/${ref_us}.${qry_us}.cds" -query "$tmpfile")
            local _id; _id=$(echo "$_blastn" | cut -f 2 | head -n 1)

            ## Not exist
            if [ -z "$_id" ]; then
                outflag_li[i]=true
                continue
            fi

            ## No match
            if [[ "$_id" != "${ref_li[i]}" && "$_id" != "${qry1_li[i]}" && "$_id" != "${qry2_li[i]}" ]]; then
                outflag_li[i]=true
                continue
            fi

            ## 外群追加
            python3 -m biotp rename_header \
                "$tmpfile" \
                "$tmpfile" \
                "${out}" \
                "" \
                ""

            cat "$tmpfile" >> "${taskdir}/${ref_li[i]}.cds"
            rm "$tmpfile"
        done
    }

    function construct_msa_microsynteny() {
        for ((i=1; i<${#ref_li[@]}+1; i++)); do
            ## No out
            if ${outflag_li[i]}; then
                continue
            fi
            mafft "${taskdir}/${ref_li[i]}.cds" > "${taskdir}/${ref_li[i]}.cds.fasta"
        done
    }

    function merge_msa_microsynteny() {
        touch "${taskdir}/merged.fasta"
        for ((i=1; i<${#ref_li[@]}+1; i++)); do
            ## No out
            if ${outflag_li[i]}; then
                continue
            fi
            cat "${taskdir}/${ref_li[i]}.cds.fasta" >> "${taskdir}/merged.fasta"
        done
        python3 -m biotp merge_seqs_by_seqids \
            "${taskdir}/merged.fasta" \
            "${taskdir}/merged.fasta"
    } 

    function estimate_mltree() {
        raxml-ng \
            --msa "${taskdir}/merged.fasta" \
            --all \
            --model GTR \
            --bs-trees 500 \
            --threads 8 \
            --redo
        Rscript ${SCRIPT}/make_tree.R \
            "${taskdir}/merged.fasta.raxml.support" \
            "${taskdir}/merged.fasta.raxml.png"
    }

    function main() {
        parse_args "$@"
        make_taskdir
        to_bed
        to_cds
        search_microsynteny
        makeblastdb_microsynteny
        declare_microsynteny
        retrieve_microsynteny
        construct_msa_microsynteny
        merge_msa_microsynteny
        estimate_mltree
    }
    main "$@"
}

function compare_besthit_one_to_two_microsynteny() {
    function usage() {
        cat <<EOS
Usage:  compare_besthit_one_to_two_microsynteny <arg1> <arg2> <arg3> <arg4> <arg5> <arg6> <arg7>

    arg1: genus
    arg2: ref
    arg3: qry
    arg4: feature
    arg5: seqid_ref
    arg6: seqid_qry1
    arg7: seqid_qry2

EOS
        exit 1
    }

    function parse_args() {
        if [[ $# != 7 ]]; then
            usage
        fi
        genus="$1"
        ref="$2"
        qry="$3"
        feat="$4"
        seqid_ref="$5"
        seqid_qry1="$6"
        seqid_qry2="$7"
        ref_us=${ref/ /_}
        qry_us=${qry/ /_}
    }

    function to_bed() {
        cd $taskdir
        for org in "$ref" "$qry"; do
            org_us=${org/ /_}
            echo "${org_us}.genome.gff -> ${org_us}.bed"
            python3 -m biotp slice_lines_by_seqids \
                "${DATA}/${genus}/${org_us}.genome.gff" \
                "${taskdir}/${org_us}.genome.sliced.gff" \
                "${seqid_ref}" "${seqid_qry1}" "${seqid_qry2}"
            python3 -m jcvi.formats.gff bed --type=gene --key=${feat} "${taskdir}/${org_us}.genome.sliced.gff" -o "${org_us}.bed"
        done
        cd $ROOT
    }
    
    function to_cds() {
        cd $taskdir
        for org in "$ref" "$qry"; do
            org_us=${org/ /_}
            echo "${org_us}.cds.all.fasta -> ${org_us}.cds"
            python3 -m biotp slice_headers_by_ids \
                "${DATA}/${genus}/${org_us}.cds.all.fasta" \
                "${taskdir}/${org_us}.cds.sliced.fasta" \
                "${seqid_ref}" "${seqid_qry1}" "${seqid_qry2}"
            python3 -m biotp rename_headers_to_features \
                "${taskdir}/${org_us}.cds.sliced.fasta" \
                "${taskdir}/${org_us}.cds.sliced.${feat}.fasta" \
                "$feat" 
            python3 -m jcvi.formats.fasta format "${taskdir}/${org_us}.cds.sliced.${feat}.fasta" "${org_us}.cds"
        done
        cd $ROOT
    }

    function search_microsynteny() {
        cd $taskdir
        python3 -m jcvi.compara.catalog ortholog "${ref_us}" "${qry_us}" --no_strip_names
        python3 -m jcvi.compara.synteny mcscan "${ref_us}.bed" "${ref_us}.${qry_us}.lifted.anchors" --iter=2 -o "${ref_us}.${qry_us}.i2.blocks"
        cd $ROOT
    }

    function output_besthit() {
        python3 -m biotp output_besthit_one_to_two_microsynteny \
            "${taskdir}/${ref_us}.bed" \
            "${taskdir}/${qry_us}.bed" \
            "${taskdir}/${ref_us}.${qry_us}.i2.blocks" \
            "${taskdir}/${ref_us}.${qry_us}.anchors" \
            "${taskdir}/${seqid_ref}.${seqid_qry1}.${seqid_qry2}.csv"
    }

    function make_besthit_bar() {
        Rscript ${SCRIPT}/count_besthit.R \
            "${taskdir}/${seqid_ref}.${seqid_qry1}.${seqid_qry2}.csv" \
            "${taskdir}/${seqid_ref}.${seqid_qry1}.${seqid_qry2}.png" \
            "${seqid_ref}" \
            "${seqid_qry1}" \
            "${seqid_qry2}"
    }

    function main() {
        parse_args "$@"
        make_taskdir
        to_bed
        to_cds
        search_microsynteny
        output_besthit
        make_besthit_bar
    }

    main "$@"
}


function visualize_one_to_two_synteny() {
    function usage() {
        cat <<EOS
Usage:  visualize_one_to_two_synteny <arg1> <arg2> <arg3> <arg4> <arg5> <arg6>

    arg1: ref_org
    arg2: qry_org
    arg3: qry_feat
    arg4: qry_feat
    arg5: num_seqs
    arg6: seqs

EOS
        exit 1
    }

    function parse_args() {
        if [[ $# -lt 6 ]]; then
            usage
        fi
        ref="${1// /_}"
        ref_genus=${ref%%_*}
        shift
        qry="${1// /_}"
        qry_genus=${qry%%_*}
        shift
        qry_feat="$1"
        shift
        ref_feat="$1"
        shift
        num_seqs=$1
        shift
        for ((i=1; i<=num_seqs; i++)); do
            seqs+=("$1")
            shift
        done 
    }

    function to_bed() {
        cd $taskdir
        python3 -m biotp slice_lines_by_seqids \
            "${DATA}/${ref_genus}/${ref}.genome.gff" \
            "${taskdir}/${ref}.genome.slice.gff" \
            "${seqs[@]}"
        python3 -m jcvi.formats.gff bed \
            --type=gene \
            --key=${ref_feat} \
            "${taskdir}/${ref}.genome.slice.gff" \
            -o "${ref}.bed"

        python3 -m biotp slice_lines_by_seqids \
            "${DATA}/${qry_genus}/${qry}.genome.gff" \
            "${taskdir}/${qry}.genome.slice.gff" \
            "${seqs[@]}"
        python3 -m jcvi.formats.gff bed \
            --type=gene \
            --key=${qry_feat} \
            "${taskdir}/${qry}.genome.slice.gff" \
            -o "${qry}.bed"
        cd $ROOT
    }
    
    function to_cds() {
        cd $taskdir
        python3 -m biotp slice_headers_by_ids \
            "${DATA}/${ref_genus}/${ref}.cds.all.fasta" \
            "${taskdir}/${ref}.cds.slice.fasta" \
            "${seqs[@]}"
        python3 -m biotp rename_headers_to_features \
            "${taskdir}/${ref}.cds.slice.fasta" \
            "${taskdir}/${ref}.cds.slice.${ref_feat}.fasta" \
            "$ref_feat" 
        python3 -m jcvi.formats.fasta format \
            "${taskdir}/${ref}.cds.slice.${ref_feat}.fasta" \
            "${ref}.cds"

        python3 -m biotp slice_headers_by_ids \
            "${DATA}/${qry_genus}/${qry}.cds.all.fasta" \
            "${taskdir}/${qry}.cds.slice.fasta" \
            "${seqs[@]}"
        python3 -m biotp rename_headers_to_features \
            "${taskdir}/${qry}.cds.slice.fasta" \
            "${taskdir}/${qry}.cds.slice.${qry_feat}.fasta" \
            "$qry_feat" 
        python3 -m jcvi.formats.fasta format \
            "${taskdir}/${qry}.cds.slice.${qry_feat}.fasta" \
            "${qry}.cds"
        cd $ROOT
    }

    function search_microsynteny() {
        cd $taskdir
        python3 -m jcvi.compara.catalog ortholog "$ref" "$qry" --no_strip_names
        python3 -m jcvi.compara.synteny mcscan "${ref}.bed" "${ref}.${qry}.lifted.anchors" --iter=2 -o "${ref}.${qry}.i2.blocks"
        touch "blocks"
        cat "${ref}.${qry}.i2.blocks" > "blocks"
        touch "blocks.layout"
        echo "# x, y, rotation, ha, va, color, ratio, label" >> "blocks.layout"
        echo "0.5, 0.5, 0, left, center, #009E73, 1, ${seqs[1]}" >> "blocks.layout"
        echo "0.5, 0.3, 0, left, center, #E69F00, 1, ${seqs[2]}" >> "blocks.layout"
        echo "0.5, 0.7, 0, left, center, #56B4E9, 1, ${seqs[3]}" >> "blocks.layout"
        echo "# edges" >> "blocks.layout"
        echo "e, 0, 1" >> "blocks.layout"
        echo "e, 0, 2" >> "blocks.layout"
        python3 -m jcvi.formats.bed merge "${ref}.bed" "${qry}.bed" -o "${ref}.${qry}.bed"
        python3 -m jcvi.graphics.synteny blocks "${ref}.${qry}.bed" "blocks.layout" --glyphstyle=arrow --shadestyle=line
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


function compare_besthit_one_to_two_synteny() {
    function usage() {
        cat <<EOS
Usage:  compare_besthit_one_to_two_synteny <arg1> <arg2> <arg3> <arg4> <arg5>

    arg1: ref_org
    arg2: qry_org
    arg3: qry_feat
    arg4: qry_feat
    arg5: num_seqs
    arg6: seqs

EOS
        exit 1
    }

    function parse_args() {
        if [[ $# -lt 6 ]]; then
            usage
        fi
        ref="${1// /_}"
        ref_genus=${ref%%_*}
        shift
        qry="${1// /_}"
        qry_genus=${qry%%_*}
        shift
        ref_feat="$1"
        shift
        qry_feat="$1"
        shift
        num_seqs=$1
        shift
        for ((i=1; i<=num_seqs; i++)); do
            seqs+=("$1")
            shift
        done 
    }

    function to_bed() {
        cd $taskdir
        python3 -m biotp slice_lines_by_seqids \
            "${DATA}/${ref_genus}/${ref}.genome.gff" \
            "${taskdir}/${ref}.genome.slice.gff" \
            "${seqs[@]}"
        python3 -m jcvi.formats.gff bed \
            --type=gene \
            --key=${ref_feat} \
            "${taskdir}/${ref}.genome.slice.gff" \
            -o "${ref}.bed"

        python3 -m biotp slice_lines_by_seqids \
            "${DATA}/${qry_genus}/${qry}.genome.gff" \
            "${taskdir}/${qry}.genome.slice.gff" \
            "${seqs[@]}"
        python3 -m jcvi.formats.gff bed \
            --type=gene \
            --key=${qry_feat} \
            "${taskdir}/${qry}.genome.slice.gff" \
            -o "${qry}.bed"
        cd $ROOT
    }
    
    function to_cds() {
        cd $taskdir
        python3 -m biotp slice_headers_by_ids \
            "${DATA}/${ref_genus}/${ref}.cds.all.fasta" \
            "${taskdir}/${ref}.cds.slice.fasta" \
            "${seqs[@]}"
        python3 -m biotp rename_headers_to_features \
            "${taskdir}/${ref}.cds.slice.fasta" \
            "${taskdir}/${ref}.cds.slice.${ref_feat}.fasta" \
            "$ref_feat" 
        python3 -m jcvi.formats.fasta format \
            "${taskdir}/${ref}.cds.slice.${ref_feat}.fasta" \
            "${ref}.cds"

        python3 -m biotp slice_headers_by_ids \
            "${DATA}/${qry_genus}/${qry}.cds.all.fasta" \
            "${taskdir}/${qry}.cds.slice.fasta" \
            "${seqs[@]}"
        python3 -m biotp rename_headers_to_features \
            "${taskdir}/${qry}.cds.slice.fasta" \
            "${taskdir}/${qry}.cds.slice.${qry_feat}.fasta" \
            "$qry_feat" 
        python3 -m jcvi.formats.fasta format \
            "${taskdir}/${qry}.cds.slice.${qry_feat}.fasta" \
            "${qry}.cds"
        cd $ROOT
    }

    function search_microsynteny() {
        cd $taskdir
        python3 -m jcvi.compara.catalog ortholog "$ref" "$qry" --no_strip_names
        python3 -m jcvi.compara.synteny mcscan "${ref}.bed" "${ref}.${qry}.lifted.anchors" --iter=2 -o "${ref}.${qry}.i2.blocks"
        cd $ROOT
    }

    function output_besthit() {
        python3 -m biotp output_besthit_one_to_two_synteny \
            "${taskdir}/${ref}.bed" \
            "${taskdir}/${qry}.bed" \
            "${taskdir}/${ref}.${qry}.i2.blocks" \
            "${taskdir}/${ref}.${qry}.anchors" \
            "${taskdir}/hits.tsv"
    }

    function make_besthit_bar() {
        Rscript ${SCRIPT}/count_besthit.R \
            "${taskdir}/${seqid_ref}.${seqid_qry1}.${seqid_qry2}.csv" \
            "${taskdir}/${seqid_ref}.${seqid_qry1}.${seqid_qry2}.png" \
            "${seqid_ref}" \
            "${seqid_qry1}" \
            "${seqid_qry2}"
    }

    function main() {
        parse_args "$@"
        make_taskdir
        to_bed
        to_cds
        search_microsynteny
        output_besthit
        # make_besthit_bar
    }

    main "$@"
}