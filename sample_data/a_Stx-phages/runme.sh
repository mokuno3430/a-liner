#!/bin/sh

prefix=Stx-prophage-regions 
a-liner --xlsx input/sequence_config.xlsx \
        --seq_thickness 0.75 --seq_color black --scale tick \
        --blastn input/blastn_outfmt6.txt \
        --min_identity 80 --alignment_alpha 0.35 --colormap 0 \
        --gb input/*.gb --gff3 input/*.gff --gene_thickness 6 \
        --feature_color_map input/color_legend_map.txt \
        --feature CDS pseudogene --gene_edge_color black \
        --scatter input/GCcontent_files/*_GCcontent_w1000-s100.txt \
        --scatter_min 30 --scatter_max 70 --scatter_space 0.6 --scatter_ylines 30 50 70 \
        --marker_color "#c33399" --background_color "white" \
        --highlight input/highlights.txt --h_alpha 0.5 --h_thickness 8 \
	    --out ${prefix} 1> ${prefix}.log
