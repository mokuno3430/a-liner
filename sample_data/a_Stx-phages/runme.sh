#!/bin/sh

prefix=Stx-phage_loci
a-liner --xlsx seq_info.xlsx \
        --blastn blastn_outfmt6.txt \
        --seq_color black \
        --gff3 gff_files/*.gff3 \
        --gff_xlsx gff_files/cladeI_10290_gff.xlsx \
        --feature CDS pseudogene \
        --gene_edge_color black \
        --scale tick \
        --highlight highlights.txt \
        --h_alpha 0.5 \
        --h_thickness 4 \
        --colormap 0 \
        --left_margin 0.25 \
	    --out ${prefix} 1> log_${prefix}
