#!/bin/sh

prefix=ostrich-emu_sex-chromosomes
a-liner -i seq_info/lyla_emu-chrW.txt seq_info/lyla_emu-chrZ.txt \
        seq_info/lyla_ostrich-chrZ.txt seq_info/lyla_ostrich-chrW.txt \
        --seq_thickness 0.5 --seq_color '#3b2a2a' \
        --figure_size 6 2.8 \
        --scale both --tick_width 5000000 \
        --minimap2 minimap2/*_shift.paf \
	    --min_alignment_len 10000 --min_identity 85 \
        --highlight highlights.txt \
        --h_alpha 1 --h_thickness 1 \
        --sp_highlight highlights.txt \
        --scatter est-copy/*_est-copy.txt \
        --marker_color '#515d66' \
        --background_color '#fbfcfb' \
        --scatter_ylines 0 1 2 3 --scatter_max 3 \
        --out ${prefix} 1> log_${prefix}
