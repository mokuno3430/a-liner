import sys
import os.path
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from . import alignment, common, genes, scatterplot, highlight
from . import seqs as scf


def func_set_axes( ax, size ):
    ##set data range
    ax.set_xlim(0, size.xlim_max)
    ax.set_ylim(0, size.ylim_max)

    ##non-display axis
    for spine in ["right", "left", "bottom", "top"]:
        ax.spines[ spine ].set_color("none")
    ax.set_yticklabels([])
    ax.set_xticklabels([])
    ax.tick_params( length=0 )


def func_print_messages( ):
    print( ' '.join( sys.argv ))

    
def main():
    args = common.get_args()
    func_print_messages()
    if args.input:
        seqs = scf.input_scaffold_tsv( args.input, args.seq_color, args.seq_thickness )
    elif args.xlsx:
        seqs = scf.input_scaffold_xlsx( args.xlsx, args.xlsx_sheet, args.seq_color, args.seq_thickness )
    if not seqs:
        print( f"Error: no seqs!", file=sys.stderr )
        exit()
        
    size=common.Size( seqs, args )

    fig = plt.figure( figsize=size.figsize_inch )
    ax = fig.add_subplot(111)
    func_set_axes( ax, size )
    
    scf.plot_scaffolds( ax, seqs )
    est_left_space =  scf.plot_scaffold_names( ax, seqs, args, size )
    if args.left_margin != -1 :
        est_left_space = args.left_margin
    size.output_parameters( est_left_space )
    scf.output_parameters( args )
        
    ##plot scalebar or ticks
    scalebar = scf.Scalebar( size, seqs, args )
    ratio = args.gene_thickness if args.gff3 or args.gff_xlsx or args.gb else 1
    scalebar.plot_scale( ax, seqs, ratio, size.histograms )
    scalebar.output_parameters()

    ##set colormap
    heatmap = alignment.Colormap( args.min_identity, 100, args.colormap, args.alignment_alpha )

    ##plot alignment
    input_formats = [ args.alignment, args.blastn, args.lastz, args.mummer, args.minimap2 ]
    func_plot_alignmment = [ alignment.load_original, alignment.load_blastn, alignment.load_lastz, alignment.load_mummer, alignment.load_minimap2 ]
    flag = False
    for files, func_plot in zip( input_formats, func_plot_alignmment ):
        for fn in files:
            flag = func_plot( seqs, ax, heatmap, size, fn, args.min_identity, args.min_alignment_len, args.include_nonadjacent ) or flag

    ##set and plot colormap legend
    if flag :
        heatmap.output_parameters( args.min_alignment_len, args.include_nonadjacent )
        heatmap_legend = alignment.Colorbox( size, heatmap, args.scale, scalebar.width )
        heatmap_legend.plot( ax, heatmap )
            
    ##plot genes
    input_formats = [ args.gff3, args.gff_xlsx, args.gb ]
    func_plot_genes = [ genes.plot_genes_from_gff, genes.plot_genes_from_gff_excel, genes.plot_genes_from_gb ]
    flag = False
    for files, func_plot in zip( input_formats, func_plot_genes ):
        for fn in files:
            flag = func_plot( seqs, ax, size, fn, args ) or flag
    genes.output_genes_parameters( args, flag )

    ## plot highlight
    flag = False
    for fn in args.highlight:
        flag = highlight.plot_highlight( seqs, ax, size, fn, args.h_alpha, args.h_thickness, 0 ) or flag
    highlight.output_parameters( flag, args.h_alpha, args.h_thickness )

    flag = False
    for fn in args.sp_highlight:
        flag = highlight.plot_highlight( seqs, ax, size, fn, args.sp_h_alpha, args.h_thickness, 1 ) or flag
    highlight.output_parameters4sp( flag, args.sp_h_alpha )

    
    ##plot scatter
    flag = False
    scatterplot.plot_background( seqs, ax, size, args )
    for fn in args.scatter:
        flag = scatterplot.plot_scatterplot( seqs, ax, size, fn, args ) or flag
    scatterplot.output_parameters( flag, args )


    pdf_file = args.out + '.pdf'
    with PdfPages(pdf_file) as pp:
        if est_left_space < 0.5:
            fig.subplots_adjust(left=est_left_space, right=0.95, bottom=0.05, top=0.95)
            bbox = None
        else:
            print("Note: sequence name is long or font size is large → switched to auto-estimated margin", file=sys.stderr)
            print("      Figure width may exceed the specified --figure_size.", file=sys.stderr)
            bbox = 'tight'
        pp.savefig(fig, bbox_inches=bbox)
    plt.clf()
    print(f"\n{pdf_file} has been generated.")
    
