import argparse
from matplotlib import colors as mcolors
from matplotlib import markers
import matplotlib.colors as clr
import os.path
import sys

def get_args( argv=None ):
    valid_marker_list = get_marker_styles_list()
    parser = argparse.ArgumentParser( formatter_class=argparse.MetavarTypeHelpFormatter, allow_abbrev=False )

    group_i = parser.add_mutually_exclusive_group(required=True)
    group_i.add_argument( '-i', '--input', help='File(s): sequence info for display. '
                          'Format: tab-delimited with columns [n, seq_ID, start(1-based), end(1-based), strand(+ or -), display name] ',
                          type=str, metavar='file')
    group_i.add_argument( '--xlsx', help='File: seq_info.xlsx', type=str, metavar='Excel')
    parser.add_argument( '--xlsx_sheet', help='sheet name of seq_info.xlsx', type=str, default=0 )
    group_i.add_argument( '--sep_input', help=argparse.SUPPRESS, type=str, nargs='*', metavar='file')
    parser.add_argument( '--out', help='Optional: prefix of PDF file (default: %(default)s).', type=str, default='out' )
    parser.add_argument( '--figure_size', help='Optional: figure size as [width height](inch). If width=0 → set to 6. If height=0 → auto(default: %(default)s).',
                         type=validate_float, nargs=2, default=[6,0], metavar=( 'width', 'height') )

    layout_group = parser.add_argument_group("Sequence layout options")
    layout_group.add_argument( '--seq_layout', help='Optional: sequence layout (default: %(default)s).',
                               choices=['left', 'center', 'right'], default='left', type=str )
    layout_group.add_argument( '--margin_bw_seqs', help='Optional: vertical margin between adjacent sequences. Default -1 means auto-adjust.',
                               type=float, default=-1, metavar='float')
    layout_group.add_argument( '--xlim_max', help='Optional: maximum x-axis coordinate for plotting (bp). Default -1 means auto-adjust.',
                               type=int, default=-1, metavar='int')
    layout_group.add_argument( '--left_margin', help='Optional: left side margin of the figure. Default -1 means auto-adjust. range 0.05-0.50.',
                               type=validate_left_margin, default=-1, metavar='float')

    seq_group = parser.add_argument_group("Sequence drawing options")
    seq_group.add_argument( '--seq_color', help='Optional: color of sequences (default %(default)s).', default='grey', type=validate_color, metavar='str' )
    seq_group.add_argument( '--seq_font_size', help='Optional: font size of sequence names (pt). Default 6.',
                            type=float, default=6, metavar='float')
    seq_group.add_argument( '--seq_thickness', help='Optional: thickness of sequence lines (pt) (default: %(default)s).', type=float, default=1.5, metavar='float' )

    scale_group = parser.add_argument_group("Sequence scale options")
    scale_group.add_argument( '--scale', help='Optional: how to display scale. '
                              '"legend" = show scale bar, "tick" = show axis ticks on sequences, "both" = show both (default: %(default)s).',
                              type=str, choices=['legend', 'tick', 'both'], default='legend')
    scale_group.add_argument( '--tick_width', help='Optional: scale width of axis (bp) (default -1 means auto).', type=int, default=-1, metavar='int' )
    scale_group.add_argument( '--tick_font_size', help='Optional: font size of ticks (pt) (default: %(default)s).', type=float, default=3, metavar='float' )

    alignment_files = parser.add_argument_group("Sequence alignment files")
    alignment_files.add_argument( '-a', '--alignment', help='File(s): custom alignment data. Format: tab-delimited with columns '
                                  '[seq_ID1, start1(1-based), end1(1-based), seq_ID2, start2(1-based), end2(1-based), identity(%%)].',
                                  type=str, nargs='*', default=[], metavar='file')
    alignment_files.add_argument( '--blastn', help="File(s): blastn output. Example: blastn -db ref.fa -query query.fa -out blastn.txt -outfmt 6",
                                  type=str, nargs='*', default=[], metavar='file')
    alignment_files.add_argument( '--lastz', help='File(s): lastz output. Example: lastz ref.fa query.fa --format=general --output=lastz.txt',
                                  type=str, nargs='*', default=[], metavar='file')
    alignment_files.add_argument( '--mummer', help='File(s): MUMmer show-coords output. Example: show-coords -H out.delta > show-coords.tsv',
                                  type=str, nargs='*', default=[], metavar='file')
    alignment_files.add_argument( '--minimap2', help='File(s): minimap2 PAF output. Example: minimap2 -c ref.fa query.fa > out.paf', type=str, nargs='*', default=[], metavar='file')

    alignment_group = parser.add_argument_group("Sequence alignment options")
    alignment_group.add_argument( '--min_identity', help='Optional: minimum sequence identity (%%). '
                                  'Alignments below this threshold will be ignored (default: %(default)s).',
                                  type=int, default=70, metavar='int')
    alignment_group.add_argument( '--min_alignment_len', help='Optional: minimum alignment length (bp). '
                                  'Alignments shorter than this will be ignored (default: %(default)s).',
                                  type=int, default=0, metavar='int')
    alignment_group.add_argument( '--alignment_alpha', help='Optional: transparency (alpha) of alignment coloring, '
                                  'range 0–1 (0 = fully transparent, 1 = opaque) (default: %(default)s).',
                                  type=validate_range, default=0.5, metavar='float')
    alignment_group.add_argument( '--colormap', help='Optional: colormap for sequence identity. '
                                  '0 = bone_r, 1 = hot_r, 2 = BuPu, 3 = YlOrRd, 4 = YlGnBu, 5 = rainbow (original) (default: %(default)s).',
                                  choices=[0, 1, 2, 3, 4, 5], default=5, type=int )
    alignment_group.add_argument( '--include_nonadjacent', help='Include alignments between non-adjacent sequences (default: only adjacent).',
                                  action='store_true', default=False )
    
    gene_files = parser.add_argument_group("Gene annotation files")
    gene_files.add_argument( '--gff3', help='File(s): gene annotation in GFF format.', type=str, nargs='*', default=[], metavar='gff3')
    gene_files.add_argument( '--gff_xlsx', help='File(s): GFF format in Excel files', type=str, nargs='*', default=[], metavar='Excel')
    gene_files.add_argument( '--gb', help='File(s): genbank format.', type=str, nargs='*', default=[], metavar='genbank')

    feature_legend_group = parser.add_argument_group("Gene / feature legend options")
    feature_legend_group.add_argument( '--feature_color_map', help='TSV file specifying feature labels, matching keywords, and colors. '
                             'Features matching the keywords will be colored accordingly and shown in the legend.', type=str, metavar='FILE')
    feature_legend_group.add_argument( '--feature_color_legend_font_size', help='Optional: font size of gene names (pt) (default: %(default)s).', type=float, default=5, metavar='float' )
    feature_legend_group.add_argument( '--feature_color_legend_marker_size', help='Optional: marker size (default: %(default)s).', default=5, type=float, metavar='float')
    feature_legend_group.add_argument( '--feature_color_legend_ncol', help='Optional: 0 means auto.', default=0, type=int, metavar='int')
    
    gene_group = parser.add_argument_group("Gene / feature drawing options")
    gene_group.add_argument( '--feature', help="Optional: GFF/GenBank feature types to draw (space-separated)(default: gene).", type=str, nargs='*', default={ 'gene' } )
    gene_group.add_argument( '--gene_thickness', help='Optional: relative thickness of gene arrows compared to seq_thickness (default: %(default)s).',
                             type=float, default=3, metavar='float' )
    gene_group.add_argument( '--gene_label_attr', help='Optional: attribute key used for feature labels (default: %(default)s).',
                             type=str, default='gene', metavar='str' )
    gene_group.add_argument( '--gene_font_size', help='Optional: font size of gene names (pt) (default: %(default)s).', type=float, default=3, metavar='float' )
    gene_group.add_argument( '--gene_font_rotation', help='Optional: rotation angle of gene names (degrees) (default: %(default)s).', type=float, default=75, metavar='float')
    gene_group.add_argument( '--gene_color', help='Optional: fill color of gene arrows (default: %(default)s).', type=validate_color, default='white', metavar='str')
    gene_group.add_argument( '--gene_edge_color', help='Optional: edge (outline) color of gene arrows (default: %(default)s) (no outline).', type=validate_color, default=None, metavar='str')

    highlight_group = parser.add_argument_group("Highlight options")
    highlight_group.add_argument( '--highlight', help='File(s): highlight regions. Format: tab-delimited with columns [seq_ID, start(1-based), end(1-based), color]',
                                  type=str, nargs='*', default=[], metavar='file')
    highlight_group.add_argument( '--h_alpha', help='Optional: transparency of highlights (0=transparent, 1=opaque) (default: %(default)s).', type=validate_range, default=0.3, metavar="float")
    highlight_group.add_argument( '--h_thickness', help='Optional: relative thickness of highlights compared to sequence thickness (default: %(default)s).',
                                  type=float, default=3.5, metavar='float')
    
    scatter_group = parser.add_argument_group("Scatter plot options")
    scatter_group.add_argument( '--scatter', help='File(s): scatterplot data. Format: tab-delimited with columns [seq_ID, position(1-based), value]',
                                type=str, nargs='*', default=[], metavar='file')
    scatter_group.add_argument( '--marker_color', help='Optional: marker color (default: %(default)s).', default='deeppink', type=validate_color, metavar='str' )
    scatter_group.add_argument( '--marker_size', help='Optional: marker size (default: %(default)s).', default=3, type=float, metavar='float')
    scatter_group.add_argument( '--marker_style', help=f"Optional: marker style. Valid choices: {', '.join(valid_marker_list)} (default: %(default)s).",
                                default='.', type=validate_marker_style, metavar='str')
    scatter_group.add_argument( '--scatter_space', help='Optional: relative height of scatterplot compared to alignment space (default: %(default)s). '
                                'For example, 0.8 means 80%% of the alignment height.', default=0.8, type=float, metavar='float')
    scatter_group.add_argument( '--scatter_min', help='Optional: minimum value of y-axis (default: %(default)s).', type=float, default=0, metavar='float')
    scatter_group.add_argument( '--scatter_max', help='Optional: maximum value of y-axis (default: %(default)s).', type=float, default=4, metavar='float')
    scatter_group.add_argument( '--scatter_ylines', help='Optional: add horizontal reference lines at the given y values (list of floats).',
                                type=validate_float, metavar="float", nargs='*', default=[])
    scatter_group.add_argument( '--background_color', help='Optional: background color of scatter plot (default: %(default)s).', default='whitesmoke', type=validate_color, metavar='str')
    scatter_group.add_argument( '--sp_highlight', help='File(s): highlight regions for scatter plot. Format: tab-delimited with columns [seq_ID, start(1-based), end(1-based), color]',
                                  type=str, nargs='*', default=[], metavar='file')
    scatter_group.add_argument( '--sp_h_alpha', help='Optional: transparency of highlights for scatter plot(default: %(default)s).', type=validate_range, default=0.3, metavar="float")

    
    # Version info
    parser.add_argument( '-v', '--version', action='version', version='%(prog)s v.1.1.0', default=True )

    if argv is None:
        return parser.parse_args()
    else:
        return parser.parse_args(argv)


class Size:
    def __init__( self, seqs, args ):
        length = self._cal_total_length( seqs )
        self.alignment_height = 0.8 if args.scale == 'legend' else 0.95
        self.alignments = [ self.alignment_height ] * len( seqs )
        self.alignments[-1] = 0
        self.histogram_height = self.alignment_height * args.scatter_space
        self.histograms = self._set_histogram_space( seqs, args.scatter )
        self.margin_bw_seqs = int( max( length )/200 ) if args.margin_bw_seqs == -1 else args.margin_bw_seqs
        self.xlim_max = self.cal_xlim_max( length, seqs, args.xlim_max )
        self.bottom_margin = self.__set_bottom_margin( args )
        self.top_margin = 0.5 if args.scale == 'legend' else 0.65
        self.ylim_max = self.__set_ymax( seqs )
        self.figsize_inch = []
        self.figsize_inch.append( args.figure_size[0] if args.figure_size[0] != 0 else 6 )
        self.figsize_inch.append( self.ylim_max * 0.5 if args.figure_size[1] == 0 else args.figure_size[1] )
        self.PT2INCH4Y = self.ylim_max /( self.figsize_inch[1] * 0.9 * 72 )
        self.PT2INCH4X = self.xlim_max /( self.figsize_inch[0] * 72 )
        self.seq_thickness = self._set_seq_thickness( args )
        self.margin_bw_gspace_aspace = 1.5 * self.PT2INCH4Y
        self.margin_bw_seq_aspace = 1.5 * self.PT2INCH4Y
        self.seq_layout = args.seq_layout
        if args.gene_thickness >= 1 and ( args.gff3 or args.gff_xlsx or args.gb ):
            self.margin_bw_seq_aspace = ( self.seq_thickness * ( args.gene_thickness - 1 )* 0.5 ) + 1.5 * self.PT2INCH4Y               
        if args.h_thickness >= 1 and args.highlight :
            if self.margin_bw_seq_aspace < ( self.seq_thickness * ( args.h_thickness - 1 )* 0.5 ) + 1.5 * self.PT2INCH4Y :
                self.margin_bw_seq_aspace = ( self.seq_thickness * ( args.h_thickness - 1 )* 0.5 ) + 1.5 * self.PT2INCH4Y
        self._set_scaffold_layout( seqs )


    def _set_seq_thickness( self, args ):
        default_value = self.PT2INCH4Y * 1.5
        seq_thickness = self.PT2INCH4Y * args.seq_thickness
        if seq_thickness >= 0.5 :
            print( f"Warnig: The value of '{args.seq_thickness}' specified for --seq_thickness option is too large. The default value of 1.5 is used instead.",
                   file=sys.stderr )
            return default_value
        if ( args.gff3 + args.gff_xlsx + args.gb ) and seq_thickness * args.gene_thickness >= 0.5 :
            print( f"Warnig: The value of '{args.seq_thickness}' specified for --seq_thickness option is too large. The default value of 1.5 is used instead.",
                   file=sys.stderr )
            return default_value
        if args.highlight and seq_thickness * args.h_thickness >= 0.5 :
            print( f"Warnig: The value of '{args.seq_thickness}' specified for --seq_thickness option is too large. The default value of 1.5 is used instead.",
                   file=sys.stderr )
            return default_value
        return seq_thickness

    def __set_bottom_margin( self, args ):
        bottom = 1.4
        if not( args.blastn + args.alignment + args.minimap2 + args.lastz + args.mummer ) :
            bottom = bottom - 0.2 if args.scale == 'tick' else bottom - 0.1
        if not( args.gff3 + args.gff_xlsx + args.gb ) or args.gene_font_size == 0:
            bottom -= 0.4
        return bottom
    
    def __set_ymax( self, seqs ):
        y = self.bottom_margin
        for i in range( len( seqs )):
            y += self.alignments[i] + self.histograms[i]
        y += self.top_margin
        return y

    def _cal_total_length( self, seqs ):
        length = []
        for i in range( len( seqs )):
            total = 0
            for j in range( len( seqs[i] )):
                total += seqs[i][j].length
            length.append( total )
        return length

    def cal_xlim_max( self, length, seqs, xlim_max ):
        sum_len = []
        for i in range( len( length ) ):
            sum_len.append( length[i] + self.margin_bw_seqs * ( len( seqs[i] ) - 1 ) )
        return max( sum_len ) if max( sum_len ) > xlim_max else xlim_max

    def _set_scaffold_layout( self, seqs ):
        length = self._cal_total_length( seqs )
        for i in range( len( length ) ):
            length[i] += self.margin_bw_seqs * ( len( seqs[i] ) - 1 )

        y = self.bottom_margin
        for i in range( len( seqs )):
            x = self.__set_xvalue( length[i] )
            for j in range( len( seqs[i] )):
                seqs[i][j].set_origins( x, y, self.seq_thickness )
                x += seqs[i][j].length + self.margin_bw_seqs
            y += self.alignments[i] + self.histograms[i]
        y += self.top_margin

    def __set_xvalue( self, length ):
        if( self.seq_layout == 'center' ):
            return int( ( self.xlim_max - length )/2 )
        if( self.seq_layout == 'left' ):
            return 0
        return self.xlim_max - length

    def _set_histogram_space( self, seqs, files ):
        histograms = [0] * len( seqs )
        for fn in files:
            if not os.path.isfile( fn ):
                print( f"Error: '{fn}' is not found.", file=sys.stderr )
                continue
            with open( fn , 'r' ) as file:
                for line in file:
                    buf = line.rstrip( '\n' ).split( '\t' )
                    i_list = detect_index_update( buf[0], int( buf[1] ), int( buf[1] ), seqs )
                    for i, j in i_list:
                        if( i == -1 ):
                            continue
                        histograms[ i ] = self.histogram_height
        return histograms

    def output_parameters( self, left_margin ):
        print( '' )
        print( '## Size paramenters' )
        print( '--figure_size %.2f %.2f' % ( self.figsize_inch[0], self.figsize_inch[1] ))
        print( '--seq_layout %s' % ( self.seq_layout ))
        print( '--margin_bw_seqs %d' % ( self.margin_bw_seqs ))
        print( '--xlim_max %d' % ( self.xlim_max ))
        print( '--left_margin %.2f' % ( left_margin ))


class Text:
    def __init__( self, label, size, x, y, ha, va, color='black' ):
        self.label = label
        self.size = size
        self.origin_x = x
        self.origin_y = y
        self.color = color
        self.ha = ha
        self.va = va

    def output( self, ax ):
        if self.label == 'BLANK':
            return
        if self.size == 0:
            return
        ax.text( self.origin_x, self.origin_y, self.label, fontsize = self.size, color = self.color, ha=self.ha, va=self.va )
        

class Color:
    replaced_color = 'grey'
    colors = dict(mcolors.BASE_COLORS, **mcolors.CSS4_COLORS)

    def __init__( self, color, alpha ):
        self.color, self.alpha = self.__evaluate_input( color, alpha )

    def __evaluate_input( self, color, alpha ):
        if clr.is_color_like( color ):
            return color, alpha
        else:
            print( f"Error: '{color}' is an invalid color.", file=sys.stderr )
            return Color.replaced_color, 0


def detect_index_update( ID, start, end, seqs ):
    i_list = []    
    for i in range( len( seqs )):
        for j in range( len( seqs[i] )):
            if str( ID ) != seqs[i][j].ID:
                continue
            if( seqs[i][j].name == 'BLANK' ):
                continue
            if( seqs[i][j].start <= start and end <= seqs[i][j].end ):
                i_list.append( [i, j] )
    return i_list


def validate_range( value ):
    try:
        value = float(value)
    except ValueError:
        raise argparse.ArgumentTypeError(f"{value} is not a valid float.")

    if not ( 0 <= value <= 1 ):
        raise argparse.ArgumentTypeError(f"Value must be in the range [0, 1]. You provided {value}." )
    return value

def validate_left_margin( value ):
    try:
        value = float(value)
    except ValueError:
        raise argparse.ArgumentTypeError(f"{value} is not a valid float.")

    if not ( 0.05 <= value <= 0.5 ):
        raise argparse.ArgumentTypeError(f"Value must be in the range [0.05, 0.50]. You provided {value}." )
    return value


def validate_float(value):
    try:
        return float(value)
    except ValueError:
        raise argparse.ArgumentTypeError(f"Invalid float value: '{value}'")


def validate_color( value ):
    try:
        mcolors.to_rgba( value )
        return value
    except ValueError:
        raise argparse.ArgumentTypeError(f"'{value}' is not a valid color. Use a named color or hex code like '#RRGGBB'.")

    
def validate_marker_style( value ):
    valid_markers = set( get_marker_styles_list() )
    if value not in valid_markers:
        raise argparse.ArgumentTypeError(
            f"Invalid marker style '{value}'. Valid options are: {', '.join(sorted(valid_markers))}"
        )
    return value


def get_marker_styles_list():
    valid_markers = ['.', ',', 'o', 'v', '^', '<', '>', '8', 's', 'p', '*', 'h', 'H', 'D', 'd', 'P', 'X' ]
    return sorted(set(valid_markers))
