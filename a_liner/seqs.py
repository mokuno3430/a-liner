import matplotlib
import sys
import os.path
import csv
import math
import pandas as pd
from . import common

class Scaffold:
    def __init__( self, order, ID, start, end, strand, name, color, thickness ):
        self.order = order
        self.ID = ID
        if start < end :
            self.start, self.end = start, end
        else :
            self.start, self.end = end, start
        self.strand = strand
        self.name = name if name != 'nan' else ''
        self.length = self.end - self.start + 1
        self.origin_x = 0
        self.origin_y = 0
        self.height = 0
        self.color = common.Color( color, 1 )
        
    def set_origins( self, x, y, height ):
        self.origin_x = x
        self.origin_y = y
        self.height = height

    def convert_position2xcoord( self, pos ):
        if( self.strand == '+' ):
            x = self.origin_x + pos - self.start
        else:
            x = self.origin_x + self.length - pos + self.start
        return x

    def convert_position2ycoord( self, pos ):
        return self.origin_y + pos

    def plot_line( self, ax ):
        if self.height == 0 :
            return
        if self.name == 'BLANK':
            return
        x = [self.origin_x, self.origin_x + self.length, self.origin_x + self.length, self.origin_x ]
        y = [self.origin_y, self.origin_y, self.origin_y + self.height , self.origin_y + self.height ]
        ax.fill( x, y, color=self.color.color, lw=0, alpha=self.color.alpha, zorder=2 )

    def plot_line_on_histogram( self, ax, h_height ):
        if self.height == 0 :
            return
        if self.name == 'BLANK':
            return
        x = [self.origin_x, self.origin_x + self.length, self.origin_x + self.length, self.origin_x ]
        y = [self.origin_y + h_height, self.origin_y + h_height, self.origin_y + self.height + h_height, self.origin_y + self.height + h_height ]
        ax.fill( x, y, color=self.color.color, lw=0, alpha=self.color.alpha, zorder=2 )


def input_scaffold_tsv( scf_files, color, thickness ):
    seqs = []
    for fn in scf_files:
        if not os.path.isfile( fn ):
            print( f"Error: '{fn}' is not found.", file=sys.stderr )
            continue
        with open(fn, 'r') as file:
            buf = csv.reader(file, delimiter='\t')
            n = 0
            seqlist = []
            for lineno, row in enumerate(buf, start=1):
                if len(row) != 5:
                    print(f"Error {fn} Line {lineno}: Expected 5 columns, got {len(row)} columns.")
                    sys.exit(1)
                scaffold, scaf_start, scaf_end, scaf_strand, scaf_name = row
                # int check
                try:
                    scaf_start = int(scaf_start)
                    scaf_end = int(scaf_end)
                except ValueError:
                    print(f"Error {fn} Line {lineno}: seq_start or seq_end is not an integer: start='{scaf_start}', end='{scaf_end}'", file=sys.stderr)
                    sys.exit(1)
                # strand check
                if scaf_strand not in ('+', '-'):
                    print(f"Error {fn} Line {lineno}: seq_strand must be '+' or '-': got '{scaf_strand}'", file=sys.stderr )
                    sys.exit(1)
                aninstance = Scaffold( n, scaffold, scaf_start, scaf_end, scaf_strand, scaf_name, color, thickness )
                seqlist.append( aninstance )
                n += 1
            seqs.append( seqlist )
    return seqs


def validate_n_column(df, xlsx):
    try:
        n_values = df['n'].astype(int).tolist()
    except ValueError:
        print(f"Error {xlsx}: column A must contain only integers (0,1,2,...).", file=sys.stderr)
        sys.exit(1)

    unique_n = sorted(df['n'].unique())
    expected = list(range(len(unique_n)))

    if unique_n != expected:
        print(f"Error {xlsx}: column A must contain consecutive integers starting from 0.\n"
              f"  Found: {unique_n}\n"
              f"  Expected: {expected}\n"
              f"Please fix the 'n' column in your Excel file.", file=sys.stderr)
        sys.exit(1)
    return pd.Series(n_values, index=df.index)


def input_scaffold_xlsx( xlsx, sheet_name, color, thickness ):
    try:
        df = pd.read_excel(xlsx, sheet_name=sheet_name)
    except ValueError:
        print(f"Error: Sheet '{sheet_name}' not found in {xlsx}.", file=sys.stderr)
        sys.exit(1)
    
    df['n'] = validate_n_column(df, xlsx)

    max_value = df['n'].max()
    counts = [ 0 ] * ( max_value + 1 )
    seqs = [ [] for _ in range(max_value + 1) ]

    for row in df.itertuples(index=False):
        # int check
        try:
            scaf_start = int( row.start )
            scaf_end = int( row.end )
        except ValueError:
            print(f"Error {xlsx}: seq_start (column C) or seq_end (column D) is not an integer: start='{row.start}', end='{row.end}'", file=sys.stderr)
            sys.exit(1)
        # strand check
        if row.strand not in ('+', '-'):
            print(f"Error {xlsx}: seq_strand (column E) must be '+' or '-': got '{row.strand}'", file=sys.stderr )
            sys.exit(1)
        aninstance = Scaffold( counts[row.n], str(row.ID), row.start, row.end, row.strand, str(row.name), color, thickness )
        counts[ row.n ] += 1
        seqs[ row.n ].append( aninstance )
    return seqs


def plot_scaffolds( ax, seqs ):
    for i in range( len( seqs )):
        for j in range( len( seqs[i] )):
            seqs[i][j].plot_line( ax )

            
def plot_scaffold_names( ax, seqs, args, size ):
    histograms, PT2INCH4Y, PT2INCH4X = size.histograms, size.PT2INCH4Y, size.PT2INCH4X
    h_flag = 0
    for i in range( len( histograms ) ):
        h_flag += histograms[i]
    
    ratio = args.gene_thickness if args.gff3 or args.gff_xlsx or args.gb else 1 
    s_flag = 0
    for i in range( len( seqs )):
        s_flag += len( seqs[i] ) - 1

    est_left_space = 0.05
    for i in range( len( seqs )):
        for j in range( len( seqs[i] )):
            left_space = 0.05
            if args.seq_font_size == 0 :
                continue
            if seqs[i][j].name == 'BLANK' :
                continue
            if s_flag != 0 :
                x = seqs[i][j].origin_x
                if args.scale == 'legend' :
                    shift = 2 * PT2INCH4Y
                else:
                    shift = ( args.tick_font_size + 3 ) * PT2INCH4Y
                y = seqs[i][j].origin_y + histograms[i] + seqs[i][j].height * ( 1 + ratio ) * 0.5 + shift
                aninstance = common.Text( seqs[i][j].name, args.seq_font_size, x, y, 'left', 'bottom' )
            else :
                x = seqs[i][j].origin_x - 0.5 * args.seq_font_size * 0.55 * PT2INCH4X
                if h_flag != 0 : 
                    x -= 5 * args.tick_font_size * 0.55 * PT2INCH4X
                tmp = (( len( seqs[i][j].name ) + 0.5 )* args.seq_font_size * 0.55 )/( 72 * 0.9 * size.figsize_inch[0] )
                tmp2 = ( 0.9 - tmp )* seqs[i][j].origin_x / size.xlim_max
                left_space += tmp - tmp2
                y = seqs[i][j].origin_y + args.seq_font_size * 0.08 * PT2INCH4Y + histograms[i] * 0.5
                aninstance = common.Text( seqs[i][j].name, args.seq_font_size, x, y, 'right', 'center' )
            aninstance.output( ax )    
            est_left_space = left_space if est_left_space < left_space else est_left_space
    return est_left_space
    

def output_parameters( args ):
    print( '' )
    print( '## Sequence lines paramenters' )    
    print( '--seq_color "%s"' % ( args.seq_color ))
    print( '--seq_font_size %.1f' % ( args.seq_font_size ))
    print( '--seq_thickness %.2f' % ( args.seq_thickness ))


class Scalebar:
    def __init__( self, size, seqs, args ):
        self.scale_type = args.scale
        self.width = self._set_width( size.xlim_max, seqs, args.tick_width )
        self.mr, self.unit = self._set_unit_for_tick( seqs )
        self.bar_color = 'black'
        self.legend = self._convert_for_display_w_unit( self.width )
        self.tick_font_size = args.tick_font_size
        self.legend_font_size = 6
        self.PT2INCH4Y = size.PT2INCH4Y
        self.origin_x, self.origin_y = self._set_origin( size )
        self.label = common.Text( self.legend, self.legend_font_size, self.origin_x + self.width/2, self.origin_y - 3 * self.PT2INCH4Y, 'center', 'top' )
        
    def _set_origin( self, size ):
        x = size.xlim_max - self.width
        y = ( self.legend_font_size + 3 )* self.PT2INCH4Y
        return x, y

    def _set_width( self, xlim_max, seqs, user_width ):
        max_seq_len = self._set_max_seq_len( seqs )
        i = xlim_max
        if user_width < xlim_max * 0.5 and user_width != -1 :
            return user_width
        if user_width >= xlim_max * 0.5:
            print( f"Warnig: The value of '{user_width}' in the --tick-width option is too large. The optimal value is estimated automatically. ",
                   file=sys.stderr )
        width = self._cal_width_cand( i )
        while width > max_seq_len:
            i /= 2
            width = self._cal_width_cand( i )
        return width

    def _cal_width_cand( self, i ):
        if math.pow( 10, int( math.log10( i )))/ i <= 0.2 :
            return math.pow( 10, int( math.log10( i )))
        if math.pow( 10, int( math.log10( i )))/ i <= 0.4 :
            return math.pow( 10, int( math.log10( i )))/2
        return math.pow( 10, int( math.log10( i )))/5

    def _set_max_seq_len( self, seqs ):
        max_seq_len = 0
        for i in range( len( seqs )):            
            for j in range( len(seqs[i] )):
                if seqs[i][j].name == 'BLANK' :
                    continue
                if max_seq_len >= seqs[i][j].length :
                    continue
                max_seq_len = seqs[i][j].length                
        return max_seq_len

    def _set_seq_max_end_value( self, seqs ):
        max_value = 0
        for i in range( len( seqs )):            
            for j in range( len(seqs[i] )):
                if seqs[i][j].name == 'BLANK' :
                    continue
                if max_value >= seqs[i][j].end :
                    continue
                max_value = seqs[i][j].end                    
        return max_value
                
    def _set_unit_for_tick( self, seqs ):
        max_value = self.width
        if 9 <= math.log10( max_value ):
            mr = 10 ** 9
            unit = 'Gbp'
        elif 5 <= math.log10( max_value ):
            mr = 10 ** 6
            unit = 'Mbp'
        elif 2 <= math.log10( max_value ):
            mr = 10 ** 3
            unit = 'kbp'
        else:
            mr = 1
            unit = 'bp'
        return mr, unit

    def _convert_for_display_w_unit( self, i ):
        if self.width /self.mr < 1 :
            text = "{0:0.1f} {1}".format( i/(self.mr ), self.unit )
        else :
            text="{0:,} {1}".format( int( i/self.mr ), self.unit )
        return text
        
    def _convert_for_display_wo_unit( self, i ):
        if self.width /self.mr < 1 :
            text = "{0:0.1f}".format( i/(self.mr ) )
        else :
            text="{0:,}".format( int( i/self.mr ) )
        return text

    def plot_scale( self, ax, seqs, ratio, histograms ):
        if self.scale_type == 'legend' or self.scale_type == 'both' :
            self._plot_scale_legend( ax )
        if self.scale_type == 'tick' or self.scale_type == 'both' :
            self._plot_scale_tick( ax, seqs, ratio, histograms )

    def _plot_scale_legend( self, ax ):
        ax.plot( [ self.origin_x, self.origin_x + self.width ], [ self.origin_y, self.origin_y ], color=self.bar_color, lw=0.5 )
        self.label.output( ax )
    
    def _plot_scale_tick( self, ax, seqs, ratio, histograms ):
        func = {'tick' : self._convert_for_display_w_unit, 'both' : self._convert_for_display_wo_unit }
        for i in range( len( seqs )):
            for j in range( len( seqs[i] )):
                if seqs[i][j].name == 'BLANK' :
                    continue
                start = int( math.ceil( seqs[i][j].start / self.width ) * self.width )
                if histograms[i] == 0 :
                    baseline = ( 1 + ratio )/2 * seqs[i][j].height
                else :
                    baseline = seqs[i][j].height + histograms[i]
                y1 = seqs[i][j].convert_position2ycoord( baseline + 0.7 * self.PT2INCH4Y )
                y2 = seqs[i][j].convert_position2ycoord( baseline + 2.0 * self.PT2INCH4Y )
                y3 = seqs[i][j].convert_position2ycoord( baseline + 2.1 * self.PT2INCH4Y )
                for m in range(start, seqs[i][j].end + 1, int( self.width ) ):
                    x = seqs[i][j].convert_position2xcoord( m + 0.5 )
                    ax.plot( [x, x], [y1, y2], color=self.bar_color, lw=0.25 )
                    text = func[ self.scale_type ]( m )
                    aninstance = common.Text( text, self.tick_font_size, x, y3, 'center', 'bottom' )
                    aninstance.output( ax )

    def output_parameters( self ):
        print( '' )
        print( '## Scalebar paramenters' )
        print( '--scale %s' % ( self.scale_type ))
        print( '--tick_width %d' % ( self.width ))
        print( '--tick_font_size %.2f' % ( self.tick_font_size ))
