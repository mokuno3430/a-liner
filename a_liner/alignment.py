import sys
import re
import csv
import os.path
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.colors as clr
from . import common

class Colormap:
    pallet = {
        'red':((0.00, 0.10, 0.10),
                (0.25, 0.10, 0.10),
                (0.40, 0.30, 0.30),
                (0.60, 1.00, 1.00),
                (0.80, 0.90, 0.90),
                (1.00, 0.70, 0.70)),
        'green':((0.00, 0.10, 0.10),
                (0.25, 0.60, 0.60),
                (0.40, 0.80, 0.80),
                (0.60, 0.75, 0.75),
                (0.80, 0.30, 0.30),
                (1.00, 0.15, 0.15)),
        'blue':((0.00, 0.40, 0.40),
                (0.25, 1.00, 1.00),
                (0.40, 0.25, 0.25),
                (0.60, 0.00, 0.00),
                (0.80, 0.05, 0.05),
                (1.00, 0.20, 0.20))
    }
    mycm = clr.LinearSegmentedColormap('original', pallet )

    cmaps = [ cm.bone_r, cm.hot_r, cm.BuPu, cm.YlOrRd, cm.YlGnBu, mycm ]
    cmap_list = [ 'bone_r', 'hot_r', 'BuPu', 'YlOrRd', 'YlGnBu', 'original' ]
    
    def __init__( self, min_identity, max_identity, cm, alpha ):
        self.min_identity = min_identity
        self.max_identity = max_identity
        self.alpha = alpha
        self.cm = cm

    def convert_identity2color( self, identity ):
        color_value=( identity - self.min_identity )/( self.max_identity - self.min_identity )
        align_color=Colormap.cmaps[ self.cm ]( color_value )
        return align_color

    def output_parameters( self, min_alignment_len, include_nonadjacent ):
        print( '' )
        print( '## Sequence alignment paramenters' )
        print( '--min_identity %d' % ( self.min_identity))
        print( '--min_alignment_len %d' % ( min_alignment_len ))
        print( '--alignment_alpha %.2f' % ( self.alpha ))
        print( '--colormap %d (%s)' % ( self.cm, Colormap.cmap_list[self.cm]  ))
        if include_nonadjacent :
            print( '--include_nonadjacent' )


class Colorbox:
    def __init__( self, size, heatmap, scale_flag, scalebar_length ):
        self.PT2INCH4X = size.PT2INCH4X
        self.PT2INCH4Y = size.PT2INCH4Y
        self.height = 8 * self.PT2INCH4Y
        self.width = size.xlim_max * 0.2
        self.tick_font_size = 4
        self.legend_font_size = 6
        self.origin_x = size.xlim_max * 0.75 - scalebar_length if scale_flag != 'tick' else size.xlim_max * 0.8
        self.origin_y = float( ( self.tick_font_size + 6 ) * self.PT2INCH4Y - self.height / 2 )
        self.min_identity = heatmap.min_identity
        self.max_identity = heatmap.max_identity
        identity_range = self.max_identity - self.min_identity
        self.SCALE = 5
        self.BIN = 0.2
        if identity_range <= 1:
            self.SCALE = 0.25
            self.BIN = 0.01
        elif identity_range <= 2:
            self.SCALE = 0.5
            self.BIN = 0.05
        elif identity_range <= 5:
            self.SCALE = 1
            self.BIN = 0.1
        elif identity_range <= 10:
            self.SCALE = 2
            self.BIN = 0.1
        self.cell_width = (self.width / ( identity_range ))*self.BIN

        
    def plot( self, ax, heatmap ):
        text_x = self.origin_x - 0.8 * self.legend_font_size * 0.55 * self.PT2INCH4X
        text_legend = common.Text( 'identity (%)', self.legend_font_size, text_x, self.origin_y + self.height / 2,  'right', 'center' )
        text_legend.output( ax )
        cell_originx = self.origin_x
        sum = self.min_identity
        for i in np.arange( self.min_identity, self.max_identity + self.BIN, self.BIN ):
            align_color=heatmap.convert_identity2color( i )
            colorbox_x = [ cell_originx, cell_originx + self.cell_width - 1, cell_originx + self.cell_width - 1, cell_originx ]
            colorbox_y = [ self.origin_y, self.origin_y, self.origin_y + self.height, self.origin_y + self.height ]
            cell_originx += self.cell_width
            ax.fill( colorbox_x, colorbox_y, color=align_color, alpha=heatmap.alpha, linewidth=0 )

        text_x = self.origin_x
        y = self.origin_y - 2 * self.PT2INCH4Y
        for i in np.arange( self.min_identity, self.max_identity + self.BIN, self.SCALE ):
            x = text_x + self.cell_width/2
            text_x += self.cell_width * self.SCALE / self.BIN
            if self.SCALE == 0.25 :
                num_legend = common.Text( '{:.02f}'.format( i ), self.tick_font_size, x, y, 'center', 'top' )
            elif self.SCALE == 0.5 :
                num_legend = common.Text( '{:.01f}'.format( i ), self.tick_font_size, x, y, 'center', 'top' )
            else :
                num_legend = common.Text( int( i ), self.tick_font_size, x, y, 'center', 'top' )
            num_legend.output( ax )

    def output_parameters( self ):
        print( '' )
        print( '##Colorbox paramenters:' )
        print( '  width: %.2f' % ( self.width ))
        print( '  height: %.2f' % ( self.height ))
        print( '  origin_x %.2f' % ( self.origin_x ))
        print( '  origin_y: %.2f' % ( self.origin_y ))


def convert_position2coord( A, B, A_start, A_end, B_start, B_end, margin, margin2, A_h_height, B_h_height ):
    x1 = A.convert_position2xcoord( A_start )
    x2 = A.convert_position2xcoord( A_end )
    x3 = B.convert_position2xcoord( B_end )
    x4 = B.convert_position2xcoord( B_start )
    x = [ x1, x2, x3, x4 ]

    if( A.origin_y < B.origin_y ):
        if A_h_height == 0 :
            y1 = A.convert_position2ycoord( margin + A.height )
        else :
            y1 = A.convert_position2ycoord( A.height + A_h_height + margin2 )
        y2 = B.convert_position2ycoord( margin * -1 )
    else:
        y1 = A.convert_position2ycoord( margin * -1 )
        if B_h_height == 0 :
            y2 = B.convert_position2ycoord( margin + B.height )
        else :
            y2 = B.convert_position2ycoord( B.height + B_h_height + margin2 )
    y = [ y1, y1, y2, y2 ]

    return x, y


def plot_alignments( ax, seqs, i_list, m_list, size, heatmap, qstart, qend, dstart, dend, color, include_nonadjacent ):
    flag = False
    for i, j in i_list :
        if i == -1 :
            continue
        for m, n in m_list:
            if m == -1 :
                continue
            if i == m :
                continue
            if not include_nonadjacent and abs( i - m ) != 1 :
                continue
            x, y = convert_position2coord( seqs[i][j], seqs[m][n], qstart, qend, dstart,dend,
                                           size.margin_bw_seq_aspace, size.margin_bw_gspace_aspace, size.histograms[i], size.histograms[m] )
            ax.fill( x, y, color=color, alpha=heatmap.alpha, lw=0 )
            flag = True
    return flag


def load_original( seqs, ax, heatmap, size, fn, min_identity, min_alignment_len, include_nonadjacent ):
    flag = False
    if not os.path.isfile( fn ):
        print( f"Error: --original: '{fn}' is not found.", file=sys.stderr )
        return flag
    with open( fn , 'r' ) as file:
        for line in file:
            buf = line.rstrip( '\n' ).split( '\t' )
            if( float( buf[6] ) < min_identity ):
                continue
            if( int( buf[2] ) - int( buf[1] ) + 1 < min_alignment_len ):
                continue
            color = heatmap.convert_identity2color( float( buf[6] ))
            i_list = common.detect_index_update( buf[0], int( buf[1] ), int( buf[2] ), seqs )
            if( int( buf[4] ) < int( buf[5] ) ):
                m_list = common.detect_index_update( buf[3], int( buf[4] ), int( buf[5] ), seqs )
                flag = plot_alignments( ax, seqs, i_list, m_list, size, heatmap, int( buf[1] ), int( buf[2] ), int( buf[4] ), int( buf[5] )+1, color, include_nonadjacent ) or flag
            else:
                m_list = common.detect_index_update( buf[3], int( buf[5] ), int( buf[4] ), seqs )
                flag = plot_alignments( ax, seqs, i_list, m_list, size, heatmap, int( buf[1] ), int( buf[2] ), int( buf[4] )+1, int( buf[5] ), color, include_nonadjacent ) or flag
    return flag

    
def load_blastn( seqs, ax, heatmap, size, fn, min_identity, min_alignment_len, include_nonadjacent ):
    flag = False
    if not os.path.isfile( fn ):
        print( f"Error: --blastn: '{fn}' is not found.", file=sys.stderr )
        return flag
    with open( fn , 'r' ) as file:
        for line in file:
            buf = line.rstrip( '\n' ).split( '\t' )
            if( float( buf[2] ) < min_identity ):
                continue
            if( int( buf[7] ) - int( buf[6] ) + 1 < min_alignment_len ):
                continue
            color = heatmap.convert_identity2color( float( buf[2] ))
            i_list = common.detect_index_update( buf[0], int( buf[6] ), int( buf[7] ), seqs )
            if( int( buf[8] ) < int( buf[9] ) ):
                m_list = common.detect_index_update( buf[1], int( buf[8] ), int( buf[9] ), seqs )
                flag = plot_alignments( ax, seqs, i_list, m_list, size, heatmap, int( buf[6] ), int( buf[7] )+1, int( buf[8] ), int( buf[9] )+1, color, include_nonadjacent ) or flag
            else:
                m_list = common.detect_index_update( buf[1], int( buf[9] ), int( buf[8] ), seqs )
                flag = plot_alignments( ax, seqs, i_list, m_list, size, heatmap, int( buf[6] ), int( buf[7] )+1, int( buf[8] )+1, int( buf[9] ), color, include_nonadjacent ) or flag
    return flag


def load_lastz( seqs, ax, heatmap, size, fn, min_identity, min_alignment_len, include_nonadjacent ):
    flag = False
    if not os.path.isfile( fn ):
        print( f"Error: --lastz: '{fn}' is not found.", file=sys.stderr )
        return flag
    with open( fn , 'r' ) as file:
        for line in file:
            if( line[0:1] == "#" ):
                continue
            buf = line.rstrip( '\n' ).split( '\t' )
            s_start = int( buf[4] ) + 1
            if( float( buf[12][:-1] ) < min_identity ):
                continue
            if( int( buf[5] ) - s_start + 1 < min_alignment_len ):
                continue
            i_list = common.detect_index_update( buf[1], s_start, int( buf[5] ), seqs )
            if( buf[7] == '+' ):
                s_pos = int( buf[9] ) + 1
                e_pos = int( buf[10] )
                m_list = common.detect_index_update( buf[6], s_pos, e_pos, seqs )
            else:
                s_pos = int( buf[8] ) - int( buf[9] )
                e_pos = int( buf[8] ) - int( buf[10] )
                m_list = common.detect_index_update( buf[6], e_pos, s_pos, seqs )
                s_pos += 1
            color = heatmap.convert_identity2color( float( buf[12][:-1] ))
            flag = plot_alignments( ax, seqs, i_list, m_list, size, heatmap, s_start, int( buf[5] )+1, s_pos, e_pos + 1, color, include_nonadjacent ) or flag 
    return flag


def load_mummer( seqs, ax, heatmap, size, fn, min_identity, min_alignment_len, include_nonadjacent ):
    flag = False
    if not os.path.isfile( fn ):
        print( f"Error: --mummer: '{fn}' is not found.", file=sys.stderr )
        return flag
    with open( fn , 'r' ) as file:
        for line in file:
            buf = line.rstrip( '\n' ).split( )
            if( float( buf[9] ) < min_identity ):
                continue
            if( int( buf[1] ) - int( buf[0] ) + 1 < min_alignment_len ):
                continue
            i_list = common.detect_index_update( buf[11], int( buf[0] ), int( buf[1] ), seqs )
            if int( buf[3] ) < int( buf[4] ):
                q_start = int( buf[3] )
                q_end = int( buf[4] )
            else :
                q_start = int( buf[4] )
                q_end = int( buf[3] )
            m_list = common.detect_index_update( buf[12], q_start, q_end, seqs )
            color = heatmap.convert_identity2color( float( buf[9] ))
            if int( buf[3] ) < int( buf[4] ) :
                flag = plot_alignments( ax, seqs, i_list, m_list, size, heatmap, int( buf[0] ), int( buf[1] )+1, q_start, q_end + 1, color, include_nonadjacent ) or flag
            else:
                flag = plot_alignments( ax, seqs, i_list, m_list, size, heatmap, int( buf[0] ), int( buf[1] )+1, q_end + 1, q_start, color, include_nonadjacent ) or flag
    return flag


def func_cal_identity_minimap2( fields ):
    # Parse CIGAR from optional fields    
    cigar = ''
    for f in fields[12:]:
        if f.startswith("cg:Z:"):
            cigar = f[5:]        
            break

    # Calculate alignment match length and number of gap openings from CIGAR string.
    alignment_match = 0
    gapopen = 0
    prev_op = None
    identity = 0
    # Parse CIGAR string
    for length, op in re.findall(r'(\d+)([MIDNSHP=X])', cigar):
        length = int(length)
        if op in ('M', '=', 'X'):
            alignment_match += length
        elif op in ('I', 'D'):
            if prev_op != op:
                gapopen += 1
        prev_op = op

    matches = int( fields[9] )
    identity = ( matches / alignment_match ) * 100 if alignment_match > 0 else 0
    return identity


def load_minimap2( seqs, ax, heatmap, size, fn, min_identity, min_alignment_len, include_nonadjacent ):
    flag = False
    if not os.path.isfile( fn ):
        print( f"Error: --minimap2: '{fn}' is not found.", file=sys.stderr )
        return flag
    with open( fn , 'r' ) as file:
        for line in file:
            buf = line.rstrip( '\n' ).split( )
            # Extract basic fields
            qseqid, qlen, qstart, qend, strand, sseqid, slen, sstart, send, matches, alnlen = buf[0:11]
            qstart = int( qstart ) + 1
            sstart = int( sstart ) + 1
            if( int( qend ) - int( qstart ) + 1 < min_alignment_len ):
                continue
            if( int( send ) - int( sstart ) + 1 < min_alignment_len ):
                continue            
            identity = func_cal_identity_minimap2( buf )
            if identity < min_identity:
                continue
            color = heatmap.convert_identity2color( identity )
            i_list = common.detect_index_update( qseqid, int( qstart ), int( qend ), seqs )
            m_list = common.detect_index_update( sseqid, int( sstart ), int( send ), seqs )
            if strand == '+' :
                flag = plot_alignments( ax, seqs, i_list, m_list, size, heatmap, int( qstart ), int( qend )+1, int( sstart ), int( send )+1, color, include_nonadjacent ) or flag 
            else :
                flag = plot_alignments( ax, seqs, i_list, m_list, size, heatmap, int( qstart ), int( qend )+1, int( send )+1, int( sstart ), color, include_nonadjacent ) or flag
    return flag
