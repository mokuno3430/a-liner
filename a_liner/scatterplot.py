import os.path
import sys
import csv
from . import common

def plot_scatterplot( seqs, ax, size, fn, args ):
    color = args.marker_color
    msize = args.marker_size
    m_style = args.marker_style
    min_y = args.scatter_min
    max_y = args.scatter_max
    flag = False
    if not os.path.isfile( fn ):
        return flag
    with open( fn , 'r' ) as file:
        x_list = []
        y_list = []
        for line in file:
            buf = line.rstrip( '\n' ).split( '\t' )
            i_list = common.detect_index_update( buf[0], int( buf[1] ), int( buf[1] ), seqs )
            for i, j in i_list:
                if( i == -1 ):
                    continue
                if( float( buf[2] ) < min_y or max_y < float( buf[2] )):
                    continue
                x = seqs[i][j].convert_position2xcoord( int(buf[1] ))
                y0 = seqs[i][j].origin_y + seqs[i][j].height + size.margin_bw_seq_aspace
                y1 = seqs[i][j].origin_y + size.histograms[i] - size.margin_bw_gspace_aspace
                y = ( (y1 - y0 ) * (float( buf[2] ) - min_y ) )/(max_y - min_y) + y0
                x_list.append( x )
                y_list.append( y )
    if( len( x_list ) != 0 ):
        ax.scatter( x_list, y_list, color=color, marker=m_style, s=msize, zorder=4, lw=0 )
        flag = True
    return flag

def plot_background( seqs, ax, size, args ):
    min_y = args.scatter_min
    max_y = args.scatter_max
    color = common.Color( args.background_color, 1 )    
    for i in range( len( seqs )):
        if( size.histograms[i] == 0 ):
            continue
        n = 0
        for j in range( len( seqs[i] )):
            if seqs[i][j].name == 'BLANK':
                continue
            x_start = seqs[i][j].origin_x
            x_end = seqs[i][j].origin_x + seqs[i][j].length
            y_bottom = seqs[i][j].origin_y + size.margin_bw_seq_aspace + seqs[i][j].height
            y_top = seqs[i][j].origin_y + size.histogram_height - size.margin_bw_gspace_aspace
            x = [ x_start, x_end, x_end, x_start ]
            y = [ y_bottom, y_bottom, y_top, y_top ]
            seqs[i][j].plot_line_on_histogram( ax, size.histograms[i] )
            ax.fill( x, y, color=color.color, lw=0, alpha=color.alpha, zorder=0 )

            for k in args.scatter_ylines:
                if k < min_y :
                    continue
                if max_y < k :
                    continue
                y = ( (y_top - y_bottom ) * ( k - min_y ) )/(max_y - min_y ) + y_bottom
                ax.hlines( y=y, xmin=seqs[i][j].origin_x, xmax=seqs[i][j].origin_x + seqs[i][j].length, lw=0.25, ls='--', color='grey', zorder=2 )

            tick_set = set( args.scatter_ylines )
            tick_set.update( [args.scatter_min, args.scatter_max] )
            x_pos = seqs[i][j].origin_x - 1 * 0.55 * size.PT2INCH4X
            if n != 0 :
                continue
            for t in tick_set :
                y = ( (y_top - y_bottom ) * ( t - min_y ) )/(max_y - min_y ) + y_bottom
                ax.text( x_pos, y, t, fontsize = args.tick_font_size, ha='right', va='center' )
            n += 1


def output_parameters( flag, args ):
    if not flag :
        return
    print( '' )
    print( '## Scatter plot paramenters:' )
    print( '--marker_color "%s"' % ( args.marker_color ))
    print( '--marker_size %.1f' % ( args.marker_size ))
    print( '--marker_style "%s"' % ( args.marker_style ))
    print( '--scatter_space %.1f' % ( args.scatter_space ))
    print( '--scatter_min %.1f' % ( args.scatter_min ))
    print( '--scatter_max %.1f' % ( args.scatter_max ))
    print( '--scatter_ylines %s' % ( " ".join(map( str, args.scatter_ylines))))
    print( '--background_color "%s"' % ( args.background_color ))
