import os
import sys
import csv
from . import common

class Highlight:

    def __init__( self, seq, start, end, height, color, y, zorder ):
        self.height = height
        self.color = color
        self.pos = self.__convert_position2coord( seq, start, end, y )
        self.zorder = zorder
        
    def __convert_position2coord( self, seq, h_start, h_end, y ):
        x1 = seq.convert_position2xcoord( h_start )
        x2 = seq.convert_position2xcoord( h_end + 1 )
        y1 = seq.convert_position2ycoord( seq.height /2 + self.height / 2 + y )
        y2 = seq.convert_position2ycoord( seq.height /2 - self.height / 2 + y )
        x = [ x1, x2, x2, x1 ]
        y = [ y1, y1, y2, y2 ]
        return [ x, y ]

    def plot( self, ax ):
        if not self.pos:
            return
        ax.fill( self.pos[0], self.pos[1], color=self.color.color, alpha=self.color.alpha, lw=0, zorder=self.zorder )


def plot_highlight( seqs, ax, size, fn, alpha, thickness, n ):
    flag = False
    try:
        if not os.path.exists(fn):
            raise FileNotFoundError(f"Error: The file '{fn}' does not exist. Please check the file path.")
        
        if os.path.getsize(fn) == 0:
            raise ValueError(f"Error: The file '{fn}' is empty. Please provide a valid input file.")

        with open( fn , 'r' ) as file:
            buf = csv.reader(file, delimiter='\t')
            for lineno, row in enumerate(buf, start=1):
                if len(row) != 4:
                    print(f"[ERROR] Line {lineno}: Expected 4 columns, got {len(row)} columns.", file=sys.stderr )
                    sys.exit(1)
                seqid, start_pos, end_pos, color = row
                
                try:
                    start_pos = int(start_pos)
                    end_pos = int(end_pos)
                except ValueError:
                    print(f"[ERROR] Line {lineno}: start or end is not an integer: start='{start_pos}', end='{end_pos}'", file=sys.stderr)
                    sys.exit(1)

                color1 = common.Color( color, alpha )                
                if start_pos < end_pos :
                    h_lists = detect_index_for_highlights( seqid, start_pos, end_pos, seqs )
                else :
                    h_lists = detect_index_for_highlights( seqid, end_pos, start_pos, seqs )
                if n == 0 :
                    for i, j, start, end in h_lists:
                        if size.histograms[i] == 0:
                            continue
                        height = size.histograms[i] - seqs[i][j].height - size.margin_bw_seq_aspace - size.margin_bw_gspace_aspace
                        aninstance = Highlight( seqs[i][j], start, end, height, color1, size.histograms[i]/2, 1 )
                        aninstance.plot( ax )
                        flag = True
                else :
                    for i, j, start, end in h_lists:
                        height = seqs[i][j].height * thickness                        
                        zorder = 3 if ( thickness == 1 or alpha == 1 ) else 1
                        aninstance = Highlight( seqs[i][j], start, end, height, color1, 0, zorder )
                        aninstance.plot( ax )
                        if size.histograms[i] != 0 :
                            height = seqs[i][j].height
                            white_base = Highlight( seqs[i][j], start, end, height * 1.1, common.Color( "white", 1 ), size.histograms[i], 3 )
                            # white_base.plot( ax )
                            aninstance = Highlight( seqs[i][j], start, end, height, color1, size.histograms[i], 3 )
                            # aninstance.plot( ax )
                        flag = True

    except (FileNotFoundError, ValueError) as e:
        print(e, file=sys.stderr)
        return flag
    
    return flag


def output_parameters( flag, alpha, thickness ):
    if not flag:
        return
    print( '' )
    print( '## Highlights paramenters:' )
    print( '--h_thickness %.2f' % ( thickness ))
    print( '--h_alpha %s' % ( alpha ))


def output_parameters4sp( flag, alpha ):
    if not flag:
        return
    print( '' )
    print( '## Scatter plot Highlights paramenters:' )
    print( '--sp_h_alpha %s' % ( alpha ))

    
def detect_index_for_highlights( ID, start, end, seqs ):
    h_list = []
    for i in range( len( seqs )):
        for j in range( len( seqs[i] )):
            if str( ID ) != seqs[i][j].ID :
                continue
            if( seqs[i][j].name == 'BLANK' ):
                continue
            if( seqs[i][j].start <= start and end <= seqs[i][j].end ):
                h_list.append( [ i, j, start, end ] )
            elif( seqs[i][j].start < start and start <= seqs[i][j].end and seqs[i][j].end < end ):
                h_list.append( [ i, j, start, seqs[i][j].end ] )
            elif( start < seqs[i][j].start and seqs[i][j].start < end and end <= seqs[i][j].end ):
                h_list.append( [ i, j, seqs[i][j].start, end ] )
            elif( start < seqs[i][j].start and seqs[i][j].end < end ):
                h_list.append( [ i, j, seqs[i][j].start, seqs[i][j].end ] )
    return h_list
