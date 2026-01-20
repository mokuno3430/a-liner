import os.path
import sys
import csv
from matplotlib.colors import to_rgb
import pandas as pd
from Bio import SeqIO
from . import common


class Gene:
    def __init__( self, scaffold, g_start, g_end, g_strand, height, tail_length, color, edge_color ):
        self.start, self.end, self.strand = self.__convert_position2coord( scaffold, g_start, g_end, g_strand )
        self.y_mid_coord = scaffold.convert_position2ycoord( scaffold.height /2 )
        self.width = self.end - self.start
        self.height = height
        self.color = common.Color( color, 1 )
        self.edge_color = edge_color
        self.x, self.y  = self.__set_gene_coord( tail_length )

    def __set_gene_coord( self, tail_length ):
        x1 = self.start
        x3 = self.end

        y1 = self.y_mid_coord - self.height / 2
        y2 = self.y_mid_coord
        y3 = self.y_mid_coord + self.height / 2

        if( self.strand == '+' ): # => or >
            x2 = self.end - tail_length if self.width > tail_length else self.start
            x = [ x1, x2, x3, x2, x1 ]
            y = [ y1, y1, y2, y3, y3 ]
        elif( self.strand == '-' ): # <= or <
            x2 = self.start + tail_length if self.width > tail_length else self.end
            x = [ x1, x2, x3, x3, x2 ]
            y = [ y2, y1, y1, y3, y3 ]
        else:
            x = [ x1, x3, x3, x1 ]
            y = [ y1, y1, y3, y3 ]
        return x, y

    def __convert_position2coord( self, scaffold, g_start, g_end, g_strand):
        if( scaffold.strand == '+' ):
            start = scaffold.convert_position2xcoord( g_start )
            end = scaffold.convert_position2xcoord( g_end )
            return start, end, g_strand
        else:
            start = scaffold.convert_position2xcoord( g_end )
            end = scaffold.convert_position2xcoord( g_start )
            if( g_strand == '+' ):
                return  start, end, '-'
            elif( g_strand == '-' ):
                return  start, end, '+'
            else:
                return start, end, g_strand

    def plot( self, ax ):        
        if self.edge_color is not None :
            edge_color, lw = self.edge_color, 0.15
        elif to_rgb(self.color.color) == (1.0, 1.0, 1.0):
            edge_color, lw = 'black', 0.15
        else:
            edge_color, lw = None, 0
        ax.fill( self.x, self.y, color=self.color.color, alpha=self.color.alpha, lw=lw, edgecolor=edge_color, zorder=4 )


def func_get_gene_name_from_gff( attributes, gene_label_attr ):
    gene_name = ''
    for attr in attributes.split(';'):
        if '=' in attr:
            key, value = attr.split('=', 1)
            if key == gene_label_attr:
                gene_name = value
                break
    return gene_name


def func_plot_gene_name( seq, ax, start, end, gene_name, gene_font_size, height, rotation, PT2INCH4Y ):
    if gene_font_size == 0 or gene_name == '' :
        return
    text_x = seq.convert_position2xcoord( ( start + end )/2 )
    text_y = seq.convert_position2ycoord( -1 * height/2 - PT2INCH4Y )
    ax.text( text_x, text_y, gene_name, fontsize = gene_font_size, color = 'black', ha='center', va='top', rotation=rotation, style='italic' )


def deal_gff( seqs, ax, size, args, buf, color ):
    feature_set = args.feature
    # Only extract 'gene' or desired features    
    if len( buf ) == 10:
        seqid, source, type_, start, end, score, strand, phase, attributes, color = buf
    elif len( buf ) == 9:
        seqid, source, type_, start, end, score, strand, phase, attributes = buf
    else:
        return False

    if type_ not in feature_set:
        return False

    i_list = common.detect_index_update( seqid, int(start), int(end), seqs )    
    flag = False
    for i, j in i_list:        
        if( i == -1 ):
            continue            
        height = seqs[i][j].height * args.gene_thickness
        tail_length = ( height / 2 ) * ( size.xlim_max / size.ylim_max ) * ( size.figsize_inch[1] / size.figsize_inch[0] )
        aninstance = Gene( seqs[i][j], int(start), int(end)+1, strand, height, tail_length, color, args.gene_edge_color )
        aninstance.plot( ax )

        gene_name = func_get_gene_name_from_gff( attributes, args.gene_label_attr )
        func_plot_gene_name( seqs[i][j], ax, int(start), int(end), gene_name, args.gene_font_size, height, args.gene_font_rotation, size.PT2INCH4Y )
        flag = True
    return flag


def plot_genes_from_gff( seqs, ax, size, fn, args ):
    color = args.gene_color
    flag = False
    if not os.path.isfile( fn ):
        print( f"Error: '{fn}' is not found.", file=sys.stderr )
        return flag
    with open( fn, 'r') as infile:
        for line in infile:
            line = line.strip()
            if not line or line.startswith("#"):
                if line.startswith("##FASTA"):
                    break  # Stop parsing at ##FASTA
                continue
            buf = line.split('\t')
            flag = deal_gff( seqs, ax, size, args, buf, color ) or flag
    return flag


def plot_genes_from_gff_excel( seqs, ax, size, fn, args ):
    color = args.gene_color
    flag = False
    if not os.path.isfile( fn ):
        print( f"Error: '{fn}' is not found.", file=sys.stderr )
        return flag
    df = pd.read_excel( fn, header=None )
    for index, row in df.iterrows():
        buf = row.tolist()
        flag = deal_gff( seqs, ax, size, args, buf, color ) or flag
    return flag

        
def is_genbank_format( file_path, max_lines=100, threshold=1 ):
    genbank_keywords = ['LOCUS', 'DEFINITION', 'ACCESSION', 'VERSION', 'FEATURES', 'ORIGIN']
    found = 0
    try:
        with open(file_path, 'r') as f:
            for i, line in enumerate(f):
                if any(line.startswith(k) for k in genbank_keywords):
                    found += 1
                    if found >= threshold:
                        return True
                if i >= max_lines:
                    break
        return False
    except Exception:
        return False

    
def plot_genes_from_gb( seqs, ax, size, fn, args ):
    color = args.gene_color
    edge_color = args.gene_edge_color
    feature_set = args.feature
    flag = False
    if not os.path.isfile( fn ):
        print( f"Error: '{fn}' is not found.", file=sys.stderr )
        return flag
    if not is_genbank_format( fn ) :
        print( f"Error: '{fn}' is not genbank format.", file=sys.stderr )
        return flag
    records = SeqIO.parse(fn, "genbank")
    for record in records:
        seq_id = record.id
        for feature in record.features:
            if feature.type not in feature_set:
                continue
            start = int(feature.location.start) + 1
            end = int(feature.location.end)
            strand = "+" if feature.location.strand == 1 else "-" if feature.location.strand == -1 else "."
            gene_name = feature.qualifiers.get( args.gene_label_attr, [""])[0]

            i_list = common.detect_index_update( seq_id, int(start), int(end), seqs )
            for i, j in i_list:
                if( i == -1 ):
                    continue
                height = seqs[i][j].height * args.gene_thickness
                tail_length = ( height / 2 ) * ( size.xlim_max / size.ylim_max ) * ( size.figsize_inch[1] / size.figsize_inch[0] )
                aninstance = Gene( seqs[i][j], start, end+1, strand, height, tail_length, color, args.gene_edge_color )
                aninstance.plot( ax )

                func_plot_gene_name( seqs[i][j], ax, start, end+1, gene_name, args.gene_font_size, height, args.gene_font_rotation, size.PT2INCH4Y )            
                flag = True
    return flag

def output_genes_parameters( args, flag ):
    if not flag:
        return
    print( '' )
    print( '## Genes paramenters' )
    print( '--gene_thickness %.2f' % ( args.gene_thickness ))
    print( '--gene_font_size %.1f' % ( args.gene_font_size ))
    print( '--gene_font_rotation %d' % ( args.gene_font_rotation ))
    print( '--gene_color "%s"' % ( args.gene_color ))
    print( '--gene_edge_color "%s"' % ( args.gene_edge_color ))
    print( '--feature %s' % (  " ".join(map( str, args.feature ))))
    print( '--gene_label_attr %s' % ( args.gene_label_attr ))
