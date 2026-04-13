import os.path
import re
import sys
import csv
import math
from matplotlib.path import Path
from matplotlib.lines import Line2D
from matplotlib.colors import to_rgb
from matplotlib.colors import to_rgba
import matplotlib.pyplot as plt
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
            x = [ x1, x2, x3, x2, x1, x1 ]
            y = [ y1, y1, y2, y3, y3, y1 ]
        elif( self.strand == '-' ): # <= or <
            x2 = self.start + tail_length if self.width > tail_length else self.end
            x = [ x1, x2, x3, x3, x2, x1 ]
            y = [ y2, y1, y1, y3, y3, y2 ]
        else:
            x = [ x1, x3, x3, x1, x1 ]
            y = [ y1, y1, y3, y3, y1 ]
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
        elif to_rgb(self.color.color) == (1.0, 1.0, 1.0) or self.color.color == "None":
            edge_color, lw = 'black', 0.15
        else:
            edge_color, lw = None, 0
        ax.fill( self.x, self.y, color=self.color.color, alpha=self.color.alpha, zorder=4, lw=0 )
        if edge_color is not None :        
            ax.plot( self.x, self.y, zorder=6, lw=lw, color=edge_color )

        
class Feature_color_legend:
    def __init__( self, args ):
        self.font_size = args.feature_color_legend_font_size
        self.padding = 2
        self.marker_size = args.feature_color_legend_marker_size
        self.legend_dict, self.pattern_dict, self.pattern2legend = self._read_file( args.feature_color_map )
        self._flag = {}
        self.ncol = args.feature_color_legend_ncol
        self.nrow = 0
        self.loc = 'upper center'
        self.anchor = ( 0.5, 0.97 )
        s_len = 0.5
        self.gene_marker = Path( [[-1*s_len, -1*s_len], [0, -1*s_len], [s_len, 0], [0, s_len], [-1*s_len, s_len], [-1*s_len, -1*s_len] ] )
        self.edge_color = args.gene_edge_color
        
    def _read_file( self, fn ):
        legend_dict = {}
        pattern_dict = {}
        pattern2legend = {}
        if fn is None:
            return legend_dict, pattern_dict, pattern2legend
        try:
            if not os.path.exists( fn ):
                raise FileNotFoundError(f"Error: {fn} does not exist.")
            if os.path.getsize(fn) == 0 :
                raise ValueError(f"Error: {fn} is empty.")
            df = pd.read_csv( fn, header=None, sep='\t' )

            if df.shape[1] < 3:
                raise ValueError( f"Error: {fn} must have three columns (label, pattern, color)." )
            if df.iloc[:, :3].isnull().any().any():
                raise ValueError( f"Error: {fn} contains missing values in required columns." )

            df = df.iloc[:, :3]
            df.columns = ["label", "pattern", "color"]

            df_valid = df[df["color"].apply( self._is_valid_color )].copy()
            legend_dict = dict(zip(df_valid.iloc[:, 0], df_valid.iloc[:, 2]))
            pattern_dict = dict(zip(df_valid.iloc[:, 1], df_valid.iloc[:, 2]))
            pattern2legend = dict(zip(df_valid.iloc[:, 1], df_valid.iloc[:, 0]))
            return( legend_dict, pattern_dict, pattern2legend )
        
        except (FileNotFoundError, ValueError) as e:
            print( e, file=sys.stderr )
            sys.exit( 1 )

    def set_gene_color( self, name, desc, color ):
        if not self.pattern_dict:
            return color
        for text, set_color in self.pattern_dict.items():
            buf = text.split("/")
            for s_text in buf :
                if re.compile( s_text ).search( name ) or re.compile( s_text ).search( desc ):
                    color = set_color
                    self._flag[ text ] = 1
                    break
        return color

    def output( self, fig ):
        if not self.legend_dict:
            return
        legend_handles = []
        flag = {}
        for l, c in self.pattern_dict.items():
            if l not in self._flag:
                continue
            if self.pattern2legend[l] in flag:
                continue
            if self.edge_color is None :
                handle = Line2D( [], [], marker=self.gene_marker, linestyle="",
                                 markersize=self.marker_size, markerfacecolor=c, markeredgecolor="none", markeredgewidth=0.0, label=self.pattern2legend[ l ] )
            else:
                handle = Line2D( [], [], marker=self.gene_marker, linestyle="",
                                 markersize=self.marker_size, markerfacecolor=c, markeredgecolor=self.edge_color, markeredgewidth=0.15, label=self.pattern2legend[ l ] )
            legend_handles.append(handle)
            flag[ self.pattern2legend[ l ] ] = 1

        if not legend_handles:
            return

        self.ncol, self.nrow = self._estimate_ncol( fig )
        fig.legend(
            frameon=False,
            handles=legend_handles,
            # bbox_transform=fig.transFigure,
            loc=self.loc,
            ncol=self.ncol,
            edgecolor="white",
            fontsize=self.font_size,
            handlelength=0.8,
            handletextpad=0.4,
            columnspacing=1.4,
            bbox_to_anchor=self.anchor )
        return

    def _estimate_ncol(self, fig ):
        if not self.legend_dict:
            return 0, 0
        if not self._flag:
            return 0, 0
        fig.canvas.draw()
        renderer = fig.canvas.get_renderer()
        text_widths = []
        for label in self._flag:
            t = fig.text(0, 0, self.pattern2legend[label], fontsize=self.font_size )
            bbox = t.get_window_extent( renderer=renderer )
            t.remove()
            text_widths.append(bbox.width)

        fig_width_px = fig.get_size_inches()[0] * fig.dpi
        avg_width = max(text_widths) + self.padding
        if self.ncol != 0 :
            ncol = self.ncol
        else:
            ncol = max(1, int(fig_width_px // avg_width))
            ncol = min(ncol, len(self._flag))            
        nrow = math.ceil( len( self._flag ) / ncol )
        return ncol, nrow

    def _is_valid_color(self, c ):
        try:
            to_rgba(c)
            return True
        except ValueError:
            print(f"'{c}' is not a valid color. Use a named color or hex code like '#RRGGBB'.", file=sys.stderr)
            return False

class RNA:
    def __init__( self ):
        self.start = -1
        self.end = -1
        self.strand = "."
        
    def set_RNA_position( self, type_, start, end, strand ):
        if re.compile( "gene" ).search( type_ ):
            self.start = -1
            self.end = -1
            self.strand = "."
        #if re.compile( "tRNA" ).search( type_ ):
            #return
        if re.compile( "RNA" ).search( type_ ):
            self.start = start
            self.end = end
            self.strand = strand
        
    def set_exon_strand( self, start, end, strand ):
        #print( self.start, self.end, start, end )
        if self.start == -1 :
            return strand
        elif self.strand == "+" and self.end == end :            
            return strand
        elif self.strand == "-" and self.start == start:            
            return strand
        return "."


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


def deal_gff( seqs, ax, size, args, buf, color, gene_legend, mRNA ):
    feature_set = args.feature
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
        gene_name = func_get_gene_name_from_gff( attributes, args.gene_label_attr )
        if len( buf ) == 9:
            color = gene_legend.set_gene_color( gene_name, attributes, color )
        tail_length = ( height / 2 ) * ( size.xlim_max / size.ylim_max ) * ( size.figsize_inch[1] / size.figsize_inch[0] )
        if type_ == "exon" or type_ == "CDS" :
            strand = mRNA.set_exon_strand( int(start), int(end), strand )
        aninstance = Gene( seqs[i][j], int(start), int(end)+1, strand, height, tail_length, color, args.gene_edge_color )
        aninstance.plot( ax )

        func_plot_gene_name( seqs[i][j], ax, int(start), int(end), gene_name, args.gene_font_size, height, args.gene_font_rotation, size.PT2INCH4Y )
        flag = True
    return flag


def plot_genes_from_gff( seqs, ax, size, fn, args, gene_legend ):
    color = args.gene_color
    flag = False
    if not os.path.isfile( fn ):
        print( f"Error: '{fn}' is not found.", file=sys.stderr )
        return flag
    mRNA = RNA( )
    with open( fn, 'r') as infile:
        for line in infile:
            line = line.strip()
            if not line or line.startswith("#"):
                if line.startswith("##FASTA"):
                    break  # Stop parsing at ##FASTA
                continue
            buf = line.split('\t')
            mRNA.set_RNA_position( buf[2], int(buf[3]), int(buf[4]), buf[6] )
            flag = deal_gff( seqs, ax, size, args, buf, color, gene_legend, mRNA ) or flag
    return flag


def plot_genes_from_gff_excel( seqs, ax, size, fn, args, gene_legend ):
    color = args.gene_color
    flag = False
    if not os.path.isfile( fn ):
        print( f"Error: '{fn}' is not found.", file=sys.stderr )
        return flag
    mRNA = RNA()
    df = pd.read_excel( fn, header=None )
    for index, row in df.iterrows():
        buf = row.tolist()
        if buf[0].startswith("#"):
            continue
        mRNA.set_RNA_position( buf[2], int(buf[3]), int(buf[4]), buf[6] )
        flag = deal_gff( seqs, ax, size, args, buf, color, gene_legend, mRNA ) or flag
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

    
def plot_genes_from_gb( seqs, ax, size, fn, args, gene_legend ):
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
    mRNA = RNA( )
    for record in records:
        seq_id = record.id
        for feature in record.features:
            start = int(feature.location.start) + 1
            end = int(feature.location.end)
            strand = "+" if feature.location.strand == 1 else "-" if feature.location.strand == -1 else "."
            mRNA.set_RNA_position( feature.type, start, end, strand )
            if feature.type not in feature_set:
                continue
            i_list = common.detect_index_update( seq_id, int(start), int(end), seqs )
            for i, j in i_list:
                if( i == -1 ):
                    continue
                gene_name = feature.qualifiers.get( args.gene_label_attr, [""])[0]
                product = feature.qualifiers.get( "product", [""])[0]
                color = gene_legend.set_gene_color( gene_name, product, args.gene_color )
                height = seqs[i][j].height * args.gene_thickness
                tail_length = ( height / 2 ) * ( size.xlim_max / size.ylim_max ) * ( size.figsize_inch[1] / size.figsize_inch[0] )
                mRNA.set_RNA_position( feature.type, start, end, strand )
                if feature.type == "exon" or feature.type == "CDS" :
                    strand = mRNA.set_exon_strand( start, end, strand )
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
