# a-liner
This repository contains the `a-liner` script and sample data.

## Contents
- [Requirements](#requirements)
- [Installation](#installation)
- [Sample Data](#sample-data)
- [Options](#options)
- [Prepare input files](#prepare-input-files)
- [Citation](#citation)

## Requirements
This project has been tested with the following environment:
- `python 3.13.5`
- `matplotlib 3.10.3`
- `pandas 2.3.1`
- `numpy 2.3.2`
- `biopython 1.85`
- `bcbio-gff 0.7.1`
- `openpyxl 3.1.5`


## Installation
### Install via Bioconda (Recommended)
You can install a-liner directly from Bioconda:

```
conda install -c bioconda a-liner
```
> Note: If you haven't configured Bioconda before, follow the instructions [here](https://bioconda.github.io/) to set up the necessary channels.

### Optional: Create a dedicated conda environment
```
conda create -n a-liner-env a-liner
conda activate a-liner-env
conda install -c bioconda a-liner
```

### Install manually from source (alternative method)
If you prefer not to use Bioconda:
```
git clone https://github.com/mokuno3430/a-liner.git
cd a-liner
chmod u+x a-liner
```

Required dependencies:
```
conda install python=3.13.5 matplotlib=3.10.3 numpy pandas biopython bcbio-gff openpyxl 
```

## Sample Data

## Options
### General options
```
  -h, --help            show this help message and exit
  -i, --input [file ...]
                        File(s): sequence info for display.
                        Format: tab-delimited with columns [seq_ID, start(1-based), end(1-based), strand(+ or -), display name]
                        Sequences are arranged from bottom to top in the figure.
  --xlsx Excel          File: seq_info.xlsx
  --xlsx_sheet str      sheet name of seq_info.xlsx
  --out str             Optional: prefix of PDF file (default: out).
  --figure_size width height
                        Optional: figure size as [width height](inch).
                        If width=0 → set to 6. If height=0 → auto(default: [6, 0]).
  -v, --version         show program's version number and exit
```

### Sequence layout options
```
  --seq_layout {left,center,right}
                        Optional: sequence layout (default: left).
  --margin_bw_seqs float
                        Optional: vertical margin between adjacent sequences. Default -1 means auto-adjust.
  --xlim_max int        Optional: maximum x-axis coordinate for plotting (bp). Default -1 means auto-adjust.
  --left_margin float   Optional: left side margin of the figure.
                        Default -1 means auto-adjust. range 0.05-0.50.
```

### Sequence drawing options
```
  --seq_color str       Optional: color of sequences (default grey).
  --seq_font_size float Optional: font size of sequence names (pt). Default 6.
  --seq_thickness float Optional: thickness of sequence lines (pt) (default: 1.5).
```

### Sequence scale options
```
  --scale {legend,tick,both}
                        Optional: how to display scale.
                        "legend" = show scale bar, "tick" = show axis ticks on sequences,
                        "both" = show both (default: legend).
  --tick_width int      Optional: scale width of axis (bp) (default -1 means auto).
  --tick_font_size float
                        Optional: font size of ticks (pt) (default: 3).
```

### Sequence alignment files
```
  -a, --alignment [file ...]
                        File(s): custom alignment data.
                        Format: tab-delimited with columns [seq_ID1, start1, end1, seq_ID2, start2, end2, identity(%)].
  --blastn [file ...]   File(s): blastn output.
                        Example: blastn -db ref.fa -query query.fa -out blastn.txt -outfmt 6
  --lastz [file ...]    File(s): lastz output.
                        Example: lastz ref.fa query.fa --format=general --output=lastz.txt
  --mummer [file ...]   File(s): MUMmer show-coords output.
                        Example: show-coords -H out.delta > show-coords.tsv
  --minimap2 [file ...]
                        File(s): minimap2 PAF output.
                        Example: minimap2 -c ref.fa query.fa > out.paf
```

### Sequence alignment options
```
  --min_identity int    Optional: minimum sequence identity (%).
                        Alignments below this threshold will be ignored (default: 70).
  --min_alignment_len int
                        Optional: minimum alignment length (bp).
                        Alignments shorter than this will be ignored (default: 0).
  --alignment_alpha float
                        Optional: transparency (alpha) of alignment coloring, range 0–1
                        (0 = fully transparent, 1 = opaque) (default: 0.5).
  --colormap {0,1,2,3,4,5}
                        Optional: colormap for sequence identity.
                        0 = bone_r, 1 = hot_r, 2 = BuPu, 3 = YlOrRd, 4 = YlGnBu, 5 = rainbow (original) (default: 5).
  --include_nonadjacent Include alignments between non-adjacent sequences (default: only adjacent).
```

### Gene annotation files
```
  --gff3 [gff3 ...]     File(s): gene annotation in GFF format.
  --gff_xlsx [Excel ...]
                        File(s): GFF format in Excel files
  --gb [genbank ...]    File(s): genbank format.
```

### Gene / feature drawing options
```
  --feature [str ...]   Optional: GFF/GenBank feature types to draw (space-separated)(default: gene).
  --gene_thickness float
                        Optional: relative thickness of gene arrows compared to seq_thickness (default: 3).
  --gene_label_attr str
                        Optional: attribute key used for feature labels (default: Name).
  --gene_font_size float
                        Optional: font size of gene names (pt) (default: 3).
  --gene_font_rotation float
                        Optional: rotation angle of gene names (degrees) (default: 75).
  --gene_color str      Optional: fill color of gene arrows (default: black).
  --gene_edge_color str
                        Optional: edge (outline) color of gene arrows (default: None) (no outline).
```

### Highlight options
```
  --highlight [file ...]
                        File(s): highlight regions.
                        Format: tab-delimited with columns [seq_ID, start(1-based), end(1-based), color]
  --h_alpha float       Optional: transparency of highlights (0=transparent, 1=opaque) (default: 0.3).
  --h_thickness float   Optional: relative thickness of highlights compared to sequence thickness (default: 3.5).
```

### Scatter plot options
```
  --scatter [file ...]  File(s): scatterplot data.
                        Format: tab-delimited with columns [seq_ID, position(1-based), value]
  --marker_color str    Optional: marker color (default: deeppink).
  --marker_size float   Optional: marker size (default: 3).
  --marker_style str    Optional: marker style.
                        Valid choices: *, ,, ., 8, <, >, D, H, P, X, ^, d, h, o, p, s, v (default: .).
  --scatter_space float
                        Optional: relative height of scatterplot compared to alignment space (default: 0.8).
                        For example, 0.8 means 80% of the alignment height.
  --scatter_min float   Optional: minimum value of y-axis (default: 0).
  --scatter_max float   Optional: maximum value of y-axis (default: 4).
  --scatter_ylines [float ...]
                        Optional: add horizontal reference lines at the given y values (list of floats).
  --background_color str
                        Optional: background color of scatter plot (default: whitesmoke).
  --sp_highlight [file ...]
                        File(s): highlight regions for scatter plot.
                        Format: tab-delimited with columns [seq_ID, start(1-based), end(1-based), color]
  --sp_h_alpha float    Optional: transparency of highlights for scatter plot(default: 0.3).
```

## Prepare input files
### Sequence Configuration Files (required)
a-liner supports two ways to specify how sequences are arranged in the figure:  
**TSV-based configuration** and **Excel-based configuration**.

**Option 1: TSV-based configuration** (`--input`)  
Prepare one TSV file **per sequence track (row)**.  
- TSV files are drawn **from bottom to top** in the order they are given to `--input`.
- Within each TSV file, sequences are drawn **from left to right** in the order they appear in the file.

Each TSV file must contain the following five columns:
- **Column 1:** Sequence ID
- **Column 2:** Start position (1-based)
- **Column 3:** End position (1-based)
- **Column 4:** Strand (`+` or `-`)
- **Column 5:** Display name

This format is suitable when the number of sequence tracks is small or when each track is prepared independently.  

**Option 2: Excel-based configuration** (`--xlsx`)

Alternatively, all sequence configuration information can be provided in a single Excel file using the `--xlsx` option.
- All sequence tracks are described in one sheet.
- Sequences assigned to the same track are drawn **from left to right** in the order they appear in the sheet.
- Track indices are specified explicitly, allowing flexible control of layout.

The Excel file must contain the following columns:
- **Column A:** Track index (0-based; `0` corresponds to the bottom track)
- **Column B:** Sequence ID
- **Column C:** Start position (1-based)
- **Column D:** End position (1-based)
- **Column E:** Strand (`+` or `-`)
- **Column F:** Display name

This format is recommended when many tracks are used or when layout needs to be adjusted interactively.


### Alignment Input Files (optional)
a-liner supports alignment results generated by major alignment tools, including **BLASTN**, **minimap2**, **LASTZ**, and **MUMmer**.  
Below are minimal example commands to produce alignment outputs in formats compatible with a-liner.

#### BLASTN
Use the `-outfmt 6` option to generate tabular BLASTN output.
```
makeblastdb -in seq1.fa -dbtype nucl
blastn -query seq2.fa -db seq1.fa -outfmt 6 -out output_blastn.txt
```

#### minimap2
Use the `-c` option to generate PAF output that includes CIGAR strings in the `cg` tag.
```
minimap2 -c seq1.fa seq2.fa > output_minimap2.paf
```

#### LASTZ
Use `--format=general` to produce a general tabular format.

```
lastz seq1.fa[multiple] seq2.fa --format=general --output=output_lastz.txt
```

#### MUMmer
Run `show-coords` with `-H` to generate a headerless coordinate table.
```
nucmer --prefix output_nucmer seq1.fa seq2.fa
show-coords -H output_nucmer.delta > output_nucmer.mcoords
```
  

### Annotation Files (optional)
a-liner supports several annotation formats.
GFF3 is the recommended format, but genbank format is also available.

Supported formats:
- **GFF3 files:** Load standard GFF3 annotations using the `--gff3` option.
- **GFF3 with embedded FASTA:** GFF3 files that include an embedded FASTA section are also supported and can be loaded using the `--gff3` option.
- **Excel files converted from GFF:** Excel files created from a GFF file.  
  (e.g., by manually importing or editing the annotation table) can be loaded with the `--gff_xlsx` option.
- **GenBank flat files:** Load GenBank annotations using the `--gb` option.  
  (Color information in GenBank files is NOT currently supported.)

**Coloring features**  
If your annotation file includes an additional column appended to the right of the standard 9-column GFF format, and this column contains a valid color name (e.g., `red`, `#FF0000`), a-liner will use this color when drawing the corresponding gene feature.

This applies to both:  
- plain GFF files with an extra color column
- Excel files derived from such GFF files
  

### Highlight Regions Files (optional)
To highlight specific regions, prepare a tab-delimited file and use either the `--highlight` or `--sp_highlight` option.  
The format is the same for both:

- **Column 1:** Sequence ID  
- **Column 2:** Start position (1-based)
- **Column 3:** End position (1-based)
- **Column 4:** Color specification (hex code, e.g., `#FF0000`, or a Matplotlib color name, e.g., `red`)

No header or index is required.  

**Example (highlights.txt):**

```
plasmid1    5000    15000    #FF9999
```
  

### Scatter Plot Data Files (optional)
To add scatter plots, provide a tab-delimited file using the `--scatter` option.  
Each row represents one data point.
- **Column 1:** Sequence ID
- **Column 2:** Position (1-based)
- **Column 3:** Value (e.g., SNP density, read depth, GC content)

Additional notes:
- Multiple scatter plot files can be specified at once by listing them after the `--scatter` option.
- Scatter plots are drawn above the corresponding sequence track.
- The y-axis scale is shared across all sequences. The minimum and maximum values can be controlled using `--scatter_min` and `--scatter_max`.
- Horizontal reference lines can be added at specified values using the `--scatter_ylines` option.
- Rows do not need to be sorted.
  

## Citation
Please cite the tool as follows until the formal publication becomes available:

> Okuno M. A-liner: linear alignment visualizer for genome comparisons.  
> Available at: https://github.com/mokuno3430/a-liner
