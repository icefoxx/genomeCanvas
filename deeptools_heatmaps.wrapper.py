import os, re
import io
from sys import argv, stdout
from pathlib import Path
import argparse
import logging
import inspect
import subprocess
import time
import datetime
import gzip
from itertools import islice, combinations
from collections import OrderedDict
import json ## https://appdividend.com/2022/01/29/how-to-convert-python-string-to-dictionary/

import pandas as pd
import numpy as np
from sklearn.decomposition import PCA

import seaborn as sns
import matplotlib.pyplot as plt
import plotly.graph_objects as go
import plotly.express as px
import matplotlib
import matplotlib.cm
from scipy.stats import gaussian_kde, pearsonr, spearmanr, mannwhitneyu
matplotlib.use('Agg')

from pybedtools import BedTool

from pdb import set_trace

import random

from sklearn.cluster import KMeans, MiniBatchKMeans
from sklearn.manifold import TSNE
from sklearn.preprocessing import StandardScaler
from sklearn.neighbors import NearestNeighbors

# Optional deps
try:
    import fastcluster
    import scipy.cluster.hierarchy as sch
except Exception:
    fastcluster = None
    sch = None

try:
    import igraph as ig
    import leidenalg as la
except Exception:
    ig = None
    la = None

# UMAP (for clustering)
try:
    import umap  # module presence check
    from umap import UMAP
except Exception:
    umap = None
    UMAP = None

class CustomFormatter(logging.Formatter):

    grey = "\x1b[38;20m" # test how it looks in bash with e.g.: echo -e "\x1b[34;20m text\e[0m"
    yellow = "\x1b[33;20m"
    blue = "\x1b[34;20m"
    pink = "\x1b[35;20m"
    red = "\x1b[31;20m"
    bold_red = "\x1b[31;1m"
    reset = "\x1b[0m"
    format = "###\t[%(asctime)s] %(filename)s:%(lineno)d: %(name)s %(levelname)s: %(message)s"

    FORMATS = {
        logging.DEBUG: blue + format + reset,
        logging.INFO: grey + format + reset,
        logging.WARNING: yellow + format + reset,
        logging.ERROR: red + format + reset,
        logging.CRITICAL: bold_red + format + reset
    }

    def format(self, record):
        log_fmt = self.FORMATS.get(record.levelno)
        formatter = logging.Formatter(log_fmt)
        return formatter.format(record)

def configureLogging(analysisPrefix=(os.path.basename(__file__)).replace(".py","")):
    logger = logging.getLogger()
    logger.disabled = True
    logger.setLevel(logging.INFO)
    
    streamhdlr = logging.StreamHandler(stdout)
    filehdlr  = logging.FileHandler(f"{analysisPrefix}.log")
    
    logger.addHandler(streamhdlr)
    logger.addHandler(filehdlr)
    
    streamhdlr.setLevel(logging.INFO)
    filehdlr.setLevel(logging.INFO)
    
    lgrPlainFormat = logging.Formatter('###\t[%(asctime)s] %(filename)s:%(lineno)d: %(name)s %(levelname)s: %(message)s')
    filehdlr.setFormatter(lgrPlainFormat)
    streamhdlr.setFormatter(CustomFormatter())

def str2bool(v):
    if str(v).lower() in ('yes', 'true', 't', 'y', '1'):
        return True
    elif str(v).lower() in ('no', 'false', 'f', 'n', '0'):
        return False
    else:
        logging = logging.getLogger(inspect.currentframe().f_code.co_name)
        logging.error("Unrecognized parameter was set for '{}'. Program was aborted.".format(v))
        exit()

def peak_center_bin_flank_type(v):
    """Argparse type: number of bins on each side of peak center (minimum 1)."""
    try:
        n = int(v)
    except (TypeError, ValueError):
        raise argparse.ArgumentTypeError(f"peakCenterBinFlank must be an integer >= 1 (got {v!r})")
    if n < 1:
        raise argparse.ArgumentTypeError(f"peakCenterBinFlank must be >= 1 (got {n})")
    return n

def testFileGzippness(inputFileName):
    with gzip.open(inputFileName, 'r') as fh:
        try:
            fh.read(1)
            logging.info(f"Gzipness True: Input file {inputFileName} is a valid gzip file by OSError")
            return True
        except OSError:
            logging.info(f"Gzipness False: Input file {inputFileName} is NOT a valid gzip file by OSError")
            return False

def rgb2hex(r,g,b): #https://stackoverflow.com/questions/3380726/converting-a-rgb-color-tuple-to-a-six-digit-code
    return "#{:02x}{:02x}{:02x}".format(r,g,b)

def hex2rgb(hexcode):
    return tuple(map(ord,hexcode[1:].decode('hex')))

def plotCorrMatrix(dataDf, outfilePrefix, samples, conditionGrp, metaGrp, samples2condGrp, samples2metaGrp):
    for cType in ["pearson", "kendall", "spearman"]:
        corrDf = dataDf.corr(method=cType)
        ax = sns.clustermap( corrDf,
                        method="ward",
                        cmap='Spectral_r', 
                        annot=True, 
                        annot_kws={"size": 7}, 
                        linewidths=.5,
                        vmin=-1, 
                        vmax=1, 
                        figsize=(15,15),
                        cbar_kws={'label': f"{cType} correlation"})
        ax.ax_heatmap.set_xticklabels(ax.ax_heatmap.get_xticklabels(), rotation=45, horizontalalignment='right')
        plt.savefig(f"{outfilePrefix}.corr.{cType}.pdf", bbox_inches='tight', dpi=300)
        corrDf.to_csv(f"{outfilePrefix}.corr.{cType}.tsv", sep="\t")

def plotPCA(dfPlot, outfilePrefix, samples, conditionGrp, metaGrp, samples2condGrp, samples2metaGrp, colorPalette = "Spectral"):
    
    dfPlot.dropna(inplace=True)
    sample_names_hover = list(dfPlot)
    df = dfPlot.T

    ### Compute data for max components number and no more than 10:
    maxComponentsNum = np.min([10, len(list(dfPlot))])
    pca = PCA(n_components=maxComponentsNum)
    components = pca.fit_transform(df)
    componentsExpl = list(pca.explained_variance_ratio_)

    ### Extract stats for first 2 components and save projection related files:
    PCsLabel = [f"PC{i+1}" for i in range(0, maxComponentsNum)]

    # Additionally generate a simple skree plot:
    PC_values = np.arange(pca.n_components_) + 1
    screePlot = go.Figure()

    screePlot.add_trace(
        go.Scatter(
            x=PC_values,
            y=pca.explained_variance_ratio_,
            line=dict(color="#e87d72")
        ))

    screePlot.add_trace(
        go.Bar(
            x=PC_values,
            y=pca.explained_variance_ratio_,
            marker_color="#ededed",
        ))

    screePlot.update_layout(template='plotly_white',
                      xaxis_title='Principal Component',
                      yaxis_title="Proportion of Variance Explained",
                      title="Scree Plot",
                      showlegend=False)
    
    
    total_var = np.sum(componentsExpl[:2]) * 100
    labels = {lab: f"PC {i} ({pcap*100:.2f}%)" for i, lab, pcap in zip([1,2], ["x","y"], componentsExpl[:2])}
    plottitle=f'Principal Component Analysis<br>Total Explained Variance by first 2 PCs: {total_var:.2f}%'

    used = set()
    conditionGrp2Color = OrderedDict()
    conditionGrpUnique = [x for x in conditionGrp if x not in used and (used.add(x) or True)]

    cmap = matplotlib.cm.get_cmap(colorPalette) # color selections based on https://stackoverflow.com/questions/25408393/getting-individual-colors-from-a-color-map-in-matplotlib and https://stackoverflow.com/questions/52251582/matplotlib-has-no-attribute-cm-when-deploying-an-app
    plotVal2cmapVal = matplotlib.colors.Normalize(vmin=0, vmax=len(conditionGrpUnique)-1)
    for grpIdx, grp in enumerate(conditionGrpUnique):
        color = cmap(plotVal2cmapVal(grpIdx))
        conditionGrp2Color[grp] = ['rgb({},{},{})'.format(int(255*color[0]), int(255*color[1]), int(255*color[2])), rgb2hex(int(255*color[0]), int(255*color[1]), int(255*color[2]))]
    
    color_discrete_map = {}
    for sample in samples:
        color_discrete_map[sample] = conditionGrp2Color[samples2condGrp[sample]][1]
    
    pcaPlot = px.scatter(
            x=components[:,0], y=components[:,1],
            labels=labels,
            title=plottitle,
            template = "ggplot2",
            color=samples,
            color_discrete_map=color_discrete_map,
#             width=900, 
#             height=850,
            hover_name=sample_names_hover
        )
    pcaPlot.update_traces(marker=dict(size=12,
                                      line=dict(width=2,
                                      color='DarkSlateGrey')),
                          selector=dict(mode='markers'))
    # pcaPlot = px.scatter(
#             x=components[:,0], y=components[:,1],
#             labels=labels,
#             title=plottitle,
#             template = "ggplot2",
# #             width=900, 
# #             height=850,
#             hover_name=sample_names_hover
#         )

    pcaPlot.write_html(f"{outfilePrefix}.pca.html")
    screePlot.write_html(f"{outfilePrefix}.scree.html")

def plotRibbon(df, outfilePrefix, samples, conditionGrp, metaGrp, samples2condGrp, samples2metaGrp, colorPalette = "Spectral"):
    samples2condGrp = OrderedDict()
    samples2metaGrp = OrderedDict()
    for s, c, m in zip(samples, conditionGrp, metaGrp):
        samples2condGrp[s] = c
        samples2metaGrp[s] = m

    tracks = []
    x_placement = 0
    ribbonWidth = 0.5
    sample_size = len(list(df))#df.shape[1]-1
    genomicPosition = list(df.index)

    used = set()
    conditionGrp2Color = OrderedDict()
    conditionGrpUnique = [x for x in conditionGrp if x not in used and (used.add(x) or True)]

    cmap = matplotlib.cm.get_cmap(colorPalette) # color selections based on https://stackoverflow.com/questions/25408393/getting-individual-colors-from-a-color-map-in-matplotlib and https://stackoverflow.com/questions/52251582/matplotlib-has-no-attribute-cm-when-deploying-an-app
    plotVal2cmapVal = matplotlib.colors.Normalize(vmin=0, vmax=len(conditionGrpUnique)-1)
    for grpIdx, grp in enumerate(conditionGrpUnique):
        color = cmap(plotVal2cmapVal(grpIdx))
        conditionGrp2Color[grp] = ['rgb({},{},{})'.format(int(255*color[0]), int(255*color[1]), int(255*color[2])), rgb2hex(int(255*color[0]), int(255*color[1]), int(255*color[2]))]

    for idx, sample in enumerate(list(df)):
        x_placement = idx+1
        z_raw = list(df[sample])#spectra[:, i]
        x = [] # Relative placement of individual "BigWig" files
        y = [] # Genomic bins
        z = [] # Signal enrichment
        ci = int(255/sample_size*x_placement) # ci = "color index"
        for j in range(0, len(z_raw)):
            z.append([z_raw[j], z_raw[j]])
            y.append([genomicPosition[j], genomicPosition[j]])
            x.append([x_placement, x_placement+ribbonWidth])
        color = cmap(plotVal2cmapVal(x_placement))
        tracks.append(dict(
            z=z,
            x=x,
            y=y,
            name=sample,
            colorscale=[ [i, conditionGrp2Color[samples2condGrp[sample]][0]] for i in np.arange(0,1.1,0.1) ],
            showscale=False,
            type='surface',
        ))

    ribbonPlot = go.Figure(data=tracks)#[go.Surface(x=x, y=y, z=z)])

    camera = dict(eye=dict(x=-1, y=-2, z=.5)) # relative starting position, based on: https://plotly.com/python/3d-camera-controls/
    ribbonPlot.update_layout(scene = dict(
                          xaxis_title='Trace<br>(x-axis)',
                          yaxis_title='<br>Position relative<br>to peak center<br>(y-axis)',
                          zaxis_title='Average signal enrichment<br>(z-axis)'),
                         #width=700,
                          height=700,
                          margin=dict(r=20, b=10, l=10, t=10),
                          font=dict(
                            family="Arial",
                            size=12,
                            color="black"
                          ),
                      scene_camera=camera,
                      template="plotly_white"
                     )

    ribbonPlot.write_html(f"{outfilePrefix}.ribbonPlot.html")

BW_META_COLUMNS = ["bigWig", "bigWigLabel", "conditionGrp", "metaGrp", "colorMap", "profileColor"]


def _normalize_bw_meta_header_field(field):
    field = str(field).strip()
    if field == "#bigWig":
        return "bigWig"
    return field


def _bw_meta_header_match(fields):
    """Return (is_header, normalized_fields) if fields are a valid metadata header row."""
    fields = [str(f).strip() for f in fields]
    if len(fields) != len(BW_META_COLUMNS):
        return False, None
    normalized = [_normalize_bw_meta_header_field(f) for f in fields]
    if normalized == BW_META_COLUMNS:
        return True, normalized
    return False, None


def _log_bw_meta_column_help(lgr):
    lgr.warning("The input TSV metadata file columns should be:")
    lgr.warning("1. 'bigWig' - full path to bigwig file")
    lgr.warning("2. 'bigWigLabel' - displayed name for the big wig file")
    lgr.warning("3. 'conditionGrp' - currently only used for ribbon plot to assign the same color for the traces from the same condition.")
    lgr.warning("4. 'metaGrp' - condition (keep the same name is you dont have e.g. KO and WT for each). This is used to re-normalzie the enrichment within the group")
    lgr.warning("5. 'colorMap' - name of the color map or customized colors if useColorList option is used")
    lgr.warning("6. 'profileColor' - hex code for color for profile plot")


def _parse_bw_meta_file_lines(metaDataFile, lgr):
    """
    Read metadata TSV: skip blank lines and '#' comment lines.
    Recognize optional header (plain or '#bigWig' first column). Return data lines only.
    """
    data_lines = []
    has_header = False

    with open(metaDataFile, "r") as fh:
        for raw in fh:
            line = raw.rstrip("\n\r")
            if not line.strip():
                continue

            if line.lstrip().startswith("#"):
                is_hdr, _ = _bw_meta_header_match(line.split("\t"))
                if is_hdr:
                    if has_header:
                        lgr.warning("Duplicate metadata header line ignored.")
                    else:
                        has_header = True
                        lgr.info("Metadata file: recognized header line (#bigWig ...).")
                    continue
                lgr.debug("Skipping metadata comment line: %s", line[:120])
                continue

            is_hdr, _ = _bw_meta_header_match(line.split("\t"))
            if is_hdr and not has_header and not data_lines:
                has_header = True
                lgr.info("Metadata file: recognized header line (bigWig ...).")
                continue

            data_lines.append(line)

    return data_lines, has_header


def _pad_bw_meta_columns(df, defcolor="Spectral_r", defprofilecolor="auto"):
    if len(list(df)) == 2:
        df.insert(df.shape[1], 2, "def")
    if len(list(df)) == 3:
        df.insert(df.shape[1], 3, "def")
    if len(list(df)) == 4:
        df.insert(df.shape[1], 4, defcolor)
    if len(list(df)) == 5:
        df.insert(df.shape[1], 5, defprofilecolor)
    return df


def readbwMeta(metaDataFile, lgr, defcolor="Spectral_r", defprofilecolor="auto"):
    if not os.path.exists(metaDataFile):
        msg = "Metadata file does not exist... {}".format(metaDataFile)
        lgr.critical(msg)
        raise FileNotFoundError(msg)

    data_lines, has_header = _parse_bw_meta_file_lines(metaDataFile, lgr)
    if not data_lines:
        msg = "Metadata file has no data rows (only header/comments?). Program was aborted."
        lgr.critical(msg)
        raise ValueError(msg)

    tsv_body = "\n".join(data_lines)
    if has_header:
        header_line = "\t".join(BW_META_COLUMNS)
        df = pd.read_csv(
            io.StringIO(header_line + "\n" + tsv_body),
            sep="\t",
            header=0,
            index_col=False,
        )
        df.columns = [_normalize_bw_meta_header_field(c) for c in df.columns]
        if list(df.columns) != BW_META_COLUMNS:
            msg = (
                "Metadata header columns are invalid after parsing. "
                f"Expected {BW_META_COLUMNS}, got {list(df.columns)}."
            )
            lgr.critical(msg)
            _log_bw_meta_column_help(lgr)
            raise ValueError(msg)
    else:
        df = pd.read_csv(io.StringIO(tsv_body), sep="\t", header=None, index_col=False)
        df = _pad_bw_meta_columns(df, defcolor=defcolor, defprofilecolor=defprofilecolor)
        if len(list(df)) != 6:
            msg = "Metadata file does not have all the required columns. Please correct that. Program was aborted."
            lgr.critical(msg)
            _log_bw_meta_column_help(lgr)
            raise ValueError(msg)
        df.columns = BW_META_COLUMNS

    return df.copy()


def parseArgs():
    parser = argparse.ArgumentParser()
    #group = parser.add_mutually_exclusive_group(required=True)
    #group.add_argument("-p", "--prepareMode", help="[REQUIRED with -p or -c] -p will invoke blank config file generation, -c will run the analysis.", action="store_true", default=False, required=True, dest="prepareMode")
    #group.add_argument("-c", "--configFileName", help="[REQUIRED with -p or -c] Name of the config file.", default=False, action="store", type=str, required=False, dest="configFileName")
    #parser.add_argument("-q", "--queueType", help="(optional) Queue type. 'local' will run locally without submitting to HPC.", action="store", type=str, required=False, dest="queueType")
    
    requiredParams = parser.add_argument_group('REQUIRED parameters')
    requiredParams.add_argument("-r", "--regions", help="comma separated or *.lst file with a list of input BED files. If the regions file has '*.vout.anno.Ranks.tsv' ending, it means its part of the differential peak calling pipeline, thus the regions will be extracted from that annotated file. Note that this parameter is REQUIRED unless the script works in BIN mode", default="Wait, what? This should never be here, did someone forgot to specify regions in the regions mode?", action="store", type=str, required=False, dest="regions")
    requiredParams.add_argument("-p", "--analysisPrefix", help="Analysis prefix part of which can contain output directory, but this would have to be full path to work properly.", action="store", type=str, required=True, dest="analysisPrefix")
    requiredParams.add_argument("-m", "--metaDataFile", help="Tab-separated metadata file. Optional header row: bigWig, bigWigLabel, conditionGrp, metaGrp, colorMap, profileColor (or same with '#bigWig' as first column). Lines starting with '#' are comments unless that header form. Without a header, columns are assigned in that order (2-6 columns allowed; missing fields get defaults).", action="store", type=str, required=True, dest="metaDataFile")

    optionalParams = parser.add_argument_group('OPTIONAL parameters')
    optionalParams.add_argument("-l", "--regionsLabels", help="comma separated or *.lst file with a list of input BED files LABELS. Leave this as 'auto' if you wish the automated selection.", default='auto', action="store", type=str, required=False, dest="regionsLabels")
    optionalParams.add_argument("-f", "--flank", help="Flanking region size in bp. By default='2000", default=2000, action="store", type=int, required=False, dest="flank")
    optionalParams.add_argument("-@", "--n_proc", help="Number of processors/threads. By default = '16'.", default=16, action="store", type=int, required=False, dest="n_proc")
    optionalParams.add_argument("--refPointLabel", help="refPointLabel. By default = ' '.", default=" ", action="store", type=str, required=False, dest="refPointLabel")
    optionalParams.add_argument("--xAxisLabel", help="xAxisLabel. By default = 'peak center'.", default="peak center", action="store", type=str, required=False, dest="xAxisLabel")
    optionalParams.add_argument("--heatmapHeightBase", help="heatmapHeightBase. By default = '10'.", default=10, action="store", type=int, required=False, dest="heatmapHeightBase")
    optionalParams.add_argument("--sortRegions", help="sortRegions. By default = 'descend'. ", default="descend", action="store", type=str, required=False, dest="sortRegions", choices=["descend", "ascend", "no", "keep"])
    optionalParams.add_argument("--sortUsing", help="sortUsing (method to average the BIN-based enrichment data). By default = 'mean'. ", default="mean", action="store", type=str, required=False, dest="sortUsing", choices=["mean", "median", "max", "min", "sum", "region_length"])
    optionalParams.add_argument("--sortUsingSamples", help="sortUsingSamples; Comma separated list of samples, given as relateive positions, counted starting from 1, from left to right. Example: --sortUsingSamples 1,3. By default = ' ', which has i.e. no value and that use all samples (this is also default mode for deeptools).", default=" ", action="store", type=str, required=False, dest="sortUsingSamples")
    optionalParams.add_argument("--whatToShow", help="whatToShow. By default = 'plot, heatmap and colorbar'.", default="plot, heatmap and colorbar", action="store", type=str, required=False, dest="whatToShow")
    optionalParams.add_argument("--useColorList", help="useColorList; True or False, which indicates that the colors provided in input metadata file are specified as 'colorList' (i.e. custom colors) instead of 'colorMap' (which use also custom colors, but from predefined sets of heatmaps listed here: https://matplotlib.org/2.0.2/users/colormaps.html). Be defaults = False.", default=False, action="store", type=str2bool, required=False, dest="useColorList")
    optionalParams.add_argument("--extras", help="Plot extras? i.e. correlation plots, etc. By default = False, setting to False will essentially only plot stuff related with deeptoolsHeatmaps.", default=True, action="store", type=str2bool, required=False, dest="extras")
    optionalParams.add_argument("--overwrite", help="Should the matrix files be overwritten? By default = False.", default=False, action="store", type=str2bool, required=False, dest="overwrite")
    optionalParams.add_argument("--diffPeakThresholds", help="If the script is run as a part of differential peak calling pipeline, specify here which threshold of peaks should be extracted from supplied VOOM file. By default = 'Up2_Down2'.", default='Up2_Down2', action="store", type=str, required=False, dest="diffPeakThresholds", choices=["Up2FC_Down2FC", "Up2_Down2", "Up2NoFDR_Down2NoFDR", "Up_Down", "UpNoFDR_DownNoFDR", "Up2FC_Down2FC_Control", "Up2_Down2_Control", "Up2NoFDR_Down2NoFDR_Control", "Up_Down_Control", "UpNoFDR_DownNoFDR_Control", "Control"])
    optionalParams.add_argument("--maxRegionsPG", help="Max number of regions to be included into the PairGrid plot. This is to prevent overwhelming from too many regions, as the RAM consumption can grow out of hands at some point. If the numebr of regions is higher than specificed here, the data is ranomly subsampled. By default = 5000, which is reasonably high to keep the I think be quite confident of the output stats.", default=5000, action="store", type=int, required=False, dest="maxRegionsPG")


    optionalParams.add_argument(
        "--inputType",
        help="Input mode: 'region' (default) for user-supplied peaks or 'bin' for fixed-width genomic bins; heat-map steps are skipped in 'bin' mode.",
        default="region",
        type=str,
        choices=["region", "bin"],
        dest="inputType",
    )

    optionalParams.add_argument(
        "--clusterRegions",
        help="Cluster regions for heatmap. 'none' (default) disables clustering; 'kmeans' adds --kmeans, 'hclust' adds --hclust for using the deeptools internal clustering methods from the plotHeatmap function. External methods: hclustExt (hierarchal clustering using Ward), kmeansExt (k-means from scikit learn, preferred for data < 50k datapoints), miniBatchKMeansExt (MiniBatchKMeans from scikit learn, preferred for data >50k points), leidenExt (using Leiden graph-based clustering from leidenalg package). The usage of the clustering is valid only when --inputType is set to `region` and exactly one region BED file is provided at the input. Note, using hierarchal clustering is only possible for 'small' sample sizes, ≤2-5k regions feels OK; 5-10k is doable on with more RAM added; >10k gets painful; >20k is often infeasible.",
        default="none",
        type=str,
        choices=["none", "kmeans", "hclust", "hclustExt", "kmeansExt", "miniBatchKMeansExt", "leidenExt", "leidenKmeansExt"],
        dest="clusterRegions",
    )

    optionalParams.add_argument(
        "--gBinSize",
        help="Genomic bin size in bp. By default = 1000. This is only used when '--inputType bin' is specified.",
        default=1000,
        type=int,
        dest="gBinSize",
    )

    optionalParams.add_argument(
        "--clusterNumber",
        help="Number of clusters to use for clustering regions. By default = 5. This is only used when '--inputType region' and '--clusterRegions' are specified.",
        default=5,
        type=int,
        dest="clusterNumber",
    )

    optionalParams.add_argument(
        "--rndSeed",
        help="Random seed for reproducibility. Use an integer (default: '7201') or 'auto' to disable seeding.",
        default="7201",
        type=str,
        dest="rndSeed",
    )

    optionalParams.add_argument(
        "--reuseRenormedMatrix",
        help="If True (default), reuse an existing *.renormMatrix.mat.gz in clustering mode instead of recomputing it. Set to False to always overwrite.",
        default=True,
        action="store",  # keep as string to pass through str2bool
        type=str2bool,
        required=False,
        dest="reuseRenormedMatrix",
    )

    optionalParams.add_argument(
        "--clusteringDimReduct",
        help="Space on which to perform clustering: 'raw' matrix (default), PCA-reduced, or t-SNE embedding.",
        default="raw",
        type=str,
        choices=["raw", "pca", "tsne", "umap"],
        dest="clusteringDimReduct",
    )

    optionalParams.add_argument(
        "--leidenResolutionParameter",
        help=("Resolution parameter for Leiden (leidenalg). Practical ranges: "
            "~0.2-0.8 -> fewer/larger clusters; 1.0 (default); 1.5-5.0 -> more/smaller clusters. "
            "Exact effect depends on graph construction and data scale."),
        default=1.0,              # leidenalg default
        type=float,
        dest="leidenResolutionParameter",
    )

    optionalParams.add_argument(
        "--renormApproach",
        help="How to renormalize the matrix before plotting heatmaps/profiles. 'profileMinMax' (default) uses per-metaGrp min/max based on profile curves; 'qnorm' performs quantile normalization of matrix values per metaGrp across samples.",
        default="profileMinMax",
        type=str,
        choices=["profileMinMax", "qnorm"],
        dest="renormApproach",
    )

    optionalParams.add_argument(
        "--clusterUsingSamples",
        help="1-based comma-separated sample indices to use for clustering (e.g., '1,3,5'). Defaults to all samples byt setting .",
        default="",
        type=str,
        required=False,
        dest="clusterUsingSamples",
    )

    optionalParams.add_argument(
        "-s", "--start_point",
        help="referencePoint start point for deeptools heatmap plotting; 'center' (default), 'TSS', 'TES'",
        default="center",
        type=str,
        choices=["center", "TSS", "TES"],
        dest="start_point",
    )

    optionalParams.add_argument(
        '--bodyLength',
        help="Length of the body region when --start_point is set to TSS or TES. Default: 5000 (5kb).",
        default=5000,
        type=int,
        dest="bodyLength",
    )

    optionalParams.add_argument(
        "--peakCenterBinFlank",
        help="Number of bins on each side of the peak-center bin used for BED-group enrichment comparisons (window = center +/- N bins). Minimum = 1. Default = 5.",
        default=5,
        type=peak_center_bin_flank_type,
        required=False,
        dest="peakCenterBinFlank",
    )

    optionalParams.add_argument(
        "--peakCenterSignifThreshold",
        help="P-value threshold for peak-center BED comparison clustermaps: cells with p >= threshold are shown in grey (clustering still uses all values). Default = 0.05.",
        default=0.05,
        type=float,
        required=False,
        dest="peakCenterSignifThreshold",
    )

    args = parser.parse_args()
    
    lgr = logging.getLogger(inspect.currentframe().f_code.co_name)
    lgr.info("Current working directory: {}".format(os.getcwd()))
    lgr.info("Command used to run the program: python {}".format(' '.join(str(x) for x in argv)))
    #logging.info("Prepare mode (-p flag): {}".format(args.prepareMode))
    #logging.info("ConfigFileName (-c flag): {}".format(args.configFileName))
    errors = False
    
    if args.regions == "Wait, what? This should never be here, did someone forgot to specify regions in the regions mode?" and args.inputType == "region":
        msg = "Regions file (-r flag) is required in 'region' input mode. Please correct that. Program was aborted."
        lgr.critical(msg)
        errors = True
    
    if (args.regions).endswith(".lst") and args.inputType == "region":
        if os.path.exists(args.regions) == True:
            lgr.info("Loading in {} file...".format(args.regions))
            regions = list(pd.read_csv(args.regions, header=None)[0])
            diffPeakMode = False
            
            for region in regions:
                if os.path.exists(region) == False:
                    lgr.critical("Specified region file from *.lst file does not exist... {}".format(region))
                    errors = True
        else:
            lgr.critical("Regions file does not exist... {}".format(args.regions))
            errors = True
    elif (args.regions).endswith(".vout.anno.Ranks.tsv") and args.inputType == "region":
        if os.path.exists(args.regions) == True:
            lgr.info("Loading in {} file...".format(args.regions))
            regions = pd.read_csv(args.regions, sep="\t")
            diffPeakMode = True
        else:
            lgr.critical("Regions file does not exist... {}".format(args.regions))
    else:
        if args.inputType == "region":
            regions = (args.regions).split(",")
            for region in regions:
                if os.path.exists(region) == False:
                    lgr.critical("Specified region file from comma-separated list does not exist... {}".format(region))
                    errors = True
            diffPeakMode = False
        else:
            lgr.info("Input type is 'bin', thus no regions file is required, but the genomic bin size (--gBinSize flag) is used instead.")
            regions = []
            diffPeakMode = False
    lgr.info("Regions (-r flag): {}".format(args.regions))
    
    if args.inputType == "region":
        if (args.regionsLabels).endswith(".lst"):
            lgr.info("Loading in {} file...".format(args.regionsLabels))
            regionsLabels = list(pd.read_csv(args.regionsLabels, header=None)[0])
        else:
            regionsLabels = (args.regionsLabels).split(",")
        lgr.info("Regions' labels (-l flag): {}".format(regionsLabels))
        
        if (args.regionsLabels != "auto" and len(regionsLabels) != len(regions)) and not (args.regions).endswith(".vout.anno.Ranks.tsv"):
            msg = "List of regions of interest (i.e. BED files), is not equal to the lenght of regions labels. Please correct that. Program was aborted."
            lgr.critical(msg)
            errors = True
    else:
        lgr.info("Input type is 'bin', thus no regions labels are required, but the genomic bin size (--gBinSize flag) is used instead.")
        regionsLabels = []
    
    if args.clusterNumber < 1:
        msg = "Cluster number (--clusterNumber flag) must be at least 1. Please correct that. Program was aborted."
        lgr.critical(msg)
        errors = True
    
    if len(regions) > 1 and args.clusterRegions != "none":
        msg = "Clustering of regions (--clusterRegions flag) is only supported for single region file. Please correct that. Program was aborted."
        lgr.critical(msg)
        errors = True
    
    lgr.info("Analysis prefix (-p flag): {}".format(args.analysisPrefix))
    lgr.info("Metadata file (-m flag): {}".format(args.metaDataFile))
    tmp = readbwMeta(args.metaDataFile, lgr)

    lgr.info("Flank (-p flag): {}".format(args.flank))
    lgr.info("n_proc (-@ flag): {}".format(args.n_proc))
    lgr.info("refPointLabel (--refPointLabel flag): '{}'".format(args.refPointLabel))
    lgr.info("xAxisLabel (--xAxisLabel flag): {}".format(args.xAxisLabel))
    lgr.info("heatmapHeightBase (--heatmapHeightBase flag): {}".format(args.heatmapHeightBase))
    lgr.info("sortRegions (--sortRegions flag): {}".format(args.sortRegions))
    lgr.info("sortUsing (--sortUsing flag): {}".format(args.sortUsing))
    lgr.info("sortUsingSamples (--sortUsingSamples flag): {}".format(args.sortUsingSamples))
    lgr.info("whatToShow (--whatToShow flag): {}".format(args.whatToShow))
    lgr.info("use colorList instead of colorMap (--useColorList flag): {}".format(args.useColorList))
    lgr.info("Plot Extras? (--extras flag): {}".format(args.extras))
    lgr.info("overwrite matrix file? (--overwrite flag): {}".format(args.overwrite))
    lgr.info("Specify differential thresholds to plot (--diffPeakThresholds flag): {}".format(args.diffPeakThresholds))
    lgr.info("Max number of regions for PairGrid plot (--maxRegionsPG flag): {}".format(args.maxRegionsPG))
    lgr.info("Input type (--inputType flag): {}".format(args.inputType))
    lgr.info("Cluster regions (--clusterRegions flag): {}".format(args.clusterRegions))
    lgr.info("Genomic bin size (--gBinSize flag): {}".format(args.gBinSize))
    lgr.info("Cluster number (--clusterNumber flag): {}".format(args.clusterNumber))
    lgr.info("Random seed (--rndSeed flag): {}".format(args.rndSeed))
    lgr.info("Reuse renormalized matrix if available (--reuseRenormedMatrix flag): {}".format(args.reuseRenormedMatrix))
    lgr.info("Clustering dimension reduction (--clusteringDimReduct flag): {}".format(args.clusteringDimReduct))
    lgr.info("Leiden resolution parameter (--leidenResolutionParameter flag): {}".format(args.leidenResolutionParameter))
    lgr.info("Renormalization approach (--renormApproach flag): {}".format(args.renormApproach))
    lgr.info("Cluster using samples (--clusterUsingSamples flag): {}".format(args.clusterUsingSamples if args.clusterUsingSamples != "" else "all samples"))
    lgr.info("Peak-center bin flank (--peakCenterBinFlank flag): {}".format(args.peakCenterBinFlank))
    lgr.info("Peak-center significance threshold (--peakCenterSignifThreshold flag): {}".format(args.peakCenterSignifThreshold))


    if errors == True:
        msg = f"Critical errors encountered during argiments parsing (more details above). Program was aborted."
        lgr.critical(msg)
        raise Exception(msg)

    return regions, args.regions, regionsLabels, args.analysisPrefix, args.metaDataFile, args.flank, args.n_proc, args.refPointLabel, args.xAxisLabel, args.heatmapHeightBase, args.sortRegions, args.sortUsing, args.sortUsingSamples, args.whatToShow, args.useColorList, args.extras, args.overwrite, diffPeakMode, args.diffPeakThresholds, args.maxRegionsPG, args.inputType, args.clusterRegions, args.gBinSize, args.clusterNumber, args.rndSeed, args.reuseRenormedMatrix, args.clusteringDimReduct, args.leidenResolutionParameter, args.renormApproach, args.clusterUsingSamples, args.start_point, args.bodyLength, args.peakCenterBinFlank, args.peakCenterSignifThreshold

def corrfuncR(x, y, **kws):
    if len(x) > 2:
        r, _ = pearsonr(x, y)
        rho, _ = spearmanr(x, y)
        ax = plt.gca()
        ax.annotate("r = {:.3f}\nRho = {:.3f}".format(r, rho),
                    xy=(.9, .4), 
                    xycoords=ax.transAxes,
                    fontsize=17, ha="right")
    else:
        ax = plt.gca()
        ax.annotate("r = NA\nRho = NA",
                    xy=(.9, .4), 
                    xycoords=ax.transAxes,
                    fontsize=17, ha="right")

def densityScatter(xM, yM, **kws):
    x = np.array(xM)
    y = np.array(yM)

    if len(x) > 10:
        # Calculate the point density
        xy = np.vstack([x,y])
        try:
            z = gaussian_kde(xy)(xy)
            ax = sns.scatterplot(x=xM, y=yM, hue=z, s=3, edgecolor=None, palette="plasma", legend=None)
        except np.linalg.LinAlgError:
            ax = sns.scatterplot(x=xM, y=yM, s=3, edgecolor=None, color="#70777d", legend=None)
    else:
        ax = sns.scatterplot(x=xM, y=yM, s=3, edgecolor=None, color="#70777d", legend=None)

def _load_pairgrid_tsv_sampled(tsv_path: str, max_rows: int) -> pd.DataFrame:
    """
    Stream-read a multiBigwigSummary TSV and reservoir-sample up to max_rows rows
    without loading the entire file into memory.
    """
    k = max(1, int(max_rows))
    reservoir = []
    n_seen = 0
    col_names = None

    for chunk in pd.read_csv(tsv_path, sep="\t", chunksize=50000, low_memory=False):
        if col_names is None:
            chunk = chunk.copy()
            chunk.columns = [c.replace("'", "").replace("#", "") for c in chunk.columns]
            col_names = list(chunk.columns)
        for row in chunk.itertuples(index=False, name=None):
            n_seen += 1
            if len(reservoir) < k:
                reservoir.append(row)
            else:
                j = np.random.randint(0, n_seen)
                if j < k:
                    reservoir[j] = row

    if not reservoir:
        return pd.DataFrame(columns=col_names or [])
    return pd.DataFrame(reservoir, columns=col_names)


def plotPairGrid(npzTsvFile, pariGridFilePrefix, maxRegionsPG=5000, overwrite=False, lgr=None):
    out_png = f"{pariGridFilePrefix}.png"
    out_pdf = f"{pariGridFilePrefix}.pdf"
    if os.path.exists(out_png) and os.path.exists(out_pdf) and not overwrite:
        if lgr:
            lgr.info(f"PairGrid outputs exist, skipping: {out_png}")
        return

    ### Preprocess the input data GSE231362s.FOXA1__VS__FH.Up_Down.multiBigwigSummary.tsv
    df = _load_pairgrid_tsv_sampled(npzTsvFile, maxRegionsPG)
    if df.empty:
        if lgr:
            lgr.warning(f"PairGrid skipped: no rows in {npzTsvFile}")
        return

    df["Region"] = (
        df["chr"].astype(str) + ":" + df["start"].astype(str) + "-" + df["end"].astype(str)
    )
    df.set_index("Region", inplace=True)
    df.drop(columns=["chr", "start", "end"], inplace=True)
    df.columns = [autobreaklinetitle(col) for col in list(df)]
    sample_cols = list(df.columns)
    df[sample_cols] = np.log10(df[sample_cols].astype(np.float32) + 1.0)
    df.replace([np.inf, -np.inf], np.nan, inplace=True)

    # df.dropna(inplace=True)
    df = df.replace(np.nan, 0)

    ### double check if all columns have at least 5 (arbitrary number) unique values in them, if not, add a small random number to each value to make them unique. This is to protect from the scenario where one sample had only NaN/0 values or their equivalent, and in effect the result plot was not correclty generated and pipeline was failing. We however want to see the result because it will immediately look like annomally.
    for sample in sample_cols:
        if len(df[sample].unique()) < 5:
            df[sample] = df[sample].apply(lambda x: x + np.random.uniform(0, 0.0001))
    
    # Plot
    g = sns.PairGrid(df)
    g.map_diag(sns.histplot, color="#dfe2e6", stat='density', edgecolor="#65737e", linewidth=1, bins=15)
    if len(df) >= 10:
        g.map_diag(sns.kdeplot, color="#343d46")
    g.map_lower(densityScatter)
    g.map_upper(corrfuncR)
    g.fig.subplots_adjust(top=0.95)
    g.fig.suptitle(f'PairGrid plot for the log10(enrichment+1) of {len(df)} genomic regions')
    
    g.savefig(out_png, bbox_inches='tight', dpi=300)
    g.savefig(out_pdf, bbox_inches='tight', dpi=300)
    plt.close(g.fig)
    if lgr:
        lgr.info(f"PairGrid plot saved: {out_png}, {out_pdf}")

def plotProfile(dataDf, prefix, samplesList, samples2color, numberOfPeaks):
    """
    This function is used to generate an equivalent of `plotProfile` function from deeptools, but using our precalculated data from dataDf
    """

    lgr = logging.getLogger(inspect.currentframe().f_code.co_name)
    # print(dataDf)
    # print(dataDf.max())
    # open a new figure:
    fig = plt.figure(figsize=(10, 10))
    ax = fig.add_subplot(111)

    ### Generate scatter plot, with hidden data points and only displayed as lineplot, for each sample and add it to the figure, using idx column (index) as x axis, and sample column as y axis:
    maxes = []
    for sample in samplesList:
        #ax.scatter(dataDf.index, dataDf[sample], s=0.1, color=samples2color[sample])
        ax.plot(dataDf["idx"], dataDf[sample], label=sample, linewidth=2, color=samples2color[sample])
        maxes.append(dataDf[sample].max())

    ### set y axis limits:
    # print(np.max(maxes))
    # ax.set_ylim(0, np.max(maxes))

    ### Add title with the number of peaks from the heatmap (i.e. the dataDf length)
    plt.title(f"Number of peaks: {numberOfPeaks}")

    ### Add legend that will be located outside of the plot box:
    ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), fontsize=12)

    ### Save the figure:
    plt.savefig(f"{prefix}.profile.pdf", bbox_inches='tight', dpi=300)

    lgr.info(f"Plotting profile for {len(samplesList)} samples completed.")

def autobreaklinetitle(title,sep="[-_,\.:]",lmax=15):
    outsep = "-"
    sss = [ re.sub("['# ]","",rr) for rr in re.split(sep,title) if len(rr) ]
    newtitle, tmp = "", ""
    for ss in sss:
        tmp0 = tmp
        tmp += ss
        if len(tmp) > lmax:
            newtitle += tmp0.strip(outsep) + "\n"
            tmp = ss
        tmp += outsep
    newtitle += tmp.strip(outsep)
    newtitle = "\n" + newtitle
    return newtitle

def setRandomSeed(rndSeed: str, lgr=None):
    """
    Set reproducible random seeds across numpy and random.
    rndSeed: 'auto' or an integer-like string.
    """
    if rndSeed == "auto":
        if lgr: lgr.info("rndSeed='auto' -> seeding disabled (stochastic).")
        return None
    try:
        seed = int(rndSeed)
    except Exception:
        seed = 7201
        if lgr: lgr.warning("rndSeed not int/auto; falling back to 7201.")
    random.seed(seed)
    np.random.seed(seed)
    if lgr: lgr.info(f"Random seed set to {seed}.")
    return seed


def loadDeepToolsMatrix(matPath: str, lgr=None):
    """
    Load a deepTools computeMatrix .mat(.gz) into:
      header: dict parsed from JSON line
      coords: DataFrame with columns [chrm,start,end,name,score,strand]
      X: np.ndarray (n_regions x n_bins_total) float32
    """
    if not os.path.exists(matPath):
        raise FileNotFoundError(f"Matrix file not found: {matPath}")
    # Read first line to parse header
    opener = gzip.open if matPath.endswith(".gz") else open
    with opener(matPath, "rb") as fh:
        first = fh.readline().decode("utf-8").strip()
    if first.startswith("@"):
        header = json.loads(first[1:])
        skiprows = 1
    else:
        # rare older format; try to parse JSON later if it appears in line 2
        header = None
        skiprows = 0

    # If header is None, attempt to read again
    if header is None:
        with opener(matPath, "rb") as fh:
            lines = [fh.readline().decode("utf-8") for _ in range(2)]
        if lines and lines[0].startswith("@"):
            header = json.loads(lines[0][1:])
            skiprows = 1
        elif len(lines) > 1 and lines[1].startswith("@"):
            header = json.loads(lines[1][1:])
            skiprows = 1
        else:
            raise RuntimeError("Could not locate deepTools JSON header in matrix file.")

    totalBins = int(header["sample_boundaries"][-1])
    usecols = list(range(0, 6 + totalBins))
    names = ["chrm", "start", "end", "name", "score", "strand"] + list(range(totalBins))

    df = pd.read_csv(
        matPath, sep="\t", header=None, skiprows=skiprows,
        usecols=usecols, compression="infer", low_memory=False, dtype={0:str,1:int,2:int,3:str,4:str,5:str}
    )
    df.columns = names

    coords = df[["chrm","start","end","name","score","strand"]].copy()
    X = df.iloc[:, 6:].to_numpy(dtype=np.float32)

    # Handle NaNs
    if np.isnan(X).any():
        col_means = np.nanmean(X, axis=0)
        inds = np.where(np.isnan(X))
        X[inds] = np.take(col_means, inds[1])
    if lgr:
        lgr.info(f"Loaded matrix {matPath}: regions={X.shape[0]}, bins={X.shape[1]}")
    return header, coords, X


def sanitize_label_for_filename(label: str) -> str:
    """Return a filesystem-safe token derived from a sample or group label."""
    safe = re.sub(r"[^\w\-.]+", "_", str(label).strip())
    safe = re.sub(r"_+", "_", safe).strip("_.")
    return safe or "sample"


def peak_center_combo_label(bed_label: str, sample_label: str) -> str:
    """Build a BED x sample group label for peak-center comparison matrices."""
    return f"BED_{bed_label}._.sample_{sample_label}"


def compare_peak_center_between_bed_groups(matFile,
                                           outPrefix,
                                           samplesList,
                                           metaGrp,
                                           regionsLabels,
                                           lgr=None,
                                           center_bin_flank=5,
                                           signif_threshold=0.05,
                                           overwrite=False):
    """
    Test whether peak-center enrichment differs between BED x bigWig groups.

    Uses the original (non-renormalized) computeMatrix values supplied to
    plotHeatmap. Comparisons are performed separately for each metaGrp (target)
    from the metadata file. Within a metaGrp, every combination of BED region
    group and bigWig sample defines a group (e.g.
    BED_DNTTIP1._.sample_DNTTIP1_DMSO_D0). For each group, per-region mean signal
    across the center bin +/- flank is computed; pairwise Mann-Whitney U tests
    (two-sided) yield a symmetric p-value matrix. A ward-clustered seaborn
    clustermap of -log10(p) is saved alongside each TSV; non-significant cells
    (p >= signif_threshold) are shown in grey while clustering uses all values.
    """
    if not os.path.exists(matFile):
        if lgr:
            lgr.warning(f"Peak-center BED comparison skipped: matrix not found ({matFile})")
        return

    header, coords, X = loadDeepToolsMatrix(matFile, lgr=lgr)

    bed_group_boundaries = header["group_boundaries"]
    n_bed_groups = len(bed_group_boundaries) - 1
    if n_bed_groups < 2:
        if lgr:
            lgr.info("Peak-center BED comparison skipped: matrix has fewer than 2 region groups.")
        return

    sample_boundaries = header["sample_boundaries"]
    segmentLength = sample_boundaries[1] - sample_boundaries[0]
    n_samples = len(sample_boundaries) - 1

    if regionsLabels and regionsLabels != ["auto"]:
        bed_labels = [str(l).replace("\n", " ").strip() for l in regionsLabels[:n_bed_groups]]
    elif header.get("group_labels"):
        bed_labels = [str(l).replace("\n", " ").strip() for l in header["group_labels"][:n_bed_groups]]
    else:
        bed_labels = [f"bedGroup{i}" for i in range(n_bed_groups)]

    center_offset = segmentLength // 2
    lo = max(0, center_offset - center_bin_flank)
    hi = min(segmentLength - 1, center_offset + center_bin_flank)

    metaGrp2sample_idx = OrderedDict()
    for si in range(n_samples):
        mg = metaGrp[si] if si < len(metaGrp) else header["sample_labels"][si]
        metaGrp2sample_idx.setdefault(mg, []).append(si)

    for mg, sample_indices in metaGrp2sample_idx.items():
        combo_keys = []
        combo_values = OrderedDict()

        for gi in range(n_bed_groups):
            row_start, row_end = bed_group_boundaries[gi], bed_group_boundaries[gi + 1]
            for si in sample_indices:
                sample_label = (
                    samplesList[si] if si < len(samplesList) else header["sample_labels"][si]
                ).replace("\n", " ").strip()
                combo_key = peak_center_combo_label(bed_labels[gi], sample_label)
                combo_keys.append(combo_key)

                s, e = sample_boundaries[si], sample_boundaries[si + 1]
                center_cols = list(range(s + lo, s + hi + 1))

                if row_end <= row_start:
                    combo_values[combo_key] = np.array([], dtype=np.float32)
                else:
                    block = X[row_start:row_end, center_cols]
                    region_means = np.nanmean(block, axis=1)
                    combo_values[combo_key] = region_means[~np.isnan(region_means)]

        if len(combo_keys) < 2:
            if lgr:
                lgr.info(
                    f"Peak-center BED comparison skipped for metaGrp '{mg}': "
                    f"need at least 2 BED x sample groups (found {len(combo_keys)})."
                )
            continue

        safe_mg = sanitize_label_for_filename(mg)
        out_tsv = f"{outPrefix}.peakCenterBEDcomparison.binFlank_{center_bin_flank}.{safe_mg}.pvalues.tsv"
        out_pdf = f"{outPrefix}.peakCenterBEDcomparison.binFlank_{center_bin_flank}.{safe_mg}.pvalues.clustermap.pdf"

        if os.path.exists(out_tsv) and os.path.exists(out_pdf) and not overwrite:
            if lgr:
                lgr.info(f"Peak-center BED comparison exists, skipping metaGrp '{mg}': {out_tsv}")
            continue

        pval_mat = pd.DataFrame(np.nan, index=combo_keys, columns=combo_keys)
        for ki in range(len(combo_keys)):
            pval_mat.iloc[ki, ki] = 1.0

        for i, j in combinations(range(len(combo_keys)), 2):
            ki, kj = combo_keys[i], combo_keys[j]
            v1, v2 = combo_values[ki], combo_values[kj]
            if len(v1) == 0 or len(v2) == 0:
                pval = np.nan
            else:
                try:
                    _, pval = mannwhitneyu(v1, v2, alternative="two-sided")
                except ValueError:
                    pval = np.nan
            pval_mat.iloc[i, j] = pval
            pval_mat.iloc[j, i] = pval

        pval_mat.to_csv(out_tsv, sep="\t", float_format="%.6g")

        neglog = -np.log10(pval_mat.astype(float).clip(lower=1e-300, upper=1.0))
        np.fill_diagonal(neglog.values, 0.0)
        neglog = neglog.fillna(0.0)

        pval_arr = pval_mat.astype(float).values
        significant = (pval_arr < signif_threshold) & ~np.isnan(pval_arr)
        np.fill_diagonal(significant, False)
        color_mask = ~significant

        cmap = matplotlib.cm.get_cmap("viridis_r").copy()
        cmap.set_bad(color="#d9d9d9")

        n_combos = len(combo_keys)
        figsize = max(8.0, min(24.0, 0.45 * n_combos + 4.0))
        annot = n_combos <= 20
        g = sns.clustermap(
            neglog,
            mask=color_mask,
            method="ward",
            cmap=cmap,
            annot=annot,
            fmt=".2f",
            annot_kws={"size": max(4, int(8 - n_combos / 15))},
            linewidths=0.5,
            figsize=(figsize, figsize),
            cbar_kws={"label": f"-log10(p-value); grey if p >= {signif_threshold}"},
        )
        g.ax_heatmap.set_xticklabels(
            g.ax_heatmap.get_xticklabels(), rotation=45, horizontalalignment="right"
        )
        g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_yticklabels(), rotation=0)
        plt.savefig(out_pdf, bbox_inches="tight", dpi=300)
        plt.close(g.fig)

        if lgr:
            counts = ", ".join(f"{k}={len(combo_values[k])}" for k in combo_keys)
            lgr.info(
                f"Peak-center BED comparison for metaGrp '{mg}' "
                f"({len(combo_keys)} BED x sample groups, center bins {lo}-{hi} of {segmentLength}; "
                f"regions: {counts}) -> {out_tsv}, {out_pdf}"
            )


def writeClusterBEDs(coordsDf: pd.DataFrame, labels: np.ndarray, outPrefix: str, lgr=None):
    """
    Write one BED per cluster. Returns list of BED paths in cluster order (C1..Ck).
    """
    uniq = np.unique(labels)
    bedPaths = []
    outdir = os.path.dirname(outPrefix)
    if outdir and not os.path.exists(outdir):
        os.makedirs(outdir, exist_ok=True)
    for i, clu in enumerate(uniq, start=1):
        sub = coordsDf[labels == clu][["chrm","start","end"]]
        bed = f"{outPrefix}.cluster{i}.bed"
        sub.to_csv(bed, sep="\t", header=False, index=False)
        bedPaths.append(bed)
        if lgr: lgr.info(f"Wrote {len(sub)} regions to {bed}")
    return bedPaths


def plotDendrogram(Z, outPrefix: str, lgr=None):
    """Save dendrogram as PNG and PDF."""
    fig = plt.figure(figsize=(10, 6))
    sch.dendrogram(Z, no_labels=True, count_sort=True)
    plt.tight_layout()
    png = f"{outPrefix}.dendrogram.png"
    pdf = f"{outPrefix}.dendrogram.pdf"
    plt.savefig(png, dpi=200, bbox_inches="tight")
    plt.savefig(pdf, dpi=200, bbox_inches="tight")
    plt.close(fig)
    if lgr: lgr.info(f"Dendrogram saved: {png}, {pdf}")

def runPCAandTSNE(X: np.ndarray, labels: np.ndarray, outPrefix: str, seed: int, lgr=None, embedding_cache: dict | None = None):
    """
    Visualization only. If embedding_cache is provided, reuse:
      - 'pca_scores' + 'pca_evr' (+ 'pca_components') for PCA artifacts
      - 'tsne' for t-SNE scatter
    Otherwise compute as before.
    """
    def _auto_marker_size(n: int) -> int:
        # sizes are in points^2 (matplotlib scatter semantics)
        if n < 250:
            return 6
        elif n < 500:
            return 5
        elif n < 1000:
            return 4
        elif n < 2500:
            return 3
        elif n < 5000:
            return 2
        else:
            return 1
    
    def _pick_colors(n_clusters: int):
        # <=10 -> tab10, 11..20 -> tab20, >20 -> evenly-sampled nipy_spectral # more at https://matplotlib.org/stable/users/explain/colors/colormaps.html
        if n_clusters <= 10:
            cmap = matplotlib.cm.get_cmap("tab10", n_clusters)
        elif n_clusters <= 20:
            cmap = matplotlib.cm.get_cmap("tab20", n_clusters)
        else:
            cmap = matplotlib.cm.get_cmap("nipy_spectral", n_clusters)
        return [cmap(i) for i in range(n_clusters)]

    # Use cached PCA if available; else compute
    use_cached_pca = embedding_cache is not None and "pca_scores" in embedding_cache and "pca_evr" in embedding_cache
    comps = None

    if use_cached_pca:
        scores = embedding_cache["pca_scores"]
        evr = embedding_cache["pca_evr"]
        comps = embedding_cache.get("pca_components", None)
    else:
        Xs = StandardScaler(with_mean=True, with_std=True).fit_transform(X)
        npc = min(50, Xs.shape[1])
        pca = PCA(n_components=npc, svd_solver="auto", random_state=seed if isinstance(seed, int) else None)
        scores = pca.fit_transform(Xs)
        evr = pca.explained_variance_ratio_
        comps = pca.components_

    # Save PCA results (scores + scree)
    pd.DataFrame(scores[:, :2], columns=["PC1","PC2"]).assign(cluster=labels).to_csv(f"{outPrefix}.pca.scores.tsv", sep="\t", index=False)
    pd.DataFrame({"PC": np.arange(1, len(evr)+1), "explainedVariance": evr}).to_csv(f"{outPrefix}.pca.scree.tsv", sep="\t", index=False)

    # Optionally save PCA loadings if available
    if comps is not None and comps.ndim == 2 and comps.shape[0] >= 2:
        pd.DataFrame(comps[:2, :], index=["PC1","PC2"]).to_csv(f"{outPrefix}.pca.loadings.tsv", sep="\t", index=True, header=False)

    # Scree plot
    fig = plt.figure(figsize=(6,4))
    plt.bar(np.arange(1, len(evr)+1), evr)
    plt.plot(np.arange(1, len(evr)+1), evr, marker="o", linewidth=1)
    plt.xlabel("Principal Component")
    plt.ylabel("Proportion of variance")
    plt.tight_layout()
    plt.savefig(f"{outPrefix}.pca.scree.png", dpi=200)
    plt.savefig(f"{outPrefix}.pca.scree.pdf", dpi=200)
    plt.close(fig)

    # PCA scatter
    fig = plt.figure(figsize=(6, 6))
    ax = fig.add_subplot(111)

    uniq = np.unique(labels)
    colors = _pick_colors(len(uniq))
    ms = _auto_marker_size(scores.shape[0])  # dynamic marker size

    for i, clu in enumerate(uniq):
        sel = labels == clu
        ax.scatter(
            scores[sel, 0], scores[sel, 1],
            s=ms, marker="o", linewidths=0,
            label=f"C{int(i+1)}", alpha=0.8, color=colors[i]
        )

    ax.set_xlabel(f"PC1 ({evr[0]*100:.1f}%)")
    ax.set_ylabel(f"PC2 ({evr[1]*100:.1f}%)")

    # Legend outside (right)
    ax.legend(
        title="Clusters",
        bbox_to_anchor=(1.02, 1), loc="upper left", borderaxespad=0.0,
        frameon=False, fontsize=8, markerscale=1.5
    )

    fig.tight_layout()
    fig.savefig(f"{outPrefix}.pca.scatter.png", dpi=200, bbox_inches="tight")
    fig.savefig(f"{outPrefix}.pca.scatter.pdf", dpi=200, bbox_inches="tight")
    plt.close(fig)

    # t-SNE: reuse if present; else compute on PCA-50
    if embedding_cache is not None and "tsne" in embedding_cache:
        emb = embedding_cache["tsne"]
    else:
        X50 = scores[:, :min(50, scores.shape[1])]
        tsne = TSNE(n_components=2, init="pca", learning_rate="auto",
                    random_state=seed if isinstance(seed, int) else None, perplexity=30)
        emb = tsne.fit_transform(X50)
    pd.DataFrame(emb, columns=["TSNE1","TSNE2"]).assign(cluster=labels).to_csv(f"{outPrefix}.tsne.tsv", sep="\t", index=False)

    fig = plt.figure(figsize=(6, 6))
    ax = fig.add_subplot(111)

    ms_tsne = _auto_marker_size(emb.shape[0])  # dynamic marker size
    # reuse the same cluster colors so PCA and t-SNE match
    for i, clu in enumerate(uniq):
        sel = labels == clu
        ax.scatter(
            emb[sel, 0], emb[sel, 1],
            s=ms_tsne, marker="o", linewidths=0,
            label=f"C{int(i+1)}", alpha=0.8, color=colors[i]
        )

    ax.set_xlabel("t-SNE 1")
    ax.set_ylabel("t-SNE 2")

    # Legend outside (right)
    ax.legend(
        title="Clusters",
        bbox_to_anchor=(1.02, 1), loc="upper left", borderaxespad=0.0,
        frameon=False, fontsize=8, markerscale=1.5
    )

    fig.tight_layout()
    fig.savefig(f"{outPrefix}.tsne.scatter.png", dpi=200, bbox_inches="tight")
    fig.savefig(f"{outPrefix}.tsne.scatter.pdf", dpi=200, bbox_inches="tight")
    plt.close(fig)

    if lgr: lgr.info(f"PCA/t-SNE artifacts saved with prefix: {outPrefix}")

def runUMAP(labels: np.ndarray, outPrefix: str, lgr=None, embedding_cache: dict | None = None):
    """
    Generate UMAP scatter + artifacts if an embedding was computed.
    """
    if embedding_cache is None or "umap" not in embedding_cache:
        if lgr: lgr.warning("UMAP embedding not available in cache; skipping UMAP plots.")
        return
    emb = embedding_cache["umap"]
    uniq = np.unique(labels)

    # Save coordinates
    pd.DataFrame(emb, columns=["UMAP1", "UMAP2"]).assign(cluster=labels).to_csv(
        f"{outPrefix}.umap.tsv", sep="\t", index=False
    )

    # Save parameters for reproducibility
    params = embedding_cache.get("umap_params", None)
    if params is not None:
        with open(f"{outPrefix}.umap.params.json", "w") as fh:
            json.dump(params, fh, indent=2)

    # Color selection consistent with PCA/t-SNE
    def _auto_marker_size(n: int) -> int:
        if n < 250: return 6
        elif n < 500: return 5
        elif n < 1000: return 4
        elif n < 2500: return 3
        elif n < 5000: return 2
        else: return 1

    def _pick_colors(n_clusters: int):
        if n_clusters <= 10:
            cmap = matplotlib.cm.get_cmap("tab10", n_clusters)
        elif n_clusters <= 20:
            cmap = matplotlib.cm.get_cmap("tab20", n_clusters)
        else:
            cmap = matplotlib.cm.get_cmap("nipy_spectral", n_clusters)
        return [cmap(i) for i in range(n_clusters)]

    colors = _pick_colors(len(uniq))
    ms = _auto_marker_size(emb.shape[0])

    fig = plt.figure(figsize=(6, 6))
    ax = fig.add_subplot(111)

    for i, clu in enumerate(uniq):
        sel = labels == clu
        ax.scatter(emb[sel, 0], emb[sel, 1], s=ms, label=f"C{clu}", alpha=0.85, edgecolors="none", c=[colors[i]])

    ax.set_xlabel("UMAP1", fontsize=12)
    ax.set_ylabel("UMAP2", fontsize=12)
    ax.set_title("UMAP embedding (clusters)", fontsize=13)
    ax.legend(title="Clusters", bbox_to_anchor=(1.02, 1), loc="upper left", borderaxespad=0.0, frameon=False, fontsize=8, markerscale=1.5)
    fig.tight_layout()
    fig.savefig(f"{outPrefix}.umap.scatter.png", dpi=200, bbox_inches="tight")
    fig.savefig(f"{outPrefix}.umap.scatter.pdf", dpi=200, bbox_inches="tight")
    plt.close(fig)
    if lgr: lgr.info(f"Wrote UMAP artifacts: {outPrefix}.umap.*")

def compute_clustering_embedding(X: np.ndarray, mode: str, seed, lgr=None):
    """
    Returns:
      X_for_cluster: np.ndarray used for clustering
      emb_cache: dict with precomputed embeddings for plotting:
        - 'pca_scores': PCA scores (n x k)
        - 'pca_evr': explained variance ratio (len k)
        - 'pca_components': PCA loadings (k x d)
        - 'tsne': t-SNE (n x 2), when mode == 'tsne'
    """
    emb_cache = {}

    if mode == "raw":
        if lgr: lgr.info("Clustering space: RAW matrix")
        return X, None

    # Standardize first
    Xs = StandardScaler(with_mean=True, with_std=True).fit_transform(X)

    # PCA (min(50, d))
    npc = min(50, Xs.shape[1])
    pca = PCA(n_components=npc, svd_solver="auto", random_state=seed if isinstance(seed, int) else None)
    Xp = pca.fit_transform(Xs)
    evr = pca.explained_variance_ratio_
    emb_cache["pca_scores"] = Xp
    emb_cache["pca_evr"] = evr
    emb_cache["pca_components"] = pca.components_

    if mode == "pca":
        if lgr: lgr.info(f"Clustering space: PCA ({npc} comps), EVR first 3: {np.round(evr[:3], 4)}")
        return Xp, emb_cache

    elif mode == "tsne":

        # t-SNE(2D) on PCA-50
        tsne = TSNE(n_components=2, init="pca", learning_rate="auto",
                    random_state=seed if isinstance(seed, int) else None, perplexity=30)
        Xt = tsne.fit_transform(Xp)
        emb_cache["tsne"] = Xt
        if lgr: lgr.info("Clustering space: t-SNE (2D on PCA-50)")
        return Xt, emb_cache

    elif mode == "umap":
        if UMAP is None:
            raise ImportError("umap-learn is not available; install it and re-run (e.g., conda install -c conda-forge umap-learn).")
        # Reuse cached PCA if available, else compute like other branches did
        if "pca_scores" in emb_cache:
            Xp = emb_cache["pca_scores"]
        else:
            Xs = StandardScaler(with_mean=True, with_std=True).fit_transform(X)
            npc = min(50, Xs.shape[1])
            pca = PCA(n_components=npc, random_state=seed if isinstance(seed, int) else None)
            Xp = pca.fit_transform(Xs)
            evr = pca.explained_variance_ratio_
            emb_cache["pca_scores"] = Xp
            emb_cache["pca_evr"] = evr
            emb_cache["pca_components"] = pca.components_
        reducer = UMAP(
            n_components=2,
            n_neighbors=15,
            min_dist=0.1,
            metric="euclidean",
            random_state=seed if isinstance(seed, int) else None,
        )
        Xu = reducer.fit_transform(Xp)
        emb_cache["umap"] = Xu
        emb_cache["umap_params"] = {
            "n_components": 2,
            "n_neighbors": 15,
            "min_dist": 0.1,
            "metric": "euclidean",
            "random_state": int(seed) if isinstance(seed, int) else None,
        }
        if lgr: lgr.info("Clustering space: UMAP (2D on PCA-50)")
        return Xu, emb_cache



def externalCluster_and_Replot(matPath: str,
                               basePrefix: str,
                               method: str,
                               clusterNumber: int,
                               rndSeed: str,
                               heatmapHeight: int,
                               refPointLabel: str,
                               xAxisLabel: str,
                               colorsFlag: str,
                               colorMaps: str,
                               sortUsingSamplesString: str,
                               samples_for_recompute, labelsHeatmap_for_recompute, labelsHeatmapPrefix_for_recompute, samplesList, conditionGrp, metaGrp, profileColors,
                               lgr,
                               reuseRenormedMatrix: bool,
                               clusteringDimReduct: str,
                               leidenResolutionParameter: float,
                               renormApproach: str,
                               clusterUsingSamples: str = "",
                               start_point: str = "center",
                               bodyLength: int = 1000):
    """
    Driver for external clustering + deepTools reorder + plotting.
    method in {'hclustExt','kmeansExt','miniBatchKMeansExt','leidenExt'}.
    Produces:
      - cluster BEDs
      - reordered/labeled matrix
      - heatmap with --sortRegions keep
      - PCA/t-SNE for kmeans*/leiden; dendrogram for hclustExt
    """
    header, coords, X = loadDeepToolsMatrix(matPath, lgr)
    seed = setRandomSeed(rndSeed, lgr)

    # Decide which bins (columns) are used in clustering
    sample_labels = header.get("sample_labels", [])
    boundaries = header.get("sample_boundaries", [])
    n_samples = len(sample_labels)

    def _log_range(i1, start, end, include=True):
        lab = sample_labels[i1 - 1] if 1 <= i1 <= n_samples else f"idx{i1}"
        msg = ("include" if include else "exclude")
        # half-open [start:end), show both range and size
        if lgr: lgr.info(f"{msg} bins [{start}:{end - 1}] (n={end - start}) -> sample '{lab}' (index {i1})")

    if clusterUsingSamples:
        # Parse 1-based indices "1,3,5"
        raw_parts = [x.strip() for x in clusterUsingSamples.split(",") if x.strip()]
        try:
            sel_idx = sorted({int(x) for x in raw_parts})
        except Exception:
            raise ValueError(f"--clusterUsingSamples must be a comma-separated list of 1-based integers, got: {clusterUsingSamples}")

        bad = [i for i in sel_idx if i < 1 or i > n_samples]
        if bad:
            raise ValueError(f"--clusterUsingSamples contains out-of-range indices (1..{n_samples}): {bad}")

        selected_cols = []
        if lgr: lgr.info(f"Clustering on subset of samples (1-based): {sel_idx}")
        for i1 in sel_idx:
            s, e = boundaries[i1 - 1], boundaries[i1]
            selected_cols.extend(range(s, e))
            _log_range(i1, s, e, include=True)

        # Log the excluded ranges for transparency
        excluded = [i for i in range(1, n_samples + 1) if i not in sel_idx]
        for i1 in excluded:
            s, e = boundaries[i1 - 1], boundaries[i1]
            _log_range(i1, s, e, include=False)

        X_used = X[:, selected_cols]
        if lgr: lgr.info(f"Total bins selected for clustering: {X_used.shape[1]} out of {X.shape[1]}")
    else:
        # Default: use all samples, but still log the mapping
        if lgr: lgr.info("Clustering on ALL samples (default). Included bin ranges per sample:")
        for i1 in range(1, n_samples + 1):
            s, e = boundaries[i1 - 1], boundaries[i1]
            _log_range(i1, s, e, include=True)
        X_used = X

    # Making sure there are no NaNs or Infs in the matrix:
    X_used = np.nan_to_num(X_used, nan=0.0, posinf=0.0, neginf=0.0)


    # Choose feature space for clustering
    X_for_cluster, emb_cache = compute_clustering_embedding(X_used, clusteringDimReduct, seed, lgr)

    # Compute labels
    if method == "hclustExt":
        if fastcluster is None or sch is None:
            raise ImportError("fastcluster/scipy not available; install to use hclustExt.")
        if clusterNumber < 1:
            raise ValueError("clusterNumber must be >=1 for hclustExt.")
        Z = fastcluster.linkage_vector(X_for_cluster, method="ward")  # Euclidean
        labels = sch.fcluster(Z, t=clusterNumber, criterion="maxclust").astype(int)
    elif method == "kmeansExt":
        if clusterNumber < 1:
            raise ValueError("clusterNumber must be >=1 for kmeansExt.")
        km = KMeans(n_clusters=clusterNumber, n_init=25, random_state=seed)
        labels = km.fit_predict(X_for_cluster) + 1
    elif method == "miniBatchKMeansExt":
        if clusterNumber < 1:
            raise ValueError("clusterNumber must be >=1 for miniBatchKMeansExt.")
        mb = MiniBatchKMeans(n_clusters=clusterNumber, batch_size=4096, n_init=10, random_state=seed)
        labels = mb.fit_predict(X_for_cluster) + 1
    elif method == "leidenExt":
        if ig is None or la is None:
            raise ImportError("igraph/leidenalg not available; install to use leidenExt.")
        # Build kNN graph (undirected) with weights ~ similarity
        n_neighbors = max(5, min(15, X_for_cluster.shape[0]-1))
        nn = NearestNeighbors(n_neighbors=n_neighbors, metric="euclidean")
        nn.fit(X_for_cluster)
        dists, idxs = nn.kneighbors(X_for_cluster, return_distance=True)
        edges, weights = [], []
        n = X_for_cluster.shape[0]
        for i in range(n):
            for j, dist in zip(idxs[i,1:], dists[i,1:]):  # skip self
                if i < j:
                    edges.append((i, int(j)))
                    weights.append(1.0 / (1.0 + float(dist)))
        g = ig.Graph(n=n, edges=edges, directed=False)
        g.es["weight"] = weights

        part = la.find_partition(
            g,
            la.RBConfigurationVertexPartition,
            weights="weight",
            seed=int(seed) if isinstance(seed, int) else None,
            resolution_parameter=float(leidenResolutionParameter),
        )
        labels = np.array(part.membership, dtype=int) + 1
        lgr.info(f"leidenExt produced {len(np.unique(labels))} clusters (resolution={leidenResolutionParameter}); --clusterNumber is ignored for this method.")

    elif method == "leidenKmeansExt":
        if ig is None or la is None:
            raise ImportError("igraph/leidenalg not available; install to use leidenKmeansExt.")
        if clusterNumber < 1:
            raise ValueError("clusterNumber must be >=1 for leidenKmeansExt.")

        # 1) Leiden (same kNN + weights construction as above)
        n_neighbors = max(5, min(15, X_for_cluster.shape[0]-1))
        nn = NearestNeighbors(n_neighbors=n_neighbors, metric="euclidean")
        nn.fit(X_for_cluster)
        dists, idxs = nn.kneighbors(X_for_cluster, return_distance=True)
        edges, weights = [], []
        n = X_for_cluster.shape[0]
        for i in range(n):
            for j, dist in zip(idxs[i,1:], dists[i,1:]):  # skip self
                if i < j:
                    edges.append((i, int(j)))
                    weights.append(1.0 / (1.0 + float(dist)))
        g = ig.Graph(n=n, edges=edges, directed=False)
        g.es["weight"] = weights

        part = la.find_partition(
            g,
            la.RBConfigurationVertexPartition,
            weights="weight",
            seed=int(seed) if isinstance(seed, int) else None,
            resolution_parameter=float(leidenResolutionParameter),
        )
        le_labels = np.array(part.membership, dtype=int)  # 0..(L-1)
        L = int(le_labels.max()) + 1
        lgr.info(f"leidenKmeansExt: Leiden found {L} communities at resolution={leidenResolutionParameter}")

        # 2) If requested K >= L, just keep Leiden as-is. Nothing to merge.
        if clusterNumber >= L:
            labels = le_labels + 1  # 1..L
            lgr.info(
                f"leidenKmeansExt: requested K={clusterNumber} ≥ Leiden communities L={L}; "
                f"using Leiden labels as-is (no k-means merge)."
            )
        else:
            # 3) Merge communities down to exactly K via KMeans on community centroids
            centroids = np.zeros((L, X_for_cluster.shape[1]), dtype=float)
            for c in range(L):
                sel = (le_labels == c)
                # guard if any empty (shouldn't happen, but be safe)
                if np.any(sel):
                    centroids[c, :] = X_for_cluster[sel, :].mean(axis=0)
                else:
                    centroids[c, :] = 0.0

            km = KMeans(n_clusters=clusterNumber, n_init=50, random_state=seed)
            super_assign = km.fit_predict(centroids)  # assigns each Leiden community to one of K merged clusters (0..K-1)

            # Map each point's Leiden community to its merged (kmeans) super-cluster
            merged = np.array([super_assign[c] for c in le_labels], dtype=int)
            labels = merged + 1  # 1..K

            lgr.info(
                f"leidenKmeansExt: merged L={L} Leiden communities down to K={clusterNumber} via KMeans-on-centroids."
            )

    else:
        raise ValueError(f"Unsupported external method: {method}")
    lgr.info(f"Clustering method={method} on space={clusteringDimReduct}; clusters={len(np.unique(labels))}")

    # BEDs and group labels
    outPrefix = f"{basePrefix}.{method}"
    bedPaths = writeClusterBEDs(coords, labels, outPrefix, lgr)
    groupLabels = [f"C{i+1}" for i in range(len(bedPaths))]
    regionsLabelsFlag = " --regionsLabel " + " ".join(f"'{gl}'" for gl in groupLabels)

    # Recompute a deepTools matrix for groups = cluster BEDs
    clusterMat, groupLabels = computeMatrix_for_cluster_BEDs(
        bedPaths=bedPaths,
        outPrefix=outPrefix,              # e.g. {analysisPrefix}.renorm.leidenExt
        flank=header["upstream"][0],      
        n_proc=int(header.get("proc number", 1)),
        samples=samples_for_recompute,   
        labelsHeatmap=labelsHeatmap_for_recompute,     
        labelsHeatmapPrefix=labelsHeatmapPrefix_for_recompute, 
        sortUsingSamplesString=sortUsingSamplesString,
        lgr=lgr,
        reuseRenormedMatrix=reuseRenormedMatrix,
        start_point = start_point,
        bodyLength=bodyLength
    )

    # Renormalize this NEW matrix and plot two heatmaps
    renormOutPrefix = f"{outPrefix}.clustersFromRenorm"
    renormMat = renormalize_matrix_and_plot(
        matFile=clusterMat,
        outPrefix=renormOutPrefix,
        samplesList=samplesList,
        conditionGrp=conditionGrp,
        metaGrp=metaGrp,
        profileColors=profileColors,
        heatmapHeightBase=heatmapHeight,      # we already computed H in main; pass it in
        colorsFlag=colorsFlag,
        colorMaps=colorMaps,
        refPointLabel=refPointLabel,
        xAxisLabel=xAxisLabel,
        sortUsingSamplesString=sortUsingSamplesString,
        regionsLabelsFlag=regionsLabelsFlag,
        lgr=lgr,
        reuseRenormedMatrix=reuseRenormedMatrix,
        renormApproach=renormApproach
    )

    # Extra visualizations
    if method == "hclustExt":
        plotDendrogram(Z, outPrefix, lgr)
    else:
        runPCAandTSNE(X_used, labels, outPrefix, seed, lgr, embedding_cache=emb_cache)
        if emb_cache is not None and "umap" in emb_cache:
            runUMAP(labels, outPrefix, lgr, embedding_cache=emb_cache)

    # Save cluster label table
    rid = coords[["chrm","start","end"]].astype(str).agg(lambda r: f"{r['chrm']}:{r['start']}-{r['end']}", axis=1)
    pd.DataFrame({"region": rid, "cluster": labels}).to_csv(f"{outPrefix}.clusters.tsv", sep="\t", index=False)

    return renormMat

def computeMatrix_for_cluster_BEDs(bedPaths,
                                   outPrefix,
                                   flank,
                                   n_proc,
                                   samples,                 # space-separated bigWig paths
                                   labelsHeatmap,           # "$lab_0" "$lab_1" ...
                                   labelsHeatmapPrefix,     # "lab_0=$'..' && lab_1=$'..' && "
                                   sortUsingSamplesString,  # artifact, not really used but kept for consistency
                                   lgr,
                                   reuseRenormedMatrix: bool,
                                   start_point: str = "center",
                                   bodyLength = 0):
    """
    Re-run computeMatrix using cluster BEDs as groups (order = cluster1, cluster2, ...).
    Returns path to the new matrix and the group labels we used.
    """
    # C1, C2, ...
    groupLabels = [f"C{i+1}" for i in range(len(bedPaths))]
    # regionsLabelsFlag = " --regionsLabel " + " ".join(f"'{gl}'" for gl in groupLabels)

    outMat = f"{outPrefix}.clusters.computeMatrix.mat.gz"

    # Reuse if requested and existing
    if reuseRenormedMatrix and os.path.exists(outMat):
        lgr.info(f"[reuseRenormedMatrix=True] Reusing existing clustered computeMatrix: {outMat}")
        return outMat, groupLabels
    else:
        # regionsLabelsFlag = " --regionsLabel " + " ".join(f"'{gl}'" for gl in groupLabels)
        cmd = (
            f"""{labelsHeatmapPrefix} """
            f"""computeMatrix reference-point -p {n_proc} --referencePoint {start_point} """
            f"""-S {samples} --samplesLabel {labelsHeatmap} """
            f"""-R {' '.join(bedPaths)} """
            f"""--beforeRegionStartLength {flank} --afterRegionStartLength {flank} """
            f"""--skipZeros --missingDataAsZero """
            f"""-o {outMat}""" # {regionsLabelsFlag} # we need to add labels at plotHeatmap level
        )
        if bodyLength > 0 and start_point != "center":
            cmd = (
                f"""{labelsHeatmapPrefix} """
                f"""computeMatrix scale-regions -p {n_proc} --startLabel "start" --endLabel "end" """
                f"""-S {samples} --samplesLabel {labelsHeatmap} """
                f"""-R {' '.join(bedPaths)} """
                f"""--regionBodyLength {bodyLength} """
                f"""--beforeRegionStartLength {flank} --afterRegionStartLength {flank} """
                f"""--skipZeros --missingDataAsZero """
                f"""-o {outMat}""" # {regionsLabelsFlag} # we need to add labels at plotHeatmap level
            )

        lgr.info(f"Recomputing matrix per clusters: {cmd}")
        proc = subprocess.run(cmd, shell=True, capture_output=True, text=True)
        if proc.returncode != 0 or not os.path.exists(outMat):
            lgr.error("computeMatrix for clusters failed.")
            lgr.error(f"STDOUT:\n{proc.stdout}\nSTDERR:\n{proc.stderr}")
            raise RuntimeError("computeMatrix (clusters) failed.")
        return outMat, groupLabels

def _quantile_normalize_2d(M: np.ndarray) -> np.ndarray:
    """
    Quantile-normalize columns of a 2D array (features x samples).
    NaNs are treated as 0s for stability. Returns a new array.
    """
    if M.ndim != 2:
        raise ValueError("Expected 2D array for quantile normalization.")
    X = np.array(M, dtype=np.float64, copy=True)
    # Safety for NaNs
    X[np.isnan(X)] = 0.0

    # Sort each column
    X_sorted = np.sort(X, axis=0)  # shape: (n_features, n_samples)
    # Mean of ranks across columns
    rank_means = np.mean(X_sorted, axis=1)

    # For each column, map original values to rank_means by argsort ranks
    X_norm = np.empty_like(X)
    # argsort twice gives ranks 0..n_features-1; ties get increasing order (OK in practice)
    ranks = np.argsort(np.argsort(X, axis=0), axis=0)
    # Broadcast rank_means per column
    n_samples = X.shape[1]
    for j in range(n_samples):
        X_norm[:, j] = rank_means[ranks[:, j]]
    return X_norm


def _qnorm_per_meta_group_on_matrix(enrichDf: pd.DataFrame,
                                    header: dict,
                                    samplesList: list[str],
                                    metaGrp: list[str],
                                    lgr=None) -> pd.DataFrame:
    """
    Quantile-normalize the matrix per meta-group across *samples*.
    Implementation detail:
      - For each metaGrp, collect sample blocks via header['sample_boundaries'].
      - For each sample in the group, flatten its (n_regions x n_bins) block into a 1D vector.
      - Stack vectors => (n_regions*n_bins, n_samples_in_group), quantile-normalize columns.
      - Reshape each normalized column back to (n_regions x n_bins) and write into original columns.
    Returns a NEW DataFrame with qnorm-applied values.
    """
    # Build mapping: sample index -> [start,end) in column-space
    boundaries = header["sample_boundaries"]  # e.g. [0, 400, 800, ...]
    n_samples = len(boundaries) - 1
    if n_samples != len(samplesList):
        if lgr: lgr.warning(f"Header sample count ({n_samples}) != samplesList length ({len(samplesList)}). Assuming header order is truth.")
    sample_blocks = [(boundaries[i], boundaries[i+1]) for i in range(n_samples)]

    # metaGrp -> indices of samples in this group (by header order)
    metaGrp2idxs = OrderedDict()
    for idx, mg in enumerate(metaGrp):
        metaGrp2idxs.setdefault(mg, []).append(idx)

    # We'll create a copy to avoid mutating caller's df until end
    outDf = enrichDf.copy()

    n_regions = enrichDf.shape[0]
    for grp, idxs in metaGrp2idxs.items():
        if len(idxs) <= 1:
            # nothing to normalize across: one sample -> skip
            if lgr: lgr.info(f"[qnorm] metaGrp='{grp}' has {len(idxs)} sample(s); skipping qnorm for this group.")
            continue

        # Gather flattened columns per sample (in this metaGrp)
        flat_cols = []
        blocks = []  # store (start,end) to write back
        for sidx in idxs:
            s, e = sample_blocks[sidx]
            block_vals = enrichDf.loc[:, list(range(s, e))].to_numpy(dtype=np.float64, copy=False)  # (n_regions, n_bins)
            flat_cols.append(block_vals.reshape(-1))  # (n_regions*n_bins,)
            blocks.append((s, e))
        M = np.column_stack(flat_cols)  # shape: (n_regions*n_bins, n_group_samples)

        # Quantile normalize across columns
        Mq = _quantile_normalize_2d(M)

        # Write back: reshape each column to (n_regions, n_bins) and assign
        for col_j, (s, e) in enumerate(blocks):
            n_bins = e - s
            restored = Mq[:, col_j].reshape(n_regions, n_bins)
            outDf.loc[:, list(range(s, e))] = restored

        if lgr: lgr.info(f"[qnorm] Applied quantile normalization for metaGrp='{grp}' across {len(blocks)} sample(s).")
    return outDf


def renormalize_matrix_and_plot(matFile,
                                outPrefix,
                                samplesList,
                                conditionGrp,
                                metaGrp,
                                profileColors,
                                heatmapHeightBase,
                                colorsFlag,
                                colorMaps,
                                refPointLabel,
                                xAxisLabel,
                                sortUsingSamplesString,
                                regionsLabelsFlag,
                                lgr,
                                reuseRenormedMatrix: bool,
                                renormApproach: str):
    
    # ---------- read header + prepare ----------
    topColumns = 2
    headerString = None
    if testFileGzippness(matFile):
        fh = gzip.open(matFile, "rb")
    else:
        fh = open(matFile, "rb")
    lines = list(islice(fh, topColumns))
    fh.close()

    switch = 0
    for row in lines:
        row = row.decode("ascii")
        if row.startswith("@"):
            headerString = row
            header = json.loads(row[1:].strip())
            segmentLength = header["sample_boundaries"][1]
            bwsNum = len(header["sample_labels"])
            switch = 1
        elif switch == 1:
            rowLength = len(row.strip().split("\t"))

    if headerString is None:
        raise RuntimeError("Header not found while renormalizing.")

    # ---------- load enrichment block ----------
    enrichDf = pd.read_csv(
        matFile,
        skiprows=topColumns-1,
        header=None,
        compression="infer",
        sep="\t",
        usecols=list(range(-6, rowLength))[rowLength - header["sample_boundaries"][-1]:],
        names=["chrm","start","end","name","score","strand"] + list(range(header["sample_boundaries"][-1])),
        index_col=["chrm","start","end","name","score","strand"]
    )

    # map bins -> sample
    binEdge2sample = OrderedDict()
    for be, sample in zip(header["sample_boundaries"], header["sample_labels"]):
        binEdge2sample[be] = sample
    bin2sample = OrderedDict()
    for b in range(header["sample_boundaries"][-1]):
        if b in binEdge2sample:
            currentSample = binEdge2sample[b]
        bin2sample[b] = currentSample

    # ---------- mean profiles (your current approach) ----------
    enrichMean = enrichDf.mean()
    if isinstance(header["bin size"], list):
        binSize = header["bin size"][0]
    else:
        binSize = header["bin size"]

    dataDictSimple = OrderedDict()
    for sample, enrStart, enrEnd in zip(samplesList, header['sample_boundaries'][:-1], header['sample_boundaries'][1:]):
        dataDictSimple[sample] = list(enrichMean[enrStart:enrEnd])
    dataDictSimple["idx"] = list(range(int(-(segmentLength*binSize)/2), int((segmentLength*binSize)/2), binSize))
    dataDf = pd.DataFrame.from_dict(dataDictSimple).set_index("idx")
    dataDf.to_csv(f"{outPrefix}.MeanSampleEnrichment.tsv", sep="\t")
    lgr.info(f"Saved {outPrefix}.MeanSampleEnrichment.tsv")

    # ---------- min-max per metaGrp on profile level (your current approach) ----------
    samples2metaGrp = {s: m for s, m in zip(samplesList, metaGrp)}
    metaGrp2samples = OrderedDict()
    for s, m in zip(samplesList, metaGrp):
        metaGrp2samples.setdefault(m, []).append(s)

    dataRenormDf = dataDf.copy()
    minMaxPerGrp = {}
    for grp, cols in metaGrp2samples.items():
        ddmin = dataRenormDf[cols].min().min()
        ddmax = dataRenormDf[cols].max().max()
        minMaxPerGrp[grp] = (ddmin, ddmax)
        for sample in cols:
            dataRenormDf[sample] = dataRenormDf[sample].apply(lambda x: (x - ddmin) / (ddmax - ddmin))
    dataRenormDf.to_csv(f"{outPrefix}.MinMaxPerGrpSampleEnrichment.tsv", sep="\t")

    # ---------- decide reuse vs recompute ----------
    renormMat = f"{outPrefix}.renormMatrix.mat.gz"
    if reuseRenormedMatrix and os.path.exists(renormMat):
        lgr.info(f"[reuseRenormedMatrix=True] Reusing existing renormalized matrix: {renormMat}")
        # We still want profile/plots; read header/enrichment from the *input* matFile for the profile,
        # then skip the renormalization write step and go straight to plotting using renormMat.
        skip_renorm = True
    else:
        lgr.info(f"[reuseRenormedMatrix={reuseRenormedMatrix}] Will (re)compute renormalized matrix: {renormMat}")
        skip_renorm = False
        
    # ---------- write renormalized matrix (if needed) ----------
    if not skip_renorm:
        tmpMat = renormMat.replace(".gz", "")
        with open(tmpMat, "w") as out:
            out.write(headerString)

        if renormApproach == "profileMinMax":
            # apply same min-max to every bin by its sample's metaGrp (current/default behavior)
            lgr.info("[renorm] Using profileMinMax (per-metaGrp min-max on profile level).")
            for b in bin2sample:
                mn, mx = minMaxPerGrp[samples2metaGrp[bin2sample[b]]]
                enrichDf[b] = enrichDf[b].apply(lambda x: (x - mn) / (mx - mn))

            outDf = enrichDf  # already transformed

        elif renormApproach == "qnorm":
            lgr.info("[renorm] Using qnorm (per-metaGrp quantile normalization of matrix values).")
            # Perform quantile normalization across samples within each meta-group
            outDf = _qnorm_per_meta_group_on_matrix(
                enrichDf=enrichDf,
                header=header,
                samplesList=samplesList,
                metaGrp=metaGrp,
                lgr=lgr
            )

            try:
                samples2meta = {s: m for s, m in zip(samplesList, metaGrp)}
                # Build sample blocks again
                boundaries = header["sample_boundaries"]
                n_samples_hdr = len(boundaries) - 1
                for si in range(n_samples_hdr):
                    s, e = boundaries[si], boundaries[si+1]
                    mg = samples2meta[samplesList[si]] if si < len(samplesList) else None
                    if mg is None:
                        continue
                    block = outDf.loc[:, list(range(s, e))]
                    # Min-max within metaGrp across all its samples (same min/max as profile approach)
                    # Compute per-group min/max once:
                # precompute:
                metaGrp2samples = OrderedDict()
                for s, m in zip(samplesList, metaGrp):
                    metaGrp2samples.setdefault(m, []).append(s)
                grp_minmax = {}
                for m, cols in metaGrp2samples.items():
                    # build block indices for all samples in this metaGrp
                    idxs = []
                    for si, sname in enumerate(samplesList):
                        if sname in cols:
                            idxs.extend(list(range(boundaries[si], boundaries[si+1])))
                    vals = outDf.loc[:, idxs].to_numpy()
                    grp_minmax[m] = (float(np.min(vals)), float(np.max(vals)))
                # rescale by group's min-max
                for si in range(n_samples_hdr):
                    s, e = boundaries[si], boundaries[si+1]
                    m = samples2meta[samplesList[si]] if si < len(samplesList) else None
                    if m is None:
                        continue
                    mn, mx = grp_minmax[m]
                    if mx > mn:
                        outDf.loc[:, list(range(s, e))] = (outDf.loc[:, list(range(s, e))] - mn) / (mx - mn)
            except Exception as _e:
                lgr.warning(f"[renorm][qnorm] Optional post-qnorm [0,1] rescale skipped due to: {_e}")

        else:
            raise ValueError(f"Unknown renormApproach: {renormApproach}")

        outDf.to_csv(tmpMat, mode="a", header=False, sep="\t")
        subprocess.run(f"gzip -f {tmpMat}", shell=True, check=True)
        lgr.info(f"Renormalized matFile and saved as {renormMat} (approach={renormApproach})")

    # ---------- plotting (heatmaps + profile) ----------
    # heatmap height: one row per cluster group
    n_groups = len(header["group_boundaries"]) - 1
    pxH = heatmapHeightBase * max(1, n_groups)

    # default-scale heatmap
    cmd1 = (
        f"plotHeatmap --missingDataColor '#ffffff' -m {renormMat} "
        f"--heatmapHeight {pxH} --interpolationMethod bilinear "
        f"--refPointLabel '{refPointLabel}' --xAxisLabel '{xAxisLabel}' "
        f"-out {outPrefix}.deeptools.renormMatrix.pdf "
        f"--zMin 0 "
        f"--{colorsFlag} {colorMaps} {sortUsingSamplesString} " # --sortRegions keep
        f"{regionsLabelsFlag}"
    )
    lgr.info(f"Plotting heatmap (renorm): {cmd1}")
    subprocess.Popen(cmd1, shell=True, stdin=None, stdout=None, stderr=None)

    # bounded heatmap
    cmd2 = (
        f"plotHeatmap --missingDataColor '#ffffff' -m {renormMat} "
        f"--heatmapHeight {pxH} --interpolationMethod bilinear "
        f"--refPointLabel '{refPointLabel}' --xAxisLabel '{xAxisLabel}' "
        f"-out {outPrefix}.deeptools.renormMatrix.zMax_1.zMin_0.pdf "
        f"--zMax 1.1 --zMin 0 "
        f"--{colorsFlag} {colorMaps} {sortUsingSamplesString}" # --sortRegions keep
        f" {regionsLabelsFlag}"
    )
    lgr.info(f"Plotting heatmap (renorm z-bounded): {cmd2}")
    subprocess.Popen(cmd2, shell=True, stdin=None, stdout=None, stderr=None)

    # profile
    samples2color = OrderedDict()
    if 'auto' in profileColors:
        pal = sns.color_palette("hls", len(samplesList))
        for s, c in zip(samplesList, pal):
            samples2color[s] = c
    else:
        for s, c in zip(samplesList, profileColors):
            samples2color[s] = c

    dataDf_reset = dataDf.reset_index().rename(columns={"idx": "idx"})
    plotProfile(dataDf_reset, outPrefix, samplesList, samples2color,
                numberOfPeaks=header["group_boundaries"][-1])
    return renormMat


def main():
    configureLogging()
    lgr = logging.getLogger(inspect.currentframe().f_code.co_name)
    regionsList, regionsFileName, regionsLabels, analysisPrefix, metaDataFile, flank, n_proc, refPointLabel, xAxisLabel, heatmapHeightBase, sortRegions, sortUsing, sortUsingSamples, whatToShow, useColorList, extras, overwrite, diffPeakMode, diffPeakThresholds, maxRegionsPG, inputType, clusterRegions, gBinSize, clusterNumber, rndSeed, reuseRenormedMatrix, clusteringDimReduct, leidenResolutionParameter, renormApproach, clusterUsingSamples, start_point, bodyLength, peakCenterBinFlank, peakCenterSignifThreshold = parseArgs()
    
    ### Prepare the BW list, labels and other params:
    df = readbwMeta(metaDataFile,lgr)
    lgr.info(f"Identified {len(df)} input samples (bw files).")
#     lab_0=$'H3K27Ac_WdrKO_r\nTreg_rep1'
#     lab_1=$'H3K27Ac_WdrKO_r\nTreg_rep2'

    if clusterRegions != "none":
        if inputType != "region":
            lgr.critical("--clusterRegions is only valid with '--inputType region'.")
            raise Exception("Clustering requested in 'bin' mode.")
        if len(regionsList) != 1:
            lgr.critical(f"Clustering works on ONE region BED only; you supplied {len(regionsList)} files.")
            raise Exception("Too many region files for clustering.")

    permissionIssue = 0
    samples = ""
    labels = ""
    labelsHeatmap = ""
    labelsHeatmapPrefix = ""
    colorMaps = ""
    conditionGrp = list(df["conditionGrp"])
    metaGrp = list(df["metaGrp"])
    samplesList = list(df["bigWigLabel"])
    #samplesList = [re.sub("(.{15})", "\\1"+"\n", col.replace("'","").replace("#",""), 0, re.DOTALL) for col in samplesList]
    samplesList = [ autobreaklinetitle(col) for col in samplesList ]
    profileColors = list(df["profileColor"])
    for idx, row in df.iterrows():
        if os.path.exists(row["bigWig"]) == False:
            permissionIssue = 1
            lgr.error(f"BigWig file {row['bigWig']} is not accessible >>> exiting")
        else:
            samples += f""" {row["bigWig"]}"""
            labels += f""" '{row["bigWigLabel"].replace(" ","_")}'"""
            #labelsHeatmapPrefix += "lab_{}=$'{}' && ".format(idx, re.sub("(.{15})", "\\1"+r"\\n", row["bigWigLabel"], 0, re.DOTALL))
            labelsHeatmapPrefix += "lab_{}=$'{}' && ".format(idx, autobreaklinetitle(row["bigWigLabel"]))
            labelsHeatmap += """ "$lab_{}" """.format(idx)
            colorMaps += f""" '{row["colorMap"]}'"""
    if permissionIssue == 1:
        msg = "At least one of the input BigWig files was not accessible, program was aborted."
        lgr.error(msg)
        raise Exception(msg)
        
    ### Extract peaks from voom file if working as a part of differential peak calling pipeline:
    if diffPeakMode == True:
        voomDf = regionsList.copy()
        voomDf["chrm"] = voomDf["Region"].apply(lambda x: x.split(":")[0])
        voomDf["start"] = voomDf["Region"].apply(lambda x: int(x.split(":")[1].split("-")[0]))
        voomDf["end"] = voomDf["Region"].apply(lambda x: int(x.split(":")[1].split("-")[1]))
        
        ctrlCat = "Control"
        if "Up2FC_Down2FC" in diffPeakThresholds:
            upCat = "Up2FC"
            downCat = "Down2FC"
            upColor = "#d73502"
            downColor = "#0066f2"
            statType = 'p.value'
            statVal = 1.1
            FCcut = 1
        elif "Up2_Down2" in diffPeakThresholds:
            upCat = "Up2"
            downCat = "Down2"
            upColor = "#c63254"
            downColor = "#005073"
            statType = 'q.value'
            statVal = 0.05
            FCcut = 1
        elif "Up2NoFDR_Down2NoFDR" in diffPeakThresholds:
            upCat = "Up2NoFDR"
            downCat = "Down2NoFDR"
            upColor = "#ec6449"
            downColor = "#107dac"
            statType = 'p.value'
            statVal = 0.05
            FCcut = 1
        elif "Up_Down" in diffPeakThresholds:
            upCat = "Up"
            downCat = "Down"
            upColor = "#f9994b"
            downColor = "#189ad3"
            statType = 'q.value'
            statVal = 0.05
            FCcut = 0
        elif "UpNoFDR_DownNoFDR" in diffPeakThresholds:
            upCat = "UpNoFDR"
            downCat = "DownNoFDR"
            upColor = "#fbcb5a"
            downColor = "#1ebbd7"
            statType = 'p.value'
            statVal = 0.05
            FCcut = 0
        elif diffPeakThresholds == "Control":
            upCat = downCat = "NA"
            upColor = "#fbcb5a"
            downColor = "#1ebbd7"
            statType = 'p.value'
            statVal = 1
            FCcut = 0
            
            prefix = f'{analysisPrefix}.Control'
            regionsList = []
            regionsLabels = []
        else:
            msg = "Something went wrong, currently no other options for diffPeakThresholds are supported, except of: Up2FC_Down2FC, Up2_Down2, Up2NoFDR_Down2NoFDR, Up_Down, UpNoFDR_DownNoFDR; + (_)Control. Program execution was aborted."
            lgr.error(msg)
            raise Exception(msg)
        ctrlColor = "#808080"
        
        ctrlDf = voomDf[voomDf["Regulation"] == ctrlCat].copy()
        ctrlDf = ctrlDf[["chrm", "start", "end"]].copy()
        optimalCtrl = len(ctrlDf)
        
        if "Up" in diffPeakThresholds and "Down" in diffPeakThresholds:
            upDf = voomDf[(voomDf[statType] < statVal) & (voomDf["log2FC"] > FCcut)].copy()
            upDf = upDf[["chrm", "start", "end"]].copy()
            downDf = voomDf[(voomDf[statType] < statVal) & (voomDf["log2FC"] < -FCcut)].copy()
            downDf = downDf[["chrm", "start", "end"]].copy()
        
            if len(ctrlDf) <= 2:
                ### if there is not enough control peaks, we take top 10% of the least significantly changing sites.
                ctrlDf = voomDf.copy()
                ctrlDf["absLog2FC"] = ctrlDf["log2FC"].apply(lambda x: abs(x))
                ctrlDf["negPVal"] = ctrlDf["p.value"].apply(lambda x: -x)
                ctrlDf["negFDR"] = ctrlDf["q.value"].apply(lambda x: -x)
                ctrlDf.sort_values(by=["absLog2FC", "negPVal", "negFDR"], ascending=True, inplace=True)
                ctrlDf = ctrlDf.head(int(len(ctrlDf)/10))
                ctrlDf = ctrlDf[["chrm", "start", "end"]].copy()
        
            ### Find the optimal number of control peaks to display (we dont want to simply put out all of them, because usually their volume would overwhelm the plot, thus we wish to balance them roughly equally with other traces)
            upNum = len(upDf)
            downNum = len(downDf)
            ctrlNum = len(ctrlDf)
            if upNum > 10 and downNum > 10:
                optimalCtrl = np.min([int((upNum+downNum)/2), ctrlNum])
            elif upNum > 10 and downNum <= 10:
                optimalCtrl = np.min([upNum, ctrlNum])
            elif upNum <= 10 and downNum > 10:
                optimalCtrl = np.min([downNum, ctrlNum])
            else:
                optimalCtrl = ctrlNum
        
            ### Extract BED region files
            prefix = f'{analysisPrefix}.{upCat}_{downCat}'
            regionsList = []
            regionsLabels = []
            if upNum > 10:
                bedName = f"{prefix}.{upCat}.bed"
                upDf.to_csv(bedName, index=False, header=False, sep="\t")
                regionsList.append(bedName)
                regionsLabels.append(f"{upCat}\n[{upNum} regions]")
        
            if downNum > 10:
                bedName = f"{prefix}.{downCat}.bed"
                downDf.to_csv(bedName, index=False, header=False, sep="\t")
                regionsList.append(bedName)
                regionsLabels.append(f"{downCat}\n[{downNum} regions]")
        
        # optimalCtrl = 100 ### this is an easy way to force the algorithm to work with less regions in order to test it.
        if optimalCtrl > 1 and "Control" in diffPeakThresholds:
            bedName = f"{prefix}.{ctrlCat}.bed"
            ctrlDf = ctrlDf.sample(n=optimalCtrl)
            ctrlDf.to_csv(bedName, index=False, header=False, sep="\t")
            regionsList.append(bedName)
            regionsLabels.append(f"{ctrlCat}\n[{optimalCtrl} regions]")
    
    haltFurtherExecution = False
    permissionIssue = 0
    regions = ""
    if inputType == "region":
        if regionsList != []:
            for bed in regionsList:
                if os.path.exists(bed) == False:
                    permissionIssue = 1
                    lgr.error(f"BED file {bed} is not accessible >>> exiting")
                else:
                    regions += f" {bed}"
            if permissionIssue == 1:
                msg = "At least one of the input BED files was not accessible, program was aborted."
                lgr.error(msg)
                raise Exception(msg)
        else:
            lgr.warning(f"NO REGIONS FOUND TO PLOT, if working as a part of differential peak calling pipeline, it means that from the selected threshold, no differential peaks were found.")
            haltFurtherExecution = True
            #raise Exception("No regions found to plot in regionsList")
    else:
        pass

    if haltFurtherExecution == True:
        pass
    else:
    
        if len(samplesList) != len(set(samplesList)):
            msg = "The list of sample labels contain non unique names, correct that and re-run"
            lgr.error(msg)
            raise Exception(msg)
    
        pxH = heatmapHeightBase*len(regionsList)
        if diffPeakMode == True:
            pass#prefix = f'{analysisPrefix}.{upCat}_{downCat}'
        else:
            if inputType == "region":
                prefix = f'{analysisPrefix}.flank_{flank}'
            else:
                prefix = f'{analysisPrefix}.bin_{gBinSize}'
        matFile = f"{prefix}.mat.gz"
    
        ### Generate matrix for PCA and correlation from deeptools:
        if inputType == "region":
            if len(regionsList) > 1:
                mergedBed = BedTool(regionsList[0])
                for bed in regionsList[1:]:
                    mergedBed = mergedBed.cat(bed, postmerge=True)
                mergedBed.saveas(f'{prefix}.mergedRegions.bed')
                regionsInit = f'{prefix}.mergedRegions.bed'
                lgr.info(f"Merged the regions for deeptools correlation and correlation and saved as {regionsInit}")
            else:
                regionsInit = regionsList[0]

            mergedBed = BedTool(regionsInit)
        else:
            mergedBed = []
        
        if len(mergedBed) > 10 or inputType == "bin":
            npzFile = f"{prefix}.multiBigwigSummary.npz"
            npzTsvFile = f"{prefix}.multiBigwigSummary.tsv"
            if inputType == "bin":
                command  = (f"multiBigwigSummary bins -p {n_proc} --labels {labels} -b {samples} -bs {gBinSize} -o {npzFile} --outRawCounts {npzTsvFile}")
            else:
                command = f"""multiBigwigSummary BED-file -p {n_proc} --labels {labels} -b {samples} -o {npzFile} --BED {regionsInit} --outRawCounts {npzTsvFile} """
            print("command: ", command)
            if os.path.exists(npzFile) == False or overwrite == True:
                lgr.info(f"Calculate matrix with command: {command}")
                subprocess.Popen(command, shell=True, stdout=subprocess.PIPE).stdout.read()
            else:
                lgr.info(f"{npzFile} file exists, skipping matrix generation.")
            
            # Convert the table to log10 scale:
            log_tsv = f"{prefix}.multiBigwigSummary_log10.tsv"
            if not os.path.exists(log_tsv) or overwrite:
                df_raw = pd.read_csv(npzTsvFile, sep="\t")
                # first 3 columns are genomic coords per deepTools docs
                # https://deeptools.readthedocs.io/en/develop/content/tools/multiBigwigSummary.html
                coord_cols = list(df_raw.columns[:3])
                sample_cols = list(df_raw.columns[3:])
                # transform only sample columns
                df_raw[sample_cols] = np.log10(df_raw[sample_cols].astype(float) + 1.0)
                df_raw.to_csv(log_tsv, sep="\t", index=False)
        
            ### Generate actually nice looking PairGrid plot:
            pariGridFilePrefix = f"{prefix}.PairGrid"
            plotPairGrid(npzTsvFile, pariGridFilePrefix, maxRegionsPG, overwrite=overwrite, lgr=lgr)
        
            ### Get the list of samples to be sorted by:
            sortUsingSamplesString = ""
            if sortUsingSamples != " ":
                sortUsingSamplesString = "--sortUsingSamples {}".format(' '.join(str(x) for x in sortUsingSamples.split(",")))
        
            ### Generate PCA and correlation plots:
            if os.path.exists(npzFile) == True and (extras == True or diffPeakMode == True):
                command = f"""{labelsHeatmapPrefix} plotCorrelation -in {npzFile} --labels {labelsHeatmap} --corMethod pearson --removeOutliers --skipZeros --plotTitle 'Pearson Correlation of Average Enrichment Scores\n{os.path.basename(prefix)}' --whatToPlot scatterplot --log1p -o {prefix}.deeptoolsPearsonCorr.scatter.png --plotHeight 50 --plotWidth 50 & """
                lgr.info(f"Plotting deeptools Pearson correlation scatter with command: {command}")
        #         subprocess.Popen(command, shell=True, stdout=subprocess.PIPE).stdout.read()
                subprocess.Popen(command, shell=True, stdin=None, stdout=None, stderr=None)
            
                command = f"""{labelsHeatmapPrefix} plotCorrelation -in {npzFile} --labels {labelsHeatmap} --corMethod spearman --skipZeros --plotTitle 'Spearman Correlation of Average Enrichment Scores\n{os.path.basename(prefix)}' --whatToPlot scatterplot --log1p -o {prefix}.deeptoolsSpearmanCorr.scatter.png --plotHeight 50 --plotWidth 50 & """
                lgr.info(f"Plotting deeptools Spearman correlation scatter with command: {command}")
        #         subprocess.Popen(command, shell=True, stdout=subprocess.PIPE).stdout.read()
                subprocess.Popen(command, shell=True, stdin=None, stdout=None, stderr=None)
            
                command = f"""plotPCA -in {npzFile}  -o {prefix}.deeptoolsPCA.pdf -T 'PCA - {os.path.basename(prefix)}' --plotHeight 15 --plotWidth 10  --transpose & """
                lgr.info(f"Plotting deeptools PCA with command: {command}")
        #         subprocess.Popen(command, shell=True, stdout=subprocess.PIPE).stdout.read()
                subprocess.Popen(command, shell=True, stdin=None, stdout=None, stderr=None)
            
                command = f"""plotCorrelation -in {npzFile}  --corMethod spearman --skipZeros --plotTitle 'Spearman Correlation of Average Enrichment Scores\n{os.path.basename(prefix)}' --whatToPlot heatmap --colorMap bwr --plotNumbers -o {prefix}.deeptoolsSpearmanCorr.pdf --outFileCorMatrix {prefix}.deeptoolsSpearmanCorrValues.tab --plotHeight 16 --plotWidth 12 & """
                lgr.info(f"Plotting deeptools Spearman correlation with command: {command}")
        #         subprocess.Popen(command, shell=True, stdout=subprocess.PIPE).stdout.read()
                subprocess.Popen(command, shell=True, stdin=None, stdout=None, stderr=None)
            
                command = f"""plotCorrelation -in {npzFile}  --corMethod pearson --removeOutliers --skipZeros --plotTitle 'Pearson Correlation of Average Enrichment Scores\n{os.path.basename(prefix)}' --whatToPlot heatmap --colorMap bwr --plotNumbers -o {prefix}.deeptoolsPearsonCorr.pdf --outFileCorMatrix {prefix}.deeptoolsPearsonCorrValues.tab --plotHeight 16 --plotWidth 12 & """
                lgr.info(f"Plotting deeptools Pearson correlation with command: {command}")
        #         subprocess.Popen(command, shell=True, stdout=subprocess.PIPE).stdout.read()
                subprocess.Popen(command, shell=True, stdin=None, stdout=None, stderr=None)

                ### Generate correlation heatmap for log10 transformed data (using custom plotting approach for this one only because the input is not the npz file but the tsv file):
                df_log = pd.read_csv(log_tsv, sep="\t")
                sample_cols = list(df_log.columns[3:])
                corr_log = df_log[sample_cols].corr(method="pearson")
                corr_log.to_csv(f"{prefix}.deeptoolsPearsonCorrValues.log10.tab", sep="\t")

                g = sns.clustermap(
                    corr_log, method="ward", cmap="bwr", vmin=-1, vmax=1,
                    annot=True, annot_kws={"size": 7}, linewidths=.5, figsize=(12,12)
                )
                g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xticklabels(), rotation=45, ha='right')
                plt.savefig(f"{prefix}.deeptoolsPearsonCorr.log10.pdf", bbox_inches='tight', dpi=300)
                plt.close()



        #         command = f"""plotCorrelation -in {npzFile}  --corMethod spearman --skipZeros --plotTitle 'Spearman Correlation of Average Scores\n{os.path.basename(prefix)}' --whatToPlot heatmap --colorMap bwr -min -1 -max 1 --plotNumbers -o {prefix}.deeptoolsSpearmanCorr.minMaxSet.pdf --plotHeight 60 --plotWidth 50 & """
        #         lgr.info(f"Plotting deeptools Spearman correlation with command: {command}")
        # #         subprocess.Popen(command, shell=True, stdout=subprocess.PIPE).stdout.read()
        #         subprocess.Popen(command, shell=True, stdin=None, stdout=None, stderr=None)
        #         
        #         command = f"""plotCorrelation -in {npzFile}  --corMethod pearson --skipZeros --plotTitle 'Pearson Correlation of Average Scores\n{os.path.basename(prefix)}' --whatToPlot heatmap --colorMap bwr -min -1 -max 1 --plotNumbers -o {prefix}.deeptoolsPearsonCorr.minMaxSet.pdf --plotHeight 60 --plotWidth 50 & """
        #         lgr.info(f"Plotting deeptools Pearson correlation with command: {command}")
        # #         subprocess.Popen(command, shell=True, stdout=subprocess.PIPE).stdout.read()
        #         subprocess.Popen(command, shell=True, stdin=None, stdout=None, stderr=None)
            
            else:
                if extras == True:
                    msg = f"Input Matrix file {npzFile} does not exist - PCA and correlation plots will not be generated"
                    lgr.error(msg)
                    raise Exception(msg)
        

            
        
        
            ### Generate heatmap plot:
            if inputType == "region":
                ### Generate heatmaps
                ### Calculate matrix (or skip):
                command = f"""{labelsHeatmapPrefix} computeMatrix reference-point -p {n_proc} --referencePoint {start_point} -S {samples} --samplesLabel {labelsHeatmap} -R {regions} --beforeRegionStartLength {flank} --afterRegionStartLength {flank} --skipZeros --missingDataAsZero -o {matFile}"""
                if os.path.exists(matFile) == False or overwrite == True:
                    lgr.info(f"Calculate matrix with command: {command}")
                    subprocess.Popen(command, shell=True, stdout=subprocess.PIPE).stdout.read()
                else:
                    lgr.info(f"{matFile} file exists, skipping matrix generation.")
                    
                if useColorList == True:
                    colorsFlag = "colorList"
                else:
                    colorsFlag = "colorMap"
            
                if regionsLabels == ['auto']:
                    if regions.count(" ") == 1:
                        regionsLabelsFlag = r" --regionsLabel ' ' "
                    else:
                        regionsLabelsFlag = r" --regionsLabel "
                        for i in range(regions.count(" ")):
                            regionsLabelsFlag += f"'bed{i}' "
                        #regionsLabelsFlag = r" --regionsLabel {}".format(regions.count(" ") * "' ' ")
                else:
                    regionsLabelsFlag = " --regionsLabel {}".format(' '.join(f'"{x}"' for x in regionsLabels))
                if clusterRegions != "none":
                    lgr.info("Clustering enabled (--clusterRegions); omitting --regionsLabel so deepTools names clusters automatically.")
                    regionsLabelsFlag = ""
            
                if os.path.exists(matFile) == True:
                    clusterFlag = ""
                    if clusterRegions == "kmeans":
                        clusterFlag = f"--kmeans {clusterNumber}"
                    elif clusterRegions == "hclust":
                        clusterFlag = f"--hclust {clusterNumber}"
                    
                    command = f"""plotHeatmap --missingDataColor '#ffffff' -m {matFile} --heatmapHeight {pxH} --interpolationMethod bilinear --refPointLabel '{refPointLabel}' --xAxisLabel '{xAxisLabel}' --outFileSortedRegions {prefix}.sortedRegions.txt -out {prefix}.deeptools.pdf --{colorsFlag} {colorMaps} {sortUsingSamplesString} {regionsLabelsFlag} {clusterFlag}""" # Removing from params: --outFileNameMatrix {prefix}.deeptoolsHeatmap.outMatrix.tsv.gz
                    lgr.info(f"Plotting heatmap with command: {command}")
                    subprocess.Popen(command, shell=True, stdin=None, stdout=None, stderr=None)
                
                    command = f"""plotHeatmap --missingDataColor '#ffffff' -m {matFile} --heatmapHeight {pxH} --interpolationMethod bilinear --refPointLabel '{refPointLabel}' --xAxisLabel '{xAxisLabel}' -out {prefix}.deeptools.svg --{colorsFlag} {colorMaps} {sortUsingSamplesString} {regionsLabelsFlag} {clusterFlag}""" # Removing from params: --outFileNameMatrix {prefix}.deeptoolsHeatmap.outMatrix.tsv.gz
                    lgr.info(f"Plotting heatmap with command: {command}")
                    subprocess.Popen(command, shell=True, stdin=None, stdout=None, stderr=None)
                
                    command = f"""plotHeatmap --missingDataColor '#ffffff' -m {matFile} --heatmapHeight {pxH} --interpolationMethod bilinear --refPointLabel '{refPointLabel}' --xAxisLabel '{xAxisLabel}' -out {prefix}.deeptools.png --{colorsFlag} {colorMaps} {sortUsingSamplesString} {regionsLabelsFlag} {clusterFlag}""" # Removing from params: --outFileNameMatrix {prefix}.deeptoolsHeatmap.outMatrix.tsv.gz
                    lgr.info(f"Plotting heatmap with command: {command}")
                    subprocess.Popen(command, shell=True, stdin=None, stdout=None, stderr=None)
                else:
                    msg = f"Input Matrix file {matFile} does not exist - plot will not be generated"
                    lgr.error(msg)
                    raise Exception(msg)
            else:
                lgr.info("Heatmap plotting is only supported for 'region' input type, skipping heatmap plotting.")
    
        else:
            lgr.info("Not enough regions to plot correlation and PCA")
    
    ### Convert matFile to other stuff:
    # Note, the issue and something to be resolved is that other "extras" plots are only generated while calculating the "{prefix}.renormMatrix.mat.gz" for the first time. After that they would no longer be re-calculated. this is and issue potentially.
    numberOfPeaks = -1
    if diffPeakMode == False and inputType == "region":
        if os.path.exists(matFile) and len(regionsList) >= 2 and clusterRegions == "none":
            compare_peak_center_between_bed_groups(
                matFile=matFile,
                outPrefix=prefix,
                samplesList=samplesList,
                metaGrp=metaGrp,
                regionsLabels=regionsLabels,
                lgr=lgr,
                center_bin_flank=peakCenterBinFlank,
                signif_threshold=peakCenterSignifThreshold,
                overwrite=overwrite,
            )

        ### We do not want to renormalize for the usage as a part of differential peak calling pipeline, because all the samples for that application by definition belong to the same renormalization "group"
        renormMatFileName = f"{prefix}.renormMatrix.mat.gz"
        if os.path.exists(renormMatFileName) == False or overwrite == True:
    
            topColumns=2 # Number of top rows to read in and to establish the total number of "column" in the structure of the input file. Its needed because first 3 or 6 (or more ?) columns in the file (except header(s)) will contain region coordinates, and the rest will be the enrichment signals. 

            if testFileGzippness(matFile) == True:
                infile = gzip.open(matFile, 'rb')
            else:
                infile = open(matFile, 'r')
            infileSlice = list(islice(infile, topColumns))
            infile.close()

            switch = 0
            for row in infileSlice:
                row = row.decode('ascii')
                if row[0] == "@":
                    headerString = row
                    ### The header looks e.g. like this, and if it doesn't, this means you might have used different deeptools version (in older ones the header was more complicate, there were additional column header line etc.)
                    # @{"verbose":false,"scale":1,"skip zeros":true,"nan after end":false,"sort using":"mean","unscaled 5 prime":0,"body":0,"sample_labels":["H3K27Ac::N6_WT-Rep1","H3K2
                    # 7Ac::N6_WT-Rep2","H3K27Ac::N6_WT-Rep3","H3K27Ac::N6_IKN159Y-Rep1","H3K27Ac::N6_IKN159Y-Rep2","H3K27Ac::N6_IKN159Y-Rep3","H3K7Me::WT-Rep1","H3K7Me::WT-Rep2","H3K7M
                    # e::WT-Rep3","H3K7Me::Mut-Rep1","H3K7Me::Mut-Rep2","H3K7Me::Mut-Rep3","WTIK::N6_WT-Rep1","WTIK::N6_WT-Rep2","WTIK::N6_WT-Rep3","IKN159Y::N6_IKN159Y-Rep1","IKN159Y:
                    # :N6_IKN159Y-Rep2","IKN159Y::N6_IKN159Y-Rep3"],"downstream":10000,"unscaled 3 prime":0,"group_labels":["K27ac-H3K27Ac-N6_IKN159Y.Repro_n48595"],"bin size":50,"upst
                    # ream":10000,"group_boundaries":[0,48503],"sample_boundaries":[0,400,800,1200,1600,2000,2400,2800,3200,3600,4000,4400,4800,5200,5600,6000,6400,6800,7200],"max thre
                    # shold":2000.0,"ref point":["center","center","center","center","center","center","center","center","center","center","center","center","center","center","center",
                    # "center","center","center"],"min threshold":null,"sort regions":"keep","proc number":1,"bin avg type":"mean","missing data as zero":true}
                    header = json.loads(row[1:].strip())
                    segmentLength = header["sample_boundaries"][1]
                    bwsNum = len(header['sample_labels'])
                    switch = 1
                elif switch == 1:
                    rowLength = len(row.strip().split("\t"))
            #samples = header["sample_labels"]
            numberOfPeaks = header["group_boundaries"][-1]

            ### reopen the input file using Pandas, and focusing only on the "important" rows and columns, which only contain enrichment signals:
            enrichDf = pd.read_csv( matFile,
                                    skiprows=topColumns-1,
                                    header=None,
                                    compression="infer",
                                    sep="\t",
                                    usecols=list(range(-6,rowLength))[rowLength-header["sample_boundaries"][-1]:],
                                    names=["chrm", "start", "end", "name", "score", "strand"] + list(range(0,header["sample_boundaries"][-1])),
                                    index_col=["chrm", "start", "end", "name", "score", "strand"]
                                    )
            binEdge2sample = OrderedDict()
            for be, sample in zip(header["sample_boundaries"], header["sample_labels"]):
                binEdge2sample[be] = sample
            bin2sample = OrderedDict()
            for b in list(range(0,header["sample_boundaries"][-1])):
                if b in binEdge2sample:
                    currentSample = binEdge2sample[b]
                bin2sample[b] = currentSample

            ## Get mean values per column (i.e. the aggregated line plots from the top part of deeptools heatmap); then reorganize and save as dataMatrix for future reference.
    
            enrichMean = enrichDf.mean()
            dataDictAll = OrderedDict()
            dataDictSimple = OrderedDict()
            if type(header["bin size"]) == type([]):
                binSize = header["bin size"][0]
            else:
                binSize = header["bin size"]
            for sample, condition, meta, enrStart, enrEnd in zip(samplesList, conditionGrp, metaGrp, header['sample_boundaries'][:-1], header['sample_boundaries'][1:]):
                dataDictAll[sample] = [sample, condition, meta, enrStart, enrEnd, enrichMean[enrStart:enrEnd]]
                dataDictSimple[sample] = list(enrichMean[enrStart:enrEnd])
            dataDictSimple["idx"] = list(range(int(-(segmentLength*binSize)/2), int((segmentLength*binSize)/2), binSize))

            dataDf = pd.DataFrame.from_dict(dataDictSimple)
            dataDf.set_index("idx", inplace=True)
            dataDf.to_csv(f"{prefix}.MeanSampleEnrichment.tsv", sep="\t")
            lgr.info(f"Saved {prefix}.MeanSampleEnrichment.tsv")

            dataDf = pd.read_csv(f"{prefix}.MeanSampleEnrichment.tsv", sep="\t", index_col="idx")

            ### Gather other global data:
            samples2condGrp = OrderedDict()
            samples2metaGrp = OrderedDict()
            condGrp2samples = OrderedDict()
            metaGrp2samples = OrderedDict()
            dataRenormDf = dataDf.copy()
            for s, c, m, p in zip(samplesList, conditionGrp, metaGrp, profileColors):
                samples2condGrp[s] = c
                samples2metaGrp[s] = m

                if c not in condGrp2samples:
                    condGrp2samples[c] = []
                condGrp2samples[c].append(s)

                if m not in metaGrp2samples:
                    metaGrp2samples[m] = []
                metaGrp2samples[m].append(s)

            

            ### Min-Max normalize data per metaGrp (i.e. its usually antibody) on the profile level:
            minMaxPerGrp = {}
            for grp in metaGrp2samples:
                cols_to_norm = metaGrp2samples[grp]
                ddmin = dataRenormDf[cols_to_norm].min().min()
                ddmax = dataRenormDf[cols_to_norm].max().max()
                minMaxPerGrp[grp] = (ddmin, ddmax)
                for sample in cols_to_norm:
                    dataRenormDf[sample] = dataRenormDf[sample].apply(lambda x: (x-ddmin)/(ddmax-ddmin))

            if extras == True:
                dataRenormDf.to_csv(f"{prefix}.MinMaxPerGrpSampleEnrichment.tsv", sep="\t")
                lgr.info(f"Saved {prefix}.MinMaxPerGrpSampleEnrichment.tsv")


            ### Generate plots again:
            # if extras == True:
            # # if len(samplesList) > 1:
        #         plotCorrMatrix(dataRenormDf, f'{analysisPrefix}.MinMaxPerGrp', samplesList, conditionGrp, metaGrp, samples2condGrp, samples2metaGrp)
        #         plotPCA(dataRenormDf, f'{analysisPrefix}.MinMaxPerGrp', samplesList, conditionGrp, metaGrp, samples2condGrp, samples2metaGrp)
                # plotRibbon(dataRenormDf, f'{prefix}.MinMaxPerGrp', samplesList, conditionGrp, metaGrp, samples2condGrp, samples2metaGrp)
                # lgr.info(f"Plotted CorrMatrix, PCA and Ribbon plots (renormalized)")


            ### Min-Max normalize data per metaGrp (i.e. its usually antibody) on the matrix level:
    
    
            # add header to file:
            renormMatFile = open(renormMatFileName.replace('.renormMatrix.mat.gz','.renormMatrix.mat'), 'w')
            renormMatFile.write(headerString)
            renormMatFile.close()
    
            # renorm enrichDf:

            # print(samples2metaGrp)
            for b in bin2sample:
                enrichDf[b] = enrichDf[b].apply(lambda x: (x-minMaxPerGrp[samples2metaGrp[bin2sample[b]]][0]) / (minMaxPerGrp[samples2metaGrp[bin2sample[b]]][1] - minMaxPerGrp[samples2metaGrp[bin2sample[b]]][0]))
            enrichDf.to_csv(renormMatFileName.replace(".renormMatrix.mat.gz",".renormMatrix.mat"), mode='a', header=False, sep="\t")
    
            # gzip renormalized matrix
            command = f"gzip {renormMatFileName.replace('.renormMatrix.mat.gz','.renormMatrix.mat')}"
            subprocess.Popen(command, shell=True, stdout=subprocess.PIPE).stdout.read()
            lgr.info(f"Renormalized matFile and saved as {renormMatFileName}")
        else:
            lgr.info(f"{renormMatFileName} file exists, skipping matrix renormalization.")
            with open(npzTsvFile, 'r') as fp:
                for count, line in enumerate(fp):
                    pass
            numberOfPeaks = count-1 #print('Total Lines', count + 1)
            
            

        ### Generate profilePlot:
        samples2color = OrderedDict()
        if 'auto' in profileColors:
            profileColors = sns.color_palette("hls", len(samplesList))
        for s, c in zip(samplesList, profileColors):
            samples2color[s] = c

        dataDf = pd.read_csv(f"{prefix}.MeanSampleEnrichment.tsv", sep="\t")

        ### Generate an equivalent of `plotProfile` function from deeptools, but using our precalculated data from dataDf:
        plotProfile(dataDf, prefix, samplesList, samples2color, numberOfPeaks)
        lgr.info(f"Plotted profiles (not renormalized)")

        ### Generate plots:
        if extras == True:
            ### Gather other global data:
            samples2condGrp = OrderedDict()
            samples2metaGrp = OrderedDict()
            condGrp2samples = OrderedDict()
            metaGrp2samples = OrderedDict()
            dataRenormDf = dataDf.copy()
            for s, c, m, p in zip(samplesList, conditionGrp, metaGrp, profileColors):
                samples2condGrp[s] = c
                samples2metaGrp[s] = m

                if c not in condGrp2samples:
                    condGrp2samples[c] = []
                condGrp2samples[c].append(s)

                if m not in metaGrp2samples:
                    metaGrp2samples[m] = []
                metaGrp2samples[m].append(s)
            # if len(samplesList) > 1:
        #         plotCorrMatrix(dataDf, analysisPrefix, samplesList, conditionGrp, metaGrp, samples2condGrp, samples2metaGrp)
        #         plotPCA(dataDf, analysisPrefix, samplesList, conditionGrp, metaGrp, samples2condGrp, samples2metaGrp)
            #plotRibbon(dataDf, prefix, samplesList, conditionGrp, metaGrp, samples2condGrp, samples2metaGrp)
            lgr.info(f"Plotted CorrMatrix and PCA (not renormalized)")

        
        ### External clustering branch (raw + renorm)
        extMethods = {"hclustExt", "kmeansExt", "miniBatchKMeansExt", "leidenExt", "leidenKmeansExt"}
        if clusterRegions in extMethods:
            colorsFlag = "colorList" if useColorList else "colorMap"
            if clusterRegions == "leidenExt":
                lgr.info("leidenExt ignores --clusterNumber; cluster count is determined by Leiden.")

            # heatmap height (per group)
            H = heatmapHeightBase  # this is per group; final function multiplies by n_groups

            if os.path.exists(renormMatFileName):
                if clusterUsingSamples and clusterUsingSamples.strip():
                    # Normalize e.g. "1, 5" -> "1_5" (sorted, unique, 1-based kept)
                    toks = [t for t in re.split(r"[,\s]+", clusterUsingSamples.strip()) if t]
                    sel_idx = []
                    for t in toks:
                        try:
                            sel_idx.append(int(t))
                        except Exception:
                            # ignore non-integer tokens silently
                            continue
                    sel_idx = sorted(set(sel_idx))
                    clusterUsingSuffix = f"clusterUsingSamples_{'_'.join(str(i) for i in sel_idx)}" if sel_idx else "clusterUsingSamples_all"
                else:
                    clusterUsingSuffix = "clusterUsingSamples_all"

                outDir = (
                    f"{analysisPrefix}.clustering."
                    f"clusteringDimReduct_{clusteringDimReduct}."
                    f"clusterNumber_{clusterNumber}."
                    f"clusterRegions_{clusterRegions}."
                    f"renormApproach_{renormApproach}."
                    f"{clusterUsingSuffix}."
                    f"leidenResolutionParameter_{leidenResolutionParameter}"
                )
                os.makedirs(outDir, exist_ok=True)
                if lgr: lgr.info(f"Output directory: {outDir}")

                basePrefix = os.path.join(outDir, f"{analysisPrefix}.renorm")

                lgr.info(f"External clustering [{clusterRegions}] on RENORM matrix: {renormMatFileName}")
                externalCluster_and_Replot(
                    matPath=renormMatFileName,
                    basePrefix=basePrefix,
                    method=clusterRegions,
                    clusterNumber=clusterNumber,
                    rndSeed=rndSeed,
                    heatmapHeight=H,
                    refPointLabel=refPointLabel,
                    xAxisLabel=xAxisLabel,
                    colorsFlag=colorsFlag,
                    colorMaps=colorMaps,
                    sortUsingSamplesString=sortUsingSamplesString,
                    lgr=lgr,
                    samples_for_recompute=samples,
                    labelsHeatmap_for_recompute=labelsHeatmap,
                    labelsHeatmapPrefix_for_recompute=labelsHeatmapPrefix,
                    samplesList=samplesList, 
                    conditionGrp=conditionGrp, 
                    metaGrp=metaGrp, 
                    profileColors=profileColors,
                    reuseRenormedMatrix=reuseRenormedMatrix,
                    clusteringDimReduct=clusteringDimReduct,
                    leidenResolutionParameter=leidenResolutionParameter,
                    renormApproach=renormApproach,
                    clusterUsingSamples=clusterUsingSamples,
                    start_point=start_point,
                    bodyLength=bodyLength,
                )
            else:
                lgr.warning(f"Renormalized matrix not found; skipping external clustering: {renormMatFileName}")


        ### Generate renormalized heatmap plot:

        if os.path.exists(renormMatFileName) == True:
            clusterFlag = ""
            if clusterRegions == "kmeans":
                clusterFlag = f"--kmeans {clusterNumber}"
            elif clusterRegions == "hclust":
                clusterFlag = f"--hclust {clusterNumber}"

            if clusterRegions != "none":
                lgr.info("Clustering enabled (--clusterRegions); omitting --regionsLabel for renormalized heatmaps as well.")
                regionsLabelsFlag = ""

            command = f"""plotHeatmap --missingDataColor '#ffffff' -m {renormMatFileName} --heatmapHeight {pxH} --interpolationMethod bilinear --refPointLabel '{refPointLabel}' --xAxisLabel '{xAxisLabel}' --outFileSortedRegions {prefix}.sortedRegions.renormMatrix.txt -out {prefix}.deeptools.renormMatrix.pdf --{colorsFlag} {colorMaps} {sortUsingSamplesString} {regionsLabelsFlag} {clusterFlag}""" # removing from params --outFileNameMatrix {prefix}.deeptoolsHeatmap.outMatrix.renormMatrix.tsv.gz
            lgr.info(f"Plotting heatmap with command: {command}")
    #         subprocess.Popen(command, shell=True, stdout=subprocess.PIPE).stdout.read()
            subprocess.Popen(command, shell=True, stdin=None, stdout=None, stderr=None)
    
    #         command = f"""plotHeatmap --missingDataColor '#ffffff' -m {renormMatFileName} --heatmapHeight {pxH} --interpolationMethod bilinear --refPointLabel '{refPointLabel}' --xAxisLabel '{xAxisLabel}' -out {prefix}.deeptoolsHeatmap.renormMatrix.svg --{colorsFlag} {colorMaps} {sortUsingSamplesString} {regionsLabelsFlag}""" # removing from params --outFileNameMatrix {prefix}.deeptoolsHeatmap.outMatrix.renormMatrix.tsv.gz
    #         lgr.info(f"Plotting heatmap with command: {command}")
    # #         subprocess.Popen(command, shell=True, stdout=subprocess.PIPE).stdout.read()
    #         subprocess.Popen(command, shell=True, stdin=None, stdout=None, stderr=None)
    
            command = f"""plotHeatmap --missingDataColor '#ffffff' -m {renormMatFileName} --heatmapHeight {pxH} --interpolationMethod bilinear --refPointLabel '{refPointLabel}' --xAxisLabel '{xAxisLabel}' -out {prefix}.deeptools.renormMatrix.zMax_1.zMin_0.pdf --zMax 1.1 --zMin 0 --{colorsFlag} {colorMaps} {sortUsingSamplesString} {regionsLabelsFlag} {clusterFlag}"""
            lgr.info(f"Plotting heatmap with command: {command}")
    #         subprocess.Popen(command, shell=True, stdout=subprocess.PIPE).stdout.read()
            subprocess.Popen(command, shell=True, stdin=None, stdout=None, stderr=None)

            command = f"""plotProfile -m {renormMatFileName} \
                -out {prefix}.renormMatrix.deeptoolsProfile.pdf \
                --numPlotsPerRow 2 \
                --plotTitle "Profile per sample" \
                --yMax {" ".join(["1.2"] * len(samplesList))} \
                --yMin {" ".join(["-0.2"] * len(samplesList))}"""
            lgr.info(f"Plotting renormalized profile with command: {command}")
            subprocess.Popen(command, shell=True, stdin=None, stdout=None, stderr=None)
    
    #         command = f"""plotHeatmap --missingDataColor '#ffffff' -m {renormMatFileName} --heatmapHeight {pxH} --interpolationMethod bilinear --refPointLabel '{refPointLabel}' --xAxisLabel '{xAxisLabel}' -out {prefix}.deeptoolsHeatmap.renormMatrix.zMax_1.zMin_0.svg --zMax 1.1 --zMin 0 --{colorsFlag} {colorMaps} {sortUsingSamplesString} {regionsLabelsFlag}"""
    #         lgr.info(f"Plotting heatmap with command: {command}")
    # #         subprocess.Popen(command, shell=True, stdout=subprocess.PIPE).stdout.read()
    #         subprocess.Popen(command, shell=True, stdin=None, stdout=None, stderr=None)
        else:
            msg = f"Input Renormalized Matrix file {renormMatFileName} does not exist - plot will not be generated"
            lgr.error(msg)
            raise Exception(msg)

if __name__ == "__main__":
    main()


