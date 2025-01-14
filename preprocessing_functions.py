from pybiomart import Server

def get_ensemble_labels(adata, organism, server='http://www.ensembl.org',
                        mart='ENSEMBL_MART_ENSEMBL',
                        extended_labeling=False,
                        include_cell_cycle=False):
    """
    Retrieves and labels genes based on various biological categories from Ensembl.

    Parameters:
    - adata (AnnData): Annotated Data object.
    - organism (str): Organism name for which gene labels are retrieved.
    - server (str, optional): Ensembl server URL. Default is 'http://www.ensembl.org'.
    - mart (str, optional): Ensembl Mart name. Default is 'ENSEMBL_MART_ENSEMBL'.
    - extended_labeling (bool, optional): If True, includes additional labeling categories.
    - include_cell_cycle (bool, optional): If True, includes cell cycle genes.

    Returns:
    - adata (AnnData): Annotated Data object with added labels.
    """
    # If want archive ensemble versions, change server:
    # e.g., for 109 enter server = 'http://feb2023.archive.ensembl.org/'
    ensembl_server = Server(host=server)
    dataset = ensembl_server.marts[mart].datasets[organism]

    # Mitochondrial genes
    mito_genes = dataset.query(attributes=['ensembl_gene_id', 'external_gene_name'],
                               filters={'chromosome_name': ['MT']})
    adata.var['mt'] = adata.var['gene_ids'].isin(mito_genes['Gene stable ID'])

    # Ribosome Genes
    # Cytosolic ribosome GO:0022626  ## Includes RPS and RPL proteins
    ribo_genes = dataset.query(attributes=['ensembl_gene_id', 'external_gene_name', 'gene_biotype'],
                               filters={'link_go_closure': ['GO:0022626']})
    adata.var['ribo'] = adata.var.gene_symbols.str.startswith(("RPS", "RPL")) | adata.var['gene_ids'].isin(ribo_genes['Gene stable ID'])


    # Transcription factors: GO:0006355  ## Actually: regulation of DNA-templated transcription
    tf_genes = dataset.query(attributes=['ensembl_gene_id', 'external_gene_name'],
                             filters={'link_go_closure': ['GO:0006355']})
    adata.var['tf'] = adata.var['gene_ids'].isin(tf_genes['Gene stable ID'])

    # Hemoglobin complex
    hb_genes = dataset.query(attributes=['ensembl_gene_id', 'external_gene_name'],
                             filters={'link_go_closure': ['GO:0005833']})
    adata.var['hb'] = adata.var['gene_ids'].isin(hb_genes['Gene stable ID'])

    # Platelet alpha granule
    pl_genes = dataset.query(attributes=['ensembl_gene_id', 'external_gene_name'],
                             filters={'link_go_closure': ['GO:0031091']})
    adata.var['pl'] = adata.var['gene_ids'].isin(pl_genes['Gene stable ID'])

    # Ribosomal RNA
    rRNA_genes = dataset.query(attributes=['ensembl_gene_id', 'external_gene_name'],
                               filters={'biotype': ['rRNA']})
    adata.var['rRNA'] = adata.var['gene_ids'].isin(rRNA_genes['Gene stable ID'])

    # Chromosome Location
    chromosome_loc = dataset.query(attributes=['ensembl_gene_id', 'chromosome_name'])
    chromosome_loc = chromosome_loc.set_index('Gene stable ID').to_dict()['Chromosome/scaffold name']
    adata.var['chromosome_loc'] = adata.var.gene_ids.map(chromosome_loc)

    if extended_labeling:
        # Polycomb group proteins
        pcg_genes = dataset.query(attributes=['ensembl_gene_id', 'external_gene_name'],
                                  filters={'link_go_closure': ['GO:0031519']})
        adata.var['pcg'] = adata.var['gene_ids'].isin(pcg_genes['Gene stable ID'])

        # Extracellular Matrix: GO:0031012
        exc_genes = dataset.query(attributes=['ensembl_gene_id', 'external_gene_name'],
                                  filters={'link_go_closure': ['GO:0031012']})
        adata.var['exc'] = adata.var['gene_ids'].isin(exc_genes['Gene stable ID'])

        # Cytosol: GO:0005829
        cyt_genes = dataset.query(attributes=['ensembl_gene_id', 'external_gene_name'],
                                  filters={'link_go_closure': ['GO:0005829']})
        adata.var['cyt'] = adata.var['gene_ids'].isin(cyt_genes['Gene stable ID'])

        # Membrane: GO:0016020
        mem_genes = dataset.query(attributes=['ensembl_gene_id', 'external_gene_name'],
                                  filters={'link_go_closure': ['GO:0016020']})
        adata.var['mem'] = adata.var['gene_ids'].isin(mem_genes['Gene stable ID'])

    if include_cell_cycle:
        # Cell Cycle genes as defined by cc.genes in Seurat_4.3.0
        hs_s_genes = ["MCM5", "PCNA", "TYMS", "FEN1", "MCM2", "MCM4", "RRM1", "UNG", "GINS2",
                      "MCM6", "CDCA7", "DTL", "PRIM1", "UHRF1", "MLF1IP", "HELLS", "RFC2", "RPA2", "NASP",
                      "RAD51AP1", "GMNN", "WDR76", "SLBP", "CCNE2", "UBR7", "POLD3", "MSH2", "ATAD2", "RAD51",
                      "RRM2", "CDC45", "CDC6", "EXO1", "TIPIN", "DSCC1", "BLM", "CASP8AP2", "USP1", "CLSPN",
                      "POLA1", "CHAF1B", "BRIP1", "E2F8"]
        hs_g2m_genes = ["HMGB2", "CDK1", "NUSAP1", "UBE2C", "BIRC5", "TPX2", "TOP2A", "NDC80", "CKS2",
                        "NUF2", "CKS1B", "MKI67", "TMPO", "CENPF", "TACC3", "FAM64A", "SMC4", "CCNB2", "CKAP2L", "CKAP2",
                        "AURKB", "BUB1", "KIF11", "ANP32E", "TUBB4B", "GTSE1", "KIF20B", "HJURP", "CDCA3", "HN1", "CDC20",
                        "TTK", "CDC25C", "KIF2C", "RANGAP1", "NCAPD2", "DLGAP5", "CDCA2", "CDCA8", "ECT2", "KIF23", "HMMR",
                        "AURKA", "PSRC1", "ANLN", "LBR", "CKAP5", "CENPE", "CTCF", "NEK2", "G2E3", "GAS2L3", "CBX5", "CENPA"]

        hs_homologs = dataset.query(attributes=['ensembl_gene_id', 'external_gene_name',
                                               'hsapiens_homolog_ensembl_gene', 'hsapiens_homolog_associated_gene_name',
                                               'hsapiens_homolog_orthology_type'])

        s_genes = list(hs_homologs[hs_homologs['Human gene name'].isin(hs_s_genes)]['Gene stable ID'])
        g2m_genes = list(hs_homologs[hs_homologs['Human gene name'].isin(hs_g2m_genes)]['Gene stable ID'])

        adata.var['s_genes'] = adata.var_names.isin(adata.var[adata.var['gene_ids'].isin(s_genes)].index)
        adata.var['g2m_genes'] = adata.var_names.isin(adata.var[adata.var['gene_ids'].isin(g2m_genes)].index)

    return adata.copy()


import numpy as np
from scipy.stats import median_abs_deviation


def is_outlier(adata, metric: str, nmads: int):
    """
    Detects outliers in a specified metric column of the AnnData object.

    Parameters:
    - adata (AnnData): Annotated Data object.
    - metric (str): Name of the metric column to analyze for outliers.
    - nmads (int): Number of median absolute deviations (MADs) used to define outliers.

    Returns:
    - outlier (pd.Series): Boolean Series indicating whether each data point is an outlier.
    """
    metric_values = adata.obs[metric]
    median_value = np.median(metric_values)
    mad_value = median_abs_deviation(metric_values)

    outlier = (metric_values < median_value - nmads * mad_value) | (
        median_value + nmads * mad_value < metric_values
    )

    return outlier


def mad_outlier_filtering(adata, filter_out = True, factor = 5,  plot = True, ):
    adata.obs["outlier"] = (
    is_outlier(adata, "log1p_total_counts", factor)
    | is_outlier(adata, "log1p_n_genes_by_counts", factor)
    | is_outlier(adata, "pct_counts_in_top_20_genes", factor)
    )
    if plot:
        sc.pl.scatter(adata, x="pct_counts_in_top_20_genes", y="pct_unspliced", color="outlier", legend_loc='best')
        sc.pl.scatter(adata, x="log1p_n_genes_by_counts", y="pct_counts_in_top_20_genes", legend_loc='best', color="outlier")
        sc.pl.scatter(adata, x="log1p_total_counts", y="pct_counts_in_top_20_genes", legend_loc='best', color="outlier")
        sc.pl.scatter(adata, x = "total_counts", y ="n_genes_by_counts", color="outlier")

        if filter_out:
            adata = adata[(~adata.obs.outlier)].copy()

            print('after filtering')

            sc.pl.scatter(adata, x="pct_counts_in_top_20_genes", y="pct_unspliced", color="outlier", legend_loc='best')
            sc.pl.scatter(adata, x="log1p_n_genes_by_counts", y="pct_counts_in_top_20_genes", legend_loc='best', color="outlier")
            sc.pl.scatter(adata, x="log1p_total_counts", y="pct_counts_in_top_20_genes", legend_loc='best', color="outlier")
            sc.pl.scatter(adata, "total_counts", "n_genes_by_counts", color="outlier")


    return adata



import rpy2.robjects as ro
from rpy2.robjects import pandas2ri

def cook_soup(adata, data, data_tod, genes, cells, soupx_groups,
 r_script_path = "/proj/rnaatlas/nobackup/private/EmilioTemp/macaque_sc/macaque_sc_git/3_preprocessing/make_soup.R"):
    """
    Applies the soupX method to preprocess single-cell RNA-seq data.

    Parameters:
    - adata (AnnData): Annotated Data object.
    - data (type): Description of data parameter.
    - data_tod (type): Description of data_tod parameter.
    - genes (type): Description of genes parameter.
    - cells (type): Description of cells parameter.
    - soupx_groups (type): Description of soupx_groups parameter.
    - r_script_path (str): File path to the R script containing the soupX method.

    Returns:
    - adata (AnnData): Annotated Data object with soupX-processed counts.
    """
    # Your existing Python code
    adata.layers["raw_counts"] = adata.X

    # Convert Python objects to R objects
    ro.globalenv['data'] = data
    ro.globalenv['data_tod'] = data_tod
    ro.globalenv['genes'] = genes
    ro.globalenv['cells'] = cells
    ro.globalenv['soupx_groups'] = soupx_groups

    # Load the R script containing the R function
    ro.r.source(r_script_path)

    # Get the R function into the Python environment
    make_soup = ro.globalenv['make_soup']

    # Call the R function
    out = make_soup(data, data_tod, genes, cells, soupx_groups)

    # Further processing or combining results if needed
    adata.layers["soupX_counts"] = out.T
    adata.X = adata.layers["soupX_counts"]

    return adata

import scanpy as sc
def call_soupx(adata, adata_raw, min_genes=200, plot = False):
    """
    Calls the soupX method to remove ambient RNA from single-cell RNA-seq data.

    Parameters:
    - adata (AnnData): Annotated Data object.
    - adata_raw (AnnData): Annotated Data object containing raw unfiltered reads.
    - min_genes (int, optional): Minimum number of genes for filtering cells. Default is 200.

    Returns:
    - adata (AnnData): Annotated Data object with ambient RNA removed.
    """
    print('╭( •̀ㅂ•́)و Computing ambient RNA...')

    # Filter cells based on the minimum number of genes
    sc.pp.filter_cells(adata, min_genes=min_genes)

    # Create a copy of adata for preprocessing
    adata_pp = adata.copy()

    # Preprocess data
    sc.pp.normalize_total(adata_pp)
    sc.pp.log1p(adata_pp)
    sc.pp.pca(adata_pp)
    sc.pp.neighbors(adata_pp)
    sc.tl.leiden(adata_pp, key_added="soupx_groups")

    # Extract soupx_groups information
    soupx_groups = adata_pp.obs["soupx_groups"]

    # Prepare variables for R magic
    cells = adata.obs_names
    genes = adata.var_names
    data = adata.X.T

    # Load raw unfiltered reads
    data_tod = adata_raw.X.T
    del adata_raw

    # Remove ambient RNA using soupX
    adata = cook_soup(adata, data, data_tod, genes, cells, soupx_groups)

    print('Amount of RNA counts kept:')
    print(adata.layers['soupX_counts'].sum() / adata.layers['raw_counts'].sum())
    print('(⌐■_■) Done!')
    return adata, adata_pp


import pandas as pd
import numpy as np
import scanpy as sc
import anndata

def calculate_pct_unspliced(path_to_velocito, cache=False):
    """
    Calculate the percentage of unspliced reads for each cell.

    Parameters:
    - path_to_velocito (str): Path to the directory containing Velocyto output files.
    - cache (bool): Whether to cache the read matrices. Default is False.

    Returns:
    - pd.DataFrame: A DataFrame containing barcodes and the percentage of unspliced reads.
    """
    mtx_files = ["spliced.mtx.gz", "unspliced.mtx.gz"]
    barcodes_file = 'barcodes.tsv.gz'

    # Read spliced and unspliced matrices
    adata_spliced = sc.read(path_to_velocito + mtx_files[0], cache=cache).T
    adata_unspliced = sc.read(path_to_velocito + mtx_files[1], cache=cache).T

    # Read barcodes
    barcodes = pd.read_csv(path_to_velocito + barcodes_file, sep='\t', header=None, index_col=0)
    barcodes = barcodes.rename_axis('barcodes')

    # Calculate sum of spliced and unspliced counts
    sum_spliced = np.sum(adata_spliced.X, axis=1)
    sum_unspliced = np.sum(adata_unspliced.X, axis=1)

    # Calculate percentage of unspliced reads
    pct_unspliced = 100 * (sum_unspliced / (sum_spliced + sum_unspliced))

    # Add the percentage of unspliced reads to the barcodes DataFrame
    barcodes['pct_unspliced'] = pct_unspliced

    # Delete the read matrices to save memory
    del adata_spliced, adata_unspliced

    return barcodes




import scanpy as sc
import pandas as pd
import numpy as np
import scrublet as scr

def preprocess_adata(star_out_directory: str,
                     mt_limit: float,
                     organism: str,
                     manual_scrublet_threshold: bool = False,
                     scrublet_threshold: float = 0.25,
                     unspliced_limit: float = 0.4,
                     check_unspliced: bool = True,
                     check_ambient: bool = True,
                     var_names: ['gene_symbols', 'gene_ids'] = 'gene_symbols',
                     extended_labeling: bool = False,
                     include_cell_cycle: bool = False,
                     perform_filtering: bool = False,
                     plot: bool = False,
                     cache: bool = False):
    """
    Preprocesses single-cell RNA-seq data in an AnnData object.

    Parameters:
    - star_out_directory (str): Path to the directory containing STAR output files.
    - mt_limit (float): Threshold for mitochondrial genes to define outliers.
    - organism (str): Organism identifier.
    - manual_scrublet_threshold (bool, optional): define if one should have a manual threshold to call doubulets
    - scrublet_threshold (float, optional): Sets the threshold to call doublets, only if manual_scrublet_threshold = True.
    - unspliced_limit (float, optional): Threshold for unspliced transcripts to define outliers. Default is 0.4.
    - check_unspliced (bool, optional): Whether to check and filter barcodes based on unspliced transcripts. Default is True.
    - check_ambient (bool, optional): Whether to check and remove ambient RNA. Default is True.
    - var_names (str or list, optional): The variable names to use for annotation ('gene_symbols' or 'gene_ids'). Default is 'gene_symbols'.
    - extended_labeling (bool, optional): Whether to perform extended gene labeling. Default is False.
    - include_cell_cycle (bool, optional): Whether to include cell cycle information. Default is False.
    - perform_filtering (bool, optional): Whether to perform the outlier and splice filtering, after detection and labeling
    - plot (bool, optional): Whether to generate plots during preprocessing. Default is False.
    - cache (bool, optional): Whether to cache the loaded data for future use. Default is False.

    Returns:
    - adata (AnnData): Preprocessed AnnData object.
    """
    # Read in filtered reads
    adata = sc.read_10x_mtx(star_out_directory + 'Gene/filtered/', var_names=var_names, cache=cache)
    if var_names == 'gene_symbols':
        adata.var['gene_symbols'] = adata.var.index
    elif var_names == 'gene_ids':
        adata.var['gene_ids'] = adata.var.index

    # Remove ambient RNA
    if check_ambient:
        # Read in raw reads
        adata_raw = sc.read_10x_mtx(star_out_directory + 'Gene/raw/', var_names=var_names, cache=cache)
        if plot:
            adata_raw.obs['total_genes'] = adata_raw.X.sum(axis=1)
            df = pd.DataFrame(data={
                'total_genes': adata_raw.obs['total_genes'].sort_values(ascending=False),
                'index': range(len(adata_raw.obs['total_genes'].sort_values(ascending=False)))
            })
            barcode_rankplot(df, 'index', 'total_genes', 'Droplets', 'nUMIs', 'Raw barcode rank plot', log_scale=True, vertical_line=adata.shape[0])

        # Compute and remove ambient RNA
        adata = call_soupx(adata, adata_raw, min_genes=0)

    # Remove barcodes having low amount unspliced transcripts
    if check_unspliced:
        pct_unspliced = calculate_pct_unspliced(star_out_directory + 'Velocyto/filtered/', cache=cache)
        adata.obs = pd.merge(adata.obs, pct_unspliced, left_index=True, right_index=True)
        adata.obs["splice_outlier"] = adata.obs['pct_unspliced'] > unspliced_limit

    # Gene labeling and calculate general statistics on barcodes
    adata = get_ensemble_labels(adata, organism=organism + '_gene_ensembl',
                                server='http://feb2023.archive.ensembl.org/',
                                mart='ENSEMBL_MART_ENSEMBL',
                                extended_labeling=extended_labeling,
                                include_cell_cycle=include_cell_cycle)
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt", "ribo", "hb", "pl"],
                               inplace=True, percent_top=[20], log1p=True)

    # Define outlier barcodes to filter out based on mitochondrial and unspliced limit
    adata.obs["mt_outlier"] = adata.obs["pct_counts_mt"] > mt_limit
    adata.obs['outlier'] = adata.obs['mt_outlier'] | (check_unspliced & adata.obs['splice_outlier'])

    # Plot overview before filtering
    if plot:
        sc.pl.scatter(adata, "total_counts", "n_genes_by_counts", color="outlier")
        sc.pl.violin(adata, keys=['pct_counts_mt', 'pct_counts_ribo', 'pct_counts_hb', 'pct_counts_pl'], multi_panel=True)
        if check_unspliced:
            sc.pl.scatter(adata, x="pct_counts_mt", y="pct_unspliced", color="outlier", legend_loc='best')
            sc.pl.scatter(adata, x="pct_counts_mt", y="pct_counts_in_top_20_genes", legend_loc='best', color="outlier")
            sc.pl.scatter(adata, x="pct_unspliced", y="pct_counts_in_top_20_genes", legend_loc='best', color="outlier")

        else:
            sc.pl.scatter(adata, x="n_genes_by_counts", y="pct_counts_mt", color="outlier", ax=axes[1], color_map=None, show=False)

    if perform_filtering:

        # Perform filtering
        print(f"Total number of cells: {adata.n_obs}")
        adata = adata[(~adata.obs.outlier)].copy()
        print(f"Number of cells after filtering of low-quality cells: {adata.n_obs}")

        sc.pp.calculate_qc_metrics(adata, qc_vars=["mt", "ribo", "hb", "pl"],
                                   inplace=True, percent_top=[20], log1p=True)

        # Plot again after filtering
        if plot:
            sc.pl.scatter(adata, "total_counts", "n_genes_by_counts", color="pct_counts_mt")
            sc.pl.violin(adata, keys=['pct_counts_mt', 'pct_counts_ribo', 'pct_counts_hb', 'pct_counts_pl'], multi_panel=True)

            if check_unspliced:
                sc.pl.scatter(adata, "total_counts", "n_genes_by_counts", color="pct_unspliced")
                sc.pl.scatter(adata, x="pct_counts_mt", y="pct_unspliced", color="pct_counts_in_top_20_genes", legend_loc='best')
                sc.pl.scatter(adata, x="pct_counts_mt", y="pct_counts_in_top_20_genes", color="pct_unspliced", legend_loc='best')
                sc.pl.scatter(adata, x="pct_counts_mt", y="log1p_n_genes_by_counts", color="pct_unspliced", legend_loc='best')


            else:
                sc.pl.scatter(adata, x="pct_counts_mt", y="pct_counts_in_top_20_genes", color="log1p_n_genes_by_counts", legend_loc='best')
                sc.pl.scatter(adata, x="n_genes_by_counts", y="pct_counts_mt", color="outlier", ax=axes[1], color_map=None, show=False)

    # Perform doublet check
    if check_ambient:
        adata.X = adata.layers['raw_counts'].copy()
    scrub = scr.Scrublet(adata.X, expected_doublet_rate=0.05, random_state=4)
    adata.obs['doublet_scores'], adata.obs['predicted_doublets'] = scrub.scrub_doublets()

    if manual_scrublet_threshold:
        adata.obs['predicted_doublets'] = scrub.call_doublets(threshold=scrublet_threshold)
    scrub.plot_histogram()

    adata.obs['doublet_info'] = adata.obs["predicted_doublets"].astype(str)

    if check_ambient:
        adata.X = adata.layers['soupX_counts'].copy()

    if plot:
        p6 = sc.pl.scatter(adata, "total_counts", "pct_counts_in_top_20_genes", color="doublet_info")

    if perform_filtering:
        adata = adata[adata.obs.predicted_doublets == False].copy()

        sc.pp.calculate_qc_metrics(adata, qc_vars=["mt", "ribo", "hb", "pl"],
        inplace=True, percent_top=[20], log1p=True)

    return adata


import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd

def plot_line(df, x_col, y_col, x_label, y_label, title, log_scale=True, vertical_line=None):
    """
    Plot a line graph using the provided DataFrame.

    Parameters:
    - df (pd.DataFrame): DataFrame containing the data.
    - x_col (str): Column name for the x-axis.
    - y_col (str): Column name for the y-axis.
    - x_label (str): Label for the x-axis.
    - y_label (str): Label for the y-axis.
    - title (str): Title of the plot.
    - log_scale (bool, optional): Whether to use log scale for both axes. Default is True.
    - vertical_line (float, optional): X-axis position for a vertical line. Default is None.
    """
    # Set up the figure and axes
    plt.figure(figsize=(6, 4))
    ax = sns.lineplot(x=x_col, y=y_col, data=df)

    # Set labels and title
    ax.set(xlabel=x_label, ylabel=y_label)
    ax.set_title(title)

    # Set log scale if specified
    if log_scale:
        ax.set_yscale('log')
        ax.set_xscale('log')

    # Draw a vertical line at the specified position
    if vertical_line is not None:
        ax.axvline(x=vertical_line, color='red', linestyle='--', label=f'Filtered cutoff at {vertical_line}th cell')

    # Show the legend
    ax.legend()

    # Show the plot
    plt.show()
    return ax
def barcode_rankplot(df, x_col, y_col, x_label, y_label, title, log_scale=True, vertical_line=None):
    """
    Plot a line graph using the provided DataFrame.

    Parameters:
    - df (pd.DataFrame): DataFrame containing the data.
    - x_col (str): Column name for the x-axis.
    - y_col (str): Column name for the y-axis.
    - x_label (str): Label for the x-axis.
    - y_label (str): Label for the y-axis.
    - title (str): Title of the plot.
    - log_scale (bool, optional): Whether to use log scale for both axes. Default is True.
    - vertical_line (float, optional): X-axis position for a vertical line. Default is None.

    Returns:
    - fig (matplotlib.figure.Figure): The Figure object containing the plot.
    """
    # Create a new figure
    fig, ax = plt.subplots(figsize=(6, 4))
    sns.lineplot(x=x_col, y=y_col, data=df, ax=ax)

    # Set labels and title
    ax.set(xlabel=x_label, ylabel=y_label)
    ax.set_title(title)

    # Set log scale if specified
    if log_scale:
        ax.set_yscale('log')
        ax.set_xscale('log')

    # Draw a vertical line at the specified position
    if vertical_line is not None:
        ax.axvline(x=vertical_line, color='red', linestyle='--', label=f'Filtered cutoff at {vertical_line}th cell')

    # Show the legend
    ax.legend()
    return fig


import scanpy as sc
import pandas as pd
import anndata

def read_velocyto(path, var_names: ['gene_symbols', 'gene_ids'] = 'gene_symbols', incl_ambiguous = True, make_unique: bool = True, cache: bool = False):
    """
    Read Velocyto output files into an AnnData object.

    Parameters:
    ----------
    path : str
        The path to the directory containing Velocyto output files.
    var_names : {'gene_symbols', 'gene_ids'}, optional, default: 'gene_symbols'
        The variable names to use for annotation (gene_symbols or gene_ids).
    make_unique : bool, optional, default: True
        Whether to make the variable names unique.
    cache : bool, optional, default: False
        Whether to cache the loaded data for future use.

    Returns:
    -------
    adata : AnnData
        An AnnData object containing spliced, unspliced, and ambiguous data.

    Notes:
    ------
    This function assumes that the Velocyto output files have the following naming convention:
    - "ambiguous.mtx.gz"
    - "spliced.mtx.gz"
    - "unspliced.mtx.gz"
    - "barcodes.tsv.gz"
    - "features.tsv.gz"
    """
    if incl_ambiguous:

        mtx_files = ["ambiguous.mtx.gz", "spliced.mtx.gz", "unspliced.mtx.gz"]
        barcodes_file = 'barcodes.tsv.gz'
        features_file = 'features.tsv.gz'

        # Read spliced, unspliced, and ambiguous matrices
        adata_spliced = sc.read(path + mtx_files[1], cache=cache)
        adata_unspliced = sc.read(path + mtx_files[2], cache=cache)
        adata_ambiguous = sc.read(path + mtx_files[0], cache=cache)

        # Read features file based on variable names
        if var_names == 'gene_symbols':
            features = pd.read_csv(path + features_file, sep='\t', header=None, index_col=1)
            features = features.rename_axis('gene_symbols')
            features.columns = ['gene_ids', 'feature_type']
            features['gene_symbols'] = features.index

        elif var_names == 'gene_ids':
            features = pd.read_csv(path + features_file, sep='\t', header=None, index_col=0)
            features = features.rename_axis('gene_ids')
            features.columns = ['gene_symbols', 'feature_type']
            features['gene_ids'] = features.index


        # Read barcodes
        barcodes = pd.read_csv(path + barcodes_file, sep='\t', header=None, index_col=0)
        barcodes = barcodes.rename_axis('barcodes')

        # Create AnnData object
        adata = anndata.AnnData(
            X=adata_spliced.X.T,
            obs=barcodes,
            var=features
        )

        # Make variable names unique if specified
        if make_unique:
            adata.var_names_make_unique()

        # Add layers for unspliced and ambiguous data
        adata.layers['spliced'] = adata.X
        adata.layers['unspliced'] = adata_unspliced.X.T
        adata.layers['ambiguous'] = adata_ambiguous.X.T
        adata.layers['velocyto'] = adata.layers['spliced'] + adata.layers['unspliced'] +  adata.layers['ambiguous']
        adata.X = adata.layers['velocyto'].copy()

        del adata_unspliced, adata_ambiguous, adata_spliced

        adata.obs['sum_spliced'] = np.sum(adata.layers['spliced'], axis = 1)
        adata.obs['sum_unspliced'] = np.sum(adata.layers['unspliced'], axis = 1)
        adata.obs['sum_ambiguous'] = np.sum(adata.layers['ambiguous'], axis = 1)
        adata.obs['pct_unspliced'] = 100*(adata.obs['sum_unspliced']/(adata.obs['sum_spliced'] + adata.obs['sum_unspliced']))
    else:
        mtx_files = [ "spliced.mtx.gz", "unspliced.mtx.gz"]
        barcodes_file = 'barcodes.tsv.gz'
        features_file = 'features.tsv.gz'

        # Read spliced, unspliced, and ambiguous matrices
        adata_spliced = sc.read(path + mtx_files[0], cache=cache)
        adata_unspliced = sc.read(path + mtx_files[1], cache=cache)


        # Read features file based on variable names
        if var_names == 'gene_symbols':
            features = pd.read_csv(path + features_file, sep='\t', header=None, index_col=1)
            features = features.rename_axis('gene_symbols')
            features.columns = ['gene_ids', 'feature_type']
            features['gene_symbols'] = features.index

        elif var_names == 'gene_ids':
            features = pd.read_csv(path + features_file, sep='\t', header=None, index_col=0)
            features = features.rename_axis('gene_ids')
            features.columns = ['gene_symbols', 'feature_type']
            features['gene_ids'] = features.index


        # Read barcodes
        barcodes = pd.read_csv(path + barcodes_file, sep='\t', header=None, index_col=0)
        barcodes = barcodes.rename_axis('barcodes')

        # Create AnnData object
        adata = anndata.AnnData(
            X=adata_spliced.X.T,
            obs=barcodes,
            var=features
        )

        # Make variable names unique if specified
        if make_unique:
            adata.var_names_make_unique()

        # Add layers for unspliced and ambiguous data
        adata.layers['spliced'] = adata.X
        adata.layers['unspliced'] = adata_unspliced.X.T
        adata.layers['velocyto'] = adata.layers['spliced'] + adata.layers['unspliced']
        adata.X = adata.layers['velocyto'].copy()

        del adata_unspliced, adata_spliced

        adata.obs['sum_spliced'] = np.sum(adata.layers['spliced'], axis = 1)
        adata.obs['sum_unspliced'] = np.sum(adata.layers['unspliced'], axis = 1)
        adata.obs['pct_unspliced'] = 100*(adata.obs['sum_unspliced']/(adata.obs['sum_spliced'] + adata.obs['sum_unspliced']))

    return adata

from pathlib import Path, PurePath
from typing import BinaryIO, Literal
from scanpy._utils import Empty, _empty
from scanpy import read

from anndata import (
    AnnData,
    read_csv,
    read_excel,
    read_h5ad,
    read_hdf,
    read_loom,
    read_mtx,
    read_text,
)
def read_STARsolo(
    path: Path | str,
    *,
    var_names: Literal["gene_symbols", "gene_ids"] = "gene_symbols",
    matrix_name: str = 'matrix',
    make_unique: bool = True,
    cache: bool = False,
    cache_compression: Literal["gzip", "lzf"] | None | Empty = _empty,
    gex_only: bool = True,
    prefix: str | None = None,
) -> AnnData:
    """\
    Read 10x-Genomics-formatted mtx directory.

    Parameters
    ----------
    path
        Path to directory for `.mtx` and `.tsv` files,
        e.g. './filtered_gene_bc_matrices/hg19/'.
    var_names
        The variables index.
    make_unique
        Whether to make the variables index unique by appending '-1',
        '-2' etc. or not.
    cache
        If `False`, read from source, if `True`, read from fast 'h5ad' cache.
    cache_compression
        See the h5py :ref:`dataset_compression`.
        (Default: `settings.cache_compression`)
    gex_only
        Only keep 'Gene Expression' data and ignore other feature types,
        e.g. 'Antibody Capture', 'CRISPR Guide Capture', or 'Custom'
    prefix
        Any prefix before `matrix.mtx`, `genes.tsv` and `barcodes.tsv`. For instance,
        if the files are named `patientA_matrix.mtx`, `patientA_genes.tsv` and
        `patientA_barcodes.tsv` the prefix is `patientA_`.
        (Default: no prefix)

    Returns
    -------
    An :class:`~anndata.AnnData` object
    """
    path = Path(path)
    prefix = "" if prefix is None else prefix
    is_legacy = (path / f"{prefix}genes.tsv").is_file()
    adata = _read_STARsolo(
        path,
        var_names=var_names,
        matrix_name=matrix_name,
        make_unique=make_unique,
        cache=cache,
        cache_compression=cache_compression,
        prefix=prefix,
        is_legacy=is_legacy,
    )
    if is_legacy or not gex_only:
        return adata
    gex_rows = adata.var["feature_types"] == "Gene Expression"
    return adata[:, gex_rows].copy()


def _read_STARsolo(
    path: Path,
    *,
    var_names: Literal["gene_symbols", "gene_ids"] = "gene_symbols",
    matrix_name: str = 'matrix',
    make_unique: bool = True,
    cache: bool = False,
    cache_compression: Literal["gzip", "lzf"] | None | Empty = _empty,
    prefix: str = "",
    is_legacy: bool,
) -> AnnData:
    """
    Read mex from output from Cell Ranger v2- or v3+
    """
    suffix = "" if is_legacy else ".gz"
    adata = read(
        path / f"{prefix}{matrix_name}.mtx{suffix}",
        cache=cache,
        cache_compression=cache_compression,
    ).T  # transpose the data
    genes = pd.read_csv(
        path / f"{prefix}{'genes' if is_legacy else 'features'}.tsv{suffix}",
        header=None,
        sep="\t",
    )
    if var_names == "gene_symbols":
        var_names_idx = pd.Index(genes[1].values)
        if make_unique:
            var_names_idx = anndata.utils.make_index_unique(var_names_idx)
        adata.var_names = var_names_idx
        adata.var["gene_ids"] = genes[0].values
    elif var_names == "gene_ids":
        adata.var_names = genes[0].values
        adata.var["gene_symbols"] = genes[1].values
    else:
        raise ValueError("`var_names` needs to be 'gene_symbols' or 'gene_ids'")
    if not is_legacy:
        adata.var["feature_types"] = genes[2].values
    barcodes = pd.read_csv(path / f"{prefix}barcodes.tsv{suffix}", header=None)
    adata.obs_names = barcodes[0].values
    return adata



import os
def listdir(path):
    return [file for file in os.listdir(path) if not file.startswith('.DS_Store')]



#Scrublet functions
import scrublet as scr
def call_scrublet(adata, expected_doublet_rate = 0.05,
                  min_counts=3, min_cells=3, min_gene_variability_pctl=85,
                  n_prin_comps=30, plot_histogram = True, plot_umap = True, doublet_threshold = None, plot_title = None
                 ):
    scrub = scr.Scrublet(adata.X, expected_doublet_rate=expected_doublet_rate)
    adata.obs['doublet_score'], adata.obs['predicted_doublet'] = scrub.scrub_doublets(min_counts=min_counts, min_cells=min_cells,
                                                          min_gene_variability_pctl=min_gene_variability_pctl, n_prin_comps=n_prin_comps)
    adata.uns['scrublet'] = scrub

    if doublet_threshold != None:
        adata = set_doublet_threshold(adata, doublet_threshold, plot_histogram, plot_umap)

    _scrub_plots(adata, plot_histogram, plot_umap, plot_title)

    return adata


def set_doublet_threshold(adata, doublet_threshold, plot_histogram = True, plot_umap = True, plot_title = None):
    adata.uns['scrublet'].call_doublets(threshold=doublet_threshold)
    adata.obs['doublet_score'] = adata.uns['scrublet'].doublet_scores_obs_
    adata.obs['predicted_doublet'] = adata.uns['scrublet'].predicted_doublets_
    _scrub_plots(adata, plot_histogram, plot_umap, plot_title = plot_title)
    return adata


def _scrub_plots(adata, plot_histogram, plot_umap, plot_title = None):
    if plot_histogram:
        fig, axs = adata.uns['scrublet'].plot_histogram()
        if plot_title != None:
            fig.suptitle(plot_title, weight = 'bold', y = 1.02)
            plt.show()

    if plot_umap:
        if 'UMAP' not in adata.uns['scrublet']._embeddings.keys():
            adata.uns['scrublet'].set_embedding('UMAP', scr.get_umap(adata.uns['scrublet'].manifold_obs_, 10, min_dist=0.3))
        fig, axs = adata.uns['scrublet'].plot_embedding('UMAP', order_points=True)
        if plot_title != None:
            fig.suptitle(plot_title, weight = 'bold', y = 1.02)
            plt.show()


#export functions
from contextlib import contextmanager
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt

@contextmanager
def save_plots_to_pdf(pdf_filename):
    pdf = PdfPages(pdf_filename)
    original_show = plt.show

    def save_figure(*args, **kwargs):
        fig = plt.gcf()
        pdf.savefig(fig, bbox_inches='tight')
        original_show(*args, **kwargs)

    plt.show = save_figure
    try:
        yield
    finally:
        plt.show = original_show
        pdf.close()


from pybiomart import Server

def label_hs_orthologs(adata, organism, server, mart):
    attributes = ['ensembl_gene_id',  # Macaque gene ID
                  'external_gene_name',  # Macaque gene name
                  'hsapiens_homolog_ensembl_gene',  # Human ortholog gene ID
                  'hsapiens_homolog_associated_gene_name']  # Human ortholog gene name
    ensembl_server = Server(host=server)
    dataset = ensembl_server.marts[mart].datasets[organism + '_gene_ensembl']

    response = dataset.query(attributes = attributes)

    grouped_df = response.groupby('Gene stable ID')['Human gene stable ID'].agg(lambda x: ', '.join(x.astype(str))).reset_index()

    adata.var['hsapiens_orthologs_id'] = adata.var.index.map(grouped_df.set_index('Gene stable ID').to_dict()['Human gene stable ID'])


    grouped_df = response.groupby('Gene stable ID')['Human gene name'].agg(lambda x: ', '.join(x.astype(str))).reset_index()

    adata.var['hsapiens_orthologs_name'] = adata.var.index.map(grouped_df.set_index('Gene stable ID').to_dict()['Human gene name'])

    return(adata)

from pybiomart import Server
def get_orthologs_df(organism, mart, server):
    """
    e.g.
    mart = 'ENSEMBL_MART_ENSEMBL'
    server = 'http://feb2023.archive.ensembl.org/'
    organism = 'mfascicularis'

    """
    attributes = ['ensembl_gene_id',  # Macaque gene ID
              'external_gene_name',  # Macaque gene name
              'hsapiens_homolog_ensembl_gene',  # Human ortholog gene ID
              'hsapiens_homolog_associated_gene_name']  # Human ortholog gene name
    ensembl_server = Server(host=server)
    dataset = ensembl_server.marts[mart].datasets[organism + '_gene_ensembl']

    response = dataset.query(attributes = attributes)
    response = response.dropna()
    response = response.reset_index(drop = True)

    return response
