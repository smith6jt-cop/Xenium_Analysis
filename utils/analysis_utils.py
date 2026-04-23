"""
Utility functions for Xenium spatial transcriptomics analysis
"""

from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc


def load_xenium_data(
    data_path: Path,
    sample_name: str,
    file_type: str = "h5ad"
) -> sc.AnnData:
    """
    Load Xenium data from various file formats.

    Parameters
    ----------
    data_path : Path
        Path to data directory
    sample_name : str
        Sample identifier
    file_type : str
        File format ('h5ad', 'h5', 'csv')

    Returns
    -------
    AnnData object
    """
    file_path = data_path / f"{sample_name}.{file_type}"

    if file_type == "h5ad":
        adata = sc.read_h5ad(file_path)
    elif file_type == "h5":
        adata = sc.read_10x_h5(file_path)
    elif file_type == "csv":
        adata = sc.read_csv(file_path)
    else:
        raise ValueError(f"Unsupported file type: {file_type}")

    return adata


def add_spatial_coordinates(
    adata: sc.AnnData,
    coords_file: Path | None = None,
    x_col: str = "x",
    y_col: str = "y"
) -> sc.AnnData:
    """
    Add spatial coordinates to AnnData object.

    Parameters
    ----------
    adata : AnnData
        AnnData object
    coords_file : Path, optional
        Path to CSV file with coordinates
    x_col : str
        Column name for x coordinates
    y_col : str
        Column name for y coordinates

    Returns
    -------
    AnnData with spatial coordinates in obsm['spatial']
    """
    if coords_file is not None:
        coords = pd.read_csv(coords_file, index_col=0)
        spatial_coords = coords[[x_col, y_col]].values
    else:
        # Try to get from obs
        if x_col in adata.obs and y_col in adata.obs:
            spatial_coords = adata.obs[[x_col, y_col]].values
        else:
            raise ValueError("Spatial coordinates not found")

    adata.obsm['spatial'] = spatial_coords
    return adata


def calculate_qc_metrics(
    adata: sc.AnnData,
    mt_prefix: str = "MT-",
    ribo_prefix: tuple[str, str] = ("RPS", "RPL")
) -> sc.AnnData:
    """
    Calculate comprehensive QC metrics.

    Parameters
    ----------
    adata : AnnData
        AnnData object
    mt_prefix : str
        Prefix for mitochondrial genes
    ribo_prefix : tuple
        Prefixes for ribosomal genes

    Returns
    -------
    AnnData with QC metrics in obs
    """
    # Identify gene types
    adata.var['mt'] = adata.var_names.str.startswith(mt_prefix)
    adata.var['ribo'] = adata.var_names.str.startswith(ribo_prefix)

    # Calculate metrics
    sc.pp.calculate_qc_metrics(
        adata,
        qc_vars=['mt', 'ribo'],
        percent_top=None,
        log1p=False,
        inplace=True
    )

    return adata


def filter_cells_and_genes(
    adata: sc.AnnData,
    min_genes: int = 10,
    min_cells: int = 3,
    max_mt_pct: float = 20.0,
    min_counts: int = 50
) -> sc.AnnData:
    """
    Filter low-quality cells and genes.

    Parameters
    ----------
    adata : AnnData
        AnnData object
    min_genes : int
        Minimum genes per cell
    min_cells : int
        Minimum cells per gene
    max_mt_pct : float
        Maximum mitochondrial percentage
    min_counts : int
        Minimum UMI counts per cell

    Returns
    -------
    Filtered AnnData
    """
    n_cells_before = adata.n_obs
    n_genes_before = adata.n_vars

    # Filter
    sc.pp.filter_cells(adata, min_genes=min_genes)
    sc.pp.filter_cells(adata, min_counts=min_counts)
    sc.pp.filter_genes(adata, min_cells=min_cells)

    if 'pct_counts_mt' in adata.obs:
        adata = adata[adata.obs.pct_counts_mt < max_mt_pct, :].copy()

    print(f"Filtered: {n_cells_before - adata.n_obs} cells, "
          f"{n_genes_before - adata.n_vars} genes")

    return adata


def normalize_and_hvg(
    adata: sc.AnnData,
    target_sum: float = 1e4,
    n_top_genes: int = 2000,
    flavor: str = 'seurat_v3'
) -> sc.AnnData:
    """
    Normalize data and identify highly variable genes.

    Parameters
    ----------
    adata : AnnData
        AnnData object
    target_sum : float
        Target sum for normalization
    n_top_genes : int
        Number of highly variable genes
    flavor : str
        Method for HVG selection

    Returns
    -------
    Normalized AnnData with HVG annotation
    """
    # Store raw counts
    adata.layers['counts'] = adata.X.copy()

    # Normalize
    sc.pp.normalize_total(adata, target_sum=target_sum)
    sc.pp.log1p(adata)
    adata.layers['log1p_norm'] = adata.X.copy()

    # Identify HVGs
    sc.pp.highly_variable_genes(
        adata,
        n_top_genes=n_top_genes,
        flavor=flavor,
        subset=False
    )

    return adata


def run_standard_workflow(
    adata: sc.AnnData,
    n_pcs: int = 50,
    n_neighbors: int = 15,
    resolution: float = 0.5
) -> sc.AnnData:
    """
    Run standard scanpy workflow: scaling, PCA, neighbors, UMAP, clustering.

    Parameters
    ----------
    adata : AnnData
        Preprocessed AnnData object
    n_pcs : int
        Number of principal components
    n_neighbors : int
        Number of neighbors for graph
    resolution : float
        Resolution for Leiden clustering

    Returns
    -------
    AnnData with embeddings and clusters
    """
    # Scale
    sc.pp.scale(adata, max_value=10)

    # PCA
    sc.tl.pca(adata, svd_solver='arpack', n_comps=n_pcs)

    # Neighbors and UMAP
    sc.pp.neighbors(adata, n_neighbors=n_neighbors, n_pcs=min(30, n_pcs))
    sc.tl.umap(adata)

    # Clustering
    sc.tl.leiden(adata, resolution=resolution)

    return adata


def export_to_csv(
    adata: sc.AnnData,
    output_dir: Path,
    sample_name: str,
    include_spatial: bool = True
) -> None:
    """
    Export AnnData to CSV files.

    Parameters
    ----------
    adata : AnnData
        AnnData object
    output_dir : Path
        Output directory
    sample_name : str
        Sample identifier
    include_spatial : bool
        Whether to include spatial coordinates
    """
    output_dir.mkdir(parents=True, exist_ok=True)

    # Export metadata
    adata.obs.to_csv(output_dir / f"{sample_name}_metadata.csv")

    # Export expression matrix
    expr_df = pd.DataFrame(
        adata.X.toarray() if hasattr(adata.X, 'toarray') else adata.X,
        index=adata.obs_names,
        columns=adata.var_names
    )
    expr_df.to_csv(output_dir / f"{sample_name}_expression.csv")

    # Export spatial coordinates if available
    if include_spatial and 'spatial' in adata.obsm:
        spatial_df = pd.DataFrame(
            adata.obsm['spatial'],
            index=adata.obs_names,
            columns=['x', 'y']
        )
        spatial_df.to_csv(output_dir / f"{sample_name}_spatial_coords.csv")

    print(f"Exported data to {output_dir}")


def create_summary_report(
    adata: sc.AnnData,
    sample_name: str,
    output_file: Path
) -> pd.DataFrame:
    """
    Create a summary report of the dataset.

    Parameters
    ----------
    adata : AnnData
        AnnData object
    sample_name : str
        Sample identifier
    output_file : Path
        Output file path

    Returns
    -------
    DataFrame with summary statistics
    """
    summary = {
        'Sample': sample_name,
        'Total cells': adata.n_obs,
        'Total genes': adata.n_vars,
        'Mean counts per cell': adata.obs['total_counts'].mean() if 'total_counts' in adata.obs else np.nan,
        'Median counts per cell': adata.obs['total_counts'].median() if 'total_counts' in adata.obs else np.nan,
        'Mean genes per cell': adata.obs['n_genes_by_counts'].mean() if 'n_genes_by_counts' in adata.obs else np.nan,
        'Has spatial coords': 'spatial' in adata.obsm,
        'Has cell types': 'celltype' in adata.obs,
    }

    summary_df = pd.DataFrame([summary]).T
    summary_df.columns = ['Value']

    if output_file:
        summary_df.to_csv(output_file)

    return summary_df
