import os

import numpy as np
import pandas as pd


from anndata import AnnData
import scanpy as sc
from scanpy.plotting import palettes
palette = palettes.default_102

from math import ceil

from matplotlib import gridspec
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.colors import to_rgba
from matplotlib.colors import ListedColormap
from mpl_toolkits.axes_grid1 import make_axes_locatable


from scipy.sparse import issparse

import bin2cell as b2c

from PIL import Image
# setting needed so PIL can load the large TIFFs
Image.MAX_IMAGE_PIXELS = None

# setting needed so cv2 can load the large TIFFs
os.environ["OPENCV_IO_MAX_IMAGE_PIXELS"] = str(2**40)
import cv2


from typing import List, Dict, Optional, Union, Tuple



def compute_spot_density(
    adata: AnnData, 
    row_key: str ="array_row", 
    col_key: str ="array_col", 
    n_row_bins: int =200, 
    n_col_bins: int =200
) -> AnnData:
    """
    Compute a lightweight 2D density map of spatial transcriptomics spots 
    and store it in `adata.obs['spot_density']`.

    Parameters
    ----------
    adata : AnnData
        Annotated data object with spot coordinates in `adata.obs`.
    row_key : str, default="array_row"
        Column in `adata.obs` for row coordinates.
    col_key : str, default="array_col"
        Column in `adata.obs` for column coordinates.
    n_row_bins : int, default=200
        Number of bins along the row axis.
    n_col_bins : int, default=200
        Number of bins along the column axis.

    Returns
    -------
    adata : AnnData
        Modified AnnData with `adata.obs['spot_density']`.
    """

    # Define binning
    row_bins = np.linspace(adata.obs[row_key].min(), adata.obs[row_key].max(), n_row_bins + 1)
    col_bins = np.linspace(adata.obs[col_key].min(), adata.obs[col_key].max(), n_col_bins + 1)

    # Digitize spots into bins
    row_idx = np.digitize(adata.obs[row_key], bins=row_bins) - 1
    col_idx = np.digitize(adata.obs[col_key], bins=col_bins) - 1

    row_idx = np.clip(row_idx, 0, n_row_bins - 1)
    col_idx = np.clip(col_idx, 0, n_col_bins - 1)

    # Compute density per bin
    density_map = np.zeros((n_row_bins, n_col_bins), dtype=int)
    for r, c in zip(row_idx, col_idx):
        density_map[r, c] += 1

    # Map density back to each spot
    spot_density = density_map[row_idx, col_idx]

    # Store in AnnData
    adata.obs["spot_density"] = spot_density

    return adata


def plot_image_options(
    adata: AnnData,
    library_id: str,
    plot_dict: Dict[str, str],
    plot_axes: bool = False,
    image_orientation: Optional[List[str]] = None,
    scatter_orientation: Optional[List[str]] = None,
    shared_orientation: Optional[List[str]] = None,
    figsize: Optional[Tuple[int, int]] = (18, 6)
) -> None:
    """
    Plot H&E images from AnnData using the same affine transforms as preview_spatial_zoom_region.
    """
    if image_orientation is None:
        image_orientation = []
    if scatter_orientation is None:
        scatter_orientation = []
    if shared_orientation is None:
        shared_orientation = []

    for img_key, basis_name in plot_dict.items():
        spatial_data = adata.uns['spatial'][library_id]
        img = spatial_data['images'][img_key].copy()
        scale_factor = spatial_data['scalefactors'][f'tissue_{img_key}_scalef']

        coords = adata.obsm[basis_name]
        x = coords[:, 0] * scale_factor
        y = coords[:, 1] * scale_factor

        # Crop to scatter bounds
        x_min, x_max = x.min(), x.max()
        y_min, y_max = y.min(), y.max()
        img = img[int(y_min):int(y_max), int(x_min):int(x_max)].copy()
        x -= x_min
        y -= y_min

        # -----------------------------
        # Image-only transforms
        # -----------------------------
        img, _, _, Hf, Wf = apply_affine_to_image_and_coords(
            img=img, x=None, y=None, transforms=image_orientation
        )

        # -----------------------------
        # Scatter-only transforms
        # -----------------------------
        _, x, y, Hf, Wf = apply_affine_to_image_and_coords(
            img=None, x=x, y=y, transforms=scatter_orientation
        )

        # -----------------------------
        # Shared transforms (both)
        # -----------------------------
        img_final, x_final, y_final, Hf, Wf = apply_affine_to_image_and_coords(
            img=img, x=x, y=y, transforms=shared_orientation
        )

        # Determine plotting limits
        xmin, xmax = x_final.min(), x_final.max()
        ymin, ymax = y_final.min(), y_final.max()
        xmin = max(0.0, xmin)
        xmax = min(float(Wf - 1), xmax)
        ymin = max(0.0, ymin)
        ymax = min(float(Hf - 1), ymax)

        # Plot
        fig, ax = plt.subplots(1, 1, figsize=figsize)
        ax.imshow(img_final, origin='lower')
        ax.set_title(img_key)
        ax.set_xlim(xmin, xmax)
        ax.set_ylim(ymin, ymax)  # invert y to match plotting orientation
        if not plot_axes:
            ax.axis('off')
        plt.tight_layout()
        plt.show()



def make_rectangle(x, y, width, height):
    return patches.Rectangle((x, y), width, height, linewidth=1, edgecolor="red", facecolor="none")

# -----------------------------
# Build affine matrix for transforms
# -----------------------------
def build_affine_matrix(h, w, transforms: List[str]):
    """
    Build a single 3x3 affine matrix for a list of rotate/flip transforms.
    Returns: M (3x3), new_h, new_w
    """
    M_total = np.eye(3, dtype=np.float32)
    for t in transforms:
        M = np.eye(3, dtype=np.float32)
        if t == "flip_h":
            M = np.array([[-1, 0, w-1],
                          [0, 1, 0],
                          [0, 0, 1]], dtype=np.float32)
        elif t == "flip_v":
            M = np.array([[1, 0, 0],
                          [0, -1, h-1],
                          [0, 0, 1]], dtype=np.float32)
        elif t == "rotate_cw":
            M = np.array([[0, 1, 0],
                          [-1, 0, w-1],
                          [0, 0, 1]], dtype=np.float32)
            h, w = w, h
        elif t == "rotate_ccw":
            M = np.array([[0, -1, h-1],
                          [1, 0, 0],
                          [0, 0, 1]], dtype=np.float32)
            h, w = w, h
        M_total = M @ M_total
    return M_total, h, w

# -----------------------------
# Apply affine matrix to image and/or scatter coordinates
# -----------------------------
def apply_affine_to_image_and_coords(img=None, x=None, y=None, transforms=None):
    """
    Apply a list of transforms to image and/or scatter coordinates independently.
    - img: np.ndarray or None
    - x, y: scatter coordinates or None
    Returns transformed img, x, y, new_h, new_w
    """
    if transforms is None or len(transforms) == 0:
        h = img.shape[0] if img is not None else (0 if x is None else max(y)+1)
        w = img.shape[1] if img is not None else (0 if x is None else max(x)+1)
        return img, x, y, h, w

    # Determine starting shape
    h = img.shape[0] if img is not None else (0 if x is None else max(y)+1)
    w = img.shape[1] if img is not None else (0 if x is None else max(x)+1)

    # Build affine matrix
    M, h_new, w_new = build_affine_matrix(h, w, transforms)

    # Transform image only
    if img is not None:
        img_t = cv2.warpAffine(img, M[:2, :], (w_new, h_new), flags=cv2.INTER_NEAREST)
    else:
        img_t = None

    # Transform scatter only
    if x is not None and y is not None:
        coords_h = np.stack([x, y, np.ones_like(x)], axis=0)
        coords_trans = M @ coords_h
        x_t, y_t = coords_trans[0], coords_trans[1]
    else:
        x_t, y_t = x, y

    return img_t, x_t, y_t, h_new, w_new


def preview_spatial_zoom_region(
    adata: AnnData,
    library_id: str,
    basis_name: str = "spatial",
    img_key: str = "hires",
    plot_scatter: bool = True,
    scatter_size: float = 0.1,
    scatter_color: str = "blue",
    scatter_alpha: float = 0.02,
    plot_rectangle: bool = True,
    row_start: int = None,
    row_end: int = None,
    col_start: int = None,
    col_end: int = None,
    plot_axes: bool = True,
    image_orientation: Optional[List[str]] = None,
    scatter_orientation: Optional[List[str]] = None,
    shared_orientation: Optional[List[str]] = None,
    observe_selected_area = False,
    selected_area_plot_dict = {'hires' : 'spatial'},
):
    if image_orientation is None:
        image_orientation = []
    if scatter_orientation is None:
        scatter_orientation = []
    if shared_orientation is None:
        shared_orientation = []

    # -----------------------------
    # Load image and coordinates
    # -----------------------------
    img = adata.uns["spatial"][library_id]["images"][img_key].copy()
    scale_factor = adata.uns["spatial"][library_id]["scalefactors"][f"tissue_{img_key}_scalef"]
    coords = adata.obsm[basis_name] * scale_factor
    x_coords = coords[:, 0]
    y_coords = coords[:, 1]

    # Crop image to scatter bounds
    x_min, x_max = int(x_coords.min()), int(x_coords.max())
    y_min, y_max = int(y_coords.min()), int(y_coords.max())
    img = img[y_min:y_max, x_min:x_max].copy()

    # Translate scatter to crop coordinates
    x = x_coords - x_min
    y = y_coords - y_min

    # -----------------------------
    # Image-only transforms
    # -----------------------------
    img, _, _, h, w = apply_affine_to_image_and_coords(img=img, x=None, y=None, transforms=image_orientation)

    # -----------------------------
    # Scatter-only transforms
    # -----------------------------
    _, x, y, h, w = apply_affine_to_image_and_coords(img=None, x=x, y=y, transforms=scatter_orientation)

    # -----------------------------
    # Shared transforms applied to both
    # -----------------------------
    img, x, y, h, w = apply_affine_to_image_and_coords(img=img, x=x, y=y, transforms=shared_orientation)

    # -----------------------------
    # Rectangle in final coordinates
    # -----------------------------
    rect = None
    if plot_rectangle and None not in (row_start, row_end, col_start, col_end):
        x0 = col_start
        y0 = row_start
        width = col_end - col_start
        height = row_end - row_start
        rect = make_rectangle(x0, y0, width, height)

    # -----------------------------
    # Plot
    # -----------------------------
    fig, axes = plt.subplots(1, 2, figsize=(16, 8))
    axes[0].imshow(img, origin="lower")
    axes[0].set_title("H&E")
    if rect:
        axes[0].add_patch(rect)
    if not plot_axes:
        axes[0].axis("off")

    axes[1].imshow(img, origin="lower")
    if plot_scatter:
        axes[1].scatter(x, y, s=scatter_size, c=scatter_color, alpha=scatter_alpha,
                        edgecolors="none", rasterized=True)
    if rect:
        axes[1].add_patch(make_rectangle(rect.get_x(), rect.get_y(), rect.get_width(), rect.get_height()))
    axes[1].set_title("H&E with scatter")
    if not plot_axes:
        axes[1].axis("off")

    plt.show()
        
    
    
    
    if observe_selected_area:
        missing = [name for name, val in {"row_start": row_start, "row_end": row_end, "col_start": col_start, "col_end": col_end}.items() if val is None]
        if missing: raise ValueError(f"Missing required values: {', '.join(missing)}")
            
        x, y = get_transformed_coords(adata, library_id,scatter_orientation=scatter_orientation, shared_orientation=shared_orientation)
        mask = (y >= row_start) & (y <= row_end) & (x >= col_start) & (x <= col_end)
        sub = adata[mask]
        
        plot_spatial_image(
            sub,
            library_id,
            plot_dict = selected_area_plot_dict,
            image_orientation = image_orientation,
            scatter_orientation = scatter_orientation,
            shared_orientation = shared_orientation,
            )


def get_transformed_coords(adata, library_id, basis='spatial',
                           scatter_orientation=None, shared_orientation=None):
    coords = adata.obsm[basis].copy() * adata.uns["spatial"][library_id]["scalefactors"]["tissue_hires_scalef"]
    x, y = coords[:,0], coords[:,1]
    x_min, y_min = x.min(), y.min()
    x -= x_min
    y -= y_min
    _, x, y, _, _ = apply_affine_to_image_and_coords(img=None, x=x, y=y, transforms=scatter_orientation or [])
    _, x, y, _, _ = apply_affine_to_image_and_coords(img=None, x=x, y=y, transforms=shared_orientation or [])
    return x, y


def get_values(adata, color_col: str, use_raw: bool = False) -> np.ndarray:
    if color_col in adata.obs.columns:
        return adata.obs_vector(color_col)
    elif color_col in adata.var_names:
        values = adata[:, color_col].X
        if hasattr(values, "todense"):
            return np.array(values.todense()).flatten()
        else:
            return np.array(values).flatten()
    elif use_raw and adata.raw is not None and color_col in adata.raw.var_names:
        values = adata.raw[:, color_col].X
        if hasattr(values, "todense"):
            return np.array(values.todense()).flatten()
        else:
            return np.array(values).flatten()
    else:
        raise ValueError(f"'{color_col}' not found in obs, var_names, or raw.var_names")

def plot_spatial_image(
    adata: AnnData,
    library_id: str,
    plot_dict: Dict[str, str],
    color: List[str or float] = None,
    groups=None,
    size: int =10,
    plot_axes: bool = False,
    image_orientation: Optional[List[str]] = None,
    scatter_orientation: Optional[List[str]] = None,
    shared_orientation: Optional[List[str]] = None,
    figsize: Optional[Tuple[int, int]] = (10, 10),
    use_raw: bool = False,
    ncols: int = 2,
    cmap: str = "OrRd",
    category_colors: Optional[Dict[str, str]] = None,
    show_legend: bool = True,
    save: str = False,
) -> None:

    if image_orientation is None:
        image_orientation = []
    if scatter_orientation is None:
        scatter_orientation = []
    if shared_orientation is None:
        shared_orientation = []
        
    

    for img_key, basis_name in plot_dict.items():
        spatial_data = adata.uns['spatial'][library_id]
        img = spatial_data['images'][img_key].copy()
        scale_factor = spatial_data['scalefactors'][f'tissue_{img_key}_scalef']

        coords = adata.obsm[basis_name]
        x = coords[:, 0] * scale_factor
        y = coords[:, 1] * scale_factor
        
        if len(x) == 0:
            raise ValueError(
                f"No spatial coordinates found for basis '{basis_name}' in current data object.\n"
                f"Image key: '{img_key}'\n"
                f"Library ID: '{library_id}'\n"
                "Check that the data contains tissue spots/cells."
            )
        
        
        
        # Crop to scatter bounds
        x_min, x_max = x.min(), x.max()
        y_min, y_max = y.min(), y.max()
        img = img[int(y_min):int(y_max), int(x_min):int(x_max)].copy()
        x -= x_min
        y -= y_min

        # -----------------------------
        # Image-only transforms
        # -----------------------------
        img, _, _, Hf, Wf = apply_affine_to_image_and_coords(
            img=img, x=None, y=None, transforms=image_orientation
        )

        # -----------------------------
        # Scatter-only transforms
        # -----------------------------
        _, x, y, Hf, Wf = apply_affine_to_image_and_coords(
            img=None, x=x, y=y, transforms=scatter_orientation
        )

        # -----------------------------
        # Shared transforms (both)
        # -----------------------------
        img_final, x_final, y_final, Hf, Wf = apply_affine_to_image_and_coords(
            img=img, x=x, y=y, transforms=shared_orientation
        )
        # Determine plotting limits
        xmin, xmax = x_final.min(), x_final.max()
        ymin, ymax = y_final.min(), y_final.max()

        # Normalize color to always be a list
        color_list = [None] if color is None else ([color] if isinstance(color, str) else list(color))

        # Handle subplot grid
        ncols_to_use = 1 if color is None or len(color_list) == 1 else ncols
        nrows = ceil(len(color_list) / ncols_to_use)

        figsize_to_use = (6 * ncols_to_use, 6 * nrows) if figsize is None else figsize

        fig, axes = plt.subplots(nrows, ncols_to_use, figsize=figsize_to_use)
        axes = np.atleast_1d(axes).flatten()

        for ax, color_col in zip(axes, color_list):
            ax.imshow(img_final, origin='lower')

            if color_col is None:
                ax.set_title("Image only")
                ax.set_xlim(xmin, xmax)
                ax.set_ylim(ymin, ymax)
                if not plot_axes:
                    ax.axis("off")
                continue

            values = get_values(adata, color_col, use_raw=use_raw)
            values = np.ravel(values)

            # Sort for plotting
            if np.issubdtype(values.dtype, np.number):
                order = np.argsort(values)[::-1]
            else:
                if groups is not None:
                    cat_map = {cat: i for i, cat in enumerate(groups)}
                    mapped = [cat_map[v] for v in values]
                    order = np.argsort(mapped)
                else:
                    order = np.argsort(values)

            values = values[order]
            xp_sorted, yp_sorted = x_final[order], y_final[order]

            # Handle categorical
            if isinstance(values.dtype, pd.CategoricalDtype) or np.issubdtype(values.dtype, np.object_):
                if category_colors is not None:
                    color_dict = category_colors
                else:
                    unique_vals = np.unique(values)
                    palette = plt.cm.tab20.colors
                    color_dict = {cat: palette[i % len(palette)] for i, cat in enumerate(unique_vals)}

                colors_mapped = [color_dict[v] for v in values]
                scatter = ax.scatter(xp_sorted, yp_sorted, s=size, c=colors_mapped, edgecolor='none')

                if show_legend:
                    import matplotlib.patches as mpatches
                    handles = [mpatches.Patch(color=c, label=l) for l, c in color_dict.items()]
                    ax.legend(handles=handles, title="Legend", loc='center left', bbox_to_anchor=(1, 0.5))

            # Continuous
            elif np.issubdtype(values.dtype, np.number):
                scatter = ax.scatter(xp_sorted, yp_sorted, s=size, c=values, cmap=cmap, edgecolor='none')

                divider = make_axes_locatable(ax)
                cax = divider.append_axes("right", size="5%", pad=0.05)
                cbar = fig.colorbar(scatter, cax=cax)
                cbar.set_label(color_col, rotation=270, labelpad=10)

            else:
                raise TypeError("Unsupported data type for plotting")

            ax.set_title(f"Plotting {color_col}")
            ax.set_xlim(xmin, xmax)
            ax.set_ylim(ymin, ymax)
            if not plot_axes:
                ax.axis("off")

        # Hide extra axes
        for ax in axes[len(color_list):]:
            ax.axis("off")

        plt.tight_layout()

        if save:
            plt.savefig(save)
            plt.show()
        else:
            plt.show()

      
        
def has_non_integer(matrix) -> bool:
    """
    Check if a matrix contains any non-integer values.

    Parameters
    ----------
    matrix : np.ndarray or sparse matrix
        Matrix to check.

    Returns
    -------
    bool
        True if any non-integer values are present.
    """
    if issparse(matrix):
        # Only check the stored non-zero values
        return not np.all(np.mod(matrix.data, 1) == 0)
    else:
        # Dense matrix: check in chunks if very large
        chunk_size = 1000
        n_rows = matrix.shape[0]
        for start in range(0, n_rows, chunk_size):
            end = min(start + chunk_size, n_rows)
            chunk = matrix[start:end]
            if not np.all(np.mod(chunk, 1) == 0):
                return True
        return False
    

def convert_to_int(adata: AnnData) -> None:
    """
    Round entries of `adata.X` to nearest integer in-place.

    Handles dense or sparse matrices.

    Parameters
    ----------
    adata : AnnData
        Annotated data object with `X` matrix.

    Returns
    -------
    None
    """
    if issparse(adata.X):
        # Only round the non-zero entries
        adata.X.data = np.rint(adata.X.data)
    else:
        # Dense matrix: round in-place, optionally in chunks for very large matrices
        try:
            # Attempt in-place rounding (fast)
            np.rint(adata.X, out=adata.X)
        except MemoryError:
            # Chunked rounding if memory is tight
            chunk_size = 1000
            n_rows = adata.X.shape[0]
            for start in range(0, n_rows, chunk_size):
                end = min(start + chunk_size, n_rows)
                np.rint(adata.X[start:end], out=adata.X[start:end])

                
