import tifffile
import numpy as np
import cv2 as cv
import matplotlib.pyplot as plt
import math
from pathlib import Path
from multiprocessing import Pool
from tqdm import tqdm

# --- CONFIGURATION ---
MIN_DYNAMIC_RANGE = 30  # Minimum intensity difference (max-min) to process an image
MONTAGE_COLS = 5  # Number of columns in QC montage


def imread(path):
    return np.squeeze(tifffile.imread(path))


def imwrite(path, image):
    path.parent.mkdir(parents=True, exist_ok=True)
    tifffile.imwrite(path, image)


def get_foreground(image, mask):
    """Returns image with background pixels zeroed out."""
    return np.where(mask != 0, image, 0)


def get_background(image, mask):
    """Returns image with foreground (cells/network) zeroed out."""
    return np.where(mask == 0, image, 0)


def create_overlay(raw, mask, color=(0, 255, 255)):
    """Creates a cyan overlay of the mask onto the raw image."""
    if raw.dtype != np.uint8:
        raw = cv.normalize(raw, None, 0, 255, cv.NORM_MINMAX, dtype=cv.CV_8U)
    rgb = cv.cvtColor(raw, cv.COLOR_GRAY2RGB)
    rgb[mask > 0] = color
    return rgb


def create_montage(trial_name, image_paths, output_path):
    """Generates a grid montage sorted by genotype (WT vs Mutant)."""
    if not image_paths: return

    wts = sorted([p for p in image_paths if 'wt' in get_genotype_from_name(p.name).lower()])
    muts = sorted([p for p in image_paths if 'mutant' in get_genotype_from_name(p.name).lower()])

    rows_wt = math.ceil(len(wts) / MONTAGE_COLS) if wts else 0
    rows_mut = math.ceil(len(muts) / MONTAGE_COLS) if muts else 0
    total_rows = rows_wt + rows_mut

    if total_rows == 0: return

    fig, axes = plt.subplots(total_rows, MONTAGE_COLS, figsize=(MONTAGE_COLS * 3, total_rows * 3.5))
    fig.suptitle(f"{trial_name}", fontsize=20, fontweight='bold')

    axes_flat = axes.flatten() if total_rows * MONTAGE_COLS > 1 else [axes]
    current_idx = 0

    def plot_group(paths, label_color):
        nonlocal current_idx
        for p in paths:
            ax = axes_flat[current_idx]
            try:
                ax.imshow(imread(p))
                ax.set_title(p.stem.replace("MAX_chan1_", "")[-15:], color=label_color, fontsize=8, fontweight='bold')
            except Exception:
                ax.text(0.5, 0.5, "Error", ha='center')
            ax.axis('off')
            current_idx += 1

        while current_idx % MONTAGE_COLS != 0:
            axes_flat[current_idx].axis('off')
            current_idx += 1

    if wts: plot_group(wts, '#555555')
    if muts: plot_group(muts, '#E63946')

    plt.tight_layout()
    plt.savefig(output_path, dpi=600)
    plt.close(fig)


def get_genotype_from_name(filename):
    clean = filename.lower().replace('_', ' ').replace('-', ' ')
    if any(x in clean.split() for x in ['v', 'v1', 'v2']):
        return 'Mutant'
    return 'WT'


def _worker(args):
    paths, func = args
    img_path, mask_path = paths

    image = imread(img_path)
    mask = imread(mask_path)

    results = func(image, mask)
    if not results: return

    # Place outputs next to the folder containing MIPs (works for Trial/MIPs and Trial/Genotype/MIPs)
    results_root = None
    for parent in img_path.parents:
        if parent.name == "MIPs":
            # Write alongside the owner of MIPs (genotype folder if present, otherwise trial)
            results_root = parent.parent / "analysis_results"
            break
    if results_root is None:
        return

    for key, data in results.items():
        if hasattr(data, 'to_csv'):
            save_dir = results_root / 'quantification'
            save_dir.mkdir(parents=True, exist_ok=True)
            data.to_csv(save_dir / f"{img_path.stem}_{key}", index=False)
        else:
            save_path = results_root / key / img_path.name
            imwrite(save_path, data)


def build_tasks(parent_dir):
    tasks = []
    for img in parent_dir.rglob("MIPs/chan1/*.tif"):
        mask = img.parent.parent / "chan1_masks" / img.name
        if mask.exists():
            tasks.append((img, mask))
    return tasks


def run_pipeline(tasks, process_func, cores=8):
    if not tasks: return
    print(f"Processing {len(tasks)} files using {cores} cores...")

    pool_args = [(t, process_func) for t in tasks]

    with Pool(processes=int(cores)) as pool:
        for _ in tqdm(pool.imap_unordered(_worker, pool_args), total=len(tasks)):
            # tqdm consumes the iterator; results are written inside _worker
            pass
