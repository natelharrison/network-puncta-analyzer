import argparse
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.spatial import cKDTree
from skimage import feature
import pipeline_utils

# --- CONFIGURATION ---
LOG_MIN_SIGMA = 1.5  # Min blob size
LOG_MAX_SIGMA = 2.5  # Max blob size (Increased from 2)
LOG_NUM_SIGMA = 5  # Steps (Increased for better sizing)
LOG_THRESHOLD = 0.025  # Threshold (Stricter: 0.02 vs 0.002)
LOG_OVERLAP = 1  # Allowed overlap fraction (tuned)
EDGE_MARGIN = 5  # Pixels to ignore near image borders to prevent edge artifacts


def analyze_puncta(image, mask):
    # 1. Normalize to float32 for stable LoG detection
    img_float = image.astype(np.float32, copy=False)
    if img_float.max() > 0:
        img_norm = img_float / img_float.max()
    else:
        img_norm = img_float

    # 2. Blob Detection
    blobs = feature.blob_log(
        img_norm,
        min_sigma=LOG_MIN_SIGMA,
        max_sigma=LOG_MAX_SIGMA,
        num_sigma=LOG_NUM_SIGMA,
        threshold=LOG_THRESHOLD,
        overlap=LOG_OVERLAP,
        exclude_border=5
    )

    # 3. Filter by Mask
    if len(blobs) > 0:
        ys = blobs[:, 0]
        xs = blobs[:, 1]
        sigmas = blobs[:, 2]

        # Drop detections whose extent (approx radius ~ sqrt(2)*sigma) touches borders
        H_img, W_img = img_norm.shape
        radii = np.sqrt(2) * sigmas
        inside_y = (ys - radii >= EDGE_MARGIN) & (ys + radii < H_img - EDGE_MARGIN)
        inside_x = (xs - radii >= EDGE_MARGIN) & (xs + radii < W_img - EDGE_MARGIN)
        keep = inside_y & inside_x
        blobs = blobs[keep]

        # Convert to ints after filtering
        ys = blobs[:, 0].astype(int)
        xs = blobs[:, 1].astype(int)

        H, W = mask.shape
        ys = np.clip(ys, 0, H - 1)
        xs = np.clip(xs, 0, W - 1)

        is_background = mask[ys, xs] == 0
        blobs = blobs[is_background]

    # 4. Metrics
    count = len(blobs)
    bg_area = np.count_nonzero(mask == 0)
    avg_nnd = 0.0

    if count > 1:
        coords = blobs[:, :2]
        tree = cKDTree(coords)
        dists, _ = tree.query(coords, k=2)
        avg_nnd = np.mean(dists[:, 1])

    stats = pd.DataFrame({
        'puncta_count': [count],
        'avg_nnd': [avg_nnd],
        'puncta_density': [count / bg_area if bg_area > 0 else 0],
        'search_area': [bg_area]
    })

    # 5. Generate Label Mask
    puncta_labels = np.zeros_like(image, dtype=np.uint16)

    if count > 0:
        ys = blobs[:, 0].astype(int)
        xs = blobs[:, 1].astype(int)

        # Create IDs from 1 to Count
        ids = np.arange(1, count + 1, dtype=np.uint16)

        # Assign IDs to the coordinates
        # If duplicates exist due to rounding, the last one wins (counts are preserved in CSV)
        puncta_labels[ys, xs] = ids


    return {
        'puncta_labels': puncta_labels,  # Return the explicit label map
        'puncta_coords.csv': pd.DataFrame(blobs, columns=['y', 'x', 'sigma']),
        'puncta_stats.csv': stats
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("parent_dir", type=Path)
    parser.add_argument("--cores", type=int, default=8)
    args = parser.parse_args()

    tasks = pipeline_utils.build_tasks(args.parent_dir)
    pipeline_utils.run_pipeline(tasks, analyze_puncta, args.cores)
