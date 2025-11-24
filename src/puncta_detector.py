import argparse
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.spatial import cKDTree
from skimage import feature, morphology
import pipeline_utils

# --- CONFIGURATION ---
# Laplacian of Gaussian (LoG) settings
LOG_MIN_SIGMA = 1  # Min blob size
LOG_MAX_SIGMA = 2  # Max blob size
LOG_NUM_SIGMA = 5  # Steps between min/max
LOG_THRESHOLD = 0.002  # Intensity threshold for detection
LOG_OVERLAP = 0.5  # Allowed overlap fraction


def analyze_puncta(image, mask):
    # 1. Masking: Analyze BACKGROUND only
    im_background = pipeline_utils.get_background(image, mask)
    bg_area = np.count_nonzero(mask == 0)

    # 2. Blob Detection
    blobs = feature.blob_log(
        im_background,
        min_sigma=LOG_MIN_SIGMA,
        max_sigma=LOG_MAX_SIGMA,
        num_sigma=LOG_NUM_SIGMA,
        threshold=LOG_THRESHOLD,
        overlap=LOG_OVERLAP
    )

    # 3. Metrics
    count = len(blobs)
    avg_nnd = 0.0

    if count > 1:
        # Nearest Neighbor Distance
        coords = blobs[:, :2]
        tree = cKDTree(coords)
        dists, _ = tree.query(coords, k=2)  # k=2 gets the nearest neighbor (k=1 is self)
        avg_nnd = np.mean(dists[:, 1])

    stats = pd.DataFrame({
        'puncta_count': [count],
        'avg_nnd': [avg_nnd],
        'puncta_density': [count / bg_area if bg_area > 0 else 0],
        'search_area': [bg_area]
    })

    # 4. Generate visual label mask
    puncta_mask = np.zeros_like(image, dtype=np.uint16)
    if count > 0:
        # Mark blob centers
        puncta_mask[blobs[:, 0].astype(int), blobs[:, 1].astype(int)] = 1
        # Dilate for visibility (radius=1)
        puncta_mask = morphology.dilation(puncta_mask, morphology.disk(1))

    return {
        'puncta_labels': morphology.label(puncta_mask).astype(np.uint16),
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