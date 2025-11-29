import argparse
from pathlib import Path
import numpy as np
import pandas as pd
import cv2 as cv
from skimage import measure, morphology
from ridge_detector import RidgeDetector
import pipeline_utils

# --- CONFIGURATION ---
# Steger's Algorithm Parameters
RIDGE_WIDTHS = [1, 2, 3, 4]  # Range of line widths to detect (in pixels)
RIDGE_CONTRAST_LOW = 30  # Lower threshold for hysteresis
RIDGE_CONTRAST_HIGH = 120  # Upper threshold for hysteresis

# Reconstruction Limits
MAX_FILAMENT_THICKNESS = 2  # Cap to prevent artifacts from becoming blobs
MIN_FILAMENT_THICKNESS = 2  # Minimum visibility
WIDTH_SCALE_FACTOR = 1  # Multiplier if we want to thin the lines (1.0 = true width)


def normalize_to_8bit(image):
    """Robust 8-bit conversion."""
    if image.dtype == np.uint8: return image

    min_val, max_val = float(image.min()), float(image.max())
    # Prevent noise amplification in empty images
    if (max_val - min_val) < pipeline_utils.MIN_DYNAMIC_RANGE:
        return np.zeros_like(image, dtype=np.uint8)

    scaled = (image - min_val) / (max_val - min_val) * 255.0
    return np.clip(scaled, 0, 255).astype(np.uint8)


def draw_ridges(detector, shape):
    """Reconstructs binary mask from detected ridge objects."""
    output = np.zeros(shape, dtype=np.uint8)

    for contour in detector.contours:
        # Default thickness
        thickness = MIN_FILAMENT_THICKNESS

        # Calculate thickness based on Steger's width estimation
        if hasattr(contour, 'width_l') and hasattr(contour, 'width_r'):
            # Average total width along the contour segment
            avg_raw_width = np.mean(contour.width_l + contour.width_r)

            # Apply scaling and clamping
            calculated_thickness = round(avg_raw_width * WIDTH_SCALE_FACTOR)
            thickness = int(np.clip(calculated_thickness, MIN_FILAMENT_THICKNESS, MAX_FILAMENT_THICKNESS))
            thickness = 1

        # Draw the line segment
        pts = np.stack([contour.col, contour.row], axis=1).astype(np.int32)
        cv.polylines(output, [pts], isClosed=False, color=255, thickness=thickness)

    return output


def analyze_networks(image, mask):
    # 1. Normalize
    im_8bit = normalize_to_8bit(image)

    # 2. Run Steger's Ridge Detection
    det = RidgeDetector(
        line_widths=RIDGE_WIDTHS,
        low_contrast=RIDGE_CONTRAST_LOW,
        high_contrast=RIDGE_CONTRAST_HIGH,
        min_len=5,  # Ignore ridges shorter than this length
        max_len=0,  # Ignore ridges longer than this length, set to 0 for no limit
        dark_line=False,  # Set to True if detecting black ridges in white background, False otherwise
        estimate_width=True,  # Estimate width for each detected ridge point
        extend_line=True,  # Tend to preserve ridges near junctions if set to True
        correct_pos=False,  # Correct ridge positions with asymmetric widths if set to True
    )
    det.detect_lines(im_8bit)

    # 3. Create Binary Mask
    im_ridges = draw_ridges(det, im_8bit.shape)

    # 4. Strict Masking: Remove ridges found outside tissue mask
    im_ridges = pipeline_utils.get_foreground(im_ridges, mask)

    # 5. Quantification
    im_bool = im_ridges > 0
    labels = morphology.label(im_bool)

    props = pd.DataFrame(measure.regionprops_table(
        labels, intensity_image=im_bool, properties=('label', 'area', 'mean_intensity')
    ))

    # Calculate metrics
    tissue_area = np.count_nonzero(mask)
    total_pixels = 0
    structure_index = 0

    if tissue_area > 0 and not props.empty:
        total_pixels = props['area'].sum()
        # Structure index = Integrated Density / Tissue Area
        structure_index = (props['area'] * props['mean_intensity']).sum() / tissue_area

    stats = pd.DataFrame({
        'total_ridge_pixels': [total_pixels],
        'structure_index': [structure_index],
        'object_count': [len(np.unique(mask)) - 1],
        'tissue_area': [tissue_area]
    })

    return {
        'network_binary': im_ridges,
        'network_overlay': pipeline_utils.create_overlay(im_8bit, im_ridges),
        'network_stats.csv': stats
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("parent_dir", type=Path)
    parser.add_argument("--cores", type=int, default=8)
    args = parser.parse_args()

    tasks = pipeline_utils.build_tasks(args.parent_dir)
    pipeline_utils.run_pipeline(tasks, analyze_networks, args.cores)