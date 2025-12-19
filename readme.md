# Network & Puncta Quantification Pipeline

A high-throughput pipeline for quantifying **filamentous networks** and **puncta** in 2D fluorescence microscopy images.  
All metrics are **normalized to internal WT controls per trial**, enabling robust cross-experiment comparisons.

The pipeline is designed for batch processing, reproducibility, and QC-driven analysis, with optional interactive tissue masking via **Micro-SAM**.

---

## Features

- End-to-end processing from raw OIB files to quantified CSVs
- Micro-SAM–assisted tissue masking (CPU or GPU)
- Parallelized batch analysis
- Automated QC overlays and montages
- Per-image and per-trial summary statistics

---

## Installation

### 1. Clone the repository
```text
git clone <repo_url>
cd abboud_project
```

### 2. Create the Conda environment
```text
conda env create -f environment.yml
```

### 3. Activate the environment
```text
conda activate mucin_env
```

### 4. (Optional) Install Micro-SAM support

Required only if you plan to generate tissue masks interactively.

**CPU**
```text
conda install -c conda-forge micro_sam
```

**GPU**
```text
conda install pytorch torchvision torchaudio pytorch-cuda=11.8 -c pytorch -c nvidia
```

---

## Data Organization

### Raw data (input)

Original `.oib` files are organized by trial:

```text
<Project_Root>/
├── Trial_01/    # original .oib files
├── Trial_02/
└── Trial_03/
```

---

### Analysis-ready layout

After conversion, each trial should follow this structure:

```text
<Project_Root>/
├── Trial_01/
│   └── MIPs/
│       ├── chan0/          # secondary channel (unused)
│       ├── chan1/          # primary channel TIFs
│       └── chan1_masks/    # tissue masks (manual or Micro-SAM)
├── Trial_02/
│   └── MIPs/...
└── Trial_03/
    └── MIPs/...
```

Notes:
- Mask filenames **must exactly match** the corresponding images in `chan1`
- Masks may be binary or instance-labeled
- Any non-zero pixel is treated as tissue

---

## Workflow

### Step 1 — Convert OIB files to max-projection TIFs (optional)

If starting from raw `.oib` files:

```text
python src/convert_oibs.py "<project_root>"
```

This step:
- Recursively scans for `.oib` files
- Generates max-projection TIFs
- Writes outputs to `MIPs/chan0` and `MIPs/chan1` per trial

Skip this step if MIPs already exist.

---

### Step 2 — Generate tissue masks with Micro-SAM

Masks are created interactively using napari.

1. Launch napari from the active environment:
```text
napari
```

2. Open:
```text
Plugins > Segment Anything for Microscopy > Image Series Annotator
```

3. Configure per trial:
- **Input folder**  
  `<Project_Root>/Trial_X/MIPs/chan1`
- **Output folder**  
  `<Project_Root>/Trial_X/MIPs/chan1_masks`

4. In **Advanced settings**:
- Set **Custom weights path** to your SAM checkpoint  
  (e.g. `sam_vit_b_01ec64.pth`)
- Select a model appropriate for your CPU/GPU

5. Annotate tissue using boxes, points, or scribbles  
6. Refine masks with napari paint/erase tools as needed  
7. Save masks to disk

---

### Step 3 — Run automated analysis

Execute the full analysis pipeline:

```text
python run_analysis.py "<project_root>" --cores 8
```

This step will:
- Validate that all images have matching masks
- Detect puncta and filamentous networks
- Compile per-image and per-trial statistics
- Generate QC overlays and summary plots

Adjust `--cores` to control multiprocessing.

---

## Outputs

Each trial produces an `analysis_results/` directory:

```text
analysis_results/
├── network_binary/          # binary ridge/network masks
├── network_overlay/         # QC overlays
├── puncta_labels/           # labeled puncta detections
├── quantification/
│   ├── per_image_*.csv
│   ├── DATA_Networks_Master.csv
│   └── DATA_Puncta_Master.csv
├── Montage_Network.png
└── PLOT_Summary.png
```

---

## Notes & Best Practices

- Ensure tissue masks fully cover valid regions to avoid biased normalization
- Large images benefit from GPU-accelerated Micro-SAM
- Inspect QC overlays to verify segmentation quality
- WT normalization is performed **within each trial**, not globally

---

## Citation

If you use this pipeline in a publication, please cite the associated study and acknowledge Micro-SAM where applicable.

---

## Extensions

The pipeline is modular and can be extended to support:
- Alternative network detectors
- Additional per-cell or per-region metrics
- 3D or time-series data

Contact the maintainer for guidance.
