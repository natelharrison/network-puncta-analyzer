# Network & Puncta Quantification Pipeline

A high-throughput pipeline for quantifying filamentous networks and puncta in 2D fluorescence microscopy images.
Metrics are normalized to internal WT controls per trial. Includes optional Micro-SAM-assisted masking and batch QC.

---
## Methodology
- Networks: Steger ridge detector (via `ridge_detector`) with configurable widths/contrast; masked to non-zero pixels in `chan1_masks`.
- Puncta: Laplacian of Gaussian blob detector (`skimage.feature.blob_log`) with configurable sigmas/threshold/overlap; masked to non-zero pixels in `chan1_masks`.
- Normalization: Within each comparison group, metrics are divided by the WT mean; mapping of WT/Mutant comes from folder names or filenames.
- Reproducibility: Run parameters and environment details are recorded in `analysis_runs/<run_id>/run_meta.json`.

---
## Features
- End-to-end from OIB → MIPs → masks → metrics/plots
- Micro-SAM tissue masking (CPU or GPU)
- Parallelized batch processing
- Per-image stats plus combined CSVs and plots with run metadata

---
## Installation
1) Clone/download and enter the repo
```
git clone <repo_url>
cd abboud_project
```
2) Create env
```
conda env create -f environment.yml
```
3) Activate
```
conda activate mucin_env
```
4) (Optional) Micro-SAM + GPU
```
conda install -c conda-forge micro_sam           # CPU
conda install pytorch torchvision torchaudio pytorch-cuda=11.8 -c pytorch -c nvidia  # GPU (optional)
```

---
## Data Organization

### Raw OIBs
```
<Project_Root>/
├─ Trial_01/    # original .oib files
├─ Trial_02/
└─ Trial_03/
```

### Analysis-ready layout (per-genotype folders under each trial)
```
<Project_Root>/
├─ Trial_01/
│  ├─ Rosa/               # WT example
│  │  ├─ MIPs/
│  │  │  ├─ chan0/
│  │  │  ├─ chan1/
│  │  │  └─ chan1_masks/
│  │  └─ analysis_results/
│  └─ V1/                 # Mutant example
│     ├─ MIPs/
│     │  ├─ chan0/
│     │  ├─ chan1/
│     │  └─ chan1_masks/
│     └─ analysis_results/
└─ Trial_02/
   └─ ...
```
Notes:
- Mask filenames must exactly match the images in `chan1`
- Masks may be binary or instance-labeled; any non-zero pixel is treated as tissue

---
## Workflow

### 1) Convert OIBs to max-projection TIFs (optional)
```
python src/convert_oibs.py "<project_root>" --skip-existing
```
- Scans for `.oib` and writes max-projections to `MIPs/chan0` and `MIPs/chan1`.
- `--skip-existing` avoids overwriting existing TIFs.

### 2) Create tissue masks with Micro-SAM in napari
- Launch: `napari`
- Plugins > Segment Anything for Microscopy > Image Series Annotator
- Per genotype folder:
  - Input: `<Project_Root>/Trial_X/<Genotype>/MIPs/chan1`
  - Output: `<Project_Root>/Trial_X/<Genotype>/MIPs/chan1_masks`
- Advanced: set Custom weights path to your SAM checkpoint (e.g., `sam_vit_b_01ec64.pth`), pick CPU/GPU model.
- Annotate (boxes/points/scribbles), refine with paint/erase, save masks with matching filenames.

### 3) Run automated analysis
```
python run_analysis.py "<project_root>" --cores 8 --run-id <optional_name>
```
- Finds every `MIPs/chan1/*.tif` (trial/genotype aware), checks for matching masks.
- Writes per-image outputs to each genotype’s `analysis_results/`.
- Writes combined stats/plots/metadata to `<project_root>/analysis_runs/<run_id>/` (timestamp if not set).

---
## Outputs
- Per genotype `analysis_results/`:
  - `network_binary/`, `network_overlay/`, `puncta_labels/`
  - `quantification/` per-image CSVs
- Project-level `analysis_runs/<run_id>/`:
  - `DATA_Networks_Combined.csv`, `DATA_Puncta_Combined.csv`
  - `PLOT_Per_Trial_Breakdown.png`, `PLOT_Global_Summary.png`
  - `run_meta.json` (parameters, cores, platform, run id)

---
## Notes
- WT/Mutant mapping: folder names (`Rosa`, `Rosa 1/2`, `BL6`, `V`, `V1`, `V2`) or filenames. To change mapping, edit `FOLDER_GENOTYPE_MAP` and `detect_strain_info` in `src/compile_stats.py`.
- Normalization is within each comparison group; ensure WT controls exist per group.
- Inspect QC overlays and montages to confirm masking and detections.
