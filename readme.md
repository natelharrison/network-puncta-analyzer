# Network & Puncta Quantification Pipeline

This software provides a high-throughput method for quantifying network and puncta phenotypes in fluorescence microscopy images. The pipeline distinguishes between **Filamentous Networks** (webs/ridges) and **Discrete Puncta** (dots). Originally designed for Mucin analysis, this pipeline is applicable to any stain exhibiting ridge-like or blob-like structures. To account for batch effects and experimental variation, all metrics are normalized against internal Wildtype controls within each trial (reported as "Relative to WT").

## Pipeline Summary

The pipeline executes the following logic:
1.  **Preprocessing:** Converts Olympus (`.oib`) files to TIFF format (optional) and applies robust normalization.
2.  **Network Detection:** Utilizes Steger’s Ridge Detection to identify curvilinear structures based on intensity gradients.
3.  **Puncta Detection:** Utilizes Laplacian of Gaussian (LoG) blob detection on background regions (areas excluding the network).
4.  **Statistical Normalization:** Aggregates data by trial and calculates metrics relative to the mean of the Wildtype control group.

---

## Directory Structure & Inputs

The pipeline requires a specific directory hierarchy. The structure differs slightly depending on whether you start with raw `.oib` files or pre-processed `.tif` images.

### 1. Initial State (Raw OIBs)
If starting with Olympus files, place them directly inside your Trial folders.

```text
E:\...\Project_Name\
├── Trial_R13_Rosa_V1\
│   ├── WT_Control_Image_01.oib
│   ├── WT_Control_Image_02.oib
│   ├── Mutant_Image_01.oib
│   └── Mutant_Image_02.oib
└── Trial_R14_Rosa_V2\
    ├── (Corresponding .oib files...)
```

### 2. Final State (Ready for Analysis)
After running the conversion script (Step 1) and performing manual segmentation (Step 2), the directory will contain the necessary `chan1` and `chan1_masks` folders.

**Legend:**
* `[Generated]`: Created automatically by `convert_oibs.py`.
* `[Manual]`: Created by the user during segmentation.
* `[Output]`: Created automatically by `run_analysis.py`.

```text
E:\...\Project_Name\
├── Trial_R13_Rosa_V1\
│   │
│   ├── (Original .oib files)
│   │
│   ├── max_z_projections/
│   │   ├── chan0/                    # Secondary channel (unused)
│   │   ├── chan1/                    # [Generated] Primary Mucin Staining (Gray TIFs)
│   │   │   ├── MAX_chan1_WT_Image_01.tif
│   │   │   └── MAX_chan1_Mutant_Image_01.tif
│   │   │
│   │   └── chan1_masks/              # [Manual] Instanced Network Segmentation Masks
│   │       ├── MAX_chan1_WT_Image_01.tif
│   │       └── MAX_chan1_Mutant_Image_01.tif
│   │
│   └── analysis_results/             # [Output]
│       ├── Montage_Network.png       # QC Visualization
│       ├── network_binary/           # Detected ridges
│       ├── puncta_labels/            # Detected dots
│       └── quantification/           # Raw CSV data
│
└── Trial_R14_Rosa_V2\
    └── ...
```

---

## Installation

### 1. Python Environment (Miniforge)
We recommend **Miniforge** for package management.
1.  Download and install [Miniforge3](https://github.com/conda-forge/miniforge#miniforge3).
2.  Open your terminal (or Miniforge Prompt on Windows).

```bash
# Create environment
conda create -n mucin_env python=3.9
conda activate mucin_env

# Install dependencies
pip install -r requirements.txt
```

### 2. GPU Support (Required for Micro-SAM)
To run the segmentation tools efficiently, install the PyTorch version compatible with your hardware.

* **Windows/Linux (NVIDIA GPU):**
    ```bash
    pip3 install torch torchvision torchaudio --index-url [https://download.pytorch.org/whl/cu118](https://download.pytorch.org/whl/cu118)
    ```
* **MacOS (M1/M2):** Standard installation.

---

## Workflow

### Step 1: Convert Images (Optional)
Extracts Max Z-Projections from Olympus `.oib` files. This generates the `chan1` folder.

```bash
python src/convert_oibs.py "E:\Mic\Project_Name"
```

### Step 2: Segmentation (Manual)
Generate masks to exclude background noise.
1.  Launch Napari: `napari` in miniforge terminal with environment open.
2.  Load the `chan1` folder.
3.  Use the **Micro-SAM** plugin (Batch Segment) to annotate tissue.
4.  **Save Output:** Ensure masks are saved in a folder named `chan1_masks` at the same level as `chan1`. File names must match exactly.

### Step 3: Automated Analysis
Executes detection algorithms and compiles statistics.

**Windows Tip: Copying Paths**
1.  Hold **Shift** and **Right-Click** the experiment folder.
2.  Select **"Copy as path"**.
3.  Paste into the command line.

```bash
python run_analysis.py "E:\Mic\Project_Name"
```

---

## Output Description

Results are saved in `analysis_results` within each trial folder, plus summary files in the root directory.

* **`Montage_Network.png`:** Visual quality control showing network segmentation overlaid on raw images.
* **`DATA_Networks_Master.csv`:** Aggregated metrics for network structures (Total Pixels, Structure Index).
* **`DATA_Puncta_Master.csv`:** Aggregated metrics for puncta (Count, Density, Nearest Neighbor Distance).
* **`PLOT_Summary.png`:** Statistical visualization (Strip plots with Mean ± SE).

## Methodology Reference

| Metric | Description | Algorithm |
| :--- | :--- | :--- |
| **Network Area** | Total pixels belonging to curvilinear structures. | Steger's Ridge Detection |
| **Structure Index** | Measures network density/clumping. | $\sum (Area \times Intensity) / TissueArea$ |
| **Puncta Density** | Count of discrete dots per unit area. | Laplacian of Gaussian (LoG) |
| **Avg NND** | Average Nearest Neighbor Distance (Clustering). | k-d tree spatial query |
