# Network & Puncta Quantification Pipeline

**DISCLAIMER:** Currently, the genotype distinction logic (Mutant vs. Wildtype) is hardcoded for specific projects. The analysis modules rely on filename patterns (e.g., containing 'v1'/'v2' for mutants and 'Rosa'/'BL6' for wildtypes) to assign groups. While the detection scripts should work with the correct file directory structure,  statistical compilation requires these specific naming conventions in this version. I plan on addressing this in the near future.

---

This software provides a high-throughput method for quantifying network and puncta phenotypes in 2D fluorescence microscopy images (or 2D intensity projection). The pipeline distinguishes between two distinct structural classes:

1.  **Filamentous Networks:** Curvilinear webs or ridge-like structures.
2.  **Discrete Puncta:** Distinct blobs or dots.

While originally designed for mucin analysis, the detection pipelines are agnostic to the biological target and can be applied to any images these morphological features. To account for batch effects and experimental variation across imaging sessions, all metrics are normalized against the mean of the Wildtype control group within each trial.

## Pipeline Summary

The pipeline executes the following logic:
1.  **Preprocessing:** Converts Olympus (`.oib`) files to TIFF format (optional) and applies.
2.  **Network Detection:** Utilizes **Steger’s Ridge Detection** to identify curvilinear structures.
3.  **Puncta Detection:** Utilizes **Laplacian of Gaussian (LoG)** blob detection, specifically restricted to background regions (areas excluding the detected network).
4.  **Statistical Normalization:** Aggregates data by trial and calculates metrics relative to the mean of the Wildtype control group.

---

## Normalization Logic: "Relative to WT"

We calculate a **Relative** value for every metric:

$$ \text{Relative Value} = \frac{\text{Raw Sample Value}}{\text{Mean of Internal WT Controls}} $$

* **WT Control Mean:** Always equals **1.0**.
* **Interpretation:** A value of **1.5** indicates a 50% increase relative to the baseline for that specific trial.

---

## Directory Structure & Inputs

The pipeline requires a specific directory hierarchy.

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
│   │   ├── chan1/                    # [Generated] Primary Channel (Gray TIFs)
│   │   │   ├── MAX_chan1_WT_Image_01.tif
│   │   │   └── MAX_chan1_Mutant_Image_01.tif
│   │   │
│   │   └── chan1_masks/              # [Manual] Tissue Segmentation Masks
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

### 1. Set Up Miniforge
We recommend **Miniforge** for package management.
1.  Download the [Miniforge3 installer](https://github.com/conda-forge/miniforge#miniforge3) for your OS.
2.  Run the installer (double-click on Windows, or `bash Miniforge3-*.sh` on macOS/Linux) and finish the prompts.
3.  Launch the **Miniforge Prompt/terminal**.

### 2. Create the Base Environment (conda + pip)
`environment.yml` includes the base stack (Napari, PyQt5, ridge-detector via pip, etc.). Create and activate it:

```bash
conda env create -f environment.yml
conda activate mucin_env
```

### 3. Optional: Micro-SAM
Needed in order to create segmentation masks for the network objects
* **Micro-SAM:**
    ```bash
    conda install -c conda-forge micro_sam
    ```
* **NVIDIA GPU (CUDA 11.8):**
    ```bash
    conda install pytorch torchvision torchaudio pytorch-cuda=11.8 -c pytorch -c nvidia
    ```


---

## Workflow

### Step 1: Convert Images (Optional)
Extracts Max Z-Projections from Olympus `.oib` files. This generates the `chan1` folder.

```bash
python src/convert_oibs.py "E:\Mic\Project_Name"
```

### Step 2: Segmentation (Manual)
Generate masks to exclude background noise.
1.  Launch Napari: `napari` (in Miniforge terminal with environment active).
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
