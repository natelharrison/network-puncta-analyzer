# Network & Puncta Quantification Pipeline

This project provides high-throughput quantification of filamentous networks and puncta in 2D fluorescence images. Metrics are normalized to internal WT controls per trial.

## Pipeline Summary
1) Preprocess: convert Olympus .oib to max-projection TIFs (optional).
2) Network detection: Steger ridge detector.
3) Puncta detection: Laplacian of Gaussian.
4) Normalize metrics vs WT per trial.

## Directory Structure & Inputs

### Initial (raw OIBs)
```
E:\...\Project_Name\
??? Trial_R13_Rosa_V1\ (original .oib files)
??? Trial_R14_Rosa_V2\ (original .oib files)
??? Trial_R15_BL6_V1\ (original .oib files)
```

### Ready for analysis
```
E:\...\Project_Name\
??? Trial_R13_Rosa_V1\
?   ??? MIPs/
?       ??? chan0/          # Secondary channel (unused)
?       ??? chan1/          # [Generated] primary channel TIFs
?       ??? chan1_masks/    # [Manual] tissue masks (instance labels)
??? Trial_R14_Rosa_V2\
?   ??? MIPs/...
??? Trial_R15_BL6_V1\
    ??? MIPs/...
```

## Installation
1) Install Miniforge and open its prompt.
2) Create env:
```
conda env create -f environment.yml
conda activate mucin_env
```
3) Optional: add Micro-SAM and CUDA PyTorch if needed for mask creation:
```
conda install -c conda-forge micro_sam
conda install pytorch torchvision torchaudio pytorch-cuda=11.8 -c pytorch -c nvidia
```

## Workflow
### Step 1: Convert OIBs (optional)
```
python src/convert_oibs.py "E:\Mic\Project_Name"
```
### Step 2: Segmentation (manual)
- Launch `napari`
- Load `MIPs/chan1`
- Use Micro-SAM plugin to segment tissue
- Save masks to `MIPs/chan1_masks` with matching filenames

### Step 3: Automated analysis
```
python run_analysis.py "E:\Mic\Project_Name"
```

## Outputs
- `analysis_results/` under each trial
  - `network_binary/`, `puncta_labels/`, `quantification/`
  - `Montage_Network.png`, `DATA_Networks_Master.csv`, `DATA_Puncta_Master.csv`, `PLOT_Summary.png`

## Notes
- Instance masks are assumed for `chan1_masks`.
- Structure index currently uses binary ridge intensity; adjust if you need grayscale-weighted metrics.
