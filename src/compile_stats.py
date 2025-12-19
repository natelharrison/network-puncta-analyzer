import argparse
from pathlib import Path
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib
matplotlib.use("Agg")  # headless backend for script usage
import matplotlib.pyplot as plt
from scipy.stats import ttest_ind
from tqdm import tqdm
import pipeline_utils
import traceback
import json
import sys
import platform
from datetime import datetime
import network_detector
import puncta_detector


def get_genotype(filename):
    return pipeline_utils.get_genotype_from_name(filename)


def detect_strain_info(filename):
    clean = filename.lower().replace('_', ' ').replace('-', ' ').replace('.', ' ')
    if 'rosa 1' in clean: return 'Rosa 1', 'WT'
    if 'rosa 2' in clean: return 'Rosa 2', 'WT'
    if 'rosa' in clean:   return 'Rosa', 'WT'
    if 'bl6' in clean:    return 'BL6', 'WT'
    if 'v1' in clean:     return 'V1', 'Mutant'
    if 'v2' in clean:     return 'V2', 'Mutant'
    return 'Unknown', 'Unknown'


FOLDER_GENOTYPE_MAP = {
    # Wildtype keys
    'rosa': ('Rosa', 'WT'),
    'rosa 1': ('Rosa 1', 'WT'),
    'rosa 2': ('Rosa 2', 'WT'),
    'bl6': ('BL6', 'WT'),
    # Mutant keys
    'v': ('V', 'Mutant'),
    'v1': ('V1', 'Mutant'),
    'v2': ('V2', 'Mutant'),
}


def assign_comparison_group(trial_name, strain):
    trial_clean = trial_name.upper()
    if 'R14' in trial_clean:
        if strain in ['Rosa 1', 'V1']:
            return "R14 (Rosa1 vs V1)"
        if strain in ['Rosa 2', 'V2']:
            return "R14 (Rosa2 vs V2)"
    if 'R13' in trial_clean:
        return "R13 (Rosa vs V1)"
    if 'R15' in trial_clean:
        return "R15 (BL6 vs V1)"
    return trial_name


def normalize_by_group(df, metric_col, group_col='comparison_group', ctrl_label='WT'):
    """Normalizes data so Control Mean = 1.0 per group."""
    norm_col = f'{metric_col}_norm'
    df[norm_col] = np.nan
    for group in df[group_col].unique():
        mask_group = df[group_col] == group
        mask_ctrl = mask_group & (df['genotype'] == ctrl_label)
        ctrl_mean = df.loc[mask_ctrl, metric_col].mean()
        if ctrl_mean and ctrl_mean > 0:
            df.loc[mask_group, norm_col] = df.loc[mask_group, metric_col] / ctrl_mean
    return norm_col


def load_and_combine(data_dir, suffix):
    files = list(data_dir.rglob(f"*{suffix}"))
    if not files:
        return None

    dfs = []
    for p in files:
        if p.name.startswith("DATA_") or p.name.startswith("SUMMARY_") or "Combined" in p.name:
            continue

        try:
            df = pd.read_csv(p)

            # Determine trial/genotype from filesystem
            analysis_root = None
            for parent in p.parents:
                if parent.name == "analysis_results":
                    analysis_root = parent
                    break

            trial_name = "Unknown"
            genotype_folder = "NA"
            if analysis_root:
                mips_owner = analysis_root.parent  # could be Trial or Genotype
                has_mips_here = (mips_owner / "MIPs").exists()
                has_mips_parent = mips_owner.parent and (mips_owner.parent / "MIPs").exists()

                if has_mips_here and not has_mips_parent:
                    # Genotype-level: Trial/Genotype/{MIPs,analysis_results}
                    trial_name = mips_owner.parent.name if mips_owner.parent else "Unknown"
                    genotype_folder = mips_owner.name
                elif has_mips_here and has_mips_parent:
                    # Uncommon nested case: prefer parent as trial
                    trial_name = mips_owner.parent.name
                    genotype_folder = mips_owner.name
                elif not has_mips_here and has_mips_parent:
                    # Trial-level: Trial/{MIPs,analysis_results}
                    trial_name = mips_owner.name
                    genotype_folder = "NA"
                else:
                    # Fallback: treat analysis_results parent as trial
                    trial_name = mips_owner.name
                    genotype_folder = "NA"

            image_name = p.name.replace(f"_{suffix}", "")
            strain, genotype = detect_strain_info(p.name)

            # Override strain/genotype if folder provides a mapping
            folder_key = genotype_folder.lower()
            if folder_key in FOLDER_GENOTYPE_MAP:
                strain, genotype = FOLDER_GENOTYPE_MAP[folder_key]

            group_name = assign_comparison_group(trial_name, strain)

            df.insert(0, 'image_name', image_name)
            df['trial'] = trial_name
            df['folder_group'] = genotype_folder
            df['strain'] = strain
            df['genotype'] = genotype
            df['comparison_group'] = group_name
            dfs.append(df)
        except Exception as e:
            print(f"[!] Skipped file {p}: {e}")

    return pd.concat(dfs, ignore_index=True) if dfs else None


def plot_metric_polished(ax, data, x_col, y_col, hue_col, title, ylabel, show_legend=False):
    colors = {"WT": "#444444", "Mutant": "#D62828"}

    # Filter out Unknowns to keep plots clean
    data = data[data['genotype'].isin(["WT", "Mutant"])]

    order = sorted(data[x_col].unique())
    hue_order = ["WT", "Mutant"]

    if 'Norm' in title:
        ax.axhline(1.0, linestyle=':', color='black', alpha=0.4, linewidth=1)

    sns.stripplot(
        data=data,
        x=x_col, y=y_col, hue=hue_col,
        palette=colors, order=order, hue_order=hue_order,
        dodge=True, jitter=True, size=6, alpha=0.4,
        ax=ax, legend=False, zorder=0
    )

    sns.pointplot(
        data=data,
        x=x_col, y=y_col, hue=hue_col,
        palette=colors, order=order, hue_order=hue_order,
        dodge=0.4, capsize=0.1,
        errorbar='se', markers="D", linestyles="",
        markersize=6,
        ax=ax, legend=False, zorder=10
    )

    ax.set_title(title, fontweight='bold', fontsize=14, pad=10)
    ax.set_ylabel(ylabel, fontsize=12)
    ax.set_xlabel("")

    # Statistics
    for i, group in enumerate(order):
        sub = data[data[x_col] == group]
        wt = sub[sub['genotype'] == 'WT'][y_col].dropna()
        mut = sub[sub['genotype'] == 'Mutant'][y_col].dropna()

        if len(wt) > 1 and len(mut) > 1:
            stat, p = ttest_ind(wt, mut, equal_var=False)
            sig = "***" if p < 0.001 else "**" if p < 0.01 else "*" if p < 0.05 else "ns"

            y_max = sub[y_col].max()
            # Dynamic height for sig stars
            y_pos = y_max + (data[y_col].max() * 0.05)

            ax.text(i, y_pos, f"{sig}\n(p={p:.3f})",
                    ha='center', va='bottom', fontsize=9, fontweight='bold')

    if show_legend:
        from matplotlib.lines import Line2D
        custom_lines = [
            Line2D([0], [0], color=colors['WT'], marker='D', linestyle='', label='WT (Mean±SE)'),
            Line2D([0], [0], color=colors['Mutant'], marker='D', linestyle='', label='Mutant (Mean±SE)'),
        ]
        ax.legend(handles=custom_lines, loc='upper right', frameon=False)


def make_nice_plots(net_df, puncta_df, output_dir):
    sns.set_theme(style="ticks", context="talk")

    plots_to_make = []

    if net_df is not None:
        plots_to_make.append((net_df, 'total_ridge_pixels_norm', 'Network Amount (Norm)', 'Fold Change vs WT'))

    if puncta_df is not None:
        plots_to_make.append((puncta_df, 'puncta_count_norm', 'Puncta Count (Norm)', 'Fold Change vs WT'))
        plots_to_make.append((puncta_df, 'puncta_density_norm', 'Puncta Density (Norm)', 'Fold Change vs WT'))

    if not plots_to_make: return

    # --- PLOT 1: Per-Trial Breakdown ---
    fig1, axes = plt.subplots(len(plots_to_make), 1, figsize=(10, 5 * len(plots_to_make)), sharex=True)
    if len(plots_to_make) == 1: axes = [axes]

    for i, (data, col, title, ylab) in enumerate(plots_to_make):
        plot_metric_polished(axes[i], data, 'comparison_group', col, 'genotype', title, ylab, show_legend=(i == 0))

    plt.setp(axes[-1].get_xticklabels(), rotation=25, ha='right')
    sns.despine(offset=10, trim=True)
    plt.tight_layout()
    plt.savefig(output_dir / "PLOT_Per_Trial_Breakdown.png", dpi=300)

    # --- PLOT 2: Global Summary ---
    fig2, axes2 = plt.subplots(1, len(plots_to_make), figsize=(5 * len(plots_to_make), 6))
    if len(plots_to_make) == 1: axes2 = [axes2]

    for i, (data, col, title, ylab) in enumerate(plots_to_make):
        plot_metric_polished(axes2[i], data, 'genotype', col, 'genotype', f"Global: {title}", ylab, show_legend=False)

    sns.despine(offset=10, trim=True)
    plt.tight_layout()
    plt.savefig(output_dir / "PLOT_Global_Summary.png", dpi=300)
    print("[OK] Plots saved.")


def save_network_montages(data_dir: Path):
    """Create per-trial network QC montages from saved overlays."""
    print(f"\n--- Generating Network Montages ---")
    for trial_dir in sorted([p for p in data_dir.iterdir() if p.is_dir()]):
        overlay_dir = trial_dir / "analysis_results" / "network_overlay"
        if not overlay_dir.exists(): continue

        images = sorted(overlay_dir.glob("*.tif"))
        if not images: continue

        output_path = overlay_dir.parent / "Montage_Network.png"
        try:
            pipeline_utils.create_montage(trial_dir.name, images, output_path)
            print(f"[OK] {trial_dir.name}: {output_path.name}")
        except Exception:
            print(f"[!] Failed montage for {trial_dir.name}")
            traceback.print_exc()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("data_dir", type=Path)
    parser.add_argument("--run-id", type=str, default=None,
                        help="Run identifier (default: timestamp)")
    parser.add_argument("--cores", type=str, default=None,
                        help="Cores used (for metadata only)")
    args = parser.parse_args()

    run_id = args.run_id or datetime.now().strftime("%Y%m%d_%H%M%S")
    run_dir = args.data_dir / "analysis_runs" / run_id
    run_dir.mkdir(parents=True, exist_ok=True)

    print(f"--- Compiling Stats & Cleaning Data in {args.data_dir} ---")
    print(f"[i] Run ID: {run_id}")

    net_df = load_and_combine(args.data_dir, "network_stats.csv")
    puncta_df = load_and_combine(args.data_dir, "puncta_stats.csv")

    # Clean & Save CSVs
    if net_df is not None:
        normalize_by_group(net_df, 'total_ridge_pixels')
        normalize_by_group(net_df, 'structure_index')
        if 'tissue_area' in net_df.columns:
            net_df.rename(columns={'tissue_area': 'mask_area'}, inplace=True)

        # Save
        cols_to_save = [c for c in net_df.columns if 'Unnamed' not in c]
        net_df[cols_to_save].to_csv(run_dir / "DATA_Networks_Combined.csv", index=False)

    if puncta_df is not None:
        normalize_by_group(puncta_df, 'puncta_count')
        normalize_by_group(puncta_df, 'puncta_density')
        normalize_by_group(puncta_df, 'avg_nnd')
        if 'tissue_area' in puncta_df.columns:
            puncta_df.rename(columns={'tissue_area': 'mask_area'}, inplace=True)

        cols_to_save = [c for c in puncta_df.columns if 'Unnamed' not in c]
        puncta_df[cols_to_save].to_csv(run_dir / "DATA_Puncta_Combined.csv", index=False)

    # Save run metadata
    meta = {
        "run_id": run_id,
        "data_dir": str(args.data_dir.resolve()),
        "timestamp": datetime.now().isoformat(),
        "python_version": sys.version.split()[0],
        "platform": platform.platform(),
        "cores": args.cores,
        "network_params": {
            "ridge_widths": network_detector.RIDGE_WIDTHS,
            "ridge_low_contrast": network_detector.RIDGE_CONTRAST_LOW,
            "ridge_high_contrast": network_detector.RIDGE_CONTRAST_HIGH,
            "width_scale_factor": network_detector.WIDTH_SCALE_FACTOR,
        },
        "puncta_params": {
            "log_min_sigma": puncta_detector.LOG_MIN_SIGMA,
            "log_max_sigma": puncta_detector.LOG_MAX_SIGMA,
            "log_num_sigma": puncta_detector.LOG_NUM_SIGMA,
            "log_threshold": puncta_detector.LOG_THRESHOLD,
            "log_overlap": puncta_detector.LOG_OVERLAP,
            "edge_margin": puncta_detector.EDGE_MARGIN,
        },
    }
    with open(run_dir / "run_meta.json", "w") as f:
        json.dump(meta, f, indent=2)

    save_network_montages(args.data_dir)
    make_nice_plots(net_df, puncta_df, run_dir)


if __name__ == "__main__":
    main()
