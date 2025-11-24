import subprocess
import sys
import argparse
from pathlib import Path


def validate_inputs(data_dir: Path) -> bool:
    print(f"--- Checking inputs in: {data_dir} ---")
    images = list(data_dir.rglob("chan1/*.tif"))
    masks = list(data_dir.rglob("chan1_masks/*.tif"))

    if not images or not masks:
        print("[!] ERROR: Missing images or masks.")
        return False

    print(f"[✓] Found {len(images)} images and {len(masks)} masks.")
    return True


def run_step(script_name, args):
    script_path = Path("src") / script_name
    print(f"\n>>> Running {script_name}...")
    cmd = [sys.executable, str(script_path)] + args
    try:
        subprocess.run(cmd, check=True)
    except subprocess.CalledProcessError:
        print(f"\n[!] Pipeline failed at {script_name}")
        sys.exit(1)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("dir", type=Path)
    parser.add_argument("--cores", type=str, default="8")
    args = parser.parse_args()

    if not validate_inputs(args.dir): sys.exit(1)
    path_str = str(args.dir)

    run_step("puncta_detector.py", [path_str, "--cores", args.cores])
    run_step("network_detector.py", [path_str, "--cores", args.cores])
    run_step("compile_stats.py", [path_str])  # Now includes montage generation

    print(f"\n[✓] Pipeline Complete.")


if __name__ == "__main__":
    main()