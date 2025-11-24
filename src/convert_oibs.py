import argparse
from pathlib import Path
import numpy as np
import oiffile
import tifffile
from tqdm import tqdm


def convert_single_oib(image_path: Path):
    # Use context manager to safely handle file opening
    with oiffile.OifFile(image_path) as oib:
        image = oib.asarray()

    # Process first 2 channels
    for i, chan in enumerate(image[:2]):
        chan_name = f"chan{i}"

        # OUTPUT: Trial_Folder / max_z_projections / chanX / Image.tif
        save_dir = image_path.parent / 'max_z_projections' / chan_name
        save_dir.mkdir(parents=True, exist_ok=True)

        save_path = save_dir / f'MAX_{chan_name}_{image_path.stem}.tif'

        # Create Max Projection
        max_z = np.max(chan, axis=0)
        tifffile.imwrite(save_path, max_z)


def main():
    parser = argparse.ArgumentParser(description="Step 1: Convert OIBs to Max-Z TIFs")
    parser.add_argument('dir', type=Path, help="Root folder containing .oib files")
    args = parser.parse_args()

    oib_files = list(args.dir.rglob('*.oib'))

    if not oib_files:
        print(f"No .oib files found in {args.dir}")
        return

    print(f"Found {len(oib_files)} OIB files. Converting...")

    for image_path in tqdm(oib_files, desc="Converting"):
        try:
            convert_single_oib(image_path)
        except Exception as e:
            print(f"Error converting {image_path.name}: {e}")

    print("\n[!] Conversion Complete.")
    print("[!] Now please generate your masks in the 'chan1_masks' folders before running analysis.")


if __name__ == "__main__":
    main()