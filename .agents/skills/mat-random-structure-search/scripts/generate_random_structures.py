"""
Generate random crystal structures for a given composition (AIRSS-style).

Generates random crystal structures using pymatgen's Structure.from_spacegroup()
with randomized lattice parameters and Wyckoff positions. Structures are filtered
for minimum interatomic distances to avoid unphysical configurations.

This implements the core idea of Ab Initio Random Structure Searching (AIRSS)
by Pickard & Needs, but uses MLIP relaxation instead of DFT for efficiency.

Usage:
    python generate_random_structures.py --composition NaCl --num_structures 100
    python generate_random_structures.py --composition Li2ZrCl6 --num_structures 50 --spacegroups 12,14,62,166

Requirements:
    - Conda environment: base-agent
    - Required packages: pymatgen, numpy
"""

import sys
from pathlib import Path

if __package__ in (None, ""):
    sys.path.insert(0, str(Path(__file__).resolve().parents[4]))

import argparse
import json

import numpy as np
from pymatgen.core import Composition


# Common space groups for inorganic crystals, weighted by frequency in ICSD
COMMON_SPACEGROUPS = [
    # Cubic
    225,
    227,
    221,
    229,
    216,
    226,
    # Hexagonal
    194,
    186,
    167,
    166,
    193,
    164,
    # Trigonal
    148,
    155,
    160,
    161,
    # Tetragonal
    139,
    140,
    136,
    129,
    141,
    127,
    # Orthorhombic
    62,
    63,
    64,
    58,
    55,
    57,
    59,
    61,
    33,
    36,
    # Monoclinic
    14,
    12,
    15,
    13,
    11,
    # Triclinic
    2,
    1,
]


from src.utils.analysis.random_structures import estimate_volume_per_atom as estimate_volume_per_atom


from src.utils.analysis.random_structures import get_min_distance as get_min_distance


from src.utils.analysis.random_structures import check_min_distances as check_min_distances


from src.utils.analysis.random_structures import generate_random_structure as generate_random_structure


from src.utils.analysis.random_structures import _random_lattice_params as _random_lattice_params


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Generate random crystal structures (AIRSS-style)"
    )
    parser.add_argument(
        "--composition",
        required=True,
        help="Chemical composition (e.g., NaCl, Li2ZrCl6, SrTiO3)",
    )
    parser.add_argument(
        "--num_structures",
        type=int,
        default=100,
        help="Number of structures to generate (default: 100)",
    )
    parser.add_argument(
        "--spacegroups",
        type=str,
        default=None,
        help="Comma-separated space group numbers to use (default: common SGs)",
    )
    parser.add_argument(
        "--volume_min",
        type=float,
        default=0.6,
        help="Minimum volume scale factor (default: 0.6)",
    )
    parser.add_argument(
        "--volume_max",
        type=float,
        default=1.8,
        help="Maximum volume scale factor (default: 1.8)",
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=None,
        help="Random seed for reproducibility",
    )
    parser.add_argument(
        "--output_dir",
        required=True,
        help="Directory to save generated structures",
    )
    args = parser.parse_args()

    if args.seed is not None:
        np.random.seed(args.seed)

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    composition = Composition(args.composition)
    reduced = composition.reduced_composition
    print(f"Composition: {reduced.reduced_formula}")
    print(f"Num atoms per formula unit: {int(reduced.num_atoms)}")

    # Determine space groups to use
    if args.spacegroups:
        spacegroups = [int(sg) for sg in args.spacegroups.split(",")]
    else:
        spacegroups = COMMON_SPACEGROUPS

    print(f"Using {len(spacegroups)} space groups")
    print(f"Target: {args.num_structures} structures")

    # Generate structures
    generated = []
    attempts = 0
    max_total_attempts = args.num_structures * 20  # safety limit

    while len(generated) < args.num_structures and attempts < max_total_attempts:
        # Pick random space group
        sg = int(np.random.choice(spacegroups))
        # Random volume scale
        vol_scale = np.random.uniform(args.volume_min, args.volume_max)

        structure = generate_random_structure(
            reduced, sg, volume_scale=vol_scale, max_attempts=10
        )
        attempts += 1

        if structure is not None:
            generated.append(
                {
                    "structure": structure,
                    "spacegroup": sg,
                    "volume_scale": vol_scale,
                }
            )
            if len(generated) % 10 == 0:
                print(f"  Generated {len(generated)}/{args.num_structures}...")

    print(f"Successfully generated {len(generated)} structures ({attempts} attempts)")

    # Save structures
    manifest_entries = []
    for i, item in enumerate(generated):
        formula = item["structure"].composition.reduced_formula
        cif_name = f"{i:04d}_{formula}_sg{item['spacegroup']}.cif"
        cif_path = output_dir / cif_name
        item["structure"].to(filename=str(cif_path))

        manifest_entries.append(
            {
                "index": i,
                "formula": formula,
                "spacegroup": item["spacegroup"],
                "volume_scale": round(item["volume_scale"], 3),
                "num_atoms": len(item["structure"]),
                "volume": round(item["structure"].volume, 2),
                "cif_file": cif_name,
            }
        )

    # Save manifest
    manifest_path = output_dir / "generation_manifest.json"
    manifest = {
        "composition": reduced.reduced_formula,
        "num_generated": len(generated),
        "total_attempts": attempts,
        "spacegroups_used": list(set(e["spacegroup"] for e in manifest_entries)),
        "volume_range": [args.volume_min, args.volume_max],
        "seed": args.seed,
        "structures": manifest_entries,
    }
    with open(manifest_path, "w") as f:
        json.dump(manifest, f, indent=2)

    # Print summary
    print(f"\n{'='*60}")
    print(f"Random Structure Generation: {reduced.reduced_formula}")
    print(f"{'='*60}")
    print(f"  Generated:    {len(generated)}")
    print(f"  Attempts:     {attempts}")
    print(f"  Success rate: {len(generated)/max(attempts,1)*100:.1f}%")
    print(f"  Output:       {output_dir}")
    print(f"  Manifest:     {manifest_path}")

    # Space group distribution
    sg_counts: dict[int, int] = {}
    for e in manifest_entries:
        sg_counts[e["spacegroup"]] = sg_counts.get(e["spacegroup"], 0) + 1
    print("\n  Space group distribution:")
    for sg, count in sorted(sg_counts.items(), key=lambda x: -x[1])[:10]:
        print(f"    SG {sg:3d}: {count} structures")

    print("\n  Next step: Relax all structures with an MLIP, then rank by energy.")

    # Save input configs for reproducibility
    from src.utils.config_utils import save_skill_inputs

    save_skill_inputs(args, args.output_dir)


if __name__ == "__main__":
    main()
