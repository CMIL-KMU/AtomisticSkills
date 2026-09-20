"""
Drug Discovery Utilities

Core cheminformatics operations for drug discovery workflows.
Provides standalone functions for SMILES parsing, molecule standardization,
3D conformer generation, descriptor calculation, and fingerprinting.
"""

import logging
from typing import List, Tuple, Optional

from rdkit import Chem
from atomistic_analysis.molecules import mol_from_smiles as mol_from_smiles
from atomistic_analysis.molecules import compute_descriptors as compute_descriptors
from atomistic_analysis.molecules import compute_fingerprints as compute_fingerprints

logger = logging.getLogger(__name__)


# ==================== SMILES Parsing ====================


def parse_smiles_from_string(
    smiles: str, name: Optional[str] = None
) -> List[Tuple[str, Optional[str]]]:
    """Parse a single SMILES string into [(smiles, name)] format."""
    return [(smiles, name)]


def parse_smiles_from_file(filepath: str) -> List[Tuple[str, Optional[str]]]:
    """
    Parse SMILES from file. Supports multiple formats:
    - SMILES (no name)
    - SMILES\\tNAME (tab-separated)
    - SMILES NAME (space-separated)

    Lines starting with # are ignored.
    """
    records = []
    with open(filepath, "r") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue

            # Try tab-separated first, then space-separated
            if "\t" in line:
                parts = line.split("\t", 1)
            else:
                parts = line.split(None, 1)

            smiles = parts[0]
            name = parts[1].strip() if len(parts) > 1 else None
            records.append((smiles, name))

    return records


# ==================== Molecule Operations ====================




def standardize_mol(mol: Chem.Mol, mode: str = "cleanup") -> Optional[Chem.Mol]:
    """
    Apply RDKit MolStandardize transformations.

    Args:
        mol: RDKit molecule
        mode: Standardization mode
            - "none": No standardization
            - "cleanup": Basic cleanup (sanitization, charge parent)
            - "parent": Largest fragment + cleanup
            - "uncharged": Parent + neutralize
            - "tautomer": Parent + canonical tautomer

    Returns:
        Standardized molecule
    """
    if mode == "none":
        return mol

    try:
        from rdkit.Chem.MolStandardize import rdMolStandardize

        if mode in ["cleanup", "parent", "uncharged", "tautomer"]:
            # Always start with cleanup
            mol = rdMolStandardize.Cleanup(mol)

        if mode in ["parent", "uncharged", "tautomer"]:
            # Get largest fragment
            mol = rdMolStandardize.FragmentParent(mol)

        if mode == "uncharged":
            # Neutralize charges
            uncharger = rdMolStandardize.Uncharger()
            mol = uncharger.uncharge(mol)

        if mode == "tautomer":
            # Get canonical tautomer
            te = rdMolStandardize.TautomerEnumerator()
            mol = te.Canonicalize(mol)

        return mol

    except ImportError:
        logger.warning(
            "rdkit.Chem.MolStandardize not available, returning original molecule"
        )
        return mol
    except Exception as e:
        logger.error(f"Standardization failed: {e}")
        return None


def enumerate_protomers(
    smiles: str, ph_min: float = 6.4, ph_max: float = 8.4, max_variants: int = 10
) -> List[str]:
    """
    Enumerate pH-dependent protonation states using Dimorphite-DL.

    Returns:
        List of SMILES strings for different protomers
    """
    try:
        from dimorphite_dl import DimorphiteDL

        dimorphite = DimorphiteDL(
            min_ph=ph_min, max_ph=ph_max, max_variants=max_variants, silent=True
        )
        protomers = dimorphite.protonate(smiles)
        return protomers if protomers else [smiles]
    except ImportError:
        logger.warning("dimorphite_dl not available, returning original SMILES")
        return [smiles]
    except Exception as e:
        logger.error(f"Protomer enumeration failed: {e}")
        return [smiles]


# ==================== 3D Conformer Generation ====================


def embed_conformers(mol: Chem.Mol, num_confs: int = 10) -> List[int]:
    """
    Generate 3D conformers using ETKDG.

    Returns:
        List of conformer IDs
    """
    from rdkit.Chem import AllChem

    try:
        conf_ids = AllChem.EmbedMultipleConfs(
            mol, numConfs=num_confs, params=AllChem.ETKDGv3()
        )
        return list(conf_ids)
    except Exception as e:
        logger.error(f"Conformer embedding failed: {e}")
        return []


def optimize_conformers(
    mol: Chem.Mol, max_iters: int = 500, force_field: str = "MMFF94"
) -> List[Tuple[int, float, bool]]:
    """
    Optimize all conformers with specified force field.

    Returns:
        List of (conformer_id, energy, converged) tuples
    """
    from rdkit.Chem import AllChem

    results = []

    for conf_id in range(mol.GetNumConformers()):
        try:
            if force_field.upper().startswith("MMFF"):
                props = AllChem.MMFFGetMoleculeProperties(mol)
                ff = AllChem.MMFFGetMoleculeForceField(mol, props, confId=conf_id)
                converged = ff.Minimize(maxIts=max_iters) == 0
                energy = ff.CalcEnergy()
            else:  # UFF
                ff = AllChem.UFFGetMoleculeForceField(mol, confId=conf_id)
                converged = ff.Minimize(maxIts=max_iters) == 0
                energy = ff.CalcEnergy()

            results.append((conf_id, energy, converged))
        except Exception as e:
            logger.error(f"Optimization failed for conformer {conf_id}: {e}")
            results.append((conf_id, float("inf"), False))

    return results


def select_best_conformer(
    opt_results: List[Tuple[int, float, bool]],
) -> Tuple[int, float, bool]:
    """Select conformer with lowest energy."""
    if not opt_results:
        return (-1, float("inf"), False)
    return min(opt_results, key=lambda x: x[1])


# ==================== PDBQT Conversion ====================


def mol_to_pdbqt(mol: Chem.Mol, output_path: str, conf_id: int = -1) -> None:
    """Convert RDKit molecule to PDBQT format using Meeko."""
    from meeko import MoleculePreparation, PDBQTWriterLegacy

    preparator = MoleculePreparation()
    mol_setups = preparator.prepare(mol)

    for setup in mol_setups:
        pdbqt_string = PDBQTWriterLegacy.write_string(setup)[0]
        with open(output_path, "w") as f:
            f.write(pdbqt_string)
        break  # Only write first setup


# ==================== Molecular Descriptors ====================




# ==================== Fingerprints ====================
