"""Existing drug discovery utilities shared with the installed analysis package."""

import logging
from typing import List, Tuple, Optional, Dict, Any
from rdkit import Chem
from rdkit.Chem import rdMolAlign
from rdkit.Chem import Descriptors, Crippen, Lipinski, QED, rdMolDescriptors

logger = logging.getLogger(__name__)


def mol_from_smiles(smiles: str) -> Optional[Chem.Mol]:
    """Convert SMILES string to RDKit molecule."""
    try:
        mol = Chem.MolFromSmiles(smiles)
        return mol
    except Exception as e:
        logger.error(f"Failed to parse SMILES '{smiles}': {e}")
        return None


def compute_descriptors(
    smiles: str, name: Optional[str] = None, include_sandp_tpsa: bool = False
) -> Dict[str, Any]:
    """
    Compute molecular descriptors and drug-likeness heuristics.

    Returns:
        Dictionary with descriptors and validity flags
    """
    mol = mol_from_smiles(smiles)
    if mol is None:
        return {
            "smiles": smiles,
            "name": name or smiles,
            "valid": False,
            "error": "Invalid SMILES",
        }

    try:
        # Descriptors
        mw = round(Descriptors.MolWt(mol), 2)
        logp = round(Crippen.MolLogP(mol), 2)
        tpsa = round(Descriptors.TPSA(mol, includeSandP=include_sandp_tpsa), 2)
        hbd = rdMolDescriptors.CalcNumHBD(mol)
        # Two HBA definitions, both reported explicitly. `Lipinski.NumHAcceptors` is an
        # alias whose meaning CHANGED between rdkit 2025.09.4 and 2025.09.6 (caffeine:
        # 6 -> 3), so it is never called here.
        #   hba          - strict SMARTS acceptor count: excludes amide and pyrrole-type
        #                  N whose lone pair is delocalised. Chemically correct count.
        #   hba_lipinski - raw N+O count. The crude surrogate Lipinski 1997 actually
        #                  specified for the Rule of Five, and what Ro5 is scored on.
        hba = rdMolDescriptors.CalcNumHBA(mol)
        hba_lipinski = rdMolDescriptors.CalcNumLipinskiHBA(mol)
        rotatable = Lipinski.NumRotatableBonds(mol)
        num_atoms = mol.GetNumAtoms()
        num_heavy = mol.GetNumHeavyAtoms()
        formal_charge = Chem.GetFormalCharge(mol)
        qed = round(QED.qed(mol), 2)

        # Lipinski Ro5 -- scored on the N+O surrogate, per Lipinski 1997
        lipinski_violations = sum([mw > 500, logp > 5, hbd > 5, hba_lipinski > 10])

        # Veber
        veber_pass = rotatable <= 10 and tpsa <= 140

        return {
            "smiles": Chem.MolToSmiles(mol),
            "name": name or Chem.MolToSmiles(mol),
            "valid": True,
            "molecular_weight": mw,
            "logp": logp,
            "tpsa": tpsa,
            "hbd": hbd,
            "hba": hba,
            "hba_lipinski": hba_lipinski,
            "hba_definitions": {
                "hba": "rdMolDescriptors.CalcNumHBA -- strict SMARTS acceptor count",
                "hba_lipinski": "rdMolDescriptors.CalcNumLipinskiHBA -- N+O count, scores Ro5",
            },
            "rotatable_bonds": rotatable,
            "num_atoms": num_atoms,
            "num_heavy_atoms": num_heavy,
            "formal_charge": formal_charge,
            "qed": qed,
            "lipinski_violations": lipinski_violations,
            "lipinski_ro5_pass": lipinski_violations <= 1,
            "veber_pass": veber_pass,
        }
    except Exception as e:
        return {
            "smiles": smiles,
            "name": name or smiles,
            "valid": False,
            "error": str(e),
        }


def compute_fingerprints(
    records: List[Tuple[str, Optional[str]]],
    radius: int = 2,
    fp_size: int = 2048,
    use_chirality: bool = False,
    use_features: bool = False,
    compute_similarity: bool = True,
) -> Dict[str, Any]:
    """
    Compute Morgan fingerprints and optional similarity matrix.

    Returns:
        Dictionary with fingerprints, compounds info, and optional similarity matrix
    """
    from rdkit.Chem import rdFingerprintGenerator
    from rdkit import DataStructs

    # Create fingerprint generator
    if use_features:
        invgen = rdFingerprintGenerator.GetMorganFeatureAtomInvGen()
        fpgen = rdFingerprintGenerator.GetMorganGenerator(
            radius=radius,
            fpSize=fp_size,
            includeChirality=use_chirality,
            atomInvariantsGenerator=invgen,
        )
    else:
        fpgen = rdFingerprintGenerator.GetMorganGenerator(
            radius=radius, fpSize=fp_size, includeChirality=use_chirality
        )

    # Compute fingerprints
    compounds = []
    fps = []

    for smiles, name in records:
        mol = mol_from_smiles(smiles)
        if mol is None:
            compounds.append({"smiles": smiles, "name": name or smiles, "valid": False})
            fps.append(None)
            continue

        try:
            fp = fpgen.GetFingerprint(mol)
            canonical_smiles = Chem.MolToSmiles(mol)

            compounds.append(
                {
                    "smiles": canonical_smiles,
                    "name": name or canonical_smiles,
                    "valid": True,
                    "fingerprint": fp.ToBitString(),
                }
            )
            fps.append(fp)
        except Exception as e:
            compounds.append(
                {
                    "smiles": smiles,
                    "name": name or smiles,
                    "valid": False,
                    "error": str(e),
                }
            )
            fps.append(None)

    result = {
        "n_compounds": len(compounds),
        "n_valid": sum(1 for c in compounds if c.get("valid")),
        "compounds": compounds,
    }

    # Compute similarity matrix if requested
    if compute_similarity:
        n = len(fps)
        matrix = [[None] * n for _ in range(n)]

        valid_indices = [i for i, fp in enumerate(fps) if fp is not None]
        valid_fps = [fps[i] for i in valid_indices]

        for i, idx_i in enumerate(valid_indices):
            for j, idx_j in enumerate(valid_indices):
                if idx_i <= idx_j:
                    sim = DataStructs.TanimotoSimilarity(valid_fps[i], valid_fps[j])
                    matrix[idx_i][idx_j] = round(float(sim), 4)
                    matrix[idx_j][idx_i] = matrix[idx_i][idx_j]

        result["similarity_matrix"] = matrix

    return result


def check_molecule_identity(ref: Chem.Mol, probe: Chem.Mol) -> None:
    """
    Verify that the reference and probe molecules represent the same compound.

    Compares InChIKeys (first 14 characters, the connectivity-only block) to
    catch the silent-failure mode where a user passes a docked pose of
    compound A against a reference of compound B: CalcRMS will happily produce
    a plausible-looking number if the atom counts line up, but the number is
    meaningless. We compare only the connectivity block so that protonation
    or tautomer differences between the docked form and the crystal form do
    not trigger a false mismatch.

    Args:
        ref: Reference molecule (crystal ligand).
        probe: Docked pose molecule.

    Raises:
        ValueError: If the molecules appear to be different compounds.
    """
    try:
        ref_key = Chem.inchi.MolToInchiKey(ref)
        probe_key = Chem.inchi.MolToInchiKey(probe)
    except Exception as e:
        raise ValueError(f"Cannot compute InChIKey for identity check: {e}")

    ref_connectivity = ref_key.split("-")[0]
    probe_connectivity = probe_key.split("-")[0]

    if ref_connectivity != probe_connectivity:
        raise ValueError(
            f"Reference and docked pose appear to be different molecules. "
            f"Reference InChIKey: {ref_key}, "
            f"docked pose InChIKey: {probe_key}. "
            f"Check that --reference and --docked correspond to the same compound."
        )


def symmetry_corrected_rmsd(
    ref: Chem.Mol,
    probe: Chem.Mol,
) -> float:
    """
    Compute the minimum heavy-atom RMSD over all valid atom mappings.

    Uses RDKit's CalcRMS, which enumerates molecular automorphisms (symmetry
    mappings) and returns the minimum RMSD *without* performing rigid-body
    alignment between probe and reference. This is the correct comparison
    for docking validation: the docked pose and the crystal reference are
    expected to share the receptor's coordinate frame, so any alignment step
    would artificially deflate the RMSD and mask a failing protocol.

    `symmetrizeConjugatedTerminalGroups=True` (default since RDKit 2022.09)
    correctly handles carboxylates, nitro groups, amidinium groups, and other
    conjugated terminal groups where the oxygens or nitrogens are chemically
    equivalent but formally labelled differently. Passed explicitly here for
    provenance and version independence.

    Args:
        ref: Reference molecule with 3D coordinates.
        probe: Docked pose with 3D coordinates.

    Returns:
        Minimum in-place heavy-atom RMSD in Angstroms.
    """
    return float(
        rdMolAlign.CalcRMS(
            probe,
            ref,
            symmetrizeConjugatedTerminalGroups=True,
        )
    )
