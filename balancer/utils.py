from __future__ import annotations

from typing import Literal

from rdkit import Chem
from rdkit.Chem import Mol, Atom


class ReactionError(Exception):
    """ use this for errors raised related to `Reaction` class """
    pass


class FragmentError(Exception):
    pass


def get_neighbors(mol: Mol, atom: Atom | int, map_type: Literal['all', 'mapped_only', 'unmapped_only']) -> list[Atom]:
    if isinstance(atom, int):
        atom = mol.GetAtomWithIdx(atom)
    else:
        assert isinstance(atom, Atom)
    if map_type == 'mapped_only':
        return [n for n in atom.GetNeighbors() if n.GetAtomMapNum() > 0]
    elif map_type == 'unmapped_only':
        return [n for n in atom.GetNeighbors() if n.GetAtomMapNum() == 0]
    else:
        return list(atom.GetNeighbors())


def smiles_to_mol_list(chemicals_smiles: str) -> [Chem.Mol]:
    """ convert a SMILES string to a list of `Chem.Mol` objects """
    chemicals_mols = []
    for individual_molecule_smiles in sorted(chemicals_smiles.split(".")):
        mol = Chem.MolFromSmiles(individual_molecule_smiles)
        if mol is None:
            raise ReactionError(f"Failed to create molecule from smiles: {individual_molecule_smiles}")
        chemicals_mols.append(mol)
    return chemicals_mols


def has_unmapped_atom(m: Chem.Mol) -> bool:
    """ if this molecule has unmapped atom (atom map number == 0) """
    return any(a.GetAtomMapNum() == 0 for a in m.GetAtoms())
