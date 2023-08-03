from __future__ import annotations

from collections import defaultdict

from rdkit import Chem
from rdkit.Chem import Mol, Atom, RWMol


class ReactionError(Exception):
    """ use this for errors raised related to `Reaction` class """
    pass


def get_neighbors(mol: Mol, atom: Atom | int, mapped=True) -> list[Atom]:
    if isinstance(atom, int):
        atom = mol.GetAtomWithIdx(atom)
    else:
        assert isinstance(atom, Atom)
    if mapped:
        return [n for n in atom.GetNeighbors() if n.GetAtomMapNum() > 0]
    else:
        return [n for n in atom.GetNeighbors() if n.GetAtomMapNum() == 0]


def get_unmapped_fragments(mol: Mol, use_original_mol=False, auto_add_dummies=False) -> list[Mol]:
    # identifier the bonds that connect a mapped atom with an unmapped atom
    if not use_original_mol:
        mol = Mol(mol)
    bonds_of_interest = []

    for a in mol.GetAtoms():
        a.SetProp("atom_type", "non-terminal")
        a.SetProp("original_idx", str(a.GetIdx()))

    for bond in mol.GetBonds():
        a1 = bond.GetBeginAtom()
        a2 = bond.GetEndAtom()
        a1_mn = a1.GetAtomMapNum()
        a2_mn = a2.GetAtomMapNum()
        if a1_mn * a2_mn == 0 and a1_mn != a2_mn:
            if a1_mn == 0:
                terminal_atom = a1
            else:
                terminal_atom = a2
            # a terminal atom is an unmapped atom has at least one mapped neighboring atom
            terminal_atom.SetProp("atom_type", "terminal")
            bonds_of_interest.append(bond)

    if auto_add_dummies:
        mol_f = Chem.FragmentOnBonds(mol, (b.GetIdx() for b in bonds_of_interest), addDummies=True)
        fragments = Chem.GetMolFrags(mol_f, asMols=True)
        unmapped_fragments = [f for f in fragments if all(a.GetAtomMapNum() == 0 for a in f.GetAtoms())]
        return unmapped_fragments
    else:
        mol_f = Chem.FragmentOnBonds(mol, (b.GetIdx() for b in bonds_of_interest), addDummies=False)
        fragments = Chem.GetMolFrags(mol_f, asMols=True)
        unmapped_fragments = [f for f in fragments if all(a.GetAtomMapNum() == 0 for a in f.GetAtoms())]

        unmapped_fragments_with_dummies = []
        for f in unmapped_fragments:
            fwd = RWMol(f)
            for a in f.GetAtoms():
                if a.GetProp("atom_type") == "terminal":
                    original_idx = int(a.GetProp("original_idx"))
                    mapped_neighbor_atoms = get_neighbors(mol, original_idx, mapped=True)
                    for mna in mapped_neighbor_atoms:
                        dummy = Atom("*")
                        dummy.SetProp("symbol_for_dummy", mna.GetSymbol())
                        dummy.SetProp("original_idx", str(mna.GetIdx()))
                        dummy_index = fwd.AddAtom(dummy)
                        dummy_bond_type = mol.GetBondBetweenAtoms(original_idx, mna.GetIdx()).GetBondType()
                        fwd.AddBond(a.GetIdx(), dummy_index, order=dummy_bond_type)
            fwd = fwd.GetMol()
            # Chem.SanitizeMol(fwd)  # TODO this should remove excess implicit hydrogen atoms
            unmapped_fragments_with_dummies.append(fwd)
        return unmapped_fragments_with_dummies


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


def get_explicit_fragment(frag: Mol) -> Mol:
    rwf = Chem.RWMol(frag)
    for a in frag.GetAtoms():
        if a.GetSymbol() == "*":
            actual_symbol = a.GetProp("symbol_for_dummy")
            new_atom = Atom(actual_symbol)
            rwf.ReplaceAtom(a.GetIdx(), new_atom, preserveProps=True)
    return rwf.GetMol()


def get_unique_fragments(fragments: list[Mol], consider_terminal_neighbors=False) -> dict[str, Mol | dict[str, Mol]]:
    if consider_terminal_neighbors:
        # d = defaultdict(dict)
        d = dict()
        for f in fragments:
            f_smi = Chem.MolToSmiles(f)
            exf_smi = Chem.MolToSmiles(get_explicit_fragment(f))
            d[f"{f_smi}|{exf_smi}|"] = f
        return d
    else:
        return {Chem.CanonSmiles(Chem.MolToSmiles(f)): f for f in fragments}


def get_frag_svg(frag: Mol, show_terminal_neighbors=False):
    from rdkit.Chem.Draw import MolDraw2DSVG
    d2d = MolDraw2DSVG(-1, -1)
    if show_terminal_neighbors:
        exp_frag = get_explicit_fragment(frag)
        hl_aids = []
        for a in exp_frag.GetAtoms():
            try:
                assert len(a.GetProp("symbol_for_dummy")) > 0
                hl_aids.append(a.GetIdx())
            except (AssertionError, KeyError) as e:
                continue
        d2d.DrawMolecule(exp_frag, highlightAtoms=hl_aids)
    else:
        d2d.DrawMolecule(frag)
    d2d.FinishDrawing()
    text = d2d.GetDrawingText()
    return text
