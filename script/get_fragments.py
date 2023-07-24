"""
get unmapped molecular fragments from a reaction smiles, the fragments can come from either reactants or products
"""
from balancer.schema import Reaction
from rdkit import Chem

# def get_unmapped_fragment_from_molecule(mol: Chem.Mol):
#


def get_unmapped_fragments_from_reaction_smiles(
        reaction_meta:dict,
        reaction_identifier:str = "dummy",
        reaction_smiles: str = "C[CH:1](C)[c:2]1[cH:3][cH:4][c:5]([cH:6][cH:7]1)[C:8](=[O:9])[CH2:10][CH2:11][CH2:12]Cl.O=C1CCC(N1[Br:13])=O>ClCCl>[O:9]=[C:8]([c:5]1[cH:6][cH:7][c:2]([cH:3][cH:4]1)[CH2:1][Br:13])[CH:10]1[CH2:11][CH2:12]1",
        from_reactants: bool = True,
):

    from rdkit.Chem import Draw

    reaction = Reaction.from_reaction_smiles(reaction_smiles=reaction_smiles, identifier=reaction_identifier, meta=reaction_meta)
    reaction.get_molecules()
    if from_reactants:
        molecules = reaction.reactants
    else:
        molecules = reaction.products


    for i_mol, mol in enumerate(molecules):
        img = Draw.MolToImage(mol)
        img.save(f"mol-{i_mol}.png")

        # identifier the indices of unmapped atoms
        unmapped_atoms = []
        atoms = mol.GetAtoms()
        for atom in atoms:
            atom_map_number = atom.GetAtomMapNum()
            if atom_map_number == 0:
                unmapped_atoms.append(atom)
                # print("unmapped:", atom.GetSymbol(), atom.GetIdx())

        # identifier the bonds that connect a mapped atom with an unmapped atom
        unmapped_atom_indices = [a.GetIdx() for a in unmapped_atoms]
        mapped_atom_indices = [i for i in range(len(atoms)) if i not in unmapped_atom_indices]
        bonds_of_interest = []
        for atom in unmapped_atoms:
            neighbor_atoms = atom.GetNeighbors()
            neighbor_atoms_mapped = [a for a in neighbor_atoms if a.GetIdx() in mapped_atom_indices]
            for neighbor_atom_mapped in neighbor_atoms_mapped:
                # print("these are mapped neighbors of:", atom, atom.GetSymbol(), atom.GetIdx())
                # print(neighbor_atom_mapped.GetSymbol(), neighbor_atom_mapped.GetIdx())
                bond = mol.GetBondBetweenAtoms(atom.GetIdx(), neighbor_atom_mapped.GetIdx())
                bonds_of_interest.append(bond)

        # fragment the molecule based on the above bonds
        mol_f = Chem.FragmentOnBonds(mol, (b.GetIdx() for b in bonds_of_interest))
        fragments = Chem.GetMolFrags(mol_f, asMols=True)
        grid = Draw.MolsToGridImage(fragments, useSVG=True)
        with open(f'molfrags-{i_mol}.svg', 'w') as f:
            f.write(grid)
        # collect all fragments that only have unmapped atoms, there should be only two cases:
        # all atoms in this frag are mapped, or all atoms are unmapped
        # in the end return all unmapped fragments as rdkit.Mol

        # print(Chem.MolToSmiles(mol))
        # print("="*12)


get_unmapped_fragments_from_reaction_smiles({})