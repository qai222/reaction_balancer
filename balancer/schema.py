from __future__ import annotations

import pandas as pd
from collections import defaultdict
from pandas._typing import FilePath
from pydantic import BaseModel
from rdkit.Chem import SanitizeFlags, SanitizeMol
from rdkit.Chem.AllChem import ReactionFromSmarts
from rdkit.Chem.Draw import MolDraw2DSVG

from balancer.utils import *


class Reaction:

    def __init__(
            self,
            identifier: str,
            reaction_smiles: str,
            reactants_smiles: str,
            agents_smiles: str,
            products_smiles: str,
            reactants: tuple[Chem.Mol] = (),
            agents: tuple[Chem.Mol] = (),
            products: tuple[Chem.Mol] = (),
            meta: dict[str, str | int | float] = None,
    ):
        """
        describing a reaction from pistachio
        
        :param identifier: the string identifier, e.g. 'US20010001799A1_0719'
        :param reaction_smiles: the reaction SMILES without ChemAxon annotations
        :param reactants_smiles: 
        :param agents_smiles: 
        :param products_smiles: 
        :param reactants: a list of `Chem.Mol` objects
        :param agents: a list of `Chem.Mol` objects
        :param products: a list of `Chem.Mol` objects
        :param meta: metadata about this reaction, e.g. {"line_number": 2, "reaction_class": "1.7.5	Hydroxy to triflyloxy"}
        """
        if meta is None:
            meta = dict()
        self.meta = meta
        self.identifier = identifier
        self.products = products
        self.agents = agents
        self.reactants = reactants
        self.products_smiles = products_smiles
        self.agents_smiles = agents_smiles
        self.reactants_smiles = reactants_smiles
        self.reaction_smiles = reaction_smiles
        # self.get_molecules()

    def __repr__(self):
        return f"{self.identifier}: {self.reaction_smiles}"

    def __hash__(self):
        return hash(self.reaction_smiles)

    def __eq__(self, other: Reaction):
        return self.reaction_smiles == other.reaction_smiles

    @property
    def has_molecules(self):
        return len(self.reactants) > 0 and len(self.products) > 0

    def get_elements(self, components: Literal['r', 'a', 'p']) -> set[str]:
        """ get the element set of reactants/agens/products """
        assert self.has_molecules
        if components == 'r':
            ms = self.reactants
        elif components == 'a':
            ms = self.agents
        else:
            ms = self.products
        elements = []
        for m in ms:
            elements += [a.GetSymbol() for a in m.GetAtoms()]
        return set(elements)

    def inspect(self):
        """
        update meta with a "tag" dictionary for the reaction
        
        use this to produce a csv file like
                                 line_number_in_pistachio,ra_ge_p,r_ge_p,r_unmapped,p_unmapped
        <reaction identifier>
        """
        d = {
            # if r+a >= p
            "ra_ge_p": self.get_elements('r').union(self.get_elements('a')).issuperset(self.get_elements('p')),
            # if r>=p
            "r_ge_p": self.get_elements('r').issuperset(self.get_elements('p')),
            # any unmapped atoms in r
            "r_unmapped": any(has_unmapped_atom(m) for m in self.reactants),
            # any unmapped atoms in p
            "p_unmapped": any(has_unmapped_atom(m) for m in self.products),
        }
        self.meta.update(d)

    def as_dict(self) -> dict:
        """ for serialization """
        d = {
            "identifier": self.identifier,
            "reaction_smiles": self.reaction_smiles,
            "reactants_smiles": self.reactants_smiles,
            "agents_smiles": self.agents_smiles,
            "products_smiles": self.products_smiles,
        }
        d_meta = {f"meta__{k}": v for k, v in self.meta.items()}
        d.update(d_meta)
        return d

    def get_molecules(self):
        """ run rdkit smiles to mol function """
        for role in ['reactants', 'agents', 'products']:
            mols = smiles_to_mol_list(getattr(self, role + "_smiles"))
            setattr(self, role, mols)

    @classmethod
    def from_reaction_smiles(cls, reaction_smiles: str, identifier: str = None, meta: dict[str, str] = None):
        """ construction method from a reaction SMILES """
        assert reaction_smiles.count(">") == 2, f"cannot find 2 `>` in the reaction smiles: {reaction_smiles}"

        reactants_smiles, agents_smiles, products_smiles = reaction_smiles.split(">")
        if identifier is None:
            identifier = reaction_smiles
        return cls(identifier, reaction_smiles, reactants_smiles, agents_smiles, products_smiles, meta=meta)

    def get_unmapped_fragments(self, where: Literal['r', 'p', 'a']) -> tuple[list[Mol], list[list[Fragment]]]:
        if not self.has_molecules:
            self.get_molecules()
        if where.startswith('r'):
            molecules = list(self.reactants)
        elif where.startswith('p'):
            molecules = list(self.products)
        else:
            molecules = list(self.reactants) + list(self.products)

        frag_lol = []
        for mol in molecules:
            frag_lol.append(Fragment.get_unmapped_fragments(mol))
        return molecules, frag_lol

    def get_unmapped_fragments_unique_flat(self, where: Literal['r', 'p', 'a'], exclude_whole=True) -> list[Fragment]:
        _, frag_lol = self.get_unmapped_fragments(where=where)
        unique = []
        for frag_lst in frag_lol:
            for frag in frag_lst:
                if frag not in unique:
                    if exclude_whole and "*" not in frag.smiles:
                        continue
                    else:
                        unique.append(frag)
        unique = sorted(unique, key=lambda x: x.as_tuple())
        return unique

    @staticmethod
    def load_reactions_from_csv(csv_file: FilePath) -> list[Reaction]:
        df = pd.read_csv(csv_file)
        reactions = []
        for record in df.to_dict(orient="records"):
            meta = dict()
            for k, v in record.items():
                if k.startswith("meta__"):
                    meta[k.lstrip("meta__")] = v
            reaction = Reaction.from_reaction_smiles(
                reaction_smiles=record['reaction_smiles'],
                identifier=record['identifier'],
                meta=meta
            )
            reactions.append(reaction)
        return reactions

    def get_reaction_svg(self, width=-1, height=-1) -> str:
        rxn = ReactionFromSmarts(self.reaction_smiles, useSmiles=True)
        d2d = MolDraw2DSVG(width, height)
        d2d.DrawReaction(rxn, highlightByReactant=True)
        d2d.FinishDrawing()
        text = d2d.GetDrawingText()
        return text


class Fragment(BaseModel):
    smiles: str
    """ canonical smiles with dummies """

    dummies_map_to: list[str]
    """ symbols that dummies in canonical smiles are mapped to """

    mol_json: str
    """ json dump of the fragment Mol object, which may be not canonical """

    atom_properties: dict[int, dict[str, str]] = dict()
    """ atomic properties not included in json dump """

    @property
    def identifier(self):
        return str(self.as_tuple())

    def as_tuple(self):
        return self.smiles, tuple(self.dummies_map_to)

    def as_dict(self):
        return {"smiles": self.smiles, "dummies_map_to": self.dummies_map_to}

    def __eq__(self, other: Fragment):
        return self.as_tuple() == other.as_tuple()

    def __hash__(self):
        return hash(self.as_tuple())

    def __repr__(self):
        return f"SMILES: {self.smiles} | DUMMIES: {self.dummies_map_to}"

    @staticmethod
    def get_unmapped_fragments(mol: Mol) -> list[Fragment]:
        """
        - Given a molecule, find the molecular fragments in which all atoms are unmapped.
        - A terminal atom of an unmapped molecular fragment is an atom has at least one mapped atom neighbor
        in the original molecule.
        - The mapped nearest neighbors of a terminal atom are noted with symbol `*` (`add_Dummies=True`)
        - The `*`s in the canonical smiles of a molecular fragment are assigned the following properties

        :param mol:
        :return:
        """
        # first creat a copy
        mol = Mol(mol)
        for a in mol.GetAtoms():
            a.SetProp("frag_tag", "non-terminal")
            a.SetProp("original_index", str(a.GetIdx()))

        # identifier the bonds that connect a mapped atom with an unmapped atom
        bonds_of_interest = []

        # dict[<terminal_atom_original_index>]->[(<BondType>, <nn_symbol>, <nn_idx>, <nn_atomic_number>), ]
        # ex. lookup_table[3] -> [(Chem.rdchem.BondType.SINGLE, "C", 2, 6), ]
        lookup_table = defaultdict(list)
        for bond in mol.GetBonds():
            a1 = bond.GetBeginAtom()
            a2 = bond.GetEndAtom()
            a1: Atom
            a2: Atom
            a1_mn = a1.GetAtomMapNum()
            a2_mn = a2.GetAtomMapNum()
            if a1_mn * a2_mn == 0 and a1_mn != a2_mn:
                if a1_mn == 0:
                    terminal_atom = a1
                    terminal_nn = a2
                else:
                    terminal_atom = a2
                    terminal_nn = a1
                # a terminal atom is an unmapped atom has at least one mapped neighboring atom
                terminal_atom.SetProp("frag_tag", "terminal")
                bonds_of_interest.append(bond)
                lookup_table[terminal_atom.GetIdx()].append(
                    [bond.GetBondType(), terminal_nn.GetSymbol(), terminal_nn.GetIdx(), terminal_nn.GetAtomicNum()])
                # [bond.GetBondType()].append((terminal_nn.GetSymbol(), terminal_nn.GetIdx()))

        mol_f = Chem.FragmentOnBonds(mol, (b.GetIdx() for b in bonds_of_interest), addDummies=True)
        fragments = Chem.GetMolFrags(mol_f, asMols=True, sanitizeFrags=False)
        for fragment in fragments:
            # see https://github.com/rdkit/rdkit/issues/46 and https://github.com/rdkit/rdkit/discussions/3901
            try:
                SanitizeMol(fragment, sanitizeOps=SanitizeFlags.SANITIZE_ALL ^ SanitizeFlags.SANITIZE_KEKULIZE)
            except Exception as e:
                raise FragmentError(e.__str__())
        unmapped_fragments = [f for f in fragments if all(a.GetAtomMapNum() == 0 for a in f.GetAtoms())]

        for uf in unmapped_fragments:
            uf: Chem.Mol
            terminal_atoms = []
            for a in uf.GetAtoms():
                try:
                    frag_tag = a.GetProp('frag_tag')
                except KeyError:
                    continue
                if frag_tag == 'terminal':
                    terminal_atoms.append(a)

            for terminal_atom in terminal_atoms:
                dummies_and_bonds = [(a, uf.GetBondBetweenAtoms(a.GetIdx(), terminal_atom.GetIdx())) for a in
                                     terminal_atom.GetNeighbors() if a.GetSymbol() == "*"]
                assignments = lookup_table[int(terminal_atom.GetProp("original_index"))]
                n_assignments = len(assignments)
                assert n_assignments == len(dummies_and_bonds)

                already_assigned_dummies = set()
                for assign in assignments:
                    nn_bt, nn_sym, nn_idx, nn_atomic_number = assign
                    for i_dummy in range(n_assignments):
                        if i_dummy in already_assigned_dummies:
                            continue
                        dummy_atom, bond = dummies_and_bonds[i_dummy]
                        dummy_atom.SetIsotope(0)
                        if bond.GetBondType() == nn_bt:
                            dummy_atom.SetProp("original_index", str(nn_idx))
                            dummy_atom.SetProp("original_symbol", nn_sym)
                            dummy_atom.SetProp("original_atomic_number", str(nn_atomic_number))
                            dummy_atom.SetProp("frag_tag", f"terminal-nn")
                            already_assigned_dummies.add(i_dummy)
                assert len(already_assigned_dummies) == n_assignments

        result = []
        for uf in unmapped_fragments:
            can_smi = Chem.MolToSmiles(uf, canonical=True)
            # get a tuple of atomic symbols for dummies in can_smi
            dummies_map_to = []
            for can_atom_order in Chem.CanonicalRankAtoms(uf):
                a = uf.GetAtomWithIdx(can_atom_order)
                if a.GetSymbol() == "*":
                    dummies_map_to.append(a.GetProp("original_symbol"))
            atom_properties = dict()
            for a in uf.GetAtoms():
                atom_properties[a.GetIdx()] = a.GetPropsAsDict(includePrivate=False, includeComputed=False)
            frag = Fragment(
                smiles=can_smi,
                dummies_map_to=dummies_map_to,
                mol_json=Chem.MolToJSON(uf),
                atom_properties=atom_properties
            )
            result.append(frag)
        return result

    def get_mol(self):
        frag_mol = Chem.JSONToMols(self.mol_json)[0]
        frag_mol: Mol
        for aid, props in self.atom_properties.items():
            a = frag_mol.GetAtomWithIdx(aid)
            a: Atom
            for k, v in props.items():
                a.SetProp(k, v)
        return frag_mol

    def get_svg(self, width=-1, height=-1, show_terminal_nn=True):
        frag_mol = self.get_mol()
        d2d = MolDraw2DSVG(width, height)
        if show_terminal_nn:
            for a in frag_mol.GetAtoms():
                if a.GetSymbol() == "*":
                    a: Atom
                    a.SetIsotope(int(a.GetProp("original_atomic_number")))
        d2d.DrawMolecule(frag_mol)
        d2d.FinishDrawing()
        text = d2d.GetDrawingText()
        return text
