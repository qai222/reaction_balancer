from __future__ import annotations

from typing import Literal

from rdkit import Chem


class ReactionError(Exception):
    """ use this for errors raised related to `Reaction` class """
    pass


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
    """ if this molecue has unmapped atom (atom map number == 0) """
    return any(a.GetAtomMapNum() == 0 for a in m.GetAtoms())


class Reaction:

    def __init__(
            self,
            identifier: str,
            reaction_smiles: str,
            reactants_smiles: [str],
            agents_smiles: [str],
            products_smiles: [str],
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
