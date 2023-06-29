from balancer.schema import Reaction

"""
inspect the reaction smiles from `pistachio.smi`, write out results to a `csv` file

1. you may encounter SMILES strings that violate rdkit SMILES grammar (valence constraints most likely), 
    capture such errors in `Reaction.meta`, they will be included in `Reaction.inspect`
2. some reaction SMILES have cxsmiles (ChemAxon) annotations, some do not
3. `Chem.Mol` objects will not be created in `Reaction.reactants/agents/products` until you run `Reaction.get_molecules`
4. use `pandas` for creating `csv` file

some dummy comments
"""
