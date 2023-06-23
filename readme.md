# Project title: 
Byproduct inference and reaction balancer

# Motivation:
Organic reaction equations stored in databases are often
1. unbalanced;
2. incomplete with missing/omitted byproducts. 

Such equations 
1. cannot support a bijective atom mapping between reactants and products;
2. cannot be used for thermochemistry calculations.

[//]: # (why does Alex think this is important?)

# Formal description
```
Given 
1. a complete set of reactant molecular graphs {R_i}, and
2. an incomplete set of products {P_j}, 

find
1. nonnegative stoichiometric coefficients {x_i}, {x_j} and {x_k}, and
2. a set of *reasonable* byproduct molecular graphs {B_k}

such that the following is balanced:
\sum_i x_iR_i = \sum_j x_jP_j + \sum_k x_kB_k.
```

- what to do if we have multiple solutions? e.g. `HCOOH` and `HH.COO`.
  - maybe reaction conditions should also be included
- do we count hydrogen?
- what is *reasonable*?
    - equivalent to what actually happened
    - can be found in databases, e.g. pubchem
    - heuristics
        - valence constraints
        - num of bond rearrangement is small
        - energies
        - mechanism?

# Notes
1. Reaction SMILES
   1. While in [theory](https://www.daylight.com/meetings/summerschool01/course/basics/smirks.html), the agents (things between `>`) of 
   a reaction SMILES cannot contribute to its products (things after the 2nd `>`),
   in practice (Pistachio, CAS, and Reaxys) this is allowed.
   2. Quick visualization of reaction SMILES is available through [ASKCOS](https://askcos.mit.edu/).
2. Atom mapping
   1. Atom mappings in USPTO came from `Indigo` which is not very reliable. 
   2. Atom mappings in Pistachio came from `NameRxn` which is more reliable when the reaction class is not `Unrecognized`.
   3. A set of atom mappers are wrapped in [ASKCOS](https://gitlab.com/mlpds_mit/ASKCOS/atom-mapping-services/-/tree/main/)
3. Side product vs byproduct
    ```
    intended reaction:
    A + B + C + D > E > AB + CD
    AB: desired product
    CD: byproduct
    
    side reaction
    A + C > AC
    AC: side product
    ```

4. Workflow
   - reaction smi -> atom mapped reaction smi -> unmapped atoms -> byproducts
5. Reaction classification
   - [askcos](https://gitlab.com/mlpds_mit/ASKCOS/askcos-core/-/tree/dev/askcos/synthetic/reaction_classification)
    and [api](https://askcos-demo.mit.edu/api/v2/?format=api)