import json

from rdkit import Chem

from balancer.schema import Reaction
from balancer.utils import get_unique_fragments, get_frag_svg


def accept_reaction(reaction: Reaction) -> bool:
    try:
        assert reaction.meta['ra_ge_p']
        assert reaction.meta['r_ge_p']
        assert reaction.meta['r_unmapped']
        assert not reaction.meta['p_unmapped']
        return True
    except AssertionError:
        return False


if __name__ == '__main__':
    from collections import defaultdict

    reactions = Reaction.load_reactions_from_csv("../data/Pistachio_2022Q4/reactions_by_class/1.1.3.csv")
    reactions = [r for r in reactions if accept_reaction(r)]
    reaction_class_name = reactions[0].meta['nextmove_class']
    reaction_class_id = reactions[0].meta['nextmove_class_identifier']

    frag_dict = dict()
    frag_dict: dict[str, Chem.Mol]
    frag_comb_to_reactions = defaultdict(set)
    frag_comb_to_reactions: dict[str, set[Reaction]]
    for r in reactions:
        reactants, unmapped_frags_list = r.get_unmapped_fragments('r')
        unmapped_frags_unique_dict = get_unique_fragments([f for fl in unmapped_frags_list for f in fl],
                                                          consider_terminal_neighbors=True)
        for k, v in unmapped_frags_unique_dict.items():
            frag_dict[k] = v
        frag_comb = " + ".join(sorted(unmapped_frags_unique_dict.keys()))
        frag_comb_to_reactions[frag_comb].add(r)
    frag_svg_dict = {k: get_frag_svg(v, show_terminal_neighbors=True) for k, v in frag_dict.items()}
    frag_comb_to_reaction_info = {k: [(vv.reaction_smiles, vv.identifier) for vv in v] for k, v in
                                  frag_comb_to_reactions.items()}
    results = {
        "reaction_class_name": reaction_class_name,
        "reaction_class_identifier": reaction_class_id,
        "frag_comb_to_reaction_info": frag_comb_to_reaction_info,
        "frag_svg_dict": frag_svg_dict,
    }

    with open(f"{reaction_class_id}.json", "w") as f:
        json.dump(results, f)
