import glob
import random
from collections import defaultdict
from os.path import basename
from typing import Literal
from tqdm import tqdm
from loguru import logger

from balancer.schema import Reaction

MAX_REACTION_PER_FRAGMENT_TUPLE = 50

def get_line_count(fn):
    with open(fn, "rb") as f:
        num_lines = sum(1 for _ in f)
    return num_lines

REACTION_CLASS_CSVS = [
    fn for fn in glob.glob("../data/Pistachio_2022Q4/reactions_by_class/*.csv")
    if not basename(fn).startswith("0.0")
]
REACTION_CLASS_CSVS = sorted(REACTION_CLASS_CSVS, key=lambda x: get_line_count(x), reverse=True)


def accept_reaction(reaction: Reaction) -> bool:
    try:
        assert reaction.meta['ra_ge_p']
        assert reaction.meta['r_ge_p']
        assert reaction.meta['r_unmapped']
        assert not reaction.meta['p_unmapped']
        return True
    except AssertionError:
        return False



def inspect_reactions(
        reactions: list[Reaction],
        inspect_fragments_from: Literal["r", "p", "a"] = "r",
        exclude_whole_molecules_in_fragmentation=True,
        max_reactions_per_fragment_tuple=MAX_REACTION_PER_FRAGMENT_TUPLE,
):
    reactions = {r.identifier: r for r in reactions if accept_reaction(r)}
    logger.info(f"inspecting # of reactions: {len(reactions)}")

    fragment_dictionary = dict()  # id -> example fragment object, occurrence
    fragment_tuple_dictionary = defaultdict(list)  # frag tuple -> list[reaction_id]
    for reaction in tqdm(reactions.values()):
        unique_fragments = reaction.get_unmapped_fragments_unique_flat(
            where=inspect_fragments_from,
            exclude_whole=exclude_whole_molecules_in_fragmentation
        )
        for uf in unique_fragments:
            if uf.identifier not in fragment_dictionary:
                fragment_dictionary[uf.identifier] = [uf, 1]
            else:
                fragment_dictionary[uf.identifier][1] += 1
        fragment_tuple = tuple(sorted([ff.identifier for ff in unique_fragments]))
        fragment_tuple_dictionary[fragment_tuple].append(reaction.identifier)

    reaction_dictionary = dict()
    fragment_tuple_dictionary_slim = dict()
    random.seed(42)
    for fragment_tuple in fragment_tuple_dictionary:
        reaction_ids = fragment_tuple_dictionary[fragment_tuple]
        if max_reactions_per_fragment_tuple < len(reaction_ids):
            reaction_ids_slim = random.sample(reaction_ids, max_reactions_per_fragment_tuple)
        else:
            reaction_ids_slim = reaction_ids
        fragment_tuple_dictionary_slim[fragment_tuple] = len(reaction_ids), reaction_ids_slim
        for rid in reaction_ids_slim:
            reaction_dictionary[rid] = reactions[rid]

    frag_to_frag_tuples = defaultdict(list)
    for frag in fragment_dictionary:
        fts = [ft for ft in fragment_tuple_dictionary if frag in ft]
        frag_to_frag_tuples[frag] = fts
    return fragment_dictionary, fragment_tuple_dictionary_slim, frag_to_frag_tuples, reaction_dictionary


def inspect_all(n_class=5, dump_pkl_as=None, from_most_popular=True):
    """
    inspect molecular fragments

    :param n_class: how many reaction classes to go through
    :param dump_pkl_as:
    :param from_most_popular: start from most popular reaction class or randomly select `n_class` classes
    :return:
    """
    reactions = []
    if from_most_popular:
        csvs = REACTION_CLASS_CSVS
    else:
        random.seed(42)
        csvs = random.sample(REACTION_CLASS_CSVS, n_class)
    i = 0
    for csv in csvs:
        r_cls_id = basename(csv).rstrip(".csv")
        reactions += Reaction.load_reactions_from_csv(csv)
        logger.info(f"loading: {r_cls_id}, loaded classes: {i}")
        i += 1
        if i >= n_class:
            break
    return_data = inspect_reactions(reactions)
    if dump_pkl_as:
        with open(dump_pkl_as, "wb") as f:
            pickle.dump(return_data, f)
    return return_data


if __name__ == '__main__':
    import pickle

    # data = inspect_one("1.1.1")
    # with open("1.1.1.pkl", "wb") as f:
    #     pickle.dump(data, f)

    data = inspect_all(5, dump_pkl_as="inspect.pkl", from_most_popular=False)
