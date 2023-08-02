import os
import pickle
import sys

import more_itertools as mit
import pandas as pd
from loguru import logger
from pandas._typing import FilePath
from rdkit import RDLogger
from tqdm import tqdm

from balancer.schema import Reaction, ReactionError

"""
inspect the reaction smiles from `pistachio.smi`, write out results to a `csv` file
"""

RDLogger.DisableLog('rdApp.*')  # https://github.com/rdkit/rdkit/issues/2683#issuecomment-538872880

PISTACHIO_VERSION = "2022Q4"


class PistachioError(Exception):
    pass


def load_pistachio_mapped_smi(smi_file: FilePath, version: str, batch_size=10000, wdir: FilePath = "./",
                              dump_reactions=False):
    os.makedirs(wdir, exist_ok=True)
    logfile = f"{wdir}/load_pistachio_mapped_smi-{version}.log"
    logger_id = logger.add(logfile)
    logger.remove(0)
    logger.info(f"loading: {smi_file}")
    with open(smi_file, "r") as f:
        lines = f.readlines()
    logger.info(f"loaded # of lines: {len(lines)}")

    i_line = 0
    i_chunk = 0
    for line_chunk in tqdm(list(mit.chunked(lines, batch_size))):
        pistachio_reactions = []
        reaction_dicts = []
        for line in line_chunk:
            items = line.strip().split("\t")
            reaction_smiles, document_identifier, compound_identifier, nextmove_class_identifier, nextmove_class = items
            if reaction_smiles.endswith("|") and reaction_smiles.count("|") >= 2:
                reaction_smiles = reaction_smiles.split()[0]
            reaction_identifier = f"{version}--{i_line}"
            r = Reaction.from_reaction_smiles(
                reaction_smiles, reaction_identifier,
                meta={
                    "document_identifier": document_identifier,
                    "compound_identifier": compound_identifier,
                    "nextmove_class_identifier": nextmove_class_identifier,
                    "nextmove_class": nextmove_class,
                }
            )
            try:
                r.get_molecules()
                r.inspect()
            except ReactionError as e:
                logger.warning(
                    f"Error in parsing molecules for reaction {reaction_identifier}: {reaction_smiles}\n{str(e)}")
            reaction_dicts.append(r.as_dict())
            pistachio_reactions.append(r)
            i_line += 1
        i_chunk += 1
        df_meta = pd.DataFrame.from_records(reaction_dicts)
        df_meta.to_csv(f"{wdir}/meta-{i_chunk:06}.csv", index=False)
        if dump_reactions:
            with open(f"{wdir}/reactions-{i_chunk:06}.pkl", "wb") as f:
                pickle.dump(pistachio_reactions, f)
    logger.remove(logger_id)
    logger.add(sys.stdout)


def group_reactions():
    """
    need to first combine csvs
    `head -n 1 meta-000001.csv > combined.out && tail -n+2 -q meta-*.csv >> combined.out`
    """
    df = pd.read_csv("../data/Pistachio_2022Q4/reactions/combined.csv")
    df = df.drop_duplicates(subset=['reaction_smiles'])
    for class_id, df_class in tqdm(df.groupby("meta__nextmove_class_identifier")):
        df_class.to_csv(f"../data/Pistachio_2022Q4/reactions_by_class/{class_id}.csv", index=False)


if __name__ == '__main__':
    load_pistachio_mapped_smi("../data/Pistachio_2022Q4/pistachio.smi", PISTACHIO_VERSION,
                              wdir="../data/Pistachio_2022Q4/reactions")
    group_reactions()
