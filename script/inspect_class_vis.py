import pickle

import dash_bootstrap_components as dbc
import pandas as pd
from dash import Dash, Input, Output, callback, dash_table, html

from balancer.schema import Fragment, Reaction

_ROW_1_HEIGHT = "300px"
_REACTION_SVG_WIDTH = 1200
_REACTION_SVG_HEIGHT = 300

with open("inspect.pkl", "rb") as f:
    FRAGMENT_DICTIONARY, FRAGMENT_TUPLE_DICTIONARY, FRAG_TO_FRAG_TUPLES, REACTION_DICTIONARY = pickle.load(f)


def get_fragment_tuple_list(fragment_identifier: str):
    return FRAG_TO_FRAG_TUPLES[fragment_identifier]


def get_n_reactions(fragment_tuple: tuple[str, ...]):
    return FRAGMENT_TUPLE_DICTIONARY[fragment_tuple][0]


def get_sample_reactions(fragment_tuple: tuple[str, ...]):
    return [REACTION_DICTIONARY[i] for i in FRAGMENT_TUPLE_DICTIONARY[fragment_tuple][1]]


def get_fragment_df():
    records = []
    for fid, (frag, num) in FRAGMENT_DICTIONARY.items():
        frag: Fragment
        r = {
            "identifier": fid,
            "SMILES": frag.smiles,
            "dummies": ",".join(frag.dummies_map_to),
            "occurrence": num,
        }
        records.append(r)
    records = sorted(records, key=lambda x: x["occurrence"], reverse=True)
    df = pd.DataFrame.from_records(records)
    return df


FRAGMENT_DF = get_fragment_df()


def get_fragment_tuple_df(fragment_id: str):
    ft_list = get_fragment_tuple_list(fragment_id)
    records = []
    for ft in ft_list:
        r = {
            'tuple': " + ".join(ft),
            '# of reactions': get_n_reactions(ft),
        }
        records.append(r)
    records = sorted(records, key=lambda x: r['# of reactions'], reverse=True)
    df = pd.DataFrame.from_records(records)
    return df


app = Dash(external_stylesheets=[dbc.themes.BOOTSTRAP])
TABLE1 = dash_table.DataTable(
    FRAGMENT_DF.to_dict('records'), [{"name": i, "id": i} for i in FRAGMENT_DF.columns], id='tbl1',
    style_data={'whiteSpace': 'normal', 'height': 'auto'},
    style_table={'height': _ROW_1_HEIGHT, 'overflowY': 'auto'}
)

TABLE2 = dash_table.DataTable([], id='tbl2', style_data={'whiteSpace': 'normal', 'height': 'auto'},
                              style_table={'height': _ROW_1_HEIGHT, 'overflowY': 'auto'})

columns_format = [
    dict(id='svg', name='svg', presentation='markdown'),
    dict(id='id', name='reaction_id'),
    dict(id='smiles', name='smiles'),
    dict(id='class_id', name='class_id'),
    dict(id='class_name', name='class_name'),
]
TABLE3 = dash_table.DataTable([], id='tbl3', style_data={'whiteSpace': 'normal', 'height': 'auto'},
                              style_table={'height': '400px', 'overflowY': 'auto'}, page_size=1, columns=columns_format,

                              markdown_options={
                                  'html': True
                              },

                              )

ROW1 = dbc.Row(
    [
        dbc.Col(
            TABLE1,
            className="p-3",
            width=6,
        ),
        dbc.Col(
            TABLE2,
            className="p-3",
            width=6
        ),
    ], className="px-3")

ROW2 = dbc.Row(
    dbc.Col(
        TABLE3,
        className="px-6",
        width=12
    ),
)

app.layout = html.Div(
    [ROW1, ROW2]
)


@callback(Output('tbl2', 'data'), Input('tbl1', 'active_cell'))
def update_tbl2(active_cell):
    if not active_cell:
        return []
    i_row = active_cell['row']
    frag_id = FRAGMENT_DF.iloc[i_row, 0]

    df = get_fragment_tuple_df(frag_id)
    return df.to_dict('records')


@callback(Output('tbl3', 'data'), Input('tbl2', 'active_cell'), Input('tbl2', 'data'))
def update_reaction(active_cell, tbl2_data):
    if not active_cell:
        return []
    i_row = active_cell['row']
    tbl2_record = tbl2_data[i_row]
    frag_tuple_repr = tbl2_record['tuple']
    frag_tuple = tuple([*frag_tuple_repr.split(" + ")])
    reactions = get_sample_reactions(frag_tuple)
    assert get_n_reactions(frag_tuple) == tbl2_record["# of reactions"], f"{get_n_reactions(frag_tuple)}, {tbl2_record['# of reactions']}"
    records = []
    for r in reactions:
        r: Reaction
        record = {
            "id": r.identifier,
            "smiles": r.reaction_smiles,
            "class_id": r.meta['nextmove_class_identifier'],
            "class_name": r.meta['nextmove_class'],
            "svg": r.get_reaction_svg(width=_REACTION_SVG_WIDTH, height=_REACTION_SVG_HEIGHT)
        }
        records.append(record)
    return records


if __name__ == "__main__":
    app.run(debug=False)
