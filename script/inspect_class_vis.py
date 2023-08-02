import json
import random

import dash_bootstrap_components as dbc
import pandas as pd
from dash import Dash, Input, Output, callback, dash_table
from rdkit.Chem.AllChem import ReactionFromSmarts
from rdkit.Chem.Draw import MolDraw2DSVG

MAX_REACTIONS_PER_FRAG = 10

with open("1.1.1.json", "r") as f:
    data = json.load(f)

reaction_class_name = data["reaction_class_name"]
reaction_class_identifier = data['reaction_class_identifier']
frag_comb_to_reaction_info = data["frag_comb_to_reaction_info"]

app = Dash(external_stylesheets=[dbc.themes.BOOTSTRAP])

records = []
for f, info_list in frag_comb_to_reaction_info.items():
    r = {"fragment set": f, "# of reactions": len(info_list)}
    records.append(r)
df = pd.DataFrame.from_records(records)


def smi2svg(reaction_smiles):
    rxn = ReactionFromSmarts(reaction_smiles, useSmiles=True)
    d2d = MolDraw2DSVG(1200, 300)
    # d2d = MolDraw2DSVG(-1, -1)
    # dopts = d2d.drawOptions()
    d2d.DrawReaction(rxn, highlightByReactant=True)
    # d2d.DrawReaction(rxn, highlightByReactant=False)
    d2d.FinishDrawing()
    svg = d2d.GetDrawingText().replace('svg:', '')
    return svg


table1 = dash_table.DataTable(df.to_dict('records'), [{"name": i, "id": i} for i in df.columns], id='tbl',
                              style_data={'whiteSpace': 'normal', 'height': 'auto'}, )

app.layout = dbc.Row(
    [
        dbc.Label('Inspection', className="col-12"),
        dbc.Col(
            dbc.Container([
                table1,
            ]), className="col-2"
        ),
        dbc.Col(id="tbl_out", className="col-9"),
    ]
)


@callback(Output('tbl_out', 'children'), Input('tbl', 'active_cell'))
def update_graphs(active_cell):
    if not active_cell:
        return "select a set of fragments"
    i_row = active_cell['row']
    i_column = active_cell['column']
    frag = df.iloc[i_row, 0]
    svgs = frag_comb_to_reaction_info[frag]
    if len(svgs) > MAX_REACTIONS_PER_FRAG:
        random.seed(42)
        svgs = random.sample(svgs, MAX_REACTIONS_PER_FRAG)

    columns_format = [
        dict(id='reaction', name='reaction', presentation='markdown'),
        dict(id='pistachio_id', name='pistachio_id'),
    ]
    records = []
    for smi, rid in svgs:
        record = {
            "pistachio_id": rid,
            "reaction": smi2svg(smi),
        }
        records.append(record)
    return dash_table.DataTable(
        markdown_options={
            'html': True
        },
        data=pd.DataFrame.from_records(records).to_dict("records"),
        columns=columns_format,
    )


if __name__ == "__main__":
    app.run(debug=False)
