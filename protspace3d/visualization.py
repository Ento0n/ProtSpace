#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import pandas as pd

SYMBOLS = [
    "circle",
    "square",
    "diamond",
    "cross",
    "x",
    "circle-open",
    "square-open",
    "diamond-open"
]


# https://github.com/sacdallago/bio_embeddings/blob/develop/bio_embeddings/visualize/plotly_plots.py
def render(df: pd.DataFrame, selected_column: str):
    import plotly.graph_objects as go

    col_groups = df[selected_column].unique().tolist()

    df["class_index"] = np.ones(len(df)) * -100

<<<<<<< HEAD
    data = []
    for group_idx, group_value in enumerate(col_groups):
        df_group = df[df[selected_column] == group_value]
        trace = go.Scatter3d(
            x=df_group["x"],
            y=df_group["y"],
            z=df_group["z"],
            mode="markers",
            name=group_value,
            # TODO: figure something out to deal with the colors
            # 10 colors are available; once those are used, pick different symbol
            marker=dict(symbol=SYMBOLS[group_idx % 8]),
=======
        app.layout = dbc.Container(
            [
                # Header
                dbc.Row(
                    [
                        dbc.Col(
                            html.H1("ProtSpace3D", style={"color": "white"}),
                            width=4,
                            style={"background-color": "black"},
                        ),
                        dbc.Col(width=7, style={"background-color": "black"}),
                        dbc.Col(
                            html.Img(
                                src=app.get_asset_url("logo.png"),
                                alt="Rostlab-logo",
                                style={"height": "60px", "width": "60px"},
                            ),
                            style={"background-color": "black"},
                            width=1,
                        ),
                    ]
                ),
                # space between header and content below
                dbc.Row([html.Br()]),
                # graph and controls
                dbc.Row(
                    [
                        dbc.Col(
                            [
                                dcc.Dropdown(
                                    self.csv_header,
                                    self.csv_header[0],
                                    id="dd_menu",
                                    searchable=False,
                                    clearable=False,
                                ),
                                dcc.Graph(
                                    id="graph",
                                    figure=self.fig,
                                    clear_on_unhover=True,
                                    style={
                                        "width": "90%",
                                        "height": "90vh",
                                    },
                                ),
                            ],
                            width=8,
                        ),
                        dbc.Col(
                            [
                                dashbio.NglMoleculeViewer(id="ngl_molecule_viewer"),
                            ],
                            width=4,
                        ),
                    ]
                ),
            ],
            fluid=True,
>>>>>>> 3dd30fe (Working ngl molecule viewer)
        )
        data.append(trace)
        df.loc[df[selected_column] == group_value, "class_index"] = group_idx

    fig = go.Figure(data=data)
    fig.update_layout(
        # Remove axes ticks and labels as they are usually not informative
        scene=dict(
            xaxis=dict(showticklabels=False, showspikes=False, title=""),
            yaxis=dict(showticklabels=False, showspikes=False, title=""),
            zaxis=dict(showticklabels=False, showspikes=False, title=""),
        ),
    )
    return fig
