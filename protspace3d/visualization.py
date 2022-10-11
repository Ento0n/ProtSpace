#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import pandas as pd

import dash
from dash import Dash, dcc, html
from dash import Input, Output
from dash.exceptions import PreventUpdate

import plotly.graph_objects as go

SYMBOLS = [
    "circle",
    "square",
    "diamond",
    "cross",
    "x",
    "circle-open",
    "square-open",
    "diamond-open",
]

styles = {
    "pre": {
        "border": "thin lightgrey solid",
        "overflowX": "scroll",
        "overflowY": "scroll",
        "width": "100%",
        "height": "90vh",
    }
}


def init_app(fig: go.Figure, csv_header: list):

    app = Dash(__name__)

    app.layout = html.Div(
        [
            dcc.Markdown("""**3D scatter plot of embedding space.**"""),
            dcc.Dropdown(
                csv_header,
                csv_header[0],
                id="dd_menu",
                searchable=False,
                clearable=False,
            ),
            dcc.Graph(
                id="graph", figure=fig, clear_on_unhover=True, style=styles["pre"]
            ),
            dcc.Store(id="store_data", storage_type="memory"),
        ]
    )

    return app


# https://github.com/sacdallago/bio_embeddings/blob/develop/bio_embeddings/visualize/plotly_plots.py
def render(df: pd.DataFrame, selected_column: str):
    col_groups = df[selected_column].unique().tolist()

    df["class_index"] = np.ones(len(df)) * -100

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


@dash.callback(
    Output("graph", "figure"), Input("store_data", "data"), Input("dd_menu", "value")
)
def update_graph(df, xaxis_column_name):
    # Check whether an input is triggered
    ctx = dash.callback_context
    if not ctx.triggered:
        raise PreventUpdate

    fig = render(pd.DataFrame.from_records(df), selected_column=xaxis_column_name)
    return fig
