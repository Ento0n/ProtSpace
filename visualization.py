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
    "diamond-open",
]


# https://github.com/sacdallago/bio_embeddings/blob/develop/bio_embeddings/visualize/plotly_plots.py
def render(df: pd.DataFrame, selected_column: str):
    import plotly.graph_objects as go

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
            marker=dict(symbol=SYMBOLS[group_idx % 10]),
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
