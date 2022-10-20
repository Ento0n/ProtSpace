#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
from pandas import DataFrame

import dash
from dash import Dash, dcc, html
from dash import Input, Output
from dash.exceptions import PreventUpdate
import dash_bootstrap_components as dbc

import plotly.graph_objects as go


class Visualizator:
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

    def __init__(self, fig: go.Figure, csv_header: list[str]):
        self.fig = fig
        self.csv_header = csv_header

    def init_app(self):
        """
        Initializes app & Builds html layout for Dash
        :return: layout
        """
        app = Dash(__name__, external_stylesheets=[dbc.themes.BOOTSTRAP])

        app.layout = dbc.Container(
            [
                dbc.Row(
                    [
                        dbc.Col(html.H1("ProtSpace3D"), width=4),
                        dbc.Col(width=4),
                        dbc.Col(html.Img(src=app.get_asset_url("logo.png")), width=4),
                    ]
                ),
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
                                    # style=self.styles["pre"],
                                    style={"width": "90%", "height": "90vh"},
                                ),
                                dcc.Store(id="store_data", storage_type="memory"),
                            ],
                            width=12,
                        )
                    ]
                ),
            ]
        )

        return app

    @staticmethod
    # https://github.com/sacdallago/bio_embeddings/blob/develop/bio_embeddings/visualize/plotly_plots.py
    def render(df: DataFrame, selected_column: str):
        """
        Renders the plotly graph with the selected column in the dataframe df
        :param df: dataframe
        :param selected_column: column of the dataframe
        :return: plotly graphical object
        """
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
                # 10 colors are available; once those are used, pick different symbol
                marker=dict(symbol=Visualizator.SYMBOLS[group_idx % 8]),
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
        fig.update_traces(hoverinfo="name", hoverlabel=dict(namelength=-1))
        return fig

    @staticmethod
    @dash.callback(
        Output("graph", "figure"),
        Input("store_data", "data"),
        Input("dd_menu", "value"),
    )
    def update_graph(df: DataFrame, selected_value: str):
        """
        Renders new graph for selected drop down menu value
        :param df: dataframe
        :param selected_value: selected column of dropdown menu
        :return: graph to be displayed
        """
        # Check whether an input is triggered
        ctx = dash.callback_context
        if not ctx.triggered:
            raise PreventUpdate

        # Convert records df form to origin
        df = DataFrame.from_records(df)

        fig = Visualizator.render(df, selected_column=selected_value)
        return fig
