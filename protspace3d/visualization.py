#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
from pandas import DataFrame

from dash import Dash, dcc, html
import dash_bootstrap_components as dbc

import dash_bio as dashbio

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

    representation_options = [
        {"label": "backbone", "value": "backbone"},
        {"label": "ball+stick", "value": "ball+stick"},
        {"label": "cartoon", "value": "cartoon"},
        {"label": "hyperball", "value": "hyperball"},
        {"label": "licorice", "value": "licorice"},
        {"label": "axes+box", "value": "axes+box"},
        {"label": "helixorient", "value": "helixorient"},
    ]

    def __init__(self, fig: go.Figure, csv_header: list[str]):
        self.fig = fig
        self.csv_header = csv_header

    @staticmethod
    def get_header(app: Dash):
        header = dbc.Row(
            [
                dbc.Col(
                    html.H1("ProtSpace3D", style={"color": "white"}),
                    width=4,
                    style={"background-color": "black"},
                ),
                dbc.Col(width=7, style={"background-color": "black"}),
                dbc.Col(
                    html.A(
                        [
                            html.Img(
                                src=app.get_asset_url("logo.png"),
                                alt="Rostlab-logo",
                                style={"height": "60px", "width": "60px"},
                            ),
                        ],
                        href="https://rostlab.org/",
                        target="_blank",
                    ),
                    style={"background-color": "black"},
                    width=1,
                ),
            ]
        )

        return header

    def init_app_pdb(self):
        """
        Initializes app & Builds html layout for Dash
        :return: layout
        """
        app = Dash(__name__, external_stylesheets=[dbc.themes.BOOTSTRAP])

        app.layout = dbc.Container(
            [
                # Header
                self.get_header(app),
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
                                        "width": "110%",
                                        "height": "90vh",
                                    },
                                ),
                            ],
                            width=8,
                        ),
                        dbc.Col(
                            [
                                dcc.Dropdown(
                                    id="representation_dropdown",
                                    options=self.representation_options,
                                    multi=True,
                                    value=["cartoon"],
                                ),
                                dashbio.NglMoleculeViewer(id="ngl_molecule_viewer"),
                            ],
                            width=4,
                        ),
                    ]
                ),
            ],
            fluid=True,
        )

        return app

    def init_app(self):
        """
        Initializes app & Builds html layout for Dash
        :return: layout
        """
        app = Dash(__name__, external_stylesheets=[dbc.themes.BOOTSTRAP])

        app.layout = dbc.Container(
            [
                # Header
                self.get_header(app),
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
                            width=12,
                        ),
                    ]
                ),
            ],
            fluid=True,
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
        # iterate over different values of the selected column
        for group_idx, group_value in enumerate(col_groups):
            # extract df with only group value
            df_group = df[df[selected_column] == group_value]
            trace = go.Scatter3d(
                x=df_group["x"],
                y=df_group["y"],
                z=df_group["z"],
                mode="markers",
                name=group_value,
                # 10 colors are available; once those are used, pick different symbol
                marker=dict(symbol=Visualizator.SYMBOLS[group_idx % 8]),
                text=df_group.index.to_list(),
            )
            data.append(trace)
            # Give the different group values a number
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

        # Set hoverinfo
        fig.update_traces(
            hoverinfo=["name", "text"],
            hoverlabel=dict(namelength=-1),
            hovertemplate="%{text}",
        )

        # Set legend in right upper corner of the plot
        fig.update_layout(legend=dict(yanchor="top", y=0.99, xanchor="right", x=0.99))

        return fig
