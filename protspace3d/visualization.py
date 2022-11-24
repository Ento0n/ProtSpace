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

    def get_header(self, app: Dash):
        header = dbc.Row(
            [
                dbc.Col(
                    html.H1("ProtSpace3D", style={"color": "white"}),
                    style={"background-color": "black"},
                    xxl=10,
                    xl=10,
                    lg=9,
                    md=9,
                    sm=7,
                    xs=7,
                ),
                dbc.Col(
                    [
                        dbc.Button(
                            "",
                            id="help_button",
                            class_name="bi bi-question-lg",
                            color="dark",
                            style={
                                "margin-top": "10px",
                                "margin-bottom": "5px",
                                "margin-right": "20px",
                                "background-color": "black",
                            },
                        ),
                        self.get_help_modal(),
                    ],
                    xxl=1,
                    xl=1,
                    lg=2,
                    md=2,
                    sm=3,
                    xs=3,
                    style={"background-color": "black"},
                ),
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
                    xxl=1,
                    xl=1,
                    lg=1,
                    md=1,
                    sm=2,
                    xs=2,
                ),
            ]
        )

        return header

    @staticmethod
    def get_disclaimer_modal():
        modal = dbc.Modal(
            [
                dbc.ModalHeader(dbc.ModalTitle("Disclaimer"), close_button=False),
                dbc.ModalBody("Here the disclaimer text is shown!"),
                dbc.ModalFooter(dbc.Button("Agree", id="disclaimer_modal_button")),
            ],
            id="disclaimer_modal",
            size="xl",
            is_open=True,
            backdrop="static",
        )

        return modal

    @staticmethod
    def get_help_modal():
        modal = dbc.Modal(
            [
                dbc.ModalHeader(dbc.ModalTitle("Help"), close_button=True),
                dbc.ModalBody("Here the help text is shown!"),
            ],
            id="help_modal",
            size="xl",
            is_open=False,
            backdrop=True,
        )

        return modal

    @staticmethod
    def get_app():
        app = Dash(
            __name__, external_stylesheets=[dbc.themes.BOOTSTRAP, dbc.icons.BOOTSTRAP]
        )

        return app

    def init_app_pdb(self, original_id_col: list):
        """
        Initializes app & Builds html layout for Dash
        :return: layout
        """
        app = self.get_app()

        app.layout = dbc.Container(
            [
                # Header
                self.get_header(app),
                # sizing of the molecule viewer
                dcc.Location(id="url"),
                html.Div(id="molviewer_sizing_div", hidden=True),
                # storage to save the selected molecules
                # Needed for image download name
                dcc.Store(id="mol_name_storage"),
                # graph and controls
                dbc.Row(
                    [
                        dbc.Col(
                            [
                                dcc.Dropdown(
                                    self.csv_header,
                                    self.csv_header[0],
                                    id="dd_menu",
                                    style={"margin-top": "5px"},
                                ),
                                dcc.Graph(
                                    id="graph",
                                    figure=self.fig,
                                    clear_on_unhover=True,
                                    style={
                                        "width": "100%",
                                        "height": "80vh",
                                    },
                                    responsive=True,
                                ),
                            ],
                            id="left_col",
                            width=6,
                            style={"border-right": "solid black 1px"},
                        ),
                        dbc.Col(
                            [
                                dbc.Row(
                                    [
                                        dbc.Col(
                                            [
                                                dcc.Markdown(
                                                    "Molecules:",
                                                    style={
                                                        "margin-top": "20px",
                                                        "margin-bottom": "0px",
                                                        "padding-top": "0px",
                                                        "padding-bottom": "0px",
                                                        "height": "30px",
                                                    },
                                                ),
                                            ],
                                            xxl=9,
                                            xl=8,
                                            lg=7,
                                            md=6,
                                            sm=5,
                                            xs=4,
                                        ),
                                        dbc.Col(
                                            [
                                                dbc.Button(
                                                    "",
                                                    id="reset_view_button",
                                                    class_name="bi bi-arrow-counterclockwise",
                                                    color="dark",
                                                    outline=True,
                                                    style={
                                                        "margin-top": "5px",
                                                        "margin-bottom": "5px",
                                                        "margin-right": "20px",
                                                    },
                                                ),
                                                dbc.Button(
                                                    "Settings",
                                                    id="settings_button",
                                                    outline=True,
                                                    color="dark",
                                                    style={
                                                        "margin-top": "5px",
                                                        "margin-bottom": "5px",
                                                    },
                                                ),
                                            ],
                                            xxl=3,
                                            xl=4,
                                            lg=5,
                                            md=6,
                                            sm=7,
                                            xs=8,
                                        ),
                                    ]
                                ),
                                dcc.Dropdown(
                                    id="molecules_dropdown",
                                    options=original_id_col,
                                    multi=True,
                                    style={"margin-bottom": "5px"},
                                ),
                                html.Div(
                                    [
                                        dashbio.NglMoleculeViewer(
                                            id="ngl_molecule_viewer",
                                        ),
                                    ],
                                    id="moleculeviewer_div",
                                    style={
                                        "border-bottom": "1px solid grey",
                                        "border-right": "1px solid grey",
                                        "margin-left": "0px",
                                    },
                                ),
                                dbc.Offcanvas(
                                    id="offcanvas",
                                    title="Settings",
                                    is_open=False,
                                    children=[
                                        dcc.Markdown("Representations:"),
                                        dcc.Dropdown(
                                            id="representation_dropdown",
                                            options=self.representation_options,
                                            multi=True,
                                            value=["cartoon"],
                                        ),
                                        html.Br(),
                                        dbc.Row(
                                            [
                                                dbc.Col(
                                                    [
                                                        dcc.Markdown("Start:"),
                                                        dcc.Dropdown(
                                                            id="range_start",
                                                            disabled=True,
                                                        ),
                                                    ]
                                                ),
                                                dbc.Col(
                                                    [
                                                        dcc.Markdown("End:"),
                                                        dcc.Dropdown(
                                                            id="range_end",
                                                            disabled=True,
                                                        ),
                                                    ]
                                                ),
                                                dbc.Col(
                                                    [
                                                        dcc.Markdown(
                                                            "Highlighted atoms:",
                                                        ),
                                                        dcc.Dropdown(
                                                            id="selected_atoms",
                                                            multi=True,
                                                            disabled=True,
                                                        ),
                                                    ]
                                                ),
                                            ]
                                        ),
                                        html.Br(),
                                        dcc.Markdown("Spacing:"),
                                        dcc.Slider(
                                            id="spacing_slider",
                                            min=10,
                                            max=200,
                                            value=50,
                                            marks=None,
                                            tooltip={
                                                "placement": "bottom",
                                                "always_visible": False,
                                            },
                                        ),
                                        dcc.Markdown("Space distribution:"),
                                        dcc.Slider(
                                            id="distribution_slider",
                                            min=0,
                                            max=12,
                                            value=6,
                                            step=1,
                                            marks=None,
                                        ),
                                        dbc.Button(
                                            "Recalculate molecule viewing size",
                                            id="recal_size_button",
                                            class_name="d-grid mx-auto",
                                            color="dark",
                                            outline=True,
                                            style={"margin-bottom": "10px"},
                                        ),
                                        dcc.Markdown("Height:"),
                                        dcc.Slider(
                                            id="height_slider",
                                            min=200,
                                            max=5000,
                                            value=300,
                                            marks=None,
                                            tooltip={
                                                "placement": "bottom",
                                                "always_visible": False,
                                            },
                                        ),
                                        dcc.Markdown("Width:"),
                                        dcc.Slider(
                                            id="width_slider",
                                            min=200,
                                            max=5000,
                                            value=500,
                                            marks=None,
                                            tooltip={
                                                "placement": "bottom",
                                                "always_visible": False,
                                            },
                                        ),
                                        dbc.Row(
                                            [
                                                dbc.Col(
                                                    dcc.Input(
                                                        id="filename_input",
                                                        type="text",
                                                        placeholder="filename",
                                                        style={
                                                            "height": "38px",
                                                            "margin-right": "20px",
                                                        },
                                                    ),
                                                    width=6,
                                                ),
                                                dbc.Col(
                                                    dbc.Button(
                                                        "Download image",
                                                        id="download_molecule_button",
                                                        color="dark",
                                                        outline=True,
                                                        disabled=True,
                                                    ),
                                                    width=6,
                                                ),
                                            ]
                                        ),
                                    ],
                                    style={"width": "50%"},
                                ),
                            ],
                            id="right_col",
                            width=6,
                        ),
                    ]
                ),
                # modal with disclaimer that opens on startup
                # has to be at the end, otherwise automatic sizing doesn't work...
                self.get_disclaimer_modal(),
            ],
            fluid=True,
        )

        return app

    def init_app(self):
        """
        Initializes app & Builds html layout for Dash
        :return: layout
        """
        app = self.get_app()

        app.layout = dbc.Container(
            [
                # Header
                self.get_header(app),
                # model with disclaimer that opens on startup
                self.get_disclaimer_modal(),
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
                                ),
                                dcc.Graph(
                                    id="graph",
                                    figure=self.fig,
                                    clear_on_unhover=True,
                                    style={
                                        "width": "100%",
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
    def render(df: DataFrame, selected_column: str, original_id_col: object):
        """
        Renders the plotly graph with the selected column in the dataframe df
        :param df: dataframe
        :param selected_column: column of the dataframe
        :param original_id_col: the colum "original id" of the mapped csv file
        :return: plotly graphical object
        """

        # custom separator to sort str, int and float (str case-insensitive)
        # order: 1. int and float 2. str 3. rest 4. NA
        def my_comparator(item):
            if isinstance(item, float) or isinstance(item, int):
                return 0, item
            elif item == "NA":
                return 3, item
            elif isinstance(item, str):
                item = item.lower()
                return 1, item
            else:
                return 2, item

        mapped_index = None
        if original_id_col is not None:
            # swap index
            mapped_index = df.index
            df.index = original_id_col

        col_groups = df[selected_column].unique().tolist()

        col_groups.sort(key=my_comparator)

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
        fig.update_layout(legend=dict(yanchor="top", y=0.97, xanchor="right", x=0.99))

        # change margins of the graph
        fig.update_layout(margin=dict(l=1, r=1, t=1, b=1))

        # swap index again
        if original_id_col is not None:
            df.index = mapped_index

        return fig
