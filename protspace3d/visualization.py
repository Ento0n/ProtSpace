#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import sys
from colorsys import hls_to_rgb

import dash_bio as dashbio
import dash_bootstrap_components as dbc
import numpy as np
import plotly.graph_objects as go
from dash import Dash, dcc, html
from pandas import DataFrame


class Visualizator:
    SYMBOLS = [
        "circle",
        "square",
        "diamond",
        "cross",
        "circle-open",
        "diamond-open",
        "square-open",
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
    def n_samples_equation(n, max_out, min_val, reverse: bool = False) -> float:
        # the underlying equation is (x-10)/(x+10)
        numerator = n - 10
        denominator = n + 10
        res = numerator / denominator

        # less than 10 for n should be 0
        if res < 0:
            res = 0

        if reverse:
            res = 1 - res

        # out is the min value that can be picked
        out = min_val + max_out * res

        return out

    @staticmethod
    def gen_distinct_colors(n, sort: bool = True):
        color_list = list()
        np.random.seed(42)
        hues = np.arange(0, 360, 360 / n)
        hues = hues[np.random.permutation(hues.size)]
        for hue in hues:
            min_sat = Visualizator.n_samples_equation(n, 50, 30, reverse=True)
            saturation = min_sat + np.random.ranf() * (80 - min_sat)
            saturation = 100
            min_lum = Visualizator.n_samples_equation(n, 40, 40, reverse=True)
            luminosity = min_lum + np.random.ranf() * (80 - min_lum)
            luminosity = 50
            color_list.append(tuple([hue / 360, luminosity / 100, saturation / 100]))
        if sort:
            color_list.sort()

        # round small values, otherwise dash has difficulties to display
        rgb_list = []
        for h, l, s in color_list:
            rgb = hls_to_rgb(round(h, 2), round(l, 2), round(s, 2))

            rgb = list(rgb)
            for idx, value in enumerate(rgb):
                rgb[idx] = value * 255
            rgb = tuple(rgb)

            rgb_list.append(rgb)
        return rgb_list

    @staticmethod
    # https://github.com/sacdallago/bio_embeddings/blob/develop/bio_embeddings/visualize/plotly_plots.py
    def render(
        df: DataFrame,
        selected_column: str,
        original_id_col: object,
        umap_paras: dict,
        umap_flag: bool = True,
    ):
        """
        Renders the plotly graph with the selected column in the dataframe df
        :param df: dataframe
        :param selected_column: column of the dataframe
        :param original_id_col: the colum "original id" of the mapped csv file
        :param umap_flag: flag is set if umap calculations are used.
        :return: plotly graphical object
        """

        # custom separator to sort str, int and float (str case-insensitive)
        # order: 1. int and float 2. str 3. rest 4. NA
        def my_comparator(val):
            if isinstance(val, float) or isinstance(val, int):
                return 0, val
            elif val == "NA":
                return 3, val
            elif isinstance(val, str):
                val = val.lower()
                return 1, val
            else:
                return 2, val

        mapped_index = None
        if original_id_col is not None:
            # swap index
            mapped_index = df.index
            df.index = original_id_col

        col_groups = df[selected_column].unique().tolist()

        col_groups.sort(key=my_comparator)

        color_list = Visualizator.gen_distinct_colors(n=len(col_groups))

        fig = go.Figure()

        numeric_flag = False
        if all(
            [
                isinstance(item, int) or isinstance(item, float) or item == "NA"
                for item in col_groups
            ]
        ):
            # Only numeric values and "NA" in the group
            numeric_flag = True
            colorscale = list()

            # Find min and max value of the group, excluding "NA"
            min_val = sys.float_info.max
            max_val = sys.float_info.min
            for item in col_groups:
                if isinstance(item, int) or isinstance(item, float):
                    if min_val > item:
                        min_val = item
                    if max_val < item:
                        max_val = item

            # Fill the colorscale with the used colors for the groups
            for idx, item in enumerate(col_groups):
                if item != "NA":
                    colorscale.append(f"rgb{color_list[idx]}")

            # create colorbar
            colorbar = dict(
                title="Colorbar",
                lenmode="fraction",
                len=0.5,
                yanchor="bottom",
                ypad=50,
            )

            # create and add a dummy trace that holds the colorbar
            color_trace = go.Scatter3d(
                x=[None],
                y=[None],
                z=[None],
                mode="markers",
                marker=dict(
                    colorscale=colorscale,
                    showscale=True,
                    colorbar=colorbar,
                    cmin=min_val,
                    cmax=max_val,
                ),
                showlegend=False,
            )

            fig.add_trace(color_trace)

        # Figure out how many symbols to use depending on number of column groups
        n_symbols = int(
            Visualizator.n_samples_equation(
                n=len(col_groups), max_out=len(Visualizator.SYMBOLS) - 3, min_val=3
            )
        )

        df["class_index"] = np.ones(len(df)) * -100

        if umap_flag:
            x = "x_umap"
            y = "y_umap"
            z = "z_umap"
        else:
            x = "x_pca"
            y = "y_pca"
            z = "z_pca"

        # iterate over different values of the selected column
        for group_idx, group_value in enumerate(col_groups):
            # Only show NA in legend if colorbar is shown
            show_legend = True
            if numeric_flag:
                if group_value != "NA":
                    show_legend = False

            # extract df with only group value
            df_group = df[df[selected_column] == group_value]
            trace = go.Scatter3d(
                x=df_group[x],
                y=df_group[y],
                z=df_group[z],
                mode="markers",
                name=group_value,
                marker=dict(
                    color=f"rgb{color_list[group_idx]}",
                    symbol=Visualizator.SYMBOLS[group_idx % n_symbols],
                    line=dict(color="black", width=1),
                ),
                text=df_group.index.to_list(),
                showlegend=show_legend,
            )
            fig.add_trace(trace)
            # Give the different group values a number
            df.loc[df[selected_column] == group_value, "class_index"] = group_idx

        if umap_flag:
            fig.update_layout(
                # Remove axes ticks and labels as they are usually not informative
                scene=dict(
                    xaxis=dict(showticklabels=False, showspikes=False, title=""),
                    yaxis=dict(showticklabels=False, showspikes=False, title=""),
                    zaxis=dict(showticklabels=False, showspikes=False, title=""),
                ),
            )
        else:
            # extract variance column
            pca_variance = df["variance"].to_list()

            fig.update_layout(
                # Remove axes ticks and labels as they are usually not informative
                scene=dict(
                    xaxis=dict(
                        showticklabels=False,
                        showspikes=False,
                        title=f"PC1 ({float(pca_variance[0]):.1f}%)",
                    ),
                    yaxis=dict(
                        showticklabels=False,
                        showspikes=False,
                        title=f"PC2 ({float(pca_variance[1]):.1f}%)",
                    ),
                    zaxis=dict(
                        showticklabels=False,
                        showspikes=False,
                        title=f"PC3 ({float(pca_variance[2]):.1f}%)",
                    ),
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

        # Update title
        fig.update_layout(
            title={
                "text": "UMAP"
                + f"<br>n_neighbours: {umap_paras['n_neighbours']},"
                + f" min_dist: {umap_paras['min_dist']}, metric: {umap_paras['metric']}",
                "y": 0.98,
                "x": 0.4,
            }
        )

        # swap index again
        if original_id_col is not None:
            df.index = mapped_index

        return fig

    @staticmethod
    def get_app():
        app = Dash(
            __name__, external_stylesheets=[dbc.themes.BOOTSTRAP, dbc.icons.BOOTSTRAP]
        )

        return app

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
    def get_download_toast():
        toast = dbc.Toast(
            "Html file successfully saved in output folder!",
            header="HTML created",
            id="html_download_toast",
            is_open=False,
            dismissable=True,
            duration=4000,
            style={
                "position": "fixed",
                "top": 66,
                "left": 10,
                "width": 350,
            },
        )

        return toast

    def get_graph_container(self):
        graph_container = (
            dbc.Offcanvas(
                id="graph_offcanvas",
                is_open=False,
                title="Graph settings",
                children=[
                    dcc.Markdown("Dimensionality reduction"),
                    dbc.RadioItems(
                        options=[
                            {"label": "PCA", "value": "PCA"},
                            {"label": "UMAP", "value": "UMAP"},
                        ],
                        value="UMAP",
                        id="dim_red_radio",
                        inline="True",
                    ),
                    dcc.Markdown("HTML"),
                    dbc.Row(
                        [
                            dbc.Col(
                                [
                                    dcc.Dropdown(
                                        self.csv_header,
                                        self.csv_header[0],
                                        id="html_dd",
                                    )
                                ],
                                width=9,
                            ),
                            dbc.Col(
                                [
                                    dbc.Button(
                                        "Download",
                                        id="html_download_button",
                                        outline=True,
                                        color="dark",
                                    )
                                ]
                            ),
                        ]
                    ),
                ],
                style={"width": "50%"},
                placement="end",
            ),
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
                        ],
                        width=10,
                    ),
                    dbc.Col(
                        [
                            dbc.Button(
                                "",
                                id="graph_settings_button",
                                class_name="bi bi-gear-wide-connected",
                                outline=True,
                                color="dark",
                                style={
                                    "margin-top": "5px",
                                    "margin-bottom": "5px",
                                    "margin-right": "15px",
                                },
                            ),
                        ],
                        width=2,
                    ),
                ]
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
        )

        return graph_container

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
                # Toast that is displayed if a html file is created
                self.get_download_toast(),
                # graph and controls
                dbc.Row(
                    [
                        dbc.Col(
                            self.get_graph_container(),
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
                                                    "",
                                                    id="molecules_settings_button",
                                                    class_name="bi bi-gear-wide-connected",
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
                                    id="molecules_offcanvas",
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
                                            min=3,
                                            max=9,
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
                # toast that is displayed if a html file is created
                self.get_download_toast(),
                # space between header and content below
                dbc.Row([html.Br()]),
                # graph and controls
                dbc.Row(
                    [
                        dbc.Col(
                            self.get_graph_container(),
                            width=12,
                        ),
                    ]
                ),
            ],
            fluid=True,
        )

        return app
