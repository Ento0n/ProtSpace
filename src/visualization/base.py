#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import dash_bootstrap_components as dbc
from dash import Dash, dcc, html
import plotly.graph_objects as go

metric_options = [
    "euclidean",
    "cosine",
    "manhattan",
    "chebyshev",
    "minkowski",
    "canberra",
    "braycurtis",
    "haversine",
    "mahalanobis",
    "wminkowski",
    "seuclidean",
    "correlation",
    "hamming",
    "jaccard",
    "dice",
    "russellrao",
    "kulsinski",
    "rogerstanimoto",
    "sokalmichener",
    "sokalneath",
    "yule",
]

main_button_style = {
    "margin-top": "5px",
    "margin-bottom": "5px",
    "margin-right": "15px",
}


def get_app():
    """
    Initializes dash application
    :return: application
    """
    app = Dash(
        __name__,
        external_stylesheets=[dbc.themes.BOOTSTRAP, dbc.icons.BOOTSTRAP],
    )

    return app


def init_app(umap_paras: dict, csv_header: list[str], fig: go.Figure):
    """
    Set up the layout of the application
    :return: application
    """
    app = get_app()

    app.layout = dbc.Container(
        [
            # Header
            get_header(app),
            # model with disclaimer that opens on startup
            get_disclaimer_modal(),
            # toast that is displayed if a html file is created
            get_download_toast(),
            # toast that displays the information of a selected protein
            get_info_toast(),
            # graph and controls
            dbc.Row(
                [
                    dbc.Col(
                        get_graph_container(umap_paras, False, csv_header, fig),
                        width=12,
                    ),
                ]
            ),
        ],
        fluid=True,
    )

    return app


def get_graph_offcanvas(umap_paras: dict, umap_paras_string: str):
    """
    Creates layout of the offcanvas for the graph.
    :param umap_paras: Parameters of the UMAP calculation
    :param umap_paras_string: UMAP parameters in string format.
    :return: graph offcanvas layout
    """
    offcanvas = dbc.Offcanvas(
        id="graph_offcanvas",
        is_open=False,
        title="Graph settings",
        style={"width": "50%"},
        placement="end",
        children=[
            dcc.Markdown("Dimensionality reduction"),
            dbc.RadioItems(
                options=[
                    {"label": "UMAP", "value": "UMAP"},
                    {"label": "PCA", "value": "PCA"},
                ],
                value="UMAP",
                id="dim_red_radio",
                inline="True",
            ),
            html.Br(),
            dcc.Markdown("HTML"),
            dbc.Button(
                "Download all files",
                id="button_html_all",
                color="dark",
                outline=True,
            ),
            html.Br(),
            html.Br(),
            dcc.Markdown("UMAP parameters"),
            dbc.Row(
                [
                    dbc.Col(
                        [
                            dcc.Markdown("n_neighbours:"),
                            dbc.Input(
                                id="n_neighbours_input",
                                type="number",
                                min=0,
                                step=1,
                                value=umap_paras["n_neighbours"],
                            ),
                        ],
                        width=4,
                    ),
                    dbc.Col(
                        [
                            dcc.Markdown("min_dist:"),
                            dbc.Input(
                                id="min_dist_input",
                                type="number",
                                min=0,
                                # Set max to 1 since min_dist must be less than or equal to spread,
                                # which is 1.0 by default and not changed
                                max=1,
                                step=0.1,
                                value=umap_paras["min_dist"],
                            ),
                        ],
                        width=4,
                    ),
                    dbc.Col(
                        [
                            dcc.Markdown("metric:"),
                            dcc.Dropdown(
                                id="metric_input",
                                options=metric_options,
                                value=umap_paras["metric"],
                            ),
                        ],
                        width=4,
                    ),
                ],
            ),
            html.Br(),
            dbc.Button(
                "Recalculate UMAP",
                id="umap_recalculation_button",
                color="dark",
                outline=True,
            ),
            html.Br(),
            html.Br(),
            dcc.Dropdown(
                id="last_umap_paras_dd",
                value=umap_paras_string,
                options=[umap_paras_string],
                clearable=False,
                searchable=False,
            ),
        ],
    )

    return offcanvas


def get_graph_container(
    umap_paras: dict, pdb: bool, csv_header: list[str], fig: go.Figure
):
    """
    Creates the layout for the graph Row
    :param umap_paras: umap parameters
    :param pdb: flag whether pdb layout is needed or not.
    :param csv_header: headers of the csv file
    :param fig: graph Figure
    :return: Layout of the offcanvas
    """
    # UMAP parameters in string format
    umap_paras_string = str(
        str(umap_paras["n_neighbours"])
        + " ; "
        + str(umap_paras["min_dist"])
        + " ; "
        + umap_paras["metric"],
    )

    # width sizing of the dropdown menu column
    if pdb:
        xs = 6
        sm = 7
        md = 8
        xxl = 9
    else:
        xs = 8
        sm = 9
        md = 10
        xxl = 11

    graph_container = (
        # Storage to save whether Highlighting circle is already displayed or not
        dcc.Store(id="highlighting_bool", storage_type="memory", data=False),
        # Storage to save last camera data (relayoutData)
        dcc.Store(id="relayoutData_save", storage_type="memory", data={}),
        get_graph_offcanvas(umap_paras, umap_paras_string),
        dbc.Row(
            children=[
                dbc.Col(
                    [
                        dcc.Dropdown(
                            csv_header,
                            csv_header[0],
                            id="dd_menu",
                            style={"margin-top": "5px"},
                        ),
                    ],
                    xs=xs,
                    sm=sm,
                    md=md,
                    xxl=xxl,
                ),
                dbc.Col(
                    xs=12 - xs,
                    sm=12 - sm,
                    md=12 - md,
                    xxl=12 - xxl,
                    children=[
                        dbc.Stack(
                            direction="horizontal",
                            children=[
                                dbc.Button(
                                    "",
                                    class_name="bi bi-download",
                                    id="html_download_button",
                                    outline=True,
                                    color="dark",
                                    style=main_button_style,
                                ),
                                dbc.Button(
                                    "",
                                    id="graph_settings_button",
                                    class_name="bi bi-gear-wide-connected",
                                    outline=True,
                                    color="dark",
                                    style=main_button_style,
                                ),
                            ],
                        ),
                    ],
                ),
            ],
        ),
        dcc.Graph(
            id="graph",
            figure=fig,
            clear_on_unhover=True,
            style={
                "width": "100%",
                "height": "80vh",
            },
            responsive=True,
        ),
    )

    return graph_container


def get_disclaimer_modal():
    """
    Layout for the modal (window) with the disclaimer displayed at startup of the website.
    :return: Disclaimer layout modal
    """
    modal = dbc.Modal(
        [
            dbc.ModalHeader(dbc.ModalTitle("Disclaimer"), close_button=False),
            dbc.ModalBody(
                children=[
                    dbc.Alert(
                        color="warning",
                        children="All fonts and representations used un this app are based on dash, "
                        "dash bootstrap components and dash_bio. No liability is accepted by the creators "
                        "of this website.",
                    )
                ]
            ),
            dbc.ModalFooter(dbc.Button("Agree", id="disclaimer_modal_button")),
        ],
        id="disclaimer_modal",
        size="xl",
        is_open=True,
        backdrop="static",
    )

    return modal


def get_help_modal():
    """
    Layout for the modal (window) with the help text.
    :return: Help layout modal
    """
    modal = dbc.Modal(
        [
            dbc.ModalHeader(dbc.ModalTitle("Help"), close_button=True),
            dbc.ModalBody(
                children=[
                    dcc.Markdown(
                        """
                        Graph
                                Navigation
                                    Orbital rotation by clicking and dragging with left mouse button
                                    Pan by clicking and dragging with right mouse button
                                Selection
                                    By clicking on a molecule-dot in the graph, the infobox is shown, 
                                    and the selected molecule is indicated by a black circle around it.
                    """
                    ),
                ]
            ),
        ],
        id="help_modal",
        size="xl",
        is_open=False,
        backdrop=True,
    )

    return modal


def get_download_toast():
    """
    Layout for the toast (infobox) shown when downloading an html file.
    :return: download toast layout
    """
    toast = dbc.Toast(
        "Html file(s) successfully saved in output folder!",
        header="Download HTML",
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


def get_info_toast():
    """
    Layout for the toast (infobox) shown when clicking on a molecule in the graph.
    :return: info toast layout
    """
    toast = dbc.Toast(
        children=[
            html.Div(id="expanded_seq_div"),
            html.Div(id="collapsed_seq_div"),
            html.Button(id="expand_seq_button"),
            html.Button(id="collapse_seq_button"),
            html.Div(id="group_info_expanded_div"),
            html.Div(id="group_info_collapsed_div"),
            dbc.Button(id="group_info_expand_button"),
            dbc.Button(id="group_info_collapse_button"),
        ],
        id="info_toast",
        is_open=False,
        dismissable=True,
        body_style={
            "max-height": "50vh",
            "overflow": "auto",
        },
        style={
            "position": "fixed",
            "top": 166,
            "left": 10,
            "width": 200,
        },
    )

    return toast


def get_header(app: Dash):
    """
    Layout for the black header of the application
    :param app: the application itself
    :return: layout of the header
    """
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
                    get_help_modal(),
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
