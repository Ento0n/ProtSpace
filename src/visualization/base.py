#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import dash_bootstrap_components as dbc
import plotly.graph_objects as go
from dash import Dash, dcc, html

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

help_modal_icon_style = {
    "position": "relative",
    "bottom": "8px",
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


def init_app(
    umap_paras: dict,
    csv_header: list[str],
    fig: go.Figure,
    dim_red: str,
    tsne_paras: dict,
    original_id_col: list,
):
    """
    Set up the layout of the application
    :return: application
    """
    app = get_app()

    app.layout = dbc.Container(
        [
            # get all side components like header, toasts etc
            get_side_components(app),
            # modal with disclaimer that opens on startup
            get_disclaimer_modal(),
            # graph and controls
            dbc.Row(
                [
                    dbc.Col(
                        get_graph_container(
                            umap_paras,
                            False,
                            csv_header,
                            fig,
                            dim_red,
                            tsne_paras,
                            original_id_col,
                        ),
                        width=12,
                    ),
                ]
            ),
        ],
        fluid=True,
    )

    return app


def get_side_components(app: Dash):
    """
    Collection of side components like header, toasts etc. shared by normal and pdb app
    :param app: the Dash application
    :return: side components in a html Div object
    """

    side_components = html.Div(
        children=[
            # Header
            get_header(app),
            # Storage to save the clicked on molecule in the graph,
            # needed for replacing the clicked molecule in the list
            dcc.Store(id="clicked_mol_storage"),
            # toast that is displayed if a html file is created
            get_download_toast(),
            # Toasts in container, so they stack below each other...
            dbc.Container(
                [
                    # toast that displays the information of a selected protein
                    get_info_toast(),
                    html.Br(),
                    # toast that displays the nearest neighbours of the selected points
                    get_neighbour_toast(),
                ],
                style={
                    "position": "fixed",
                    "top": 166,
                    "left": 10,
                    "width": 350,
                },
            ),
        ]
    )

    return side_components


def get_graph_offcanvas(
    umap_paras: dict,
    umap_paras_string: str,
    dim_red: str,
    tsne_paras: dict,
    tsne_paras_string: str,
):
    """
    Creates layout of the offcanvas for the graph.
    :param umap_paras: Parameters of the UMAP calculation
    :param umap_paras_string: UMAP parameters in string format.
    :param dim_red: the initial dimensionality reduction
    :param tsne_paras: Parameters of the TSNE calculation
    :return: graph offcanvas layout
    """
    offcanvas = dbc.Offcanvas(
        id="graph_offcanvas",
        is_open=False,
        title="Graph settings",
        style={"width": "50%", "max-width": "600px"},
        placement="end",
        children=[
            dbc.Spinner(
                children=[
                    html.Div(
                        id="load_umap_spinner",
                        hidden=True,
                    )
                ],
                fullscreen=True,
                delay_show=80,
            ),
            dcc.Markdown("Download"),
            dbc.Button(
                "Download all files",
                id="button_graph_all",
                color="dark",
                outline=True,
            ),
            html.Br(),
            html.Br(),
            dcc.Markdown("Dimensions"),
            dbc.RadioItems(
                options=[
                    {"label": "3D", "value": "3D"},
                    {"label": "2D", "value": "2D"},
                ],
                value="3D",
                id="dim_radio",
                inline=True,
            ),
            html.Br(),
            dbc.Switch(
                id="nearest_neighbours_switch",
                label="Display nearest neighbours",
                value=False,
            ),
            html.Br(),
            dbc.Switch(
              id="correlation_collapse_switch",
              label="Display correlation scores",
              value=False,
            ),
            html.Br(),
            dcc.Markdown("Dimensionality reduction"),
            dbc.Tabs(
                id="dim_red_tabs",
                children=[
                    dbc.Tab(
                        label="UMAP",
                        tab_id="UMAP",
                        children=[
                            html.Br(),
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
                                                value=umap_paras[
                                                    "n_neighbours"
                                                ],
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
                    ),
                    dbc.Tab(label="PCA", tab_id="PCA"),
                    dbc.Tab(
                        label="t-SNE",
                        tab_id="TSNE",
                        children=[
                            html.Br(),
                            dbc.Row(
                                children=[
                                    dbc.Col(
                                        children=[
                                            dcc.Markdown("Iterations:"),
                                            dbc.Input(
                                                id="iterations_input",
                                                type="number",
                                                min=250,
                                                step=10,
                                                value=tsne_paras["iterations"],
                                            ),
                                        ],
                                        width=4,
                                    ),
                                    dbc.Col(
                                        children=[
                                            dcc.Markdown("Perplexity:"),
                                            dbc.Input(
                                                id="perplexity_input",
                                                type="number",
                                                min=1,
                                                step=1,
                                                value=tsne_paras["perplexity"],
                                            ),
                                        ],
                                        width=4,
                                    ),
                                    dbc.Col(
                                        children=[
                                            dcc.Markdown("Learning rate:"),
                                            dbc.Input(
                                                id="learning_rate_input",
                                                min=1,
                                                step=1,
                                                value=tsne_paras[
                                                    "learning_rate"
                                                ],
                                            ),
                                        ],
                                        width=4,
                                    ),
                                ]
                            ),
                            html.Br(),
                            dcc.Markdown("Metric:"),
                            dcc.Dropdown(
                                id="tsne_metric_input",
                                options=metric_options,
                                value=tsne_paras["tsne_metric"],
                            ),
                            html.Br(),
                            dbc.Button(
                                "Recalculate t-SNE",
                                id="tsne_recalculation_button",
                                color="dark",
                                outline=True,
                            ),
                            html.Br(),
                            html.Br(),
                            dcc.Dropdown(
                                id="last_tsne_paras_dd",
                                value=tsne_paras_string,
                                options=[tsne_paras_string],
                                clearable=False,
                                searchable=False,
                            ),
                        ],
                    ),
                ],
                active_tab=dim_red,
            ),
        ],
    )

    return offcanvas


def get_settings_button_tooltip(button_id: str):
    """
    Returns the tooltip for the settings buttons
    :param button_id: target button id the tooltip is to be displayed
    :return: tooltip
    """
    tooltip = dbc.Tooltip(
        "Settings",
        target=button_id,
        placement="bottom",
    )
    return tooltip


def get_graph_download_button_tooltip(button_id: str):
    """
    Returns the tooltip for the download graph button
    :param button_id: target button id the tooltip is to be displayed
    :return: tooltip
    """
    tooltip = dbc.Tooltip(
        "Graph download", target=button_id, placement="bottom"
    )
    return tooltip


def get_graph_container(
    umap_paras: dict,
    pdb: bool,
    csv_header: list[str],
    fig: go.Figure,
    dim_red: str,
    tsne_paras: dict,
    original_id_col: list,
):
    """
    Creates the layout for the graph Row
    :param umap_paras: umap parameters
    :param pdb: flag whether pdb layout is needed or not.
    :param csv_header: headers of the csv file
    :param fig: graph Figure
    :param dim_red: initial dimensionality reduction
    :param tsne_paras: TSNE parameters
    :param original_id_col: list with the original IDs
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

    # TSNE parameters in string format
    tsne_paras_string = str(
        str(tsne_paras["iterations"])
        + " ; "
        + str(tsne_paras["perplexity"])
        + " ; "
        + str(tsne_paras["learning_rate"])
        + " ; "
        + str(tsne_paras["tsne_metric"])
    )

    # width sizing of the dropdown menu column and whether 2 dropdowns are above the graph or not
    if pdb:
        xs = 6
        sm = 7
        md = 8
        xxl = 9

        dropdown_s = dcc.Dropdown(
            csv_header,
            csv_header[0],
            id="dd_menu",
            style={"margin-top": "5px"},
        )
    else:
        xs = 8
        sm = 9
        md = 10
        xxl = 11

        dropdown_s = dbc.Row(
            children=[
                dbc.Col(
                    children=[
                        dcc.Dropdown(
                            csv_header,
                            csv_header[0],
                            id="dd_menu",
                            style={"margin-top": "5px"},
                        ),
                    ],
                    width=6,
                ),
                dbc.Col(
                    children=[
                        dcc.Dropdown(
                            id="molecules_dropdown",
                            options=original_id_col,
                            multi=True,
                            style={"margin-top": "5px"},
                        ),
                    ],
                    width=6,
                ),
            ],
        )

    graph_container = (
        # Storage to save whether Highlighting circle is already displayed or not
        dcc.Store(id="highlighting_bool", storage_type="memory", data=False),
        # Storage to save last camera data (relayoutData)
        dcc.Store(id="relayoutData_save", storage_type="memory", data={}),
        get_graph_offcanvas(
            umap_paras,
            umap_paras_string,
            dim_red,
            tsne_paras,
            tsne_paras_string,
        ),
        get_settings_button_tooltip(button_id="graph_settings_button"),
        get_graph_download_button_tooltip(button_id="graph_download_button"),
        get_correlation_scores_collapse(),
        dbc.Row(
            children=[
                dbc.Col(
                    [
                        dropdown_s,
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
                                    id="graph_download_button",
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


def get_correlation_scores_collapse():
    collapse = dbc.Collapse(
        dbc.Card(dbc.CardBody("This is a test!")),
        id="correlation_collapse",
        is_open=False,
        style={"margin-top": "5px"}
    )

    return collapse


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
                        children=[
                            dcc.Markdown(
                                """
                            This is the ProtSpace3D tool. There is no usage summary statistics provided by
                            the developers.
                            """
                            ),
                            dcc.Markdown(
                                """
                            Considering the data collection in the background and the displayed in the web application,
                            the Python library Dash is used and liable.
                            """
                            ),
                            dcc.Markdown(
                                """
                            Please refer to the policy of its developers Plotly: https://plotly.com/privacy/
                            for more information.
                            """
                            ),
                        ],
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
                    dbc.ListGroup(
                        children=[
                            dbc.ListGroupItem(
                                children=[
                                    html.H4("Graph"),
                                    html.Br(),
                                    html.H5("Group Selection"),
                                    html.P(
                                        "Group selection is done by clicking on"
                                        " the dropdown menu above the graph and"
                                        " selecting the wanted group."
                                    ),
                                    html.H5("Buttons"),
                                    dbc.Stack(
                                        direction="horizontal",
                                        gap=3,
                                        children=[
                                            html.I(
                                                className="bi bi-download",
                                                style=help_modal_icon_style,
                                            ),
                                            html.P(
                                                "Download a html or png file of"
                                                " the selected group."
                                            ),
                                        ],
                                    ),
                                    dbc.Stack(
                                        direction="horizontal",
                                        gap=3,
                                        children=[
                                            html.I(
                                                className=(
                                                    "bi bi-gear-wide-connected"
                                                ),
                                                style=help_modal_icon_style,
                                            ),
                                            html.P(
                                                "Open the settings to the"
                                                " graph."
                                            ),
                                        ],
                                    ),
                                    html.H5("Navigation"),
                                    html.B("Orbital rotation:"),
                                    html.P(
                                        "Click and hold the left mouse button."
                                    ),
                                    html.B("Pan:"),
                                    html.P(
                                        "Click and hold the right mouse button."
                                    ),
                                    html.B("Zoom:"),
                                    html.P(
                                        "Scrolling with the mouse wheel zooms"
                                        " in and out while in the graph with"
                                        " the cursor"
                                    ),
                                    html.H5("Molecule selection"),
                                    html.P(
                                        "A molecule is selected by clicking on"
                                        " the corresponding dot in the graph"
                                    ),
                                    html.H5("Legend"),
                                    html.P(
                                        "By clicking on a group in the legend,"
                                        " it is hidden. Clicking on a hidden"
                                        " group shows it again."
                                    ),
                                    html.P(
                                        "Double click on a displayed group to"
                                        " only show the selected group. Double"
                                        " click again on it displays all"
                                        " groups."
                                    ),
                                ]
                            ),
                            dbc.ListGroupItem(
                                children=[
                                    html.H4("Molecule viewer"),
                                    html.Br(),
                                    html.H5("Molecule selection"),
                                    html.P(
                                        "In the dropdown menu above the"
                                        " molecule viewer the selected"
                                        " molecule(s) are displayed. More can"
                                        " be selected by opening the"
                                        " dropdown-menu and seledting more."
                                        " Less can be selected by clicking on"
                                        " the x next to the molecule ID"
                                    ),
                                    html.H5("Buttons"),
                                    dbc.Stack(
                                        direction="horizontal",
                                        gap=3,
                                        children=[
                                            html.I(
                                                className=(
                                                    "bi bi-arrow-counterclockwise"
                                                ),
                                                style=help_modal_icon_style,
                                            ),
                                            html.P(
                                                "Reset the view of the molecule"
                                                " viewer."
                                            ),
                                        ],
                                    ),
                                    dbc.Stack(
                                        direction="horizontal",
                                        gap=3,
                                        children=[
                                            html.I(
                                                className=(
                                                    "bi bi-gear-wide-connected"
                                                ),
                                                style=help_modal_icon_style,
                                            ),
                                            html.P(
                                                "Open the settings to the"
                                                " molecule viewer."
                                            ),
                                        ],
                                    ),
                                    html.H5("Navigation"),
                                    html.B("Orbital rotation:"),
                                    html.P(
                                        "Click and hold the left mouse button."
                                    ),
                                    html.B("Pan:"),
                                    html.P(
                                        "Click and hold the right mouse button."
                                    ),
                                    html.B("Zoom:"),
                                    html.P(
                                        "Scrolling with the mouse wheel zooms"
                                        " in and out while in the graph with"
                                        " the cursor"
                                    ),
                                ]
                            ),
                        ]
                    )
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
    Layout for the toast (infobox) shown when downloading the graph as file.
    :return: download toast layout
    """
    toast = dbc.Toast(
        "File(s) successfully saved in output folder!",
        header="Download graph",
        id="graph_download_toast",
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
            "max-height": "35vh",
            "overflow": "auto",
        },
        style={"width": 200},
        header_style={"overflow": "auto"},
    )

    return toast


def get_neighbour_toast():
    """
    Layout for the toast (window showing nearest neighbours) shown when clicking on a molecule in the graph.
    :return: neihgbour toast layout
    """
    toast = dbc.Toast(
        id="neighbour_toast",
        is_open=False,
        dismissable=True,
        body_style={
            "max-height": "35vh",
            "overflow": "auto",
        },
    )

    return toast


def get_help_button_tooltip(button_id: str):
    """
    Returns the tooltip for the help button
    :param button_id: target button id the tooltip is to be displayed
    :return: tooltip
    """
    tooltip = dbc.Tooltip(
        "Help",
        target=button_id,
        placement="left",
    )
    return tooltip


def get_header(app: Dash):
    """
    Layout for the black header of the application
    :param app: the application itself
    :return: layout of the header
    """
    header = dbc.Row(
        [
            dbc.Col(
                dbc.Stack(
                    direction="horizontal",
                    gap=5,
                    children=[
                        html.H1("RostSpace", style={"color": "white"}),
                        dcc.Loading(
                            color="white",
                            style={},
                            children=[
                                html.Div(
                                    id="load_graph_spinner",
                                    hidden=True,
                                )
                            ],
                        ),
                    ],
                ),
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
                    get_help_button_tooltip(button_id="help_button"),
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
