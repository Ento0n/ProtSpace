#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
from pathlib import Path
from dash import Dash, dcc, html

import dash
from dash import Input, Output
from dash.exceptions import PreventUpdate

from preprocessing import data_preprocessing
from visualization import render, init_app

@dash.callback(Output("graph", "figure"), Input("dd_menu", "value"))
def update_graph(xaxis_column_name):
    ctx = dash.callback_context
    if not ctx.triggered:
        raise PreventUpdate

    fig = render(df, selected_column=xaxis_column_name)
    print("Update call!")
    fig.update_traces(hoverinfo="none", hovertemplate=None)
    fig.update_layout(clickmode="event+select")
    fig.update_traces(
        marker=dict(size=6, line=dict(width=1, color="DarkSlateGrey")),
        selector=dict(mode="markers"),
    )
    return fig


def parse_args():
    """Creates and returns the ArgumentParser object

    Run example: python processing.py -d data -b VA --sep , --uid_col 0
    """

    # Instantiate the parser
    parser = argparse.ArgumentParser(
        description="ProtSpace3D"
    )  # TODO add description
    # Required positional argument
    parser.add_argument(
        "-d",
        "--data",
        required=True,
        type=str,
        help=(
            "A path to the data directory containing containing as a stem"
            " basename with different extensions for different files (.csv,"
            " .h5, .fasta)"
        ),
    )
    # Required positional argument
    parser.add_argument(
        "-b",
        "--basename",
        required=True,
        type=str,
        help="Base filename for data in directory",
    )
    # Required positional argument
    parser.add_argument(
        "--sep",
        required=False,
        type=str,
        default=",",
        help="Separator for CSV file",
    )
    # Optional argument
    parser.add_argument(
        "--uid_col",
        type=int,
        default=0,
        help="CSV column index containing the unique identifieres",
    )

    args = parser.parse_args()
    data_dir_path = Path(args.data)
    basename = Path(args.basename)
    csv_sep = args.sep
    uid_col = args.uid_col

    return data_dir_path, basename, csv_sep, uid_col


# Parse arguments
data_dir_path, basename, csv_sep, uid_col = parse_args()


# Preprocessing
df, fig, csv_header = data_preprocessing(data_dir_path, basename, csv_sep, uid_col)


# --- APP creation ---
# TODO separate code in a different script
app = init_app(fig, csv_header)


if __name__ == "__main__":
    app.run_server(debug=True)
