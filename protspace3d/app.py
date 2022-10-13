#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
from pathlib import Path

from preprocessing import data_preprocessing
from visualization import init_app


def parse_args():
    """Creates and returns the ArgumentParser object

    Run example: python protspace3d/app.py -d data -b VA --sep , --uid_col 0
    """

    # Instantiate the parser
    parser = argparse.ArgumentParser(description="ProtSpace3D")  # TODO add description
    # Required argument
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
    # Required argument
    parser.add_argument(
        "-b",
        "--basename",
        required=True,
        type=str,
        help="Base filename for data in directory",
    )
    # Optional argument
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
        required=False,
        type=int,
        default=0,
        help="CSV column index containing the unique identifieres",
    )
    # Optional argument
    parser.add_argument(
        "--html_cols",
        required=False,
        type=int,
        help="CSV column data to be saved as html",
        nargs="+",
    )

    args = parser.parse_args()
    data_dir_path = Path(args.data)
    basename = Path(args.basename)
    csv_sep = args.sep
    uid_col = args.uid_col
    html_cols = args.html_cols

    return data_dir_path, basename, csv_sep, uid_col, html_cols


def main():
    # Parse arguments
    data_dir_path, basename, csv_sep, uid_col, html_cols = parse_args()

    # Preprocessing
    df, fig, csv_header = data_preprocessing(
        data_dir_path, basename, csv_sep, uid_col, html_cols
    )

    # --- APP creation ---
    application = init_app(fig, csv_header)

    return application, html_cols


if __name__ == "__main__":
    app, html = main()

    if html is None:
        app.run_server(debug=True)
