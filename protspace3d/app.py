#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
from pathlib import Path

from preprocessing import DataPreprocessor
from visualization import Visualizator

from callbacks import get_callbacks

from dash import Input, Output
from dash.exceptions import PreventUpdate
import dash
from pandas import DataFrame
import dash_bio.utils.ngl_parser as ngl_parser
import pandas as pd


class Parser:
    def __init__(self):
        (
            self.data_dir_path,
            self.basename,
            self.csv_sep,
            self.uid_col,
            self.html_cols,
            self.pdb_flag,
        ) = self._parse_args()

    def get_params(self):
        """
        :return: class parameters
        """
        return (
            self.data_dir_path,
            self.basename,
            self.csv_sep,
            self.uid_col,
            self.html_cols,
            self.pdb_flag,
        )

    @staticmethod
    def _parse_args():
        """
        Creates and returns the ArgumentParser object

        Run example: python protspace3d/app.py -d data/ex1 -b VA --sep , --uid_col 0
        """

        # Instantiate the parser
        parser = argparse.ArgumentParser(description="ProtSpace3D")
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
            help="CSV column index containing the unique identifiers",
        )
        # Optional argument
        parser.add_argument(
            "--html_cols",
            required=False,
            type=int,
            help="CSV columns to be saved as html",
            nargs="+",
        )
        # Optional argument
        parser.add_argument(
            "--pdb",
            required=False,
            action="store_true",
            help="PDB flag, if set 3D structures can be viewed.",
        )

        args = parser.parse_args()
        data_dir_path = Path(args.data)
        basename = Path(args.basename)
        csv_sep = args.sep
        uid_col = args.uid_col
        html_cols = args.html_cols
        pdb_flag = args.pdb

        return data_dir_path, basename, csv_sep, uid_col, html_cols, pdb_flag


def main():
    """
    Handles the process of the application
    :return: app & html_flag
    """
    # Create Application object
    parser = Parser()

    # Parse arguments
    data_dir_path, basename, csv_sep, uid_col, html_cols, pdb_flag = parser.get_params()

    # Create data preprocessor object
    data_preprocessor = DataPreprocessor(
        data_dir_path, basename, csv_sep, uid_col, html_cols
    )

    # Preprocessing
    df, fig, csv_header = data_preprocessor.data_preprocessing()

    # initialize structure container if flag set
    struct_container = None
    if pdb_flag:
        struct_container = data_preprocessor.init_structure_container()

    # Create visualization object
    visualizator = Visualizator(fig, csv_header, df)

    # --- APP creation ---
    application = visualizator.init_app()

    # html flag set or not
    html_flag = True if html_cols is not None else False

    return application, html_flag, df, struct_container


if __name__ == "__main__":
    app, html, df, struct_container = main()

    get_callbacks(app, df, struct_container)

    if not html:
        app.run_server(debug=True)
