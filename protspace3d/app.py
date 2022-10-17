#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
from pathlib import Path

from preprocessing import DataPreprocessor
from visualization import Visualizator


class Parser:
    def __init__(self):
        (
            self.data_dir_path,
            self.basename,
            self.csv_sep,
            self.uid_col,
            self.html_cols,
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

        args = parser.parse_args()
        data_dir_path = Path(args.data)
        basename = Path(args.basename)
        csv_sep = args.sep
        uid_col = args.uid_col
        html_cols = args.html_cols

        return data_dir_path, basename, csv_sep, uid_col, html_cols


def main():
    """
    Handles the process of the application
    :return: app & html_flag
    """
    # Create Application object
    parser = Parser()

    # Parse arguments
    data_dir_path, basename, csv_sep, uid_col, html_cols = parser.get_params()

    # Create data preprocessor object
    data_preprocessor = DataPreprocessor(
        data_dir_path, basename, csv_sep, uid_col, html_cols
    )

    # Preprocessing
    df, fig, csv_header = data_preprocessor.data_preprocessing()

    # Create visualization object
    visualizator = Visualizator(fig, csv_header)

    # --- APP creation ---
    application = visualizator.init_app()

    # html flag set or not
    html_flag = True if html_cols is not None else False

    return application, html_flag


if __name__ == "__main__":
    app, html = main()

    if not html:
        app.run_server(debug=True)
