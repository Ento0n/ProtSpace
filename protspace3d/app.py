#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
from pathlib import Path

from preprocessing import DataPreprocessor
from visualization import Visualizator

from callbacks import get_callbacks, get_callbacks_pdb


class Parser:
    def __init__(self):
        (
            self.output_d,
            self.hdf_path,
            self.csv_path,
            self.csv_sep,
            self.uid_col,
            self.html_cols,
            self.pdb_flag,
            self.reset,
        ) = self._parse_args()

    def get_params(self):
        """
        :return: class parameters
        """
        return (
            self.output_d,
            self.hdf_path,
            self.csv_path,
            self.csv_sep,
            self.uid_col,
            self.html_cols,
            self.pdb_flag,
            self.reset,
        )

    @staticmethod
    def _parse_args():
        """
        Creates and returns the ArgumentParser object

        Pla2g2:                 python protspace3d/app.py -o data/Pla2g2 --hdf data/Pla2g2/emb_prott5.h5 --csv data/Pla2g2/Pla2g2.csv --sep ";" --pdb data/Pla2g2/colabfold/pdb
                                python protspace3d/app.py -o data/Pla2g2 --hdf data/Pla2g2/emb_esm2.h5 --csv data/Pla2g2/Pla2g2.csv --sep ";" --html 1
        conotoxins_swiss_prot:  python protspace3d/app.py -o data/conotoxins_swiss_prot --hdf data/conotoxins_swiss_prot/Swiss_Prot_Conotoxins_mapped.h5 --csv data/conotoxins_swiss_prot/Swiss_Prot_Conotoxins_mapped.csv
        conotoxins:             python protspace3d/app.py -o data/conotoxins --hdf data/conotoxins/Conotoxins_mapped.h5 --csv data/conotoxins/Conotoxins_mapped.csv --pdb data/conotoxins/colabfold/pdb


        """

        # Instantiate the parser
        parser = argparse.ArgumentParser(description="ProtSpace3D")
        # Required argument
        parser.add_argument(
            "-o",
            "--output",
            required=True,
            type=str,
            help=(
                "Name of the output folder in the project directory where generated files will be stored."
            ),
        )
        # Required argument
        parser.add_argument(
            "--hdf",
            required=True,
            type=str,
            help="Path to HDF5-file containing the per protein embeddings as a key-value pair, aka UID-embedding",
        )
        # Required argument
        parser.add_argument(
            "--csv",
            required=False,
            type=str,
            help="Path to CSV-file containg groups/features by which the dots in the 3D-plot are colored",
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
            type=str,
            help="Path the directory that holds the .pdb files.",
        )
        # Optional argument
        parser.add_argument(
            "--reset",
            required=False,
            action="store_true",
            help="Generated df.csv is deleted and recalculated.",
        )

        args = parser.parse_args()
        output_d = Path(args.output)
        hdf_path = Path(args.hdf)
        csv_path = args.csv
        csv_sep = args.sep
        uid_col = args.uid_col
        html_cols = args.html_cols
        pdb_d = args.pdb
        reset = args.reset

        return output_d, hdf_path, csv_path, csv_sep, uid_col, html_cols, pdb_d, reset


def setup():
    """
    Handles the process of the application
    :return: app & html_flag
    """
    # Create Application object
    parser = Parser()

    # Parse arguments
    (
        output_d,
        hdf_path,
        csv_path,
        csv_sep,
        uid_col,
        html_cols,
        pdb_d,
        reset,
    ) = parser.get_params()

    # Set up csv argument in path, even if not set
    if csv_path is not None:
        csv_path = Path(csv_path)
    else:
        csv_path = Path("")

    # pdb directory given or not
    pdb_flag = True if pdb_d is not None else False

    # Create data preprocessor object
    data_preprocessor = DataPreprocessor(
        output_d, hdf_path, csv_path, csv_sep, uid_col, html_cols, reset
    )

    # Preprocessing
    (
        df,
        fig,
        csv_header,
        original_id_col,
    ) = data_preprocessor.data_preprocessing()

    # initialize structure container if flag set
    structure_container = None
    if pdb_flag:
        pdb_d = Path(pdb_d)
        structure_container = data_preprocessor.init_structure_container(pdb_d)

    # Create visualization object
    visualizator = Visualizator(fig, csv_header)

    # get ids of the proteins
    if original_id_col is not None:
        ids = original_id_col
    else:
        ids = df.index.to_list()

    # --- APP creation ---
    if pdb_flag:
        application = visualizator.init_app_pdb(ids)
    else:
        application = visualizator.init_app()

    # html cols set or not
    html_flag = True if html_cols is not None else False

    return (
        application,
        html_flag,
        df,
        structure_container,
        pdb_flag,
        original_id_col,
    )


def main():
    app, html, df, struct_container, pdb, orig_id_col = setup()

    # don't start server if html is needed
    if not html:
        # different callbacks for different layout
        if pdb:
            get_callbacks(app, df, orig_id_col)
            get_callbacks_pdb(app, df, struct_container, orig_id_col)
        else:
            get_callbacks(app, df, orig_id_col)

        app.run_server(debug=True)


if __name__ == "__main__":
    main()
