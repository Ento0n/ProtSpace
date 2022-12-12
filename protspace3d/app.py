#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import sys
from pathlib import Path

from preprocessing import DataPreprocessor
from visualization import Visualizator

from callbacks import get_callbacks, get_callbacks_pdb

import yaml


class LoadConfFile(argparse.Action):
    def __call__(self, parser, namespace, values, option_string=None):
        first_arg = sys.argv[1]

        if first_arg == "-conf" or first_arg == "--configuration":
            with open(values, "r") as f:
                dictionary = yaml.full_load(f)

            arguments = self.yaml_to_parser(dictionary)

            parser.parse_args(arguments, namespace)
        else:
            raise Exception(
                "Argument -conf or --configuration must be the first argument!"
            )

    @staticmethod
    def yaml_to_parser(dictionary: dict):
        arguments = list()
        if "o" in dictionary.keys():
            arguments.append("-o")
            arguments.append(str(dictionary["o"]))
        if "hdf" in dictionary.keys():
            arguments.append("--hdf")
            arguments.append(str(dictionary["hdf"]))
        if "csv" in dictionary.keys():
            arguments.append("--csv")
            arguments.append(str(dictionary["csv"]))
        if "sep" in dictionary.keys():
            arguments.append("--sep")
            arguments.append(str(dictionary["sep"]))
        if "uid_col" in dictionary.keys():
            arguments.append("--uid_col")
            arguments.append(str(dictionary["uid_col"]))
        if "html_cols" in dictionary.keys():
            arguments.append("--html_cols")
            arguments.append(str(dictionary["html_cols"]))
        if "pdb" in dictionary.keys():
            arguments.append("--pdb")
            arguments.append(str(dictionary["pdb"]))
        if "reset" in dictionary.keys():
            if dictionary["reset"]:
                arguments.append("--reset")
        if "n_neighbours" in dictionary.keys():
            arguments.append("--n_neighbours")
            arguments.append(str(dictionary["n_neighbours"]))
        if "min_dist" in dictionary.keys():
            arguments.append("--min_dist")
            arguments.append(str(dictionary["min_dist"]))
        if "metric" in dictionary.keys():
            arguments.append("--metric")
            arguments.append(str(dictionary["metric"]))

        return arguments


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
            self.conf,
            self.n_neighbours,
            self.min_dist,
            self.metric,
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
            self.conf,
            self.n_neighbours,
            self.min_dist,
            self.metric,
        )

    @staticmethod
    def _parse_args():
        """
        Creates and returns the ArgumentParser object

        Pla2g2:                 python protspace3d/app.py -o data/Pla2g2 --hdf data/Pla2g2/emb_prott5.h5 --csv data/Pla2g2/Pla2g2.csv --sep ";" --pdb data/Pla2g2/colabfold/pdb
                                python protspace3d/app.py -o data/Pla2g2 --hdf data/Pla2g2/emb_esm2.h5 --csv data/Pla2g2/Pla2g2.csv --sep ";" --html 1
        conotoxins_swiss_prot:  python protspace3d/app.py -o data/conotoxins_swiss_prot --hdf data/conotoxins_swiss_prot/Swiss_Prot_Conotoxins_mapped.h5 --csv data/conotoxins_swiss_prot/Swiss_Prot_Conotoxins_mapped.csv
        conotoxins:             python protspace3d/app.py -o data/conotoxins --hdf data/conotoxins/Conotoxins.h5 --csv data/conotoxins/Conotoxins_mapped.csv --pdb data/conotoxins/colabfold/pdb


        """

        # Instantiate the parser
        parser = argparse.ArgumentParser(description="ProtSpace3D")
        # Required argument
        parser.add_argument(
            "-o",
            "--output",
            required=False,
            type=str,
            help=(
                "Name of the output folder in the project directory where generated files will be stored."
            ),
        )
        # Required argument
        parser.add_argument(
            "--hdf",
            required=False,
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
            "-conf",
            "--configuration",
            required=False,
            type=str,
            action=LoadConfFile,
            help="Path to configuration file.",
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
        # Optional argument
        parser.add_argument(
            "--n_neighbours",
            required=False,
            default=25,
            type=int,
            help="UMAP parameter n_neighbours, default: 25",
        )
        parser.add_argument(
            "--min_dist",
            required=False,
            default=0.5,
            type=float,
            help="UMAP parameter min_dist, default: 0.5",
        )
        parser.add_argument(
            "--metric",
            required=False,
            default="euclidean",
            help="Metric used for UMAP calculation, default: euclidean",
        )

        args = parser.parse_args()
        output_d = Path(args.output) if args.output is not None else None
        hdf_path = Path(args.hdf) if args.hdf is not None else None
        csv_path = Path(args.csv) if args.csv is not None else None
        csv_sep = args.sep
        uid_col = args.uid_col
        html_cols = args.html_cols
        pdb_d = Path(args.pdb) if args.pdb is not None else None
        reset = args.reset
        conf_file = args.configuration
        n_neighbours = args.n_neighbours
        min_dist = args.min_dist
        metric = args.metric

        return (
            output_d,
            hdf_path,
            csv_path,
            csv_sep,
            uid_col,
            html_cols,
            pdb_d,
            reset,
            conf_file,
            n_neighbours,
            min_dist,
            metric,
        )


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
        conf,
        n_neighbours,
        min_dist,
        metric,
    ) = parser.get_params()

    # Check if h5 file is present
    if hdf_path is None or output_d is None:
        if hdf_path is None and output_d is None:
            raise Exception(
                "hdf path and output directory must be given either in config file or as argument!"
            )
        elif hdf_path is None:
            raise Exception(
                "hdf path must be given either in config file or as argument!"
            )
        else:
            raise Exception(
                "output directory must be given either in config file or as argument!"
            )

    # pdb directory given or not
    pdb_flag = True if pdb_d is not None else False

    # put UMAP parameters in dictionary
    umap_paras = dict()
    umap_paras["n_neighbours"] = n_neighbours
    umap_paras["min_dist"] = min_dist
    umap_paras["metric"] = metric

    # Create data preprocessor object
    data_preprocessor = DataPreprocessor(
        output_d,
        hdf_path,
        csv_path,
        csv_sep,
        uid_col,
        html_cols,
        reset,
        umap_paras,
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
        umap_paras,
        output_d,
    )


def main():
    app, html, df, struct_container, pdb, orig_id_col, umap_paras, output_d = setup()

    # don't start server if html is needed
    if not html:
        # different callbacks for different layout
        if pdb:
            get_callbacks(app, df, orig_id_col, umap_paras, output_d)
            get_callbacks_pdb(app, df, struct_container, orig_id_col)
        else:
            get_callbacks(app, df, orig_id_col, umap_paras, output_d)

        app.run_server(debug=True)


if __name__ == "__main__":
    main()
