#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import sys
from pathlib import Path

import yaml

from src.callbacks import get_callbacks, get_callbacks_pdb
from src.preprocessing import DataPreprocessor
from src.structurecontainer import StructureContainer
from src.visualization.visualizator import Visualizator


class LoadConfFile(argparse.Action):
    def __call__(
        self,
        parser: argparse.ArgumentParser,
        namespace: argparse.Namespace,
        values: str,
        option_string: str = None,
    ):
        """
        Handling of conf file
        :param parser: the argument parser itself
        :param namespace: argument parser specific variable
        :param values: value of conf argument
        :param option_string: argument parser specific variable
        :return: arguments in parser list format
        """
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
        """
        Transform the arguments from the yaml file format to the argument parser format
        :param dictionary: arguments given in conf file
        :return: arguemnts in parser format
        """
        arguments = list()
        if "o" in dictionary.keys():
            arguments.append("-o")
            arguments.append(str(dictionary["o"]))
        if "output" in dictionary.keys():
            arguments.append("--output")
            arguments.append(str(dictionary["output"]))
        if "hdf" in dictionary.keys():
            arguments.append("--hdf")
            arguments.append(str(dictionary["hdf"]))
        if "csv" in dictionary.keys():
            arguments.append("--csv")
            arguments.append(str(dictionary["csv"]))
        if "fasta" in dictionary.keys():
            arguments.append("--fasta")
            arguments.append(str(dictionary["fasta"]))
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
        if "json" in dictionary.keys():
            arguments.append("--json")
            arguments.append(str(dictionary["json"]))
        if "reset" in dictionary.keys():
            if dictionary["reset"]:
                arguments.append("--reset")
        if "pca" in dictionary.keys():
            arguments.append("--pca")
        if "tsne" in dictionary.keys():
            arguments.append("--tsne")
        if "iterations" in dictionary.keys():
            arguments.append("--iterations")
            arguments.append(str(dictionary["iterations"]))
        if "perplexity" in dictionary.keys():
            arguments.append("--perplexity")
            arguments.append(str(dictionary["perplexity"]))
        if "learning_rate" in dictionary.keys():
            arguments.append("--learning_rate")
            arguments.append(str(dictionary["learning_rate"]))
        if "tsne_metric" in dictionary.keys():
            arguments.append("--tsne_metric")
            arguments.append(str(dictionary["tsne_metric"]))
        if "n_neighbours" in dictionary.keys():
            arguments.append("--n_neighbours")
            arguments.append(str(dictionary["n_neighbours"]))
        if "min_dist" in dictionary.keys():
            arguments.append("--min_dist")
            arguments.append(str(dictionary["min_dist"]))
        if "metric" in dictionary.keys():
            arguments.append("--metric")
            arguments.append(str(dictionary["metric"]))
        if "port" in dictionary.keys():
            arguments.append("--port")
            arguments.append(str(dictionary["port"]))
        if "verbose" in dictionary.keys():
            arguments.append("--verbose")
            arguments.append(str(dictionary["verbose"]))
        if "v" in dictionary.keys():
            arguments.append("-v")

        return arguments


class Parser:
    def __init__(self):
        (
            self.output_d,
            self.hdf_path,
            self.csv_path,
            self.fasta_path,
            self.csv_sep,
            self.uid_col,
            self.html_cols,
            self.pdb_d,
            self.json_d,
            self.reset,
            self.conf,
            self.pca_flag,
            self.tsne_flag,
            self.iterations,
            self.perplexity,
            self.learning_rate,
            self.tsne_metric,
            self.n_neighbours,
            self.min_dist,
            self.metric,
            self.port,
            self.verbose,
        ) = self._parse_args()

    def get_params(self):
        """
        :return: class parameters
        """
        return (
            self.output_d,
            self.hdf_path,
            self.csv_path,
            self.fasta_path,
            self.csv_sep,
            self.uid_col,
            self.html_cols,
            self.pdb_d,
            self.json_d,
            self.reset,
            self.conf,
            self.pca_flag,
            self.tsne_flag,
            self.iterations,
            self.perplexity,
            self.learning_rate,
            self.tsne_metric,
            self.n_neighbours,
            self.min_dist,
            self.metric,
            self.port,
            self.verbose,
        )

    @staticmethod
    def _parse_args():
        """
        Creates and returns the ArgumentParser object
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
                "Name of the output folder where generated files will be"
                " stored."
            ),
        )
        # Required argument
        parser.add_argument(
            "--hdf",
            required=False,
            type=str,
            help=(
                "Path to HDF5-file containing the per protein embeddings as a"
                " key-value pair, aka UID-embedding"
            ),
        )
        # Required argument
        parser.add_argument(
            "--csv",
            required=False,
            type=str,
            help=(
                "Path to CSV-file containing groups/features by which the dots"
                " in the 3D-plot are colored"
            ),
        )
        parser.add_argument(
            "-f",
            "--fasta",
            required=False,
            type=str,
            help=(
                "Path to fasta file containing the ID and the according"
                " sequence."
            ),
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
            help="Separator of CSV file",
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
            help=(
                "CSV columns to be saved as html, either the column index(es)"
                " or column name(s)."
            ),
            nargs="+",
        )
        # Optional argument
        parser.add_argument(
            "--pdb",
            required=False,
            type=str,
            help="Path the directory that holds the .pdb files.",
        )
        parser.add_argument(
            "--json",
            required=False,
            type=str,
            help="Path the directory that holds the .json files.",
        )
        # Optional argument
        parser.add_argument(
            "--reset",
            required=False,
            action="store_true",
            help="Precomputed file df.csv is deleted and recalculated.",
        )
        # Optional argument
        parser.add_argument(
            "--pca",
            required=False,
            action="store_true",
            help="PCA is initially used as dimensionality reduction",
        )
        # Optional argument
        parser.add_argument(
            "--tsne",
            required=False,
            action="store_true",
            help="t-SNE is initially used as dimensionality reduction",
        )
        # Optional argument
        parser.add_argument(
            "--iterations",
            required=False,
            default=1000,
            type=int,
            help="t-SNE parameter n_iters, default: 1000",
        )
        # Optional argument
        parser.add_argument(
            "--perplexity",
            required=False,
            default=30.0,
            type=float,
            help="t-SNE parameter perplexity, default: 30.0",
        )
        # Optional argument
        parser.add_argument(
            "--learning_rate",
            required=False,
            default="auto",
            help="t-SNE parameter learning_rate, default: auto",
        )
        # Optional argument
        parser.add_argument(
            "--tsne_metric",
            required=False,
            default="euclidean",
            help="t-SNE parameter metric, default: euclidean",
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
        parser.add_argument(
            "--port",
            required=False,
            type=int,
            default=8050,
            help="Port on which the website is locally hosted.",
        )
        parser.add_argument(
            "-v",
            "--verbose",
            required=False,
            action="store_true",
            help=(
                "By setting the verbose parameter, the script prints its"
                " internal operations."
            ),
        )

        args = parser.parse_args()
        output_d = Path(args.output) if args.output is not None else None
        hdf_path = Path(args.hdf) if args.hdf is not None else None
        csv_path = Path(args.csv) if args.csv is not None else None
        fasta_path = Path(args.fasta) if args.fasta is not None else None
        csv_sep = args.sep
        uid_col = args.uid_col
        html_cols = args.html_cols
        pdb_d = Path(args.pdb) if args.pdb is not None else None
        json_d = Path(args.json) if args.json is not None else None
        reset = args.reset
        conf_file = args.configuration
        pca_flag = args.pca
        tsne_flag = args.tsne
        iterations = args.iterations
        perplexity = args.perplexity
        learning_rate = args.learning_rate
        tsne_metric = args.tsne_metric
        n_neighbours = args.n_neighbours
        min_dist = args.min_dist
        metric = args.metric
        port = args.port
        verbose = args.verbose

        return (
            output_d,
            hdf_path,
            csv_path,
            fasta_path,
            csv_sep,
            uid_col,
            html_cols,
            pdb_d,
            json_d,
            reset,
            conf_file,
            pca_flag,
            tsne_flag,
            iterations,
            perplexity,
            learning_rate,
            tsne_metric,
            n_neighbours,
            min_dist,
            metric,
            port,
            verbose,
        )


def required_arguments_check(hdf_path: Path, output_d: Path):
    """
    Checks whether the required arguments for the script are given after processing conf file and arguments
    :param hdf_path:    Path to the hdf file
    :param output_d:    Path to the output directory
    :return: None
    """
    if hdf_path is None or output_d is None:
        if hdf_path is None and output_d is None:
            raise Exception(
                "hdf path and output directory must be given either in config"
                " file or as argument!"
            )
        elif hdf_path is None:
            raise Exception(
                "hdf path must be given either in config file or as argument!"
            )
        else:
            raise Exception(
                "output directory must be given either in config file or as"
                " argument!"
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
        fasta_path,
        csv_sep,
        uid_col,
        html_cols,
        pdb_d,
        json_d,
        reset,
        conf,
        pca_flag,
        tsne_flag,
        iterations,
        perplexity,
        learning_rate,
        tsne_metric,
        n_neighbours,
        min_dist,
        metric,
        port,
        verbose,
    ) = parser.get_params()

    required_arguments_check(hdf_path, output_d)

    dim_red = "UMAP"
    if pca_flag:
        dim_red = "PCA"
    elif tsne_flag:
        dim_red = "TSNE"

    # put UMAP parameters in dictionary
    umap_paras = dict()
    umap_paras["n_neighbours"] = n_neighbours
    umap_paras["min_dist"] = min_dist
    umap_paras["metric"] = metric

    # Put TSNE parameters in dictionary
    tsne_paras = dict()
    tsne_paras["iterations"] = iterations
    tsne_paras["perplexity"] = perplexity
    tsne_paras["learning_rate"] = learning_rate
    tsne_paras["tsne_metric"] = tsne_metric

    # Create data preprocessor object
    data_preprocessor = DataPreprocessor(
        output_d,
        hdf_path,
        csv_path,
        fasta_path,
        csv_sep,
        uid_col,
        html_cols,
        reset,
        dim_red,
        umap_paras,
        tsne_paras,
        verbose,
    )

    # Preprocessing
    (
        df,
        fig,
        csv_header,
        original_id_col,
        embeddings,
        embedding_uids,
        distance_dic,
        fasta_dict,
    ) = data_preprocessor.data_preprocessing()

    # initialize structure container if flag set
    structure_container = StructureContainer(pdb_d, json_d)

    # Create visualization object
    visualizator = Visualizator(fig, csv_header, dim_red)

    # get ids of the proteins
    if original_id_col is not None:
        ids = original_id_col
    else:
        ids = df.index.to_list()

    umap_paras_dict = data_preprocessor.get_umap_paras_dict(df)
    tsne_paras_dict = data_preprocessor.get_tsne_paras_dict(df)

    # --- APP creation ---
    if structure_container.pdb_flag:
        application = visualizator.get_pdb_app(ids, umap_paras, tsne_paras)
    else:
        application = visualizator.get_base_app(umap_paras, tsne_paras, ids)

    return (
        application,
        True if html_cols is not None else False,
        df,
        structure_container,
        original_id_col,
        dim_red,
        umap_paras,
        umap_paras_dict,
        tsne_paras,
        tsne_paras_dict,
        output_d,
        csv_header,
        port,
        embeddings,
        embedding_uids,
        distance_dic,
        fasta_dict,
    )


def main():
    """
    Most general processing of the script
    :return: None
    """
    (
        app,
        html,
        df,
        struct_container,
        orig_id_col,
        dim_red,
        umap_paras,
        umap_paras_dict,
        tsne_paras,
        tsne_paras_dict,
        output_d,
        csv_header,
        port,
        embeddings,
        embedding_uids,
        distance_dic,
        fasta_dict,
    ) = setup()

    # don't start server if html is needed
    if not html:
        # different callbacks for different layout
        if struct_container.pdb_flag:
            get_callbacks(
                app,
                df,
                orig_id_col,
                umap_paras,
                tsne_paras,
                output_d,
                csv_header,
                embeddings,
                embedding_uids,
                distance_dic,
                umap_paras_dict,
                tsne_paras_dict,
                fasta_dict,
                struct_container,
            )
            get_callbacks_pdb(app, df, struct_container, orig_id_col)
        else:
            get_callbacks(
                app,
                df,
                orig_id_col,
                umap_paras,
                tsne_paras,
                output_d,
                csv_header,
                embeddings,
                embedding_uids,
                distance_dic,
                umap_paras_dict,
                tsne_paras_dict,
                fasta_dict,
                struct_container,
            )

        app.run_server(debug=True, port=port)


if __name__ == "__main__":
    main()
