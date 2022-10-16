#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from pathlib import Path
import h5py
from scipy.spatial.distance import pdist, squareform
import numpy as np
import pandas as pd
import re
from pandas import DataFrame

import dash
from dash import Input, Output
from dash.exceptions import PreventUpdate

from visualization import Visualizator


class DataPreprocessor:

    NON_WORD_RE = re.compile(r"\W")
    AXIS_NAMES = ["x", "y", "z"]

    def __init__(
        self,
        data_dir_path: Path,
        basename: Path,
        csv_separator: str,
        uid_col: int,
        html_cols: list[int],
    ):
        self.data_dir_path = data_dir_path
        self.basename = basename
        self.csv_separator = csv_separator
        self.uid_col = uid_col
        self.html_cols = html_cols

    def data_preprocessing(self):
        # root directory that holds, proteins.fasta, embeddings.h5, labels.csv and some output_file.html
        root = Path(self.data_dir_path)
        emb_h5file = root / f"{self.basename}.h5"
        label_csv_p = root / f"{self.basename}.csv"

        self.check_files(emb_h5file, label_csv_p)

        df_csv = pd.read_csv(
            label_csv_p, sep=self.csv_separator, index_col=self.uid_col
        )
        df_csv.fillna(" <NA> ", inplace=True)

        # save index name for df.csv
        index_name = df_csv.index.name

        # Unify notation of the UIDs
        csv_uids = self.unify_seq_uids(df_csv.index)
        df_csv.index = csv_uids

        df_embeddings, csv_header = self.read_df_csv(
            root, self.uid_col, df_csv, emb_h5file, csv_uids, index_name
        )

        # save df_embeddings in global variable df
        global df
        df = df_embeddings

        # generate initial figure
        fig = Visualizator.render(df, selected_column=csv_header[0])

        return df_embeddings, fig, csv_header

    @staticmethod
    def check_files(emb_h5file: Path, label_csv_p: Path):
        # Check whether all files are present
        files = [emb_h5file, label_csv_p]
        for file in files:
            if not file.is_file():
                raise FileNotFoundError(
                    f"The {file} file is missing! Check data directory and basename."
                )

    @staticmethod
    def handle_html(
        html_cols: list[int],
        csv_header: list[str],
        data_dir_path: Path,
    ):
        # save html figures if argument is set
        if html_cols is not None:
            # -1 indicates all columns to be saved
            if html_cols == [-1]:
                for col in csv_header:
                    fig = Visualizator.render(df, selected_column=col)
                    fig.write_html(data_dir_path / f"3Dspace_{col}.html")

            else:
                # Sort given column indexes
                html_cols = sorted(html_cols)

                # Given parameters existing columns?
                for col in html_cols:
                    if col not in range(len(csv_header)):
                        raise Exception(f"Column no. {col} is not valid!")

                    fig = Visualizator.render(df, selected_column=csv_header[col])
                    fig.write_html(data_dir_path / f"3Dspace_{csv_header[col]}.html")

    def read_df_csv(
        self,
        root: Path,
        uid_col: int,
        df_csv: DataFrame,
        emb_h5file: Path,
        csv_uids: list[str],
        index_name: str,
    ):
        # load & read df.csv if present
        pres_df = root / "df.csv"
        if pres_df.is_file():
            print("Pre computed dataframe file df.csv is loaded.")
            pres_df_csv = pd.read_csv(pres_df, index_col=uid_col)

            # Check whether no. of rows equals data
            if len(pres_df_csv) != len(df_csv):
                print("# of rows doesn't match data!\nStart recalculation!")
                df_embeddings, csv_header = self.create_df(
                    root, emb_h5file, csv_uids, df_csv, index_name
                )
            else:
                # Check each column x, y & z for incompleteness
                if not self.check_coordinates(pres_df_csv):
                    print("Start recalculation!")
                    df_embeddings, csv_header = self.create_df(
                        root, emb_h5file, csv_uids, df_csv, index_name
                    )
                # columns x, y & z are fine
                else:
                    # Update df in case new columns were added to the csv
                    if (
                        len(df_csv.columns)
                        - (len(pres_df_csv.columns) - len(self.AXIS_NAMES))
                        > 0
                    ):
                        print("New column(s) were found and will be added.")
                        pres_df_csv = self.update_df(df_csv, pres_df_csv)

                        # save the new obtained df
                        pres_df_csv.to_csv(root / "df.csv", index_label=index_name)

                    # save final column names
                    csv_header = [
                        header
                        for header in pres_df_csv.columns
                        if header not in self.AXIS_NAMES
                    ]

                    # Unify df name
                    df_embeddings = pres_df_csv

        # create dataframe from files
        else:
            df_embeddings, csv_header = self.create_df(
                root, emb_h5file, csv_uids, df_csv, index_name
            )

        return df_embeddings, csv_header

    def create_df(
        self,
        root: Path,
        emb_h5file: Path,
        csv_uids: list[str],
        df_csv: DataFrame,
        index_name: str,
    ):
        # read embeddings from HDF5 format
        embeddings = self.get_embeddings(emb_h5file=emb_h5file, csv_uids=csv_uids)

        # check for proteins in csv but not in h5 file
        self.check_csv_uids(embeddings=embeddings, csv_uids=csv_uids)

        # matrix of values (protein-embeddings); n_proteins x embedding_dim
        uids, embs = zip(*embeddings.items())
        embs = np.vstack(embs)

        # data should be n_proteins x 1024 (ProtT5) OR n_proteins x 128 (ProtTucker)
        print(f"Shape of embeddings (num_proteins x embedding dim): {embs.shape}")

        # get pairwise distances; should be n_proteins x n_proteins
        pairwise_dist = squareform(self.pairwise_distances(embs))
        print(
            "Shape of pairwise distance matrix (num_proteins x num_proteins):"
            f" {pairwise_dist.shape}"
        )

        # generate umap components and merge it to CSV DataFrame
        df_umap = self.generate_umap(embs)
        df_umap.index = uids

        df_embeddings = df_csv.join(df_umap, how="right")
        csv_header = [
            header for header in df_embeddings.columns if header not in self.AXIS_NAMES
        ]

        # save dataframe
        df_embeddings.to_csv(root / "df.csv", index_label=index_name)

        return df_embeddings, csv_header

    @staticmethod
    def get_embeddings(emb_h5file: Path, csv_uids: list[str]) -> dict[str, np.ndarray]:
        """load pre-computed embeddings in .h5 file format

        Args:
            emb_h5file (str): path to hdf5 file containing embeddings
            csv_uids (list[str]): list of identifieres extracted from CSV file

        Returns:
            dict[str, np.ndarray]: dictionary with fasta headers as keys and a
                single vector (embeddings) per protein as value. Values have
                1024-dimensions for ProtT5 and 128-dimensions for ProtTucker
        """

        embeddings = dict()
        missing = list()
        print(f"Loading pre-computed embeddings from: {emb_h5file}")
        with h5py.File(emb_h5file, "r") as hdf:
            for identifier, embd in hdf.items():
                if identifier in csv_uids:
                    embeddings[identifier] = embd[:]
                else:
                    missing.append(identifier)

        # Check whether any UIDs matched with the embedding.keys
        if len(embeddings.keys()) == 0:
            raise Exception(
                "None of the Unique IDs of the h5 and the csv file matched."
            )

        print(f"Example: {next(iter(embeddings.keys()))}")
        print(f"Number of embeddings: {len(embeddings)}")
        if (nr_missed := len(missing)) > 0:
            print(f"{nr_missed} protein(s) in h5 but not in csv file:")
            print(", ".join(missing))

        return embeddings

    @staticmethod
    def pairwise_distances(data, metric="euclidean"):
        # usually euclidean or cosine distance worked best
        return pdist(data, metric=metric)

    def unify_seq_uids(self, uids: list[str]) -> list[str]:
        return list(map(lambda uid: self.NON_WORD_RE.sub("_", uid), uids))

    def generate_umap(self, data: np.ndarray) -> pd.DataFrame:
        # visualize high-dimensional embeddings with dimensionality reduction (here: umap)
        # Tutorial: https://umap-learn.readthedocs.io/en/latest/basic_usage.html
        # Parameters: https://umap-learn.readthedocs.io/en/latest/parameters.html
        import umap

        fit = umap.UMAP(
            n_neighbors=25, min_dist=0.5, random_state=42, n_components=3
        )  # initialize umap; use random_state=42 for reproducibility
        umap_fit = fit.fit_transform(data)  # fit umap to our embeddings
        df_umap = DataFrame(data=umap_fit, columns=self.AXIS_NAMES)
        return df_umap

    def check_coordinates(self, data_frame: DataFrame) -> bool:
        # Do the columns x, y and z exist?
        if not all(x in list(data_frame.columns) for x in self.AXIS_NAMES):
            print("Df corrupted as not x,y & z columns are present!")
            return False

        # Is the corresponding data complete ?
        for col in self.AXIS_NAMES:
            for value in list(data_frame[col]):
                if not isinstance(value, float):
                    # Value is corrupted
                    print(f"At least one value of the {col} column is corrupted!")
                    return False

        # All values of the x,y & z column are correct
        return True

    def update_df(self, df_csv: DataFrame, pres_df_csv: DataFrame):
        # extract column names
        df_cols = set(df_csv.columns)
        pres_df_cols = set(pres_df_csv.columns)

        # get missing columns in present df
        missing_cols = df_cols - pres_df_cols

        # add missing columns to the present df
        for col in missing_cols:
            pres_df_csv.insert(
                len(pres_df_cols) - len(self.AXIS_NAMES), col, df_csv[col]
            )
            print(
                f"Missing column {col} from the .csv file has been added to the present df.csv file."
            )

        # return updated df
        return pres_df_csv

    @staticmethod
    def check_csv_uids(embeddings: dict[str, np.ndarray], csv_uids: list[str]):
        missing = list()

        # iterate over csv uids
        for uid in csv_uids:
            if uid not in embeddings.keys():
                missing.append(uid)

        if (nr_missed := (len(missing))) > 0:
            print(f"{nr_missed} protein(s) in csv but not in h5 file:")
            print(", ".join(missing))

    @staticmethod
    @dash.callback(
        Output("store_data", "data"),
        Output("dd_menu", "value"),
        Input("dd_menu", "value"),
    )
    def store_data(value):
        # Check whether an input is triggered
        ctx = dash.callback_context
        if not ctx.triggered:
            raise PreventUpdate

        return df.to_dict("records"), value
