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

import dash_bio.utils.ngl_parser as ngl_parser

from visualization import Visualizator


class StructureContainer(object):
    def __init__(self, pdb_d):
        self.pdb_d = pdb_d
        self.point_number = None
        self.curve_number = None
        self.public_seq_id = None

    def __call__(self):
        return self.public_seq_id

    def set_focus_point(self, curve_number, point_number):
        self.curve_number = curve_number
        self.point_number = point_number
        return None

    def get_focus_point(self):
        return self.curve_number, self.point_number

    def get_structure_dir(self):
        return self.pdb_d

    def set_structure_ids(self, seq_ids):
        if isinstance(seq_ids, list):
            self.public_seq_id = [seq_id.replace(".", "_", 1) for seq_id in seq_ids]
        else:
            self.public_seq_id = [seq_ids.replace(".", "_", 1)]
        return None


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
        """
        reads & processes the files
        :return: dataframe wit collected information, graph & column headers
        """
        # root directory that holds, proteins.fasta, embeddings.h5, labels.csv and some output_file.html
        root = Path(self.data_dir_path)
        emb_h5file = root / f"{self.basename}.h5"
        label_csv_p = root / f"{self.basename}.csv"

        self._check_files(emb_h5file, label_csv_p)

        df_csv = pd.read_csv(
            label_csv_p, sep=self.csv_separator, index_col=self.uid_col
        )
        df_csv.fillna(" <NA> ", inplace=True)

        # save index name for df.csv
        index_name = df_csv.index.name

        # Unify notation of the UIDs
        csv_uids = self._unify_seq_uids(df_csv.index)
        df_csv.index = csv_uids

        df_embeddings, csv_header = self._read_df_csv(
            root, df_csv, emb_h5file, csv_uids, index_name
        )

        # save df_embeddings in global variable df
        global df
        df = df_embeddings

        # handle html saving
        DataPreprocessor._handle_html(self.html_cols, csv_header, self.data_dir_path)

        # generate initial figure
        fig = Visualizator.render(df, selected_column=csv_header[0])

        return df_embeddings, fig, csv_header

    @staticmethod
    def _check_files(emb_h5file: Path, label_csv_p: Path):
        """
        Checks whether all files are present
        :param emb_h5file: h5 file Path
        :param label_csv_p: csv file Path
        """
        # Check whether all files are present
        files = [emb_h5file, label_csv_p]
        for file in files:
            if not file.is_file():
                raise FileNotFoundError(
                    f"The {file} file is missing! Check data directory and basename."
                )

    @staticmethod
    def _handle_html(
        html_cols: list[int],
        csv_header: list[str],
        data_dir_path: Path,
    ):
        """
        Saves the html files as given by the --html_cols argument
        :param html_cols: List of given columns
        :param csv_header: List of dataframe headers
        :param data_dir_path: Path to data directory
        """
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

                # Edit input to required index numbers
                for i, num in enumerate(html_cols):
                    html_cols[i] = num - 1

                # Given parameters existing columns?
                for col in html_cols:
                    if col not in range(len(csv_header)):
                        raise Exception(f"Column no. {col} is not valid!")

                    fig = Visualizator.render(df, selected_column=csv_header[col])
                    fig.write_html(data_dir_path / f"3Dspace_{csv_header[col]}.html")

    def _read_df_csv(
        self,
        root: Path,
        df_csv: DataFrame,
        emb_h5file: Path,
        csv_uids: list[str],
        index_name: str,
    ):
        """
        If present, read df.csv and check for completion
        :param root: Path to data directory
        :param df_csv: dataframe of given csv file
        :param emb_h5file: Path to h5 file
        :param csv_uids: unique IDs of the csv file
        :param index_name: header of index row
        :return: final dataframe and its column headers
        """
        # load & read df.csv if present
        pres_df = root / "df.csv"
        if pres_df.is_file():
            print("Pre computed dataframe file df.csv is loaded.")
            pres_df_csv = pd.read_csv(pres_df, index_col=0)

            # Check whether no. of rows equals data
            if len(pres_df_csv) != len(df_csv):
                print("# of rows doesn't match data!\nStart recalculation!")
                df_embeddings, csv_header = self._create_df(
                    root, emb_h5file, csv_uids, df_csv, index_name
                )
            else:
                # Check each column x, y & z for incompleteness
                if not self._check_coordinates(pres_df_csv):
                    print("Start recalculation!")
                    df_embeddings, csv_header = self._create_df(
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
                        pres_df_csv = self._update_df(df_csv, pres_df_csv)

                        # save the new obtained df
                        pres_df_csv.to_csv(root / "df.csv", index_label=index_name)

                    # save final column names
                    csv_header = [
                        header
                        for header in pres_df_csv.columns
                        if header not in self.AXIS_NAMES or header is index_name
                    ]

                    # Unify df name
                    df_embeddings = pres_df_csv

        # create dataframe from files
        else:
            df_embeddings, csv_header = self._create_df(
                root, emb_h5file, csv_uids, df_csv, index_name
            )

        return df_embeddings, csv_header

    def _create_df(
        self,
        root: Path,
        emb_h5file: Path,
        csv_uids: list[str],
        df_csv: DataFrame,
        index_name: str,
    ):
        """
        Use data and create corresponding dataframe
        :param root: Path to data directory
        :param emb_h5file: Path to h5 file
        :param csv_uids: unique IDs of csv file
        :param df_csv: dataframe of csv file
        :param index_name: header of index column
        :return: complete dataframe and list of its headers
        """
        # read embeddings from HDF5 format
        embeddings = self._get_embeddings(emb_h5file=emb_h5file, csv_uids=csv_uids)

        # check for proteins in csv but not in h5 file
        self._check_csv_uids(embeddings=embeddings, csv_uids=csv_uids)

        # matrix of values (protein-embeddings); n_proteins x embedding_dim
        uids, embs = zip(*embeddings.items())
        embs = np.vstack(embs)

        # data should be n_proteins x 1024 (ProtT5) OR n_proteins x 128 (ProtTucker)
        print(f"Shape of embeddings (num_proteins x embedding dim): {embs.shape}")

        # get pairwise distances; should be n_proteins x n_proteins
        pairwise_dist = squareform(self._pairwise_distances(embs))
        print(
            "Shape of pairwise distance matrix (num_proteins x num_proteins):"
            f" {pairwise_dist.shape}"
        )

        # generate umap components and merge it to CSV DataFrame
        df_umap = self._generate_umap(embs)
        df_umap.index = uids

        df_embeddings = df_csv.join(df_umap, how="right")
        csv_header = [
            header for header in df_embeddings.columns if header not in self.AXIS_NAMES
        ]

        # save dataframe
        df_embeddings.to_csv(root / "df.csv", index_label=index_name)

        return df_embeddings, csv_header

    @staticmethod
    def _get_embeddings(emb_h5file: Path, csv_uids: list[str]) -> dict[str, np.ndarray]:
        """load pre-computed embeddings in .h5 file format

        Args:
            emb_h5file (str): path to hdf5 file containing embeddings
            csv_uids (list[str]): list of identifiers extracted from CSV file

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
    def _pairwise_distances(data, metric="euclidean"):
        """
        Calculate pairwise distance for given data
        :param data: embedding data
        :param metric: metric used for calculation
        :return: calculated pairwise distance
        """
        # usually euclidean or cosine distance worked best
        return pdist(data, metric=metric)

    def _unify_seq_uids(self, uids: list[str]) -> list[str]:
        """
        Replaces all special characters with _
        :param uids: unique IDs of the dataframe
        :return: new list of unique IDs
        """
        return list(map(lambda uid: self.NON_WORD_RE.sub("_", uid), uids))

    def _generate_umap(self, data: np.ndarray) -> pd.DataFrame:
        """
        generated umap for given data
        :param data: embeddings data
        :return: dataframe of the umap coordinates
        """
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

    def _check_coordinates(self, data_frame: DataFrame) -> bool:
        """
        Checks whether x, y & z column are present and complete in given dataframe
        :param data_frame: given dataframe
        :return: False if corrupted, True if not
        """
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

    def _update_df(self, df_csv: DataFrame, pres_df_csv: DataFrame):
        """
        new columns in data compared to present df.csv is added to df.csv
        :param df_csv: dataframe of data
        :param pres_df_csv: dataframe of df.csv
        :return: updated dataframe of df.csv
        """
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
    def _check_csv_uids(embeddings: dict[str, np.ndarray], csv_uids: list[str]):
        """
        Check unique IDs in csv but not in h5 file
        :param embeddings: data of the h5 file
        :param csv_uids: unique IDs of the csv file
        """
        missing = list()

        # iterate over csv uids
        for uid in csv_uids:
            if uid not in embeddings.keys():
                missing.append(uid)

        if (nr_missed := (len(missing))) > 0:
            print(f"{nr_missed} protein(s) in csv but not in h5 file:")
            print(", ".join(missing))

    def init_structure_container(self):
        root = Path(self.data_dir_path)

        structure_container = StructureContainer(root / "pdb")

        global struct_container
        struct_container = structure_container

    @staticmethod
    @dash.callback(
        Output("store_data", "data"),
        Output("dd_menu", "value"),
        Input("dd_menu", "value"),
    )
    def store_data(value):
        """
        Stores df to be used in other callback
        :param value: global value df
        :return: df in dictionary form
        """
        # Check whether an input is triggered
        ctx = dash.callback_context
        if not ctx.triggered:
            raise PreventUpdate

        return df.to_dict("records"), value

    @staticmethod
    @dash.callback(
        Output("ngl_molecule_viewer", "data"),
        Input("graph", "clickData"),
    )
    def display_molecule(clickdata):
        """
        callback function to handle the displaying of the molecule
        :param clickdata: given data by clicking on a datapoint in the 3D plot
        :return:
        """

        ctx = dash.callback_context
        if not ctx.triggered:
            raise PreventUpdate

        # dict with data of clickdata
        points = clickdata["points"][0]
        # class_index value and it's index number
        index_num = int(points["pointNumber"])
        class_index = points["curveNumber"]

        print(df)
        print(points)

        # extract df_row of selected protein
        class_df = df[df["class_index"] == class_index]
        print(f"class df: {class_df}")
        print(f"type: {type(class_df)}")
        print(f"index_num: {index_num}")
        df_row = class_df.iloc[index_num]
        # add missing name to series
        name = pd.Series(class_df.index[index_num])
        name.index = ["Name"]
        df_row = pd.concat([name, df_row])

        # extract sequence ID
        seq_id = df_row["Name"]

        # set structure container ID accordingly
        struct_container.set_structure_ids(seq_id)

        # path to .pdb file
        struct_path = str(struct_container.get_structure_dir()) + "/"

        # data format for molecule viewer
        data_list = ngl_parser.get_data(
            data_path=struct_path,
            pdb_id=seq_id,
            color="blue",
            reset_view=True,
            local=True,
        )

        return data_list
