#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
from pathlib import Path
import h5py
import pandas
from scipy.spatial.distance import pdist, squareform
import numpy as np
import pandas as pd
from pandas import DataFrame

import re

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

    def get_range(self, uid: str):
        # add .pdb file type to ID
        uid = uid + ".pdb"

        mol_range = set()
        strand = None
        with open(self.pdb_d / uid, "r") as f:
            lines = f.readlines()

            for line in lines:
                if line.startswith("ATOM"):
                    pieces = re.split("\\s+", line)

                    mol_range.add(int(pieces[5]))
                    strand = pieces[4]

        return sorted(list(mol_range)), strand


class DataPreprocessor:
    AXIS_NAMES = ["x", "y", "z"]

    def __init__(
        self,
        output_d: Path,
        hdf_path: Path,
        csv_path: Path,
        csv_separator: str,
        uid_col: int,
        html_cols: list[int],
        reset: bool,
        umap_flag: bool,
        umap_paras: dict,
    ):
        self.output_d = output_d
        self.hdf_path = hdf_path
        self.csv_path = csv_path
        self.csv_separator = csv_separator
        self.uid_col = uid_col
        self.html_cols = html_cols
        self.reset = reset
        self.umap_flag = umap_flag
        self.umap_paras = umap_paras

    def data_preprocessing(self):
        """
        reads & processes the files
        :return: dataframe wit collected information, graph & column headers
        """
        # root directory that holds, proteins.fasta, embeddings.h5, labels.csv and some output_file.html
        emb_h5file = self.hdf_path
        label_csv_p = self.csv_path

        # delete df.csv
        if self.reset:
            df_csv_path = self.output_d / "df.csv"
            if df_csv_path.is_file():
                os.remove(df_csv_path)

        csv_less_flag = self._check_files(emb_h5file, label_csv_p, self.csv_separator)

        mapped_flag = self._check_mapped(label_csv_p, self.csv_separator)

        # processed by fasta_mapper.py?
        original_id_col = None
        if mapped_flag and csv_less_flag is False:
            # UID col from fasta_mapper.py is always 0
            df_csv = pd.read_csv(label_csv_p, index_col=0)

            # Extract original ID column
            original_id_col = df_csv["original_id"].to_list()
            df_csv.drop(columns=["original_id"], inplace=True)
        # no csv given or csv only exists for mapping
        elif csv_less_flag:
            df_csv, original_id_col = self._create_csv_less_df(
                emb_h5file, mapped_flag, label_csv_p
            )
        # not mapped and csv given
        else:
            df_csv = pd.read_csv(
                label_csv_p, sep=self.csv_separator, index_col=self.uid_col
            )

        # replace empty values with NA
        df_csv.fillna("NA", inplace=True)

        # save index name for df.csv
        index_name = df_csv.index.name

        # get UIDs
        csv_uids = df_csv.index.to_list()

        df_embeddings, csv_header = self._read_df_csv(
            self.output_d,
            df_csv,
            emb_h5file,
            csv_uids,
            index_name,
            self.umap_flag,
        )

        # sort csv header alphabetically
        csv_header.sort(key=str.lower)

        # handle html saving
        DataPreprocessor._handle_html(
            self,
            self.html_cols,
            csv_header,
            self.output_d,
            original_id_col=original_id_col,
            df=df_embeddings,
        )

        # generate initial figure
        fig = Visualizator.render(
            df_embeddings,
            selected_column=csv_header[0],
            original_id_col=original_id_col,
            umap_flag=self.umap_flag,
        )

        return df_embeddings, fig, csv_header, original_id_col

    def _check_files(
        self,
        hdf_path: Path,
        csv_path: Path,
        separator: str,
    ):
        """
        Checks whether all files are present
        :param hdf_path: h5 file Path
        :param csv_path: csv file Path
        :return boolean flag whether csv less mode is activated or not
        """
        # Check whether the h5 files is present (mandatory)
        if not hdf_path.is_file():
            raise FileNotFoundError(
                f"The {hdf_path} file is missing! Check data directory and basename."
            )

        # Check whether csv file is present (csv less mode activated if not)
        csv_less_flag = False
        if not csv_path.is_file():
            csv_less_flag = True

            print(
                "No csv file found!\nActivate csv less mode, no groups or features are visualized."
            )
        # get csv headers to check whether csv only exists for mapping
        else:
            headers = self._get_headers(csv_path, separator)
            mapped_headers = ["mapped_id", "original_id"]

            mapped_headers.sort()
            headers.sort()

            if headers == mapped_headers:
                csv_less_flag = True

                print("No groups/features in csv file!")

        return csv_less_flag

    def _check_mapped(self, csv_path: Path, separator: str):
        mapped_flag = False
        if csv_path.is_file():
            # get headers of csv file and create headers to be expected by a mapped file
            headers = self._get_headers(csv_path, separator)
            mapped_headers = ["mapped_id", "original_id"]

            # sort the headers that comparison is possible
            headers.sort()
            mapped_headers.sort()

            check = all(item in headers for item in mapped_headers)
            if check or headers == mapped_headers:
                mapped_flag = True

        return mapped_flag

    @staticmethod
    def _get_headers(csv_path: Path, separator: str):
        with open(csv_path, "r", encoding="utf-8-sig") as f:
            first_line = f.readline()
            headers = first_line.split(separator)

        last_item = headers.pop()
        headers.append(last_item.strip())

        return headers

    @staticmethod
    def _create_csv_less_df(hdf_path: Path, mapped_flag: bool, csv_path: Path):
        original_id_col = None
        if mapped_flag:
            # read csv file
            df = pd.read_csv(csv_path, sep=",", index_col=0)

            # extract original ID column
            original_id_col = df["original_id"].to_list()
            df.drop(columns=["original_id"], inplace=True)

            # create empty column no_group
            df["no_group"] = ""
        else:
            # get all ids of the h5 file
            h5_uids = list()
            with h5py.File(hdf_path, "r") as hdf:
                for identifier, embd in hdf.items():
                    h5_uids.append(identifier)

            # create new dataframe with collected identifiers
            df = pd.DataFrame(index=h5_uids, columns=["no group"])
            df.index.name = "ID"

        return df, original_id_col

    def _handle_html(
        self,
        html_cols: list[int],
        csv_header: list[str],
        output_d: Path,
        original_id_col: list,
        df: DataFrame,
    ):
        """
        Saves the html files as given by the --html_cols argument
        :param html_cols: List of given columns
        :param csv_header: List of dataframe headers
        :param output_d: Path to data directory
        """
        # save html figures if argument is set
        if html_cols is not None:
            # -1 indicates all columns to be saved
            if html_cols == [-1]:
                for col in csv_header:
                    fig = Visualizator.render(
                        df,
                        selected_column=col,
                        original_id_col=original_id_col,
                        umap_flag=self.umap_flag,
                    )
                    fig.write_html(output_d / f"3Dspace_{col}.html")

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

                    fig = Visualizator.render(
                        df,
                        selected_column=csv_header[col],
                        original_id_col=original_id_col,
                        umap_flag=self.umap_flag,
                    )
                    fig.write_html(output_d / f"3Dspace_{csv_header[col]}.html")

    def _read_df_csv(
        self,
        output_d: Path,
        df_csv: DataFrame,
        hdf_path: Path,
        csv_uids: list[str],
        index_name: str,
        umap_flag: bool,
    ):
        """
        If present, read df.csv and check for completion
        :param output_d: Path to data directory
        :param df_csv: dataframe of given csv file
        :param hdf_path: Path to h5 file
        :param csv_uids: unique IDs of the csv file
        :param index_name: header of index row
        :return: final dataframe and its column headers
        """
        # load & read df.csv if present
        pres_df = output_d / "df.csv"
        if pres_df.is_file():
            print("Pre computed dataframe file df.csv is loaded.")
            pres_df_csv = pd.read_csv(pres_df, index_col=0, na_filter=False)

            # Check whether no. of rows equals data
            if len(pres_df_csv) != len(df_csv):
                print("# of rows doesn't match data!\nStart recalculation!")
                df_embeddings, csv_header = self._create_df(
                    output_d,
                    hdf_path,
                    csv_uids,
                    df_csv,
                    index_name,
                    umap_flag,
                )
            else:
                # Check each column x, y & z for incompleteness
                if not self._check_coordinates(pres_df_csv):
                    print("Start recalculation!")
                    df_embeddings, csv_header = self._create_df(
                        output_d,
                        hdf_path,
                        csv_uids,
                        df_csv,
                        index_name,
                        umap_flag,
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
                        pres_df_csv.to_csv(output_d / "df.csv", index_label=index_name)

                    # save final column names
                    csv_header = [
                        header
                        for header in pres_df_csv.columns
                        if header not in self.AXIS_NAMES
                        and header != index_name
                        and header != "variance"
                    ]

                    # Unify df name
                    df_embeddings = pres_df_csv

        # create dataframe from files
        else:
            df_embeddings, csv_header = self._create_df(
                output_d,
                hdf_path,
                csv_uids,
                df_csv,
                index_name,
                umap_flag,
            )

        return df_embeddings, csv_header

    def _create_df(
        self,
        output_d: Path,
        hdf_path: Path,
        csv_uids: list[str],
        df_csv: DataFrame,
        index_name: str,
        umap_flag: bool,
    ):
        """
        Use data and create corresponding dataframe
        :param output_d: Path to data directory
        :param hdf_path: Path to h5 file
        :param csv_uids: unique IDs of csv file
        :param df_csv: dataframe of csv file
        :param index_name: header of index column
        :return: complete dataframe and list of its headers
        """
        # read embeddings from HDF5 format
        embeddings = self._get_embeddings(emb_h5file=hdf_path, csv_uids=csv_uids)

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

        # generate dimensionality reduction components and merge it to CSV DataFrame
        if umap_flag:
            df_dim_red = self._generate_umap(embs)
            df_dim_red.index = uids
        else:
            df_dim_red = self._generate_pca(embs)
            df_dim_red.index = uids

        df_embeddings = df_csv.join(df_dim_red, how="right")
        csv_header = [
            header
            for header in df_embeddings.columns
            if header not in self.AXIS_NAMES and header != "variance"
        ]

        # save dataframe
        df_embeddings.to_csv(output_d / "df.csv", index_label=index_name)

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
            n_neighbors=self.umap_paras["n_neighbours"],
            min_dist=self.umap_paras["min_dist"],
            random_state=42,
            n_components=3,
            metric=self.umap_paras["metric"],
        )  # initialize umap; use random_state=42 for reproducibility
        umap_fit = fit.fit_transform(data)  # fit umap to our embeddings
        df_umap = DataFrame(data=umap_fit, columns=self.AXIS_NAMES)
        return df_umap

    def _generate_pca(self, data: np.ndarray):
        from sklearn.decomposition import PCA

        fit = PCA(n_components=3, random_state=42)
        pca_fit = fit.fit_transform(data)
        df_pca = DataFrame(data=pca_fit, columns=self.AXIS_NAMES)

        # extract variance information from pca
        pca_variance = list()
        for variance in fit.explained_variance_ratio_:
            pca_variance.append(variance * 100)

        variance_df = DataFrame({"variance": pca_variance})
        df_pca = pandas.concat([df_pca, variance_df], axis=1)

        return df_pca

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

    @staticmethod
    def init_structure_container(pdb_d: Path):
        structure_container = StructureContainer(pdb_d)

        return structure_container
