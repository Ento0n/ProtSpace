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
from Bio import SeqIO

from visualization.visualizator import Visualizator
from structurecontainer import StructureContainer


class DataPreprocessor:
    UMAP_AXIS_NAMES = ["x_umap", "y_umap", "z_umap"]
    PCA_AXIS_NAMES = ["x_pca", "y_pca", "z_pca"]
    AXIS_NAMES = ["x_umap", "y_umap", "z_umap", "x_pca", "y_pca", "z_pca"]

    def __init__(
        self,
        output_d: Path,
        hdf_path: Path,
        csv_path: Path,
        fasta_path: Path,
        csv_separator: str,
        uid_col: int,
        html_cols: list[str],
        reset: bool,
        umap_paras: dict,
    ):
        self.output_d = output_d
        self.hdf_path = hdf_path
        self.csv_path = csv_path
        self.fasta_path = fasta_path
        self.csv_separator = csv_separator
        self.uid_col = uid_col
        self.html_cols = html_cols
        self.reset = reset
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
            df_csv_path = self.output_d / f"df_{self.hdf_path.stem}.csv"
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

        # replace "None" values with NA
        df_csv.replace(to_replace="None", value="NA", inplace=True)

        # save index name for df.csv
        index_name = df_csv.index.name

        # get UIDs
        csv_uids = df_csv.index.to_list()

        # get sequences if fasta path is given
        fasta_dict = None
        if self.fasta_path is not None:
            fasta_dict = self._read_fasta()

        df_embeddings, csv_header, embeddings, embedding_uids = self._read_df_csv(
            self.output_d,
            df_csv,
            emb_h5file,
            csv_uids,
            index_name,
        )

        # handle html saving
        DataPreprocessor._handle_html(
            self,
            self.html_cols,
            csv_header,
            self.output_d,
            original_id_col=original_id_col,
            df=df_embeddings,
        )

        # sort csv header alphabetically
        csv_header.sort(key=str.lower)

        # generate initial figure
        fig = Visualizator.render(
            df_embeddings,
            selected_column=csv_header[0],
            original_id_col=original_id_col,
            umap_paras=self.umap_paras,
        )

        return (
            df_embeddings,
            fig,
            csv_header,
            original_id_col,
            embeddings,
            embedding_uids,
            fasta_dict,
        )

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

    def _read_fasta(self):
        fasta_sequences = SeqIO.parse(open(self.fasta_path), "fasta")
        fasta_dict = SeqIO.to_dict(fasta_sequences)
        return fasta_dict

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
        html_cols: list,
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
            # parse html_cols array from str to int if all items are numeric
            for idx, item in enumerate(html_cols):
                if item.isnumeric():
                    html_cols[idx] = int(item)

            # -1 indicates all columns to be saved
            if html_cols == [-1]:
                for col in csv_header:
                    fig = Visualizator.render(
                        df,
                        selected_column=col,
                        original_id_col=original_id_col,
                        umap_paras=self.umap_paras,
                    )
                    fig.write_html(output_d / f"3Dspace_{col}.html")

            else:
                # differentiate between given column names and indexes
                if all([isinstance(item, str) for item in html_cols]):
                    try:
                        for item in html_cols:
                            fig = Visualizator.render(
                                df,
                                selected_column=item,
                                original_id_col=original_id_col,
                                umap_paras=self.umap_paras,
                            )
                            fig.write_html(output_d / f"3Dspace_{item}.html")
                    except Exception:
                        raise Exception(
                            f"Given column name <{item}> for html output don't match column names in csv file!"
                            + f"\npossible selection: {csv_header}"
                        )

                elif all([isinstance(item, int) for item in html_cols]):
                    # Sort given column indexes
                    html_cols = sorted(html_cols)

                    possible_selection = list(range(1, len(csv_header) + 1))
                    for col in html_cols:
                        if col not in possible_selection:
                            raise Exception(
                                f"Given html column <{col}> is not valid!"
                                + f"\npossible range: {possible_selection}"
                            )

                    # Edit input to required index numbers
                    for i, num in enumerate(html_cols):
                        html_cols[i] = num - 1

                    # Given parameters existing columns?
                    for col in html_cols:
                        fig = Visualizator.render(
                            df,
                            selected_column=csv_header[col],
                            original_id_col=original_id_col,
                            umap_paras=self.umap_paras,
                        )
                        fig.write_html(output_d / f"3Dspace_{csv_header[col]}.html")
                # Mixed types in the html columns list, can't process this
                else:
                    raise Exception(
                        "Mixed input for html columns!\nEither column name or index is a valid input!"
                    )

    def _read_df_csv(
        self,
        output_d: Path,
        df_csv: DataFrame,
        hdf_path: Path,
        csv_uids: list[str],
        index_name: str,
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
        # Create embeddings
        embs = self._get_embeddings(hdf_path, csv_uids)
        embedding_uids, embs = zip(*embs.items())
        embeddings = np.vstack(embs)

        # load & read df.csv if present
        pres_df = output_d / f"df_{hdf_path.stem}.csv"
        if pres_df.is_file():
            print(f"Pre computed dataframe file df_{hdf_path.stem}.csv is loaded.")
            pres_df_csv = pd.read_csv(pres_df, index_col=0)
            pres_df_csv.fillna("NA", inplace=True)

            # Check whether no. of rows equals data
            if len(pres_df_csv) != len(df_csv):
                print("# of rows doesn't match data!\nStart recalculation!")
                df_embeddings, csv_header, embeddings = self._create_df(
                    output_d,
                    hdf_path,
                    csv_uids,
                    df_csv,
                    index_name,
                    embeddings,
                    embedding_uids,
                )
            else:
                # Check each column x, y & z for incompleteness
                if not self._check_coordinates(pres_df_csv):
                    print("Start recalculation!")
                    df_embeddings, csv_header, embeddings = self._create_df(
                        output_d,
                        hdf_path,
                        csv_uids,
                        df_csv,
                        index_name,
                        embeddings,
                        embedding_uids,
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
                        pres_df_csv.to_csv(
                            output_d / f"df_{hdf_path.stem}.csv", index_label=index_name
                        )

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
            df_embeddings, csv_header, embeddings = self._create_df(
                output_d,
                hdf_path,
                csv_uids,
                df_csv,
                index_name,
                embeddings,
                embedding_uids,
            )

        return df_embeddings, csv_header, embeddings, embedding_uids

    def _create_df(
        self,
        output_d: Path,
        hdf_path: Path,
        csv_uids: list[str],
        df_csv: DataFrame,
        index_name: str,
        embs: np.ndarray,
        embs_uids: list[str],
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

        # check for proteins in csv but not in h5 file
        self._check_csv_uids(embeddings_uids=embs_uids, csv_uids=csv_uids)

        # data should be n_proteins x 1024 (ProtT5) OR n_proteins x 128 (ProtTucker)
        print(f"Shape of embeddings (num_proteins x embedding dim): {embs.shape}")

        # get pairwise distances; should be n_proteins x n_proteins
        pairwise_dist = squareform(self._pairwise_distances(embs))
        print(
            "Shape of pairwise distance matrix (num_proteins x num_proteins):"
            f" {pairwise_dist.shape}"
        )

        # generate dimensionality reduction components and merge it to CSV DataFrame
        df_dim_red_umap = self.generate_umap(embs, self.umap_paras)
        df_dim_red_umap.index = embs_uids
        df_dim_red_pca = self._generate_pca(embs)
        df_dim_red_pca.index = embs_uids

        df_embeddings = df_csv.join([df_dim_red_umap, df_dim_red_pca], how="outer")
        csv_header = [
            header
            for header in df_embeddings.columns
            if header not in self.AXIS_NAMES and header != "variance"
        ]

        # save dataframe
        df_embeddings.to_csv(
            output_d / f"df_{hdf_path.stem}.csv", index_label=index_name
        )

        return df_embeddings, csv_header, embs

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

    @staticmethod
    def generate_umap(data: np.ndarray, umap_paras: dict) -> pd.DataFrame:
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
            n_neighbors=umap_paras["n_neighbours"],
            min_dist=umap_paras["min_dist"],
            random_state=42,
            n_components=3,
            metric=umap_paras["metric"],
        )  # initialize umap; use random_state=42 for reproducibility
        umap_fit = fit.fit_transform(data)  # fit umap to our embeddings
        df_umap = DataFrame(data=umap_fit, columns=["x_umap", "y_umap", "z_umap"])
        return df_umap

    def _generate_pca(self, data: np.ndarray):
        from sklearn.decomposition import PCA

        fit = PCA(n_components=3, random_state=42)
        pca_fit = fit.fit_transform(data)
        df_pca = DataFrame(data=pca_fit, columns=self.PCA_AXIS_NAMES)

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
            print("Df corrupted as not all x,y & z columns are present!")
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
    def _check_csv_uids(embeddings_uids: list[str], csv_uids: list[str]):
        """
        Check unique IDs in csv but not in h5 file
        :param embeddings: data of the h5 file
        :param csv_uids: unique IDs of the csv file
        """
        missing = list()

        # iterate over csv uids
        for uid in csv_uids:

            if uid not in embeddings_uids:
                missing.append(uid)

        if (nr_missed := (len(missing))) > 0:
            print(f"{nr_missed} protein(s) in csv but not in h5 file:")
            print(", ".join(missing))

    def get_umap_paras_dict(self, df: DataFrame):
        umap_paras_dict = dict()
        n_neighbours = self.umap_paras["n_neighbours"]
        min_dist = self.umap_paras["min_dist"]
        metric = self.umap_paras["metric"]
        umap_paras_string = str(n_neighbours) + " ; " + str(min_dist) + " ; " + metric
        coords_dict = dict(
            x_umap=df["x_umap"].to_list(),
            y_umap=df["y_umap"].to_list(),
            z_umap=df["z_umap"].to_list(),
        )
        umap_paras_dict[umap_paras_string] = coords_dict

        return umap_paras_dict
