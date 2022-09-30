#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys
from pathlib import Path
import h5py
from scipy.spatial.distance import pdist, squareform
import numpy as np
import pandas as pd
import re
import umap
from pandas import DataFrame

import dash
from dash import Input, Output
from dash.exceptions import PreventUpdate

from visualization import render

NON_WORD_RE = re.compile(r"\W")
AXIS_NAMES = ["x", "y", "z"]

df = pd.DataFrame()


def data_preprocessing(data_dir_path, basename, csv_separator, uid_col, html_col):
    # Assign html argument accordingly
    if html_col == -1:
        html_set = False
        html_col = 0
    else:
        html_set = True

    # root directory that holds, proteins.fasta, embeddings.h5, labels.csv and some output_file.html
    root = Path(data_dir_path)

    # load df.csv if present
    pres_df = root / "df.csv"

    if pres_df.is_file():
        pres_df_csv = pd.read_csv(pres_df, sep=csv_separator, index_col=uid_col)

    rep_seqs = root / f"{basename}.fasta"
    emb_h5file = root / f"{basename}.h5"
    label_csv_p = root / f"{basename}.csv"

    fig_3D_p = root / f"{basename}.html"
    fig_2D_p = root / f"{basename}.pdf"

    df_csv = pd.read_csv(label_csv_p, sep=csv_separator, index_col=uid_col)
    csv_uids = unify_seq_uids(
        df_csv.index
    )  # TODO: where else does this operation need to be performed?
    df_csv.index = csv_uids

    # read embeddings from HDF5 format
    embeddings = get_embeddings(emb_h5file=emb_h5file, csv_uids=csv_uids)

    # matrix of values (protein-embeddings); n_proteins x embedding_dim
    uids, embs = zip(*embeddings.items())
    embs = np.vstack(embs)

    # data should be n_proteins x 1024 (ProtT5) OR n_proteins x 128 (ProtTucker)
    print(f"Shape of embeddings (num_proteins x embedding dim): {embs.shape}")

    # get pairwise distances; should be n_proteins x n_proteins
    pdist = squareform(pairwise_distances(embs))
    print(
        "Shape of pairwise distance matrix (num_proteins x num_proteins):"
        f" {pdist.shape}"
    )

    # generate umap components and merge it to CSV DataFrame
    df_umap = generate_umap(embs)
    df_umap.index = uids

    # # --- DUMMY data ---
    # data = np.random.random((len(df_csv), 3))
    # df_umap = pd.DataFrame(data=data, columns=AXIS_NAMES)
    # df_umap.index = csv_uids
    # # --- END ---

    df_embeddings = df_csv.join(df_umap, how="right")
    csv_header = [
        header for header in df_embeddings.columns if header not in AXIS_NAMES
    ]

    # save dataframe
    df_embeddings.to_csv(root / "df.csv")

    # save df_embeddings in global variable df
    global df
    df = df_embeddings

    # generate initial figure
    fig = render(df=df_embeddings, selected_column=csv_header[html_col])

    # save_plotly_figure_to_html(fig, str(fig_3D_p))
    if html_set:
        fig.write_html(root / "plot.html")

    return df_embeddings, fig, csv_header


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
        sys.exit("None of the Unique IDs of the h5 and the csv file matched.")

    print(f"Example: {next(iter(embeddings.keys()))}")
    print(f"Number of embeddings: {len(embeddings)}")
    if (nr_missed := len(missing)) > 0:
        print(f"Lost {nr_missed} proteins due to ID-missmatch b/w CSV & FASTA")
        if nr_missed >= 10:
            print(", ".join(missing))

    return embeddings


def pairwise_distances(data, metric="euclidean"):
    # usually euclidean or cosine distance worked best
    return pdist(data, metric=metric)


def unify_seq_uids(uids: list[str]) -> list[str]:
    return list(map(lambda uid: NON_WORD_RE.sub("_", uid), uids))


def generate_umap(data: np.ndarray) -> pd.DataFrame:
    # visualize high-dimensional embeddings with dimensionality reduction (here: umap)
    # Tutorial: https://umap-learn.readthedocs.io/en/latest/basic_usage.html
    # Parameters: https://umap-learn.readthedocs.io/en/latest/parameters.html
    fit = umap.UMAP(
        n_neighbors=25, min_dist=0.5, random_state=42, n_components=3
    )  # initialize umap; use random_state=42 for reproducability
    umap_fit = fit.fit_transform(data)  # fit umap to our embeddings
    df_umap = DataFrame(data=umap_fit, columns=AXIS_NAMES)
    return df_umap


@dash.callback(
    Output("store_data", "data"), Output("dd_menu", "value"), Input("dd_menu", "value")
)
def store_data(value):
    # Check whether an input is triggered
    ctx = dash.callback_context
    if not ctx.triggered:
        raise PreventUpdate

    return df.to_dict("records"), value
