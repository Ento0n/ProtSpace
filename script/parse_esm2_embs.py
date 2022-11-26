# !/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on: Thu 24 Nov 2022 12:04:01
Description: Parse ESM2 embeddings from multiple .pt files to a single .h5 file
Usage:       python extract_esm2_embs.py

@author: tsenoner
"""
import argparse
from pathlib import Path

import h5py
import torch


def setup_arguments() -> argparse.Namespace:
    """Defines and parses required and optional arguments for the script"""
    parser = argparse.ArgumentParser(
        description="Merge & reshape PPIHP files output to a single CSV file"
    )

    parser.add_argument(
        "-esm",
        "--esm_dir",
        required=True,
        type=str,
        help="Path to the directory containing the ESM2 embeddings",
    )

    args = parser.parse_args()
    esm_dir = Path(args.esm_dir)
    return esm_dir


def extract_embeddings(esm_dir: Path) -> None:
    """Extract per-protein embeddings from `esm_dir` and merges in an .h5 file
    """
    hdf_file = esm_dir.parent / "emb_esm2.h5"
    with h5py.File(hdf_file, "w") as hdf:
        for pt_file in esm_dir.glob("*.pt"):
            data = torch.load(pt_file)
            label = data["label"]
            emb = data["mean_representations"][36]
            hdf.create_dataset(name=label, data=emb)


def main():
    esm_dir = setup_arguments()
    extract_embeddings(esm_dir)


if __name__ == "__main__":
    main()
