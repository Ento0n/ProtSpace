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
    parser.add_argument(
        "-o",
        "--out_hdf",
        required=False,
        type=str,
        help="Path to the HDF5 output file",
    )

    args = parser.parse_args()
    esm_dir = Path(args.esm_dir)
    out_hdf = Path(args.out_hdf) if args.out_hdf is not None else esm_dir.parent / "emb_esm2.h5"
    return esm_dir, out_hdf


def extract_embeddings(esm_dir: Path, out_hdf: Path) -> None:
    """Extract per-protein embeddings from `esm_dir` and merges in an .h5 file
    """
    with h5py.File(out_hdf, "w") as hdf:
        for pt_file in esm_dir.glob("*.pt"):
            data = torch.load(pt_file)
            label = data["label"]
            emb = data["mean_representations"][36]
            hdf.create_dataset(name=label, data=emb)


def main():
    esm_dir, out_hdf = setup_arguments()
    extract_embeddings(esm_dir=esm_dir, out_hdf=out_hdf)


if __name__ == "__main__":
    main()
