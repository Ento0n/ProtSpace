# !/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on: Mon 17 Oct 2022 19:40:19
Description: extract best predicted PDB models from ColabFold output
Usage:       python get_best_pdbs.py

@author: tsenoner
"""
import argparse
import shutil
from pathlib import Path, PosixPath


def setup_arguments() -> argparse.Namespace:
    """Defines required and optional arguments for the script"""
    # Instantiate the parser
    parser = argparse.ArgumentParser(
        description="Extract PDB models and their score files"
    )

    # Required argument
    parser.add_argument(
        "-cf",
        "--colabfold",
        required=True,
        type=str,
        help="A path to the data directory output generated from ColbFold",
    )

    args = parser.parse_args()
    colabfold_dir = Path(args.colabfold)
    return colabfold_dir


def get_proten_header(pdb_fn: PosixPath) -> str:
    """Get the protein header which can be found in the `.a3m` file"""
    file_nr = pdb_fn.stem.split("_")[0]
    a3m_fn = pdb_fn.parent / f"{file_nr}.a3m"
    with open(a3m_fn, "r") as a3m_handler:
        # print(a3m_handler.readlines())
        for line in a3m_handler:
            line = line.strip()
            if line.startswith("#"):
                continue
            if line.startswith(">"):
                header = line[1:]
                break
        else:
            raise ValueError(f"File {a3m_fn} is missformatted.")
    return header


def extract_rank1_models(from_dir: PosixPath, pdb_dir: PosixPath, json_dir: PosixPath) -> None:
    """Copy rank_1 files from `from_dir` to `to_dir` with header from `.a3m`

    *rank_1* files contain the best of the 5 predicted models. The `.pbd` file
    contains the 3D structure of the model itself whereas the `.json` file
    contains the scores (PAE, pLDDT, pTM).
    """
    for directory in [pdb_dir, json_dir]:
        if not directory.is_dir():
            directory.mkdir()

    rank1_files = from_dir.glob("*rank_1*")
    for pdb_file in rank1_files:
        header = get_proten_header(pdb_file)
        suffix = pdb_file.suffix
        if suffix == ".pdb":
            new_fn = pdb_dir / f"{header}{pdb_file.suffix}"
        elif suffix == ".json":
            new_fn = json_dir / f"{header}{pdb_file.suffix}"
        else:
            raise ValueError(f"Suffix type `{suffix}` is not valide.")
        # print(new_pdb_fn)
        if not new_fn.is_file():
            shutil.copy(pdb_file, new_fn)


def main():
    colabfold_dir = setup_arguments()
    predictions = colabfold_dir / "predictions"
    pdb_dir = colabfold_dir.parent / "pdb"
    json_dir = colabfold_dir.parent / "json"

    extract_rank1_models(from_dir=predictions, pdb_dir=pdb_dir, json_dir=json_dir)


if __name__ == "__main__":
    main()
