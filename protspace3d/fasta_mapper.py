# !/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on:  Thu 15 Sep 2022 17:08:06
Description: Create new fasta file with short unique ids and a csv file
             containing the mapping back to the original file
Usage:       python main.py

@author: tsenoner
"""
import argparse
from itertools import product
from math import log, ceil
from pathlib import Path
from typing import Generator

import pandas as pd
from pyfaidx import Fasta

# import typer


class FastaMapper:
    """Maps the fasta file headers to alphabet like headers (AA AB ... ZY ZZ)"""

    seqs: Fasta
    df_features: pd.DataFrame = None

    def __init__(
        self,
        fasta_input: Path,
        fasta_output: Path = None,
        csv_mapping: Path = None,
        csv_features: Path = None,
    ) -> None:
        fasta_suffixes = [".fasta", ".fna", ".ffn", ".faa", ".frn", ".fa"]
        assert fasta_input.suffix in fasta_suffixes, (
            "File does not end with common FASTA suffixes:"
            f" {' '.join(fasta_suffixes)}"
        )
        self.fasta_input = fasta_input
        self.fasta_output = (
            fasta_output
            if fasta_output is not None
            else fasta_input.with_stem(f"{fasta_input.stem}_mapped")
        )
        self.csv_mapping = (
            csv_mapping
            if csv_mapping is not None
            else self.fasta_output.with_suffix(".csv")
        )
        if csv_features is not None:
            self.df_features = pd.read_csv(csv_features, sep=";")
            self.df_features.iloc[:, 0] = self.df_features.iloc[:, 0].astype(
                "str"
            )
        self.cache = self.fasta_input.with_suffix(".fasta.fai")

    def _fasta_parser(self) -> None:
        # efficient pythonic random access to fasta subsequences
        self.seqs = Fasta(
            self.fasta_input,
            read_long_names=True,
            sequence_always_upper=True,
        )

    def _get_new_seqs_ids(self) -> Generator[str, None, None]:
        # get number of sequences and calculate how many letters are needed to
        # uniquely identify all sequences
        nr_seqs = len(self.seqs.keys())
        alphabet = "".join(map(chr, range(ord("A"), ord("Z") + 1)))
        nr_letters_needed = ceil(log(nr_seqs, len(alphabet)))

        # new identifier generator (e.g. AA, AB, AC, ...)
        new_seq_ids = (
            "".join(combi)
            for combi in product(alphabet, repeat=nr_letters_needed)
        )
        return new_seq_ids

    def create_mapping(self) -> None:
        """Create the new FASTA and CSV mapping file

        - FASTA file contains the new headers (uids) and the same sequences
        - CSV mapping file contains the columns: `mapped_id` & `original_id`
        """
        self._fasta_parser()
        new_seq_ids = self._get_new_seqs_ids()

        # create new FASTA and CSV
        data = {"mapped_id": [], "original_id": []}
        with (open(self.fasta_output, "w", encoding="utf-8") as fasta_handler,):
            for seq, new_seq_id in zip(self.seqs, new_seq_ids):
                fasta_handler.write(f">{new_seq_id}\n")
                fasta_handler.write(f"{seq}\n")
                data["mapped_id"].append(new_seq_id)
                data["original_id"].append(seq.long_name)
        df = pd.DataFrame(data)
        if self.df_features is not None:
            first_column_name = self.df_features.columns[0]
            df = df.merge(
                right=self.df_features,
                left_on="original_id",
                right_on=first_column_name,
                how="outer",
                copy=False,
            )
            df = df.drop(columns=[first_column_name])
        df.to_csv(self.csv_mapping, index=False)
        print(f"Mapped FASTA file saved at: {self.fasta_output}")
        print(f"CSV mapping to old file: {self.csv_mapping}")

        # remove cache
        self.cache.unlink()


def parse_args():
    """Defines required and optional arguments to run the script

    Run example: python protspace3d/fasta_mapper.py -fi data/ex3/Conotoxins_try1.fasta -ff data/ex3/Conotoxins_try1.csv
    """

    # Instantiate the parser
    parser = argparse.ArgumentParser(
        description="Map fasta file to unique identifiers e.g. AA AB AC BA ..."
    )
    parser.add_argument(
        "-fi",
        "--fasta_input",
        required=True,
        type=str,
        help="The .fasta file to be mapped",
    )
    parser.add_argument(
        "-fo",
        "--fasta_output",
        required=False,
        type=str,
        help=(
            "The .fasta output path of the file (default: same as fasta_input"
            " with `_mapped` added to the name)"
        ),
    )
    parser.add_argument(
        "-m",
        "--mapping",
        required=False,
        type=str,
        help=(
            "CSV path to store mapping. (default: same as fasta_output with"
            " .csv as suffix). Column names: mapped_id, original_id"
        ),
    )
    parser.add_argument(
        "-ff",
        "--feature_file",
        required=False,
        type=str,
        help=(
            "CSV file containing additional information to group the"
            " embeddings. If given, the features will be added to the CSV file."
        ),
    )

    args = parser.parse_args()
    fasta_input = Path(args.fasta_input)
    fasta_output = (
        Path(args.fasta_output) if args.fasta_output is not None else None
    )
    csv_mapping = Path(args.mapping) if args.mapping is not None else None
    csv_features = (
        Path(args.feature_file) if args.feature_file is not None else None
    )

    return fasta_input, fasta_output, csv_mapping, csv_features


def main():
    """Run main script"""
    fasta_input, fasta_output, csv_mapping, csv_features = parse_args()
    fasta_mapper = FastaMapper(
        fasta_input=fasta_input,
        fasta_output=fasta_output,
        csv_mapping=csv_mapping,
        csv_features=csv_features,
    )
    fasta_mapper.create_mapping()


if __name__ == "__main__":
    main()
