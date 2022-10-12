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
from math import ceil
from pathlib import Path, PosixPath
from typing import Generator

from pyfaidx import Fasta

# import typer


class FastaMapper:
    """Maps the fasta file headers to alphabet like headers (AA AB ... ZY ZZ)"""

    seqs: Fasta

    def __init__(
        self,
        fasta_input: Path,
        fasta_output: Path = None,
        csv_mapping: Path = None,
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
        self.cache = self.fasta_input.with_suffix(".fasta.fai")

    def _fasta_parser(self) -> None:
        # efficient pythonic random access to fasta subsequences
        self.seqs = Fasta(self.fasta_input)

    def _get_new_seqs_ids(self) -> Generator[str, None, None]:
        # get number of sequences and calculate how many letters are needed to
        # uniquely identify all sequences
        nr_seqs = len(self.seqs.keys())
        alphabet = "".join(map(chr, range(ord("A"), ord("Z") + 1)))
        nr_letters_needed = ceil(nr_seqs ** (1 / len(alphabet)))

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

        # create new FASTA and CSV file
        with (
            open(self.fasta_output, "w", encoding="utf-8") as fasta_handler,
            open(self.csv_mapping, "w", encoding="utf-8") as map_handler,
        ):
            map_handler.write("mapped_id,original_id\n")
            for seq, new_seq_id in zip(self.seqs, new_seq_ids):
                fasta_handler.write(f">{new_seq_id}\n")
                fasta_handler.write(f"{seq}\n")
                map_handler.write(f"{new_seq_id},{seq.long_name}\n")
        print(f"Mapped FASTA file saved at: {self.fasta_output}")
        print(f"CSV mapping to old file: {self.csv_mapping}")

        # remove cache
        self.cache.unlink()


def parse_args():
    """Defines required and optional arguments to run the script

    Run example: python fasta_mapping/main.py -fi example/example.fasta
    """

    # Instantiate the parser
    parser = argparse.ArgumentParser(
        description="Map fasta file to unique identifiers e.g. AA AB AC BA ..."
    )
    # Required positional argument
    parser.add_argument(
        "-fi",
        "--fasta_input",
        required=True,
        type=str,
        help="The .fasta file to be mapped",
    )
    # Required positional argument
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
    # Required positional argument
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

    args = parser.parse_args()
    fasta_input = Path(args.fasta_input)
    fasta_output = (
        Path(args.fasta_output) if args.fasta_output is not None else None
    )
    csv_mapping = Path(args.mapping) if args.mapping is not None else None

    return fasta_input, fasta_output, csv_mapping


def main():
    """Run main script"""
    fasta_input, fasta_output, csv_mapping = parse_args()

    fasta_mapper = FastaMapper(
        fasta_input=fasta_input,
        fasta_output=fasta_output,
        csv_mapping=csv_mapping,
    )
    fasta_mapper.create_mapping()


if __name__ == "__main__":
    main()
