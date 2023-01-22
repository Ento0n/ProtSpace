# !/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on: Wed 18 Jan 2023 02:27:10
Description: Parse 3FTx dset and get uniprot ID by extraction or BLAST
Usage:       python 3FTx_dset.py

@author: tsenoner
"""
import argparse
import json
import re
from pathlib import Path

import pandas as pd
from pyfaidx import Fasta
from uniprot_helper import BatchBlaster, get_tax_id


def setup_arguments() -> argparse.Namespace:
    """Defines and parses required and optional arguments for the script"""
    parser = argparse.ArgumentParser(
        description="Merge & reshape PPIHP files output to a single CSV file"
    )

    parser.add_argument(
        "-ci",
        "--csv_file_in",
        required=True,
        type=str,
        help="Path to CSV file",
    )
    parser.add_argument(
        "-fi",
        "--fasta_file_in",
        required=True,
        type=str,
        help="Path to FASTA file",
    )
    parser.add_argument(
        "-co",
        "--csv_file_out",
        required=True,
        type=str,
        help="Path to CSV output file",
    )
    parser.add_argument(
        "-fo",
        "--fasta_file_out",
        required=True,
        type=str,
        help="Path to FASTA output file",
    )
    parser.add_argument(
        "-t",
        "--taxon_mapper",
        required=True,
        type=str,
        help="Path to JSON file which maps species to taxon ids",
    )
    parser.add_argument(
        "-b",
        "--blast_dir",
        required=True,
        type=str,
        help="Path to directory where BLASTed JSON files are saved",
    )

    args = parser.parse_args()
    csv_in = Path(args.csv_file_in)
    fasta_in = Path(args.fasta_file_in)
    csv_out = Path(args.csv_file_out)
    fasta_out = Path(args.fasta_file_out)
    taxon_mapper_file = Path(args.taxon_mapper)
    blast_dir = Path(args.out_dir)
    return csv_in, fasta_in, csv_out, fasta_out, taxon_mapper_file, blast_dir


def construct_df(csv_path, fasta_path):
    df = pd.read_csv(csv_path, sep=",")
    # remove rows with `snake genomic` as `Major group``
    # df = df.loc[~df["Major group"].isin(["snake genomic"]), :]

    # rename column
    df = df.rename(
        columns={
            "Embedding ID": "embedding_id",
            "Evolutionary order": "evolutionary_order",
            "Major group": "major_group",
            "Major taxon for the purposes of this study": "major_taxon",
            "Name in fasta": "fasta_id",
            "Name (existing or suggested)": "name",
            "Original fasta header": "original_fasta_header",
            "Preliminary cysteine group": "cystein_group",
            "Species": "species",
        }
    )

    # select columns to keep
    cols2keep = [
        "cystein_group",
        # "embedding_id",
        "evolutionary_order",
        "fasta_id",
        "major_group",
        "major_taxon",
        "name",
        "original_fasta_header",
        "species",
    ]
    df = df[cols2keep]

    # read FASTA file and add sequences to DataFrame
    seqs = (
        (header, str(seq).replace("-", ""))
        for header, seq in Fasta(str(fasta_path)).items()
    )
    df_seq = pd.DataFrame(seqs, columns=["fasta_id", "seq"])
    df = pd.merge(left=df, right=df_seq, on="fasta_id")

    # remove entry with same identifier but different/updated sequence
    df = df.drop(df[df["fasta_id"] == "Walterinnesia_aegyptia_C0HKZ8"].index)

    return df


def create_taxon_mapper(taxas):
    taxon_mapper = {}
    for taxa in taxas:
        if taxa not in taxon_mapper:
            taxa_id = get_tax_id(taxa)
            if taxa_id is None:
                raise Exception(f"'{taxa}' not found")
            taxon_mapper[taxa] = taxa_id
    return taxon_mapper


def add_taxon_id(df: pd.DataFrame, taxon_mapper_file: Path) -> pd.DataFrame:
    # correct wrong taxa
    df.loc[df["species"] == "Micrurus_ tener", "species"] = "Micrurus tener"

    # read in already mapped taxon
    if taxon_mapper_file.is_file():
        with open(taxon_mapper_file, "r") as json_file:
            taxon_mapper = json.load(json_file)
    else:
        taxon_mapper = dict()
    unknown_taxon_ids = df.loc[
        ~df["species"].isin(taxon_mapper.keys()), "species"
    ].unique()

    # get taxon_id that arn't in `taxon_mapper` yet
    novel_taxon_mapper = create_taxon_mapper(taxas=unknown_taxon_ids)
    taxon_mapper.update(novel_taxon_mapper)

    # update json file
    with open(taxon_mapper_file, "w") as json_file:
        json.dump(taxon_mapper, fp=json_file, indent=4)

    # add taxon ids to DataFrame
    df["taxon_id"] = df["species"].map(taxon_mapper)
    return df


def get_uniprot_acc_ids(df) -> tuple[list, dict]:
    # get UniProt ids
    df["acc_id"] = None
    for idx, row in df.iterrows():
        uid = "blank"
        header = row["original_fasta_header"]
        # 1. check for patterns that have no UniProt ID
        # - `3FTx_\d{3}`. E.g.: Cbivi_3FTx_000
        match1 = re.match(pattern=r".*(3FTx_\d{2,3})", string=header)
        # -`unigene\d*`. E.g.: Heterodon_nasicus_unigene14895
        match2 = re.match(pattern=r".*(unigene\d{4,6})", string=header)
        if match1 or match2:
            uid = None
        # 2. check for regular UniProt entry
        elif header.startswith("sp") or header.startswith("tr"):
            uid = header.split("|")[1]
        # 3. check for odd manually descriptive naming
        elif " " in header:
            uid = None
        # 4. check if last element of sequence is UID
        else:
            header_last_part = header.split("_")[-1]
            # UniProt regular expression for accession nummers: https://www.uniprot.org/help/accession_numbers
            uniprot_accession_pattern = r"[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}"
            match_accession = re.match(
                pattern=uniprot_accession_pattern, string=header_last_part
            )
            if match_accession:
                uid = header_last_part
            if not match_accession:
                uid = None
        df.loc[idx, "acc_id"] = uid
    print(f"- found UIDs: {df['acc_id'].count()}")
    print(f"- unknown UIDs: {df['acc_id'].isna().sum()}")
    return df


def add_uniprot_acc_ids(df: pd.DataFrame, blast_dir: Path) -> pd.DataFrame:
    """Get Accession ID from BLASTp with 100% match. Favor swissprot entries.

    Gets Accession ID for all entries that have None as a value in `acc_id`.

    Args:
        df (pd.DataFrame): columns: fasta_id, acc_id
        out_dir (Path): directory containing JSON files from BLASTp results

    Returns:
        pd.DataFrame: updated accession number
    """
    df["db"] = None
    for idx, row in df.iterrows():
        # if row["acc_id"] is not None:
        #     continue

        # read JSON file
        json_file = blast_dir / f"{row['fasta_id']}.json"
        with open(json_file, "r") as json_handler:
            json_data = json.load(json_handler)

        # get accession IDs
        acc_ids = []
        dbs = []
        query_len = json_data["query_len"]
        for hits in json_data["hits"]:
            hit_acc = hits["hit_acc"]
            hit_db = hits["hit_db"]
            for hit in hits["hit_hsps"]:
                if (
                    hit["hsp_identity"] == 100.0
                    and hit["hsp_align_len"] == query_len
                ):
                    acc_ids.append(hit_acc)
                    dbs.append(hit_db)

        # if SP entry remove all other entries
        if "SP" in dbs:
            indices = [idx for idx, db in enumerate(dbs) if db == "SP"]
            acc_ids = [acc_ids[idx] for idx in indices]
            db = "SP"
        elif "TR" in dbs:
            db = "TR"
        else:
            db = None
        # add previous extracted Acc_ID to list
        if (row["acc_id"] is not None) and (row["acc_id"] not in acc_ids):
            acc_ids.insert(0, row["acc_id"])

        # add accession numbers and DB
        df.loc[idx, "acc_id"] = ",".join(acc_ids) if acc_ids else None
        df.loc[idx, "db"] = db
    return df

def run_blast(df: pd.DataFrame, blast_dir: Path) -> pd.DataFrame:
    # entries2blast = df.loc[
    #     df["acc_id"].isna(), ["fasta_id", "seq", "taxon_id"]
    # ].to_dict("records")
    entries2blast = df.loc[:, ["fasta_id", "seq", "taxon_id"]].to_dict(
        "records"
    )
    ncbi_blaster = BatchBlaster(entries=entries2blast, out_dir=blast_dir)
    ncbi_blaster.run_batch()
    df = add_uniprot_acc_ids(df=df, blast_dir=blast_dir)
    return df

def manual_curation(df: pd.DataFrame) -> pd.DataFrame:
    # remove unneded column
    df = df.drop(columns="original_fasta_header")
    # remove duplicates
    df = df.drop_duplicates(subset=["acc_id", "species", "seq"])

    # find sequences that come from genomic sequences
    # pattern: ^\d+[A-Za-z]{4}[_\d]
    genomic_id_cond = df["fasta_id"].str.match(r"^\d+[A-Za-z]{4}[_\d]") == True
    df.loc[genomic_id_cond, "data_origin"] = "genomic"
    # df.loc[genomic_id_cond, "fasta_id"] = df.loc[genomic_id_cond, "fasta_id"].str.lstrip("0123456789")
    return df


def save_data(df: pd.DataFrame, csv_file: Path, fasta_file: Path) -> None:
    # move column `fasta_id` to the front
    col_data = df.pop("fasta_id")
    df.insert(loc=0, column="fasta_id", value=col_data)

    df.to_csv(csv_file, index=False)
    with open(fasta_file, "w") as fasta_handler:
        for _, row in df.iterrows():
            fasta_handler.write(f">{row['fasta_id']}\n")
            fasta_handler.write(f"{row['seq']}\n")



def main():
    csv_in, fasta_in, csv_out, fasta_out, taxon_mapper_file, blast_dir = setup_arguments()

    df = construct_df(csv_path=csv_in, fasta_path=fasta_in)
    df = add_taxon_id(df=df, taxon_mapper_file=taxon_mapper_file)
    df = get_uniprot_acc_ids(df=df)
    df = run_blast(df=df, blast_dir=blast_dir)
    df = manual_curation(df=df)
    save_data(df=df, csv_file=csv_out, fasta_file=fasta_out)


if __name__ == "__main__":
    main()
