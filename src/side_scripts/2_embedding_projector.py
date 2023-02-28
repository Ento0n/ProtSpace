#!/usr/bin/env python3
import argparse
from pathlib import Path
import h5py
import pandas as pd

# Arguemnt parser #

parser = argparse.ArgumentParser()
parser.add_argument("-hdf", required=True, type=str, help="Path to the h5 file.")
parser.add_argument("-o", required=True, type=str, help="Output directory for the tsv file(s)")
parser.add_argument("-csv", required=False, type=str, help="Path to the csv file (Labels).")
args = parser.parse_args()

hdf_path = Path(args.hdf) if args.hdf is not None else None
out_path = Path(args.o) if args.o is not None else None
csv_path = Path(args.csv) if args.csv is not None else None

# Read Files #

embd_list = list()
uid_list = list()
with h5py.File(hdf_path, "r") as hdf:
    for identifier, embeddings in hdf.items():
        embd_list.append(embeddings[:])
        uid_list.append(identifier)

df_emb = pd.DataFrame(embd_list, index=uid_list)

df_csv = pd.read_csv(csv_path, index_col=0)

# Fit files to IDs #

df_csv, df_emb = df_csv.align(df_emb, join="left", axis=0)

df_emb.dropna(how="any", axis=0, inplace=True)

df_emb, df_csv = df_emb.align(df_csv, join="left", axis=0)

# Save files #

embd_tsv_path = out_path / hdf_path.stem
embd_tsv_path = embd_tsv_path.with_suffix(".tsv")
df_emb.to_csv(embd_tsv_path, sep="\t", header=False, index=False)

label_tsv_path = out_path / csv_path.stem
label_tsv_path = label_tsv_path.with_suffix(".tsv")
df_csv.to_csv(label_tsv_path, sep="\t", index_label="index_name")
