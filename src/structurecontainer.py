#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from pathlib import Path
import re


class StructureContainer(object):
    def __init__(self, pdb_d):
        self.pdb_d = pdb_d
        self.json_d = None
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

    def has_json_dir(self):
        if self.json_d is None:
            return False
        else:
            return True

    def set_json_dir(self, json_d: Path):
        self.json_d = json_d

    def get_json_dir(self) -> Path:
        return self.json_d

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
