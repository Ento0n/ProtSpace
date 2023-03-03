#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import re
from pathlib import Path


class StructureContainer:
    def __init__(self, pdb_d: object, json_d: object):
        self.pdb_d = pdb_d
        self.json_d = json_d
        if pdb_d is not None:
            self.pdb_flag = True
        else:
            self.pdb_flag = False
        if json_d is not None:
            self.json_flag = True
        else:
            self.json_flag = False

    def get_structure_dir(self):
        return self.pdb_d

    def set_json_dir(self, json_d: Path):
        self.json_d = json_d

    def get_json_dir(self) -> Path:
        if isinstance(self.json_d, Path):
            return self.json_d
        else:
            return Path("")

    def get_range(self, uid: str):
        if isinstance(self.pdb_d, Path):
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
